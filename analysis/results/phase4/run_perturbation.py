"""
Phase 4: Proper Embedding-Based Perturbation Analysis (FIXED)
Uses Geneformer transformer embeddings with correct gene symbol -> Ensembl -> token mapping
"""
import os, sys, gc, pickle, json
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.spatial.distance import cosine
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

MODELS_DIR = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
OUTPUT_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase4"
DATA_PATH = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
MODEL_PATH = os.path.join(MODELS_DIR, "Geneformer-V1-10M")

TARGET_GENES = [
    'CADM2','NRXN1','NLGN1','NEGR1','TCF4','MEF2C','APP','GRIN2A','GRIN2B',
    'MAPT','SCN1A','DISC1','FOXO3','CAMK2A','HOMER1','CREB1','COMT','SHANK3',
    'FMR1','SYN1','BDNF','SLC6A4'
]

MAX_CELLS = 1000
BATCH_SIZE = 25
N_PERMUTATIONS = 30
SKIP_CELLS = 0  # processing fresh 1000, will combine with existing later
MAX_TOKENS = 2048

def load_dictionaries():
    """Load token and gene name dictionaries"""
    dict_dir = os.path.join(MODELS_DIR, "geneformer")
    with open(os.path.join(dict_dir, "token_dictionary_gc104M.pkl"), 'rb') as f:
        token_dict = pickle.load(f)
    with open(os.path.join(dict_dir, "gene_name_id_dict_gc104M.pkl"), 'rb') as f:
        name_dict = pickle.load(f)
    with open(os.path.join(dict_dir, "gene_median_dictionary_gc104M.pkl"), 'rb') as f:
        median_dict = pickle.load(f)
    return token_dict, name_dict, median_dict

def build_gene_mapping(name_dict, token_dict):
    """Build symbol -> token_id mapping for target genes"""
    mapping = {}
    for gene in TARGET_GENES:
        ensembl = name_dict.get(gene)
        if ensembl and ensembl in token_dict:
            mapping[gene] = {'ensembl': ensembl, 'token_id': token_dict[ensembl]}
    print(f"Mapped {len(mapping)}/{len(TARGET_GENES)} target genes to tokens")
    return mapping

def tokenize_cell(gene_expr_dict, ensembl_expr, token_dict, median_dict, max_tokens=2048):
    """
    Tokenize a cell: rank genes by expression/median ratio, convert to tokens
    gene_expr_dict: {ensembl_id: expression_value}
    """
    # Filter to genes in token dictionary with nonzero expression
    valid = {g: e for g, e in ensembl_expr.items() if g in token_dict and e > 0}
    if not valid:
        return []
    
    # Rank by expression normalized by median (Geneformer approach)
    ranked = []
    for gene, expr in valid.items():
        median = median_dict.get(gene, 1.0)
        if median > 0:
            ranked.append((gene, expr / median))
        else:
            ranked.append((gene, expr))
    
    ranked.sort(key=lambda x: x[1], reverse=True)
    
    # Convert to token IDs (top genes first, up to max_tokens - 2 for special tokens)
    token_ids = [token_dict['<cls>']]  # CLS token
    for gene, _ in ranked[:max_tokens - 2]:
        token_ids.append(token_dict[gene])
    token_ids.append(token_dict['<eos>'])  # EOS token
    
    return token_ids

def get_embedding(model, token_ids, device):
    """Get CLS embedding from model"""
    if len(token_ids) < 3:
        return None
    input_ids = torch.tensor(token_ids, dtype=torch.long).unsqueeze(0).to(device)
    with torch.no_grad():
        outputs = model(input_ids)
        cls_emb = outputs.last_hidden_state[0, 0, :].cpu().numpy()
    return cls_emb

def perturb_gene(model, token_ids, target_token_id, device, n_perms=50):
    """Remove target gene token, measure cosine distance, build null distribution"""
    if target_token_id not in token_ids:
        return None
    
    # Baseline embedding
    baseline = get_embedding(model, token_ids, device)
    if baseline is None:
        return None
    
    # Perturbed: remove target gene
    perturbed_tokens = [t for t in token_ids if t != target_token_id]
    perturbed_emb = get_embedding(model, perturbed_tokens, device)
    if perturbed_emb is None:
        return None
    
    cos_dist = cosine(baseline, perturbed_emb)
    
    # Null distribution: remove random non-target, non-special tokens
    special = {token_ids[0], token_ids[-1], target_token_id}  # CLS, EOS, target
    available = [t for t in token_ids if t not in special]
    
    if len(available) == 0:
        return {'distance': cos_dist, 'null_distances': [], 'p_value': np.nan}
    
    n_sample = min(n_perms, len(available))
    random_tokens = np.random.choice(available, n_sample, replace=False)
    
    null_dists = []
    for rt in random_tokens:
        null_tokens = [t for t in token_ids if t != rt]
        null_emb = get_embedding(model, null_tokens, device)
        if null_emb is not None:
            null_dists.append(cosine(baseline, null_emb))
    
    if null_dists:
        p_value = (np.sum(np.array(null_dists) >= cos_dist) + 1) / (len(null_dists) + 1)
    else:
        p_value = np.nan
    
    return {'distance': cos_dist, 'null_distances': null_dists, 'p_value': p_value, 'null_mean': np.mean(null_dists) if null_dists else np.nan}

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # Load dictionaries
    token_dict, name_dict, median_dict = load_dictionaries()
    gene_mapping = build_gene_mapping(name_dict, token_dict)
    
    # Build reverse: symbol -> ensembl for all genes (for tokenization)
    symbol_to_ensembl = {v: k for k, v in name_dict.items()}  # ensembl -> symbol, we need reverse
    # Actually name_dict is symbol -> ensembl
    
    # Load model
    from transformers import BertModel
    model = BertModel.from_pretrained(MODEL_PATH)
    model.to(device)
    model.eval()
    print(f"Model loaded: {MODEL_PATH}")
    
    # Load data
    adata = sc.read_h5ad(DATA_PATH)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    
    if adata.n_obs > MAX_CELLS:
        sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=99)
        print(f"Subsampled to {adata.n_obs} cells")
    
    # Skip already-processed cells (use random_state=99 for new cells not overlapping with previous runs)
    if SKIP_CELLS > 0:
        adata = adata[SKIP_CELLS:]
        print(f"Skipping first {SKIP_CELLS} cells (already processed), {adata.n_obs} remaining")
    
    # Get gene symbols from data
    if 'feature_name' in adata.var.columns:
        data_gene_symbols = adata.var['feature_name'].values
    else:
        data_gene_symbols = adata.var.index.values
    
    # Build symbol -> ensembl for data genes
    data_symbol_to_ensembl = {}
    for sym in data_gene_symbols:
        ens = name_dict.get(sym)
        if ens:
            data_symbol_to_ensembl[sym] = ens
    print(f"Mapped {len(data_symbol_to_ensembl)}/{len(data_gene_symbols)} data genes to Ensembl IDs")
    
    # Process cells
    all_results = []
    
    for cell_idx in tqdm(range(adata.n_obs), desc="Cells"):
        cell = adata[cell_idx]
        expr_values = cell.X.toarray().flatten()
        
        # Build ensembl expression dict for tokenization
        ensembl_expr = {}
        for sym, val in zip(data_gene_symbols, expr_values):
            ens = data_symbol_to_ensembl.get(sym)
            if ens and val > 0:
                ensembl_expr[ens] = float(val)
        
        # Tokenize
        token_ids = tokenize_cell(None, ensembl_expr, token_dict, median_dict, MAX_TOKENS)
        if len(token_ids) < 10:
            continue
        
        # Get cell type
        cell_type = str(cell.obs['cell_type'].values[0]) if 'cell_type' in cell.obs.columns else 'Unknown'
        
        # Perturb each target gene
        for gene_sym, gene_info in gene_mapping.items():
            result = perturb_gene(model, token_ids, gene_info['token_id'], device, N_PERMUTATIONS)
            if result is not None:
                all_results.append({
                    'cell_idx': cell_idx,
                    'cell_type': cell_type,
                    'gene': gene_sym,
                    'cosine_distance': result['distance'],
                    'p_value': result['p_value'],
                    'null_mean': result.get('null_mean', np.nan),
                    'gene_in_cell': result['distance'] is not None,
                    'n_tokens': len(token_ids)
                })
        
        # Clear cache periodically
        if (cell_idx + 1) % BATCH_SIZE == 0:
            torch.cuda.empty_cache()
            gc.collect()
            # Save intermediate
            pd.DataFrame(all_results).to_csv(os.path.join(OUTPUT_DIR, 'intermediate_results_batch3.csv'), index=False)
            print(f"  Saved {len(all_results)} results so far")
    
    if not all_results:
        print("ERROR: No results collected!")
        return
    
    # Save raw results
    df = pd.DataFrame(all_results)
    df.to_csv(os.path.join(OUTPUT_DIR, 'cell_level_results_batch3.csv'), index=False)
    print(f"\nTotal results: {len(df)}")
    
    # Aggregate by gene
    gene_stats = []
    for gene in gene_mapping.keys():
        gd = df[df['gene'] == gene]
        if len(gd) == 0:
            continue
        dists = gd['cosine_distance'].dropna()
        nulls = gd['null_mean'].dropna()
        
        # Effect size: mean perturbation distance vs mean null distance
        mean_dist = dists.mean()
        mean_null = nulls.mean() if len(nulls) > 0 else 0
        effect = mean_dist - mean_null
        
        # Combined p-value (Fisher's method or simple median)
        pvals = gd['p_value'].dropna()
        median_p = pvals.median() if len(pvals) > 0 else np.nan
        
        gene_stats.append({
            'gene': gene,
            'n_cells': len(dists),
            'mean_cosine_distance': mean_dist,
            'std_cosine_distance': dists.std(),
            'mean_null_distance': mean_null,
            'effect_size': effect,
            'median_empirical_p': median_p,
            'min_p': pvals.min() if len(pvals) > 0 else np.nan
        })
    
    gene_df = pd.DataFrame(gene_stats).sort_values('effect_size', ascending=False)
    
    # FDR correction
    valid_p = gene_df['median_empirical_p'].dropna()
    if len(valid_p) > 0:
        _, fdr, _, _ = multipletests(valid_p, method='fdr_bh')
        gene_df.loc[valid_p.index, 'fdr'] = fdr
    
    gene_df.to_csv(os.path.join(OUTPUT_DIR, 'gene_statistics_phase4.csv'), index=False)
    print("\n=== TOP GENES BY EFFECT SIZE ===")
    print(gene_df[['gene','n_cells','mean_cosine_distance','effect_size','median_empirical_p','fdr']].to_string(index=False))
    
    # Aggregate by gene x cell type
    celltype_stats = []
    for gene in gene_mapping.keys():
        for ct in df['cell_type'].unique():
            subset = df[(df['gene'] == gene) & (df['cell_type'] == ct)]
            if len(subset) == 0:
                continue
            dists = subset['cosine_distance'].dropna()
            if len(dists) == 0:
                continue
            celltype_stats.append({
                'gene': gene, 'cell_type': ct,
                'n_cells': len(dists),
                'mean_distance': dists.mean(),
                'std_distance': dists.std(),
                'median_p': subset['p_value'].dropna().median()
            })
    
    ct_df = pd.DataFrame(celltype_stats)
    ct_df.to_csv(os.path.join(OUTPUT_DIR, 'gene_celltype_stats_phase4.csv'), index=False)
    
    print(f"\nAll results saved to {OUTPUT_DIR}")
    print("DONE")

if __name__ == '__main__':
    main()
