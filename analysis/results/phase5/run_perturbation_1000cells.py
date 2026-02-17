"""
Phase 5: Scale Up to 1000 Cells - Virtual Perturbation Screening
Extends Phase 4 to 1000 total cells by running 800 additional cells and combining results
NO EMOJI in print statements - Windows encoding issue
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
OUTPUT_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase5"
PHASE4_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\phase4"
DATA_PATH = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
MODEL_PATH = os.path.join(MODELS_DIR, "Geneformer-V1-10M")

TARGET_GENES = [
    'CADM2','NRXN1','NLGN1','NEGR1','TCF4','MEF2C','APP','GRIN2A','GRIN2B',
    'MAPT','SCN1A','DISC1','FOXO3','CAMK2A','HOMER1','CREB1','COMT','SHANK3',
    'FMR1','SYN1','BDNF','SLC6A4'
]

TARGET_CELLS = 1000  # Total target (including existing 200)
NEW_CELLS = 800      # Additional cells to process
BATCH_SIZE = 25      # Process in small batches for memory management
N_PERMUTATIONS = 30  # Null distribution samples
MAX_TOKENS = 2048
SKIP_EXISTING = 200  # Skip first 200 cells (already processed in Phase 4)

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
    
    return {'distance': cos_dist, 'null_distances': null_dists, 'p_value': p_value, 
            'null_mean': np.mean(null_dists) if null_dists else np.nan}

def load_existing_results():
    """Load Phase 4 results to combine with new data"""
    phase4_path = os.path.join(PHASE4_DIR, 'cell_level_results_combined.csv')
    if os.path.exists(phase4_path):
        print(f"Loading existing Phase 4 results from: {phase4_path}")
        df = pd.read_csv(phase4_path)
        print(f"Loaded {len(df)} existing results from {df['cell_idx'].nunique()} cells")
        return df
    else:
        print("No existing Phase 4 results found - starting fresh")
        return pd.DataFrame()

def save_intermediate_results(results, batch_num):
    """Save intermediate results to avoid data loss"""
    if results:
        df = pd.DataFrame(results)
        filepath = os.path.join(OUTPUT_DIR, f'intermediate_results_batch{batch_num}.csv')
        df.to_csv(filepath, index=False)
        print(f"Saved batch {batch_num}: {len(results)} results")

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # Load dictionaries
    token_dict, name_dict, median_dict = load_dictionaries()
    gene_mapping = build_gene_mapping(name_dict, token_dict)
    
    # Load model
    from transformers import BertModel
    model = BertModel.from_pretrained(MODEL_PATH)
    model.to(device)
    model.eval()
    print(f"Model loaded: {MODEL_PATH}")
    
    # Load existing Phase 4 results
    existing_results = load_existing_results()
    
    # Load data
    adata = sc.read_h5ad(DATA_PATH)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Skip already-processed cells and limit to target
    if SKIP_EXISTING > 0:
        adata_new = adata[SKIP_EXISTING:SKIP_EXISTING + NEW_CELLS]
        print(f"Processing cells {SKIP_EXISTING} to {SKIP_EXISTING + NEW_CELLS - 1}")
    else:
        if adata.n_obs > TARGET_CELLS:
            sc.pp.subsample(adata, n_obs=TARGET_CELLS, random_state=42)
        adata_new = adata
    
    print(f"Will process {adata_new.n_obs} new cells")
    
    # Get gene symbols from data
    if 'feature_name' in adata_new.var.columns:
        data_gene_symbols = adata_new.var['feature_name'].values
    else:
        data_gene_symbols = adata_new.var.index.values
    
    # Build symbol -> ensembl for data genes
    data_symbol_to_ensembl = {}
    for sym in data_gene_symbols:
        ens = name_dict.get(sym)
        if ens:
            data_symbol_to_ensembl[sym] = ens
    print(f"Mapped {len(data_symbol_to_ensembl)}/{len(data_gene_symbols)} data genes to Ensembl IDs")
    
    # Process cells in batches
    all_results = []
    batch_num = 1
    
    for cell_idx in tqdm(range(adata_new.n_obs), desc="Processing cells"):
        cell = adata_new[cell_idx]
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
        
        # Adjust cell_idx to account for skipped cells
        actual_cell_idx = cell_idx + SKIP_EXISTING
        
        # Perturb each target gene
        for gene_sym, gene_info in gene_mapping.items():
            result = perturb_gene(model, token_ids, gene_info['token_id'], device, N_PERMUTATIONS)
            if result is not None:
                all_results.append({
                    'cell_idx': actual_cell_idx,
                    'cell_type': cell_type,
                    'gene': gene_sym,
                    'cosine_distance': result['distance'],
                    'p_value': result['p_value'],
                    'null_mean': result.get('null_mean', np.nan),
                    'gene_in_cell': result['distance'] is not None,
                    'n_tokens': len(token_ids)
                })
        
        # Save intermediate results and clear cache periodically
        if (cell_idx + 1) % BATCH_SIZE == 0:
            save_intermediate_results(all_results[-BATCH_SIZE*len(gene_mapping):], batch_num)
            batch_num += 1
            torch.cuda.empty_cache()
            gc.collect()
            print(f"  Processed {cell_idx + 1}/{adata_new.n_obs} cells, {len(all_results)} total results")
    
    if not all_results:
        print("ERROR: No new results collected!")
        return
    
    # Save final batch if needed
    if len(all_results) % (BATCH_SIZE * len(gene_mapping)) != 0:
        save_intermediate_results(all_results[-(len(all_results) % (BATCH_SIZE * len(gene_mapping))):], batch_num)
    
    # Save new results
    new_df = pd.DataFrame(all_results)
    new_df.to_csv(os.path.join(OUTPUT_DIR, 'new_cell_results_800.csv'), index=False)
    print(f"\nNew results: {len(new_df)} from {new_df['cell_idx'].nunique()} cells")
    
    # Combine with existing results
    if len(existing_results) > 0:
        combined_df = pd.concat([existing_results, new_df], ignore_index=True)
        print(f"Combined with Phase 4: {len(combined_df)} total results from {combined_df['cell_idx'].nunique()} cells")
    else:
        combined_df = new_df
        print(f"Fresh run: {len(combined_df)} total results")
    
    combined_df.to_csv(os.path.join(OUTPUT_DIR, 'cell_level_results_1000.csv'), index=False)
    
    # Aggregate by gene
    gene_stats = []
    for gene in gene_mapping.keys():
        gd = combined_df[combined_df['gene'] == gene]
        if len(gd) == 0:
            continue
        dists = gd['cosine_distance'].dropna()
        nulls = gd['null_mean'].dropna()
        
        # Effect size: mean perturbation distance vs mean null distance
        mean_dist = dists.mean()
        mean_null = nulls.mean() if len(nulls) > 0 else 0
        effect = mean_dist - mean_null
        
        # Statistics
        pvals = gd['p_value'].dropna()
        median_p = pvals.median() if len(pvals) > 0 else np.nan
        n_sig_005 = len(pvals[pvals < 0.05]) if len(pvals) > 0 else 0
        pct_sig = (n_sig_005 / len(pvals) * 100) if len(pvals) > 0 else 0
        
        gene_stats.append({
            'gene': gene,
            'n_cells': len(dists),
            'mean_cosine_distance': mean_dist,
            'std_cosine_distance': dists.std(),
            'mean_null_distance': mean_null,
            'effect_size': effect,
            'median_empirical_p': median_p,
            'min_p': pvals.min() if len(pvals) > 0 else np.nan,
            'n_sig_005': n_sig_005,
            'pct_sig': pct_sig
        })
    
    gene_df = pd.DataFrame(gene_stats).sort_values('effect_size', ascending=False)
    
    # FDR correction
    valid_p = gene_df['median_empirical_p'].dropna()
    if len(valid_p) > 0:
        _, fdr, _, _ = multipletests(valid_p, method='fdr_bh')
        gene_df.loc[valid_p.index, 'fdr'] = fdr
    
    gene_df.to_csv(os.path.join(OUTPUT_DIR, 'gene_statistics_1000cells.csv'), index=False)
    
    print("\n=== TOP GENES BY EFFECT SIZE (1000 cells) ===")
    top_cols = ['gene','n_cells','mean_cosine_distance','effect_size','median_empirical_p','fdr','pct_sig']
    print(gene_df[top_cols].head(15).to_string(index=False))
    
    # Cell type analysis
    celltype_stats = []
    for gene in gene_mapping.keys():
        for ct in combined_df['cell_type'].unique():
            subset = combined_df[(combined_df['gene'] == gene) & (combined_df['cell_type'] == ct)]
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
    ct_df.to_csv(os.path.join(OUTPUT_DIR, 'gene_celltype_stats_1000cells.csv'), index=False)
    
    print(f"\nPhase 5 results saved to {OUTPUT_DIR}")
    print(f"Total cells processed: {combined_df['cell_idx'].nunique()}")
    print(f"Final FDR for top genes: {gene_df['fdr'].iloc[0]:.3f}")
    print("TASK 1 COMPLETE: Scaled to 1000 cells")

if __name__ == '__main__':
    main()