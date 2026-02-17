"""
Intelligence Augmentation - Scaled Perturbation Analysis for 2000+ Cells
Target: Process 2000+ cells from Siletti_DLPFC_113k.h5ad dataset
Strategy: Batch processing (150 cells/batch) with explicit memory management
Author: Claude Code (OpenClaw Subagent)
Date: 2026-02-14
"""

import os
import sys
import gc
import pickle
import json
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.spatial.distance import cosine
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import warnings
import time
from pathlib import Path

warnings.filterwarnings('ignore')

# Paths and configuration
MODELS_DIR = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
OUTPUT_DIR = r"D:\openclaw\intelligence-augmentation\analysis\results\scaled"
DATA_PATH = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\Siletti_DLPFC_113k.h5ad"
MODEL_PATH = os.path.join(MODELS_DIR, "Geneformer-V1-10M")

# Target genes from GWAS intelligence studies
TARGET_GENES = [
    'CADM2','NRXN1','NLGN1','NEGR1','TCF4','MEF2C','APP','GRIN2A','GRIN2B',
    'MAPT','SCN1A','DISC1','FOXO3','CAMK2A','HOMER1','CREB1','COMT','SHANK3',
    'FMR1','SYN1','BDNF','SLC6A4'
]

# Processing parameters
TARGET_CELLS = 2500       # Target total cells (aim for 2000+ with some buffer)
BATCH_SIZE = 150          # Cells per batch (balance between memory and efficiency)
N_PERMUTATIONS = 30       # Null distribution samples
MAX_TOKENS = 2048         # Maximum tokens per cell
RANDOM_SEED = 42          # For reproducibility

def setup_output_dir():
    """Create output directory structure"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, 'batches'), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, 'logs'), exist_ok=True)
    return True

def log_message(message, log_file="progress.log"):
    """Log progress messages with timestamp"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_path = os.path.join(OUTPUT_DIR, 'logs', log_file)
    
    # Print to console
    print(f"[{timestamp}] {message}")
    
    # Write to log file
    with open(log_path, 'a', encoding='utf-8') as f:
        f.write(f"[{timestamp}] {message}\n")

def load_dictionaries():
    """Load Geneformer token and gene dictionaries"""
    try:
        dict_dir = os.path.join(MODELS_DIR, "geneformer")
        
        with open(os.path.join(dict_dir, "token_dictionary_gc104M.pkl"), 'rb') as f:
            token_dict = pickle.load(f)
        
        with open(os.path.join(dict_dir, "gene_name_id_dict_gc104M.pkl"), 'rb') as f:
            name_dict = pickle.load(f)
        
        with open(os.path.join(dict_dir, "gene_median_dictionary_gc104M.pkl"), 'rb') as f:
            median_dict = pickle.load(f)
        
        log_message(f"Loaded dictionaries: {len(token_dict)} tokens, {len(name_dict)} genes")
        return token_dict, name_dict, median_dict
    
    except Exception as e:
        log_message(f"ERROR loading dictionaries: {e}")
        raise

def build_gene_mapping(name_dict, token_dict):
    """Build gene symbol to token mapping"""
    mapping = {}
    missing = []
    
    for gene in TARGET_GENES:
        ensembl = name_dict.get(gene)
        if ensembl and ensembl in token_dict:
            mapping[gene] = {
                'ensembl': ensembl, 
                'token_id': token_dict[ensembl]
            }
        else:
            missing.append(gene)
    
    log_message(f"Gene mapping: {len(mapping)}/{len(TARGET_GENES)} genes mapped")
    if missing:
        log_message(f"Missing genes: {missing}")
    
    return mapping

def tokenize_cell(ensembl_expr, token_dict, median_dict, max_tokens=2048):
    """
    Convert single cell expression to Geneformer tokens
    
    Args:
        ensembl_expr: {ensembl_id: expression_value}
        token_dict: Ensembl ID to token mapping
        median_dict: Gene median expression for normalization
        max_tokens: Maximum tokens per cell
    
    Returns:
        List of token IDs
    """
    # Filter to genes with expression > 0 and in token dictionary
    valid_genes = {g: e for g, e in ensembl_expr.items() 
                   if g in token_dict and e > 0}
    
    if not valid_genes:
        return []
    
    # Rank genes by expression/median ratio (Geneformer approach)
    ranked = []
    for gene, expr in valid_genes.items():
        median = median_dict.get(gene, 1.0)
        normalized = expr / median if median > 0 else expr
        ranked.append((gene, normalized))
    
    # Sort by normalized expression (descending)
    ranked.sort(key=lambda x: x[1], reverse=True)
    
    # Convert to token IDs
    token_ids = [token_dict['<cls>']]  # Start with CLS token
    
    for gene, _ in ranked[:max_tokens - 2]:  # Reserve space for CLS and EOS
        token_ids.append(token_dict[gene])
    
    token_ids.append(token_dict['<eos>'])  # End with EOS token
    
    return token_ids

def get_embedding(model, token_ids, device):
    """Get CLS embedding from Geneformer model"""
    if len(token_ids) < 3:  # Need at least CLS, one gene, EOS
        return None
    
    try:
        input_ids = torch.tensor(token_ids, dtype=torch.long).unsqueeze(0).to(device)
        
        with torch.no_grad():
            outputs = model(input_ids)
            cls_embedding = outputs.last_hidden_state[0, 0, :].cpu().numpy()
        
        return cls_embedding
    
    except Exception as e:
        log_message(f"ERROR in get_embedding: {e}")
        return None

def perturb_gene_analysis(model, token_ids, target_token_id, device, n_perms=30):
    """
    Perform gene perturbation analysis
    
    Args:
        model: Geneformer model
        token_ids: Original token sequence
        target_token_id: Token to remove (perturb)
        device: torch device
        n_perms: Number of null permutations
    
    Returns:
        Dict with perturbation results
    """
    if target_token_id not in token_ids:
        return None
    
    # Get baseline embedding (unperturbed)
    baseline_emb = get_embedding(model, token_ids, device)
    if baseline_emb is None:
        return None
    
    # Get perturbed embedding (target gene removed)
    perturbed_tokens = [t for t in token_ids if t != target_token_id]
    perturbed_emb = get_embedding(model, perturbed_tokens, device)
    if perturbed_emb is None:
        return None
    
    # Calculate cosine distance
    cos_distance = cosine(baseline_emb, perturbed_emb)
    
    # Build null distribution: remove random non-target genes
    special_tokens = {token_ids[0], token_ids[-1], target_token_id}  # CLS, EOS, target
    available_for_null = [t for t in token_ids if t not in special_tokens]
    
    null_distances = []
    if available_for_null:
        n_sample = min(n_perms, len(available_for_null))
        random_tokens = np.random.choice(available_for_null, n_sample, replace=False)
        
        for random_token in random_tokens:
            null_tokens = [t for t in token_ids if t != random_token]
            null_emb = get_embedding(model, null_tokens, device)
            
            if null_emb is not None:
                null_dist = cosine(baseline_emb, null_emb)
                null_distances.append(null_dist)
    
    # Calculate empirical p-value
    if null_distances:
        p_value = (np.sum(np.array(null_distances) >= cos_distance) + 1) / (len(null_distances) + 1)
        null_mean = np.mean(null_distances)
        null_std = np.std(null_distances)
    else:
        p_value = np.nan
        null_mean = np.nan
        null_std = np.nan
    
    return {
        'cosine_distance': cos_distance,
        'p_value': p_value,
        'null_mean': null_mean,
        'null_std': null_std,
        'n_null_samples': len(null_distances)
    }

def process_cell_batch(adata_batch, model, gene_mapping, token_dict, median_dict, 
                      device, batch_num, data_symbol_to_ensembl):
    """Process a batch of cells"""
    
    log_message(f"Processing batch {batch_num}: {adata_batch.n_obs} cells")
    
    batch_results = []
    
    for cell_idx in range(adata_batch.n_obs):
        cell = adata_batch[cell_idx]
        
        # Get expression values
        if hasattr(cell.X, 'toarray'):
            expr_values = cell.X.toarray().flatten()
        else:
            expr_values = cell.X.flatten()
        
        # Get gene symbols
        if 'feature_name' in adata_batch.var.columns:
            gene_symbols = adata_batch.var['feature_name'].values
        elif 'gene_symbols' in adata_batch.var.columns:
            gene_symbols = adata_batch.var['gene_symbols'].values
        else:
            gene_symbols = adata_batch.var.index.values
        
        # Build Ensembl expression dict
        ensembl_expr = {}
        for sym, val in zip(gene_symbols, expr_values):
            ensembl = data_symbol_to_ensembl.get(sym)
            if ensembl and val > 0:
                ensembl_expr[ensembl] = float(val)
        
        # Tokenize cell
        token_ids = tokenize_cell(ensembl_expr, token_dict, median_dict, MAX_TOKENS)
        
        if len(token_ids) < 10:  # Skip cells with too few genes
            continue
        
        # Get cell metadata
        cell_type = str(cell.obs['cell_type'].values[0]) if 'cell_type' in cell.obs.columns else 'Unknown'
        actual_cell_idx = cell_idx + (batch_num - 1) * BATCH_SIZE
        
        # Test each target gene
        for gene_symbol, gene_info in gene_mapping.items():
            result = perturb_gene_analysis(
                model, token_ids, gene_info['token_id'], device, N_PERMUTATIONS
            )
            
            if result is not None:
                batch_results.append({
                    'batch_num': batch_num,
                    'cell_idx': actual_cell_idx,
                    'cell_type': cell_type,
                    'gene': gene_symbol,
                    'cosine_distance': result['cosine_distance'],
                    'p_value': result['p_value'],
                    'null_mean': result['null_mean'],
                    'null_std': result['null_std'],
                    'n_null_samples': result['n_null_samples'],
                    'n_tokens': len(token_ids)
                })
    
    return batch_results

def save_batch_results(results, batch_num):
    """Save batch results to intermediate file"""
    if results:
        df = pd.DataFrame(results)
        filepath = os.path.join(OUTPUT_DIR, 'batches', f'batch_{batch_num:03d}_results.csv')
        df.to_csv(filepath, index=False)
        log_message(f"Saved batch {batch_num}: {len(results)} perturbation results")
        return filepath
    return None

def memory_cleanup():
    """Explicit memory cleanup"""
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()

def combine_all_batches():
    """Combine all batch results into final dataset"""
    log_message("Combining all batch results...")
    
    batch_dir = os.path.join(OUTPUT_DIR, 'batches')
    batch_files = sorted([f for f in os.listdir(batch_dir) if f.startswith('batch_') and f.endswith('.csv')])
    
    all_results = []
    for batch_file in batch_files:
        filepath = os.path.join(batch_dir, batch_file)
        try:
            df = pd.read_csv(filepath)
            all_results.append(df)
            log_message(f"Loaded {batch_file}: {len(df)} results")
        except Exception as e:
            log_message(f"ERROR loading {batch_file}: {e}")
    
    if not all_results:
        log_message("ERROR: No batch results found!")
        return None
    
    # Combine all batches
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Remove duplicates if any
    initial_len = len(combined_df)
    combined_df = combined_df.drop_duplicates(subset=['cell_idx', 'gene'])
    final_len = len(combined_df)
    
    if initial_len != final_len:
        log_message(f"Removed {initial_len - final_len} duplicate measurements")
    
    # Save combined results
    combined_path = os.path.join(OUTPUT_DIR, 'combined_perturbation_results.csv')
    combined_df.to_csv(combined_path, index=False)
    
    log_message(f"Combined results: {len(combined_df)} measurements from {combined_df['cell_idx'].nunique()} cells")
    
    return combined_df

def calculate_gene_statistics(combined_df):
    """Calculate gene-level statistics with FDR correction"""
    log_message("Calculating gene-level statistics...")
    
    gene_stats = []
    
    for gene in TARGET_GENES:
        gene_data = combined_df[combined_df['gene'] == gene].copy()
        
        if len(gene_data) == 0:
            continue
        
        # Basic statistics
        distances = gene_data['cosine_distance'].dropna()
        null_means = gene_data['null_mean'].dropna()
        p_values = gene_data['p_value'].dropna()
        
        if len(distances) == 0:
            continue
        
        # Effect size (perturbation vs null)
        mean_distance = distances.mean()
        mean_null = null_means.mean() if len(null_means) > 0 else 0
        effect_size = mean_distance - mean_null
        
        # Statistical measures
        median_p = p_values.median() if len(p_values) > 0 else np.nan
        min_p = p_values.min() if len(p_values) > 0 else np.nan
        n_sig_005 = len(p_values[p_values < 0.05]) if len(p_values) > 0 else 0
        pct_sig = (n_sig_005 / len(p_values) * 100) if len(p_values) > 0 else 0
        
        gene_stats.append({
            'gene': gene,
            'n_cells': len(distances),
            'mean_cosine_distance': mean_distance,
            'std_cosine_distance': distances.std(),
            'mean_null_distance': mean_null,
            'effect_size': effect_size,
            'median_empirical_p': median_p,
            'min_p': min_p,
            'n_significant_005': n_sig_005,
            'percent_significant': pct_sig,
            'mean_null_samples': gene_data['n_null_samples'].mean()
        })
    
    gene_df = pd.DataFrame(gene_stats)
    gene_df = gene_df.sort_values('effect_size', ascending=False)
    
    # FDR correction on median p-values
    valid_p = gene_df['median_empirical_p'].dropna()
    if len(valid_p) > 0:
        _, fdr_corrected, _, _ = multipletests(valid_p, method='fdr_bh')
        gene_df.loc[valid_p.index, 'fdr_corrected'] = fdr_corrected
        
        # Calculate overall FDR
        overall_fdr = fdr_corrected.min() if len(fdr_corrected) > 0 else 1.0
        log_message(f"FDR correction applied. Minimum FDR: {overall_fdr:.4f}")
    else:
        gene_df['fdr_corrected'] = np.nan
        overall_fdr = 1.0
        log_message("No valid p-values for FDR correction")
    
    # Save gene statistics
    gene_stats_path = os.path.join(OUTPUT_DIR, 'gene_statistics_scaled.csv')
    gene_df.to_csv(gene_stats_path, index=False)
    
    return gene_df, overall_fdr

def generate_summary_report(combined_df, gene_df, overall_fdr):
    """Generate comprehensive analysis summary"""
    
    summary_path = os.path.join(OUTPUT_DIR, 'SCALED_ANALYSIS_SUMMARY.md')
    
    # Summary statistics
    n_cells = combined_df['cell_idx'].nunique()
    n_measurements = len(combined_df)
    n_cell_types = combined_df['cell_type'].nunique()
    n_genes = len(gene_df)
    
    # Top genes
    top_genes = gene_df.head(10)
    
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write("# Intelligence Augmentation - Scaled Perturbation Analysis\n\n")
        f.write(f"**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**Analysis:** Scaled perturbation screening for cognitive enhancement\n")
        f.write(f"**Model:** Geneformer V1-10M\n\n")
        
        f.write("## Executive Summary\n\n")
        f.write(f"- **Total cells analyzed:** {n_cells:,}\n")
        f.write(f"- **Total measurements:** {n_measurements:,}\n")
        f.write(f"- **Cell types represented:** {n_cell_types}\n")
        f.write(f"- **Genes analyzed:** {n_genes}/22 target genes\n")
        f.write(f"- **Overall FDR:** {overall_fdr:.4f}\n\n")
        
        f.write("## Top Genes by Effect Size\n\n")
        f.write("| Rank | Gene | Effect Size | n_cells | FDR p-value | % Significant |\n")
        f.write("|------|------|-------------|---------|-------------|---------------|\n")
        
        for i, row in top_genes.iterrows():
            rank = top_genes.index.get_loc(i) + 1
            fdr = row.get('fdr_corrected', np.nan)
            fdr_str = f"{fdr:.3f}" if not pd.isna(fdr) else "N/A"
            f.write(f"| {rank} | **{row['gene']}** | {row['effect_size']:.4f} | {row['n_cells']} | {fdr_str} | {row['percent_significant']:.1f}% |\n")
        
        f.write(f"\n## Cell Type Distribution\n\n")
        cell_type_counts = combined_df['cell_type'].value_counts().head(10)
        for ct, count in cell_type_counts.items():
            pct = (count / n_measurements) * 100
            f.write(f"- **{ct}:** {count:,} measurements ({pct:.1f}%)\n")
        
        f.write(f"\n## Analysis Complete\n\n")
        f.write(f"Successfully scaled perturbation analysis to {n_cells:,} cells.\n")
        f.write(f"Results saved to: {OUTPUT_DIR}\n")
    
    log_message(f"Summary report saved: {summary_path}")

def main():
    """Main execution function"""
    
    # Setup
    setup_output_dir()
    log_message("=== Intelligence Augmentation - Scaled Perturbation Analysis ===")
    log_message(f"Target: {TARGET_CELLS} cells, Batch size: {BATCH_SIZE}")
    
    # Check device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    log_message(f"Using device: {device}")
    
    # Set random seed
    np.random.seed(RANDOM_SEED)
    torch.manual_seed(RANDOM_SEED)
    
    try:
        # Load dictionaries and build gene mapping
        log_message("Loading Geneformer dictionaries...")
        token_dict, name_dict, median_dict = load_dictionaries()
        gene_mapping = build_gene_mapping(name_dict, token_dict)
        
        if not gene_mapping:
            log_message("ERROR: No gene mappings found!")
            return
        
        # Load model
        log_message("Loading Geneformer model...")
        from transformers import BertModel
        model = BertModel.from_pretrained(MODEL_PATH)
        model.to(device)
        model.eval()
        log_message(f"Model loaded and moved to {device}")
        
        # Load dataset
        log_message(f"Loading dataset: {DATA_PATH}")
        adata = sc.read_h5ad(DATA_PATH)
        log_message(f"Dataset loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        
        # Limit to target number of cells
        if adata.n_obs > TARGET_CELLS:
            log_message(f"Subsampling to {TARGET_CELLS} cells...")
            sc.pp.subsample(adata, n_obs=TARGET_CELLS, random_state=RANDOM_SEED)
        
        log_message(f"Processing {adata.n_obs:,} cells")
        
        # Build gene symbol mapping
        if 'feature_name' in adata.var.columns:
            data_gene_symbols = adata.var['feature_name'].values
        elif 'gene_symbols' in adata.var.columns:
            data_gene_symbols = adata.var['gene_symbols'].values
        else:
            data_gene_symbols = adata.var.index.values
        
        data_symbol_to_ensembl = {}
        for sym in data_gene_symbols:
            ens = name_dict.get(sym)
            if ens:
                data_symbol_to_ensembl[sym] = ens
        
        log_message(f"Mapped {len(data_symbol_to_ensembl):,}/{len(data_gene_symbols):,} data genes to Ensembl IDs")
        
        # Process in batches
        n_batches = (adata.n_obs + BATCH_SIZE - 1) // BATCH_SIZE
        log_message(f"Processing {n_batches} batches of {BATCH_SIZE} cells each...")
        
        for batch_num in range(1, n_batches + 1):
            start_idx = (batch_num - 1) * BATCH_SIZE
            end_idx = min(start_idx + BATCH_SIZE, adata.n_obs)
            
            log_message(f"Batch {batch_num}/{n_batches}: cells {start_idx}-{end_idx-1}")
            
            # Extract batch
            adata_batch = adata[start_idx:end_idx].copy()
            
            # Process batch
            batch_results = process_cell_batch(
                adata_batch, model, gene_mapping, token_dict, median_dict,
                device, batch_num, data_symbol_to_ensembl
            )
            
            # Save batch results
            if batch_results:
                save_batch_results(batch_results, batch_num)
            else:
                log_message(f"WARNING: No results for batch {batch_num}")
            
            # Memory cleanup
            del adata_batch
            memory_cleanup()
            
            log_message(f"Batch {batch_num} complete: {len(batch_results) if batch_results else 0} results")
        
        # Combine all results
        log_message("Combining all batch results...")
        combined_df = combine_all_batches()
        
        if combined_df is None:
            log_message("ERROR: Failed to combine batch results")
            return
        
        # Calculate gene statistics
        gene_df, overall_fdr = calculate_gene_statistics(combined_df)
        
        # Generate summary report
        generate_summary_report(combined_df, gene_df, overall_fdr)
        
        # Final summary
        n_cells_final = combined_df['cell_idx'].nunique()
        n_measurements_final = len(combined_df)
        
        log_message("=== ANALYSIS COMPLETE ===")
        log_message(f"Final results:")
        log_message(f"  - Cells processed: {n_cells_final:,}")
        log_message(f"  - Total measurements: {n_measurements_final:,}")
        log_message(f"  - Overall FDR: {overall_fdr:.4f}")
        log_message(f"  - Results saved to: {OUTPUT_DIR}")
        
        if overall_fdr < 0.65:
            log_message(f"SUCCESS: Improved FDR from ~0.65 to {overall_fdr:.4f}")
        else:
            log_message(f"FDR not significantly improved. May need more cells or different approach.")
        
        # Top gene summary
        top_gene = gene_df.iloc[0]
        log_message(f"Top gene: {top_gene['gene']} (effect size: {top_gene['effect_size']:.4f}, {top_gene['n_cells']} cells)")
        
    except Exception as e:
        log_message(f"CRITICAL ERROR: {e}")
        import traceback
        log_message(f"Traceback: {traceback.format_exc()}")
        raise

if __name__ == '__main__':
    main()