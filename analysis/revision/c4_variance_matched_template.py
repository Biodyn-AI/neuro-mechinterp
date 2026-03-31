"""
Template: Full variance-matched control gene selection for permutation null.

This script requires:
  1. The DLPFC_11k.h5ad brain scRNA-seq dataset
  2. The Geneformer token dictionary
  3. GPU access for running Geneformer in-silico perturbations

Run on a machine with GPU after installing geneformer.
"""
import os, sys, json, pickle, gc, random
import numpy as np
import scipy.sparse as sp
from collections import defaultdict

# ===== CONFIGURATION =====
BRAIN_H5AD = "/path/to/DLPFC_11k.h5ad"  # UPDATE THIS
TOKEN_DICT_PATH = "/path/to/token_dictionary_gc30M.pkl"  # UPDATE THIS
OUTPUT_DIR = "/path/to/output"  # UPDATE THIS
N_CONTROL_GENES = 50
EXPRESSION_TOLERANCE = 0.20
VARIANCE_TOLERANCE = 0.30  # Allow +/-30% tolerance on variance
DROPOUT_TOLERANCE = 0.15   # Allow +/-15% tolerance on dropout rate
RANDOM_SEED = 42
N_BINS = 5

INTELLIGENCE_GENES = {
    "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759",
    "GRIN2B": "ENSG00000273079", "CADM2": "ENSG00000175161",
    "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
    "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322",
    "HOMER1": "ENSG00000152413", "APP": "ENSG00000142192",
    "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
    "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868",
    "FOXO3": "ENSG00000118689", "CREB1": "ENSG00000118260",
    "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
    "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576",
    "COMT": "ENSG00000093010",
}


def compute_gene_statistics(h5ad_path):
    """
    Compute per-gene expression statistics from the h5ad file.

    Returns dict: ensembl_id -> {symbol, frequency, mean_expr, variance,
                                  dropout_rate, cv}
    """
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path)
    n_cells = adata.shape[0]

    var = adata.var
    if 'feature_name' in var.columns:
        gene_symbols = var['feature_name'].values
    else:
        gene_symbols = var.index.values
    ensembl_ids = var.index.values

    X = adata.X
    if sp.issparse(X):
        X_dense = None  # Process in sparse mode for memory
        n_expressing = np.array((X > 0).sum(axis=0)).flatten()
        gene_means = np.array(X.mean(axis=0)).flatten()
        # Variance: E[X^2] - E[X]^2
        X_sq = X.copy()
        X_sq.data **= 2
        gene_var = np.array(X_sq.mean(axis=0)).flatten() - gene_means ** 2
    else:
        n_expressing = np.array((X > 0).sum(axis=0)).flatten()
        gene_means = np.mean(X, axis=0)
        gene_var = np.var(X, axis=0)

    expr_freq = n_expressing / n_cells
    dropout_rate = 1.0 - expr_freq
    cv = np.where(gene_means > 0, np.sqrt(gene_var) / gene_means, 0.0)

    gene_stats = {}
    for i, (ens, sym) in enumerate(zip(ensembl_ids, gene_symbols)):
        gene_stats[ens] = {
            'symbol': sym,
            'frequency': float(expr_freq[i]),
            'mean_expr': float(gene_means[i]),
            'variance': float(gene_var[i]),
            'dropout_rate': float(dropout_rate[i]),
            'cv': float(cv[i]),
        }

    del adata
    gc.collect()
    return gene_stats


def select_variance_matched_controls(gene_stats, token_dict_path):
    """
    Select 50 control genes matched on expression frequency, variance,
    and dropout rate.
    """
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    with open(token_dict_path, 'rb') as f:
        token_dict = pickle.load(f)

    intel_ensembls = set(INTELLIGENCE_GENES.values())
    exclude = intel_ensembls  # Add other exclusions as needed

    # For each intelligence gene, find candidates matching on multiple axes
    all_candidates = {}
    for symbol, ensembl in INTELLIGENCE_GENES.items():
        if ensembl not in gene_stats:
            continue
        target = gene_stats[ensembl]

        for ens, info in gene_stats.items():
            if ens in exclude or ens not in token_dict:
                continue

            # Check expression frequency match
            freq_ratio = info['frequency'] / target['frequency'] if target['frequency'] > 0 else float('inf')
            if not (1 - EXPRESSION_TOLERANCE <= freq_ratio <= 1 + EXPRESSION_TOLERANCE):
                continue

            # Check variance match
            if target['variance'] > 0:
                var_ratio = info['variance'] / target['variance']
                if not (1 - VARIANCE_TOLERANCE <= var_ratio <= 1 + VARIANCE_TOLERANCE):
                    continue

            # Check dropout rate match
            dropout_diff = abs(info['dropout_rate'] - target['dropout_rate'])
            if dropout_diff > DROPOUT_TOLERANCE:
                continue

            if ens not in all_candidates:
                all_candidates[ens] = {
                    'symbol': info['symbol'],
                    'frequency': info['frequency'],
                    'variance': info['variance'],
                    'dropout_rate': info['dropout_rate'],
                    'cv': info['cv'],
                    'matched_to': [],
                }
            all_candidates[ens]['matched_to'].append(symbol)

    print(f"Candidates passing mean+variance+dropout filter: {len(all_candidates)}")

    # Proceed with bin-stratified sampling (same as original)
    # ... (Use the same binning logic from permutation_null_v2_select.py)

    return all_candidates


if __name__ == "__main__":
    print("Step 1: Computing per-gene statistics from h5ad...")
    gene_stats = compute_gene_statistics(BRAIN_H5AD)

    print("Step 2: Selecting variance-matched controls...")
    candidates = select_variance_matched_controls(gene_stats, TOKEN_DICT_PATH)

    print(f"Step 3: Save selection and run Geneformer perturbations (GPU required)")
    # ... Save and run perturbations per gene ...
    print("TEMPLATE COMPLETE - implement perturbation runner as needed")
