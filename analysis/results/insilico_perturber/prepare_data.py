#!/usr/bin/env python3
"""
Prepare brain scRNA-seq data for Geneformer InSilicoPerturber analysis.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import pickle
import os

def prepare_brain_data(input_path, output_path):
    """Prepare brain data for Geneformer tokenization."""
    print(f"Loading data from: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    print(f"Original data shape: {adata.shape}")
    
    # Add n_counts column (required by Geneformer tokenizer)
    if 'sum' in adata.obs.columns:
        adata.obs['n_counts'] = adata.obs['sum']
        print("[OK] Added n_counts column from sum")
    else:
        # Calculate total counts if not available
        if hasattr(adata.X, 'toarray'):
            adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        else:
            adata.obs['n_counts'] = adata.X.sum(axis=1)
        print("[OK] Calculated n_counts from expression matrix")
    
    # Map intelligence genes to Ensembl IDs
    intelligence_genes = [
        "MEF2C", "BDNF", "GRIN2B", "CADM2", "NRXN1", "CAMK2A", 
        "GRIN2A", "SHANK3", "HOMER1", "APP", "NEGR1", "NLGN1", 
        "TCF4", "MAPT", "FOXO3", "CREB1", "FMR1", "SYN1", 
        "SCN1A", "SLC6A4", "COMT"
    ]
    
    # Create mapping from gene symbols to Ensembl IDs
    gene_symbol_to_ensembl = {}
    if 'feature_name' in adata.var.columns:
        symbol_to_ensembl = dict(zip(adata.var['feature_name'], adata.var.index))
        
        for gene in intelligence_genes:
            if gene in symbol_to_ensembl:
                gene_symbol_to_ensembl[gene] = symbol_to_ensembl[gene]
                print(f"  {gene} -> {symbol_to_ensembl[gene]}")
    
    print(f"Found {len(gene_symbol_to_ensembl)} intelligence genes with Ensembl IDs")
    
    # Save the mapping
    mapping_path = os.path.join(os.path.dirname(output_path), 'intelligence_genes_mapping.pkl')
    with open(mapping_path, 'wb') as f:
        pickle.dump(gene_symbol_to_ensembl, f)
    print(f"Saved gene mapping to: {mapping_path}")
    
    # Add ensembl_id column (required by Geneformer)
    adata.var['ensembl_id'] = adata.var.index
    
    # Ensure raw counts (Geneformer requirement)
    if adata.X.dtype != 'int':
        print("Converting expression matrix to integer counts...")
        # Round and convert to int if needed
        adata.X = np.round(adata.X).astype('int32')
    
    # Basic filtering - keep cells with reasonable count numbers
    min_counts = 500
    max_counts = 50000
    min_genes = 200
    
    print(f"Filtering cells: {min_counts} <= counts <= {max_counts}, genes >= {min_genes}")
    
    # Filter cells
    cell_filter = (adata.obs['n_counts'] >= min_counts) & (adata.obs['n_counts'] <= max_counts)
    
    # Filter genes (expressed in at least 3 cells)
    gene_filter = (adata.X > 0).sum(axis=0) >= 3
    if hasattr(gene_filter, 'A1'):
        gene_filter = gene_filter.A1
    
    print(f"Before filtering: {adata.shape}")
    adata = adata[cell_filter, gene_filter]
    print(f"After filtering: {adata.shape}")
    
    # Save prepared data
    adata.write_h5ad(output_path)
    print(f"Saved prepared data to: {output_path}")
    
    return adata, gene_symbol_to_ensembl

def attempt_simple_perturbation(adata, gene_mapping):
    """
    Attempt a simplified perturbation analysis without full Geneformer.
    This will serve as a fallback if the main approach fails.
    """
    print("\n=== FALLBACK: Simple perturbation analysis ===")
    
    # Get intelligence gene Ensembl IDs that are present
    intelligence_ensembl_ids = list(gene_mapping.values())
    present_genes = [ensembl_id for ensembl_id in intelligence_ensembl_ids if ensembl_id in adata.var.index]
    
    print(f"Found {len(present_genes)} intelligence genes in the dataset")
    
    # Simple approach: Calculate expression statistics for each gene
    results = []
    for ensembl_id in present_genes:
        gene_symbol = [k for k, v in gene_mapping.items() if v == ensembl_id][0]
        gene_idx = adata.var.index.get_loc(ensembl_id)
        
        if hasattr(adata.X, 'toarray'):
            expression = adata.X[:, gene_idx].toarray().flatten()
        else:
            expression = adata.X[:, gene_idx]
        
        mean_expr = np.mean(expression)
        std_expr = np.std(expression)
        pct_expressed = (expression > 0).sum() / len(expression) * 100
        
        results.append({
            'gene_symbol': gene_symbol,
            'ensembl_id': ensembl_id,
            'mean_expression': mean_expr,
            'std_expression': std_expr,
            'percent_expressed': pct_expressed,
            'max_expression': np.max(expression)
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('mean_expression', ascending=False)
    
    print("\nTop intelligence genes by mean expression:")
    print(results_df[['gene_symbol', 'mean_expression', 'percent_expressed']].head(10))
    
    return results_df

if __name__ == "__main__":
    # Prepare the smaller dataset first
    input_path = "D:/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
    output_path = "D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/DLPFC_11k_prepared.h5ad"
    
    adata, gene_mapping = prepare_brain_data(input_path, output_path)
    
    # Save basic statistics
    stats_df = attempt_simple_perturbation(adata, gene_mapping)
    stats_path = "D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/intelligence_genes_stats.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"\nSaved gene statistics to: {stats_path}")