#!/usr/bin/env python3
"""
Examine brain scRNA-seq data structure and prepare for Geneformer analysis.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import sys

def examine_brain_data(data_path):
    """Examine the structure of brain scRNA-seq data."""
    print(f"Loading data from: {data_path}")
    
    # Load the data
    adata = sc.read_h5ad(data_path)
    
    print(f"Data shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    
    print("\n=== Cell metadata (obs) ===")
    print(f"Column names: {list(adata.obs.columns)}")
    print(adata.obs.head())
    
    print("\n=== Gene metadata (var) ===")
    print(f"Column names: {list(adata.var.columns)}")
    print(adata.var.head())
    
    # Check if we have Ensembl IDs
    if 'ensembl_id' in adata.var.columns:
        print("\n[OK] Found ensembl_id in gene metadata")
        print(f"Number of genes with Ensembl IDs: {adata.var['ensembl_id'].notna().sum()}")
    else:
        print("\n[WARNING] No 'ensembl_id' column found in gene metadata")
        print("Available gene identifier columns:")
        for col in adata.var.columns:
            if any(keyword in col.lower() for keyword in ['gene', 'symbol', 'ensembl', 'id']):
                print(f"  - {col}: {adata.var[col].head().tolist()}")
    
    # Check for cell counts
    if 'n_counts' in adata.obs.columns:
        print("\n[OK] Found n_counts in cell metadata")
    elif 'total_counts' in adata.obs.columns:
        print("\n[OK] Found total_counts in cell metadata (can be renamed to n_counts)")
    else:
        print("\n[WARNING] No 'n_counts' or 'total_counts' column found")
        print("Calculating total counts per cell...")
        total_counts = adata.X.sum(axis=1)
        if hasattr(total_counts, 'A1'):  # If sparse matrix
            total_counts = total_counts.A1
        print(f"Total counts range: {total_counts.min():.0f} - {total_counts.max():.0f}")
    
    # Check data matrix type
    print(f"\n=== Data matrix ===")
    print(f"Type: {type(adata.X)}")
    print(f"Data type: {adata.X.dtype}")
    
    # Check if data is raw counts or normalized
    data_sample = adata.X[:100, :100]
    if hasattr(data_sample, 'toarray'):
        data_sample = data_sample.toarray()
    
    is_integer = np.allclose(data_sample, np.round(data_sample))
    has_large_values = np.any(data_sample > 10)
    
    print(f"Contains integer values: {is_integer}")
    print(f"Contains large values (>10): {has_large_values}")
    print(f"Sample values: {data_sample.flatten()[:10]}")
    
    if is_integer and has_large_values:
        print("[OK] Likely raw count data (suitable for Geneformer)")
    else:
        print("[WARNING] May be normalized data (Geneformer needs raw counts)")
    
    return adata

def check_intelligence_genes(adata):
    """Check which intelligence genes are present in the dataset."""
    intelligence_genes = [
        "MEF2C", "BDNF", "GRIN2B", "CADM2", "NRXN1", "CAMK2A", 
        "GRIN2A", "SHANK3", "HOMER1", "APP", "NEGR1", "NLGN1", 
        "TCF4", "MAPT", "FOXO3", "CREB1", "FMR1", "SYN1", 
        "SCN1A", "SLC6A4", "COMT", "NRXN1"
    ]
    
    print(f"\n=== Intelligence genes check ===")
    
    # Try different gene identifier columns
    gene_cols = []
    if 'gene_symbols' in adata.var.columns:
        gene_cols.append('gene_symbols')
    if 'gene_symbol' in adata.var.columns:
        gene_cols.append('gene_symbol')
    if hasattr(adata.var, 'index'):
        gene_cols.append('index')
    
    found_genes = {}
    for col in gene_cols:
        if col == 'index':
            gene_names = adata.var.index.tolist()
        else:
            gene_names = adata.var[col].tolist()
        
        found_in_col = []
        for gene in intelligence_genes:
            if gene in gene_names:
                found_in_col.append(gene)
        
        if found_in_col:
            found_genes[col] = found_in_col
            print(f"In {col}: {found_in_col}")
    
    if not found_genes:
        print("[WARNING] No intelligence genes found with exact symbol matching")
        print("Trying case-insensitive search...")
        
        for col in gene_cols:
            if col == 'index':
                gene_names = [g.upper() for g in adata.var.index.tolist()]
            else:
                gene_names = [str(g).upper() for g in adata.var[col].tolist()]
            
            found_in_col = []
            for gene in intelligence_genes:
                if gene.upper() in gene_names:
                    found_in_col.append(gene)
            
            if found_in_col:
                print(f"In {col} (case-insensitive): {found_in_col}")

if __name__ == "__main__":
    # Examine the smaller dataset first
    print("=" * 60)
    print("EXAMINING DLPFC_11k.h5ad")
    print("=" * 60)
    
    adata_small = examine_brain_data("D:/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")
    check_intelligence_genes(adata_small)
    
    print("\n" + "=" * 60)
    print("EXAMINING Siletti_DLPFC_113k.h5ad")
    print("=" * 60)
    
    adata_large = examine_brain_data("D:/openclaw/intelligence-augmentation/data/brain_scrna/Siletti_DLPFC_113k.h5ad")
    check_intelligence_genes(adata_large)