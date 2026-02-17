#!/usr/bin/env python3
"""
Simple data structure check for brain scRNA-seq data.
"""

import scanpy as sc
import pandas as pd
import numpy as np

print("Loading DLPFC_11k.h5ad...")
adata = sc.read_h5ad("D:/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")

print(f"Data shape: {adata.shape} (cells x genes)")
print(f"Obs columns: {list(adata.obs.columns)[:10]}...")  # First 10
print(f"Var columns: {list(adata.var.columns)[:10]}...")  # First 10

# Check gene identifiers
print(f"Gene index (first 5): {adata.var.index[:5].tolist()}")
print(f"Are gene indices Ensembl IDs? {adata.var.index[0].startswith('ENSG')}")

# Check for counts
if 'sum' in adata.obs.columns:
    print(f"Found 'sum' column in obs (likely total counts)")
    print(f"Count range: {adata.obs['sum'].min():.0f} - {adata.obs['sum'].max():.0f}")

# Check data matrix
print(f"X matrix type: {type(adata.X)}")
print(f"X dtype: {adata.X.dtype}")

# Sample values
sample = adata.X[:5, :5]
if hasattr(sample, 'toarray'):
    sample = sample.toarray()
print(f"Sample values: {sample}")

# Check intelligence genes
intelligence_genes = ["MEF2C", "BDNF", "GRIN2B", "CADM2", "NRXN1"]
print(f"\nChecking for intelligence genes in gene index...")

# Try matching gene symbols if they exist
if 'feature_name' in adata.var.columns:
    gene_symbols = adata.var['feature_name']
    found = [gene for gene in intelligence_genes if gene in gene_symbols.values]
    print(f"Found in feature_name: {found}")

# Also check var index directly
found_in_index = [gene for gene in intelligence_genes if gene in adata.var.index]
print(f"Found in var index: {found_in_index}")

print("Done!")