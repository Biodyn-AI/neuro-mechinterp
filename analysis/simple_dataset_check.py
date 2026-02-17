#!/usr/bin/env python3
import scanpy as sc
import numpy as np
import os

print("Simple dataset size check...")

try:
    print("Loading 11k dataset...")
    adata = sc.read_h5ad('../data/brain_scrna/DLPFC_11k.h5ad')
    print(f"Dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    # Test subsample to 1000 cells
    print("Testing subsample to 1000 cells...")
    adata_sub = adata.copy()
    sc.pp.subsample(adata_sub, n_obs=1000, random_state=42)
    print(f"Subsampled: {adata_sub.n_obs:,} cells")
    
    # Save test dataset
    os.makedirs('results/scaled', exist_ok=True)
    adata_sub.write_h5ad('results/scaled/brain_1000_test.h5ad')
    print("Saved test dataset: results/scaled/brain_1000_test.h5ad")
    
    print("SUCCESS - Dataset loading and subsampling works!")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()