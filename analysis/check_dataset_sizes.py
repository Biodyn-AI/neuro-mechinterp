#!/usr/bin/env python3
"""Check available dataset sizes for scaling analysis"""
import scanpy as sc
import pandas as pd
import numpy as np

print('Loading datasets to check sizes...')

datasets = [
    ('500-cell (current)', 'results/insilico_wsl/prepared_data/brain_500.h5ad'),
    ('DLPFC 11k', '../data/brain_scrna/DLPFC_11k.h5ad'), 
    ('Siletti 113k', '../data/brain_scrna/Siletti_DLPFC_113k.h5ad')
]

for name, path in datasets:
    try:
        adata = sc.read_h5ad(path)
        memory_mb = adata.X.nbytes / 1e6 if hasattr(adata.X, 'nbytes') else 0
        print(f'{name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes')
        print(f'  Memory: {memory_mb:.1f} MB')
        print(f'  Cell types: {adata.obs.columns.tolist()[:5]}...')
        del adata
    except Exception as e:
        print(f'{name}: ERROR - {e}')
    print()

print('Suggested scaling targets:')
print('- 1,000 cells: 2x current (safe)')
print('- 2,000 cells: 4x current (moderate)')  
print('- 5,000 cells: 10x current (aggressive, may OOM)')