#!/usr/bin/env python3
"""
Explore gene names in the brain dataset to find intelligence genes
"""

import scanpy as sc
import pandas as pd
import numpy as np

def main():
    print("Exploring gene names in brain dataset...")
    
    data_path = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
    
    # Load data
    adata = sc.read_h5ad(data_path)
    print(f"Dataset shape: {adata.shape}")
    
    # Intelligence genes to look for
    intelligence_genes = [
        'CADM2', 'GRIN2A', 'CAMK2A', 'MEF2C', 'APP',
        'GRIN2B', 'SCN1A', 'NRXN1', 'HOMER1', 'NEGR1',
        'SHANK3', 'BDNF', 'SYN1', 'SLC6A4', 'COMT',
        'TCF4', 'CREB1', 'FMR1', 'MAPT', 'FOXO3', 'NLGN1'
    ]
    
    print(f"Looking for {len(intelligence_genes)} intelligence genes...")
    
    # Check gene names
    gene_names = list(adata.var_names)
    print(f"Total genes in dataset: {len(gene_names)}")
    print(f"First 10 genes: {gene_names[:10]}")
    
    # Look for intelligence genes (exact match)
    found_exact = []
    for gene in intelligence_genes:
        if gene in gene_names:
            found_exact.append(gene)
            print(f"✓ Found exact: {gene}")
    
    print(f"\nFound {len(found_exact)} exact matches: {found_exact}")
    
    # Look for partial matches (case-insensitive)
    print("\nLooking for partial/case-insensitive matches...")
    found_partial = []
    
    for gene in intelligence_genes:
        if gene not in found_exact:
            # Case insensitive search
            matches = [g for g in gene_names if gene.lower() in g.lower() or g.lower() in gene.lower()]
            if matches:
                print(f"~ Partial matches for {gene}: {matches[:5]}...")  # Show first 5 matches
                found_partial.extend(matches)
    
    # Check if gene names might be Ensembl IDs
    print(f"\nFirst few gene names to check format:")
    for i in range(min(20, len(gene_names))):
        print(f"  {i+1}: {gene_names[i]}")
    
    # Check if there are any gene symbols in var
    if hasattr(adata.var, 'columns'):
        print(f"\nVariable columns available: {list(adata.var.columns)}")
        
        # Look for gene symbol columns
        symbol_cols = [col for col in adata.var.columns if 'symbol' in col.lower() or 'gene' in col.lower()]
        if symbol_cols:
            print(f"Potential gene symbol columns: {symbol_cols}")
            
            for col in symbol_cols:
                print(f"\nChecking {col} column...")
                symbols = list(adata.var[col].dropna().unique())
                print(f"First 10 symbols from {col}: {symbols[:10]}")
                
                # Check for intelligence genes in this column
                found_in_col = []
                for gene in intelligence_genes:
                    if gene in symbols:
                        found_in_col.append(gene)
                        
                print(f"Intelligence genes found in {col}: {found_in_col}")
    
    # Save gene list for reference
    gene_df = pd.DataFrame({'gene_name': gene_names})
    gene_df.to_csv('/mnt/d/openclaw/intelligence-augmentation/analysis/results/all_gene_names.csv', index=False)
    print(f"\nAll gene names saved to: all_gene_names.csv")

if __name__ == "__main__":
    main()