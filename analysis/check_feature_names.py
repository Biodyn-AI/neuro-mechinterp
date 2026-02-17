#!/usr/bin/env python3
"""
Check the feature_name column for gene symbols
"""

import scanpy as sc
import pandas as pd

def main():
    print("Checking feature_name column for gene symbols...")
    
    data_path = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
    
    # Load data
    adata = sc.read_h5ad(data_path)
    
    # Check feature_name column
    if 'feature_name' in adata.var.columns:
        feature_names = list(adata.var['feature_name'])
        print(f"Found {len(feature_names)} feature names")
        print(f"First 20 feature names: {feature_names[:20]}")
        
        # Intelligence genes to look for
        intelligence_genes = [
            'CADM2', 'GRIN2A', 'CAMK2A', 'MEF2C', 'APP',
            'GRIN2B', 'SCN1A', 'NRXN1', 'HOMER1', 'NEGR1',
            'SHANK3', 'BDNF', 'SYN1', 'SLC6A4', 'COMT',
            'TCF4', 'CREB1', 'FMR1', 'MAPT', 'FOXO3', 'NLGN1'
        ]
        
        # Check for exact matches
        found_genes = []
        for gene in intelligence_genes:
            if gene in feature_names:
                idx = feature_names.index(gene)
                ensembl_id = adata.var_names[idx]
                found_genes.append((gene, ensembl_id))
                print(f"✓ Found {gene} -> {ensembl_id}")
        
        print(f"\nTotal found: {len(found_genes)} genes")
        
        # Save mapping
        if found_genes:
            mapping_df = pd.DataFrame(found_genes, columns=['gene_symbol', 'ensembl_id'])
            mapping_df.to_csv('/mnt/d/openclaw/intelligence-augmentation/analysis/results/gene_symbol_mapping.csv', index=False)
            print(f"Gene mapping saved to: gene_symbol_mapping.csv")
        
        # Also save all feature names for reference
        feature_df = pd.DataFrame({
            'ensembl_id': adata.var_names,
            'feature_name': feature_names
        })
        feature_df.to_csv('/mnt/d/openclaw/intelligence-augmentation/analysis/results/all_gene_mapping.csv', index=False)
        print(f"All gene mapping saved to: all_gene_mapping.csv")
        
    else:
        print("No feature_name column found!")

if __name__ == "__main__":
    main()