#!/usr/bin/env python3
"""
Map intelligence genes to Ensembl IDs and create basic analysis.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import pickle

print("Loading brain data...")
adata = sc.read_h5ad("D:/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")
print(f"Loaded data: {adata.shape}")

# Intelligence genes from Savage et al. 2018
intelligence_genes = [
    "MEF2C", "BDNF", "GRIN2B", "CADM2", "NRXN1", "CAMK2A", 
    "GRIN2A", "SHANK3", "HOMER1", "APP", "NEGR1", "NLGN1", 
    "TCF4", "MAPT", "FOXO3", "CREB1", "FMR1", "SYN1", 
    "SCN1A", "SLC6A4", "COMT"
]

print(f"Looking for {len(intelligence_genes)} intelligence genes...")

# Create mapping from gene symbols to Ensembl IDs
gene_mapping = {}
gene_stats = []

if 'feature_name' in adata.var.columns:
    print("Found feature_name column in var metadata")
    
    for gene in intelligence_genes:
        # Find rows where feature_name matches
        matches = adata.var[adata.var['feature_name'] == gene]
        
        if len(matches) > 0:
            ensembl_id = matches.index[0]  # Get the Ensembl ID
            gene_mapping[gene] = ensembl_id
            
            # Get expression statistics
            gene_idx = adata.var.index.get_loc(ensembl_id)
            if hasattr(adata.X, 'toarray'):
                expression = adata.X[:, gene_idx].toarray().flatten()
            else:
                expression = adata.X[:, gene_idx].flatten()
            
            stats = {
                'gene_symbol': gene,
                'ensembl_id': ensembl_id,
                'mean_expression': float(np.mean(expression)),
                'std_expression': float(np.std(expression)),
                'percent_expressed': float((expression > 0).sum() / len(expression) * 100),
                'max_expression': float(np.max(expression)),
                'median_expression': float(np.median(expression))
            }
            
            gene_stats.append(stats)
            print(f"Found {gene}: {ensembl_id} (mean expr: {stats['mean_expression']:.2f})")
        else:
            print(f"NOT FOUND: {gene}")

print(f"\nMapped {len(gene_mapping)} out of {len(intelligence_genes)} intelligence genes")

# Save the mapping
with open("D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/intelligence_genes_mapping.pkl", 'wb') as f:
    pickle.dump(gene_mapping, f)

# Create and save statistics dataframe
stats_df = pd.DataFrame(gene_stats)
stats_df = stats_df.sort_values('mean_expression', ascending=False)

print("\nTop intelligence genes by mean expression:")
print(stats_df[['gene_symbol', 'mean_expression', 'percent_expressed']].to_string(index=False))

stats_df.to_csv("D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/intelligence_genes_stats.csv", index=False)

print(f"\nSaved mapping and statistics to insilico_perturber directory")
print(f"Total genes mapped: {len(gene_mapping)}")

# Also save a simple text summary
with open("D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/gene_mapping_summary.txt", 'w') as f:
    f.write(f"Intelligence Gene Mapping Summary\n")
    f.write(f"==================================\n\n")
    f.write(f"Total genes searched: {len(intelligence_genes)}\n")
    f.write(f"Total genes found: {len(gene_mapping)}\n")
    f.write(f"Success rate: {len(gene_mapping)/len(intelligence_genes)*100:.1f}%\n\n")
    
    f.write("Gene Symbol -> Ensembl ID Mapping:\n")
    f.write("-" * 50 + "\n")
    for gene, ensembl in gene_mapping.items():
        f.write(f"{gene:12} -> {ensembl}\n")
    
    f.write(f"\nMissing genes:\n")
    f.write("-" * 20 + "\n")
    missing = set(intelligence_genes) - set(gene_mapping.keys())
    for gene in missing:
        f.write(f"{gene}\n")

print("Analysis complete!")