#!/usr/bin/env python3
"""
Debug gene symbol mapping between GWAS genes and Geneformer tokens
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc

# Add Geneformer to path
GENEFORMER_PATH = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
sys.path.insert(0, GENEFORMER_PATH)

def investigate_gene_mapping():
    print("Gene Mapping Investigation")
    print("=" * 50)
    
    # Load token dictionary
    token_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'gene_dictionaries_30m', 'token_dictionary_gc30M.pkl')
    if not os.path.exists(token_path):
        token_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
    
    with open(token_path, 'rb') as f:
        token_dict = pickle.load(f)
    
    print(f"Token dictionary size: {len(token_dict)}")
    
    # Show sample tokens
    sample_tokens = list(token_dict.keys())[:20]
    print(f"Sample token keys: {sample_tokens}")
    
    # Load data
    data_path = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
    adata = sc.read_h5ad(data_path)
    
    # Get gene symbols from data
    if 'feature_name' in adata.var.columns:
        data_genes = adata.var['feature_name'].values
        print("Using feature_name column for data genes")
    else:
        data_genes = adata.var.index.values
        print("Using var index for data genes")
    
    print(f"Data genes count: {len(data_genes)}")
    print(f"Sample data genes: {data_genes[:20]}")
    
    # Check overlap between data genes and token dict
    overlap = set(data_genes) & set(token_dict.keys())
    print(f"Overlap between data and tokens: {len(overlap)} genes")
    
    # Our target genes
    target_genes = [
        'CADM2', 'NRXN1', 'NLGN1', 'NEGR1', 'TCF4', 'MEF2C', 'APP', 'GRIN2A', 'GRIN2B',
        'MAPT', 'SCN1A', 'DISC1', 'FOXO3', 'CAMK2A', 'HOMER1', 'CREB1', 'COMT', 'SHANK3',
        'FMR1', 'SYN1', 'BDNF', 'SLC6A4'
    ]
    
    print(f"\nAnalyzing {len(target_genes)} target genes:")
    print("-" * 50)
    
    results = []
    for gene in target_genes:
        in_data = gene in data_genes
        in_tokens = gene in token_dict
        
        # Try to find similar genes
        similar_in_data = [g for g in data_genes if gene.lower() in g.lower()][:3]
        similar_in_tokens = [g for g in token_dict.keys() if gene.lower() in g.lower()][:3]
        
        print(f"{gene:10s} | Data: {'YES' if in_data else 'NO':3s} | Tokens: {'YES' if in_tokens else 'NO':3s}")
        if similar_in_data:
            print(f"           | Similar in data: {similar_in_data}")
        if similar_in_tokens:
            print(f"           | Similar in tokens: {similar_in_tokens}")
        
        results.append({
            'gene': gene,
            'in_data': in_data,
            'in_tokens': in_tokens,
            'similar_data': similar_in_data,
            'similar_tokens': similar_in_tokens
        })
    
    # Summary stats
    in_data_count = sum(1 for r in results if r['in_data'])
    in_tokens_count = sum(1 for r in results if r['in_tokens'])
    in_both_count = sum(1 for r in results if r['in_data'] and r['in_tokens'])
    
    print(f"\nSummary:")
    print(f"  Target genes in data: {in_data_count}/{len(target_genes)}")
    print(f"  Target genes in tokens: {in_tokens_count}/{len(target_genes)}")
    print(f"  Target genes in both: {in_both_count}/{len(target_genes)}")
    
    # If genes are in data but not tokens, check for alternative symbols
    print(f"\nGenes in data but not tokens:")
    for r in results:
        if r['in_data'] and not r['in_tokens']:
            print(f"  {r['gene']} - try alternatives: {r['similar_tokens']}")
    
    # Check if we need to use Ensembl IDs instead
    print(f"\nChecking Ensembl ID format...")
    ensembl_pattern = [g for g in data_genes if g.startswith('ENSG')][:5]
    if ensembl_pattern:
        print(f"  Data contains Ensembl IDs: {ensembl_pattern}")
        
        # Check if tokens use Ensembl IDs
        token_ensembl = [g for g in token_dict.keys() if g.startswith('ENSG')][:5]
        if token_ensembl:
            print(f"  Tokens contain Ensembl IDs: {token_ensembl}")
        else:
            print(f"  Tokens do NOT contain Ensembl IDs")
    
    return results

def find_working_genes():
    """Find genes that are in both data and token dictionary for testing"""
    print(f"\nFinding working test genes...")
    
    # Load token dictionary
    token_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'gene_dictionaries_30m', 'token_dictionary_gc30M.pkl')
    if not os.path.exists(token_path):
        token_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
    
    with open(token_path, 'rb') as f:
        token_dict = pickle.load(f)
    
    # Load data
    data_path = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
    adata = sc.read_h5ad(data_path)
    
    if 'feature_name' in adata.var.columns:
        data_genes = set(adata.var['feature_name'].values)
    else:
        data_genes = set(adata.var.index.values)
    
    # Find intersection
    working_genes = list(data_genes & set(token_dict.keys()))
    
    print(f"Found {len(working_genes)} working genes")
    
    # Sample a few and check their expression
    if len(working_genes) > 0:
        print(f"Sample working genes: {working_genes[:10]}")
        
        # Check expression levels
        sample_genes = working_genes[:5]
        print(f"\nExpression check for sample genes:")
        
        for gene in sample_genes:
            gene_mask = adata.var['feature_name'] == gene if 'feature_name' in adata.var.columns else adata.var.index == gene
            if np.any(gene_mask):
                gene_idx = np.where(gene_mask)[0][0]
                expr_vals = adata[:100, gene_idx].X.toarray().flatten()  # First 100 cells
                expressing = np.sum(expr_vals > 0)
                mean_expr = np.mean(expr_vals[expr_vals > 0]) if expressing > 0 else 0
                print(f"  {gene}: {expressing}/100 cells, mean={mean_expr:.3f}")
    
    return working_genes

if __name__ == "__main__":
    results = investigate_gene_mapping()
    working_genes = find_working_genes()
    
    print(f"\nRECOMMENDATIONS:")
    print(f"1. Use working genes for initial testing: {working_genes[:5]}")
    print(f"2. Map target gene symbols to token dictionary format")
    print(f"3. Check if target genes use different symbol conventions")