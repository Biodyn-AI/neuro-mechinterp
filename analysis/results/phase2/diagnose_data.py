#!/usr/bin/env python3
"""Diagnose data structure and gene name formats."""

import pandas as pd
import scanpy as sc
from pathlib import Path

# Paths
DATA_DIR = Path(r'D:\openclaw\intelligence-augmentation\data')
BRAIN_FILE = DATA_DIR / 'brain_scrna' / 'DLPFC_11k.h5ad'
DRUG_DB_FILE = DATA_DIR / 'drugs' / 'interactions.tsv'

# GWAS genes we're looking for
GWAS_GENES = [
    'CADM2', 'APP', 'FOXO3', 'NLGN1', 'MAPT', 'TCF4', 'CREB1', 'DISC1', 
    'NEGR1', 'GRIN2A', 'FMR1', 'COMT', 'BDNF', 'CAMK2A', 'HOMER1', 
    'SHANK3', 'NRXN1', 'MEF2C', 'SLC6A4', 'SCN1A', 'GRIN2B', 'SYN1'
]

print("=" * 60)
print("PHASE 2 DATA STRUCTURE DIAGNOSIS")
print("=" * 60)

# 1. Check brain dataset structure
print("\n1. BRAIN DATASET ANALYSIS")
print("-" * 30)

try:
    adata = sc.read_h5ad(BRAIN_FILE)
    print(f"Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    
    print(f"\nObservation (cell) columns: {list(adata.obs.columns)}")
    print(f"Variable (gene) columns: {list(adata.var.columns)}")
    
    # Check gene names
    print(f"\nFirst 10 gene names: {list(adata.var.index[:10])}")
    
    # Check if gene names are symbols or IDs
    gene_names = adata.var.index.tolist()
    symbols_like = sum(1 for g in gene_names[:100] if g.isalpha() and g.isupper())
    ensembl_like = sum(1 for g in gene_names[:100] if g.startswith('ENSG'))
    
    print(f"First 100 genes analysis:")
    print(f"  Symbol-like (uppercase letters): {symbols_like}/100")
    print(f"  Ensembl-like (starts with ENSG): {ensembl_like}/100")
    
    # Check if GWAS genes are present (exact match)
    found_exact = [g for g in GWAS_GENES if g in gene_names]
    print(f"\nGWAS genes found (exact match): {len(found_exact)}/{len(GWAS_GENES)}")
    if found_exact:
        print(f"  Found: {found_exact}")
    
    # Check for case-insensitive matches
    gene_names_upper = [g.upper() for g in gene_names]
    found_case_insensitive = [g for g in GWAS_GENES if g.upper() in gene_names_upper]
    print(f"GWAS genes found (case insensitive): {len(found_case_insensitive)}/{len(GWAS_GENES)}")
    if found_case_insensitive:
        print(f"  Found: {found_case_insensitive}")
    
    # Check for partial matches
    print(f"\nPartial matches for GWAS genes:")
    for gwas_gene in GWAS_GENES[:5]:  # Check first 5
        matches = [g for g in gene_names if gwas_gene.lower() in g.lower()]
        if matches:
            print(f"  {gwas_gene}: {matches[:3]}")  # Show first 3 matches
    
    # Check if var has gene_names or gene_ids columns
    if 'gene_names' in adata.var.columns:
        print(f"\nFound 'gene_names' column in var:")
        print(f"  First 10 gene names: {list(adata.var['gene_names'][:10])}")
        
        # Check GWAS genes in gene_names column
        gene_symbols = adata.var['gene_names'].tolist()
        found_in_symbols = [g for g in GWAS_GENES if g in gene_symbols]
        print(f"  GWAS genes found in gene_names: {len(found_in_symbols)}/{len(GWAS_GENES)}")
        if found_in_symbols:
            print(f"    Found: {found_in_symbols}")
    
    # Check cell type information
    print(f"\nCell type information:")
    cell_type_cols = [col for col in adata.obs.columns if 'type' in col.lower() or 'cluster' in col.lower() or 'class' in col.lower()]
    print(f"  Potential cell type columns: {cell_type_cols}")
    
    for col in cell_type_cols[:2]:  # Check first 2 columns
        if col in adata.obs.columns:
            unique_values = adata.obs[col].unique()
            print(f"  {col}: {len(unique_values)} unique values")
            print(f"    Examples: {list(unique_values[:5])}")

except Exception as e:
    print(f"Error loading brain dataset: {e}")
    import traceback
    traceback.print_exc()

# 2. Check drug database structure
print(f"\n2. DRUG DATABASE ANALYSIS")
print("-" * 30)

try:
    drug_db = pd.read_csv(DRUG_DB_FILE, sep='\t', nrows=5)  # Load just first 5 rows to check structure
    print(f"Columns: {list(drug_db.columns)}")
    print(f"\nFirst few rows:")
    print(drug_db.head())
    
    # Load full database to check gene names
    print(f"\nChecking full database...")
    drug_db_full = pd.read_csv(DRUG_DB_FILE, sep='\t', low_memory=False)
    print(f"Total interactions: {len(drug_db_full):,}")
    
    # Find the gene column
    gene_cols = [col for col in drug_db_full.columns if 'gene' in col.lower()]
    print(f"Gene-related columns: {gene_cols}")
    
    if gene_cols:
        gene_col = gene_cols[0]
        unique_genes = drug_db_full[gene_col].nunique()
        print(f"Unique genes in {gene_col}: {unique_genes:,}")
        
        # Check for GWAS genes in drug database
        gene_names_in_db = drug_db_full[gene_col].astype(str).str.upper().unique()
        found_in_db = [g for g in GWAS_GENES if g.upper() in gene_names_in_db]
        print(f"GWAS genes found in drug DB: {len(found_in_db)}/{len(GWAS_GENES)}")
        if found_in_db:
            print(f"  Found: {found_in_db}")

except Exception as e:
    print(f"Error loading drug database: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 60)
print("DIAGNOSIS COMPLETE")
print("=" * 60)