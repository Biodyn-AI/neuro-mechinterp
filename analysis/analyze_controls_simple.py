#!/usr/bin/env python3
"""
Simple analysis to identify control genes
This version works with existing data files and doesn't require Geneformer
"""
import pandas as pd
import numpy as np
from pathlib import Path

# Load the intelligence gene results
results_file = Path("D:/openclaw/intelligence-augmentation/analysis/results/perturbation_results_streamlined.csv")
results_df = pd.read_csv(results_file)

print("=== INTELLIGENCE GENE PERTURBATION RESULTS ===")
print(f"Number of genes: {len(results_df)}")
print(f"Expression frequency range: {results_df['expression_frequency'].min():.4f} - {results_df['expression_frequency'].max():.4f}")
print(f"Perturbation impact range: {results_df['perturbation_impact'].min():.8f} - {results_df['perturbation_impact'].max():.8f}")

print("\nTop 10 genes by perturbation impact:")
print(results_df[['gene_symbol', 'expression_frequency', 'perturbation_impact']].head(10))

print(f"\nExpression frequency distribution:")
print(f"Mean: {results_df['expression_frequency'].mean():.4f}")
print(f"Median: {results_df['expression_frequency'].median():.4f}")
print(f"Std: {results_df['expression_frequency'].std():.4f}")

# Show quartiles for expression frequency
quartiles = results_df['expression_frequency'].quantile([0.25, 0.5, 0.75])
print(f"Q1: {quartiles[0.25]:.4f}, Q2: {quartiles[0.5]:.4f}, Q3: {quartiles[0.75]:.4f}")

# Look for other data files that might contain all gene information
results_dir = Path("D:/openclaw/intelligence-augmentation/analysis/results")
print(f"\nLooking for additional data files in {results_dir}...")

for file in results_dir.glob("*.csv"):
    if file.name != "perturbation_results_streamlined.csv":
        try:
            df = pd.read_csv(file, nrows=5)  # Just peek at first few rows
            print(f"- {file.name}: {len(df.columns)} columns")
            print(f"  Columns: {list(df.columns)[:5]}...")  # Show first 5 columns
        except:
            print(f"- {file.name}: Could not read")

# Let's specifically look for the all_gene files
print(f"\n=== ANALYZING ALL_GENE FILES ===")

# Try to load the all gene mapping file
try:
    all_genes_file = results_dir / "all_gene_mapping.csv"
    if all_genes_file.exists():
        # This might be large, so sample it
        all_genes_df = pd.read_csv(all_genes_file)
        print(f"All genes mapping file: {len(all_genes_df)} rows")
        print("Columns:", list(all_genes_df.columns))
        print("First 5 rows:")
        print(all_genes_df.head())
        
        # Check if it has expression frequency information
        if 'expression_frequency' in all_genes_df.columns or 'detection_rate' in all_genes_df.columns:
            print("Found expression frequency data!")
        
except Exception as e:
    print(f"Could not load all_gene_mapping.csv: {e}")

# Try the all_gene_names file
try:
    all_names_file = results_dir / "all_gene_names.csv"
    if all_names_file.exists():
        all_names_df = pd.read_csv(all_names_file)
        print(f"\nAll gene names file: {len(all_names_df)} rows")
        print("Columns:", list(all_names_df.columns))
        print("First 5 rows:")
        print(all_names_df.head())
        
except Exception as e:
    print(f"Could not load all_gene_names.csv: {e}")

print(f"\n=== SUMMARY ===")
print(f"We need to find {len(results_df)} control genes matched on expression frequency")
print("Expression frequencies to match:")
for _, row in results_df.iterrows():
    print(f"- {row['gene_symbol']}: {row['expression_frequency']:.4f}")