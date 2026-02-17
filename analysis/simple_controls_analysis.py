#!/usr/bin/env python3
"""
Simple Expression-Matched Controls Analysis
"""
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path

print("=== SIMPLE EXPRESSION-MATCHED CONTROLS ANALYSIS ===")

# Set up paths
PROJ_ROOT = Path("/mnt/d/openclaw/intelligence-augmentation")
RESULTS_DIR = PROJ_ROOT / "analysis" / "results"
BRAIN_DATA = PROJ_ROOT / "data" / "brain_scrna" / "DLPFC_11k.h5ad"
EXISTING_RESULTS = RESULTS_DIR / "perturbation_results_streamlined.csv"

print(f"Project root: {PROJ_ROOT}")
print(f"Existing results: {EXISTING_RESULTS}")

# Step 1: Load existing results
print("\n1. Loading existing intelligence gene results...")
if not EXISTING_RESULTS.exists():
    print(f"ERROR: File not found: {EXISTING_RESULTS}")
    sys.exit(1)

results_df = pd.read_csv(EXISTING_RESULTS)
print(f"Loaded {len(results_df)} intelligence genes")
print("\nTop 5 genes:")
print(results_df[['gene_symbol', 'expression_frequency', 'perturbation_impact']].head())

# Step 2: Load brain dataset to get all gene expression frequencies
print(f"\n2. Loading brain dataset: {BRAIN_DATA}")
try:
    import scanpy as sc
    print("Scanpy imported successfully")
    
    adata = sc.read_h5ad(BRAIN_DATA)
    print(f"Dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Get expression frequencies
    if hasattr(adata.X, "toarray"):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)
    
    expression_freq = (X > 0).mean(axis=0)
    
    # Get gene symbols
    if 'feature_name' in adata.var.columns:
        gene_symbols = adata.var['feature_name'].values
        ensembl_ids = adata.var_names.values
    else:
        gene_symbols = adata.var_names.values
        ensembl_ids = gene_symbols
    
    expression_data = pd.DataFrame({
        'ensembl_id': ensembl_ids,
        'gene_symbol': gene_symbols,
        'expression_frequency': expression_freq
    })
    
    expression_data = expression_data.drop_duplicates(subset='gene_symbol')
    print(f"Expression data for {len(expression_data)} unique genes")
    
except Exception as e:
    print(f"ERROR loading brain data: {e}")
    sys.exit(1)

# Step 3: Select expression-matched controls
print(f"\n3. Selecting expression-matched controls...")

intel_frequencies = results_df['expression_frequency'].values
print(f"Intelligence gene expression range: {intel_frequencies.min():.4f} - {intel_frequencies.max():.4f}")

# Get available genes (for now, just use all genes from dataset)
available_genes = expression_data.copy()

# Exclude intelligence genes
intel_genes = set(results_df['gene_symbol'].str.upper())
available_genes = available_genes[~available_genes['gene_symbol'].str.upper().isin(intel_genes)]
print(f"Available genes after excluding intelligence genes: {len(available_genes)}")

# Select controls matched on expression frequency
selected_controls = []
tolerance = 0.05

np.random.seed(42)

for _, intel_row in results_df.iterrows():
    target_freq = intel_row['expression_frequency']
    
    # Find candidates within tolerance
    candidates = available_genes[
        (available_genes['expression_frequency'] >= target_freq - tolerance) &
        (available_genes['expression_frequency'] <= target_freq + tolerance) &
        (~available_genes['gene_symbol'].isin([ctrl['gene_symbol'] for ctrl in selected_controls]))
    ]
    
    # Expand tolerance if needed
    iter_tolerance = tolerance
    while len(candidates) < 1 and iter_tolerance < 0.5:
        iter_tolerance *= 2
        candidates = available_genes[
            (available_genes['expression_frequency'] >= target_freq - iter_tolerance) &
            (available_genes['expression_frequency'] <= target_freq + iter_tolerance) &
            (~available_genes['gene_symbol'].isin([ctrl['gene_symbol'] for ctrl in selected_controls]))
        ]
    
    if len(candidates) > 0:
        selected = candidates.sample(n=1, random_state=42+len(selected_controls)).iloc[0]
        selected_controls.append({
            'gene_symbol': selected['gene_symbol'],
            'ensembl_id': selected['ensembl_id'],
            'expression_frequency': selected['expression_frequency'],
            'matched_to': intel_row['gene_symbol'],
            'matched_freq': intel_row['expression_frequency'],
            'freq_diff': abs(selected['expression_frequency'] - intel_row['expression_frequency'])
        })
        
        print(f"Matched {intel_row['gene_symbol']} ({intel_row['expression_frequency']:.4f}) -> {selected['gene_symbol']} ({selected['expression_frequency']:.4f})")
    
    if len(selected_controls) >= 21:
        break

controls_df = pd.DataFrame(selected_controls[:21])
print(f"\nSelected {len(controls_df)} control genes")
print(f"Mean frequency difference: {controls_df['freq_diff'].mean():.4f}")

# Step 4: Add housekeeping genes
print(f"\n4. Adding housekeeping genes...")
housekeeping_genes = ['ACTB', 'GAPDH', 'B2M', 'RPL13A', 'YWHAZ']
housekeeping_controls = []

for gene in housekeeping_genes:
    gene_data = expression_data[
        (expression_data['gene_symbol'].str.upper() == gene.upper()) |
        (expression_data['gene_symbol'] == gene)
    ]
    
    if len(gene_data) > 0:
        row = gene_data.iloc[0]
        housekeeping_controls.append({
            'gene_symbol': gene,
            'ensembl_id': row['ensembl_id'],
            'expression_frequency': row['expression_frequency'],
            'control_type': 'housekeeping'
        })
        print(f"Added {gene} (freq: {row['expression_frequency']:.4f})")

housekeeping_df = pd.DataFrame(housekeeping_controls)

# Step 5: Create summary comparison without running Geneformer
print(f"\n5. Creating summary comparison (without Geneformer analysis)...")

# Save control gene selections
controls_df['control_type'] = 'expression_matched'
all_controls = pd.concat([controls_df, housekeeping_df], ignore_index=True)

output_file = RESULTS_DIR / "selected_control_genes.csv"
all_controls.to_csv(output_file, index=False)
print(f"Control gene selection saved: {output_file}")

# Create basic comparison of expression frequencies
print(f"\nExpression frequency comparison:")
print(f"Intelligence genes - Mean: {results_df['expression_frequency'].mean():.4f}, Range: {results_df['expression_frequency'].min():.4f}-{results_df['expression_frequency'].max():.4f}")
if len(controls_df) > 0:
    print(f"Expression-matched controls - Mean: {controls_df['expression_frequency'].mean():.4f}, Range: {controls_df['expression_frequency'].min():.4f}-{controls_df['expression_frequency'].max():.4f}")
if len(housekeeping_df) > 0:
    print(f"Housekeeping controls - Mean: {housekeeping_df['expression_frequency'].mean():.4f}, Range: {housekeeping_df['expression_frequency'].min():.4f}-{housekeeping_df['expression_frequency'].max():.4f}")

# Generate basic report
report = f"""# Expression-Matched Controls Selection Summary

## Selected Control Genes

### Expression-Matched Controls ({len(controls_df)} genes)
"""

for _, row in controls_df.iterrows():
    report += f"- **{row['gene_symbol']}** (freq: {row['expression_frequency']:.4f}) matched to **{row['matched_to']}** (freq: {row['matched_freq']:.4f}, diff: {row['freq_diff']:.4f})\n"

if len(housekeeping_df) > 0:
    report += f"\n### Housekeeping Controls ({len(housekeeping_df)} genes)\n"
    for _, row in housekeeping_df.iterrows():
        report += f"- **{row['gene_symbol']}** (freq: {row['expression_frequency']:.4f})\n"

report += f"""

## Expression Frequency Comparison

| Gene Set | Count | Mean Freq | Min Freq | Max Freq |
|----------|-------|-----------|----------|----------|
| Intelligence Genes | {len(results_df)} | {results_df['expression_frequency'].mean():.4f} | {results_df['expression_frequency'].min():.4f} | {results_df['expression_frequency'].max():.4f} |
| Expression-Matched | {len(controls_df)} | {controls_df['expression_frequency'].mean():.4f} | {controls_df['expression_frequency'].min():.4f} | {controls_df['expression_frequency'].max():.4f} |
| Housekeeping | {len(housekeeping_df)} | {housekeeping_df['expression_frequency'].mean():.4f} | {housekeeping_df['expression_frequency'].min():.4f} | {housekeeping_df['expression_frequency'].max():.4f} |

## Quality Control

- Mean matching difference (expression-matched): {controls_df['freq_diff'].mean():.4f}
- Max matching difference (expression-matched): {controls_df['freq_diff'].max():.4f}
- All control genes are unique and distinct from intelligence genes

## Next Steps

To complete the analysis, run the Geneformer perturbation pipeline on these selected control genes using the same methodology as the intelligence genes.

The control genes have been selected to match expression frequencies, ensuring that any differences in perturbation effects are not due to expression level artifacts.
"""

report_file = RESULTS_DIR / "CONTROL_GENES_SELECTION_REPORT.md"
with open(report_file, 'w') as f:
    f.write(report)

print(f"\nReport saved: {report_file}")
print(f"\nControl gene selection completed successfully!")
print(f"- Expression-matched controls: {len(controls_df)}")
print(f"- Housekeeping controls: {len(housekeeping_df)}")
print(f"- Total controls: {len(all_controls)}")
