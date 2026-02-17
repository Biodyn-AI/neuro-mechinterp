#!/usr/bin/env python3
"""Analyze baseline 500-cell results in detail to prepare for comparison"""
import json
import pandas as pd
import numpy as np
from pathlib import Path

print("Analyzing baseline 500-cell results...")

# Load final rankings
results_dir = Path("results/insilico_wsl")
rankings_file = results_dir / "final_rankings.json"

if not rankings_file.exists():
    print(f"Rankings file not found: {rankings_file}")
    exit(1)

with open(rankings_file, 'r') as f:
    rankings = json.load(f)

print(f"Loaded rankings for {len(rankings)} genes")

# Focus on top 5 genes
TOP_GENES = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]

print("\n=== TOP 5 GENE ANALYSIS ===")
print("Gene\tShift Mean\tShift Std\tN Cells\tPathway")
print("-" * 60)

baseline_summary = {}
for gene in TOP_GENES:
    if gene in rankings:
        data = rankings[gene]
        print(f"{gene}\t{data['shift_mean']:.4f}\t{data['shift_std']:.4f}\t{data['n_cells']}\t{data['pathway']}")
        
        baseline_summary[gene] = {
            'shift_mean': data['shift_mean'],
            'shift_std': data['shift_std'],
            'n_cells': data['n_cells'],
            'cosine_sim_mean': data['cosine_sim_mean'],
            'pathway': data['pathway']
        }
    else:
        print(f"{gene}\tNOT FOUND")

# Save baseline summary for comparison
baseline_file = Path("results/baseline_500_summary.json")
with open(baseline_file, 'w') as f:
    json.dump(baseline_summary, f, indent=2)

print(f"\nBaseline summary saved: {baseline_file}")

# Calculate some statistics
shift_means = [data['shift_mean'] for data in baseline_summary.values()]
print(f"\nStatistics on top 5 genes (500 cells):")
print(f"Mean shift effect: {np.mean(shift_means):.4f}")
print(f"Std shift effect: {np.std(shift_means):.4f}")
print(f"Range: {np.min(shift_means):.4f} - {np.max(shift_means):.4f}")

# Check if we have individual gene result files
individual_results = {}
for gene in TOP_GENES:
    gene_dir = results_dir / f"perturb_{gene}"
    if gene_dir.exists():
        pickle_files = list(gene_dir.glob("*.pickle"))
        individual_results[gene] = len(pickle_files)
        print(f"{gene}: {len(pickle_files)} result files in {gene_dir}")

print(f"\nReady for comparison with larger datasets!")
print(f"Key metrics to compare:")
print("1. Ranking stability (do top genes stay on top?)")
print("2. Effect size scaling (do shift_mean values change?)")
print("3. Variance patterns (does std change with more cells?)")
print("4. Processing efficiency (time/memory requirements)")