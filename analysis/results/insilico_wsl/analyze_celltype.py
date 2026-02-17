#!/usr/bin/env python3
"""Analyze cell-type stratified results and produce summary"""
import json
import numpy as np

with open("/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/celltype_stratified_results.json") as f:
    data = json.load(f)

# Get all cell types
all_types = set()
for gene, types in data.items():
    all_types.update(types.keys())
print(f"Cell types found: {sorted(all_types)}")

# Build rankings per cell type
rankings_by_type = {}
for ct in sorted(all_types):
    genes_in_type = []
    for gene, types in data.items():
        if ct in types:
            genes_in_type.append((gene, types[ct]["mean_shift"], types[ct]["n"], types[ct].get("std", 0)))
    genes_in_type.sort(key=lambda x: x[1], reverse=True)
    rankings_by_type[ct] = genes_in_type

print("\n=== TOP 5 GENES BY CELL TYPE ===")
for ct, ranked in rankings_by_type.items():
    total_n = sum(x[2] for x in ranked)
    print(f"\n{ct} (total cells: {total_n}):")
    for i, (gene, shift, n, std) in enumerate(ranked[:5]):
        print(f"  {i+1}. {gene}: shift={shift:.4f} ± {std:.4f} (n={n})")

# Overall ranking comparison
print("\n=== RANK COMPARISON ACROSS CELL TYPES ===")
header = f"{'Gene':<10}"
for ct in sorted(all_types):
    header += f" {ct[:8]:>8}"
print(header)
print("-" * len(header))

# Get all genes
all_genes = sorted(data.keys())
for gene in all_genes:
    row = f"{gene:<10}"
    for ct in sorted(all_types):
        ranked = rankings_by_type[ct]
        rank = None
        for i, (g, s, n, sd) in enumerate(ranked):
            if g == gene:
                rank = i + 1
                break
        row += f" {rank if rank else '-':>8}"
    print(row)

# Save summary
summary = {
    "cell_types": sorted(list(all_types)),
    "rankings_by_type": {ct: [(g, s, n) for g, s, n, sd in ranked[:10]] for ct, ranked in rankings_by_type.items()},
    "top_gene_by_type": {ct: ranked[0][0] if ranked else None for ct, ranked in rankings_by_type.items()},
}
with open("/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/celltype_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nSaved celltype_summary.json")
