"""
Comprehensive analysis of Geneformer in silico perturbation results
for 21 intelligence-associated genes.
"""
import os, pickle, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

base = r"D:\openclaw\intelligence-augmentation\analysis\results\insilico_wsl"
fig_dir = r"D:\openclaw\intelligence-augmentation\analysis\figures"
os.makedirs(fig_dir, exist_ok=True)

# Pathway annotations
PATHWAYS = {
    'Synaptic Scaffolding': ['SHANK3', 'HOMER1', 'NLGN1', 'NRXN1'],
    'Glutamate Receptors': ['GRIN2A', 'GRIN2B'],
    'Synaptic Vesicle/Transmission': ['SYN1', 'SLC6A4', 'SCN1A'],
    'Neurotrophic Signaling': ['BDNF', 'CREB1', 'CAMK2A'],
    'Cell Adhesion': ['CADM2', 'NEGR1'],
    'Neurodevelopmental/Transcription': ['TCF4', 'MEF2C', 'FMR1'],
    'Neurodegeneration-Associated': ['APP', 'MAPT', 'FOXO3'],
    'Catecholamine Metabolism': ['COMT'],
}

# Gene to pathway reverse map
gene_to_pathway = {}
for pw, genes in PATHWAYS.items():
    for g in genes:
        gene_to_pathway[g] = pw

# Load all results
results = {}
raw_values = {}

for d in sorted(os.listdir(base)):
    if not d.startswith("perturb_"):
        continue
    gene = d.replace("perturb_", "")
    path = os.path.join(base, d)
    files = [f for f in os.listdir(path) if f.endswith('.pickle')]
    if not files:
        continue
    with open(os.path.join(path, files[0]), "rb") as f:
        data = pickle.load(f)
    
    for key, val in data.items():
        if isinstance(val, list):
            arr = np.array(val)
            # Values are cosine similarities between original and perturbed embeddings
            # Higher similarity = less perturbation effect
            # We compute "embedding shift" = 1 - cosine_similarity for interpretability
            shifts = 1.0 - arr
            raw_values[gene] = arr  # raw cosine similarities
            results[gene] = {
                "cosine_sim_mean": float(np.mean(arr)),
                "cosine_sim_std": float(np.std(arr)),
                "shift_mean": float(np.mean(shifts)),
                "shift_std": float(np.std(shifts)),
                "n_cells": len(arr),
                "pathway": gene_to_pathway.get(gene, "Unknown"),
            }

# Rank by cosine similarity (lower = bigger effect) or equivalently by shift (higher = bigger)
ranked = sorted(results.items(), key=lambda x: x[1]["shift_mean"], reverse=True)

print("=" * 70)
print("GENEFORMER IN SILICO PERTURBATION — FULL RANKED RESULTS")
print("=" * 70)
print(f"{'Rank':<5} {'Gene':<10} {'Cos Sim (mean±std)':<22} {'Shift':<10} {'N':>5}  {'Pathway'}")
print("-" * 70)
for i, (gene, info) in enumerate(ranked, 1):
    print(f"{i:<5} {gene:<10} {info['cosine_sim_mean']:.4f}±{info['cosine_sim_std']:.4f}    {info['shift_mean']:.4f}    {info['n_cells']:>5}  {info['pathway']}")

# Statistical tests: one-sample t-test (is shift > 0?)
print("\n=== STATISTICAL TESTS ===")
print(f"{'Gene':<10} {'t-stat':>8} {'p-value':>12} {'Significant':>12}")
for gene, vals in sorted(raw_values.items()):
    shifts = 1.0 - vals
    t_stat, p_val = stats.ttest_1samp(shifts, 0)
    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    print(f"{gene:<10} {t_stat:>8.2f} {p_val:>12.2e} {sig:>12}")

# Pairwise comparisons between top genes
print("\n=== TOP 5 PAIRWISE COMPARISONS (Mann-Whitney U) ===")
top5 = [g for g, _ in ranked[:5]]
for i in range(len(top5)):
    for j in range(i+1, len(top5)):
        g1, g2 = top5[i], top5[j]
        u_stat, p_val = stats.mannwhitneyu(1-raw_values[g1], 1-raw_values[g2], alternative='two-sided')
        print(f"  {g1} vs {g2}: U={u_stat:.0f}, p={p_val:.4f}")

# Pathway-level analysis
print("\n=== PATHWAY ANALYSIS ===")
pathway_shifts = {}
for pw, genes in PATHWAYS.items():
    pw_vals = []
    for g in genes:
        if g in raw_values:
            pw_vals.extend((1.0 - raw_values[g]).tolist())
    if pw_vals:
        pathway_shifts[pw] = np.array(pw_vals)
        print(f"{pw:<35} mean_shift={np.mean(pw_vals):.4f}  n_genes={len([g for g in genes if g in raw_values])}")

# ============ FIGURES ============

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150, 'savefig.bbox': 'tight'})
colors_pathway = {
    'Synaptic Scaffolding': '#e74c3c',
    'Glutamate Receptors': '#e67e22', 
    'Synaptic Vesicle/Transmission': '#f1c40f',
    'Neurotrophic Signaling': '#2ecc71',
    'Cell Adhesion': '#1abc9c',
    'Neurodevelopmental/Transcription': '#3498db',
    'Neurodegeneration-Associated': '#9b59b6',
    'Catecholamine Metabolism': '#95a5a6',
}

# Figure 1: Ranked bar chart
fig, ax = plt.subplots(figsize=(12, 7))
genes_ranked = [g for g, _ in ranked]
shifts_ranked = [info["shift_mean"] for _, info in ranked]
stds_ranked = [info["shift_std"] for _, info in ranked]
bar_colors = [colors_pathway.get(results[g]["pathway"], '#666') for g in genes_ranked]

bars = ax.bar(range(len(genes_ranked)), shifts_ranked, yerr=stds_ranked, 
              color=bar_colors, edgecolor='black', linewidth=0.5, capsize=3)
ax.set_xticks(range(len(genes_ranked)))
ax.set_xticklabels(genes_ranked, rotation=45, ha='right', fontweight='bold')
ax.set_ylabel('Embedding Shift (1 - cosine similarity)', fontsize=12)
ax.set_title('Geneformer In Silico Perturbation: Intelligence Gene Ranking\n(Deletion perturbation on 500 brain cells)', fontsize=13)
ax.axhline(y=np.mean(shifts_ranked), color='red', linestyle='--', alpha=0.5, label=f'Mean={np.mean(shifts_ranked):.4f}')
# Legend for pathways
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=pw) for pw, c in colors_pathway.items()]
ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=1)
ax.set_xlim(-0.5, len(genes_ranked)-0.5)
plt.savefig(os.path.join(fig_dir, 'fig1_gene_ranking.png'))
plt.close()

# Figure 2: Pathway-level box plot
fig, ax = plt.subplots(figsize=(12, 6))
pw_order = sorted(pathway_shifts.keys(), key=lambda pw: np.mean(pathway_shifts[pw]), reverse=True)
bp_data = [pathway_shifts[pw] for pw in pw_order]
bp_colors = [colors_pathway.get(pw, '#666') for pw in pw_order]
bp = ax.boxplot(bp_data, labels=[pw.replace('/', '/\n') for pw in pw_order], patch_artist=True, showmeans=True)
for patch, color in zip(bp['boxes'], bp_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_ylabel('Embedding Shift (1 - cosine similarity)')
ax.set_title('Perturbation Effect by Functional Pathway')
plt.xticks(rotation=30, ha='right', fontsize=9)
plt.savefig(os.path.join(fig_dir, 'fig2_pathway_boxplot.png'))
plt.close()

# Figure 3: Violin plot of all genes
fig, ax = plt.subplots(figsize=(14, 6))
all_shifts = [1.0 - raw_values[g] for g in genes_ranked]
parts = ax.violinplot(all_shifts, positions=range(len(genes_ranked)), showmeans=True, showmedians=True)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(colors_pathway.get(results[genes_ranked[i]]["pathway"], '#666'))
    pc.set_alpha(0.7)
ax.set_xticks(range(len(genes_ranked)))
ax.set_xticklabels(genes_ranked, rotation=45, ha='right', fontweight='bold')
ax.set_ylabel('Embedding Shift')
ax.set_title('Distribution of Perturbation Effects Across Cells')
plt.savefig(os.path.join(fig_dir, 'fig3_violin_distributions.png'))
plt.close()

# Figure 4: Heatmap - gene x summary stats
fig, ax = plt.subplots(figsize=(8, 10))
stats_matrix = np.array([[results[g]["shift_mean"], results[g]["shift_std"], results[g]["cosine_sim_mean"]] for g in genes_ranked])
# Normalize columns for heatmap
stats_norm = (stats_matrix - stats_matrix.min(axis=0)) / (stats_matrix.max(axis=0) - stats_matrix.min(axis=0))
sns.heatmap(stats_norm, annot=stats_matrix, fmt='.4f', 
            xticklabels=['Mean Shift', 'Std Dev', 'Cosine Sim'],
            yticklabels=genes_ranked, cmap='YlOrRd', ax=ax)
ax.set_title('Gene Perturbation Summary Statistics')
plt.savefig(os.path.join(fig_dir, 'fig4_heatmap.png'))
plt.close()

# Save JSON
with open(os.path.join(base, "final_rankings.json"), "w") as f:
    json.dump({g: info for g, info in ranked}, f, indent=2)

# Save for report
report_data = {
    "ranked": ranked,
    "pathway_means": {pw: float(np.mean(v)) for pw, v in pathway_shifts.items()},
}
with open(os.path.join(fig_dir, "analysis_data.json"), "w") as f:
    json.dump(report_data, f, indent=2, default=str)

print(f"\nFigures saved to {fig_dir}")
print("Analysis complete!")
