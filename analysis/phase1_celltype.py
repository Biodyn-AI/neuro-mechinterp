#!/usr/bin/env python3
"""Phase 1: Cell-type stratified re-analysis of existing perturbation results."""
import pickle, os, json, sys
import numpy as np
from datasets import load_from_disk
from collections import defaultdict
from scipy import stats

base = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
fig_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/figures'
os.makedirs(fig_dir, exist_ok=True)

ds = load_from_disk(os.path.join(base, 'brain.dataset'))
cell_types = ds['cell_type']

gene_map = {
    "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
    "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
    "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
    "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
    "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
    "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
    "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
}

# Load token dictionary to map ensembl -> token
token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
with open(token_dict_path, 'rb') as f:
    token_dict = pickle.load(f)
# Reverse: ensembl -> token
ensembl_to_token = {v: k for k, v in token_dict.items()} if isinstance(list(token_dict.keys())[0], int) else token_dict

# Check token dict structure
first_key = list(token_dict.keys())[0]
first_val = list(token_dict.values())[0]
print(f"Token dict: key type={type(first_key)}, val type={type(first_val)}")
print(f"Sample: {first_key} -> {first_val}")

# We need ensembl_id -> token_id mapping
# Let's figure out the mapping from the pickle filenames
gene_token_map = {}
for symbol in gene_map:
    pdir = os.path.join(base, f'perturb_{symbol}')
    if not os.path.exists(pdir):
        continue
    files = os.listdir(pdir)
    if files:
        # Extract token from filename: ...dict_[np.int16(XXXXX)]_raw.pickle
        fname = files[0]
        token_id = int(fname.split('(')[1].split(')')[0])
        gene_token_map[symbol] = token_id

print(f"\nGene-token mapping ({len(gene_token_map)} genes):")
for g, t in sorted(gene_token_map.items()):
    print(f"  {g}: {t}")

# For each gene, get per-cell cosine similarities and map to cell types
# Group cell types for cleaner analysis
ct_groups = {
    'Excitatory': ['glutamatergic neuron'],
    'Inhibitory': ['GABAergic interneuron', 'pvalb GABAergic cortical interneuron'],
    'Oligodendrocyte': ['oligodendrocyte'],
    'OPC': ['oligodendrocyte precursor cell'],
    'Astrocyte': ['astrocyte'],
    'Other': ['microglial cell', 'mural cell', 'macrophage'],
}
ct_to_group = {}
for group, types in ct_groups.items():
    for t in types:
        ct_to_group[t] = group

results = {}  # gene -> {celltype_group -> {shift, n, std}}

for symbol, token_id in gene_token_map.items():
    pdir = os.path.join(base, f'perturb_{symbol}')
    files = os.listdir(pdir)
    pf = os.path.join(pdir, files[0])
    with open(pf, 'rb') as f:
        data = pickle.load(f)
    
    key = (np.int16(token_id), 'cell_emb')
    cos_sims = data[key]
    
    # Find which cells express this gene
    expressing_indices = []
    for i, input_ids in enumerate(ds['input_ids']):
        if token_id in input_ids:
            expressing_indices.append(i)
    
    assert len(expressing_indices) == len(cos_sims), f"{symbol}: {len(expressing_indices)} != {len(cos_sims)}"
    
    gene_results = {}
    for idx, sim in zip(expressing_indices, cos_sims):
        ct = cell_types[idx]
        group = ct_to_group.get(ct, 'Other')
        if group not in gene_results:
            gene_results[group] = []
        gene_results[group].append(1.0 - sim)  # embedding shift
    
    results[symbol] = {}
    for group, shifts in gene_results.items():
        results[symbol][group] = {
            'mean_shift': float(np.mean(shifts)),
            'std': float(np.std(shifts)),
            'n': len(shifts),
        }

# Print results table
print("\n" + "="*100)
print("CELL-TYPE STRATIFIED PERTURBATION RESULTS")
print("="*100)

groups_order = ['Excitatory', 'Inhibitory', 'Oligodendrocyte', 'OPC', 'Astrocyte', 'Other']

# Header
header = f"{'Gene':<10}"
for g in groups_order:
    header += f" | {g:>15}"
print(header)
print("-" * len(header))

# Sort genes by overall mean shift
gene_overall = {}
for symbol in results:
    all_shifts = []
    for g in results[symbol]:
        r = results[symbol][g]
        all_shifts.extend([r['mean_shift']] * r['n'])
    gene_overall[symbol] = np.mean(all_shifts)

for symbol in sorted(gene_overall, key=gene_overall.get, reverse=True):
    row = f"{symbol:<10}"
    for g in groups_order:
        if g in results[symbol]:
            r = results[symbol][g]
            row += f" | {r['mean_shift']:.4f} (n={r['n']:>3})"
        else:
            row += f" | {'---':>15}"
    print(row)

# Rankings by cell type
print("\n" + "="*80)
print("RANKINGS BY CELL TYPE (top 10)")
print("="*80)
for group in groups_order:
    genes_in_group = [(sym, results[sym][group]['mean_shift']) 
                      for sym in results if group in results[sym] and results[sym][group]['n'] >= 5]
    genes_in_group.sort(key=lambda x: x[1], reverse=True)
    print(f"\n--- {group} ---")
    for rank, (sym, shift) in enumerate(genes_in_group[:10], 1):
        print(f"  {rank}. {sym}: {shift:.4f} (n={results[sym][group]['n']})")

# Check if rankings differ
print("\n" + "="*80)
print("RANK COMPARISON: Top 5 per cell type")
print("="*80)
for group in ['Excitatory', 'Inhibitory', 'Oligodendrocyte', 'Astrocyte']:
    genes_in_group = [(sym, results[sym][group]['mean_shift']) 
                      for sym in results if group in results[sym] and results[sym][group]['n'] >= 3]
    genes_in_group.sort(key=lambda x: x[1], reverse=True)
    top5 = [g[0] for g in genes_in_group[:5]]
    print(f"  {group}: {', '.join(top5)}")

# Save results
with open(os.path.join(base, 'celltype_stratified_results.json'), 'w') as f:
    json.dump(results, f, indent=2)

# Generate figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 3, figsize=(20, 12))
axes = axes.flatten()

for ax_idx, group in enumerate(groups_order):
    ax = axes[ax_idx]
    genes_in_group = [(sym, results[sym][group]['mean_shift'], results[sym][group]['n']) 
                      for sym in results if group in results[sym] and results[sym][group]['n'] >= 3]
    genes_in_group.sort(key=lambda x: x[1], reverse=True)
    
    if not genes_in_group:
        ax.set_title(f"{group} (insufficient data)")
        continue
    
    names = [g[0] for g in genes_in_group]
    shifts = [g[1] for g in genes_in_group]
    ns = [g[2] for g in genes_in_group]
    
    colors = ['#e74c3c' if s == max(shifts) else '#3498db' for s in shifts]
    bars = ax.barh(range(len(names)), shifts, color=colors)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('Embedding Shift')
    ax.set_title(f'{group} (n cells in group: {sum(ns)})')
    ax.invert_yaxis()

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig5_celltype_rankings.png'), dpi=150)
print(f"\nSaved fig5_celltype_rankings.png")

# Heatmap: genes x cell types
fig2, ax2 = plt.subplots(figsize=(12, 10))
top_genes = sorted(gene_overall, key=gene_overall.get, reverse=True)
main_groups = ['Excitatory', 'Inhibitory', 'Oligodendrocyte', 'Astrocyte', 'OPC']

matrix = np.zeros((len(top_genes), len(main_groups)))
for i, sym in enumerate(top_genes):
    for j, g in enumerate(main_groups):
        if g in results[sym]:
            matrix[i, j] = results[sym][g]['mean_shift']
        else:
            matrix[i, j] = np.nan

im = ax2.imshow(matrix, cmap='YlOrRd', aspect='auto')
ax2.set_xticks(range(len(main_groups)))
ax2.set_xticklabels(main_groups, rotation=45, ha='right')
ax2.set_yticks(range(len(top_genes)))
ax2.set_yticklabels(top_genes)
ax2.set_title('Embedding Shift by Gene × Cell Type')
plt.colorbar(im, ax=ax2, label='Embedding Shift')

# Add text annotations
for i in range(len(top_genes)):
    for j in range(len(main_groups)):
        if not np.isnan(matrix[i, j]):
            ax2.text(j, i, f'{matrix[i,j]:.3f}', ha='center', va='center', fontsize=7)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig6_celltype_heatmap.png'), dpi=150)
print("Saved fig6_celltype_heatmap.png")

print("\nPhase 1 complete!")
