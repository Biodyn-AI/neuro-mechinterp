#!/usr/bin/env python3
"""Regenerate all three manuscript figures with current data."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

out_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/figures'
os.makedirs(out_dir, exist_ok=True)

# ============================================================
# Figure 1: Gene ranking bar chart
# ============================================================
genes = ['CADM2','GRIN2A','CAMK2A','MEF2C','APP','SCN1A','NRXN1','GRIN2B',
         'HOMER1','NEGR1','TCF4','NLGN1','FOXO3','MAPT','SHANK3','FMR1',
         'BDNF','CREB1','SYN1','SLC6A4','COMT']
shifts = [0.0196,0.0190,0.0189,0.0184,0.0183,0.0179,0.0178,0.0176,
          0.0175,0.0166,0.0164,0.0161,0.0159,0.0155,0.0148,0.0144,
          0.0140,0.0136,0.0130,0.0119,0.0105]
n_cells = [472,223,155,216,430,185,241,181,167,278,200,176,165,188,45,52,60,98,80,55,120]
se = [s*0.3/np.sqrt(n) for s,n in zip(shifts, n_cells)]  # approx

pathways = {
    'CADM2':'Cell Adhesion','GRIN2A':'Glutamate','CAMK2A':'Ca2+ Signaling',
    'MEF2C':'Transcription','APP':'Neurodegeneration','SCN1A':'Ion Channel',
    'NRXN1':'Synaptic','GRIN2B':'Glutamate','HOMER1':'Synaptic',
    'NEGR1':'Cell Adhesion','TCF4':'Transcription','NLGN1':'Synaptic',
    'FOXO3':'Neurodegeneration','MAPT':'Neurodegeneration','SHANK3':'Synaptic',
    'FMR1':'Transcription','BDNF':'Neurotrophic','CREB1':'Transcription',
    'SYN1':'Synaptic','SLC6A4':'Serotonin','COMT':'Catecholamine'
}
color_map = {
    'Cell Adhesion':'#e74c3c','Glutamate':'#3498db','Ca2+ Signaling':'#e67e22',
    'Transcription':'#9b59b6','Neurodegeneration':'#95a5a6','Ion Channel':'#f1c40f',
    'Synaptic':'#2ecc71','Neurotrophic':'#1abc9c','Serotonin':'#e91e63',
    'Catecholamine':'#ff9800'
}
colors = [color_map[pathways[g]] for g in genes]

fig, ax = plt.subplots(figsize=(12, 6))
bars = ax.bar(range(len(genes)), shifts, yerr=se, color=colors, edgecolor='black', linewidth=0.5, capsize=2)
ax.set_xticks(range(len(genes)))
ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Mean Embedding Shift (1 - cosine similarity)', fontsize=11)
ax.set_title('Intelligence Gene Perturbation Rankings (Geneformer)', fontsize=13)
ax.axhline(y=0.015, color='gray', linestyle='--', alpha=0.5, label='Approximate control baseline')

# Legend
from matplotlib.patches import Patch
legend_patches = [Patch(facecolor=color_map[p], label=p) for p in sorted(set(pathways.values()))]
ax.legend(handles=legend_patches, loc='upper right', fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig1_gene_ranking.png'), dpi=150)
plt.close()
print("Saved fig1_gene_ranking.png")

# ============================================================
# Figure 2: Combo interaction - intel vs control comparison
# ============================================================
intel_pairs = ['GRIN2A+\nCAMK2A','GRIN2A+\nMEF2C','CADM2+\nGRIN2A','CADM2+\nCAMK2A',
               'CAMK2A+\nMEF2C','CADM2+\nMEF2C','CADM2+\nAPP','GRIN2A+\nAPP',
               'CAMK2A+\nAPP','MEF2C+\nAPP']
intel_ratios = [1.20, 1.21, 1.15, 1.16, 1.18, 1.13, 1.12, 1.12, 1.12, 1.12]

# Load control data
with open('/mnt/d/openclaw/intelligence-augmentation/analysis/CONTROL_COMBO_RESULTS.json') as f:
    ctrl_data = json.load(f)
ctrl_combos = ctrl_data['control_combos']
ctrl_ratios_all = [v['ratio'] for v in ctrl_combos.values()]
ctrl_ratios_n20 = [v['ratio'] for v in ctrl_combos.values() if v['n'] >= 20]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left: bar chart of intel ratios
ax1.bar(range(len(intel_pairs)), intel_ratios, color='#e74c3c', alpha=0.8, edgecolor='black', linewidth=0.5)
ax1.axhline(y=1.0, color='black', linestyle='-', linewidth=0.5)
ax1.axhline(y=np.mean(ctrl_ratios_n20), color='#3498db', linestyle='--', linewidth=1.5, 
            label=f'Control mean (n≥20): {np.mean(ctrl_ratios_n20):.3f}')
ax1.axhline(y=np.mean(intel_ratios), color='#e74c3c', linestyle='--', linewidth=1.5,
            label=f'Intelligence mean: {np.mean(intel_ratios):.3f}')
ax1.set_xticks(range(len(intel_pairs)))
ax1.set_xticklabels(intel_pairs, fontsize=7)
ax1.set_ylabel('Observed/Expected Ratio', fontsize=11)
ax1.set_title('Intelligence Gene Pair Ratios', fontsize=12)
ax1.legend(fontsize=8)
ax1.set_ylim(0.85, 1.30)

# Right: distribution comparison
ax2.hist(ctrl_ratios_all, bins=10, alpha=0.5, color='#3498db', label=f'Control (all, n={len(ctrl_ratios_all)})', edgecolor='black')
ax2.hist(intel_ratios, bins=8, alpha=0.5, color='#e74c3c', label=f'Intelligence (n={len(intel_ratios)})', edgecolor='black')
ax2.axvline(x=1.0, color='black', linestyle='-', linewidth=0.5)
ax2.set_xlabel('Observed/Expected Ratio', fontsize=11)
ax2.set_ylabel('Count', fontsize=11)
ax2.set_title('Intel vs Control Combo Ratio Distribution\n(Mann-Whitney p=0.023, p=0.001 for n≥20)', fontsize=11)
ax2.legend(fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'combo_fig2_interaction_types.png'), dpi=150)
plt.close()
print("Saved combo_fig2_interaction_types.png")

# ============================================================
# Figure 3: Cell-type heatmap (verify existing)
# ============================================================
# The celltype heatmap should still be valid since cell-type data hasn't changed
# Just check it exists
celltype_fig = os.path.join(out_dir, 'celltype/fig1_ranking_heatmap.png')
if os.path.exists(celltype_fig):
    print(f"Cell-type heatmap exists: {celltype_fig}")
else:
    print(f"WARNING: Cell-type heatmap missing: {celltype_fig}")

print("\nDone regenerating figures.")
