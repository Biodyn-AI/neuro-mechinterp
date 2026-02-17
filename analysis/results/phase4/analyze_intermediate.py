"""Analyze Phase 4 intermediate results (50 cells) and produce final statistics"""
import pandas as pd
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

df = pd.read_csv(r'D:\openclaw\intelligence-augmentation\analysis\results\phase4\intermediate_results.csv')
print(f"Loaded {len(df)} perturbations from {df.cell_idx.nunique()} cells")
print(f"Cell types: {df.cell_type.unique()}")

# Gene-level statistics
gene_stats = []
for gene in df.gene.unique():
    gd = df[df.gene == gene]
    dists = gd.cosine_distance.dropna()
    pvals = gd.p_value.dropna()
    nulls = gd.null_mean.dropna()
    
    mean_dist = dists.mean()
    mean_null = nulls.mean() if len(nulls) > 0 else 0
    effect = mean_dist - mean_null
    
    gene_stats.append({
        'gene': gene,
        'n_cells': len(dists),
        'mean_cosine_distance': mean_dist,
        'std': dists.std(),
        'mean_null_distance': mean_null,
        'effect_size': effect,
        'effect_ratio': mean_dist / mean_null if mean_null > 0 else np.nan,
        'median_p': pvals.median() if len(pvals) > 0 else np.nan,
        'n_significant_005': (pvals < 0.05).sum() if len(pvals) > 0 else 0,
        'pct_significant': ((pvals < 0.05).sum() / len(pvals) * 100) if len(pvals) > 0 else 0
    })

gene_df = pd.DataFrame(gene_stats).sort_values('effect_size', ascending=False)

# FDR
valid_p = gene_df['median_p'].dropna()
if len(valid_p) > 0:
    _, fdr, _, _ = multipletests(valid_p, method='fdr_bh')
    gene_df.loc[valid_p.index, 'fdr'] = fdr

gene_df.to_csv(r'D:\openclaw\intelligence-augmentation\analysis\results\phase4\gene_statistics_phase4.csv', index=False)

print("\n=== GENE RANKINGS BY EFFECT SIZE ===")
cols = ['gene','n_cells','mean_cosine_distance','effect_size','effect_ratio','median_p','fdr','n_significant_005','pct_significant']
print(gene_df[cols].to_string(index=False))

# Cell type breakdown for top genes
print("\n=== CELL TYPE BREAKDOWN (top 5 genes) ===")
top5 = gene_df.head(5).gene.values
for gene in top5:
    gd = df[df.gene == gene]
    ct_stats = gd.groupby('cell_type').agg(
        n=('cosine_distance','count'),
        mean_dist=('cosine_distance','mean'),
        sig=('p_value', lambda x: (x<0.05).sum())
    ).sort_values('mean_dist', ascending=False)
    print(f"\n{gene}:")
    print(ct_stats.to_string())

print("\nDONE")
