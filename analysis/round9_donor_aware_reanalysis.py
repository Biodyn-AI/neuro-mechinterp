import re
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

BASE = '/mnt/d/openclaw/intelligence-augmentation'
combined_path = f'{BASE}/analysis/results/final/combined_perturbation_results.csv'
adata_path = f'{BASE}/analysis/results/insilico_wsl/prepared_data/brain_500.h5ad'
out_json = f'{BASE}/analysis/round9_donor_aware_results.json'
out_csv = f'{BASE}/analysis/round9_donor_gene_stats.csv'
out_md = f'{BASE}/analysis/ROUND9_DONOR_AWARE_REPORT.md'

# Load
cdf = pd.read_csv(combined_path)
adata = ad.read_h5ad(adata_path)

# Keep rows where target gene present in cell (actual perturbation context)
cdf = cdf[cdf['target_in_cell'] == True].copy()

# Parse null list strings like [np.float32(0.00123), ...]
float_re = re.compile(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?')
def parse_null_mean(x):
    if pd.isna(x):
        return np.nan
    vals = [float(v) for v in float_re.findall(str(x))]
    return float(np.mean(vals)) if vals else np.nan

cdf['null_mean_cell'] = cdf['null_distances'].map(parse_null_mean)
cdf['delta_vs_null'] = cdf['cosine_distance'] - cdf['null_mean_cell']

# Map donor_id from cell_idx -> adata.obs
obs = adata.obs.reset_index(drop=True)
if 'donor_id' not in obs.columns:
    raise RuntimeError('donor_id column not found in adata.obs')

max_idx = obs.shape[0] - 1
cdf = cdf[(cdf['cell_idx'] >= 0) & (cdf['cell_idx'] <= max_idx)].copy()
cdf['donor_id'] = cdf['cell_idx'].astype(int).map(obs['donor_id']).astype(str)

# donor pseudobulk per gene
agg = (
    cdf.groupby(['target_gene', 'donor_id'], as_index=False, observed=True)
      .agg(
          n_cells=('delta_vs_null', 'size'),
          mean_delta=('delta_vs_null', 'mean'),
          mean_shift=('cosine_distance', 'mean'),
          mean_null=('null_mean_cell', 'mean')
      )
)

rows = []
for gene, g in agg.groupby('target_gene'):
    vals = g['mean_delta'].dropna().values
    n_donors = len(vals)
    if n_donors >= 3:
        t_stat, t_p = stats.ttest_1samp(vals, popmean=0.0)
        try:
            w_stat, w_p = stats.wilcoxon(vals)
        except Exception:
            w_stat, w_p = np.nan, np.nan
    else:
        t_stat = t_p = w_stat = w_p = np.nan

    rows.append({
        'gene': gene,
        'n_donors': n_donors,
        'n_cells_total': int(g['n_cells'].sum()),
        'donor_mean_delta': float(np.mean(vals)) if n_donors else np.nan,
        'donor_sd_delta': float(np.std(vals, ddof=1)) if n_donors > 1 else np.nan,
        't_p': float(t_p) if pd.notna(t_p) else np.nan,
        'wilcoxon_p': float(w_p) if pd.notna(w_p) else np.nan,
        'donor_min': float(np.min(vals)) if n_donors else np.nan,
        'donor_max': float(np.max(vals)) if n_donors else np.nan,
    })

res = pd.DataFrame(rows).sort_values('donor_mean_delta', ascending=False)

# BH FDR on t_p where available
p = res['t_p'].values.copy()
mask = np.isfinite(p)
q = np.full_like(p, np.nan, dtype=float)
if mask.sum() > 0:
    pm = p[mask]
    order = np.argsort(pm)
    ranked = pm[order]
    m = len(ranked)
    q_ranked = ranked * m / (np.arange(1, m + 1))
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_vals = np.empty_like(q_ranked)
    q_vals[order] = np.clip(q_ranked, 0, 1)
    q[mask] = q_vals
res['q_t_bh'] = q

# Summary
n_genes = int(res.shape[0])
nominal = int((res['t_p'] < 0.05).sum())
fdr10 = int((res['q_t_bh'] < 0.10).sum())
fdr05 = int((res['q_t_bh'] < 0.05).sum())

top = res[['gene','donor_mean_delta','t_p','q_t_bh','n_donors','n_cells_total']].head(10)

result_obj = {
    'n_rows_input': int(cdf.shape[0]),
    'n_genes': n_genes,
    'n_donors_total': int(cdf['donor_id'].nunique()),
    'nominal_t_p_lt_0_05': nominal,
    'fdr_q_lt_0_10': fdr10,
    'fdr_q_lt_0_05': fdr05,
    'top_genes_by_donor_mean_delta': top.to_dict(orient='records')
}

res.to_csv(out_csv, index=False)
with open(out_json, 'w', encoding='utf-8') as f:
    json.dump(result_obj, f, indent=2)

with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Round 9 Donor-aware Reanalysis (Pseudobulk over donors)\n\n')
    f.write('Unit-of-inference update: donor-level aggregation of per-cell matched-null deltas.\n\n')
    f.write(f"- Donors represented: **{result_obj['n_donors_total']}**\n")
    f.write(f"- Genes tested: **{n_genes}**\n")
    f.write(f"- Nominal (t-test p<0.05): **{nominal}/{n_genes}**\n")
    f.write(f"- FDR q<0.10: **{fdr10}/{n_genes}**\n")
    f.write(f"- FDR q<0.05: **{fdr05}/{n_genes}**\n\n")
    f.write('## Top genes by donor-mean delta\n\n')
    f.write(top.to_csv(sep='|', index=False))
    f.write('\n\n## Interpretation\n\n')
    f.write('- This directly addresses pseudoreplication by moving inference from cells to donors.\n')
    f.write('- If donor-level FDR significance remains absent/weak, claims should stay exploratory/set-level.\n')

print(json.dumps(result_obj, indent=2))
