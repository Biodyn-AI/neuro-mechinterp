import scanpy as sc
import os

print("Loading data...")
adata = sc.read_h5ad('/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad')
print(f"Loaded: {adata.shape}")

# Add ensembl_id column from index (index already contains ENSG IDs)
adata.var['ensembl_id'] = adata.var.index

# Add n_counts if missing
if 'n_counts' not in adata.obs.columns:
    import numpy as np
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    print(f"Added n_counts, mean={adata.obs['n_counts'].mean():.0f}")

# Ensure cell_type exists
if 'cell_type' not in adata.obs.columns:
    if 'Subset' in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs['Subset']
    elif 'subclass' in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs['subclass']
    else:
        print("obs columns:", list(adata.obs.columns[:10]))
        adata.obs['cell_type'] = 'unknown'

# Subsample to 500 cells for speed
if adata.n_obs > 500:
    sc.pp.subsample(adata, n_obs=500, random_state=42)
    print(f"Subsampled to {adata.n_obs} cells")

# Save prepared data
out_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/prepared_data'
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, 'brain_500.h5ad')
adata.write_h5ad(out_path)
print(f"Saved to {out_path}")
print(f"ensembl_id column: {adata.var['ensembl_id'][:3].tolist()}")
print(f"cell_type values: {adata.obs['cell_type'].value_counts().head().to_dict()}")
