import scanpy as sc
adata = sc.read_h5ad('/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/prepared_data/brain_500.h5ad')
print('Shape:', adata.shape)
print('Obs columns:', list(adata.obs.columns))
print()
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'category' or adata.obs[col].dtype == 'object' or adata.obs[col].nunique() < 20:
        print(f'--- {col}: {adata.obs[col].nunique()} unique ---')
        print(adata.obs[col].value_counts().head(10))
        print()
