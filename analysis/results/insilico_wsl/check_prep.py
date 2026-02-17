import scanpy as sc
a = sc.read_h5ad("/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/prepared_data/brain_500.h5ad")
print("n_counts:", "n_counts" in a.obs.columns)
print("ensembl_id:", "ensembl_id" in a.var.columns)
print("obs cols:", list(a.obs.columns[:5]))
print("shape:", a.shape)
