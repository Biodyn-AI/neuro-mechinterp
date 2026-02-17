import scanpy as sc
adata = sc.read_h5ad('/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad')
print('var columns:', list(adata.var.columns[:10]))
print('var index name:', adata.var.index.name)
print('var index[:5]:', list(adata.var.index[:5]))
print('has ensembl_id:', 'ensembl_id' in adata.var.columns)
