from datasets import load_from_disk
from collections import Counter
import pickle, os, numpy as np

base = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
ds = load_from_disk(os.path.join(base, 'brain.dataset'))

# Cell types in tokenized dataset
ct = ds['cell_type']
print("Cell types in tokenized dataset:")
for k, v in Counter(ct).most_common():
    print(f"  {k}: {v}")

# Check the pickle value - it has 472 floats for CADM2, which was expressed in 472 cells
# These must be the cosine similarities per cell
# But which cells? We need to figure out the mapping

# Let's check if the order matches the dataset order
# The pickle key is (token_id, 'cell_emb')
# The value is a list of floats (cosine similarities)
# The length matches N cells where gene was expressed

# Check which cells express CADM2 in the tokenized data
# CADM2 token = 13697
token_id = 13697
expressing_cells = []
for i, input_ids in enumerate(ds['input_ids']):
    if token_id in input_ids:
        expressing_cells.append(i)
print(f"\nCADM2 (token {token_id}) expressed in {len(expressing_cells)} cells")

# Load pickle to verify
pf = os.path.join(base, 'perturb_CADM2', 
    'in_silico_delete_perturb_CADM2_cell_embs_dict_[np.int16(13697)]_raw.pickle')
with open(pf, 'rb') as f:
    data = pickle.load(f)
vals = data[(np.int16(13697), 'cell_emb')]
print(f"Pickle has {len(vals)} values")
print(f"Match: {len(expressing_cells) == len(vals)}")
print(f"First 5 values: {vals[:5]}")
print(f"Mean cosine sim: {np.mean(vals):.6f}")
