import pickle, os
import numpy as np

# Check pickle structure
base = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
pf = os.path.join(base, 'perturb_CADM2')
files = os.listdir(pf)
print("CADM2 files:", files)

for f in files:
    fp = os.path.join(pf, f)
    with open(fp, 'rb') as fh:
        data = pickle.load(fh)
    print(f"\n{f}: type={type(data)}")
    if isinstance(data, dict):
        print(f"  keys: {list(data.keys())[:5]}")
        first_key = list(data.keys())[0]
        val = data[first_key]
        print(f"  first key: {first_key}, val type: {type(val)}")
        if isinstance(val, np.ndarray):
            print(f"  val shape: {val.shape}")
        elif isinstance(val, list):
            print(f"  val len: {len(val)}, first elem type: {type(val[0]) if val else None}")
            if val and isinstance(val[0], np.ndarray):
                print(f"  first elem shape: {val[0].shape}")

# Check tokenized dataset
from datasets import load_from_disk
ds = load_from_disk(os.path.join(base, 'brain.dataset'))
print(f"\nDataset: {ds}")
print(f"Columns: {ds.column_names}")
print(f"First example keys: {list(ds[0].keys())}")
# Check if cell_type is in dataset
if 'cell_type' in ds.column_names:
    print("cell_type present!")
else:
    print("No cell_type in dataset")
    print("Sample row:", {k: str(ds[0][k])[:100] for k in ds.column_names})
