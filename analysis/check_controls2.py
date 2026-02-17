#!/usr/bin/env python3
import pickle, numpy as np, os, glob

results_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
# Check one control pickle structure
d = glob.glob(os.path.join(results_dir, 'control_perturb_DOCK9'))[0]
pkls = glob.glob(os.path.join(d, '*.pickle'))
print(f"Pickle files: {pkls}")
with open(pkls[0], 'rb') as f:
    data = pickle.load(f)
print(f"Type: {type(data)}")
if isinstance(data, dict):
    for k, v in data.items():
        print(f"  key={k}, type={type(v)}, ", end="")
        if isinstance(v, (list, np.ndarray)):
            arr = np.array(v)
            print(f"shape={arr.shape}, mean={np.mean(arr):.4f}")
        elif isinstance(v, dict):
            print(f"dict with {len(v)} keys, first keys: {list(v.keys())[:3]}")
            for k2, v2 in list(v.items())[:1]:
                arr = np.array(v2)
                print(f"    subkey={k2}, shape={arr.shape}, mean={np.mean(arr):.4f}")
        else:
            print(f"val={v}")
