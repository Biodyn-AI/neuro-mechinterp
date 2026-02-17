import os, pickle, json
import numpy as np

base = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
results = {}

for d in sorted(os.listdir(base)):
    if d.startswith("perturb_"):
        gene = d.replace("perturb_", "")
        path = os.path.join(base, d)
        files = os.listdir(path)
        if not files:
            continue
        pkl_file = os.path.join(path, files[0])
        with open(pkl_file, "rb") as f:
            data = pickle.load(f)
        
        # Inspect structure
        if isinstance(data, dict):
            for key, val in data.items():
                print(f"{gene} key={key}, val_type={type(val)}")
                if isinstance(val, list) and len(val) > 0:
                    print(f"  list len={len(val)}, first_type={type(val[0])}")
                    if isinstance(val[0], dict):
                        print(f"  first dict keys: {list(val[0].keys())[:5]}")
                    elif isinstance(val[0], np.ndarray):
                        print(f"  first array shape: {val[0].shape}")
                    elif isinstance(val[0], (list, tuple)):
                        print(f"  first sublist len: {len(val[0])}")
                elif isinstance(val, np.ndarray):
                    print(f"  array shape: {val.shape}")
                break
        break  # just first gene to understand structure
