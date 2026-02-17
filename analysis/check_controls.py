#!/usr/bin/env python3
import pickle, numpy as np, os, glob

results_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
controls = {}
for d in sorted(glob.glob(os.path.join(results_dir, 'control_perturb_*'))):
    gene = os.path.basename(d).replace('control_perturb_', '')
    pkls = glob.glob(os.path.join(d, '*.pickle'))
    if not pkls: continue
    with open(pkls[0], 'rb') as f:
        data = pickle.load(f)
    for k, v in data.items():
        if isinstance(v, dict):
            for k2, v2 in v.items():
                arr = np.array(v2)
                shift = 1.0 - np.mean(arr)
                print(f'{gene}: n={len(arr)}, shift={shift:.4f}')
                controls[gene] = (len(arr), shift)
                break
