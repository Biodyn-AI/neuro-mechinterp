import os, pickle, json
import numpy as np

base = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
results = {}
for d in sorted(os.listdir(base)):
    if d.startswith("perturb_"):
        gene = d.replace("perturb_", "")
        path = os.path.join(base, d)
        files = os.listdir(path)
        if files:
            pkl_file = os.path.join(path, files[0])
            with open(pkl_file, "rb") as f:
                data = pickle.load(f)
            # data is a dict: {token: {"original": emb_array, "perturbed": emb_array}}
            if isinstance(data, dict):
                for key, val in data.items():
                    if isinstance(val, dict) and "original" in val:
                        orig = np.array(val["original"])
                        pert = np.array(val["perturbed"])
                        shift = np.mean(np.linalg.norm(orig - pert, axis=1))
                        results[gene] = {"shift": float(shift), "n_cells": len(orig), "token": str(key)}
                    else:
                        results[gene] = {"raw_keys": list(val.keys())[:5] if isinstance(val, dict) else str(type(val))}
                    break
            else:
                results[gene] = {"type": str(type(data))}

# Sort by shift
ranked = sorted(results.items(), key=lambda x: x[1].get("shift", 0), reverse=True)
print("=== PERTURBATION IMPACT RANKING ===")
for gene, info in ranked:
    if "shift" in info:
        print(f"  {gene}: shift={info['shift']:.4f}, n_cells={info['n_cells']}, token={info['token']}")
    else:
        print(f"  {gene}: {info}")

with open(os.path.join(base, "perturbation_rankings.json"), "w") as f:
    json.dump(dict(ranked), f, indent=2)
print(f"\nSaved to perturbation_rankings.json")
