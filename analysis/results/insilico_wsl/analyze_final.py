import os, pickle, json
import numpy as np

base = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
results = {}

for d in sorted(os.listdir(base)):
    if d.startswith("perturb_"):
        gene = d.replace("perturb_", "")
        path = os.path.join(base, d)
        files = [f for f in os.listdir(path) if f.endswith('.pickle')]
        if not files:
            continue
        pkl_file = os.path.join(path, files[0])
        with open(pkl_file, "rb") as f:
            data = pickle.load(f)
        
        for key, val in data.items():
            token_id = key[0] if isinstance(key, tuple) else key
            emb_type = key[1] if isinstance(key, tuple) else "unknown"
            
            if isinstance(val, list):
                arr = np.array(val)
                # These are cosine similarity shifts (perturbed vs original embedding)
                mean_shift = float(np.mean(arr))
                std_shift = float(np.std(arr))
                n_cells = len(arr)
                results[gene] = {
                    "mean_shift": mean_shift,
                    "std_shift": std_shift, 
                    "n_cells": n_cells,
                    "token": int(token_id),
                    "emb_type": emb_type
                }

# Rank by absolute mean shift (larger = more impact)
ranked = sorted(results.items(), key=lambda x: abs(x[1]["mean_shift"]), reverse=True)

print("=== GENEFORMER InSilicoPerturber RESULTS ===")
print(f"{'Gene':<10} {'Mean Shift':>12} {'Std':>8} {'N Cells':>8}")
print("-" * 45)
for gene, info in ranked:
    print(f"{gene:<10} {info['mean_shift']:>12.6f} {info['std_shift']:>8.6f} {info['n_cells']:>8}")

print(f"\nTotal genes analyzed: {len(results)}/21")

# Save
with open(os.path.join(base, "final_rankings.json"), "w") as f:
    json.dump(dict(ranked), f, indent=2)
print("Saved to final_rankings.json")
