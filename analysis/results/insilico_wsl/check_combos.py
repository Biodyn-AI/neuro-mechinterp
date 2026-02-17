import os, pickle, json
import numpy as np

base = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"

# Check all combo directories
combo_results = {}
for d in sorted(os.listdir(base)):
    if d.startswith("perturb_combo_"):
        combo_name = d.replace("perturb_combo_", "")
        path = os.path.join(base, d)
        files = [f for f in os.listdir(path) if f.endswith('.pickle')]
        if files:
            pkl_file = os.path.join(path, files[0])
            with open(pkl_file, "rb") as f:
                data = pickle.load(f)
            for key, val in data.items():
                if isinstance(val, list):
                    arr = np.array(val)
                    combo_results[combo_name] = {
                        "mean_cos_sim": float(np.mean(arr)),
                        "mean_shift": float(1 - np.mean(arr)),
                        "std": float(np.std(arr)),
                        "n_cells": len(arr),
                        "file": files[0]
                    }
            print(f"{combo_name}: {combo_results.get(combo_name, 'NO DATA')}")
        else:
            print(f"{combo_name}: NO PICKLE FILES")

# Also load single-gene results for synergy comparison
single_results = {}
for d in sorted(os.listdir(base)):
    if d.startswith("perturb_") and "combo" not in d:
        gene = d.replace("perturb_", "")
        path = os.path.join(base, d)
        files = [f for f in os.listdir(path) if f.endswith('.pickle')]
        if files:
            pkl_file = os.path.join(path, files[0])
            with open(pkl_file, "rb") as f:
                data = pickle.load(f)
            for key, val in data.items():
                if isinstance(val, list):
                    arr = np.array(val)
                    single_results[gene] = {
                        "mean_cos_sim": float(np.mean(arr)),
                        "mean_shift": float(1 - np.mean(arr)),
                        "std": float(np.std(arr)),
                        "n_cells": len(arr),
                    }

# Compute synergy scores
print("\n=== SYNERGY ANALYSIS ===")
print(f"{'Combo':<25} {'Combo Shift':>12} {'Sum Singles':>12} {'Synergy':>10} {'Type':>12}")
print("-" * 75)
for combo_name, cdata in sorted(combo_results.items()):
    genes = combo_name.split("_")
    if len(genes) == 2 and genes[0] in single_results and genes[1] in single_results:
        sum_singles = single_results[genes[0]]["mean_shift"] + single_results[genes[1]]["mean_shift"]
        synergy = cdata["mean_shift"] - sum_singles
        synergy_type = "SYNERGISTIC" if synergy > 0.001 else ("ANTAGONISTIC" if synergy < -0.001 else "ADDITIVE")
        print(f"{combo_name:<25} {cdata['mean_shift']:>12.6f} {sum_singles:>12.6f} {synergy:>10.6f} {synergy_type:>12}")

# Save
with open(os.path.join(base, "combo_synergy.json"), "w") as f:
    json.dump({"combo": combo_results, "singles": {k: v for k, v in single_results.items() if k in ["CADM2","GRIN2A","CAMK2A","MEF2C","APP"]}}, f, indent=2)
print("\nSaved combo_synergy.json")
