#!/usr/bin/env python3
"""Run all remaining combo perturbations sequentially as subprocesses."""
import subprocess, os, sys

base = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
pairs = [
    ("CADM2", "CAMK2A"), ("CADM2", "MEF2C"), ("CADM2", "APP"),
    ("GRIN2A", "CAMK2A"), ("GRIN2A", "MEF2C"), ("GRIN2A", "APP"),
    ("CAMK2A", "MEF2C"), ("CAMK2A", "APP"), ("MEF2C", "APP"),
]

for g1, g2 in pairs:
    pair_dir = os.path.join(base, f"perturb_combo_{g1}_{g2}")
    if os.path.exists(pair_dir):
        pickles = [f for f in os.listdir(pair_dir) if f.endswith('.pickle')]
        if pickles:
            print(f"SKIP {g1}+{g2}: already done", flush=True)
            continue
    
    print(f"\n=== Running {g1} + {g2} ===", flush=True)
    result = subprocess.run(
        ['python3', '/mnt/d/openclaw/intelligence-augmentation/analysis/run_combo_single.py', g1, g2],
        stdout=sys.stdout, stderr=sys.stdout, timeout=600
    )
    print(f"=== Done {g1}+{g2}, exit={result.returncode} ===", flush=True)

print("\nALL DONE", flush=True)
