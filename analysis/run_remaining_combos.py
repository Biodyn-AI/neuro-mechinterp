#!/usr/bin/env python3
"""Run remaining combo perturbations one at a time via subprocess."""
import subprocess, os, sys

base = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
all_pairs = [
    ("CADM2", "APP"), ("GRIN2A", "CAMK2A"), ("GRIN2A", "MEF2C"),
    ("GRIN2A", "APP"), ("CAMK2A", "MEF2C"), ("CAMK2A", "APP"), ("MEF2C", "APP"),
]

for g1, g2 in all_pairs:
    pair_dir = os.path.join(base, f"perturb_combo_{g1}_{g2}")
    if os.path.exists(pair_dir):
        pickles = [f for f in os.listdir(pair_dir) if f.endswith('.pickle')]
        if pickles:
            print(f"SKIP {g1}+{g2}: already done", flush=True)
            continue
    
    print(f"\n>>> Running {g1} + {g2} <<<", flush=True)
    proc = subprocess.Popen(
        ['python3', '/mnt/d/openclaw/intelligence-augmentation/analysis/run_combo_single.py', g1, g2],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    for line in proc.stdout:
        line = line.rstrip()
        if 'RESULT' in line or 'DONE' in line or 'Perturbing' in line or 'CUDA' in line or 'Output' in line:
            print(line, flush=True)
    proc.wait()
    print(f">>> Finished {g1}+{g2}, exit={proc.returncode} <<<", flush=True)

print("\nALL COMBINATIONS COMPLETE", flush=True)
