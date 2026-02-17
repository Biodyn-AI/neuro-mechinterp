#!/usr/bin/env python3
"""Verify Table 1 single-gene values against actual pickle data."""
import pickle, numpy as np, os, glob

results_dir = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'

# Expected from Table 1
table1 = {
    'CADM2': (0.0196, 472), 'GRIN2A': (0.0190, 223), 'CAMK2A': (0.0189, 155),
    'MEF2C': (0.0184, 216), 'APP': (0.0183, 430), 'SCN1A': (0.0179, 185),
    'NRXN1': (0.0178, 241), 'GRIN2B': (0.0176, 181), 'HOMER1': (0.0175, 167),
    'NEGR1': (0.0166, 278),
}

# Check against pickles
for gene_dir in sorted(glob.glob(os.path.join(results_dir, 'perturb_*'))):
    name = os.path.basename(gene_dir)
    if 'combo' in name or 'control' in name:
        continue
    gene = name.replace('perturb_', '')
    if gene not in table1:
        continue
    
    pkls = glob.glob(os.path.join(gene_dir, '*.pickle'))
    if not pkls:
        print(f"  {gene}: NO PICKLE")
        continue
    
    with open(pkls[0], 'rb') as f:
        data = pickle.load(f)
    
    for k, v in data.items():
        if isinstance(k, tuple) and 'cell_emb' in str(k):
            arr = np.array(v)
            shift = 1.0 - np.mean(arr)
            n = len(arr)
            expected_shift, expected_n = table1[gene]
            match_shift = abs(shift - expected_shift) < 0.001
            match_n = n == expected_n
            status = "OK" if (match_shift and match_n) else "MISMATCH"
            print(f"  {gene}: actual={shift:.4f} n={n} | table={expected_shift:.4f} n={expected_n} | {status}")
            break

# Check combo values
print("\n=== Combo Table Verification ===")
table3 = {
    'GRIN2A_CAMK2A': 0.0454, 'GRIN2A_MEF2C': 0.0452,
    'CADM2_GRIN2A': 0.0445, 'CADM2_CAMK2A': 0.0445,
    'CAMK2A_MEF2C': 0.0441, 'CADM2_MEF2C': 0.0430,
    'CADM2_APP': 0.0423, 'GRIN2A_APP': 0.0416,
    'CAMK2A_APP': 0.0416, 'MEF2C_APP': 0.0410,
}

for gene_dir in sorted(glob.glob(os.path.join(results_dir, 'perturb_combo_*'))):
    name = os.path.basename(gene_dir).replace('perturb_combo_', '')
    pkls = glob.glob(os.path.join(gene_dir, '*.pickle'))
    if not pkls:
        continue
    with open(pkls[0], 'rb') as f:
        data = pickle.load(f)
    for k, v in data.items():
        if isinstance(k, tuple) and 'cell_emb' in str(k):
            arr = np.array(v)
            shift = 1.0 - np.mean(arr)
            expected = table3.get(name, None)
            if expected:
                match = abs(shift - expected) < 0.001
                print(f"  {name}: actual={shift:.4f} | table={expected:.4f} | {'OK' if match else 'MISMATCH'}")
            break
