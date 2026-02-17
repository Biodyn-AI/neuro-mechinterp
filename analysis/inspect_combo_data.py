#!/usr/bin/env python3
"""
Inspect the raw combinatorial perturbation pickle files to understand data structure.
"""

import pickle
import numpy as np
import os
import glob

RESULTS_DIR = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"

# Find all combo pickle files
combo_dirs = glob.glob(os.path.join(RESULTS_DIR, "perturb_combo_*"))

print(f"Found {len(combo_dirs)} combination directories")
print("=" * 60)

for i, combo_dir in enumerate(sorted(combo_dirs)[:3]):  # Examine first 3
    combo_name = os.path.basename(combo_dir).replace("perturb_combo_", "")
    print(f"\n{i+1}. {combo_name}")
    print("-" * 40)
    
    # Find pickle file
    pickle_files = glob.glob(os.path.join(combo_dir, "*.pickle"))
    if not pickle_files:
        print("  No pickle files found")
        continue
    
    pickle_file = pickle_files[0]
    print(f"  Pickle file: {os.path.basename(pickle_file)}")
    
    try:
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
        
        print(f"  Data type: {type(data)}")
        
        if isinstance(data, dict):
            print(f"  Keys: {list(data.keys())}")
            
            for key, value in data.items():
                if hasattr(value, 'shape'):
                    print(f"    {key}: shape {value.shape}, dtype {value.dtype}")
                    print(f"      sample values: {value.flat[:5]}")
                else:
                    print(f"    {key}: {type(value)} - {str(value)[:100]}...")
                    
        elif hasattr(data, 'shape'):
            print(f"  Shape: {data.shape}")
            print(f"  Dtype: {data.dtype}")
            print(f"  Sample values: {data.flat[:10]}")
        else:
            print(f"  Content: {str(data)[:200]}...")
            
    except Exception as e:
        print(f"  Error loading: {e}")

print("\n" + "=" * 60)
print("SAMPLE ANALYSIS COMPLETE")