#!/usr/bin/env python3
"""Debug combo data to understand the dimension mismatch."""

import os
import pickle
import numpy as np
import pandas as pd
import glob
from scipy.spatial.distance import cosine

RESULTS_DIR = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
BASELINE_CSV = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/control_1000/control_1000.csv"

def debug_baseline():
    """Check baseline data structure."""
    print("=== BASELINE DEBUG ===")
    baseline_df = pd.read_csv(BASELINE_CSV, index_col=0)
    baseline_embeddings = baseline_df.values
    
    print(f"Baseline shape: {baseline_embeddings.shape}")
    print(f"Baseline dtype: {baseline_embeddings.dtype}")
    print(f"Baseline sample values: {baseline_embeddings[0][:5]}")
    print(f"Baseline value range: [{baseline_embeddings.min():.4f}, {baseline_embeddings.max():.4f}]")
    print(f"Baseline mean magnitude: {np.mean(np.linalg.norm(baseline_embeddings, axis=1)):.4f}")
    
    return baseline_embeddings

def debug_combo_data():
    """Check combo data structure."""
    print("\n=== COMBO DATA DEBUG ===")
    
    combo_dirs = glob.glob(os.path.join(RESULTS_DIR, "perturb_combo_*"))
    
    for combo_dir in sorted(combo_dirs)[:2]:  # Check first 2
        combo_name = os.path.basename(combo_dir).replace("perturb_combo_", "")
        print(f"\n--- {combo_name} ---")
        
        pickle_files = glob.glob(os.path.join(combo_dir, "*.pickle"))
        if not pickle_files:
            print("  No pickle files")
            continue
            
        pickle_file = pickle_files[0]
        print(f"  File: {os.path.basename(pickle_file)}")
        
        try:
            with open(pickle_file, 'rb') as f:
                data = pickle.load(f)
            
            print(f"  Data type: {type(data)}")
            
            if isinstance(data, dict):
                print(f"  Keys: {list(data.keys())}")
                
                for key, value in data.items():
                    print(f"    Key: {key}")
                    print(f"    Value type: {type(value)}")
                    
                    if isinstance(value, list):
                        arr = np.array(value)
                        print(f"    As array shape: {arr.shape}")
                        print(f"    As array dtype: {arr.dtype}")
                        print(f"    Sample values: {arr[:3] if len(arr) > 0 else 'empty'}")
                        print(f"    Value range: [{arr.min():.4f}, {arr.max():.4f}]")
                        
                        # Try computing similarity with baseline
                        if len(arr.shape) == 1:
                            print(f"    ❌ 1D array - this is the problem!")
                            print(f"    Expected: 2D array (n_cells, n_features)")
                        elif len(arr.shape) == 2:
                            print(f"    ✅ 2D array - correct format")
                            print(f"    Features per cell: {arr.shape[1]}")
                        
                    elif hasattr(value, 'shape'):
                        print(f"    Shape: {value.shape}")
                        print(f"    Dtype: {value.dtype}")
                        print(f"    Sample: {value.flat[:5]}")
                        
        except Exception as e:
            print(f"  ERROR loading: {e}")
            import traceback
            traceback.print_exc()

def test_cosine_computation():
    """Test cosine similarity computation with sample data."""
    print("\n=== COSINE SIMILARITY TEST ===")
    
    # Create test vectors
    vec1 = np.random.randn(256)
    vec2 = np.random.randn(256)
    vec3 = vec1 + 0.1 * np.random.randn(256)  # Similar to vec1
    
    print(f"Vec1 shape: {vec1.shape}")
    print(f"Vec2 shape: {vec2.shape}")
    print(f"Vec3 shape: {vec3.shape}")
    
    try:
        sim12 = 1 - cosine(vec1, vec2)
        sim13 = 1 - cosine(vec1, vec3)
        
        print(f"Similarity vec1-vec2 (random): {sim12:.6f}")
        print(f"Similarity vec1-vec3 (similar): {sim13:.6f}")
        print(f"Shift vec1-vec2: {1 - sim12:.6f}")
        print(f"Shift vec1-vec3: {1 - sim13:.6f}")
        
        print("✅ Cosine computation works")
        
    except Exception as e:
        print(f"❌ Cosine computation failed: {e}")

def main():
    baseline = debug_baseline()
    debug_combo_data()
    test_cosine_computation()

if __name__ == "__main__":
    main()