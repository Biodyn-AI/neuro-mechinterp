#!/usr/bin/env python3
"""Inspect the structure of 1000-cell perturbation pickle files"""
import pickle
import os
import numpy as np
from pathlib import Path

def inspect_pickle_file(pickle_path):
    """Inspect what's inside a pickle file"""
    print(f"🔍 Inspecting: {pickle_path}")
    
    try:
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        
        print(f"   Data type: {type(data)}")
        
        if hasattr(data, 'keys'):
            print(f"   Keys: {list(data.keys())}")
            for key in list(data.keys())[:5]:  # Show first 5 keys
                print(f"      Key: {key} -> {type(data[key])}")
                if hasattr(data[key], 'keys'):
                    print(f"         Sub-keys: {list(data[key].keys())}")
                elif hasattr(data[key], 'shape'):
                    print(f"         Shape: {data[key].shape}")
                elif hasattr(data[key], '__len__') and not isinstance(data[key], str):
                    print(f"         Length: {len(data[key])}")
        elif hasattr(data, 'shape'):
            print(f"   Shape: {data.shape}")
        elif hasattr(data, '__len__'):
            print(f"   Length: {len(data)}")
            
        return data
    except Exception as e:
        print(f"❌ Error loading {pickle_path}: {e}")
        return None

def main():
    # Check one pickle file from each gene
    base_path = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled")
    
    genes = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
    
    for gene in genes:
        gene_dir = base_path / f"perturb_{gene}_1000"
        if gene_dir.exists():
            pickle_files = list(gene_dir.glob("*.pickle"))
            if pickle_files:
                pickle_file = pickle_files[0]
                print(f"\n{'='*50}")
                print(f"Gene: {gene}")
                print('='*50)
                data = inspect_pickle_file(pickle_file)
                
                # Try to understand the structure better
                if data and hasattr(data, 'keys'):
                    for key in data.keys():
                        value = data[key]
                        print(f"   Detailed view of key '{key}':")
                        print(f"      Type: {type(value)}")
                        if hasattr(value, 'keys'):
                            for sub_key, sub_value in value.items():
                                print(f"         {sub_key}: {type(sub_value)}")
                                if hasattr(sub_value, 'shape'):
                                    print(f"            Shape: {sub_value.shape}")
                                    print(f"            Sample values: {sub_value[:3] if len(sub_value) > 3 else sub_value}")
                        elif isinstance(value, list) and len(value) > 0:
                            print(f"         List with {len(value)} elements")
                            print(f"         First element type: {type(value[0])}")
                            if hasattr(value[0], 'shape'):
                                print(f"         First element shape: {value[0].shape}")
                            if len(value) > 0 and hasattr(value[0], '__len__') and not isinstance(value[0], str):
                                print(f"         First element length: {len(value[0])}")
                                if hasattr(value[0], '__getitem__'):
                                    try:
                                        sample = value[0][:5] if len(value[0]) > 5 else value[0]
                                        print(f"         First element sample: {sample}")
                                    except:
                                        pass
                break
            else:
                print(f"❌ No pickle files found for {gene}")
        else:
            print(f"❌ Directory not found: {gene_dir}")

if __name__ == "__main__":
    main()