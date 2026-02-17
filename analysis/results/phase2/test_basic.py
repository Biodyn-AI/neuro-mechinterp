#!/usr/bin/env python3
"""Quick test to check if basic dependencies work."""

try:
    print("Testing basic imports...")
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import torch
    print("Basic imports successful!")
    
    print(f"NumPy version: {np.__version__}")
    print(f"Pandas version: {pd.__version__}")
    print(f"Scanpy version: {sc.__version__}")
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name()}")
        print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    print("Testing data loading...")
    from pathlib import Path
    data_file = Path(r'D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad')
    print(f"Data file exists: {data_file.exists()}")
    print(f"Data file size: {data_file.stat().st_size / 1e6:.1f} MB")
    
    print("All basic tests passed!")
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()