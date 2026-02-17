#!/usr/bin/env python3
import sys
import os
import torch
sys.path.append("/mnt/d/openclaw/intelligence-augmentation/models/geneformer")

print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    torch.cuda.empty_cache()
    print("GPU cache cleared")