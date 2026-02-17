#!/usr/bin/env python3
"""
Simple test to check if we can load the Geneformer model.
"""

import torch
import os

print("Testing model loading...")
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"CUDA version: {torch.version.cuda}")
    print(f"GPU: {torch.cuda.get_device_name(0)}")

# Try loading transformers
try:
    from transformers import BertForMaskedLM, BertConfig
    print("[OK] Transformers imported successfully")
except Exception as e:
    print(f"[ERROR] Failed to import transformers: {e}")
    exit(1)

# Check if model path exists
model_path = "D:/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V2-104M"
if os.path.exists(model_path):
    print(f"[OK] Model path exists: {model_path}")
    
    # List model files
    files = os.listdir(model_path)
    print(f"Model files: {files}")
    
    # Try loading config first
    try:
        config = BertConfig.from_pretrained(model_path)
        print(f"[OK] Model config loaded: vocab_size={config.vocab_size}")
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}")
        
    # Try loading just the config without the model
    try:
        print("Attempting to load model...")
        model = BertForMaskedLM.from_pretrained(model_path, config=config)
        print(f"[OK] Model loaded successfully! Type: {type(model)}")
        
        # Check model size
        total_params = sum(p.numel() for p in model.parameters())
        print(f"Total parameters: {total_params:,}")
        
    except Exception as e:
        print(f"[ERROR] Failed to load model: {e}")
        import traceback
        traceback.print_exc()
        
else:
    print(f"[ERROR] Model path does not exist: {model_path}")

print("Test completed.")