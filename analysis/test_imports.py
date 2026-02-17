#!/usr/bin/env python3
exec(open("/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py").read())

import sys
print("Testing imports...")

try:
    import scanpy as sc
    print("✓ scanpy")
except Exception as e:
    print(f"✗ scanpy: {e}")

try:
    import torch
    print("✓ torch")
except Exception as e:
    print(f"✗ torch: {e}")

try:
    import transformers
    print("✓ transformers")
except Exception as e:
    print(f"✗ transformers: {e}")

try:
    import pandas as pd
    print("✓ pandas")
except Exception as e:
    print(f"✗ pandas: {e}")

try:
    import numpy as np
    print("✓ numpy")
except Exception as e:
    print(f"✗ numpy: {e}")

print("Import testing complete")