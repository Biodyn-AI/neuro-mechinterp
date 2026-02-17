#!/usr/bin/env python3

import sys
import os

# Add scGPT to path
sys.path.append('/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT')

try:
    import scgpt
    print("✅ scGPT import successful!")
    print(f"scGPT version: {scgpt.__version__}")
    
    # Check if we can import key modules
    from scgpt.model import TransformerModel
    from scgpt.tokenizer import tokenize_and_pad_batch
    from scgpt.tasks import PerturbationTask
    print("✅ Key scGPT modules imported successfully!")
    
except ImportError as e:
    print(f"❌ Failed to import scGPT: {e}")
    sys.exit(1)