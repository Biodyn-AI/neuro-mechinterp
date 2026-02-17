#!/usr/bin/env python3
"""Test perturbation on 1000 cells with CADM2 (top gene)"""
import os
import sys
import torch
import time
from pathlib import Path

# Memory settings for larger dataset
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:512"

print("Testing 1000-cell perturbation with CADM2...")
print(f"CUDA available: {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# Setup paths
model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
data_path = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_1000_test.h5ad"
token_dict = "/mnt/d/openclaw/intelligence-augmentation/models/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/test_cadm2_1000"

print(f"Model: {model_path}")
print(f"Data: {data_path}")
print(f"Output: {output_dir}")

try:
    # Import Geneformer
    sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
    from geneformer import InSilicoPerturber
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("Initializing InSilicoPerturber...")
    start_time = time.time()
    
    # Initialize with conservative settings for larger dataset
    isp = InSilicoPerturber(
        perturb_type="delete",
        genes_to_perturb=["ENSG00000175161"],  # CADM2
        model_type="Pretrained", 
        num_classes=0,
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,
        max_ncells=1000,
        emb_layer=-1,
        forward_batch_size=8,  # Smaller batch for memory safety
        nproc=1,
        model_version="V1",
        token_dictionary_file=token_dict,
    )
    
    print(f"Initialization time: {time.time() - start_time:.1f}s")
    print("Running perturbation...")
    
    # Run perturbation
    perturb_start = time.time()
    isp.perturb_data(
        model_directory=model_path,
        input_data_file=data_path,
        output_directory=output_dir,
        output_prefix="cadm2_1000_test",
    )
    
    perturb_time = time.time() - perturb_start
    total_time = time.time() - start_time
    
    print(f"Perturbation time: {perturb_time:.1f}s")
    print(f"Total time: {total_time:.1f}s")
    
    # Check output files
    import glob
    output_files = glob.glob(f"{output_dir}/*.pickle")
    print(f"Output files: {len(output_files)}")
    for f in output_files:
        size_mb = os.path.getsize(f) / 1e6
        print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
    
    print("SUCCESS! 1000-cell perturbation completed.")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
    
    if "OutOfMemoryError" in str(e) or "CUDA out of memory" in str(e):
        print("\nOOM Error - try reducing batch size or cell count")
    elif "No such file" in str(e):
        print("\nFile not found - check paths")
    else:
        print("\nUnknown error - check logs")