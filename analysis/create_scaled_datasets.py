#!/usr/bin/env python3
"""Create properly tokenized datasets for scaling analysis"""
import scanpy as sc
import os
import sys
from pathlib import Path

def prepare_and_tokenize_dataset(n_cells=1000):
    """Prepare and tokenize dataset with specified number of cells"""
    print(f"Creating {n_cells}-cell dataset...")
    
    # Step 1: Load and prepare h5ad data
    print("Step 1: Loading and preparing data...")
    adata = sc.read_h5ad('/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad')
    print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    # Add ensembl_id column from index (index already contains ENSG IDs)
    adata.var['ensembl_id'] = adata.var.index
    
    # Add n_counts if missing
    if 'n_counts' not in adata.obs.columns:
        import numpy as np
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        print(f"Added n_counts, mean={adata.obs['n_counts'].mean():.0f}")
    
    # Ensure cell_type exists
    if 'cell_type' not in adata.obs.columns:
        if 'Subset' in adata.obs.columns:
            adata.obs['cell_type'] = adata.obs['Subset']
        elif 'subclass' in adata.obs.columns:
            adata.obs['cell_type'] = adata.obs['subclass']
        else:
            print("Available obs columns:", list(adata.obs.columns[:10]))
            adata.obs['cell_type'] = 'unknown'
    
    # Subsample to specified number of cells
    if adata.n_obs > n_cells:
        sc.pp.subsample(adata, n_obs=n_cells, random_state=42)
        print(f"Subsampled to {adata.n_obs:,} cells")
    
    # Create output directories
    prepared_dir = f'/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/prepared_{n_cells}'
    output_dir = f'/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled'
    
    os.makedirs(prepared_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Save prepared h5ad
    h5ad_path = os.path.join(prepared_dir, f'brain_{n_cells}.h5ad')
    adata.write_h5ad(h5ad_path)
    print(f"Saved prepared data: {h5ad_path}")
    
    # Step 2: Tokenize the data
    print("Step 2: Tokenizing data...")
    try:
        from geneformer import TranscriptomeTokenizer
        
        tk = TranscriptomeTokenizer(
            custom_attr_name_dict={"cell_type": "cell_type"},
            nproc=1,  # Use single process for memory efficiency
            model_version="V1",
        )
        
        # Tokenize data
        tk.tokenize_data(
            data_directory=prepared_dir,
            output_directory=output_dir,
            output_prefix=f"brain_{n_cells}",
            file_format="h5ad",
            use_generator=False
        )
        
        tokenized_path = os.path.join(output_dir, f"brain_{n_cells}.dataset")
        print(f"Tokenized dataset created: {tokenized_path}")
        
        return tokenized_path
        
    except Exception as e:
        print(f"Tokenization error: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Create scaled datasets"""
    print("Creating scaled datasets for larger dataset validation...")
    
    # Test different cell counts
    cell_counts = [1000, 2000]  # Start with conservative numbers
    
    successful_datasets = []
    
    for n_cells in cell_counts:
        try:
            tokenized_path = prepare_and_tokenize_dataset(n_cells)
            if tokenized_path and os.path.exists(tokenized_path):
                successful_datasets.append((n_cells, tokenized_path))
                print(f"✓ Successfully created {n_cells}-cell dataset")
            else:
                print(f"✗ Failed to create {n_cells}-cell dataset")
                
        except Exception as e:
            print(f"Error with {n_cells} cells: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\nSuccessfully created {len(successful_datasets)} datasets:")
    for n_cells, path in successful_datasets:
        print(f"  {n_cells:,} cells -> {path}")
    
    return successful_datasets

if __name__ == "__main__":
    main()