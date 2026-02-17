#!/usr/bin/env python3
"""Test if 5000 cells is feasible (memory check only, no tokenization)"""
import scanpy as sc
import numpy as np
import psutil
import gc

def check_memory():
    """Check current memory usage"""
    memory = psutil.virtual_memory()
    print(f"Total RAM: {memory.total / 1e9:.1f} GB")
    print(f"Available RAM: {memory.available / 1e9:.1f} GB")
    print(f"Used RAM: {memory.used / 1e9:.1f} GB ({memory.percent:.1f}%)")

def test_cell_count_feasibility(n_cells):
    """Test if a given cell count is memory-feasible"""
    print(f"\n=== Testing {n_cells:,} cells ===")
    
    check_memory()
    
    try:
        print("Loading full dataset...")
        adata = sc.read_h5ad('/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad')
        print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        
        check_memory()
        
        if adata.n_obs < n_cells:
            print(f"WARNING: Dataset only has {adata.n_obs:,} cells, cannot sample {n_cells:,}")
            return False
        
        # Test subsample
        print(f"Subsampling to {n_cells:,} cells...")
        sc.pp.subsample(adata, n_obs=n_cells, random_state=42)
        
        # Estimate memory usage
        data_memory = adata.X.nbytes if hasattr(adata.X, 'nbytes') else 0
        print(f"Subsampled to: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        print(f"Data memory: {data_memory / 1e6:.1f} MB")
        
        check_memory()
        
        # Test basic preprocessing
        print("Testing basic preprocessing...")
        adata.var['ensembl_id'] = adata.var.index
        
        if 'n_counts' not in adata.obs.columns:
            adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        
        if 'cell_type' not in adata.obs.columns:
            if 'Subset' in adata.obs.columns:
                adata.obs['cell_type'] = adata.obs['Subset']
            else:
                adata.obs['cell_type'] = 'unknown'
        
        print(f"✓ {n_cells:,} cells appears feasible for preprocessing")
        
        del adata
        gc.collect()
        
        check_memory()
        return True
        
    except MemoryError:
        print(f"✗ MEMORY ERROR: {n_cells:,} cells too large")
        return False
    except Exception as e:
        print(f"✗ ERROR: {e}")
        return False

def main():
    """Test feasibility of different cell counts"""
    print("Testing cell count feasibility for larger dataset validation")
    print("=" * 60)
    
    check_memory()
    
    # Test increasing cell counts
    cell_counts = [1000, 2000, 5000, 10000]
    feasible_counts = []
    
    for n_cells in cell_counts:
        if test_cell_count_feasibility(n_cells):
            feasible_counts.append(n_cells)
        else:
            print(f"Stopping at {n_cells:,} cells due to memory constraints")
            break
        
        print("-" * 40)
    
    print(f"\n=== SUMMARY ===")
    print(f"Feasible cell counts: {feasible_counts}")
    print(f"Maximum feasible: {max(feasible_counts) if feasible_counts else 'None'}")
    
    if 5000 in feasible_counts:
        print("✓ 5000 cells appears feasible - can attempt full pipeline")
    elif 2000 in feasible_counts:
        print("⚠ Limited to 2000 cells - still good for validation")
    else:
        print("⚠ May need to stick with 1000 cells")

if __name__ == "__main__":
    main()