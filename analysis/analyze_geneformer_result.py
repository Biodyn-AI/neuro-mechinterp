#!/usr/bin/env python3
"""Analyze Geneformer perturbation results with proper handling of the data structure"""
import pickle
import numpy as np
import json
from pathlib import Path

def analyze_geneformer_result(result_file_path, gene_name="Unknown"):
    """Analyze a Geneformer perturbation result file"""
    print(f"Analyzing {gene_name} results from {result_file_path}")
    
    if not Path(result_file_path).exists():
        print(f"File not found: {result_file_path}")
        return None
    
    try:
        # Load the pickle file
        with open(result_file_path, 'rb') as f:
            data = pickle.load(f)
        
        print(f"Data type: {type(data)}")
        
        if hasattr(data, 'keys'):
            keys = list(data.keys())
            print(f"Keys: {keys[:5]}...")  # Show first 5 keys
            
            # Geneformer typically stores results with (gene_id, emb_type) tuples as keys
            embeddings = {}
            
            for key in keys:
                if isinstance(key, tuple) and len(key) == 2:
                    gene_id, emb_type = key
                    if emb_type in embeddings:
                        embeddings[emb_type].append(data[key])
                    else:
                        embeddings[emb_type] = [data[key]]
                    
                    if len(embeddings) == 1:  # Just show first example
                        print(f"Example key: {key} -> {type(data[key])}")
                        if hasattr(data[key], 'shape'):
                            print(f"  Shape: {data[key].shape}")
            
            print(f"Embedding types found: {list(embeddings.keys())}")
            
            # For single-gene perturbation, we typically get cell embeddings
            if 'cell_emb' in embeddings:
                cell_embeddings = embeddings['cell_emb']
                print(f"Found {len(cell_embeddings)} cell embeddings")
                
                # The challenge is we need baseline embeddings to compare against
                # For now, let's just analyze the distribution of embeddings
                if cell_embeddings:
                    emb_array = np.array(cell_embeddings)
                    print(f"Embeddings array shape: {emb_array.shape}")
                    
                    # Compute some basic statistics
                    emb_norms = np.linalg.norm(emb_array, axis=1) if len(emb_array.shape) > 1 else [np.linalg.norm(emb_array)]
                    
                    result = {
                        'gene': gene_name,
                        'n_cells': len(cell_embeddings),
                        'embedding_dim': emb_array.shape[-1] if len(emb_array.shape) > 1 else len(emb_array),
                        'mean_norm': float(np.mean(emb_norms)),
                        'std_norm': float(np.std(emb_norms)),
                        'min_norm': float(np.min(emb_norms)),
                        'max_norm': float(np.max(emb_norms))
                    }
                    
                    print(f"\n=== EMBEDDING ANALYSIS FOR {gene_name} ===")
                    print(f"Cells: {result['n_cells']}")
                    print(f"Embedding dimension: {result['embedding_dim']}")
                    print(f"Mean norm: {result['mean_norm']:.4f} ± {result['std_norm']:.4f}")
                    print(f"Norm range: [{result['min_norm']:.4f}, {result['max_norm']:.4f}]")
                    
                    # Note: Without baseline embeddings, we can't compute perturbation effect
                    print(f"Note: Need baseline embeddings to compute perturbation effect")
                    
                    return result
                    
        else:
            print(f"Data structure not recognized: {type(data)}")
    
    except Exception as e:
        print(f"Error analyzing file: {e}")
        import traceback
        traceback.print_exc()
    
    return None

def check_geneformer_output_structure(result_file):
    """Inspect the structure of a Geneformer output file"""
    print(f"\n=== DETAILED STRUCTURE INSPECTION ===")
    
    try:
        with open(result_file, 'rb') as f:
            data = pickle.load(f)
        
        print(f"Root type: {type(data)}")
        
        if hasattr(data, 'keys'):
            keys = list(data.keys())
            print(f"Number of keys: {len(keys)}")
            
            # Sample a few keys to understand structure
            for i, key in enumerate(keys[:3]):
                print(f"\nKey {i}: {key}")
                print(f"  Type: {type(key)}")
                print(f"  Value type: {type(data[key])}")
                if hasattr(data[key], 'shape'):
                    print(f"  Value shape: {data[key].shape}")
                elif isinstance(data[key], (list, tuple)):
                    print(f"  Value length: {len(data[key])}")
                    if len(data[key]) > 0:
                        print(f"  First element type: {type(data[key][0])}")
        
        return data
        
    except Exception as e:
        print(f"Error inspecting structure: {e}")
        return None

def main():
    """Analyze available Geneformer result files"""
    print("Analyzing Geneformer perturbation results")
    print("="*50)
    
    results_dir = Path("results/scaled")
    
    # Look for perturbation result directories
    perturbation_dirs = [d for d in results_dir.glob("perturb_*") if d.is_dir()]
    
    print(f"Found {len(perturbation_dirs)} perturbation directories")
    
    for perturb_dir in perturbation_dirs:
        # Extract gene name and cell count from directory name
        dir_parts = perturb_dir.name.split('_')
        if len(dir_parts) >= 3:
            gene_name = dir_parts[1]
            cell_count = dir_parts[2]
            
            print(f"\n{'='*60}")
            print(f"Processing {gene_name} ({cell_count} cells)")
            print(f"Directory: {perturb_dir}")
            print(f"{'='*60}")
            
            # Find pickle files
            pickle_files = list(perturb_dir.glob("*.pickle"))
            
            if pickle_files:
                result_file = pickle_files[0]
                print(f"Analyzing file: {result_file.name}")
                print(f"File size: {result_file.stat().st_size / 1024:.1f} KB")
                
                # First, inspect the structure
                data = check_geneformer_output_structure(result_file)
                
                # Then analyze the embeddings
                result = analyze_geneformer_result(result_file, gene_name)
                
            else:
                print(f"No pickle files found in {perturb_dir}")
    
    # Additional analysis: Compare to baseline if available
    baseline_file = Path("results/baseline_500_summary.json")
    if baseline_file.exists():
        print(f"\n{'='*60}")
        print("COMPARISON WITH 500-CELL BASELINE")
        print(f"{'='*60}")
        
        with open(baseline_file, 'r') as f:
            baseline = json.load(f)
        
        print("Baseline results (500 cells):")
        for gene, data in baseline.items():
            print(f"  {gene}: shift = {data['shift_mean']:.4f}, n_cells = {data['n_cells']}")
        
        print("\nNote: Full comparison requires baseline embeddings for scaled datasets")

if __name__ == "__main__":
    main()