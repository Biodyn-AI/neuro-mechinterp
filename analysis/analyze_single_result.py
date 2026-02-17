#!/usr/bin/env python3
"""Manually analyze a single perturbation result"""
import pickle
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import json
from pathlib import Path

def analyze_perturbation_result(result_file_path, gene_name="Unknown"):
    """Analyze a single perturbation result file"""
    print(f"Analyzing {gene_name} results from {result_file_path}")
    
    if not Path(result_file_path).exists():
        print(f"File not found: {result_file_path}")
        return None
    
    try:
        # Load the pickle file
        with open(result_file_path, 'rb') as f:
            data = pickle.load(f)
        
        print(f"Data type: {type(data)}")
        
        if isinstance(data, dict):
            print(f"Dictionary keys: {list(data.keys())}")
            
            # Look for embeddings
            embeddings_key = None
            for key in data.keys():
                if 'embeddings' in key.lower():
                    embeddings_key = key
                    break
            
            if embeddings_key:
                embeddings = data[embeddings_key]
                print(f"Found embeddings under key: {embeddings_key}")
                print(f"Embeddings type: {type(embeddings)}")
                
                if isinstance(embeddings, dict):
                    print(f"Embeddings dict keys: {list(embeddings.keys())}")
                    
                    # Look for perturbed and control embeddings
                    perturbed_embs = None
                    control_embs = None
                    
                    for key, value in embeddings.items():
                        if 'perturb' in key.lower():
                            perturbed_embs = value
                            print(f"Perturbed embeddings: {key} -> {type(value)}")
                        elif 'control' in key.lower() or 'baseline' in key.lower():
                            control_embs = value
                            print(f"Control embeddings: {key} -> {type(value)}")
                    
                    if perturbed_embs is not None and control_embs is not None:
                        return compute_cosine_similarities(perturbed_embs, control_embs, gene_name)
                    else:
                        print("Could not find both perturbed and control embeddings")
                
                elif isinstance(embeddings, (list, np.ndarray)):
                    print(f"Direct embeddings array: shape {np.array(embeddings).shape}")
                    # This might be just the perturbed embeddings
                    # We need both perturbed and control to compute similarity
                    print("Need both perturbed and control embeddings for comparison")
                    
        else:
            print(f"Data is not a dictionary: {type(data)}")
            if hasattr(data, 'shape'):
                print(f"Shape: {data.shape}")
    
    except Exception as e:
        print(f"Error analyzing file: {e}")
        import traceback
        traceback.print_exc()
    
    return None

def compute_cosine_similarities(perturbed_embs, control_embs, gene_name):
    """Compute cosine similarities between perturbed and control embeddings"""
    try:
        perturbed_array = np.array(perturbed_embs)
        control_array = np.array(control_embs)
        
        print(f"Perturbed shape: {perturbed_array.shape}")
        print(f"Control shape: {control_array.shape}")
        
        # Ensure we have matching dimensions
        n_cells = min(len(perturbed_array), len(control_array))
        print(f"Using {n_cells} cells for comparison")
        
        cos_sims = []
        
        for i in range(n_cells):
            if len(perturbed_array[i].shape) > 0 and len(control_array[i].shape) > 0:
                sim = cosine_similarity([perturbed_array[i]], [control_array[i]])[0][0]
                cos_sims.append(sim)
        
        if cos_sims:
            cos_sims = np.array(cos_sims)
            
            result = {
                'gene': gene_name,
                'n_cells_processed': len(cos_sims),
                'cosine_sim_mean': float(np.mean(cos_sims)),
                'cosine_sim_std': float(np.std(cos_sims)),
                'shift_mean': float(1.0 - np.mean(cos_sims)),
                'shift_std': float(np.std(cos_sims)),
                'cosine_sim_min': float(np.min(cos_sims)),
                'cosine_sim_max': float(np.max(cos_sims))
            }
            
            print(f"\n=== RESULTS FOR {gene_name} ===")
            print(f"Cells processed: {result['n_cells_processed']}")
            print(f"Cosine similarity: {result['cosine_sim_mean']:.4f} ± {result['cosine_sim_std']:.4f}")
            print(f"Shift magnitude: {result['shift_mean']:.4f} ± {result['shift_std']:.4f}")
            print(f"Range: [{result['cosine_sim_min']:.4f}, {result['cosine_sim_max']:.4f}]")
            
            return result
        else:
            print("No valid similarities computed")
            return None
            
    except Exception as e:
        print(f"Error computing similarities: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Analyze available result files"""
    results_dir = Path("results/scaled")
    
    # Look for perturbation result directories
    perturbation_dirs = [d for d in results_dir.glob("perturb_*") if d.is_dir()]
    
    print(f"Found {len(perturbation_dirs)} perturbation directories")
    
    all_results = {}
    
    for perturb_dir in perturbation_dirs:
        # Extract gene name and cell count from directory name
        dir_parts = perturb_dir.name.split('_')
        if len(dir_parts) >= 3:
            gene_name = dir_parts[1]
            cell_count = dir_parts[2]
            
            print(f"\n{'='*50}")
            print(f"Processing {gene_name} ({cell_count} cells)")
            print(f"{'='*50}")
            
            # Find pickle files
            pickle_files = list(perturb_dir.glob("*.pickle"))
            
            if pickle_files:
                result_file = pickle_files[0]  # Use first pickle file
                result = analyze_perturbation_result(result_file, gene_name)
                
                if result:
                    key = f"{gene_name}_{cell_count}"
                    all_results[key] = result
            else:
                print(f"No pickle files found in {perturb_dir}")
    
    # Save results
    if all_results:
        results_file = results_dir / "manual_analysis_results.json"
        with open(results_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        print(f"\n{'='*50}")
        print("SUMMARY")
        print(f"{'='*50}")
        print(f"Successfully analyzed {len(all_results)} results")
        print(f"Results saved to: {results_file}")
        
        for key, result in all_results.items():
            print(f"{key}: shift = {result['shift_mean']:.4f} ± {result['shift_std']:.4f}")
    else:
        print("No results successfully analyzed")

if __name__ == "__main__":
    main()