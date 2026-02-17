#!/usr/bin/env python3
"""Analyze 1000-cell perturbation results - handles Geneformer output format"""
import os
import json
import pickle
import numpy as np
from pathlib import Path
from sklearn.metrics.pairwise import cosine_similarity

# Gene mappings
GENES = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
GENE_ENSEMBL = {
    "CADM2": 13697,
    "GRIN2A": 15103,
    "CAMK2A": 1135,
    "MEF2C": 1532,
    "APP": 8172
}

RESULTS_DIR = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled")

def load_geneformer_result(gene_symbol):
    """Load Geneformer perturbation result"""
    result_dir = RESULTS_DIR / f"perturb_{gene_symbol}_1000"
    
    if not result_dir.exists():
        print(f"❌ Directory not found: {result_dir}")
        return None
    
    pickle_files = list(result_dir.glob("*.pickle"))
    if not pickle_files:
        print(f"❌ No pickle files found for {gene_symbol}")
        return None
    
    pickle_file = pickle_files[0]
    print(f"📁 Loading {gene_symbol}: {pickle_file.name}")
    
    try:
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
        
        print(f"   Data type: {type(data)}")
        print(f"   Keys: {list(data.keys())}")
        
        return data
        
    except Exception as e:
        print(f"❌ Error loading {gene_symbol}: {e}")
        return None

def extract_embeddings_from_geneformer_data(data, gene_symbol):
    """Extract embeddings from Geneformer defaultdict structure"""
    gene_token = GENE_ENSEMBL[gene_symbol]
    
    # Look for the key with our gene token
    embedding_key = None
    for key in data.keys():
        if isinstance(key, tuple) and len(key) == 2:
            token, emb_type = key
            # Handle numpy integer types
            token_val = int(token) if hasattr(token, 'item') else token
            if token_val == gene_token and emb_type == 'cell_emb':
                embedding_key = key
                break
    
    if embedding_key is None:
        print(f"❌ Could not find embedding key for {gene_symbol} (token {gene_token})")
        print(f"   Available keys: {list(data.keys())}")
        return None, None
    
    print(f"   Found embedding key: {embedding_key}")
    
    # Extract the embedding data
    emb_data = data[embedding_key]
    print(f"   Embedding data type: {type(emb_data)}")
    
    # The data structure might be a dict with multiple keys
    if isinstance(emb_data, dict):
        print(f"   Embedding dict keys: {list(emb_data.keys())}")
        
        # Look for standard keys
        perturbed_emb = None
        control_emb = None
        
        # Try different possible key patterns
        for key, value in emb_data.items():
            key_lower = str(key).lower()
            if 'perturb' in key_lower or 'delete' in key_lower:
                perturbed_emb = value
                print(f"   Found perturbed embeddings: '{key}'")
            elif 'control' in key_lower or 'baseline' in key_lower or 'original' in key_lower:
                control_emb = value
                print(f"   Found control embeddings: '{key}'")
        
        # If not found by name, check if we have exactly 2 arrays
        if perturbed_emb is None and control_emb is None:
            values = list(emb_data.values())
            if len(values) == 2:
                perturbed_emb, control_emb = values
                print(f"   Using first two values as perturbed and control")
        
        return perturbed_emb, control_emb
    
    else:
        print(f"❌ Unexpected embedding data structure: {type(emb_data)}")
        return None, None

def compute_perturbation_metrics(perturbed_emb, control_emb, gene_symbol):
    """Compute perturbation metrics from embeddings"""
    if perturbed_emb is None or control_emb is None:
        return None
    
    try:
        # Convert to numpy arrays
        perturbed_emb = np.array(perturbed_emb)
        control_emb = np.array(control_emb)
        
        print(f"   Perturbed shape: {perturbed_emb.shape}")
        print(f"   Control shape: {control_emb.shape}")
        
        # Handle different dimensionalities
        if perturbed_emb.ndim == 1:
            perturbed_emb = perturbed_emb.reshape(1, -1)
        if control_emb.ndim == 1:
            control_emb = control_emb.reshape(1, -1)
        
        n_cells = min(len(perturbed_emb), len(control_emb))
        cos_sims = []
        
        for i in range(n_cells):
            p_emb = perturbed_emb[i].reshape(1, -1)
            c_emb = control_emb[i].reshape(1, -1)
            
            sim = cosine_similarity(p_emb, c_emb)[0][0]
            cos_sims.append(sim)
        
        cos_sims = np.array(cos_sims)
        
        # Compute statistics
        result = {
            'n_cells_processed': n_cells,
            'cosine_sim_mean': float(np.mean(cos_sims)),
            'cosine_sim_std': float(np.std(cos_sims)),
            'shift_mean': float(1.0 - np.mean(cos_sims)),
            'shift_std': float(np.std(cos_sims))
        }
        
        print(f"   ✅ {gene_symbol}: shift = {result['shift_mean']:.4f} ± {result['shift_std']:.4f} ({n_cells} cells)")
        return result
        
    except Exception as e:
        print(f"❌ Error computing metrics for {gene_symbol}: {e}")
        import traceback
        traceback.print_exc()
        return None

def debug_data_structure(data, gene_symbol):
    """Debug the data structure to understand what we're working with"""
    print(f"\n🔍 Debugging {gene_symbol} data structure:")
    print(f"Type: {type(data)}")
    
    if hasattr(data, 'keys'):
        for key in list(data.keys())[:3]:  # Show first 3 keys
            print(f"Key: {key}")
            value = data[key]
            print(f"  Value type: {type(value)}")
            
            if isinstance(value, dict):
                print(f"  Dict keys: {list(value.keys())}")
                for sub_key in list(value.keys())[:2]:  # Show first 2 sub-keys
                    sub_value = value[sub_key]
                    print(f"    {sub_key}: {type(sub_value)}")
                    if hasattr(sub_value, 'shape'):
                        print(f"      Shape: {sub_value.shape}")
            elif hasattr(value, 'shape'):
                print(f"  Shape: {value.shape}")
            elif hasattr(value, '__len__'):
                print(f"  Length: {len(value)}")

def main():
    """Analyze all 1000-cell results"""
    print("🔍 Analyzing 1000-cell perturbation results")
    print("=" * 60)
    
    # Results storage
    results_1000 = {}
    
    # Process each gene
    for gene in GENES:
        print(f"\n{'='*30}")
        print(f"Processing {gene}")
        print(f"{'='*30}")
        
        # Load data
        data = load_geneformer_result(gene)
        if data is None:
            continue
        
        # Debug data structure
        debug_data_structure(data, gene)
        
        # Extract embeddings
        perturbed_emb, control_emb = extract_embeddings_from_geneformer_data(data, gene)
        
        # Compute metrics
        if perturbed_emb is not None and control_emb is not None:
            metrics = compute_perturbation_metrics(perturbed_emb, control_emb, gene)
            if metrics is not None:
                results_1000[gene] = metrics
        else:
            print(f"❌ Could not extract embeddings for {gene}")
    
    # Save results
    output_file = RESULTS_DIR / "results_1000.json"
    with open(output_file, 'w') as f:
        json.dump(results_1000, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_file}")
    
    # Load baseline for comparison
    baseline_file = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/baseline_500_summary.json")
    baseline_results = {}
    
    if baseline_file.exists():
        with open(baseline_file, 'r') as f:
            baseline_results = json.load(f)
        print(f"📊 Loaded baseline results from: {baseline_file}")
    
    # Generate comparison report
    print("\n" + "=" * 60)
    print("📊 RESULTS COMPARISON")
    print("=" * 60)
    
    print("\n🎯 Baseline (500 cells):")
    print("| Gene   | Shift Mean | Shift Std | N Cells |")
    print("|--------|------------|-----------|---------|")
    for gene in GENES:
        if gene in baseline_results:
            data = baseline_results[gene]
            print(f"| {gene:6} | {data['shift_mean']:.4f}     | {data['shift_std']:.4f}    | {data['n_cells']:4} |")
        else:
            print(f"| {gene:6} | -          | -         | -    |")
    
    print("\n🚀 Scaled (1000 cells):")
    print("| Gene   | Shift Mean | Shift Std | N Cells |")
    print("|--------|------------|-----------|---------|")
    for gene in GENES:
        if gene in results_1000:
            data = results_1000[gene]
            print(f"| {gene:6} | {data['shift_mean']:.4f}     | {data['shift_std']:.4f}    | {data['n_cells_processed']:4} |")
        else:
            print(f"| {gene:6} | -          | -         | -    |")
    
    # Ranking comparison
    print("\n📈 GENE RANKINGS (by shift magnitude):")
    
    if baseline_results:
        baseline_ranking = sorted(
            [(gene, data['shift_mean']) for gene, data in baseline_results.items() if gene in GENES],
            key=lambda x: x[1], reverse=True
        )
        print(f"Baseline (500):  {' > '.join([f'{gene}({shift:.3f})' for gene, shift in baseline_ranking])}")
    
    if results_1000:
        scaled_ranking = sorted(
            [(gene, data['shift_mean']) for gene, data in results_1000.items() if gene in GENES],
            key=lambda x: x[1], reverse=True
        )
        print(f"Scaled (1000):   {' > '.join([f'{gene}({shift:.3f})' for gene, shift in scaled_ranking])}")
    
    print(f"\n✅ Analysis complete! Processed {len(results_1000)}/{len(GENES)} genes")

if __name__ == "__main__":
    main()