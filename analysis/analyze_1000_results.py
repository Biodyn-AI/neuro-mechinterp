#!/usr/bin/env python3
"""Analyze 1000-cell perturbation results and compare to baseline"""
import os
import json
import pickle
import numpy as np
from pathlib import Path
from sklearn.metrics.pairwise import cosine_similarity

# Genes and their paths
GENES = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
RESULTS_DIR = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled")

def load_perturbation_result(gene_symbol):
    """Load and analyze perturbation result for a gene"""
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
        if isinstance(data, dict):
            print(f"   Keys: {list(data.keys())}")
        
        return data
        
    except Exception as e:
        print(f"❌ Error loading {gene_symbol}: {e}")
        return None

def compute_perturbation_metrics(data, gene_symbol):
    """Compute perturbation metrics from loaded data"""
    if not isinstance(data, dict):
        print(f"❌ Expected dict for {gene_symbol}, got {type(data)}")
        return None
    
    # Try different possible key names for embeddings
    perturbed_keys = ['embeddings_perturbed', 'perturbed_embeddings', 'embeddings']
    baseline_keys = ['embeddings_baseline', 'baseline_embeddings', 'embeddings_control', 'control_embeddings']
    
    perturbed_emb = None
    baseline_emb = None
    
    # Find perturbed embeddings
    for key in perturbed_keys:
        if key in data:
            perturbed_emb = data[key]
            print(f"   Found perturbed embeddings: '{key}'")
            break
    
    # Find baseline embeddings 
    for key in baseline_keys:
        if key in data:
            baseline_emb = data[key]
            print(f"   Found baseline embeddings: '{key}'")
            break
    
    if perturbed_emb is None:
        print(f"❌ No perturbed embeddings found for {gene_symbol}")
        print(f"   Available keys: {list(data.keys())}")
        return None
        
    if baseline_emb is None:
        print(f"❌ No baseline embeddings found for {gene_symbol}")  
        print(f"   Available keys: {list(data.keys())}")
        return None
    
    try:
        # Convert to numpy arrays
        perturbed_emb = np.array(perturbed_emb)
        baseline_emb = np.array(baseline_emb)
        
        print(f"   Perturbed shape: {perturbed_emb.shape}")
        print(f"   Baseline shape: {baseline_emb.shape}")
        
        # Compute cosine similarities
        n_cells = min(len(perturbed_emb), len(baseline_emb))
        cos_sims = []
        
        for i in range(n_cells):
            # Ensure embeddings are 2D for cosine_similarity
            p_emb = perturbed_emb[i].reshape(1, -1)
            b_emb = baseline_emb[i].reshape(1, -1)
            
            sim = cosine_similarity(p_emb, b_emb)[0][0]
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
        return None

def main():
    """Analyze all 1000-cell results"""
    print("🔍 Analyzing 1000-cell perturbation results")
    print("=" * 60)
    
    # Results storage
    results_1000 = {}
    
    # Process each gene
    for gene in GENES:
        print(f"\nProcessing {gene}...")
        
        # Load data
        data = load_perturbation_result(gene)
        if data is None:
            continue
        
        # Compute metrics
        metrics = compute_perturbation_metrics(data, gene)
        if metrics is not None:
            results_1000[gene] = metrics
    
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
    else:
        print("⚠️  Baseline results not found")
    
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