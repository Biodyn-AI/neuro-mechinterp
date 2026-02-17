#!/usr/bin/env python3
"""Compute 1000-cell perturbation shifts from pre-computed cosine similarities"""
import pickle
import numpy as np
import json
from pathlib import Path
import statistics

# Gene mappings
GENES = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
GENE_ENSEMBL = {
    "CADM2": 13697,
    "GRIN2A": 15103, 
    "CAMK2A": 1135,
    "MEF2C": 1532,
    "APP": 8172
}

def load_cosine_similarities(gene_symbol):
    """Load pre-computed cosine similarities from pickle files"""
    pickle_path = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/perturb_{gene_symbol}_1000"
    pickle_dir = Path(pickle_path)
    
    if not pickle_dir.exists():
        print(f"❌ Directory not found: {pickle_dir}")
        return None
    
    pickle_files = list(pickle_dir.glob("*.pickle"))
    if not pickle_files:
        print(f"❌ No pickle files found for {gene_symbol}")
        return None
    
    pickle_file = pickle_files[0]
    
    try:
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
        
        # Extract cosine similarities
        gene_token = GENE_ENSEMBL[gene_symbol]
        key = (np.int16(gene_token), 'cell_emb')
        
        if key not in data:
            print(f"❌ Expected key not found: {key}")
            print(f"   Available keys: {list(data.keys())}")
            return None
        
        cosine_sims = data[key]
        
        if not isinstance(cosine_sims, list):
            print(f"❌ Expected list, got {type(cosine_sims)}")
            return None
        
        # Convert to numpy array for easier processing
        cosine_sims = np.array(cosine_sims, dtype=float)
        
        print(f"✅ {gene_symbol}: loaded {len(cosine_sims)} cosine similarities")
        print(f"   Range: [{cosine_sims.min():.4f}, {cosine_sims.max():.4f}]")
        print(f"   Mean: {cosine_sims.mean():.4f} ± {cosine_sims.std():.4f}")
        
        return cosine_sims
        
    except Exception as e:
        print(f"❌ Error loading {gene_symbol}: {e}")
        return None

def compute_perturbation_metrics(cosine_sims, gene_symbol):
    """Compute perturbation metrics from cosine similarities"""
    if cosine_sims is None or len(cosine_sims) == 0:
        return None
    
    # Convert cosine similarities to shifts
    shifts = 1.0 - cosine_sims
    
    # Compute statistics
    metrics = {
        'n_cells_processed': len(shifts),
        'cosine_sim_mean': float(np.mean(cosine_sims)),
        'cosine_sim_std': float(np.std(cosine_sims)), 
        'cosine_sim_min': float(np.min(cosine_sims)),
        'cosine_sim_max': float(np.max(cosine_sims)),
        'shift_mean': float(np.mean(shifts)),
        'shift_std': float(np.std(shifts)),
        'shift_min': float(np.min(shifts)),
        'shift_max': float(np.max(shifts))
    }
    
    return metrics

def load_baseline_results():
    """Load baseline 500-cell results for comparison"""
    baseline_file = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/baseline_500_summary.json")
    
    if not baseline_file.exists():
        print("⚠️ Baseline results not found, will extract from FINAL_REPORT.md")
        
        # Extract baseline results from the FINAL_REPORT.md
        baseline_results = {
            "CADM2": {
                "shift_mean": 0.0196,
                "shift_std": 0.0045, 
                "n_cells": 472
            },
            "GRIN2A": {
                "shift_mean": 0.0190,
                "shift_std": 0.0061,
                "n_cells": 223
            },
            "CAMK2A": {
                "shift_mean": 0.0189,
                "shift_std": 0.0059,
                "n_cells": 155
            },
            "MEF2C": {
                "shift_mean": 0.0184,
                "shift_std": 0.0065,
                "n_cells": 216
            },
            "APP": {
                "shift_mean": 0.0183,
                "shift_std": 0.0054,
                "n_cells": 430
            }
        }
        return baseline_results
    
    try:
        with open(baseline_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading baseline: {e}")
        return {}

def main():
    """Analyze 1000-cell perturbation results"""
    print("🧬 Computing 1000-cell perturbation shifts")
    print("=" * 60)
    
    # Results storage
    results_1000 = {}
    
    # Process each gene
    for gene in GENES:
        print(f"\n📊 Processing {gene}...")
        
        # Load cosine similarities
        cosine_sims = load_cosine_similarities(gene)
        if cosine_sims is None:
            continue
        
        # Compute metrics
        metrics = compute_perturbation_metrics(cosine_sims, gene)
        if metrics is not None:
            results_1000[gene] = metrics
            print(f"   Shift: {metrics['shift_mean']:.4f} ± {metrics['shift_std']:.4f} ({metrics['n_cells_processed']} cells)")
    
    # Save results
    output_file = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/results_1000.json"
    with open(output_file, 'w') as f:
        json.dump(results_1000, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_file}")
    
    # Load baseline for comparison
    baseline_results = load_baseline_results()
    
    # Generate comparison report
    print("\n" + "=" * 80)
    print("📊 RESULTS COMPARISON: 500 vs 1000 cells")
    print("=" * 80)
    
    print("\n🎯 Baseline (500 cells):")
    print("| Gene   | Shift Mean | Shift Std | N Cells | Pathway |")
    print("|--------|------------|-----------|---------|---------|")
    for gene in GENES:
        if gene in baseline_results:
            data = baseline_results[gene]
            shift_mean = data.get('shift_mean', 0)
            shift_std = data.get('shift_std', 0)
            n_cells = data.get('n_cells', data.get('n_cells_processed', 0))
            # Map to pathways from the original report
            pathways = {
                'CADM2': 'Cell Adhesion',
                'GRIN2A': 'Glutamate Receptors',
                'CAMK2A': 'Neurotrophic Signaling',
                'MEF2C': 'Neurodevelopmental/TF',
                'APP': 'Neurodegeneration'
            }
            pathway = pathways.get(gene, 'Unknown')
            print(f"| {gene:6} | {shift_mean:.4f}     | {shift_std:.4f}    | {n_cells:4} | {pathway:15} |")
    
    print("\n🚀 Scaled (1000 cells):")
    print("| Gene   | Shift Mean | Shift Std | N Cells | Change vs 500 |")
    print("|--------|------------|-----------|---------|---------------|")
    for gene in GENES:
        if gene in results_1000:
            data = results_1000[gene]
            baseline_shift = baseline_results.get(gene, {}).get('shift_mean', 0)
            change = data['shift_mean'] - baseline_shift if baseline_shift > 0 else 0
            change_pct = (change / baseline_shift * 100) if baseline_shift > 0 else 0
            print(f"| {gene:6} | {data['shift_mean']:.4f}     | {data['shift_std']:.4f}    | {data['n_cells_processed']:4} | {change:+.4f} ({change_pct:+.1f}%) |")
    
    # Ranking comparison
    print("\n📈 GENE RANKINGS (by shift magnitude):")
    
    if baseline_results:
        baseline_ranking = sorted(
            [(gene, data['shift_mean']) for gene, data in baseline_results.items() if gene in GENES],
            key=lambda x: x[1], reverse=True
        )
        print(f"Baseline (500):  {' > '.join([f'{gene}({shift:.4f})' for gene, shift in baseline_ranking])}")
    
    if results_1000:
        scaled_ranking = sorted(
            [(gene, data['shift_mean']) for gene, data in results_1000.items()],
            key=lambda x: x[1], reverse=True
        )
        print(f"Scaled (1000):   {' > '.join([f'{gene}({shift:.4f})' for gene, shift in scaled_ranking])}")
    
    # Ranking concordance analysis
    if baseline_results and results_1000:
        print("\n🎯 RANKING CONCORDANCE:")
        baseline_order = [gene for gene, _ in baseline_ranking]
        scaled_order = [gene for gene, _ in scaled_ranking]
        
        print(f"Baseline order: {' > '.join(baseline_order)}")
        print(f"Scaled order:   {' > '.join(scaled_order)}")
        
        # Calculate rank correlations
        baseline_ranks = {gene: i for i, gene in enumerate(baseline_order)}
        scaled_ranks = {gene: i for i, gene in enumerate(scaled_order)}
        
        rank_diffs = []
        for gene in GENES:
            if gene in baseline_ranks and gene in scaled_ranks:
                diff = abs(baseline_ranks[gene] - scaled_ranks[gene])
                rank_diffs.append(diff)
                print(f"  {gene}: baseline rank {baseline_ranks[gene]+1}, scaled rank {scaled_ranks[gene]+1}, diff = {diff}")
        
        mean_rank_diff = np.mean(rank_diffs) if rank_diffs else 0
        print(f"  Average rank difference: {mean_rank_diff:.1f}")
        
        if mean_rank_diff < 1:
            print("  ✅ Rankings are highly concordant")
        elif mean_rank_diff < 2:
            print("  ⚠️ Rankings are moderately concordant")
        else:
            print("  ❌ Rankings show significant differences")
    
    # Effect size analysis
    print("\n📏 EFFECT SIZE COMPARISON:")
    if baseline_results and results_1000:
        print("| Gene   | 500-cell | 1000-cell | Ratio | Interpretation |")
        print("|--------|----------|-----------|-------|----------------|")
        
        for gene in GENES:
            if gene in baseline_results and gene in results_1000:
                base_shift = baseline_results[gene]['shift_mean']
                scaled_shift = results_1000[gene]['shift_mean']
                ratio = scaled_shift / base_shift if base_shift > 0 else 0
                
                if ratio > 1.1:
                    interp = "Larger effect"
                elif ratio < 0.9:
                    interp = "Smaller effect"  
                else:
                    interp = "Similar effect"
                
                print(f"| {gene:6} | {base_shift:.4f}   | {scaled_shift:.4f}    | {ratio:.3f} | {interp:14} |")
    
    print(f"\n✅ Analysis complete! Processed {len(results_1000)}/{len(GENES)} genes")
    
    # Statistical summary
    if results_1000:
        all_shifts_1000 = [data['shift_mean'] for data in results_1000.values()]
        all_shifts_500 = [baseline_results.get(gene, {}).get('shift_mean', 0) 
                          for gene in results_1000.keys()]
        
        mean_shift_1000 = np.mean(all_shifts_1000)
        mean_shift_500 = np.mean([x for x in all_shifts_500 if x > 0])
        
        print(f"\n📊 STATISTICAL SUMMARY:")
        print(f"Mean shift (500 cells):  {mean_shift_500:.4f}")
        print(f"Mean shift (1000 cells): {mean_shift_1000:.4f}")
        print(f"Overall effect ratio:    {mean_shift_1000/mean_shift_500:.3f}")

if __name__ == "__main__":
    main()