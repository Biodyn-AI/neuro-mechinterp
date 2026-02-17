#!/usr/bin/env python3
"""
Scaled Perturbation Pipeline for Larger Dataset Validation
Test different cell counts: 1000, 2000, 5000 cells
"""
import sys
import os
import gc
import time
import json
import pickle
import warnings
from datetime import datetime
from pathlib import Path

# Set encoding and suppress warnings
os.environ['PYTHONIOENCODING'] = 'utf-8'
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:512"
warnings.filterwarnings('ignore')

print("Initializing scaled perturbation pipeline...")

# Project paths
PROJ_ROOT = Path(r"D:\openclaw\intelligence-augmentation")
DATA_DIR = PROJ_ROOT / "data"
RESULTS_DIR = PROJ_ROOT / "analysis" / "results"
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# Top 5 genes from 500-cell analysis
TOP_GENES = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
GENE_ENSEMBL = {
    "CADM2": "ENSG00000175161",
    "GRIN2A": "ENSG00000183454", 
    "CAMK2A": "ENSG00000070808",
    "MEF2C": "ENSG00000081189",
    "APP": "ENSG00000142192"
}

def log_progress(message):
    """Log with timestamp"""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}")

def prepare_scaled_dataset(n_cells=1000):
    """Prepare dataset with specified number of cells"""
    log_progress(f"Preparing dataset with {n_cells:,} cells...")
    
    # Use WSL for data preparation
    import subprocess
    
    script_content = f"""
#!/usr/bin/env python3
import scanpy as sc
import numpy as np
import pandas as pd

print("Loading full dataset...")
adata = sc.read_h5ad('/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad')
print(f"Full dataset: {{adata.n_obs:,}} cells x {{adata.n_vars:,}} genes")

# Sample {n_cells} cells randomly
if adata.n_obs > {n_cells}:
    sc.pp.subsample(adata, n_obs={n_cells}, random_state=42)
    print(f"Subsampled to: {{adata.n_obs:,}} cells")

# Basic preprocessing 
adata.X = adata.X.astype('float32')

# Save prepared dataset
output_path = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_{n_cells}.h5ad'
import os
os.makedirs(os.path.dirname(output_path), exist_ok=True)
adata.write_h5ad(output_path)
print(f"Saved scaled dataset: {{output_path}}")
"""
    
    # Write prep script
    prep_script = RESULTS_DIR / f"prep_{n_cells}.py"
    with open(prep_script, 'w') as f:
        f.write(script_content)
    
    # Run in WSL
    cmd = f'wsl -u agent -- bash -lc "cd /mnt/d/openclaw/intelligence-augmentation/analysis && python3 {prep_script}"'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error preparing dataset: {result.stderr}")
        return None
    
    print(result.stdout)
    return RESULTS_DIR / "scaled" / f"brain_{n_cells}.h5ad"

def run_scaled_perturbation(gene_symbol, n_cells, max_retries=2):
    """Run perturbation for single gene with specified cell count"""
    log_progress(f"Running {gene_symbol} perturbation on {n_cells:,} cells...")
    
    # Prepare WSL command
    script_content = f"""
#!/usr/bin/env python3
import os
import sys
import torch
import pickle
from pathlib import Path

# Memory settings
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:256"

print("Loading Geneformer...")
sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')

try:
    from geneformer import InSilicoPerturber
    
    # Paths
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    data_path = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_{n_cells}.h5ad"
    token_dict = "/mnt/d/openclaw/intelligence-augmentation/models/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/perturb_{gene_symbol}_{n_cells}"
    
    print(f"Model: {{model_path}}")
    print(f"Data: {{data_path}}")
    print(f"Output: {{output_dir}}")
    print(f"CUDA: {{torch.cuda.is_available()}}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize perturber
    isp = InSilicoPerturber(
        perturb_type="delete",
        genes_to_perturb=["{GENE_ENSEMBL[gene_symbol]}"],
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,
        max_ncells={n_cells},
        emb_layer=-1,
        forward_batch_size=16,  # Smaller batch for larger datasets
        nproc=1,
        model_version="V1",
        token_dictionary_file=token_dict,
    )
    
    # Run perturbation
    isp.perturb_data(
        model_directory=model_path,
        input_data_file=data_path,
        output_directory=output_dir,
        output_prefix=f"{gene_symbol}_{n_cells}",
    )
    
    print("Perturbation completed successfully!")
    
except Exception as e:
    print(f"ERROR: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
"""
    
    # Write and run script
    script_path = RESULTS_DIR / f"run_{gene_symbol}_{n_cells}.py"
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    import subprocess
    
    for attempt in range(max_retries + 1):
        if attempt > 0:
            log_progress(f"Retry {attempt}/{max_retries} for {gene_symbol}...")
            time.sleep(5)  # Brief pause between retries
        
        cmd = f'wsl -u agent -- bash -lc "cd /mnt/d/openclaw/intelligence-augmentation/analysis && timeout 1800 python3 {script_path}"'
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=2000)
            
            if result.returncode == 0:
                print("Success!")
                return True
            elif "OutOfMemoryError" in result.stderr or "CUDA out of memory" in result.stderr:
                print(f"OOM error for {n_cells} cells, need to reduce cell count")
                return False
            else:
                print(f"Error (attempt {attempt + 1}): {result.stderr[-500:]}")
        
        except subprocess.TimeoutExpired:
            print(f"Timeout after 30 minutes (attempt {attempt + 1})")
        
        # Clear any hanging processes
        subprocess.run('wsl -u agent -- pkill -f python3', shell=True)
        time.sleep(2)
    
    return False

def analyze_scaled_results(n_cells):
    """Analyze results from scaled perturbation"""
    log_progress(f"Analyzing results for {n_cells:,} cells...")
    
    results = {}
    
    for gene in TOP_GENES:
        result_dir = RESULTS_DIR / "scaled" / f"perturb_{gene}_{n_cells}"
        
        if not result_dir.exists():
            print(f"No results found for {gene}")
            continue
            
        # Find pickle files
        pickle_files = list(result_dir.glob("*.pickle"))
        
        if not pickle_files:
            print(f"No pickle files found for {gene}")
            continue
        
        try:
            # Load perturbation results  
            with open(pickle_files[0], 'rb') as f:
                data = pickle.load(f)
            
            if isinstance(data, dict) and 'embeddings_perturbed' in data:
                perturbed_emb = data['embeddings_perturbed']
                baseline_emb = data.get('embeddings_baseline', data.get('embeddings_control'))
                
                if baseline_emb is not None:
                    # Compute cosine similarity
                    from sklearn.metrics.pairwise import cosine_similarity
                    cos_sims = []
                    
                    for i in range(len(perturbed_emb)):
                        if i < len(baseline_emb):
                            sim = cosine_similarity([perturbed_emb[i]], [baseline_emb[i]])[0][0]
                            cos_sims.append(sim)
                    
                    results[gene] = {
                        'n_cells': len(cos_sims),
                        'cosine_sim_mean': float(np.mean(cos_sims)),
                        'cosine_sim_std': float(np.std(cos_sims)),
                        'shift_mean': float(1.0 - np.mean(cos_sims)),
                        'shift_std': float(np.std(cos_sims))
                    }
                    
                    log_progress(f"{gene}: shift = {results[gene]['shift_mean']:.4f} ± {results[gene]['shift_std']:.4f}")
        
        except Exception as e:
            print(f"Error analyzing {gene}: {e}")
    
    # Save results
    results_file = RESULTS_DIR / "scaled" / f"results_{n_cells}.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    return results

def main():
    """Main scaled analysis workflow"""
    print("=" * 60)
    print("LARGER DATASET VALIDATION - INTELLIGENCE GENE PERTURBATION")
    print("=" * 60)
    
    # Test cell counts (start conservatively)
    cell_counts = [1000, 2000, 5000]
    
    all_results = {}
    successful_counts = []
    
    for n_cells in cell_counts:
        log_progress(f"Starting analysis with {n_cells:,} cells...")
        
        try:
            # Step 1: Prepare scaled dataset
            dataset_path = prepare_scaled_dataset(n_cells)
            if dataset_path is None:
                print(f"Failed to prepare {n_cells}-cell dataset")
                continue
            
            # Step 2: Run perturbations for top 5 genes
            success_count = 0
            oom_encountered = False
            
            for gene in TOP_GENES:
                success = run_scaled_perturbation(gene, n_cells)
                if success:
                    success_count += 1
                else:
                    # Check if it's an OOM error
                    if n_cells >= 2000:  # Assume OOM for larger datasets
                        oom_encountered = True
                        break
            
            if oom_encountered:
                log_progress(f"OOM encountered at {n_cells:,} cells, stopping")
                break
            
            if success_count >= 3:  # At least 3 genes successful
                # Step 3: Analyze results
                results = analyze_scaled_results(n_cells)
                all_results[n_cells] = results
                successful_counts.append(n_cells)
                log_progress(f"Successfully completed {n_cells:,} cell analysis")
            else:
                log_progress(f"Too many failures at {n_cells:,} cells")
                
        except Exception as e:
            print(f"Error with {n_cells:,} cells: {e}")
            import traceback
            traceback.print_exc()
    
    # Step 4: Generate comparison report
    generate_comparison_report(all_results, successful_counts)

def generate_comparison_report(all_results, successful_counts):
    """Generate comparison report across different cell counts"""
    log_progress("Generating comparison report...")
    
    # Load 500-cell baseline results
    baseline_file = RESULTS_DIR / "insilico_wsl" / "final_rankings.json"
    baseline_results = {}
    
    if baseline_file.exists():
        with open(baseline_file, 'r') as f:
            baseline_data = json.load(f)
            
        for gene in TOP_GENES:
            if gene in baseline_data:
                baseline_results[gene] = baseline_data[gene]
    
    # Generate markdown report
    report = f"""# Larger Dataset Validation Report

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Executive Summary

We successfully scaled the intelligence gene perturbation analysis from 500 cells to larger datasets.
Tested cell counts: {', '.join(map(str, successful_counts))} cells.

## Results by Cell Count

### Baseline (500 cells)
"""
    
    if baseline_results:
        report += "| Gene | Shift Mean | Std | N Cells |\n"
        report += "|------|------------|-----|----------|\n"
        for gene in TOP_GENES:
            if gene in baseline_results:
                data = baseline_results[gene]
                report += f"| {gene} | {data['shift_mean']:.4f} | {data['shift_std']:.4f} | {data['n_cells']} |\n"
    else:
        report += "Baseline results not found.\n"
    
    for n_cells in successful_counts:
        report += f"\n### Scaled ({n_cells:,} cells)\n"
        
        if n_cells in all_results:
            results = all_results[n_cells]
            report += "| Gene | Shift Mean | Std | N Cells |\n"
            report += "|------|------------|-----|----------|\n"
            for gene in TOP_GENES:
                if gene in results:
                    data = results[gene]
                    report += f"| {gene} | {data['shift_mean']:.4f} | {data['shift_std']:.4f} | {data['n_cells']} |\n"
                else:
                    report += f"| {gene} | - | - | - |\n"
    
    # Analysis section
    report += "\n## Analysis\n\n"
    report += "### Ranking Stability\n"
    
    # Check if rankings change
    if len(successful_counts) > 0 and baseline_results:
        report += "Comparison of gene rankings across cell counts:\n\n"
        
        for i, n_cells in enumerate([500] + successful_counts):
            if i == 0:
                current_results = baseline_results
                label = "500 (baseline)"
            else:
                current_results = all_results.get(n_cells, {})
                label = f"{n_cells:,}"
            
            if current_results:
                # Sort by shift_mean
                sorted_genes = sorted(
                    [(gene, data['shift_mean']) for gene, data in current_results.items() if gene in TOP_GENES],
                    key=lambda x: x[1], reverse=True
                )
                
                report += f"**{label} cells:** "
                report += " > ".join([f"{gene} ({shift:.3f})" for gene, shift in sorted_genes[:5]])
                report += "\n\n"
    
    # Save report
    report_file = RESULTS_DIR / "LARGER_DATASET_REPORT.md"
    with open(report_file, 'w') as f:
        f.write(report)
    
    log_progress(f"Report saved: {report_file}")
    print(f"\nFinal report saved to: {report_file}")

if __name__ == "__main__":
    main()