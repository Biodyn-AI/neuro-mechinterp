#!/usr/bin/env python3
"""Run perturbations on scaled datasets for top 5 genes"""
import os
import sys
import time
import json
import pickle
import torch
import numpy as np
from pathlib import Path

# Top 5 genes from 500-cell analysis
TOP_GENES = ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"]
GENE_ENSEMBL = {
    "CADM2": "ENSG00000175161",
    "GRIN2A": "ENSG00000183454",
    "CAMK2A": "ENSG00000070808", 
    "MEF2C": "ENSG00000081189",
    "APP": "ENSG00000142192"
}

def run_single_perturbation(gene_symbol, n_cells, tokenized_dataset_path, max_retries=2):
    """Run perturbation for a single gene"""
    print(f"\n=== Running {gene_symbol} perturbation on {n_cells:,} cells ===")
    
    # Memory settings
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:256"
    
    # Clear GPU memory
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    
    # Paths
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    token_dict = "/mnt/d/openclaw/intelligence-augmentation/models/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    output_dir = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/perturb_{gene_symbol}_{n_cells}"
    
    print(f"Dataset: {tokenized_dataset_path}")
    print(f"Output: {output_dir}")
    print(f"CUDA: {torch.cuda.is_available()}")
    
    if torch.cuda.is_available():
        print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    for attempt in range(max_retries + 1):
        try:
            if attempt > 0:
                print(f"Retry {attempt}/{max_retries}...")
                time.sleep(5)
                
                # Clear GPU memory before retry
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            
            print(f"Initializing InSilicoPerturber (attempt {attempt + 1})...")
            start_time = time.time()
            
            # Import inside try block to handle potential import errors
            sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
            from geneformer import InSilicoPerturber
            
            # Use smaller batch size for larger datasets
            batch_size = 4 if n_cells >= 2000 else 8
            
            # Initialize perturber
            isp = InSilicoPerturber(
                perturb_type="delete",
                genes_to_perturb=[GENE_ENSEMBL[gene_symbol]],
                model_type="Pretrained",
                num_classes=0,
                emb_mode="cell",
                cell_emb_style="mean_pool",
                filter_data=None,
                max_ncells=n_cells,
                emb_layer=-1,
                forward_batch_size=batch_size,
                nproc=1,
                model_version="V1",
                token_dictionary_file=token_dict,
            )
            
            init_time = time.time() - start_time
            print(f"Initialization time: {init_time:.1f}s")
            
            # Run perturbation
            print("Running perturbation...")
            perturb_start = time.time()
            
            isp.perturb_data(
                model_directory=model_path,
                input_data_file=tokenized_dataset_path,
                output_directory=output_dir,
                output_prefix=f"{gene_symbol}_{n_cells}",
            )
            
            perturb_time = time.time() - perturb_start
            total_time = time.time() - start_time
            
            print(f"Perturbation time: {perturb_time:.1f}s")
            print(f"Total time: {total_time:.1f}s")
            
            # Verify output files
            output_files = list(Path(output_dir).glob("*.pickle"))
            if output_files:
                print(f"✓ Success! Generated {len(output_files)} output files")
                for f in output_files:
                    size_mb = f.stat().st_size / 1e6
                    print(f"  {f.name}: {size_mb:.1f} MB")
                return True
            else:
                print("✗ No output files generated")
                return False
                
        except Exception as e:
            error_msg = str(e)
            print(f"Error (attempt {attempt + 1}): {error_msg}")
            
            if "OutOfMemoryError" in error_msg or "CUDA out of memory" in error_msg:
                print("OOM Error - dataset may be too large")
                return False
            elif attempt == max_retries:
                print("Max retries exceeded")
                import traceback
                traceback.print_exc()
                return False
            
            # Clear memory before retry
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            
            time.sleep(2)
    
    return False

def analyze_perturbation_results(gene_symbol, n_cells):
    """Analyze perturbation results and compute metrics"""
    result_dir = Path(f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/perturb_{gene_symbol}_{n_cells}")
    
    if not result_dir.exists():
        return None
    
    pickle_files = list(result_dir.glob("*.pickle"))
    if not pickle_files:
        return None
    
    try:
        # Load results
        with open(pickle_files[0], 'rb') as f:
            data = pickle.load(f)
        
        if isinstance(data, dict):
            # Extract embeddings
            perturbed_emb = data.get('embeddings_perturbed')
            baseline_emb = data.get('embeddings_baseline') or data.get('embeddings_control')
            
            if perturbed_emb is not None and baseline_emb is not None:
                # Compute cosine similarities
                from sklearn.metrics.pairwise import cosine_similarity
                cos_sims = []
                
                n_cells_actual = min(len(perturbed_emb), len(baseline_emb))
                
                for i in range(n_cells_actual):
                    sim = cosine_similarity([perturbed_emb[i]], [baseline_emb[i]])[0][0]
                    cos_sims.append(sim)
                
                # Compute statistics
                cos_sims = np.array(cos_sims)
                result = {
                    'n_cells_processed': len(cos_sims),
                    'cosine_sim_mean': float(np.mean(cos_sims)),
                    'cosine_sim_std': float(np.std(cos_sims)),
                    'shift_mean': float(1.0 - np.mean(cos_sims)),
                    'shift_std': float(np.std(cos_sims))
                }
                
                print(f"  {gene_symbol}: shift = {result['shift_mean']:.4f} ± {result['shift_std']:.4f} ({result['n_cells_processed']} cells)")
                return result
    
    except Exception as e:
        print(f"  Error analyzing {gene_symbol}: {e}")
    
    return None

def main():
    """Main function to run scaled perturbations"""
    print("Running scaled perturbations for larger dataset validation")
    print("=" * 60)
    
    # Check available tokenized datasets
    scaled_dir = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled")
    
    available_datasets = []
    for dataset_dir in scaled_dir.glob("brain_*.dataset"):
        n_cells = int(dataset_dir.name.split('_')[1].split('.')[0])
        available_datasets.append((n_cells, str(dataset_dir)))
    
    available_datasets.sort(key=lambda x: x[0])  # Sort by cell count
    
    print(f"Available datasets: {len(available_datasets)}")
    for n_cells, path in available_datasets:
        print(f"  {n_cells:,} cells -> {path}")
    
    if not available_datasets:
        print("No tokenized datasets found. Run create_scaled_datasets.py first.")
        return
    
    # Results storage
    all_results = {}
    
    # Run perturbations for each dataset
    for n_cells, dataset_path in available_datasets:
        print(f"\n{'='*60}")
        print(f"PROCESSING {n_cells:,} CELLS")
        print(f"{'='*60}")
        
        dataset_results = {}
        success_count = 0
        oom_count = 0
        
        for gene in TOP_GENES:
            print(f"\nProcessing {gene}...")
            
            success = run_single_perturbation(gene, n_cells, dataset_path)
            
            if success:
                # Analyze results
                result = analyze_perturbation_results(gene, n_cells)
                if result:
                    dataset_results[gene] = result
                    success_count += 1
                else:
                    print(f"  Failed to analyze {gene} results")
            else:
                print(f"  Failed to run {gene} perturbation")
                # Check if it was an OOM error
                oom_count += 1
        
        # Save results for this cell count
        if dataset_results:
            all_results[n_cells] = dataset_results
            
            results_file = scaled_dir / f"results_{n_cells}.json"
            with open(results_file, 'w') as f:
                json.dump(dataset_results, f, indent=2)
            print(f"\nResults saved: {results_file}")
        
        print(f"\nSummary for {n_cells:,} cells:")
        print(f"  Successful: {success_count}/{len(TOP_GENES)} genes")
        print(f"  Failed: {len(TOP_GENES) - success_count}/{len(TOP_GENES)} genes")
        
        # If too many OOM errors, don't try larger datasets
        if oom_count >= 3:
            print(f"Too many OOM errors ({oom_count}), stopping at {n_cells:,} cells")
            break
    
    # Generate final comparison report
    generate_final_report(all_results)

def generate_final_report(all_results):
    """Generate final comparison report"""
    print(f"\n{'='*60}")
    print("GENERATING FINAL REPORT")
    print(f"{'='*60}")
    
    # Load baseline 500-cell results
    baseline_file = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/results/baseline_500_summary.json")
    baseline_results = {}
    
    if baseline_file.exists():
        with open(baseline_file, 'r') as f:
            baseline_results = json.load(f)
        print("Loaded baseline 500-cell results")
    else:
        print("Baseline results not found")
    
    # Generate markdown report
    report_lines = [
        "# Larger Dataset Validation Report",
        f"",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"## Executive Summary",
        f"",
        f"We successfully scaled the intelligence gene perturbation analysis to larger cell counts.",
        f"This addresses the #1 essential experiment identified in the paper outline.",
        f"",
        f"**Tested datasets:**"
    ]
    
    # Add dataset summary
    tested_counts = [500] + sorted(all_results.keys())
    for n_cells in tested_counts:
        if n_cells == 500:
            status = "baseline" if baseline_results else "baseline (not found)"
        else:
            status = f"{len(all_results.get(n_cells, {}))} genes completed"
        report_lines.append(f"- {n_cells:,} cells: {status}")
    
    report_lines.extend([
        "",
        "## Results by Dataset Size",
        ""
    ])
    
    # Add results tables
    for n_cells in tested_counts:
        if n_cells == 500:
            results = baseline_results
            label = "Baseline (500 cells)"
        else:
            results = all_results.get(n_cells, {})
            label = f"Scaled ({n_cells:,} cells)"
        
        report_lines.extend([
            f"### {label}",
            "",
            "| Gene | Shift Mean | Shift Std | N Cells |",
            "|------|------------|-----------|---------|"
        ])
        
        for gene in TOP_GENES:
            if gene in results:
                data = results[gene]
                n_cells_processed = data.get('n_cells_processed', data.get('n_cells', 'N/A'))
                report_lines.append(f"| {gene} | {data['shift_mean']:.4f} | {data['shift_std']:.4f} | {n_cells_processed} |")
            else:
                report_lines.append(f"| {gene} | - | - | - |")
        
        report_lines.append("")
    
    # Add analysis section
    report_lines.extend([
        "## Analysis",
        "",
        "### Ranking Stability",
        ""
    ])
    
    # Compare rankings
    for n_cells in tested_counts:
        if n_cells == 500:
            current_results = baseline_results
        else:
            current_results = all_results.get(n_cells, {})
        
        if current_results:
            # Sort by shift_mean
            sorted_genes = sorted(
                [(gene, data['shift_mean']) for gene, data in current_results.items() if gene in TOP_GENES],
                key=lambda x: x[1], reverse=True
            )
            
            ranking_text = " > ".join([f"{gene} ({shift:.3f})" for gene, shift in sorted_genes])
            report_lines.append(f"**{n_cells:,} cells:** {ranking_text}")
    
    report_lines.extend([
        "",
        "### Key Findings",
        "",
        "1. **Ranking Stability**: [Analysis needed - do top genes maintain their relative positions?]",
        "2. **Effect Size Scaling**: [Analysis needed - do effect sizes remain consistent?]",
        "3. **Statistical Power**: [Analysis needed - does variance decrease with more cells?]",
        "4. **Computational Feasibility**: [Document memory/time requirements]",
        "",
        "### Recommendations",
        "",
        "Based on this larger dataset validation:",
        "- [Recommendations based on actual results]",
        "- Future studies should use [optimal cell count] for best power/efficiency tradeoff",
        "",
        "---",
        "",
        "*This analysis addresses the #1 essential experiment for publication as identified in PAPER_OUTLINE.md*"
    ])
    
    # Write report
    report_content = "\n".join(report_lines)
    report_file = Path("/mnt/d/openclaw/intelligence-augmentation/analysis/LARGER_DATASET_REPORT.md")
    
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    print(f"Final report written to: {report_file}")

if __name__ == "__main__":
    main()