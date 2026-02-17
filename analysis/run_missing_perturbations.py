#!/usr/bin/env python3
"""Run missing perturbations for MEF2C and APP on 1000 cells"""
import os
import sys
import time
import torch

# Top 5 genes with their Ensembl IDs
GENE_ENSEMBL = {
    "CADM2": "ENSG00000175161",
    "GRIN2A": "ENSG00000183454", 
    "CAMK2A": "ENSG00000070808",
    "MEF2C": "ENSG00000081189",
    "APP": "ENSG00000142192"
}

MISSING_GENES = ["MEF2C", "APP"]

def run_perturbation(gene_symbol, n_cells=1000):
    """Run perturbation for a single gene"""
    print(f"\n{'='*60}")
    print(f"Running {gene_symbol} perturbation on {n_cells:,} cells")
    print(f"{'='*60}")
    
    # Memory settings for RTX 2060 (6GB)
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:256"
    
    # Clear GPU memory
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"Available memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    # Paths
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    token_dict = "/mnt/d/openclaw/intelligence-augmentation/models/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
    dataset_path = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/brain_{n_cells}.dataset"
    output_dir = f"/mnt/d/openclaw/intelligence-augmentation/analysis/results/scaled/perturb_{gene_symbol}_{n_cells}"
    
    print(f"Dataset: {dataset_path}")
    print(f"Output: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        print("Loading Geneformer...")
        start_time = time.time()
        
        # Import Geneformer
        sys.path.append('/mnt/d/openclaw/intelligence-augmentation/models/geneformer')
        from geneformer import InSilicoPerturber
        
        # Conservative settings for 6GB GPU
        batch_size = 4  # Small batch for memory efficiency
        
        print(f"Initializing InSilicoPerturber for {gene_symbol}...")
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
        
        print("Running perturbation...")
        perturb_start = time.time()
        
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=dataset_path,
            output_directory=output_dir,
            output_prefix=f"{gene_symbol}_{n_cells}",
        )
        
        perturb_time = time.time() - perturb_start
        total_time = time.time() - start_time
        
        print(f"Perturbation time: {perturb_time:.1f}s")
        print(f"Total time: {total_time:.1f}s")
        
        # Check output files
        import glob
        output_files = glob.glob(os.path.join(output_dir, "*.pickle"))
        
        if output_files:
            print(f"✓ Success! Generated {len(output_files)} output files:")
            for f in output_files:
                size_mb = os.path.getsize(f) / 1e6
                print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
            return True
        else:
            print("✗ No output files generated")
            return False
            
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run missing perturbations"""
    print("Running missing 1000-cell perturbations")
    print("="*60)
    
    success_count = 0
    
    for gene in MISSING_GENES:
        success = run_perturbation(gene, 1000)
        if success:
            success_count += 1
            print(f"✓ {gene} completed successfully")
        else:
            print(f"✗ {gene} failed")
        
        # Clear memory between runs
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        # Brief pause between runs
        time.sleep(2)
    
    print(f"\n{'='*60}")
    print(f"SUMMARY: {success_count}/{len(MISSING_GENES)} genes completed")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()