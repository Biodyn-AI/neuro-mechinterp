#!/usr/bin/env python3
"""
scGPT Perturbation Experiment - Priority 1

This script addresses the paper's #1 weakness: cross-model validation failure.
We perform perturbation analysis on the same 5 intelligence genes (CADM2, GRIN2A, 
CAMK2A, MEF2C, APP) using scGPT to compare with Geneformer results.

CRITICAL: This script MUST start with the torchtext patch for scGPT compatibility.
"""

# CRITICAL: Load torchtext shim first
exec(open("/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py").read())

import os
import sys
import json
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import warnings
warnings.filterwarnings('ignore')

# Add scGPT to path
sys.path.insert(0, '/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT/scgpt')

def main():
    print("=" * 60)
    print("scGPT Perturbation Experiment - Priority 1")
    print("=" * 60)
    
    # Configuration
    data_path = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
    checkpoint_dir = "/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT_checkpoints/brain/"
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results"
    
    # Top 5 intelligence genes for cross-validation
    intelligence_genes = ['CADM2', 'GRIN2A', 'CAMK2A', 'MEF2C', 'APP']
    
    print(f"Loading brain data from: {data_path}")
    
    # Load and prepare data
    try:
        adata = sc.read_h5ad(data_path)
        print(f"Loaded data: {adata.shape}")
        print(f"Available genes: {len(adata.var_names)}")
        
        # Check which intelligence genes are present
        available_genes = []
        for gene in intelligence_genes:
            if gene in adata.var_names:
                available_genes.append(gene)
                print(f"✓ {gene} found")
            else:
                print(f"✗ {gene} not found")
        
        if len(available_genes) == 0:
            print("ERROR: No intelligence genes found in dataset!")
            return
            
        print(f"Will analyze {len(available_genes)} genes: {available_genes}")
        
    except Exception as e:
        print(f"ERROR loading data: {e}")
        return
    
    # Load scGPT model
    print(f"\nLoading scGPT model from: {checkpoint_dir}")
    
    try:
        # Import scGPT modules
        from scgpt.model import TransformerModel
        from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
        from scgpt.tokenizer.gene_tokenizer import GeneVocab
        
        # Load vocab
        vocab_file = os.path.join(checkpoint_dir, "vocab.json")
        with open(vocab_file, 'r') as f:
            vocab_data = json.load(f)
        
        vocab = GeneVocab.from_dict(vocab_data)
        print(f"Vocab loaded: {len(vocab)} genes")
        
        # Load model checkpoint
        model_file = os.path.join(checkpoint_dir, "best_model.pt")
        ckpt = torch.load(model_file, map_location=torch.device('cpu'))
        print("Model checkpoint loaded")
        
        # Initialize model with proper parameters
        model = TransformerModel(
            ntoken=len(vocab),
            d_model=512,
            nhead=8,
            d_hid=512,
            nlayers=12,
            vocab=vocab,
            dropout=0.2,
            pad_token=vocab[vocab.pad_token],
            pad_value=0,
            do_mvc=True,
            do_dab=False,
            use_batch_labels=False,
            num_batch_labels=None,
            domain_spec_batchnorm=False,
            explicit_zero_prob=False,
        )
        
        model.load_state_dict(ckpt['model_state_dict'])
        model.eval()
        print("Model initialized and loaded")
        
    except Exception as e:
        print(f"ERROR loading scGPT model: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Prepare data for scGPT
    print("\nPreparing data for scGPT...")
    
    try:
        # Select a subset of cells for perturbation (matching Geneformer analysis)
        n_cells = min(500, adata.shape[0])
        cell_indices = np.random.choice(adata.shape[0], n_cells, replace=False)
        adata_subset = adata[cell_indices, :].copy()
        
        print(f"Selected {n_cells} cells for perturbation analysis")
        
        # Convert to count data if needed
        if 'X_raw' in adata_subset.layers:
            adata_subset.X = adata_subset.layers['X_raw']
        elif hasattr(adata_subset.X, 'todense'):
            adata_subset.X = adata_subset.X.todense()
        
        # Ensure integer counts
        adata_subset.X = np.array(adata_subset.X, dtype=int)
        
    except Exception as e:
        print(f"ERROR preparing data: {e}")
        return
    
    # Run perturbation analysis
    print(f"\nRunning perturbation analysis on {len(available_genes)} genes...")
    
    perturbation_results = {}
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model.to(device)
    
    for gene in available_genes:
        print(f"\nPerturbing {gene}...")
        
        try:
            # Get gene index in vocab
            if gene not in vocab:
                print(f"WARNING: {gene} not in model vocabulary, skipping")
                continue
                
            gene_idx = vocab[gene]
            
            # Get cells expressing this gene
            gene_col_idx = list(adata_subset.var_names).index(gene)
            expressing_cells = np.where(adata_subset.X[:, gene_col_idx] > 0)[0]
            
            if len(expressing_cells) == 0:
                print(f"WARNING: No cells express {gene}, skipping")
                continue
                
            print(f"Found {len(expressing_cells)} cells expressing {gene}")
            
            # Perform perturbation on expressing cells
            embedding_shifts = []
            
            for cell_idx in expressing_cells[:50]:  # Limit to first 50 for computational efficiency
                # Get original cell expression
                cell_expr = adata_subset.X[cell_idx, :].flatten()
                
                # Create perturbed version (set gene expression to 0)
                perturbed_expr = cell_expr.copy()
                perturbed_expr[gene_col_idx] = 0
                
                # Tokenize both versions
                original_tokens = []
                perturbed_tokens = []
                
                # Convert to gene tokens (simplified tokenization)
                for i, expr_val in enumerate(cell_expr):
                    if expr_val > 0:
                        gene_name = adata_subset.var_names[i]
                        if gene_name in vocab:
                            original_tokens.append(vocab[gene_name])
                
                for i, expr_val in enumerate(perturbed_expr):
                    if expr_val > 0:
                        gene_name = adata_subset.var_names[i]
                        if gene_name in vocab:
                            perturbed_tokens.append(vocab[gene_name])
                
                # Pad sequences to same length
                max_len = 1200  # scGPT typical max length
                if len(original_tokens) > max_len:
                    original_tokens = original_tokens[:max_len]
                if len(perturbed_tokens) > max_len:
                    perturbed_tokens = perturbed_tokens[:max_len]
                
                # Pad with pad token
                pad_token_id = vocab[vocab.pad_token]
                while len(original_tokens) < max_len:
                    original_tokens.append(pad_token_id)
                while len(perturbed_tokens) < max_len:
                    perturbed_tokens.append(pad_token_id)
                
                # Convert to tensors
                original_tensor = torch.tensor([original_tokens], dtype=torch.long).to(device)
                perturbed_tensor = torch.tensor([perturbed_tokens], dtype=torch.long).to(device)
                
                # Get embeddings
                with torch.no_grad():
                    original_output = model(original_tensor, None, None, None)
                    perturbed_output = model(perturbed_tensor, None, None, None)
                    
                    # Extract embeddings (use CLS token embedding)
                    original_emb = original_output['cell_emb'][0].cpu().numpy()
                    perturbed_emb = perturbed_output['cell_emb'][0].cpu().numpy()
                    
                    # Calculate cosine similarity
                    cos_sim = np.dot(original_emb, perturbed_emb) / (np.linalg.norm(original_emb) * np.linalg.norm(perturbed_emb))
                    embedding_shift = 1 - cos_sim
                    embedding_shifts.append(embedding_shift)
            
            # Store results
            if len(embedding_shifts) > 0:
                perturbation_results[gene] = {
                    'embedding_shifts': embedding_shifts,
                    'mean_shift': np.mean(embedding_shifts),
                    'std_shift': np.std(embedding_shifts),
                    'n_cells': len(embedding_shifts)
                }
                
                print(f"{gene}: mean_shift={np.mean(embedding_shifts):.4f} ± {np.std(embedding_shifts):.4f} (n={len(embedding_shifts)})")
            
        except Exception as e:
            print(f"ERROR perturbing {gene}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save results
    output_file = os.path.join(output_dir, "scgpt_perturbation_results.json")
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert numpy arrays to lists for JSON serialization
    json_results = {}
    for gene, data in perturbation_results.items():
        json_results[gene] = {
            'embedding_shifts': [float(x) for x in data['embedding_shifts']],
            'mean_shift': float(data['mean_shift']),
            'std_shift': float(data['std_shift']),
            'n_cells': int(data['n_cells'])
        }
    
    with open(output_file, 'w') as f:
        json.dump(json_results, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    
    # Print summary
    print("\n" + "="*60)
    print("scGPT PERTURBATION RESULTS SUMMARY")
    print("="*60)
    
    if perturbation_results:
        # Sort by mean embedding shift
        sorted_genes = sorted(perturbation_results.items(), 
                            key=lambda x: x[1]['mean_shift'], 
                            reverse=True)
        
        print(f"{'Rank':<4} {'Gene':<8} {'Mean Shift':<12} {'Std':<8} {'N Cells':<8}")
        print("-" * 50)
        
        for rank, (gene, data) in enumerate(sorted_genes, 1):
            print(f"{rank:<4} {gene:<8} {data['mean_shift']:<12.6f} "
                  f"{data['std_shift']:<8.4f} {data['n_cells']:<8}")
        
        # Save ranking for comparison
        ranking_file = os.path.join(output_dir, "scgpt_gene_ranking.csv")
        ranking_df = pd.DataFrame([
            {'rank': rank, 'gene': gene, 'mean_shift': data['mean_shift'], 
             'std_shift': data['std_shift'], 'n_cells': data['n_cells']}
            for rank, (gene, data) in enumerate(sorted_genes, 1)
        ])
        ranking_df.to_csv(ranking_file, index=False)
        print(f"\nGene ranking saved to: {ranking_file}")
        
    else:
        print("No perturbation results generated!")
    
    print("\nscGPT perturbation experiment completed!")

if __name__ == "__main__":
    main()