#!/usr/bin/env python3
"""
scGPT Perturbation Experiment - Fixed Version

This script addresses the paper's #1 weakness: cross-model validation failure.
We perform perturbation analysis on the same 5 intelligence genes using scGPT 
to compare with Geneformer results.

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
from scipy.spatial.distance import cosine
warnings.filterwarnings('ignore')

# Add scGPT to path
sys.path.insert(0, '/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT/scgpt')

def main():
    print("=" * 60)
    print("scGPT Perturbation Experiment - Fixed Version")
    print("=" * 60)
    
    # Configuration
    data_path = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
    checkpoint_dir = "/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT_checkpoints/brain/"
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results"
    
    # Top 5 intelligence genes for cross-validation (using Ensembl IDs)
    gene_mapping = {
        'CADM2': 'ENSG00000175161',
        'GRIN2A': 'ENSG00000183454', 
        'CAMK2A': 'ENSG00000070808',
        'MEF2C': 'ENSG00000081189',
        'APP': 'ENSG00000142192'
    }
    
    print(f"Loading brain data from: {data_path}")
    
    # Load and prepare data
    try:
        adata = sc.read_h5ad(data_path)
        print(f"Loaded data: {adata.shape}")
        
        # Verify all genes are present
        for symbol, ensembl_id in gene_mapping.items():
            if ensembl_id in adata.var_names:
                print(f"✓ {symbol} ({ensembl_id}) found")
            else:
                print(f"✗ {symbol} ({ensembl_id}) not found")
                return
        
        print("All intelligence genes found!")
        
    except Exception as e:
        print(f"ERROR loading data: {e}")
        return
    
    # Load scGPT model
    print(f"\nLoading scGPT model from: {checkpoint_dir}")
    
    try:
        # Import scGPT modules
        from scgpt.model import TransformerModel
        from scgpt.tokenizer.gene_tokenizer import GeneVocab
        
        # Load vocab
        vocab_file = os.path.join(checkpoint_dir, "vocab.json")
        with open(vocab_file, 'r') as f:
            vocab_data = json.load(f)
        
        vocab = GeneVocab.from_dict(vocab_data)
        print(f"Vocab loaded: {len(vocab)} genes")
        
        # Check which genes are in vocab
        genes_in_vocab = []
        for symbol, ensembl_id in gene_mapping.items():
            # Try both symbol and Ensembl ID
            if symbol in vocab:
                genes_in_vocab.append((symbol, symbol))
                print(f"✓ {symbol} found in vocab as symbol")
            elif ensembl_id in vocab:
                genes_in_vocab.append((symbol, ensembl_id))
                print(f"✓ {symbol} found in vocab as {ensembl_id}")
            else:
                print(f"✗ {symbol} not found in vocab")
        
        if len(genes_in_vocab) == 0:
            print("ERROR: No intelligence genes found in model vocabulary!")
            return
        
        print(f"Will analyze {len(genes_in_vocab)} genes")
        
        # Load model checkpoint
        model_file = os.path.join(checkpoint_dir, "best_model.pt")
        ckpt = torch.load(model_file, map_location=torch.device('cpu'))
        
        # Load model args from checkpoint
        model_args = ckpt.get('args', None)
        if model_args is None:
            print("No model args in checkpoint, using defaults")
            model_args = {
                'd_model': 512,
                'nhead': 8,
                'd_hid': 512,
                'nlayers': 12,
                'dropout': 0.2
            }
        
        # Initialize model
        model = TransformerModel(
            ntoken=len(vocab),
            d_model=getattr(model_args, 'd_model', 512),
            nhead=getattr(model_args, 'nhead', 8),
            d_hid=getattr(model_args, 'd_hid', 512),
            nlayers=getattr(model_args, 'nlayers', 12),
            vocab=vocab,
            dropout=getattr(model_args, 'dropout', 0.2),
            pad_token=vocab[vocab.pad_token],
            pad_value=0,
        )
        
        model.load_state_dict(ckpt['model_state_dict'])
        model.eval()
        print("Model initialized and loaded")
        
    except Exception as e:
        print(f"ERROR loading scGPT model: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Prepare data for perturbation
    print("\nPreparing data for perturbation...")
    
    try:
        # Select subset of cells (matching Geneformer analysis size)
        n_cells = min(500, adata.shape[0])
        np.random.seed(42)  # For reproducibility
        cell_indices = np.random.choice(adata.shape[0], n_cells, replace=False)
        adata_subset = adata[cell_indices, :].copy()
        
        print(f"Selected {n_cells} cells for analysis")
        
        # Ensure we have count data
        if hasattr(adata_subset.X, 'todense'):
            X = np.array(adata_subset.X.todense())
        else:
            X = np.array(adata_subset.X)
            
        # Make sure it's integer counts for scGPT
        X = X.astype(int)
        print(f"Data shape: {X.shape}")
        
    except Exception as e:
        print(f"ERROR preparing data: {e}")
        return
    
    # Simple perturbation approach
    print("\nRunning simplified perturbation analysis...")
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    model.to(device)
    
    perturbation_results = {}
    
    for gene_symbol, vocab_key in genes_in_vocab:
        print(f"\nPerturbing {gene_symbol} (vocab key: {vocab_key})...")
        
        try:
            # Get gene position in dataset
            ensembl_id = gene_mapping[gene_symbol]
            gene_idx = list(adata_subset.var_names).index(ensembl_id)
            
            # Find cells expressing this gene
            expressing_cells = np.where(X[:, gene_idx] > 0)[0]
            print(f"Found {len(expressing_cells)} cells expressing {gene_symbol}")
            
            if len(expressing_cells) == 0:
                print(f"No expressing cells for {gene_symbol}")
                continue
            
            # Take sample of expressing cells for computation efficiency
            max_cells = 25
            if len(expressing_cells) > max_cells:
                sample_cells = np.random.choice(expressing_cells, max_cells, replace=False)
            else:
                sample_cells = expressing_cells
            
            print(f"Analyzing {len(sample_cells)} cells")
            
            shifts = []
            
            for cell_idx in sample_cells:
                try:
                    # Get original expression
                    original_expr = X[cell_idx, :].copy()
                    
                    # Create perturbed version (knockout)
                    perturbed_expr = original_expr.copy()
                    perturbed_expr[gene_idx] = 0
                    
                    # Convert to simple token lists (top expressed genes)
                    # Sort by expression and take top genes
                    n_top_genes = min(500, np.sum(original_expr > 0))
                    
                    # Original tokens
                    original_nonzero = np.where(original_expr > 0)[0]
                    original_values = original_expr[original_nonzero]
                    original_sorted_idx = np.argsort(original_values)[::-1][:n_top_genes]
                    original_genes = original_nonzero[original_sorted_idx]
                    
                    # Perturbed tokens
                    perturbed_nonzero = np.where(perturbed_expr > 0)[0]
                    perturbed_values = perturbed_expr[perturbed_nonzero]
                    perturbed_sorted_idx = np.argsort(perturbed_values)[::-1][:n_top_genes]
                    perturbed_genes = perturbed_nonzero[perturbed_sorted_idx]
                    
                    # Convert gene indices to vocab tokens
                    original_tokens = []
                    perturbed_tokens = []
                    
                    for g_idx in original_genes:
                        gene_name = adata_subset.var['feature_name'].iloc[g_idx]
                        ensembl_name = adata_subset.var_names[g_idx]
                        
                        # Try both gene symbol and Ensembl ID
                        if gene_name in vocab:
                            original_tokens.append(vocab[gene_name])
                        elif ensembl_name in vocab:
                            original_tokens.append(vocab[ensembl_name])
                    
                    for g_idx in perturbed_genes:
                        gene_name = adata_subset.var['feature_name'].iloc[g_idx]
                        ensembl_name = adata_subset.var_names[g_idx]
                        
                        # Try both gene symbol and Ensembl ID
                        if gene_name in vocab:
                            perturbed_tokens.append(vocab[gene_name])
                        elif ensembl_name in vocab:
                            perturbed_tokens.append(vocab[ensembl_name])
                    
                    # Pad to same length
                    max_len = max(len(original_tokens), len(perturbed_tokens), 50)
                    max_len = min(max_len, 1000)  # Cap at reasonable length
                    
                    pad_token_id = vocab[vocab.pad_token]
                    
                    while len(original_tokens) < max_len:
                        original_tokens.append(pad_token_id)
                    while len(perturbed_tokens) < max_len:
                        perturbed_tokens.append(pad_token_id)
                    
                    original_tokens = original_tokens[:max_len]
                    perturbed_tokens = perturbed_tokens[:max_len]
                    
                    # Convert to tensors
                    original_tensor = torch.tensor([original_tokens], dtype=torch.long).to(device)
                    perturbed_tensor = torch.tensor([perturbed_tokens], dtype=torch.long).to(device)
                    
                    # Get embeddings
                    with torch.no_grad():
                        # Simple forward pass
                        original_output = model(original_tensor)
                        perturbed_output = model(perturbed_tensor)
                        
                        # Extract cell embeddings
                        if hasattr(original_output, 'cell_emb'):
                            original_emb = original_output.cell_emb[0].cpu().numpy()
                            perturbed_emb = perturbed_output.cell_emb[0].cpu().numpy()
                        else:
                            # Use pooled representation if cell_emb not available
                            original_emb = original_output[0].mean(dim=0).cpu().numpy()
                            perturbed_emb = perturbed_output[0].mean(dim=0).cpu().numpy()
                        
                        # Calculate embedding shift
                        cosine_sim = 1 - cosine(original_emb, perturbed_emb)
                        embedding_shift = 1 - cosine_sim
                        shifts.append(embedding_shift)
                
                except Exception as e:
                    print(f"Error processing cell {cell_idx}: {e}")
                    continue
            
            # Store results
            if len(shifts) > 0:
                perturbation_results[gene_symbol] = {
                    'embedding_shifts': shifts,
                    'mean_shift': np.mean(shifts),
                    'std_shift': np.std(shifts),
                    'n_cells': len(shifts)
                }
                
                print(f"{gene_symbol}: mean_shift={np.mean(shifts):.6f} ± {np.std(shifts):.6f} (n={len(shifts)})")
            else:
                print(f"No valid results for {gene_symbol}")
        
        except Exception as e:
            print(f"ERROR processing {gene_symbol}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save and display results
    if perturbation_results:
        # Save detailed results
        output_file = os.path.join(output_dir, "scgpt_perturbation_results.json")
        os.makedirs(output_dir, exist_ok=True)
        
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
        
        # Create ranking
        sorted_genes = sorted(perturbation_results.items(), 
                            key=lambda x: x[1]['mean_shift'], 
                            reverse=True)
        
        print("\n" + "="*60)
        print("scGPT PERTURBATION RESULTS")
        print("="*60)
        
        ranking_data = []
        for rank, (gene, data) in enumerate(sorted_genes, 1):
            print(f"{rank}. {gene}: {data['mean_shift']:.6f} ± {data['std_shift']:.6f} (n={data['n_cells']})")
            ranking_data.append({
                'rank': rank,
                'gene': gene,
                'mean_shift': data['mean_shift'],
                'std_shift': data['std_shift'],
                'n_cells': data['n_cells']
            })
        
        # Save ranking
        ranking_df = pd.DataFrame(ranking_data)
        ranking_file = os.path.join(output_dir, "scgpt_gene_ranking.csv")
        ranking_df.to_csv(ranking_file, index=False)
        
        print(f"\nResults saved to:")
        print(f"  - Detailed: {output_file}")
        print(f"  - Ranking: {ranking_file}")
        
        return ranking_data
    
    else:
        print("No successful perturbation results!")
        return None

if __name__ == "__main__":
    main()