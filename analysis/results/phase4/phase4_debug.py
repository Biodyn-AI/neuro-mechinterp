#!/usr/bin/env python3
"""
Phase 4 Debug: Diagnose embedding perturbation issues
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.spatial.distance import cosine
from tqdm import tqdm

# Add Geneformer to path
GENEFORMER_PATH = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
sys.path.insert(0, GENEFORMER_PATH)

from transformers import BertModel, BertConfig

def test_basic_functionality():
    """Test basic components step by step"""
    print("="*60)
    print("Phase 4 Diagnostic Test")
    print("="*60)
    
    # 1. Test GPU availability
    print(f"[GPU] CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"[GPU] Device: {torch.cuda.get_device_name()}")
        print(f"[GPU] Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f}GB")
    
    # 2. Test token dictionary loading
    print("\n[TOKENS] Loading token dictionary...")
    token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'gene_dictionaries_30m', 'token_dictionary_gc30M.pkl')
    if not os.path.exists(token_dict_path):
        token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
    
    try:
        with open(token_dict_path, 'rb') as f:
            token_dictionary = pickle.load(f)
        print(f"[TOKENS] ✓ Loaded {len(token_dictionary)} tokens")
        
        # Check if target genes are in dictionary
        target_genes = ['CADM2', 'NRXN1', 'NLGN1', 'NEGR1', 'TCF4']
        found_genes = [g for g in target_genes if g in token_dictionary]
        print(f"[TOKENS] Target genes found: {len(found_genes)}/{len(target_genes)}")
        print(f"[TOKENS] Found genes: {found_genes}")
        
    except Exception as e:
        print(f"[TOKENS] ✗ Failed: {e}")
        return False
    
    # 3. Test model loading
    print("\n[MODEL] Loading Geneformer V1-10M...")
    model_path = r"D:\openclaw\intelligence-augmentation\models\Geneformer\Geneformer-V1-10M"
    
    try:
        model = BertModel.from_pretrained(model_path)
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model.to(device)
        model.eval()
        print(f"[MODEL] ✓ Loaded model on {device}")
        
        # Test inference with dummy input
        dummy_input = torch.tensor([[0, 100, 200, 1]], dtype=torch.long).to(device)
        with torch.no_grad():
            outputs = model(dummy_input)
            embedding = outputs.last_hidden_state[0, 0, :].cpu().numpy()
        print(f"[MODEL] ✓ Test inference: {embedding.shape} embedding")
        
    except Exception as e:
        print(f"[MODEL] ✗ Failed: {e}")
        return False
    
    # 4. Test data loading
    print("\n[DATA] Loading brain data...")
    data_path = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
    
    try:
        adata = sc.read_h5ad(data_path)
        print(f"[DATA] ✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Check gene names
        if 'feature_name' in adata.var.columns:
            gene_symbols = adata.var['feature_name'].values
            print(f"[DATA] ✓ Using feature_name column for gene symbols")
        else:
            gene_symbols = adata.var.index.values
            print(f"[DATA] ✓ Using index for gene symbols")
            
        # Test single cell processing
        cell_0 = adata[0]
        cell_expr = dict(zip(gene_symbols, cell_0.X.toarray().flatten()))
        
        # Filter positive expression
        pos_genes = {g: e for g, e in cell_expr.items() if e > 0}
        print(f"[DATA] ✓ Cell 0 has {len(pos_genes)} genes with positive expression")
        
        # Check overlap with token dictionary
        valid_genes = {g: e for g, e in pos_genes.items() if g in token_dictionary}
        print(f"[DATA] ✓ Cell 0 has {len(valid_genes)} genes in token dictionary")
        
        if len(valid_genes) < 10:
            print(f"[DATA] ⚠ Warning: Very few valid genes. Showing first 10 positive genes:")
            for i, (g, e) in enumerate(list(pos_genes.items())[:10]):
                in_dict = "✓" if g in token_dictionary else "✗"
                print(f"[DATA]     {g}: {e:.3f} {in_dict}")
        
    except Exception as e:
        print(f"[DATA] ✗ Failed: {e}")
        return False
    
    # 5. Test tokenization and embedding
    print("\n[TOKENIZE] Testing cell tokenization...")
    
    try:
        # Rank genes by expression
        ranked_genes = sorted(valid_genes.items(), key=lambda x: x[1], reverse=True)
        print(f"[TOKENIZE] ✓ Ranked {len(ranked_genes)} genes by expression")
        
        # Convert to tokens (top 100 for test)
        token_ids = [token_dictionary[gene] for gene, _ in ranked_genes[:100]]
        token_ids = [0] + token_ids + [1]  # Add special tokens
        print(f"[TOKENIZE] ✓ Created token sequence: {len(token_ids)} tokens")
        
        # Test embedding
        input_ids = torch.tensor(token_ids, dtype=torch.long).unsqueeze(0).to(device)
        with torch.no_grad():
            outputs = model(input_ids)
            baseline_embedding = outputs.last_hidden_state[0, 0, :].cpu().numpy()
        
        print(f"[EMBED] ✓ Baseline embedding: {baseline_embedding.shape}")
        
        # Test perturbation
        if len(token_ids) > 3:  # Need at least [CLS] + 1 gene + [SEP]
            # Remove second token (first gene)
            perturbed_tokens = token_ids[:1] + token_ids[2:]  
            perturbed_input = torch.tensor(perturbed_tokens, dtype=torch.long).unsqueeze(0).to(device)
            
            with torch.no_grad():
                perturbed_outputs = model(perturbed_input)
                perturbed_embedding = perturbed_outputs.last_hidden_state[0, 0, :].cpu().numpy()
            
            # Calculate cosine distance
            cos_dist = cosine(baseline_embedding, perturbed_embedding)
            print(f"[PERTURB] ✓ Perturbation test: cosine distance = {cos_dist:.6f}")
            
            if cos_dist > 0:
                print("[PERTURB] ✓ Non-zero distance indicates perturbation is working!")
            else:
                print("[PERTURB] ⚠ Zero distance - check if perturbation is effective")
        
    except Exception as e:
        print(f"[TOKENIZE/EMBED] ✗ Failed: {e}")
        return False
    
    print("\n" + "="*60)
    print("✓ All diagnostic tests passed!")
    print("="*60)
    return True

def test_target_genes():
    """Test specific target genes in the data"""
    print("\n[TARGET GENES] Testing target gene presence...")
    
    # Load data
    data_path = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
    adata = sc.read_h5ad(data_path)
    
    # Sample small subset
    sc.pp.subsample(adata, n_obs=10, random_state=42)
    
    # Get gene names
    if 'feature_name' in adata.var.columns:
        gene_symbols = adata.var['feature_name'].values
    else:
        gene_symbols = adata.var.index.values
    
    # Load token dictionary
    token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'gene_dictionaries_30m', 'token_dictionary_gc30M.pkl')
    if not os.path.exists(token_dict_path):
        token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
    
    with open(token_dict_path, 'rb') as f:
        token_dictionary = pickle.load(f)
    
    target_genes = ['CADM2', 'NRXN1', 'NLGN1', 'NEGR1', 'TCF4', 'MEF2C', 'APP', 'GRIN2A']
    
    print(f"Testing {len(target_genes)} target genes across {adata.n_obs} cells:")
    print("-" * 60)
    
    results = []
    for gene in target_genes:
        if gene not in token_dictionary:
            print(f"{gene:8s} | NOT IN TOKEN DICT")
            continue
            
        # Count cells expressing this gene
        gene_idx = np.where(gene_symbols == gene)[0]
        if len(gene_idx) == 0:
            print(f"{gene:8s} | NOT IN DATA")
            continue
            
        gene_idx = gene_idx[0]
        expr_values = adata[:, gene_idx].X.toarray().flatten()
        expressing_cells = np.sum(expr_values > 0)
        mean_expr = np.mean(expr_values[expr_values > 0]) if expressing_cells > 0 else 0
        
        print(f"{gene:8s} | {expressing_cells:2d}/{adata.n_obs} cells | Mean expr: {mean_expr:.3f}")
        results.append({
            'gene': gene,
            'expressing_cells': expressing_cells,
            'mean_expr': mean_expr
        })
    
    return results

if __name__ == "__main__":
    # Run diagnostic tests
    success = test_basic_functionality()
    
    if success:
        # Test target genes
        gene_results = test_target_genes()
        
        print(f"\n[SUMMARY] Target gene analysis:")
        expressing_genes = [r for r in gene_results if r['expressing_cells'] > 0]
        print(f"  • {len(expressing_genes)}/{len(gene_results)} target genes found with expression")
        
        if len(expressing_genes) > 0:
            print("\n[READY] System appears ready for Phase 4 analysis!")
        else:
            print("\n[ISSUE] No target genes found with expression - check gene name mapping!")
    else:
        print("\n[FAILED] Diagnostic tests failed - fix issues before running Phase 4")