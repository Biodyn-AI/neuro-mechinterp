#!/usr/bin/env python3
"""
Simple Phase 4 test - no unicode
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import torch

# Add Geneformer to path
GENEFORMER_PATH = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
sys.path.insert(0, GENEFORMER_PATH)

from transformers import BertModel

def main():
    print("Phase 4 Simple Test")
    print("=" * 40)
    
    # Test 1: GPU
    print("GPU available:", torch.cuda.is_available())
    
    # Test 2: Token dictionary
    token_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'gene_dictionaries_30m', 'token_dictionary_gc30M.pkl')
    if not os.path.exists(token_path):
        token_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
    
    with open(token_path, 'rb') as f:
        token_dict = pickle.load(f)
    print("Tokens loaded:", len(token_dict))
    
    # Test 3: Target genes in dictionary
    targets = ['CADM2', 'NRXN1', 'NLGN1', 'TCF4', 'APP']
    found = [g for g in targets if g in token_dict]
    print("Target genes found:", found)
    
    # Test 4: Model loading
    model_path = r"D:\openclaw\intelligence-augmentation\models\Geneformer\Geneformer-V1-10M"
    model = BertModel.from_pretrained(model_path)
    device = torch.device('cuda')
    model.to(device)
    model.eval()
    print("Model loaded successfully")
    
    # Test 5: Data loading
    data_path = r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad"
    adata = sc.read_h5ad(data_path)
    print("Data shape:", adata.shape)
    
    # Test 6: Gene name mapping
    if 'feature_name' in adata.var.columns:
        gene_symbols = adata.var['feature_name'].values
        print("Using feature_name column")
    else:
        gene_symbols = adata.var.index.values
        print("Using var index")
    
    # Test 7: Single cell analysis
    cell_0 = adata[0]
    cell_expr = dict(zip(gene_symbols, cell_0.X.toarray().flatten()))
    pos_genes = {g: e for g, e in cell_expr.items() if e > 0}
    valid_genes = {g: e for g, e in pos_genes.items() if g in token_dict}
    
    print("Cell 0 stats:")
    print("  Total genes:", len(cell_expr))
    print("  Positive expression:", len(pos_genes))
    print("  Valid tokens:", len(valid_genes))
    
    # Test 8: Check target gene expression
    print("\nTarget gene expression in first 10 cells:")
    for gene in targets[:3]:  # Test first 3
        if gene in token_dict:
            gene_mask = gene_symbols == gene
            if np.any(gene_mask):
                gene_idx = np.where(gene_mask)[0][0]
                expr_vals = adata[:10, gene_idx].X.toarray().flatten()
                expressing = np.sum(expr_vals > 0)
                print(f"  {gene}: {expressing}/10 cells expressing")
            else:
                print(f"  {gene}: not found in data")
    
    # Test 9: Tokenization test
    if len(valid_genes) > 0:
        ranked = sorted(valid_genes.items(), key=lambda x: x[1], reverse=True)
        token_ids = [token_dict[g] for g, _ in ranked[:50]]
        token_ids = [0] + token_ids + [1]  # Add CLS/SEP
        
        print(f"\nTokenization test: {len(token_ids)} tokens")
        
        # Test embedding
        input_tensor = torch.tensor(token_ids).unsqueeze(0).to(device)
        with torch.no_grad():
            outputs = model(input_tensor)
            embedding = outputs.last_hidden_state[0, 0, :].cpu().numpy()
        
        print(f"Embedding shape: {embedding.shape}")
        print(f"Embedding norm: {np.linalg.norm(embedding):.3f}")
        
        # Test perturbation
        if len(token_ids) > 3:
            perturbed = token_ids[:1] + token_ids[2:]  # Remove first gene
            perturbed_tensor = torch.tensor(perturbed).unsqueeze(0).to(device)
            
            with torch.no_grad():
                perturbed_outputs = model(perturbed_tensor)
                perturbed_emb = perturbed_outputs.last_hidden_state[0, 0, :].cpu().numpy()
            
            # Cosine distance
            from scipy.spatial.distance import cosine
            distance = cosine(embedding, perturbed_emb)
            print(f"Perturbation distance: {distance:.6f}")
            
            if distance > 0.001:
                print("SUCCESS: Perturbation produces measurable distance!")
            else:
                print("WARNING: Very small perturbation distance")
    
    print("\nDiagnostic complete!")

if __name__ == "__main__":
    main()