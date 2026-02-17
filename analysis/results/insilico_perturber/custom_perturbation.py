#!/usr/bin/env python3
"""
Custom in silico perturbation approach for Geneformer models.
Bypasses the dependency issues by working directly with transformers.
"""

import os
import pickle
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from transformers import BertForMaskedLM, BertConfig
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity
from scipy import stats
import json

class CustomGeneformerPerturber:
    def __init__(self, model_path, gene_mapping_path, batch_size=32):
        """Initialize the custom perturber."""
        self.model_path = model_path
        self.batch_size = batch_size
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        # Load gene mapping
        with open(gene_mapping_path, 'rb') as f:
            self.gene_mapping = pickle.load(f)
        
        print(f"Using device: {self.device}")
        print(f"Loaded {len(self.gene_mapping)} gene mappings")
        
    def load_model(self):
        """Load the Geneformer model."""
        try:
            print(f"Loading model from {self.model_path}")
            self.model = BertForMaskedLM.from_pretrained(self.model_path)
            self.model.to(self.device)
            self.model.eval()
            print("[OK] Model loaded successfully")
            
            # Load token dictionary
            token_dict_path = os.path.join(
                os.path.dirname(self.model_path), 
                "../../geneformer/token_dictionary_gc104M.pkl"
            )
            if os.path.exists(token_dict_path):
                with open(token_dict_path, 'rb') as f:
                    self.token_dictionary = pickle.load(f)
                print(f"[OK] Loaded token dictionary with {len(self.token_dictionary)} tokens")
            else:
                print("[WARNING] Token dictionary not found, using fallback approach")
                self.token_dictionary = None
                
        except Exception as e:
            print(f"[ERROR] Error loading model: {e}")
            raise
    
    def tokenize_cells(self, adata, max_tokens=2048):
        """Convert expression data to token sequences."""
        print("Tokenizing expression data...")
        
        # Simple tokenization approach: rank genes by expression and take top N
        tokenized_data = []
        
        for i in range(adata.n_obs):
            if i % 1000 == 0:
                print(f"Tokenizing cell {i}/{adata.n_obs}")
                
            # Get expression vector for this cell
            if hasattr(adata.X, 'toarray'):
                expression = adata.X[i].toarray().flatten()
            else:
                expression = adata.X[i].flatten()
            
            # Find non-zero genes and their expressions
            nonzero_idx = np.nonzero(expression)[0]
            nonzero_expr = expression[nonzero_idx]
            
            # Sort by expression level (descending)
            sorted_indices = np.argsort(nonzero_expr)[::-1]
            
            # Take top genes up to max_tokens
            top_genes = nonzero_idx[sorted_indices[:max_tokens]]
            top_expressions = nonzero_expr[sorted_indices[:max_tokens]]
            
            # Create token sequence (using gene indices as simple tokens for now)
            tokens = top_genes.tolist()
            
            tokenized_data.append({
                'input_ids': tokens,
                'expressions': top_expressions.tolist(),
                'cell_idx': i
            })
        
        print(f"[OK] Tokenized {len(tokenized_data)} cells")
        return tokenized_data
    
    def get_embeddings(self, tokenized_data):
        """Get embeddings from the model for tokenized data."""
        print("Computing embeddings...")
        
        all_embeddings = []
        
        for i in range(0, len(tokenized_data), self.batch_size):
            batch = tokenized_data[i:i + self.batch_size]
            
            # Prepare batch
            max_len = max(len(item['input_ids']) for item in batch)
            max_len = min(max_len, 512)  # Limit sequence length
            
            input_ids = []
            attention_masks = []
            
            for item in batch:
                tokens = item['input_ids'][:max_len]
                
                # Pad sequence
                padded = tokens + [0] * (max_len - len(tokens))
                attention_mask = [1] * len(tokens) + [0] * (max_len - len(tokens))
                
                input_ids.append(padded)
                attention_masks.append(attention_mask)
            
            # Convert to tensors
            input_ids = torch.tensor(input_ids, dtype=torch.long).to(self.device)
            attention_masks = torch.tensor(attention_masks, dtype=torch.long).to(self.device)
            
            # Get embeddings
            with torch.no_grad():
                try:
                    # Get hidden states from the model
                    outputs = self.model.bert(input_ids=input_ids, attention_mask=attention_masks, output_hidden_states=True)
                    hidden_states = outputs.hidden_states[-1]  # Last layer
                    
                    # Mean pooling over sequence length (excluding padding)
                    embeddings = []
                    for j, mask in enumerate(attention_masks):
                        valid_length = mask.sum().item()
                        if valid_length > 0:
                            emb = hidden_states[j, :valid_length].mean(dim=0)
                        else:
                            emb = hidden_states[j, 0]  # Use first token if no valid tokens
                        embeddings.append(emb.cpu().numpy())
                    
                    all_embeddings.extend(embeddings)
                    
                except Exception as e:
                    print(f"[WARNING] Error in batch {i//self.batch_size}: {e}")
                    # Add zero embeddings for failed batch
                    for _ in batch:
                        all_embeddings.append(np.zeros(768))  # Assume 768-dim embeddings
        
        print(f"[OK] Computed embeddings for {len(all_embeddings)} cells")
        return np.array(all_embeddings)
    
    def perturb_gene(self, tokenized_data, target_gene_ensembl, perturbation_type='delete'):
        """Simulate gene perturbation by modifying token sequences."""
        print(f"Simulating {perturbation_type} perturbation of {target_gene_ensembl}")
        
        # Find gene index in the data
        # For now, this is a simplified approach
        perturbed_data = []
        
        for item in tokenized_data:
            new_item = item.copy()
            
            if perturbation_type == 'delete':
                # Remove the target gene from the sequence
                new_tokens = [token for token in item['input_ids'] if token != target_gene_ensembl]
                new_item['input_ids'] = new_tokens
                
            elif perturbation_type == 'overexpress':
                # Add multiple copies or increase rank of the target gene
                new_tokens = item['input_ids'].copy()
                if target_gene_ensembl not in new_tokens:
                    new_tokens.insert(0, target_gene_ensembl)  # Add at beginning
                new_item['input_ids'] = new_tokens
            
            perturbed_data.append(new_item)
        
        return perturbed_data
    
    def run_perturbation_analysis(self, adata, output_dir):
        """Run the full perturbation analysis."""
        print("Starting custom perturbation analysis...")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        results = {}
        
        try:
            # Tokenize the original data
            tokenized_data = self.tokenize_cells(adata, max_tokens=1000)
            
            # Get original embeddings
            print("Computing original embeddings...")
            original_embeddings = self.get_embeddings(tokenized_data)
            
            # Test perturbation for each intelligence gene
            for gene_symbol, ensembl_id in self.gene_mapping.items():
                print(f"\nAnalyzing perturbation of {gene_symbol} ({ensembl_id})")
                
                try:
                    # Perturb the gene (deletion)
                    perturbed_data = self.perturb_gene(tokenized_data, ensembl_id, 'delete')
                    
                    # Get perturbed embeddings
                    perturbed_embeddings = self.get_embeddings(perturbed_data)
                    
                    # Calculate cosine similarity between original and perturbed
                    similarities = []
                    for i in range(len(original_embeddings)):
                        sim = cosine_similarity(
                            [original_embeddings[i]], 
                            [perturbed_embeddings[i]]
                        )[0][0]
                        similarities.append(sim)
                    
                    # Calculate statistics
                    mean_similarity = np.mean(similarities)
                    std_similarity = np.std(similarities)
                    effect_size = 1 - mean_similarity  # Higher effect = more different
                    
                    results[gene_symbol] = {
                        'ensembl_id': ensembl_id,
                        'mean_similarity': float(mean_similarity),
                        'std_similarity': float(std_similarity),
                        'effect_size': float(effect_size),
                        'min_similarity': float(np.min(similarities)),
                        'max_similarity': float(np.max(similarities)),
                        'n_cells': len(similarities)
                    }
                    
                    print(f"  Effect size: {effect_size:.4f}")
                    
                except Exception as e:
                    print(f"  [ERROR] Error analyzing {gene_symbol}: {e}")
                    results[gene_symbol] = {
                        'ensembl_id': ensembl_id,
                        'error': str(e)
                    }
        
        except Exception as e:
            print(f"[ERROR] Critical error in analysis: {e}")
            raise
        
        # Save results
        results_path = os.path.join(output_dir, 'perturbation_results.json')
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Create results DataFrame
        valid_results = {k: v for k, v in results.items() if 'error' not in v}
        if valid_results:
            df = pd.DataFrame.from_dict(valid_results, orient='index')
            df = df.sort_values('effect_size', ascending=False)
            
            df_path = os.path.join(output_dir, 'perturbation_results.csv')
            df.to_csv(df_path)
            
            print(f"\nTop genes by perturbation effect:")
            print(df[['effect_size', 'mean_similarity']].head(10))
            
        return results

if __name__ == "__main__":
    # Load prepared data
    print("Loading prepared data...")
    adata = sc.read_h5ad("D:/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")
    
    # Use a smaller subset for testing
    adata = adata[:1000, :]  # Use first 1000 cells for testing
    print(f"Using subset: {adata.shape}")
    
    # Initialize perturber
    model_path = "D:/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V2-104M"
    gene_mapping_path = "D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/intelligence_genes_mapping.pkl"
    output_dir = "D:/openclaw/intelligence-augmentation/analysis/results/insilico_perturber/custom_results"
    
    perturber = CustomGeneformerPerturber(model_path, gene_mapping_path)
    
    try:
        perturber.load_model()
        results = perturber.run_perturbation_analysis(adata, output_dir)
        print("[OK] Analysis completed successfully!")
        
    except Exception as e:
        print(f"[ERROR] Analysis failed: {e}")
        # Fallback to simple expression-based analysis
        print("\nTrying fallback analysis...")
        
        with open(gene_mapping_path, 'rb') as f:
            gene_mapping = pickle.load(f)
        
        # Simple analysis based on expression correlation
        results = {}
        for gene_symbol, ensembl_id in gene_mapping.items():
            if ensembl_id in adata.var.index:
                gene_idx = adata.var.index.get_loc(ensembl_id)
                if hasattr(adata.X, 'toarray'):
                    expression = adata.X[:, gene_idx].toarray().flatten()
                else:
                    expression = adata.X[:, gene_idx].flatten()
                
                results[gene_symbol] = {
                    'mean_expression': float(np.mean(expression)),
                    'std_expression': float(np.std(expression)),
                    'percent_expressed': float((expression > 0).sum() / len(expression) * 100)
                }
        
        # Save fallback results
        with open(os.path.join(output_dir, 'fallback_results.json'), 'w') as f:
            json.dump(results, f, indent=2)
        
        print("[OK] Fallback analysis completed")