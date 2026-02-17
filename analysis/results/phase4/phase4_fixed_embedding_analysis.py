#!/usr/bin/env python3
"""
Phase 4: FIXED Embedding-Based Perturbation Analysis
==================================================

FIXED: Now properly maps gene symbols to Ensembl IDs for Geneformer tokens
"""

import os
import sys
import gc
import json
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F
from scipy.spatial.distance import cosine
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Add Geneformer to path
GENEFORMER_PATH = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
sys.path.insert(0, GENEFORMER_PATH)

from transformers import BertModel, BertConfig

# Configuration
CONFIG = {
    'model_path': r"D:\openclaw\intelligence-augmentation\models\Geneformer\Geneformer-V1-10M",
    'data_path': r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad",
    'output_dir': r"D:\openclaw\intelligence-augmentation\analysis\results\phase4",
    'max_cells': 200,  # Pilot run
    'batch_size': 50,  # GPU memory constraint  
    'max_tokens': 2048,  # Geneformer context length
    'device': 'cuda' if torch.cuda.is_available() else 'cpu',
    'n_permutations': 50,  # Reduced for speed
    'target_genes': [
        'CADM2', 'NRXN1', 'NLGN1', 'NEGR1', 'TCF4', 'MEF2C', 'APP', 'GRIN2A', 'GRIN2B',
        'MAPT', 'SCN1A', 'DISC1', 'FOXO3', 'CAMK2A', 'HOMER1', 'CREB1', 'COMT', 'SHANK3',
        'FMR1', 'SYN1', 'BDNF', 'SLC6A4'
    ]
}

class FixedEmbeddingPerturbationAnalyzer:
    """
    FIXED: Proper embedding-based perturbation with gene symbol->Ensembl mapping
    """
    
    def __init__(self, config):
        self.config = config
        self.device = torch.device(config['device'])
        self.results = []
        
        # Load token dictionary
        self.load_token_dictionary()
        
        # Load and map gene symbols
        self.create_gene_mapping()
        
        # Load model
        self.load_model()
        
        print(f"[INIT] Initialized on {self.device}")
        
    def load_token_dictionary(self):
        """Load Geneformer token dictionary (Ensembl IDs)"""
        token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'gene_dictionaries_30m', 'token_dictionary_gc30M.pkl')
        if not os.path.exists(token_dict_path):
            token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
            
        with open(token_dict_path, 'rb') as f:
            self.token_dictionary = pickle.load(f)
        print(f"[TOKENS] Loaded token dictionary: {len(self.token_dictionary)} tokens")
        
        # Create reverse mapping
        self.id_to_token = {v: k for k, v in self.token_dictionary.items()}
        
    def create_gene_mapping(self):
        """Create mapping from gene symbols to Ensembl IDs"""
        print("[MAPPING] Creating gene symbol -> Ensembl ID mapping...")
        
        # Load data to get gene mapping
        adata = sc.read_h5ad(self.config['data_path'])
        
        # Create mapping: gene_symbol -> ensembl_id
        self.symbol_to_ensembl = {}
        
        # Use var.index (Ensembl IDs) and feature_name (gene symbols)
        if 'feature_name' in adata.var.columns:
            for ensembl_id, gene_symbol in zip(adata.var.index, adata.var['feature_name']):
                if pd.notna(gene_symbol) and gene_symbol.strip():
                    self.symbol_to_ensembl[gene_symbol] = ensembl_id
        
        print(f"[MAPPING] Created mapping for {len(self.symbol_to_ensembl)} genes")
        
        # Map target genes to Ensembl IDs
        self.target_ensembl = {}
        missing_genes = []
        
        for gene_symbol in self.config['target_genes']:
            if gene_symbol in self.symbol_to_ensembl:
                ensembl_id = self.symbol_to_ensembl[gene_symbol]
                if ensembl_id in self.token_dictionary:
                    self.target_ensembl[gene_symbol] = ensembl_id
                    print(f"[MAPPING] {gene_symbol} -> {ensembl_id} (token: {self.token_dictionary[ensembl_id]})")
                else:
                    missing_genes.append(f"{gene_symbol} (Ensembl {ensembl_id} not in tokens)")
            else:
                missing_genes.append(f"{gene_symbol} (not in data)")
        
        print(f"[MAPPING] Successfully mapped {len(self.target_ensembl)}/{len(self.config['target_genes'])} target genes")
        
        if missing_genes:
            print(f"[MAPPING] Missing genes: {missing_genes}")
        
    def load_model(self):
        """Load Geneformer V1-10M model"""
        try:
            # Load model
            self.model = BertModel.from_pretrained(self.config['model_path'])
            self.model.to(self.device)
            self.model.eval()
            
            # Get config for embedding dimension
            config = self.model.config
            print(f"[MODEL] Loaded Geneformer V1-10M ({config.hidden_size}D embeddings)")
            
        except Exception as e:
            print(f"[ERROR] Failed to load model: {e}")
            raise
            
    def tokenize_cell(self, cell_data, max_tokens=2048):
        """
        Convert single cell to Geneformer tokens using Ensembl IDs
        
        Args:
            cell_data: AnnData cell with .X expression and .var metadata
            max_tokens: Maximum context length
            
        Returns:
            List of token IDs ranked by expression
        """
        # Get expression values and gene IDs (Ensembl)
        expr_values = cell_data.X.toarray().flatten()
        gene_ids = cell_data.var.index.values  # These are Ensembl IDs
        
        # Create expression dictionary for genes with positive expression
        gene_expr = {}
        for i, (gene_id, expr) in enumerate(zip(gene_ids, expr_values)):
            if expr > 0 and gene_id in self.token_dictionary:
                gene_expr[gene_id] = expr
        
        if not gene_expr:
            return []
            
        # Rank genes by expression (descending)
        ranked_genes = sorted(gene_expr.items(), key=lambda x: x[1], reverse=True)
        
        # Convert to token IDs
        token_ids = [self.token_dictionary[gene_id] for gene_id, _ in ranked_genes[:max_tokens-2]]
        
        # Add special tokens: [CLS] = 0, [SEP] = 1
        token_ids = [0] + token_ids + [1]
        
        return token_ids, gene_expr
        
    def get_embedding(self, token_ids):
        """Get embedding for tokenized cell"""
        if not token_ids or len(token_ids) < 2:
            return None
            
        # Convert to tensor
        input_ids = torch.tensor(token_ids, dtype=torch.long).unsqueeze(0).to(self.device)
        
        with torch.no_grad():
            outputs = self.model(input_ids)
            hidden_states = outputs.last_hidden_state
            
            # Extract [CLS] token embedding (first token)
            cls_embedding = hidden_states[0, 0, :].cpu().numpy()
            
        return cls_embedding
        
    def perturb_and_measure(self, token_ids, gene_expr, target_gene_symbol, n_permutations=50):
        """Perturb cell by removing target gene and measure embedding distance"""
        
        if target_gene_symbol not in self.target_ensembl:
            return None
            
        target_ensembl = self.target_ensembl[target_gene_symbol]
        target_token_id = self.token_dictionary[target_ensembl]
        
        # Check if target gene is in this cell
        if target_ensembl not in gene_expr:
            return {
                'target_gene': target_gene_symbol,
                'target_ensembl': target_ensembl,
                'target_in_cell': False,
                'cosine_distance': np.nan,
                'null_distances': [],
                'p_value': np.nan
            }
            
        # Get baseline embedding
        baseline_embedding = self.get_embedding(token_ids)
        if baseline_embedding is None:
            return None
            
        # Remove target gene token
        perturbed_tokens = [t for t in token_ids if t != target_token_id]
        
        # Get perturbed embedding
        perturbed_embedding = self.get_embedding(perturbed_tokens)
        if perturbed_embedding is None:
            return None
            
        # Calculate cosine distance
        cosine_dist = cosine(baseline_embedding, perturbed_embedding)
        
        # Generate null distribution by removing random non-target genes
        null_distances = []
        available_tokens = [t for t in token_ids[1:-1] if t != target_token_id]  # Exclude special tokens
        
        if len(available_tokens) >= n_permutations:
            random_tokens = np.random.choice(available_tokens, n_permutations, replace=False)
        else:
            random_tokens = available_tokens * (n_permutations // max(1, len(available_tokens)) + 1)
            random_tokens = random_tokens[:n_permutations] if len(random_tokens) > 0 else []
            
        for random_token in random_tokens:
            # Remove random gene
            null_tokens = [t for t in token_ids if t != random_token]
            null_embedding = self.get_embedding(null_tokens)
            
            if null_embedding is not None:
                null_dist = cosine(baseline_embedding, null_embedding)
                null_distances.append(null_dist)
                
        # Compute empirical p-value
        if null_distances:
            p_value = np.sum(np.array(null_distances) >= cosine_dist) / len(null_distances)
        else:
            p_value = np.nan
            
        return {
            'target_gene': target_gene_symbol,
            'target_ensembl': target_ensembl,
            'target_in_cell': True,
            'cosine_distance': cosine_dist,
            'null_distances': null_distances,
            'p_value': p_value,
            'n_null_samples': len(null_distances),
            'target_expression': gene_expr.get(target_ensembl, 0)
        }
        
    def analyze_cells(self):
        """Main analysis pipeline"""
        print("[DATA] Loading single-cell data...")
        
        # Load data
        adata = sc.read_h5ad(self.config['data_path'])
        print(f"[DATA] Loaded {adata.n_obs} cells x {adata.n_vars} genes")
        
        # Sample cells for pilot
        if adata.n_obs > self.config['max_cells']:
            sc.pp.subsample(adata, n_obs=self.config['max_cells'], random_state=42)
            print(f"[SAMPLE] Subsampled to {adata.n_obs} cells")
            
        print(f"[GENES] Processing {len(self.target_ensembl)} mapped target genes")
        
        # Process in batches
        batch_results = []
        n_batches = (adata.n_obs + self.config['batch_size'] - 1) // self.config['batch_size']
        
        for batch_idx in range(n_batches):
            print(f"\n[BATCH] Processing batch {batch_idx + 1}/{n_batches}")
            
            start_idx = batch_idx * self.config['batch_size']
            end_idx = min(start_idx + self.config['batch_size'], adata.n_obs)
            
            batch_cells = adata[start_idx:end_idx]
            
            for cell_idx in tqdm(range(batch_cells.n_obs), desc="Cells"):
                cell_data = batch_cells[cell_idx]
                
                # Tokenize cell
                token_ids, gene_expr = self.tokenize_cell(cell_data, self.config['max_tokens'])
                
                if not token_ids:
                    continue
                    
                # Get cell type if available
                cell_type = cell_data.obs.get('cell_type', 'Unknown')
                
                # Test each target gene
                for target_gene in self.target_ensembl.keys():
                    result = self.perturb_and_measure(
                        token_ids, 
                        gene_expr,
                        target_gene, 
                        self.config['n_permutations']
                    )
                    
                    if result:
                        result.update({
                            'cell_idx': start_idx + cell_idx,
                            'cell_type': cell_type,
                            'n_tokens': len(token_ids),
                            'batch': batch_idx + 1
                        })
                        batch_results.append(result)
                        
            # Clear GPU memory
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()
            
            # Save intermediate results
            if batch_results:
                batch_df = pd.DataFrame(batch_results)
                batch_file = os.path.join(self.config['output_dir'], f'batch_{batch_idx + 1}_results_fixed.csv')
                batch_df.to_csv(batch_file, index=False)
                print(f"[SAVE] Saved batch {batch_idx + 1}: {len(batch_results)} perturbations")
                
        self.results = batch_results
        
    def aggregate_and_test(self):
        """Aggregate results and perform statistical testing"""
        if not self.results:
            print("[ERROR] No results to aggregate!")
            return pd.DataFrame()
            
        df = pd.DataFrame(self.results)
        self.cell_results = df  # Store cell results
        
        # Filter out cells where target gene wasn't present
        df_valid = df[df['target_in_cell'] == True].copy()
        
        print(f"[STATS] Aggregating {len(df_valid)} valid perturbations")
        
        # Aggregate by gene
        gene_stats = []
        
        for gene_symbol in self.target_ensembl.keys():
            gene_data = df_valid[df_valid['target_gene'] == gene_symbol]
            
            if len(gene_data) == 0:
                continue
                
            # Collect distances
            distances = gene_data['cosine_distance'].dropna()
            
            if len(distances) == 0:
                continue
                
            # Collect null distribution
            all_null_distances = []
            for null_list in gene_data['null_distances']:
                if isinstance(null_list, list):
                    all_null_distances.extend(null_list)
                    
            # Statistical testing
            if all_null_distances and len(distances) > 0:
                # Wilcoxon rank-sum test
                statistic, wilcoxon_p = ranksums(distances, all_null_distances)
                
                # Empirical p-value (average across cells)
                empirical_p = gene_data['p_value'].dropna().mean()
            else:
                statistic, wilcoxon_p, empirical_p = np.nan, np.nan, np.nan
                
            gene_stats.append({
                'gene': gene_symbol,
                'ensembl_id': self.target_ensembl[gene_symbol],
                'n_cells': len(distances),
                'mean_distance': distances.mean(),
                'std_distance': distances.std(),
                'median_distance': distances.median(),
                'max_distance': distances.max(),
                'min_distance': distances.min(),
                'empirical_p_value': empirical_p,
                'wilcoxon_p_value': wilcoxon_p,
                'wilcoxon_statistic': statistic,
                'n_null_samples': len(all_null_distances),
                'mean_expression': gene_data['target_expression'].mean()
            })
            
        # Convert to DataFrame and apply FDR correction
        gene_df = pd.DataFrame(gene_stats)
        
        if len(gene_df) > 0:
            # FDR correction for empirical p-values
            valid_empirical_p = gene_df['empirical_p_value'].dropna()
            if len(valid_empirical_p) > 1:  # Need at least 2 values for FDR
                _, empirical_fdr, _, _ = multipletests(valid_empirical_p, method='fdr_bh')
                gene_df.loc[valid_empirical_p.index, 'empirical_fdr'] = empirical_fdr
                
            # FDR correction for Wilcoxon p-values  
            valid_wilcoxon_p = gene_df['wilcoxon_p_value'].dropna()
            if len(valid_wilcoxon_p) > 1:
                _, wilcoxon_fdr, _, _ = multipletests(valid_wilcoxon_p, method='fdr_bh')
                gene_df.loc[valid_wilcoxon_p.index, 'wilcoxon_fdr'] = wilcoxon_fdr
                
            # Sort by mean distance (descending)
            gene_df = gene_df.sort_values('mean_distance', ascending=False)
            
        self.gene_summary = gene_df
        
        return gene_df
        
    def save_results(self):
        """Save all results to files"""
        # Save cell-level results
        cell_file = os.path.join(self.config['output_dir'], 'cell_perturbation_results_phase4_fixed.csv')
        self.cell_results.to_csv(cell_file, index=False)
        print(f"[SAVE] Saved cell results: {cell_file}")
        
        # Save gene-level summary
        gene_file = os.path.join(self.config['output_dir'], 'gene_summary_phase4_fixed.csv')
        self.gene_summary.to_csv(gene_file, index=False)
        print(f"[SAVE] Saved gene summary: {gene_file}")
        
        # Save analysis metadata
        metadata = {
            'config': self.config,
            'n_cells_processed': len(self.cell_results),
            'n_valid_perturbations': len(self.cell_results[self.cell_results['target_in_cell'] == True]),
            'n_genes_analyzed': len(self.gene_summary),
            'target_gene_mapping': self.target_ensembl,
            'device': str(self.device),
            'model_path': self.config['model_path'],
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        metadata_file = os.path.join(self.config['output_dir'], 'analysis_metadata_phase4_fixed.json')
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
        print(f"[SAVE] Saved metadata: {metadata_file}")
        
    def print_summary(self):
        """Print analysis summary"""
        print("\n" + "="*80)
        print("PHASE 4: FIXED EMBEDDING-BASED PERTURBATION RESULTS")
        print("="*80)
        
        print(f"Analysis Summary:")
        print(f"   • Cells processed: {len(self.cell_results)}")
        print(f"   • Valid perturbations: {len(self.cell_results[self.cell_results['target_in_cell'] == True])}")
        print(f"   • Genes analyzed: {len(self.gene_summary)}")
        print(f"   • Model: Geneformer V1-10M")
        print(f"   • Method: Cosine distance in embedding space")
        print(f"   • Gene mapping: Symbol -> Ensembl ID -> Token")
        
        if len(self.gene_summary) > 0:
            print(f"\nTop 10 Genes by Mean Cosine Distance:")
            print("-"*80)
            top_genes = self.gene_summary.head(10)
            
            for idx, (_, row) in enumerate(top_genes.iterrows()):
                fdr_status = "*** SIGNIFICANT" if row.get('empirical_fdr', 1.0) < 0.05 else ""
                print(f"{idx+1:2d}. {row['gene']:8s} | Distance: {row['mean_distance']:.6f} | "
                      f"Cells: {row['n_cells']:3d} | P: {row.get('empirical_p_value', np.nan):.4f} "
                      f"{fdr_status}")
                      
            # Check for significant results
            if 'empirical_fdr' in self.gene_summary.columns:
                significant = self.gene_summary[self.gene_summary['empirical_fdr'] < 0.05]
                print(f"\nSignificant results (FDR < 0.05): {len(significant)}")
                
                if len(significant) > 0:
                    print("\nSignificant genes:")
                    for _, row in significant.iterrows():
                        print(f"   • {row['gene']}: distance={row['mean_distance']:.6f}, "
                              f"FDR={row['empirical_fdr']:.4f}")
            else:
                print("\n[WARNING] No FDR correction applied (insufficient data)")
                
        print("\n" + "="*80)

def main():
    """Run Phase 4 analysis - FIXED VERSION"""
    print("STARTING Phase 4: FIXED Embedding-Based Perturbation Analysis")
    print(f"GPU: {torch.cuda.is_available()} ({'Available' if torch.cuda.is_available() else 'Not Available'})")
    
    # Initialize analyzer
    analyzer = FixedEmbeddingPerturbationAnalyzer(CONFIG)
    
    # Check if we have any target genes mapped
    if len(analyzer.target_ensembl) == 0:
        print("[ERROR] No target genes could be mapped! Cannot proceed.")
        return None
    
    # Run analysis
    try:
        # Analyze cells
        analyzer.analyze_cells()
        
        # Aggregate and test
        gene_summary = analyzer.aggregate_and_test()
        
        # Save results
        analyzer.save_results()
        
        # Print summary
        analyzer.print_summary()
        
        print("\n[SUCCESS] Phase 4 FIXED analysis completed successfully!")
        
        return gene_summary
        
    except Exception as e:
        print(f"\n[ERROR] Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    main()