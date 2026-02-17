#!/usr/bin/env python3
"""
Phase 4: PROPER Embedding-Based Perturbation Analysis
=====================================================

This pipeline implements REAL transformer-based perturbation using Geneformer embeddings:
1. Load Geneformer V1-10M model
2. Tokenize cells using proper Geneformer tokenization
3. Get baseline embeddings from transformer
4. Perturb by removing gene tokens and re-embedding  
5. Measure COSINE DISTANCE in embedding space
6. Statistical testing with proper null distributions

Critical: This uses actual model representations, not expression statistics!
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

try:
    from geneformer import TranscriptomeTokenizer, EmbExtractor
    from transformers import BertModel, BertConfig
    GENEFORMER_AVAILABLE = True
    print("[OK] Geneformer modules loaded successfully")
except ImportError as e:
    print(f"⚠️ Geneformer import failed: {e}")
    print("Will implement manual tokenization and embedding extraction")
    GENEFORMER_AVAILABLE = False

# Configuration
CONFIG = {
    'model_path': r"D:\openclaw\intelligence-augmentation\models\Geneformer\Geneformer-V1-10M",
    'data_path': r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad",
    'output_dir': r"D:\openclaw\intelligence-augmentation\analysis\results\phase4",
    'max_cells': 200,  # Pilot run
    'batch_size': 50,  # GPU memory constraint  
    'max_tokens': 2048,  # Geneformer context length
    'device': 'cuda' if torch.cuda.is_available() else 'cpu',
    'n_permutations': 100,  # For null distribution
    'target_genes': [
        'CADM2', 'NRXN1', 'NLGN1', 'NEGR1', 'TCF4', 'MEF2C', 'APP', 'GRIN2A', 'GRIN2B',
        'MAPT', 'SCN1A', 'DISC1', 'FOXO3', 'CAMK2A', 'HOMER1', 'CREB1', 'COMT', 'SHANK3',
        'FMR1', 'SYN1', 'BDNF', 'SLC6A4'
    ]
}

class EmbeddingPerturbationAnalyzer:
    """
    Proper embedding-based perturbation using Geneformer transformer model
    """
    
    def __init__(self, config):
        self.config = config
        self.device = torch.device(config['device'])
        self.results = []
        
        # Load token dictionary
        self.load_token_dictionary()
        
        # Load model
        self.load_model()
        
        print(f"🚀 Initialized on {self.device}")
        
    def load_token_dictionary(self):
        """Load Geneformer token dictionary"""
        token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc30M.pkl')
        if not os.path.exists(token_dict_path):
            # Try 104M version
            token_dict_path = os.path.join(GENEFORMER_PATH, 'geneformer', 'token_dictionary_gc104M.pkl')
            
        if os.path.exists(token_dict_path):
            with open(token_dict_path, 'rb') as f:
                self.token_dictionary = pickle.load(f)
            print(f"📚 Loaded token dictionary: {len(self.token_dictionary)} tokens")
        else:
            raise FileNotFoundError("Token dictionary not found!")
            
        # Create reverse mapping
        self.id_to_token = {v: k for k, v in self.token_dictionary.items()}
        
    def load_model(self):
        """Load Geneformer V1-10M model"""
        try:
            from transformers import BertModel, BertConfig
            
            config_path = os.path.join(self.config['model_path'], 'config.json')
            model_path = os.path.join(self.config['model_path'], 'pytorch_model.bin')
            
            # Load config
            bert_config = BertConfig.from_pretrained(config_path)
            
            # Load model
            self.model = BertModel.from_pretrained(self.config['model_path'])
            self.model.to(self.device)
            self.model.eval()
            
            print(f"🧠 Loaded Geneformer V1-10M ({bert_config.hidden_size}D embeddings)")
            
        except Exception as e:
            print(f"❌ Failed to load model: {e}")
            raise
            
    def tokenize_cell(self, gene_expr, max_tokens=2048):
        """
        Convert single cell expression to Geneformer tokens
        
        Args:
            gene_expr: Dictionary {gene_symbol: expression_value}  
            max_tokens: Maximum context length
            
        Returns:
            List of token IDs ranked by expression
        """
        # Filter genes that exist in token dictionary
        valid_genes = {gene: expr for gene, expr in gene_expr.items() 
                      if gene in self.token_dictionary and expr > 0}
        
        if not valid_genes:
            return []
            
        # Rank genes by expression (descending)
        ranked_genes = sorted(valid_genes.items(), key=lambda x: x[1], reverse=True)
        
        # Convert to token IDs
        token_ids = [self.token_dictionary[gene] for gene, _ in ranked_genes[:max_tokens-2]]
        
        # Add special tokens: [CLS] = 0, [SEP] = 1 (standard BERT format)
        token_ids = [0] + token_ids + [1]
        
        return token_ids
        
    def get_embedding(self, token_ids):
        """
        Get embedding for tokenized cell
        
        Args:
            token_ids: List of token IDs
            
        Returns:
            numpy array: [CLS] token embedding
        """
        if not token_ids or len(token_ids) < 2:
            return None
            
        # Convert to tensor
        input_ids = torch.tensor(token_ids, dtype=torch.long).unsqueeze(0).to(self.device)
        
        with torch.no_grad():
            # Get model outputs
            outputs = self.model(input_ids)
            hidden_states = outputs.last_hidden_state
            
            # Extract [CLS] token embedding (first token)
            cls_embedding = hidden_states[0, 0, :].cpu().numpy()
            
        return cls_embedding
        
    def perturb_and_measure(self, token_ids, target_gene_symbol, n_permutations=100):
        """
        Perturb cell by removing target gene and measure embedding distance
        
        Args:
            token_ids: Original token sequence
            target_gene_symbol: Gene to remove
            n_permutations: Number of random permutations for null distribution
            
        Returns:
            dict: Perturbation results
        """
        if target_gene_symbol not in self.token_dictionary:
            return None
            
        target_token_id = self.token_dictionary[target_gene_symbol]
        
        # Check if target gene is in this cell
        if target_token_id not in token_ids:
            return {
                'target_gene': target_gene_symbol,
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
            random_tokens = available_tokens * (n_permutations // len(available_tokens) + 1)
            random_tokens = random_tokens[:n_permutations]
            
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
            'target_in_cell': True,
            'cosine_distance': cosine_dist,
            'null_distances': null_distances,
            'p_value': p_value,
            'n_null_samples': len(null_distances)
        }
        
    def analyze_cells(self):
        """Main analysis pipeline"""
        print("🔬 Loading single-cell data...")
        
        # Load data
        adata = sc.read_h5ad(self.config['data_path'])
        print(f"📊 Loaded {adata.n_obs} cells × {adata.n_vars} genes")
        
        # Sample cells for pilot
        if adata.n_obs > self.config['max_cells']:
            sc.pp.subsample(adata, n_obs=self.config['max_cells'], random_state=42)
            print(f"📝 Subsampled to {adata.n_obs} cells")
            
        # Get gene names
        if 'feature_name' in adata.var.columns:
            gene_symbols = adata.var['feature_name'].values
        else:
            gene_symbols = adata.var.index.values
            
        print(f"🧬 Processing {len(self.config['target_genes'])} target genes")
        
        # Process in batches
        batch_results = []
        n_batches = (adata.n_obs + self.config['batch_size'] - 1) // self.config['batch_size']
        
        for batch_idx in range(n_batches):
            print(f"\n🔄 Processing batch {batch_idx + 1}/{n_batches}")
            
            start_idx = batch_idx * self.config['batch_size']
            end_idx = min(start_idx + self.config['batch_size'], adata.n_obs)
            
            batch_cells = adata[start_idx:end_idx]
            
            for cell_idx in tqdm(range(batch_cells.n_obs), desc="Cells"):
                cell_data = batch_cells[cell_idx]
                
                # Get cell expression as dictionary
                cell_expr = dict(zip(gene_symbols, cell_data.X.toarray().flatten()))
                
                # Tokenize cell
                token_ids = self.tokenize_cell(cell_expr, self.config['max_tokens'])
                
                if not token_ids:
                    continue
                    
                # Get cell type if available
                cell_type = cell_data.obs.get('cell_type', 'Unknown')
                
                # Test each target gene
                for target_gene in self.config['target_genes']:
                    result = self.perturb_and_measure(
                        token_ids, 
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
            torch.cuda.empty_cache()
            gc.collect()
            
            # Save intermediate results
            if batch_results:
                batch_df = pd.DataFrame(batch_results)
                batch_file = os.path.join(self.config['output_dir'], f'batch_{batch_idx + 1}_results.csv')
                batch_df.to_csv(batch_file, index=False)
                print(f"💾 Saved batch {batch_idx + 1} results: {len(batch_results)} perturbations")
                
        self.results = batch_results
        
    def aggregate_and_test(self):
        """Aggregate results and perform statistical testing"""
        if not self.results:
            print("❌ No results to aggregate!")
            return
            
        df = pd.DataFrame(self.results)
        
        # Filter out cells where target gene wasn't present
        df_valid = df[df['target_in_cell'] == True].copy()
        
        print(f"📈 Aggregating {len(df_valid)} valid perturbations")
        
        # Aggregate by gene
        gene_stats = []
        
        for gene in self.config['target_genes']:
            gene_data = df_valid[df_valid['target_gene'] == gene]
            
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
            if all_null_distances:
                # Wilcoxon rank-sum test
                statistic, wilcoxon_p = ranksums(distances, all_null_distances)
                
                # Empirical p-value (average across cells)
                empirical_p = gene_data['p_value'].dropna().mean()
            else:
                statistic, wilcoxon_p, empirical_p = np.nan, np.nan, np.nan
                
            gene_stats.append({
                'gene': gene,
                'n_cells': len(distances),
                'mean_distance': distances.mean(),
                'std_distance': distances.std(),
                'median_distance': distances.median(),
                'max_distance': distances.max(),
                'min_distance': distances.min(),
                'empirical_p_value': empirical_p,
                'wilcoxon_p_value': wilcoxon_p,
                'wilcoxon_statistic': statistic,
                'n_null_samples': len(all_null_distances)
            })
            
        # Convert to DataFrame and apply FDR correction
        gene_df = pd.DataFrame(gene_stats)
        
        if len(gene_df) > 0:
            # FDR correction for empirical p-values
            valid_empirical_p = gene_df['empirical_p_value'].dropna()
            if len(valid_empirical_p) > 0:
                _, empirical_fdr, _, _ = multipletests(valid_empirical_p, method='fdr_bh')
                gene_df.loc[valid_empirical_p.index, 'empirical_fdr'] = empirical_fdr
                
            # FDR correction for Wilcoxon p-values  
            valid_wilcoxon_p = gene_df['wilcoxon_p_value'].dropna()
            if len(valid_wilcoxon_p) > 0:
                _, wilcoxon_fdr, _, _ = multipletests(valid_wilcoxon_p, method='fdr_bh')
                gene_df.loc[valid_wilcoxon_p.index, 'wilcoxon_fdr'] = wilcoxon_fdr
                
            # Sort by mean distance (descending)
            gene_df = gene_df.sort_values('mean_distance', ascending=False)
            
        self.gene_summary = gene_df
        self.cell_results = df
        
        return gene_df
        
    def save_results(self):
        """Save all results to files"""
        # Save cell-level results
        cell_file = os.path.join(self.config['output_dir'], 'cell_perturbation_results_phase4.csv')
        self.cell_results.to_csv(cell_file, index=False)
        print(f"💾 Saved cell results: {cell_file}")
        
        # Save gene-level summary
        gene_file = os.path.join(self.config['output_dir'], 'gene_summary_phase4.csv')
        self.gene_summary.to_csv(gene_file, index=False)
        print(f"💾 Saved gene summary: {gene_file}")
        
        # Save analysis metadata
        metadata = {
            'config': self.config,
            'n_cells_processed': len(self.cell_results),
            'n_valid_perturbations': len(self.cell_results[self.cell_results['target_in_cell'] == True]),
            'n_genes_analyzed': len(self.gene_summary),
            'device': str(self.device),
            'model_path': self.config['model_path'],
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        metadata_file = os.path.join(self.config['output_dir'], 'analysis_metadata_phase4.json')
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
        print(f"💾 Saved metadata: {metadata_file}")
        
    def print_summary(self):
        """Print analysis summary"""
        print("\n" + "="*80)
        print("🎯 PHASE 4: EMBEDDING-BASED PERTURBATION RESULTS")
        print("="*80)
        
        print(f"📊 Analysis Summary:")
        print(f"   • Cells processed: {len(self.cell_results)}")
        print(f"   • Valid perturbations: {len(self.cell_results[self.cell_results['target_in_cell'] == True])}")
        print(f"   • Genes analyzed: {len(self.gene_summary)}")
        print(f"   • Model: Geneformer V1-10M")
        print(f"   • Method: Cosine distance in embedding space")
        
        if len(self.gene_summary) > 0:
            print(f"\n🏆 Top 10 Genes by Mean Cosine Distance:")
            print("-"*80)
            top_genes = self.gene_summary.head(10)
            
            for idx, row in top_genes.iterrows():
                fdr_status = "⭐ SIGNIFICANT" if row.get('empirical_fdr', 1.0) < 0.05 else ""
                print(f"{idx+1:2d}. {row['gene']:8s} | Distance: {row['mean_distance']:.6f} | "
                      f"Cells: {row['n_cells']:3d} | P: {row.get('empirical_p_value', np.nan):.4f} "
                      f"{fdr_status}")
                      
            # Check for significant results
            if 'empirical_fdr' in self.gene_summary.columns:
                significant = self.gene_summary[self.gene_summary['empirical_fdr'] < 0.05]
                print(f"\n✨ Significant results (FDR < 0.05): {len(significant)}")
                
                if len(significant) > 0:
                    print("\nSignificant genes:")
                    for _, row in significant.iterrows():
                        print(f"   • {row['gene']}: distance={row['mean_distance']:.6f}, "
                              f"FDR={row['empirical_fdr']:.4f}")
            else:
                print("\n⚠️ No FDR correction applied (insufficient data)")
                
        print("\n" + "="*80)

def main():
    """Run Phase 4 analysis"""
    print("🚀 Starting Phase 4: Embedding-Based Perturbation Analysis")
    print(f"🎯 GPU: {torch.cuda.is_available()} ({'✅ Available' if torch.cuda.is_available() else '❌ Not Available'})")
    
    # Initialize analyzer
    analyzer = EmbeddingPerturbationAnalyzer(CONFIG)
    
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
        
        print("\n✅ Phase 4 analysis completed successfully!")
        
        return gene_summary
        
    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    main()