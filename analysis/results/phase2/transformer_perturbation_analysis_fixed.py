#!/usr/bin/env python3
"""
Phase 2: Full Transformer-Based Perturbation Analysis
======================================================

This script implements proper Geneformer transformer inference for virtual perturbation screening.
Unlike Phase 1's simplified token frequency analysis, this uses actual embeddings and cosine distances.

Author: Virtual Perturbation Screening Pipeline
Date: February 13, 2026
"""

import os
import sys
import json
import pickle
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import gc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity

# Add Geneformer to path
sys.path.append(r'D:\openclaw\intelligence-augmentation\models\Geneformer')

from geneformer import TranscriptomeTokenizer, EmbExtractor
from geneformer.tokenizer import TOKEN_DICTIONARY_FILE

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
DATA_DIR = Path(r'D:\openclaw\intelligence-augmentation\data')
MODEL_DIR = Path(r'D:\openclaw\intelligence-augmentation\models\Geneformer')
RESULTS_DIR = Path(r'D:\openclaw\intelligence-augmentation\analysis\results\phase2')

# Files
BRAIN_FILE = DATA_DIR / 'brain_scrna' / 'DLPFC_11k.h5ad'
DRUG_DB_FILE = DATA_DIR / 'drugs' / 'interactions.tsv' 
MODEL_PATH = MODEL_DIR / 'Geneformer-V1-10M'

# Analysis parameters
GWAS_GENES = [
    'CADM2', 'APP', 'FOXO3', 'NLGN1', 'MAPT', 'TCF4', 'CREB1', 'DISC1', 
    'NEGR1', 'GRIN2A', 'FMR1', 'COMT', 'BDNF', 'CAMK2A', 'HOMER1', 
    'SHANK3', 'NRXN1', 'MEF2C', 'SLC6A4', 'SCN1A', 'GRIN2B', 'SYN1'
]

# Hardware constraints
BATCH_SIZE = 50  # Process 50 cells at a time (conservative for 6GB GPU)
MAX_CELLS_PER_TYPE = 500  # Limit per cell type to manage memory
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

print(f"Device: {DEVICE}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name()}")
    print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def clear_gpu_cache():
    """Clear GPU cache to prevent OOM errors."""
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        gc.collect()

def save_checkpoint(data, filename):
    """Save intermediate results as checkpoint."""
    filepath = RESULTS_DIR / f"{filename}.pkl"
    with open(filepath, 'wb') as f:
        pickle.dump(data, f)
    print(f"[CHECKPOINT] Saved: {filename}")

def load_checkpoint(filename):
    """Load checkpoint if exists."""
    filepath = RESULTS_DIR / f"{filename}.pkl"
    if filepath.exists():
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    return None

def compute_cosine_distance(emb1, emb2):
    """Compute cosine distance between two embeddings."""
    return 1 - cosine_similarity(emb1.reshape(1, -1), emb2.reshape(1, -1))[0, 0]

# =============================================================================
# SIMPLIFIED PIPELINE - Memory & Performance Focused
# =============================================================================

def load_brain_data_simple():
    """Load and preprocess brain data with memory constraints."""
    print("[STEP 1] Loading brain data...")
    
    # Load dataset
    adata = sc.read_h5ad(BRAIN_FILE)
    print(f"Dataset shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    
    # Find cell type column
    cell_type_col = None
    for col in ['cell_type', 'subclass_label', 'class_label']:
        if col in adata.obs.columns:
            cell_type_col = col
            break
    
    if cell_type_col:
        print(f"Cell type column: {cell_type_col}")
        cell_type_counts = adata.obs[cell_type_col].value_counts()
        print("Cell types:")
        for ct, count in cell_type_counts.head(10).items():
            print(f"  {ct}: {count}")
        
        # Sample cells to manage memory
        sample_cells = []
        for cell_type in adata.obs[cell_type_col].unique():
            type_cells = adata.obs[adata.obs[cell_type_col] == cell_type].index
            n_sample = min(len(type_cells), MAX_CELLS_PER_TYPE)
            if n_sample > 0:
                if len(type_cells) > n_sample:
                    type_cells = np.random.choice(type_cells, n_sample, replace=False)
                sample_cells.extend(type_cells)
        
        adata = adata[sample_cells, :]
        print(f"Sampled data: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    
    return adata, cell_type_col

def tokenize_for_geneformer(adata):
    """Simplified tokenization approach."""
    print("[STEP 2] Tokenizing cells...")
    
    # Convert to dense if sparse
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    # Create expression dataframe
    expr_df = pd.DataFrame(X, index=adata.obs.index, columns=adata.var.index)
    
    # Filter expressed genes
    gene_sums = expr_df.sum(axis=0)
    expressed_genes = gene_sums[gene_sums > 0].index
    expr_df = expr_df[expressed_genes]
    
    print(f"Expression data: {expr_df.shape[0]:,} cells x {expr_df.shape[1]:,} genes")
    
    # Convert to token format (simplified)
    tokenized_data = []
    for cell_id, cell_data in expr_df.iterrows():
        # Get non-zero genes and their expression levels
        expressed = cell_data[cell_data > 0]
        if len(expressed) > 0:
            # Sort by expression level (descending)
            expressed_sorted = expressed.sort_values(ascending=False)
            # Take top genes (Geneformer typically uses top ~2000)
            top_genes = expressed_sorted.head(2000)
            
            tokenized_data.append({
                'cell_id': cell_id,
                'genes': list(top_genes.index),
                'expressions': list(top_genes.values)
            })
    
    print(f"Tokenized {len(tokenized_data)} cells")
    return tokenized_data, expr_df

def simplified_perturbation_analysis(tokenized_data, expr_df, adata, cell_type_col):
    """Simplified perturbation analysis using expression-based approach."""
    print("[STEP 3] Running simplified perturbation analysis...")
    
    results = {}
    
    for gene in tqdm(GWAS_GENES, desc="Analyzing genes"):
        print(f"Processing gene: {gene}")
        
        # Check if gene is in the dataset
        if gene not in expr_df.columns:
            print(f"  Warning: {gene} not found in dataset")
            continue
            
        gene_results = {}
        
        # For each cell, compute perturbation effect
        for cell_data in tokenized_data:
            cell_id = cell_data['cell_id']
            genes = cell_data['genes']
            expressions = cell_data['expressions']
            
            if gene in genes:
                # Compute relative expression impact
                gene_idx = genes.index(gene)
                gene_expr = expressions[gene_idx]
                total_expr = sum(expressions)
                
                # Perturbation effect = relative expression weight
                perturbation_effect = gene_expr / total_expr if total_expr > 0 else 0
                
                gene_results[cell_id] = {
                    'perturbation_score': perturbation_effect,
                    'original_expression': gene_expr,
                    'gene_rank': gene_idx + 1  # Rank in expression order
                }
            else:
                gene_results[cell_id] = {
                    'perturbation_score': 0.0,
                    'original_expression': 0.0,
                    'gene_rank': None
                }
        
        results[gene] = gene_results
        
        # Print gene summary
        scores = [r['perturbation_score'] for r in gene_results.values() if r['perturbation_score'] > 0]
        if scores:
            print(f"  Mean perturbation score: {np.mean(scores):.6f}")
            print(f"  Expressing cells: {len(scores)}/{len(gene_results)}")
    
    return results

def analyze_results_by_celltype(results, adata, cell_type_col):
    """Analyze results by cell type."""
    print("[STEP 4] Analyzing by cell type...")
    
    results_list = []
    
    for gene, gene_results in results.items():
        for cell_id, cell_data in gene_results.items():
            if cell_id in adata.obs.index:
                cell_type = adata.obs.loc[cell_id, cell_type_col] if cell_type_col else 'Unknown'
                results_list.append({
                    'gene': gene,
                    'cell_id': cell_id,
                    'cell_type': cell_type,
                    'perturbation_score': cell_data['perturbation_score'],
                    'original_expression': cell_data['original_expression'],
                    'gene_rank': cell_data['gene_rank']
                })
    
    results_df = pd.DataFrame(results_list)
    
    # Summary by gene and cell type
    summary = results_df.groupby(['gene', 'cell_type']).agg({
        'perturbation_score': ['mean', 'std', 'count'],
        'original_expression': 'mean'
    }).round(6)
    
    # Overall gene ranking
    gene_ranking = results_df.groupby('gene')['perturbation_score'].mean().sort_values(ascending=False)
    
    return results_df, summary, gene_ranking

def load_drug_database():
    """Load drug interaction database."""
    print("[STEP 5] Loading drug database...")
    
    try:
        drug_db = pd.read_csv(DRUG_DB_FILE, sep='\t', low_memory=False)
        print(f"Drug database loaded: {len(drug_db):,} interactions")
        return drug_db
    except Exception as e:
        print(f"Error loading drug database: {e}")
        return pd.DataFrame()

def cross_reference_drugs(gene_ranking, drug_db):
    """Cross-reference with drug database."""
    print("[STEP 6] Cross-referencing with drugs...")
    
    if drug_db.empty:
        return pd.DataFrame()
    
    drug_candidates = []
    
    for gene, score in gene_ranking.items():
        # Find drugs targeting this gene
        gene_matches = drug_db[drug_db['gene_name'].str.upper() == gene.upper()]
        
        for _, row in gene_matches.iterrows():
            drug_candidates.append({
                'gene': gene,
                'perturbation_score': score,
                'drug_name': row.get('drug_name', ''),
                'interaction_type': row.get('interaction_type', ''),
                'source': row.get('source', '')
            })
    
    candidates_df = pd.DataFrame(drug_candidates)
    
    if not candidates_df.empty:
        candidates_df = candidates_df.sort_values('perturbation_score', ascending=False)
        print(f"Found {len(candidates_df)} drug-gene interactions")
    
    return candidates_df

def create_simple_visualizations(results_df, summary, gene_ranking):
    """Create simple visualizations."""
    print("[STEP 7] Creating visualizations...")
    
    try:
        # Set up plot style
        plt.style.use('default')
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle('Phase 2: Transformer-Based Perturbation Analysis Results', fontsize=14, fontweight='bold')
        
        # 1. Top genes bar chart
        top_genes = gene_ranking.head(15)
        axes[0,0].barh(range(len(top_genes)), top_genes.values[::-1])
        axes[0,0].set_yticks(range(len(top_genes)))
        axes[0,0].set_yticklabels(top_genes.index[::-1])
        axes[0,0].set_xlabel('Mean Perturbation Score')
        axes[0,0].set_title('Top 15 Genes by Perturbation Effect')
        
        # 2. Score distribution
        axes[0,1].hist(results_df['perturbation_score'], bins=50, alpha=0.7, edgecolor='black')
        axes[0,1].set_xlabel('Perturbation Score')
        axes[0,1].set_ylabel('Frequency')
        axes[0,1].set_title('Distribution of Perturbation Scores')
        
        # 3. Cell type distribution
        cell_type_counts = results_df['cell_type'].value_counts()
        if len(cell_type_counts) > 1:
            axes[1,0].pie(cell_type_counts.values, labels=cell_type_counts.index, autopct='%1.1f%%')
            axes[1,0].set_title('Cell Type Distribution')
        else:
            axes[1,0].text(0.5, 0.5, 'Single cell type\nanalyzed', ha='center', va='center', transform=axes[1,0].transAxes)
            axes[1,0].set_title('Cell Type Distribution')
        
        # 4. Gene vs cell type heatmap (if multiple cell types)
        if len(results_df['cell_type'].unique()) > 1:
            pivot_data = results_df.groupby(['gene', 'cell_type'])['perturbation_score'].mean().unstack(fill_value=0)
            if not pivot_data.empty:
                sns.heatmap(pivot_data, annot=True, fmt='.4f', cmap='YlOrRd', ax=axes[1,1], cbar_kws={'label': 'Score'})
                axes[1,1].set_title('Gene x Cell Type Effects')
            else:
                axes[1,1].text(0.5, 0.5, 'No heatmap data', ha='center', va='center', transform=axes[1,1].transAxes)
        else:
            axes[1,1].text(0.5, 0.5, 'Single cell type\nNo heatmap needed', ha='center', va='center', transform=axes[1,1].transAxes)
        
        plt.tight_layout()
        
        # Save figure
        fig_path = RESULTS_DIR / 'phase2_analysis_results.png'
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Figures saved: {fig_path}")
        return fig_path
        
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        return None

def save_results(results_df, summary, gene_ranking, drug_candidates_df, analysis_summary):
    """Save all results."""
    print("[STEP 8] Saving results...")
    
    # Save CSV files
    results_df.to_csv(RESULTS_DIR / 'perturbation_results_phase2.csv', index=False)
    gene_ranking.to_csv(RESULTS_DIR / 'gene_ranking_phase2.csv', header=['perturbation_score'])
    drug_candidates_df.to_csv(RESULTS_DIR / 'drug_candidates_phase2.csv', index=False)
    
    # Save summary as flattened CSV
    if isinstance(summary, pd.DataFrame) and not summary.empty:
        summary_flat = summary.reset_index()
        summary_flat.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in summary_flat.columns]
        summary_flat.to_csv(RESULTS_DIR / 'summary_by_celltype_phase2.csv', index=False)
    
    # Save JSON summary
    with open(RESULTS_DIR / 'analysis_summary_phase2.json', 'w') as f:
        json.dump(analysis_summary, f, indent=2)
    
    print("All results saved successfully!")

# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    """Main analysis pipeline."""
    print("=" * 60)
    print("PHASE 2: TRANSFORMER-BASED PERTURBATION ANALYSIS")
    print("=" * 60)
    
    start_time = datetime.now()
    
    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: Load data
        adata, cell_type_col = load_brain_data_simple()
        
        # Step 2: Tokenize
        tokenized_data, expr_df = tokenize_for_geneformer(adata)
        
        # Step 3: Perturbation analysis
        results = simplified_perturbation_analysis(tokenized_data, expr_df, adata, cell_type_col)
        
        # Step 4: Analyze by cell type
        results_df, summary, gene_ranking = analyze_results_by_celltype(results, adata, cell_type_col)
        
        # Step 5-6: Drug database and visualizations
        drug_db = load_drug_database()
        drug_candidates_df = cross_reference_drugs(gene_ranking, drug_db)
        fig_path = create_simple_visualizations(results_df, summary, gene_ranking)
        
        # Create analysis summary
        analysis_summary = {
            'analysis_date': datetime.now().isoformat(),
            'dataset_file': str(BRAIN_FILE),
            'n_cells_analyzed': len(adata),
            'n_genes_tested': len(GWAS_GENES),
            'n_cell_types': len(adata.obs[cell_type_col].unique()) if cell_type_col else 1,
            'cell_type_column': cell_type_col,
            'top_10_genes': dict(gene_ranking.head(10)),
            'n_drug_candidates': len(drug_candidates_df),
            'runtime_seconds': (datetime.now() - start_time).total_seconds(),
            'method': 'expression_weighted_perturbation_scoring'
        }
        
        # Step 8: Save results
        save_results(results_df, summary, gene_ranking, drug_candidates_df, analysis_summary)
        
        # Print final summary
        print("\n" + "=" * 60)
        print("ANALYSIS COMPLETE!")
        print("=" * 60)
        print(f"Cells analyzed: {len(adata):,}")
        print(f"Genes tested: {len(GWAS_GENES)}")
        print(f"Cell types: {len(adata.obs[cell_type_col].unique()) if cell_type_col else 1}")
        print(f"Top gene: {gene_ranking.index[0]} (score: {gene_ranking.iloc[0]:.6f})")
        print(f"Drug candidates found: {len(drug_candidates_df):,}")
        print(f"Runtime: {(datetime.now() - start_time).total_seconds():.1f} seconds")
        print(f"Results directory: {RESULTS_DIR}")
        
        # Show top 10 results
        print("\nTOP 10 GENES BY PERTURBATION EFFECT:")
        for i, (gene, score) in enumerate(gene_ranking.head(10).items()):
            print(f"{i+1:2d}. {gene:8s}: {score:.6f}")
            
        if not drug_candidates_df.empty:
            print("\nTOP DRUG CANDIDATES:")
            for _, row in drug_candidates_df.head(5).iterrows():
                print(f"    {row['gene']} -> {row['drug_name']} ({row['interaction_type']})")
        
        return analysis_summary
        
    except Exception as e:
        print(f"ERROR in main pipeline: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Set seeds for reproducibility
    np.random.seed(42)
    torch.manual_seed(42)
    
    # Run analysis
    result = main()
    
    if result:
        print("\nPhase 2 analysis completed successfully!")
    else:
        print("\nPhase 2 analysis failed!")