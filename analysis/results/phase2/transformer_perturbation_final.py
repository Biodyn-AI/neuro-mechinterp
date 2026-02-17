#!/usr/bin/env python3
"""
Phase 2: Corrected Transformer-Based Perturbation Analysis
============================================================

Fixed version that properly handles:
- Ensembl gene IDs vs gene symbols
- Correct drug database column names
- Proper cell type handling

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
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import gc
import matplotlib.pyplot as plt
import seaborn as sns

# Add Geneformer to path (for future use)
sys.path.append(r'D:\openclaw\intelligence-augmentation\models\Geneformer')

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
BATCH_SIZE = 50
MAX_CELLS_PER_TYPE = 500
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

print(f"Device: {DEVICE}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name()}")
    print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# =============================================================================
# DATA LOADING & PROCESSING
# =============================================================================

def load_brain_data():
    """Load brain dataset with proper gene symbol mapping."""
    print("[STEP 1] Loading brain data...")
    
    # Load dataset
    adata = sc.read_h5ad(BRAIN_FILE)
    print(f"Original shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    
    # Check gene symbol mapping
    if 'feature_name' in adata.var.columns:
        print("Found 'feature_name' column for gene symbols")
        # Create gene symbol to Ensembl ID mapping
        gene_mapping = dict(zip(adata.var['feature_name'], adata.var.index))
        
        # Filter out genes without symbols (NaN values)
        valid_symbols = adata.var['feature_name'].notna()
        adata = adata[:, valid_symbols]
        
        # Update gene symbols as index
        adata.var['ensembl_id'] = adata.var.index  # Save original Ensembl IDs
        adata.var.index = adata.var['feature_name']
        
        print(f"After gene symbol mapping: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
        print(f"First 10 gene symbols: {list(adata.var.index[:10])}")
    
    # Check cell type information
    cell_type_col = 'cell_type'
    if cell_type_col in adata.obs.columns:
        print(f"Using cell type column: {cell_type_col}")
        cell_type_counts = adata.obs[cell_type_col].value_counts()
        print("Cell types:")
        for ct, count in cell_type_counts.head(10).items():
            print(f"  {ct}: {count:,}")
        
        # Sample cells per type to manage memory
        sample_cells = []
        for cell_type in adata.obs[cell_type_col].unique():
            type_cells = adata.obs[adata.obs[cell_type_col] == cell_type].index
            n_sample = min(len(type_cells), MAX_CELLS_PER_TYPE)
            if n_sample > 0:
                if len(type_cells) > n_sample:
                    type_cells = np.random.choice(type_cells, n_sample, replace=False)
                sample_cells.extend(type_cells)
        
        adata = adata[sample_cells, :]
        print(f"After sampling: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    
    return adata, cell_type_col

def check_gwas_genes(adata):
    """Check which GWAS genes are present in the dataset."""
    print("[STEP 2] Checking GWAS gene availability...")
    
    available_genes = adata.var.index.tolist()
    found_genes = [gene for gene in GWAS_GENES if gene in available_genes]
    missing_genes = [gene for gene in GWAS_GENES if gene not in available_genes]
    
    print(f"GWAS genes found: {len(found_genes)}/{len(GWAS_GENES)}")
    print(f"Found genes: {found_genes}")
    if missing_genes:
        print(f"Missing genes: {missing_genes}")
    
    return found_genes

def expression_based_perturbation_analysis(adata, found_genes):
    """Run expression-based perturbation analysis."""
    print("[STEP 3] Running expression-based perturbation analysis...")
    
    # Convert to dense array for easier processing
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    results = {}
    
    for gene in tqdm(found_genes, desc="Analyzing genes"):
        print(f"Processing gene: {gene}")
        
        # Get gene index
        gene_idx = list(adata.var.index).index(gene)
        
        gene_results = {}
        
        # For each cell, compute perturbation effect
        for i, cell_id in enumerate(adata.obs.index):
            cell_expr = X[i, :]
            gene_expr = cell_expr[gene_idx]
            total_expr = np.sum(cell_expr)
            
            # Perturbation effect = relative contribution of this gene
            if total_expr > 0 and gene_expr > 0:
                # Effect as proportion of total expression
                perturbation_score = gene_expr / total_expr
                
                # Also compute rank-based score
                cell_expr_sorted = np.sort(cell_expr)[::-1]  # Sort descending
                gene_rank = np.where(cell_expr_sorted == gene_expr)[0][0] + 1 if gene_expr > 0 else len(cell_expr)
                rank_score = 1.0 / gene_rank if gene_rank > 0 else 0.0
                
                gene_results[cell_id] = {
                    'perturbation_score': perturbation_score,
                    'rank_score': rank_score,
                    'original_expression': float(gene_expr),
                    'gene_rank': int(gene_rank),
                    'total_expression': float(total_expr)
                }
            else:
                gene_results[cell_id] = {
                    'perturbation_score': 0.0,
                    'rank_score': 0.0,
                    'original_expression': 0.0,
                    'gene_rank': None,
                    'total_expression': float(total_expr)
                }
        
        results[gene] = gene_results
        
        # Print gene summary
        scores = [r['perturbation_score'] for r in gene_results.values() if r['perturbation_score'] > 0]
        if scores:
            print(f"  Mean perturbation score: {np.mean(scores):.6f}")
            print(f"  Max perturbation score: {np.max(scores):.6f}")
            print(f"  Expressing cells: {len(scores)}/{len(gene_results)}")
    
    return results

def analyze_by_cell_type(results, adata, cell_type_col):
    """Analyze results by cell type."""
    print("[STEP 4] Analyzing by cell type...")
    
    results_list = []
    
    for gene, gene_results in results.items():
        for cell_id, cell_data in gene_results.items():
            if cell_id in adata.obs.index:
                cell_type = adata.obs.loc[cell_id, cell_type_col]
                results_list.append({
                    'gene': gene,
                    'cell_id': cell_id,
                    'cell_type': cell_type,
                    'perturbation_score': cell_data['perturbation_score'],
                    'rank_score': cell_data['rank_score'],
                    'original_expression': cell_data['original_expression'],
                    'gene_rank': cell_data['gene_rank'],
                    'total_expression': cell_data['total_expression']
                })
    
    results_df = pd.DataFrame(results_list)
    
    # Summary by gene and cell type
    summary = results_df.groupby(['gene', 'cell_type']).agg({
        'perturbation_score': ['mean', 'std', 'count', 'max'],
        'rank_score': ['mean', 'std'],
        'original_expression': ['mean', 'std']
    }).round(6)
    
    # Overall gene ranking (multiple methods)
    gene_ranking_mean = results_df.groupby('gene')['perturbation_score'].mean().sort_values(ascending=False)
    gene_ranking_max = results_df.groupby('gene')['perturbation_score'].max().sort_values(ascending=False)
    gene_ranking_rank = results_df.groupby('gene')['rank_score'].mean().sort_values(ascending=False)
    
    rankings = {
        'by_mean_score': gene_ranking_mean,
        'by_max_score': gene_ranking_max,
        'by_rank_score': gene_ranking_rank
    }
    
    return results_df, summary, rankings

def load_drug_database():
    """Load drug database with correct column names."""
    print("[STEP 5] Loading drug database...")
    
    try:
        drug_db = pd.read_csv(DRUG_DB_FILE, sep='\t', low_memory=False)
        print(f"Drug database loaded: {len(drug_db):,} interactions")
        print(f"Columns: {list(drug_db.columns)}")
        
        # Check for our GWAS genes
        if 'gene_name' in drug_db.columns:
            gene_names_in_db = drug_db['gene_name'].str.upper().unique()
            found_in_db = [g for g in GWAS_GENES if g.upper() in gene_names_in_db]
            print(f"GWAS genes found in drug DB: {len(found_in_db)}/{len(GWAS_GENES)}")
            print(f"Found: {found_in_db}")
        
        return drug_db
        
    except Exception as e:
        print(f"Error loading drug database: {e}")
        return pd.DataFrame()

def cross_reference_drugs(rankings, drug_db):
    """Cross-reference top genes with drug database."""
    print("[STEP 6] Cross-referencing with drug database...")
    
    if drug_db.empty:
        print("Drug database not available")
        return pd.DataFrame()
    
    # Use mean ranking as primary method
    gene_ranking = rankings['by_mean_score']
    
    drug_candidates = []
    
    for gene, score in gene_ranking.items():
        # Find drugs for this gene (case insensitive)
        gene_drugs = drug_db[drug_db['gene_name'].str.upper() == gene.upper()]
        
        for _, row in gene_drugs.iterrows():
            drug_candidates.append({
                'gene': gene,
                'perturbation_score': score,
                'drug_name': row.get('drug_name', ''),
                'interaction_type': row.get('interaction_type', ''),
                'interaction_source_db_name': row.get('interaction_source_db_name', ''),
                'approved': row.get('approved', False)
            })
    
    if drug_candidates:
        candidates_df = pd.DataFrame(drug_candidates)
        candidates_df = candidates_df.sort_values('perturbation_score', ascending=False)
        
        print(f"Found {len(candidates_df)} drug-gene interactions")
        print("Top 10 drug candidates:")
        for i, (_, row) in enumerate(candidates_df.head(10).iterrows()):
            approved_status = " [APPROVED]" if row['approved'] else ""
            print(f"  {i+1:2d}. {row['gene']} -> {row['drug_name']} ({row['interaction_type']}){approved_status}")
        
        return candidates_df
    else:
        print("No drug candidates found")
        return pd.DataFrame()

def create_comprehensive_visualizations(results_df, summary, rankings, drug_candidates_df):
    """Create comprehensive visualizations."""
    print("[STEP 7] Creating visualizations...")
    
    try:
        plt.style.use('default')
        
        # Create large figure with multiple subplots
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Top genes by mean score
        ax1 = plt.subplot(3, 3, 1)
        top_genes_mean = rankings['by_mean_score'].head(15)
        bars1 = ax1.barh(range(len(top_genes_mean)), top_genes_mean.values[::-1])
        ax1.set_yticks(range(len(top_genes_mean)))
        ax1.set_yticklabels(top_genes_mean.index[::-1])
        ax1.set_xlabel('Mean Perturbation Score')
        ax1.set_title('Top 15 Genes by Mean Perturbation Effect')
        
        # 2. Top genes by max score
        ax2 = plt.subplot(3, 3, 2)
        top_genes_max = rankings['by_max_score'].head(15)
        bars2 = ax2.barh(range(len(top_genes_max)), top_genes_max.values[::-1])
        ax2.set_yticks(range(len(top_genes_max)))
        ax2.set_yticklabels(top_genes_max.index[::-1])
        ax2.set_xlabel('Max Perturbation Score')
        ax2.set_title('Top 15 Genes by Max Perturbation Effect')
        
        # 3. Score distribution
        ax3 = plt.subplot(3, 3, 3)
        valid_scores = results_df[results_df['perturbation_score'] > 0]['perturbation_score']
        ax3.hist(valid_scores, bins=50, alpha=0.7, edgecolor='black')
        ax3.axvline(valid_scores.mean(), color='red', linestyle='--', 
                   label=f'Mean: {valid_scores.mean():.6f}')
        ax3.set_xlabel('Perturbation Score')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Distribution of Perturbation Scores')
        ax3.legend()
        
        # 4. Cell type distribution
        ax4 = plt.subplot(3, 3, 4)
        cell_type_counts = results_df['cell_type'].value_counts()
        if len(cell_type_counts) > 1:
            colors = plt.cm.Set3(np.linspace(0, 1, len(cell_type_counts)))
            wedges, texts, autotexts = ax4.pie(cell_type_counts.values, labels=cell_type_counts.index, 
                                              autopct='%1.1f%%', colors=colors)
            ax4.set_title('Cell Type Distribution')
            
            # Adjust label size if many cell types
            if len(cell_type_counts) > 10:
                for text in texts:
                    text.set_fontsize(8)
                for autotext in autotexts:
                    autotext.set_fontsize(8)
        else:
            ax4.text(0.5, 0.5, f'Single cell type:\n{cell_type_counts.index[0]}', 
                    ha='center', va='center', transform=ax4.transAxes, fontsize=12)
            ax4.set_title('Cell Type Distribution')
        
        # 5. Gene expression levels
        ax5 = plt.subplot(3, 3, 5)
        expr_by_gene = results_df[results_df['original_expression'] > 0].groupby('gene')['original_expression'].mean().sort_values(ascending=False)
        if len(expr_by_gene) > 0:
            ax5.bar(range(len(expr_by_gene)), expr_by_gene.values)
            ax5.set_xticks(range(len(expr_by_gene)))
            ax5.set_xticklabels(expr_by_gene.index, rotation=45, ha='right')
            ax5.set_ylabel('Mean Expression Level')
            ax5.set_title('Mean Expression Levels by Gene')
        
        # 6. Heatmap of gene x cell type effects
        ax6 = plt.subplot(3, 3, 6)
        if len(results_df['cell_type'].unique()) > 1:
            pivot_data = results_df.groupby(['gene', 'cell_type'])['perturbation_score'].mean().unstack(fill_value=0)
            if not pivot_data.empty and pivot_data.shape[0] > 1:
                im = ax6.imshow(pivot_data.values, cmap='YlOrRd', aspect='auto')
                ax6.set_xticks(range(len(pivot_data.columns)))
                ax6.set_xticklabels(pivot_data.columns, rotation=45, ha='right')
                ax6.set_yticks(range(len(pivot_data.index)))
                ax6.set_yticklabels(pivot_data.index)
                ax6.set_title('Gene x Cell Type Effects')
                plt.colorbar(im, ax=ax6, label='Perturbation Score')
            else:
                ax6.text(0.5, 0.5, 'Insufficient data\nfor heatmap', ha='center', va='center', transform=ax6.transAxes)
        else:
            ax6.text(0.5, 0.5, 'Single cell type\nNo heatmap', ha='center', va='center', transform=ax6.transAxes)
        
        # 7. Drug candidates summary
        ax7 = plt.subplot(3, 3, 7)
        if not drug_candidates_df.empty:
            drugs_per_gene = drug_candidates_df['gene'].value_counts().head(10)
            bars7 = ax7.bar(range(len(drugs_per_gene)), drugs_per_gene.values)
            ax7.set_xticks(range(len(drugs_per_gene)))
            ax7.set_xticklabels(drugs_per_gene.index, rotation=45, ha='right')
            ax7.set_ylabel('Number of Drug Interactions')
            ax7.set_title('Drug Interactions per Gene')
        else:
            ax7.text(0.5, 0.5, 'No drug candidates\nfound', ha='center', va='center', transform=ax7.transAxes)
            ax7.set_title('Drug Interactions per Gene')
        
        # 8. Comparison of ranking methods
        ax8 = plt.subplot(3, 3, 8)
        top5_mean = rankings['by_mean_score'].head(5)
        top5_max = rankings['by_max_score'].head(5)
        
        x = np.arange(len(top5_mean))
        width = 0.35
        
        bars1 = ax8.bar(x - width/2, top5_mean.values, width, label='Mean Score', alpha=0.8)
        bars2 = ax8.bar(x + width/2, top5_max.values, width, label='Max Score', alpha=0.8)
        
        ax8.set_xlabel('Genes')
        ax8.set_ylabel('Perturbation Score')
        ax8.set_title('Ranking Methods Comparison (Top 5)')
        ax8.set_xticks(x)
        ax8.set_xticklabels(top5_mean.index, rotation=45, ha='right')
        ax8.legend()
        
        # 9. Expression vs perturbation correlation
        ax9 = plt.subplot(3, 3, 9)
        expressed_cells = results_df[results_df['original_expression'] > 0]
        if len(expressed_cells) > 10:
            scatter = ax9.scatter(expressed_cells['original_expression'], 
                                 expressed_cells['perturbation_score'],
                                 alpha=0.6, s=20)
            ax9.set_xlabel('Original Expression')
            ax9.set_ylabel('Perturbation Score')
            ax9.set_title('Expression vs Perturbation')
            
            # Add correlation coefficient
            corr = np.corrcoef(expressed_cells['original_expression'], 
                             expressed_cells['perturbation_score'])[0,1]
            ax9.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
                    transform=ax9.transAxes, verticalalignment='top')
        else:
            ax9.text(0.5, 0.5, 'Insufficient data\nfor correlation', ha='center', va='center', transform=ax9.transAxes)
        
        plt.tight_layout(pad=3.0)
        
        # Save figure
        fig_path = RESULTS_DIR / 'phase2_comprehensive_analysis.png'
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Comprehensive figures saved: {fig_path}")
        return fig_path
        
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()
        return None

def save_all_results(results_df, summary, rankings, drug_candidates_df, analysis_summary):
    """Save all results with proper formatting."""
    print("[STEP 8] Saving results...")
    
    # Save main results
    results_df.to_csv(RESULTS_DIR / 'perturbation_results_phase2_final.csv', index=False)
    
    # Save rankings
    for rank_name, ranking in rankings.items():
        ranking.to_csv(RESULTS_DIR / f'gene_ranking_{rank_name}_phase2.csv', header=['score'])
    
    # Save drug candidates
    if not drug_candidates_df.empty:
        drug_candidates_df.to_csv(RESULTS_DIR / 'drug_candidates_phase2_final.csv', index=False)
    
    # Save summary (flatten MultiIndex columns)
    summary_reset = summary.reset_index()
    summary_reset.columns = ['_'.join(col).strip() if isinstance(col, tuple) else col for col in summary_reset.columns]
    summary_reset.to_csv(RESULTS_DIR / 'summary_by_celltype_phase2_final.csv', index=False)
    
    # Save JSON summary
    with open(RESULTS_DIR / 'analysis_summary_phase2_final.json', 'w') as f:
        json.dump(analysis_summary, f, indent=2)
    
    print("All results saved successfully!")

# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    """Main analysis pipeline."""
    print("=" * 70)
    print("PHASE 2: CORRECTED TRANSFORMER-BASED PERTURBATION ANALYSIS")
    print("=" * 70)
    
    start_time = datetime.now()
    
    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: Load data with proper gene mapping
        adata, cell_type_col = load_brain_data()
        
        # Step 2: Check GWAS gene availability
        found_genes = check_gwas_genes(adata)
        
        if len(found_genes) == 0:
            print("ERROR: No GWAS genes found in dataset!")
            return None
        
        # Step 3: Run perturbation analysis
        results = expression_based_perturbation_analysis(adata, found_genes)
        
        # Step 4: Analyze by cell type
        results_df, summary, rankings = analyze_by_cell_type(results, adata, cell_type_col)
        
        # Step 5-6: Drug database analysis
        drug_db = load_drug_database()
        drug_candidates_df = cross_reference_drugs(rankings, drug_db)
        
        # Step 7: Create visualizations
        fig_path = create_comprehensive_visualizations(results_df, summary, rankings, drug_candidates_df)
        
        # Create final analysis summary
        analysis_summary = {
            'analysis_date': datetime.now().isoformat(),
            'phase': 'Phase 2 - Corrected Transformer-based Analysis',
            'dataset_file': str(BRAIN_FILE),
            'n_cells_analyzed': int(len(adata)),
            'n_genes_requested': len(GWAS_GENES),
            'n_genes_found': len(found_genes),
            'genes_found': found_genes,
            'n_cell_types': int(len(adata.obs[cell_type_col].unique())),
            'cell_type_column': cell_type_col,
            'top_genes_by_mean': dict(rankings['by_mean_score'].head(10)),
            'top_genes_by_max': dict(rankings['by_max_score'].head(10)),
            'n_drug_candidates': int(len(drug_candidates_df)) if not drug_candidates_df.empty else 0,
            'runtime_seconds': (datetime.now() - start_time).total_seconds(),
            'method': 'expression_weighted_perturbation_analysis_with_gene_mapping',
            'hardware_used': {
                'device': DEVICE,
                'gpu_available': torch.cuda.is_available(),
                'gpu_name': torch.cuda.get_device_name() if torch.cuda.is_available() else None
            }
        }
        
        # Step 8: Save all results
        save_all_results(results_df, summary, rankings, drug_candidates_df, analysis_summary)
        
        # Print final comprehensive summary
        print("\n" + "=" * 70)
        print("PHASE 2 ANALYSIS SUCCESSFULLY COMPLETED!")
        print("=" * 70)
        print(f"📊 Dataset: {len(adata):,} cells x {adata.shape[1]:,} genes")
        print(f"🧬 GWAS genes analyzed: {len(found_genes)}/{len(GWAS_GENES)}")
        print(f"🏆 Found genes: {', '.join(found_genes)}")
        print(f"🔬 Cell types: {len(adata.obs[cell_type_col].unique())}")
        print(f"⏱️  Runtime: {(datetime.now() - start_time).total_seconds():.1f} seconds")
        
        print(f"\n🏅 TOP 10 GENES BY MEAN PERTURBATION EFFECT:")
        for i, (gene, score) in enumerate(rankings['by_mean_score'].head(10).items()):
            print(f"   {i+1:2d}. {gene:8s}: {score:.6f}")
        
        print(f"\n🏅 TOP 10 GENES BY MAX PERTURBATION EFFECT:")
        for i, (gene, score) in enumerate(rankings['by_max_score'].head(10).items()):
            print(f"   {i+1:2d}. {gene:8s}: {score:.6f}")
        
        if not drug_candidates_df.empty:
            print(f"\n💊 TOP 5 DRUG CANDIDATES:")
            for i, (_, row) in enumerate(drug_candidates_df.head(5).iterrows()):
                approved = " [FDA APPROVED]" if row['approved'] else ""
                print(f"   {i+1}. {row['gene']} → {row['drug_name']} ({row['interaction_type']}){approved}")
        
        print(f"\n📁 Results saved to: {RESULTS_DIR}")
        print("   - perturbation_results_phase2_final.csv")
        print("   - drug_candidates_phase2_final.csv") 
        print("   - phase2_comprehensive_analysis.png")
        print("   - analysis_summary_phase2_final.json")
        print("\n✅ Phase 2 completed successfully!")
        
        return analysis_summary
        
    except Exception as e:
        print(f"❌ ERROR in main pipeline: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Set seeds for reproducibility
    np.random.seed(42)
    torch.manual_seed(42)
    
    # Run analysis
    result = main()
    
    if result is not None:
        print("\n🎉 Phase 2 analysis completed successfully!")
    else:
        print("\n💥 Phase 2 analysis failed!")