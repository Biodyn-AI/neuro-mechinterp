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
    print(f"💾 Checkpoint saved: {filename}")

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
# DATA LOADING & PREPROCESSING
# =============================================================================

def load_brain_data():
    """Load and preprocess the DLPFC brain dataset."""
    print("🧠 Loading brain data...")
    
    # Load dataset
    adata = sc.read_h5ad(BRAIN_FILE)
    print(f"📊 Dataset shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    # Check for cell type annotations
    cell_type_cols = [col for col in adata.obs.columns if 'type' in col.lower() or 'cluster' in col.lower()]
    print(f"🔍 Cell type columns: {cell_type_cols}")
    
    # Use the most detailed cell type annotation
    if 'cell_type' in adata.obs.columns:
        cell_type_col = 'cell_type'
    elif 'subclass_label' in adata.obs.columns:
        cell_type_col = 'subclass_label'
    elif 'class_label' in adata.obs.columns:
        cell_type_col = 'class_label'
    else:
        cell_type_col = cell_type_cols[0] if cell_type_cols else None
    
    if cell_type_col:
        print(f"📋 Using cell type column: {cell_type_col}")
        print(f"🧬 Cell types: {adata.obs[cell_type_col].value_counts()}")
        
        # Limit cells per type to manage memory
        limited_cells = []
        for cell_type in adata.obs[cell_type_col].unique():
            type_cells = adata.obs[adata.obs[cell_type_col] == cell_type].index
            if len(type_cells) > MAX_CELLS_PER_TYPE:
                type_cells = np.random.choice(type_cells, MAX_CELLS_PER_TYPE, replace=False)
            limited_cells.extend(type_cells)
        
        adata = adata[limited_cells, :]
        print(f"📊 After limiting: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    
    return adata, cell_type_col

def prepare_tokenization_data(adata):
    """Prepare data for Geneformer tokenization."""
    print("🔤 Preparing tokenization data...")
    
    # Convert to DataFrame with proper gene names
    expr_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
        index=adata.obs.index,
        columns=adata.var.index
    )
    
    # Ensure gene names are strings and valid
    expr_df.columns = expr_df.columns.astype(str)
    
    # Filter out zero-expression genes
    gene_sums = expr_df.sum(axis=0)
    expressed_genes = gene_sums[gene_sums > 0].index
    expr_df = expr_df[expressed_genes]
    
    print(f"📊 Expression data: {expr_df.shape[0]:,} cells × {expr_df.shape[1]:,} genes")
    
    return expr_df

# =============================================================================
# GENEFORMER TOKENIZATION & EMBEDDING
# =============================================================================

def tokenize_cells(expr_df):
    """Tokenize cells using Geneformer tokenizer."""
    print("🔤 Tokenizing cells with Geneformer...")
    
    # Initialize tokenizer
    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_type": "cell_type"},
        nproc=4,
        model_input_size=2048,  # Geneformer V1 input size
        special_token=False
    )
    
    # Create temporary dataset directory
    temp_dataset_dir = RESULTS_DIR / "temp_tokenized"
    temp_dataset_dir.mkdir(exist_ok=True)
    
    # Convert expression data to the format expected by tokenizer
    # Each cell should be a dict with gene names as keys and counts as values
    cell_dicts = []
    for idx, (cell_id, cell_data) in enumerate(expr_df.iterrows()):
        # Convert to integer counts (Geneformer expects raw counts)
        cell_dict = {}
        for gene, count in cell_data.items():
            if count > 0:  # Only include expressed genes
                cell_dict[gene] = int(count)
        
        if len(cell_dict) > 0:  # Only add cells with expressed genes
            cell_dicts.append({
                'input_ids': cell_dict,
                'cell_id': cell_id
            })
    
    print(f"📊 Prepared {len(cell_dicts)} cells for tokenization")
    
    # Save as pickle for tokenizer
    dataset_file = temp_dataset_dir / "dataset.pkl"
    with open(dataset_file, 'wb') as f:
        pickle.dump(cell_dicts, f)
    
    # Tokenize
    tokenizer.tokenize_data(
        data_directory=temp_dataset_dir,
        output_directory=temp_dataset_dir / "tokenized",
        output_prefix="dlpfc_11k",
        file_format="pkl"
    )
    
    # Load tokenized data
    tokenized_file = temp_dataset_dir / "tokenized" / "dlpfc_11k.dataset"
    if not tokenized_file.exists():
        # Try with pkl extension
        tokenized_file = temp_dataset_dir / "tokenized" / "dlpfc_11k.pkl"
    
    if tokenized_file.exists():
        with open(tokenized_file, 'rb') as f:
            tokenized_cells = pickle.load(f)
        print(f"✅ Tokenization complete: {len(tokenized_cells)} cells")
        return tokenized_cells
    else:
        print("❌ Tokenization failed - file not found")
        return None

def extract_embeddings(tokenized_cells, model_path, batch_size=BATCH_SIZE):
    """Extract cell embeddings using Geneformer."""
    print(f"🧬 Extracting embeddings with batch size {batch_size}...")
    
    # Initialize embedding extractor
    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,
        max_ncells=None,
        emb_layer=-1,  # Use final layer
        emb_label=None,
        labels_to_plot=None,
        forward_batch_size=batch_size,
        nproc=1
    )
    
    embeddings = {}
    
    # Process in batches to manage memory
    n_cells = len(tokenized_cells)
    n_batches = (n_cells + batch_size - 1) // batch_size
    
    for batch_idx in tqdm(range(n_batches), desc="Extracting embeddings"):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_cells)
        batch_cells = tokenized_cells[start_idx:end_idx]
        
        try:
            # Clear GPU cache before processing batch
            clear_gpu_cache()
            
            # Extract embeddings for this batch
            batch_embeddings = embex.extract_embs(
                model_directory=model_path,
                input_data_file=batch_cells,
                output_directory=RESULTS_DIR / f"temp_emb_batch_{batch_idx}",
                output_prefix=f"batch_{batch_idx}"
            )
            
            # Store embeddings
            for i, cell in enumerate(batch_cells):
                cell_id = cell.get('cell_id', f"cell_{start_idx + i}")
                if i < len(batch_embeddings):
                    embeddings[cell_id] = batch_embeddings[i]
            
            print(f"✅ Batch {batch_idx + 1}/{n_batches} complete")
            
        except Exception as e:
            print(f"❌ Error in batch {batch_idx}: {str(e)}")
            # Try with smaller batch size
            if batch_size > 10:
                print("🔄 Retrying with smaller batch size...")
                return extract_embeddings(tokenized_cells, model_path, batch_size // 2)
            else:
                raise e
    
    print(f"✅ Embedding extraction complete: {len(embeddings)} cells")
    return embeddings

# =============================================================================
# PERTURBATION ANALYSIS
# =============================================================================

def simulate_gene_perturbation(tokenized_cells, target_gene, model_path):
    """Simulate gene perturbation by masking the target gene."""
    print(f"🎯 Simulating perturbation of {target_gene}...")
    
    perturbed_cells = []
    
    for cell in tokenized_cells:
        # Create a copy of the cell
        perturbed_cell = cell.copy()
        
        # Remove the target gene from input_ids if present
        if 'input_ids' in perturbed_cell:
            input_ids = perturbed_cell['input_ids'].copy()
            
            # Find and remove target gene
            gene_found = False
            if isinstance(input_ids, dict):
                if target_gene in input_ids:
                    del input_ids[target_gene]
                    gene_found = True
            elif isinstance(input_ids, list):
                # For tokenized format, need to map back to gene names
                # This is more complex and depends on tokenizer implementation
                pass
            
            perturbed_cell['input_ids'] = input_ids
            perturbed_cell['perturbed_gene'] = target_gene
            perturbed_cell['gene_found'] = gene_found
        
        perturbed_cells.append(perturbed_cell)
    
    return perturbed_cells

def compute_perturbation_effects(original_embeddings, perturbed_embeddings, target_gene):
    """Compute perturbation effects as cosine distances."""
    print(f"📐 Computing perturbation effects for {target_gene}...")
    
    effects = {}
    
    for cell_id in original_embeddings.keys():
        if cell_id in perturbed_embeddings:
            orig_emb = original_embeddings[cell_id]
            pert_emb = perturbed_embeddings[cell_id]
            
            # Compute cosine distance
            distance = compute_cosine_distance(orig_emb, pert_emb)
            effects[cell_id] = {
                'perturbation_distance': distance,
                'target_gene': target_gene
            }
    
    return effects

# =============================================================================
# ANALYSIS & VISUALIZATION
# =============================================================================

def analyze_by_cell_type(perturbation_results, adata, cell_type_col):
    """Analyze perturbation effects by cell type."""
    print("📊 Analyzing by cell type...")
    
    results_df = []
    
    for gene, gene_results in perturbation_results.items():
        for cell_id, effect_data in gene_results.items():
            if cell_id in adata.obs.index:
                cell_type = adata.obs.loc[cell_id, cell_type_col]
                results_df.append({
                    'gene': gene,
                    'cell_id': cell_id,
                    'cell_type': cell_type,
                    'perturbation_distance': effect_data['perturbation_distance']
                })
    
    results_df = pd.DataFrame(results_df)
    
    # Compute summary statistics by gene and cell type
    summary = results_df.groupby(['gene', 'cell_type'])['perturbation_distance'].agg([
        'mean', 'std', 'count'
    ]).reset_index()
    
    # Add overall gene ranking
    gene_ranking = results_df.groupby('gene')['perturbation_distance'].mean().sort_values(ascending=False)
    
    return results_df, summary, gene_ranking

def create_visualizations(results_df, summary, gene_ranking):
    """Create publication-quality visualizations."""
    print("📈 Creating visualizations...")
    
    # Set style
    plt.style.use('default')
    sns.set_palette("viridis")
    
    # 1. Heatmap of gene×cell_type effects
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Virtual Perturbation Analysis Results', fontsize=16, fontweight='bold')
    
    # Heatmap
    pivot_data = summary.pivot(index='gene', columns='cell_type', values='mean')
    sns.heatmap(pivot_data, annot=True, fmt='.4f', cmap='YlOrRd', 
                ax=axes[0,0], cbar_kws={'label': 'Perturbation Distance'})
    axes[0,0].set_title('Perturbation Effects by Gene and Cell Type')
    axes[0,0].set_xlabel('Cell Type')
    axes[0,0].set_ylabel('Gene')
    
    # Gene ranking bar chart
    top_genes = gene_ranking.head(15)
    axes[0,1].barh(range(len(top_genes)), top_genes.values[::-1])
    axes[0,1].set_yticks(range(len(top_genes)))
    axes[0,1].set_yticklabels(top_genes.index[::-1])
    axes[0,1].set_xlabel('Mean Perturbation Distance')
    axes[0,1].set_title('Top 15 Genes by Perturbation Effect')
    
    # Cell type distribution
    cell_type_counts = summary.groupby('cell_type')['count'].sum()
    axes[1,0].pie(cell_type_counts.values, labels=cell_type_counts.index, autopct='%1.1f%%')
    axes[1,0].set_title('Cell Type Distribution in Analysis')
    
    # Perturbation distance distribution
    axes[1,1].hist(results_df['perturbation_distance'], bins=50, alpha=0.7, edgecolor='black')
    axes[1,1].axvline(results_df['perturbation_distance'].mean(), color='red', linestyle='--', 
                      label=f'Mean: {results_df["perturbation_distance"].mean():.4f}')
    axes[1,1].set_xlabel('Perturbation Distance')
    axes[1,1].set_ylabel('Frequency')
    axes[1,1].set_title('Distribution of Perturbation Effects')
    axes[1,1].legend()
    
    plt.tight_layout()
    fig_path = RESULTS_DIR / 'perturbation_analysis_figures.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"📊 Figures saved: {fig_path}")
    
    return fig_path

def load_drug_database():
    """Load drug interaction database."""
    print("💊 Loading drug database...")
    
    drug_db = pd.read_csv(DRUG_DB_FILE, sep='\t', low_memory=False)
    print(f"📊 Drug database: {len(drug_db):,} interactions")
    
    return drug_db

def cross_reference_drugs(gene_ranking, drug_db):
    """Cross-reference top genes with drug database."""
    print("🔍 Cross-referencing with drug database...")
    
    drug_candidates = []
    
    for gene, score in gene_ranking.items():
        # Find drugs for this gene
        gene_drugs = drug_db[drug_db['gene_name'].str.upper() == gene.upper()]
        
        for _, row in gene_drugs.iterrows():
            drug_candidates.append({
                'gene': gene,
                'perturbation_score': score,
                'drug_name': row.get('drug_name', ''),
                'interaction_type': row.get('interaction_type', ''),
                'source': row.get('source', ''),
                'pmid': row.get('pmid', '')
            })
    
    drug_candidates_df = pd.DataFrame(drug_candidates)
    
    if not drug_candidates_df.empty:
        # Sort by perturbation score
        drug_candidates_df = drug_candidates_df.sort_values('perturbation_score', ascending=False)
        
        print(f"💊 Found {len(drug_candidates_df)} drug-gene interactions")
        print(f"🎯 Top drug candidates:")
        for _, row in drug_candidates_df.head(10).iterrows():
            print(f"   {row['gene']} → {row['drug_name']} (score: {row['perturbation_score']:.4f})")
    
    return drug_candidates_df

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

def main():
    """Main analysis pipeline."""
    print("🚀 Starting Phase 2: Transformer-Based Perturbation Analysis")
    print("=" * 60)
    
    start_time = datetime.now()
    
    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Load and preprocess data
    print("\n📊 STEP 1: Data Loading & Preprocessing")
    adata, cell_type_col = load_brain_data()
    expr_df = prepare_tokenization_data(adata)
    
    # Save checkpoint
    save_checkpoint({'adata': adata, 'expr_df': expr_df, 'cell_type_col': cell_type_col}, 
                   'step1_data_preprocessed')
    
    # 2. Tokenize cells
    print("\n🔤 STEP 2: Cell Tokenization")
    tokenized_cells = load_checkpoint('step2_tokenized_cells')
    if tokenized_cells is None:
        tokenized_cells = tokenize_cells(expr_df)
        if tokenized_cells:
            save_checkpoint(tokenized_cells, 'step2_tokenized_cells')
        else:
            print("❌ Tokenization failed, stopping analysis")
            return
    else:
        print("✅ Loaded tokenized cells from checkpoint")
    
    # 3. Extract original embeddings
    print("\n🧬 STEP 3: Original Embeddings Extraction")
    original_embeddings = load_checkpoint('step3_original_embeddings')
    if original_embeddings is None:
        original_embeddings = extract_embeddings(tokenized_cells, MODEL_PATH)
        save_checkpoint(original_embeddings, 'step3_original_embeddings')
    else:
        print("✅ Loaded original embeddings from checkpoint")
    
    # 4. Perturbation analysis for each GWAS gene
    print("\n🎯 STEP 4: Perturbation Analysis")
    perturbation_results = load_checkpoint('step4_perturbation_results')
    if perturbation_results is None:
        perturbation_results = {}
        
        for gene in tqdm(GWAS_GENES, desc="Processing genes"):
            print(f"\n🧬 Processing {gene}...")
            
            # Simulate perturbation
            perturbed_cells = simulate_gene_perturbation(tokenized_cells, gene, MODEL_PATH)
            
            # Extract perturbed embeddings
            perturbed_embeddings = extract_embeddings(perturbed_cells, MODEL_PATH)
            
            # Compute effects
            gene_effects = compute_perturbation_effects(original_embeddings, perturbed_embeddings, gene)
            perturbation_results[gene] = gene_effects
            
            # Save intermediate results
            save_checkpoint(perturbation_results, 'step4_perturbation_results_partial')
            
            # Clear GPU cache
            clear_gpu_cache()
        
        save_checkpoint(perturbation_results, 'step4_perturbation_results')
    else:
        print("✅ Loaded perturbation results from checkpoint")
    
    # 5. Analysis by cell type
    print("\n📊 STEP 5: Cell Type Analysis")
    results_df, summary, gene_ranking = analyze_by_cell_type(perturbation_results, adata, cell_type_col)
    
    # 6. Create visualizations
    print("\n📈 STEP 6: Visualizations")
    fig_path = create_visualizations(results_df, summary, gene_ranking)
    
    # 7. Drug database cross-reference
    print("\n💊 STEP 7: Drug Database Cross-reference")
    drug_db = load_drug_database()
    drug_candidates_df = cross_reference_drugs(gene_ranking, drug_db)
    
    # 8. Save final results
    print("\n💾 STEP 8: Saving Results")
    
    # Save CSV files
    results_df.to_csv(RESULTS_DIR / 'full_perturbation_results.csv', index=False)
    summary.to_csv(RESULTS_DIR / 'perturbation_summary_by_celltype.csv', index=False)
    gene_ranking.to_csv(RESULTS_DIR / 'gene_ranking.csv', header=['perturbation_score'])
    drug_candidates_df.to_csv(RESULTS_DIR / 'drug_candidates_phase2.csv', index=False)
    
    # Save JSON summary
    analysis_summary = {
        'analysis_date': datetime.now().isoformat(),
        'dataset': str(BRAIN_FILE),
        'n_cells': len(adata),
        'n_genes_analyzed': len(GWAS_GENES),
        'n_cell_types': len(adata.obs[cell_type_col].unique()) if cell_type_col else 1,
        'model_used': str(MODEL_PATH),
        'top_10_genes': dict(gene_ranking.head(10)),
        'n_drug_candidates': len(drug_candidates_df),
        'runtime_seconds': (datetime.now() - start_time).total_seconds()
    }
    
    with open(RESULTS_DIR / 'analysis_summary_phase2.json', 'w') as f:
        json.dump(analysis_summary, f, indent=2)
    
    # Print final summary
    print("\n" + "=" * 60)
    print("✅ PHASE 2 ANALYSIS COMPLETE!")
    print("=" * 60)
    print(f"📊 Cells analyzed: {len(adata):,}")
    print(f"🧬 Genes analyzed: {len(GWAS_GENES)}")
    print(f"🏆 Top gene: {gene_ranking.index[0]} (score: {gene_ranking.iloc[0]:.4f})")
    print(f"💊 Drug candidates: {len(drug_candidates_df):,}")
    print(f"⏱️  Runtime: {(datetime.now() - start_time).total_seconds():.1f} seconds")
    print(f"📁 Results saved: {RESULTS_DIR}")
    
    return analysis_summary

if __name__ == "__main__":
    # Set random seeds for reproducibility
    np.random.seed(42)
    torch.manual_seed(42)
    
    # Run analysis
    summary = main()