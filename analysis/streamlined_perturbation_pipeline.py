#!/usr/bin/env python3
"""
Streamlined Virtual Perturbation Screening Pipeline
Optimized for proof of concept with better progress reporting
"""
import sys
import os
import gc
import time
import json
import pickle
import warnings
from datetime import datetime
from pathlib import Path

# Set encoding and suppress warnings
os.environ['PYTHONIOENCODING'] = 'utf-8'
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

print("Initializing pipeline...")

# Project paths
PROJ_ROOT = Path(r"D:\openclaw\intelligence-augmentation")
DATA_DIR = PROJ_ROOT / "data"
MODEL_DIR = PROJ_ROOT / "models" / "Geneformer"
RESULTS_DIR = PROJ_ROOT / "analysis" / "results"
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# Data files
BRAIN_DATA = DATA_DIR / "brain_scrna" / "DLPFC_11k.h5ad"
DRUG_DATA = DATA_DIR / "drugs" / "interactions.tsv"
GENE_DICT_DIR = MODEL_DIR / "geneformer" / "gene_dictionaries_30m"

def log_progress(step, message):
    """Log progress with timestamp"""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] {step}: {message}")

def check_environment():
    """Quick environment check"""
    log_progress("SETUP", "Checking environment...")
    
    try:
        import torch
        print(f"  PyTorch: {torch.__version__}")
        print(f"  CUDA available: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"  GPU: {torch.cuda.get_device_name(0)}")
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    except ImportError as e:
        print(f"  Error importing torch: {e}")
        return None

def load_brain_data_minimal(max_cells=500):
    """Load brain data with minimal memory footprint"""
    log_progress("DATA", f"Loading brain data (max {max_cells} cells)...")
    
    try:
        import scanpy as sc
        sc.settings.verbosity = 0  # Quiet
        
        adata = sc.read_h5ad(BRAIN_DATA)
        print(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        
        # Find cell type column
        cell_type_cols = [col for col in adata.obs.columns if 'cell_type' in col.lower()]
        if cell_type_cols:
            ct_col = cell_type_cols[0]
            top_cell_type = adata.obs[ct_col].value_counts().index[0]
            filtered = adata[adata.obs[ct_col] == top_cell_type].copy()
            print(f"  Filtered to '{top_cell_type}': {filtered.n_obs:,} cells")
        else:
            filtered = adata.copy()
            print("  No cell type filtering applied")
        
        # Sample cells
        if filtered.n_obs > max_cells:
            sc.pp.subsample(filtered, n_obs=max_cells, random_state=42)
            print(f"  Subsampled to: {filtered.n_obs:,} cells")
        
        del adata
        gc.collect()
        
        print(f"  Final data: {filtered.n_obs} cells x {filtered.n_vars} genes")
        return filtered
        
    except Exception as e:
        print(f"  Error loading data: {e}")
        return None

def load_gene_dictionaries():
    """Load Geneformer gene mappings"""
    log_progress("GENES", "Loading gene dictionaries...")
    
    try:
        # Check if dictionaries exist
        files = {
            'name_to_id': GENE_DICT_DIR / "gene_name_id_dict_gc30M.pkl",
            'token_dict': GENE_DICT_DIR / "token_dictionary_gc30M.pkl", 
            'medians': GENE_DICT_DIR / "gene_median_dictionary_gc30M.pkl"
        }
        
        for name, path in files.items():
            if not path.exists():
                print(f"  ERROR: Missing {name}: {path}")
                return None, None, None
        
        # Load dictionaries
        with open(files['name_to_id'], 'rb') as f:
            gene_name_to_id = pickle.load(f)
        
        with open(files['token_dict'], 'rb') as f:
            token_dict = pickle.load(f)
        
        with open(files['medians'], 'rb') as f:
            gene_medians = pickle.load(f)
        
        print(f"  Gene name mappings: {len(gene_name_to_id):,}")
        print(f"  Token dictionary: {len(token_dict):,}")
        print(f"  Gene medians: {len(gene_medians):,}")
        
        return gene_name_to_id, token_dict, gene_medians
        
    except Exception as e:
        print(f"  Error loading gene dictionaries: {e}")
        return None, None, None

def get_intelligence_genes():
    """Get core intelligence-associated genes"""
    log_progress("TARGETS", "Defining intelligence gene targets...")
    
    # Focus on highest-confidence genes for proof of concept
    CORE_GENES = [
        "FOXO3", "MEF2C", "COMT", "BDNF", "CADM2", "NEGR1",
        "NRXN1", "SHANK3", "GRIN2A", "DRD2", "SLC6A4",
        "MAPT", "TCF4", "CACNA1C", "SCN1A", "GRIN2B",
        "NLGN1", "HOMER1", "SYN1", "CAMK2A", "CREB1",
        "APP", "APOE", "DISC1", "NRG1", "FMR1"
    ]
    
    print(f"  Core intelligence genes: {len(CORE_GENES)}")
    return CORE_GENES

def tokenize_cells_simple(adata, gene_name_to_id, token_dict, gene_medians):
    """Simple cell tokenization with progress tracking"""
    log_progress("TOKENIZE", "Converting cells to Geneformer tokens...")
    
    try:
        # Get expression matrix
        if hasattr(adata.X, "toarray"):
            X = adata.X.toarray()
        else:
            X = np.array(adata.X)
        
        gene_names = list(adata.var_names)
        print(f"  Expression matrix: {X.shape}")
        
        # Check if genes are already Ensembl IDs or gene symbols
        is_ensembl = gene_names[0].startswith("ENSG") if gene_names else False
        print(f"  Gene format: {'Ensembl IDs' if is_ensembl else 'Gene symbols'}")
        
        if is_ensembl:
            # Genes are already Ensembl IDs - filter to those in token dict
            keep_indices = [i for i, gene in enumerate(gene_names) if gene in token_dict]
            X_filtered = X[:, keep_indices]
            ensembl_ids = [gene_names[i] for i in keep_indices]
            print(f"  Genes in token dictionary: {len(ensembl_ids)}/{len(gene_names)}")
        else:
            # Map gene symbols to Ensembl IDs
            mapped_genes = {}
            for gene in gene_names:
                for variant in [gene, gene.upper()]:
                    if variant in gene_name_to_id:
                        mapped_genes[gene] = gene_name_to_id[variant]
                        break
            
            print(f"  Mapped genes: {len(mapped_genes)}/{len(gene_names)}")
            
            # Keep only mapped genes
            keep_indices = [i for i, gene in enumerate(gene_names) if gene in mapped_genes]
            X_filtered = X[:, keep_indices]
            ensembl_ids = [mapped_genes[gene_names[i]] for i in keep_indices]
        
        print(f"  Filtered expression matrix: {X_filtered.shape}")
        
        # Tokenize cells
        tokenized_cells = []
        for i in range(X_filtered.shape[0]):
            expr_vec = X_filtered[i]
            
            # Get expressed genes
            nonzero_mask = expr_vec > 0
            if not np.any(nonzero_mask):
                continue
            
            expressed_genes = np.array(ensembl_ids)[nonzero_mask]
            expressed_values = expr_vec[nonzero_mask]
            
            # Normalize and rank
            normalized_pairs = []
            for gene, value in zip(expressed_genes, expressed_values):
                median = gene_medians.get(gene, 1.0)
                if median > 0:
                    normalized_pairs.append((gene, value / median))
            
            # Sort by normalized expression
            normalized_pairs.sort(key=lambda x: x[1], reverse=True)
            
            # Convert to tokens
            tokens = []
            for gene, _ in normalized_pairs[:1024]:  # Limit tokens
                if gene in token_dict:
                    tokens.append(token_dict[gene])
            
            if tokens:
                tokenized_cells.append(tokens)
            
            if (i + 1) % 100 == 0:
                print(f"    Tokenized {i+1}/{X_filtered.shape[0]} cells")
        
        avg_tokens = np.mean([len(t) for t in tokenized_cells]) if tokenized_cells else 0
        print(f"  Successfully tokenized: {len(tokenized_cells)} cells")
        print(f"  Average tokens per cell: {avg_tokens:.0f}")
        
        return tokenized_cells, ensembl_ids
        
    except Exception as e:
        print(f"  Error in tokenization: {e}")
        return None, None

def run_simplified_perturbation(tokenized_cells, intelligence_genes, 
                               gene_name_to_id, token_dict):
    """Run simplified perturbation analysis without deep learning model"""
    log_progress("PERTURB", "Running simplified perturbation analysis...")
    
    try:
        # Map intelligence genes to tokens
        gene_tokens = {}
        for gene in intelligence_genes:
            ensembl_id = gene_name_to_id.get(gene) or gene_name_to_id.get(gene.upper())
            if ensembl_id and ensembl_id in token_dict:
                gene_tokens[gene] = token_dict[ensembl_id]
        
        print(f"  Mapped {len(gene_tokens)} intelligence genes to tokens")
        
        results = []
        
        for gene, token in gene_tokens.items():
            # Calculate perturbation impact as token frequency and position
            total_appearances = 0
            position_scores = []
            
            for tokens in tokenized_cells:
                if token in tokens:
                    total_appearances += 1
                    # Position in ranked list (lower = higher expression)
                    position = tokens.index(token)
                    position_score = 1.0 / (position + 1)  # Higher for early positions
                    position_scores.append(position_score)
            
            if position_scores:
                result = {
                    'gene_symbol': gene,
                    'token_id': token,
                    'expression_frequency': total_appearances / len(tokenized_cells),
                    'mean_position_score': np.mean(position_scores),
                    'std_position_score': np.std(position_scores),
                    'max_position_score': np.max(position_scores),
                    'cells_expressing': total_appearances,
                    'perturbation_impact': np.mean(position_scores) * (total_appearances / len(tokenized_cells))
                }
                results.append(result)
        
        # Convert to DataFrame and sort
        results_df = pd.DataFrame(results)
        if len(results_df) > 0:
            results_df = results_df.sort_values('perturbation_impact', ascending=False)
        
        print(f"  Analyzed {len(results_df)} genes")
        return results_df
        
    except Exception as e:
        print(f"  Error in perturbation analysis: {e}")
        return pd.DataFrame()

def analyze_drug_targets(results_df):
    """Cross-reference with drug database"""
    log_progress("DRUGS", "Analyzing drug repurposing candidates...")
    
    try:
        if not DRUG_DATA.exists():
            print(f"  Drug database not found: {DRUG_DATA}")
            return results_df, pd.DataFrame()
        
        drugs_df = pd.read_csv(DRUG_DATA, sep='\t')
        print(f"  Drug database: {len(drugs_df)} interactions")
        
        # Find gene column
        gene_cols = [col for col in drugs_df.columns if 'gene' in col.lower()]
        if not gene_cols:
            print("  No gene column found in drug database")
            return results_df, pd.DataFrame()
        
        gene_col = gene_cols[0]
        drug_col = [col for col in drugs_df.columns if 'drug' in col.lower() or 'compound' in col.lower()][0]
        
        # Map druggability
        drug_targets = set(drugs_df[gene_col].str.upper())
        results_df['druggable'] = results_df['gene_symbol'].str.upper().isin(drug_targets)
        
        # Get drug details for druggable targets
        drug_details = []
        for _, row in results_df[results_df['druggable']].iterrows():
            gene_drugs = drugs_df[drugs_df[gene_col].str.upper() == row['gene_symbol'].upper()]
            for _, drug_row in gene_drugs.iterrows():
                drug_details.append({
                    'gene_symbol': row['gene_symbol'],
                    'perturbation_impact': row['perturbation_impact'],
                    'drug_name': drug_row[drug_col],
                    'interaction_type': drug_row.get('interaction_types', 'unknown')
                })
        
        drug_details_df = pd.DataFrame(drug_details)
        
        n_druggable = results_df['druggable'].sum()
        print(f"  Druggable targets: {n_druggable}")
        print(f"  Drug-gene interactions: {len(drug_details_df)}")
        
        return results_df, drug_details_df
        
    except Exception as e:
        print(f"  Error in drug analysis: {e}")
        return results_df, pd.DataFrame()

def generate_summary(results_df, drug_details_df):
    """Generate final summary report"""
    log_progress("SUMMARY", "Generating summary report...")
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n_genes = len(results_df)
    n_druggable = results_df['druggable'].sum() if 'druggable' in results_df.columns else 0
    
    summary = {
        'analysis_date': timestamp,
        'genes_analyzed': int(n_genes),
        'druggable_targets': int(n_druggable),
        'drug_interactions': int(len(drug_details_df) if drug_details_df is not None else 0),
        'top_target': results_df.iloc[0]['gene_symbol'] if len(results_df) > 0 else None,
        'top_impact_score': float(results_df.iloc[0]['perturbation_impact']) if len(results_df) > 0 else 0
    }
    
    # Save results
    results_file = RESULTS_DIR / "perturbation_results_streamlined.csv"
    results_df.to_csv(results_file, index=False)
    
    if drug_details_df is not None and len(drug_details_df) > 0:
        drug_file = RESULTS_DIR / "drug_candidates_streamlined.csv"
        drug_details_df.to_csv(drug_file, index=False)
    
    summary_file = RESULTS_DIR / "analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"  Results saved to: {RESULTS_DIR}")
    print(f"  Summary: {summary}")
    
    return summary

def main():
    """Main pipeline execution"""
    print("=" * 70)
    print("STREAMLINED VIRTUAL PERTURBATION SCREENING PIPELINE")
    print("=" * 70)
    
    start_time = time.time()
    
    try:
        # Step 1: Environment check
        device = check_environment()
        if device is None:
            raise RuntimeError("Environment check failed")
        
        # Step 2: Load brain data
        adata = load_brain_data_minimal(max_cells=500)
        if adata is None:
            raise RuntimeError("Failed to load brain data")
        
        # Step 3: Load gene dictionaries
        gene_name_to_id, token_dict, gene_medians = load_gene_dictionaries()
        if gene_name_to_id is None:
            raise RuntimeError("Failed to load gene dictionaries")
        
        # Step 4: Get target genes
        intelligence_genes = get_intelligence_genes()
        
        # Step 5: Tokenize cells
        tokenized_cells, ensembl_ids = tokenize_cells_simple(
            adata, gene_name_to_id, token_dict, gene_medians)
        if tokenized_cells is None:
            raise RuntimeError("Failed to tokenize cells")
        
        # Step 6: Run perturbation analysis
        results_df = run_simplified_perturbation(
            tokenized_cells, intelligence_genes, gene_name_to_id, token_dict)
        if len(results_df) == 0:
            raise RuntimeError("No perturbation results generated")
        
        # Step 7: Drug analysis
        results_df, drug_details_df = analyze_drug_targets(results_df)
        
        # Step 8: Generate summary
        summary = generate_summary(results_df, drug_details_df)
        
        runtime_minutes = (time.time() - start_time) / 60
        
        print("\n" + "=" * 70)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 70)
        print(f"Runtime: {runtime_minutes:.1f} minutes")
        print(f"Genes analyzed: {summary['genes_analyzed']}")
        print(f"Druggable targets: {summary['druggable_targets']}")
        print(f"Top target: {summary['top_target']} (score: {summary['top_impact_score']:.4f})")
        print(f"Results directory: {RESULTS_DIR}")
        
        return results_df, drug_details_df, summary
        
    except Exception as e:
        print(f"\nPIPELINE FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None

if __name__ == "__main__":
    results_df, drug_details_df, summary = main()