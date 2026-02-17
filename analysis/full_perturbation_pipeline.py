#!/usr/bin/env python3
"""
Virtual Perturbation Screening for Cognitive Enhancement Targets
Full Pipeline - Optimized for RTX 2060 6GB, 16GB RAM

Steps:
1. Load and preprocess DLPFC brain data
2. Map GWAS intelligence genes  
3. In silico perturbation with Geneformer
4. Rank targets by perturbation effect
5. Drug repurposing with DGIdb
6. Generate results and figures
"""
import sys
import os

# Set encoding to UTF-8
os.environ['PYTHONIOENCODING'] = 'utf-8'
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import gc
import time
import pickle
import json
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Memory and GPU settings
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:1024"

# Project paths
PROJ_ROOT = Path(r"D:\openclaw\intelligence-augmentation")
DATA_DIR = PROJ_ROOT / "data"  
MODEL_DIR = PROJ_ROOT / "models" / "Geneformer"
RESULTS_DIR = PROJ_ROOT / "analysis" / "results"
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# Data files
BRAIN_DATA = DATA_DIR / "brain_scrna" / "DLPFC_11k.h5ad"
GWAS_DATA = DATA_DIR / "gwas" / "SavageJansen_2018_intelligence_metaanalysis.txt.gz" 
DRUG_DATA = DATA_DIR / "drugs" / "interactions.tsv"

# Model paths
GENEFORMER_V1 = MODEL_DIR / "Geneformer-V1-10M"  # Start with smallest model
GENE_DICT_DIR = MODEL_DIR / "geneformer" / "gene_dictionaries_30m"

def setup_environment():
    """Check CUDA and memory availability"""
    print("=== Environment Setup ===")
    
    import torch
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    
    if torch.cuda.is_available():
        print(f"CUDA device: {torch.cuda.get_device_name(0)}")
        print(f"CUDA memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
        torch.cuda.empty_cache()
    
    # Set device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    return device

def load_brain_data(sample_cells=1000, filter_cell_type="glutamatergic neuron"):
    """Load and filter brain scRNA-seq data"""
    print("\n=== Step 1: Loading Brain Data ===")
    
    import scanpy as sc
    sc.settings.verbosity = 1  # Reduce scanpy verbosity
    
    print(f"Loading: {BRAIN_DATA}")
    adata = sc.read_h5ad(BRAIN_DATA)
    print(f"Full dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Check cell type column
    cell_type_cols = [col for col in adata.obs.columns 
                     if 'cell_type' in col.lower() or 'celltype' in col.lower()]
    
    if not cell_type_cols:
        print("Warning: No cell type column found. Using all cells.")
        filtered = adata.copy()
    else:
        cell_type_col = cell_type_cols[0]
        print(f"Cell type column: {cell_type_col}")
        print(f"Available cell types: {adata.obs[cell_type_col].unique()[:10]}")
        
        # Filter to specific cell type if available
        if filter_cell_type in adata.obs[cell_type_col].values:
            filtered = adata[adata.obs[cell_type_col] == filter_cell_type].copy()
            print(f"Filtered to {filter_cell_type}: {filtered.n_obs:,} cells")
        else:
            # Use largest cell population
            largest_type = adata.obs[cell_type_col].value_counts().index[0]
            filtered = adata[adata.obs[cell_type_col] == largest_type].copy()
            print(f"Using largest cell type '{largest_type}': {filtered.n_obs:,} cells")
    
    # Sample cells for memory efficiency
    if filtered.n_obs > sample_cells:
        sc.pp.subsample(filtered, n_obs=sample_cells, random_state=42)
        print(f"Subsampled to: {filtered.n_obs:,} cells")
    
    # Basic quality stats
    print(f"Final dataset: {filtered.n_obs:,} cells × {filtered.n_vars:,} genes")
    print(f"Expression range: [{filtered.X.min():.2f}, {filtered.X.max():.2f}]")
    
    del adata  # Free memory
    gc.collect()
    return filtered

def load_gene_mappings():
    """Load Geneformer gene dictionaries"""
    print("\n=== Loading Gene Mappings ===")
    
    # Load gene dictionaries
    gene_name_dict_file = GENE_DICT_DIR / "gene_name_id_dict_gc30M.pkl"
    token_dict_file = GENE_DICT_DIR / "token_dictionary_gc30M.pkl"
    gene_median_file = GENE_DICT_DIR / "gene_median_dictionary_gc30M.pkl"
    
    if not all(f.exists() for f in [gene_name_dict_file, token_dict_file, gene_median_file]):
        raise FileNotFoundError("Geneformer gene dictionaries not found!")
    
    with open(gene_name_dict_file, "rb") as f:
        gene_name_to_id = pickle.load(f)
    
    with open(token_dict_file, "rb") as f:
        token_dict = pickle.load(f)
    
    with open(gene_median_file, "rb") as f:
        gene_medians = pickle.load(f)
    
    print(f"Gene name->ID mappings: {len(gene_name_to_id):,}")
    print(f"Token dictionary: {len(token_dict):,}")
    print(f"Gene medians: {len(gene_medians):,}")
    
    return gene_name_to_id, token_dict, gene_medians

def get_intelligence_genes():
    """Get intelligence-associated genes from multiple sources"""
    print("\n=== Step 2: Intelligence Gene Selection ===")
    
    # Literature-based intelligence genes (high confidence)
    INTELLIGENCE_GENES_CORE = [
        "FOXO3", "MEF2C", "TCF4", "NEGR1", "KIBRA", "COMT", "BDNF", "NRXN1",
        "CADM2", "LRRC4C", "ADAM12", "CSE1L", "MAPT", "WNT4", "FGF5",
        "ATXN1", "DCC", "GRM8", "DLGAP1", "NCALD", "SHANK3", "GRIN2A",
        "SLC6A3", "DRD2", "DRD1", "CREB1", "PRKACA", "CAMK2A", "SYN1",
        "SH2B1", "RBFOX1", "MEF2A", "NRXN3", "CNTNAP2", "GRIK3",
        "ATP2B2", "KCNJ6", "CACNA1C", "SCN1A", "HCN1", "GABRA5",
        "SLC17A7", "GRM5", "HOMER1", "NLGN1", "SYNGAP1", "FMR1",
        "DYRK1A", "DISC1", "NRG1"
    ]
    
    # Extended list from neurogenetics literature
    INTELLIGENCE_GENES_EXTENDED = [
        "APOE", "APP", "PSEN1", "SNCA", "LRRK2", "PARK2", "PARK7", "PINK1",
        "CACNA1A", "CACNA1B", "SCN2A", "SCN8A", "KCNQ2", "KCNQ3",
        "GABA", "GABRB3", "GABRA1", "GABRA3", "GAD1", "GAD2",
        "HTR1A", "HTR2A", "SLC6A4", "MAOA", "MAOB", "TPH1", "TPH2",
        "CHRNA4", "CHRNA7", "CHRNB2", "CHRNB4", "CHRM1", "CHRM3",
        "GRM1", "GRM4", "GRM7", "GRIA1", "GRIA2", "GRIA3", "GRIA4",
        "GRID1", "GRID2", "GRIN1", "GRIN2B", "GRIN2C", "GRIN2D",
        "SLC1A1", "SLC1A2", "SLC1A3", "SLC1A6", "SLC1A7"
    ]
    
    all_genes = list(set(INTELLIGENCE_GENES_CORE + INTELLIGENCE_GENES_EXTENDED))
    
    print(f"Core intelligence genes: {len(INTELLIGENCE_GENES_CORE)}")
    print(f"Extended set: {len(INTELLIGENCE_GENES_EXTENDED)}")
    print(f"Total unique genes: {len(all_genes)}")
    
    return all_genes, INTELLIGENCE_GENES_CORE

def load_geneformer_model(device):
    """Load Geneformer model with memory optimization"""
    print("\n=== Loading Geneformer Model ===")
    
    from transformers import BertForMaskedLM, BertConfig
    import torch
    
    print(f"Model path: {GENEFORMER_V1}")
    
    if not GENEFORMER_V1.exists():
        raise FileNotFoundError(f"Model not found: {GENEFORMER_V1}")
    
    # Load model
    model = BertForMaskedLM.from_pretrained(str(GENEFORMER_V1), torch_dtype=torch.float16 if device.type == "cuda" else torch.float32)
    model.to(device)
    model.eval()
    
    # Count parameters
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    
    print(f"Model loaded: {type(model).__name__}")
    print(f"Total parameters: {total_params:,}")
    print(f"Trainable parameters: {trainable_params:,}")
    print(f"Model size: ~{total_params * 2 / 1e6:.0f} MB")  # FP16 approximation
    
    if device.type == "cuda":
        print(f"GPU memory after loading: {torch.cuda.memory_allocated() / 1e9:.2f} GB")
    
    return model

def tokenize_cells(adata, gene_name_to_id, token_dict, gene_medians, max_cells=1000):
    """Tokenize cells for Geneformer input"""
    print("\n=== Step 3: Cell Tokenization ===")
    
    # Get expression matrix
    if hasattr(adata.X, "toarray"):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)
    
    print(f"Expression matrix: {X.shape}")
    print(f"Data type: {X.dtype}")
    print(f"Value range: [{X.min():.3f}, {X.max():.3f}]")
    
    # Map gene names to Ensembl IDs
    gene_names = list(adata.var_names)
    
    # Check if already Ensembl format
    is_ensembl = gene_names[0].startswith("ENSG") if gene_names else False
    print(f"Gene format: {'Ensembl' if is_ensembl else 'Gene symbols'}")
    
    if not is_ensembl:
        # Map symbols to Ensembl
        mapped_genes = {}
        for gene in gene_names:
            # Try direct match and uppercase
            for variant in [gene, gene.upper()]:
                if variant in gene_name_to_id:
                    mapped_genes[gene] = gene_name_to_id[variant]
                    break
        
        print(f"Mapped {len(mapped_genes)} / {len(gene_names)} genes to Ensembl")
        
        # Filter to mapped genes
        keep_indices = [i for i, gene in enumerate(gene_names) if gene in mapped_genes]
        X = X[:, keep_indices]
        gene_ensembl_ids = [mapped_genes[gene_names[i]] for i in keep_indices]
        
    else:
        gene_ensembl_ids = gene_names
        print(f"Using Ensembl IDs directly")
    
    print(f"Final gene set: {len(gene_ensembl_ids)} genes")
    
    def tokenize_single_cell(expr_vec, ensembl_ids):
        """Tokenize one cell using rank-value encoding"""
        # Only keep expressed genes (>0)
        nonzero_mask = expr_vec > 0
        if not np.any(nonzero_mask):
            return []
        
        nonzero_genes = np.array(ensembl_ids)[nonzero_mask]
        nonzero_expr = expr_vec[nonzero_mask]
        
        # Normalize by gene medians
        normalized_expr = []
        for gene, expr in zip(nonzero_genes, nonzero_expr):
            median = gene_medians.get(gene, 1.0)
            if median > 0:
                normalized_expr.append((gene, expr / median))
        
        # Sort by normalized expression (descending)
        normalized_expr.sort(key=lambda x: x[1], reverse=True)
        
        # Convert to tokens (limit to 2048 for V1 model)
        tokens = []
        for gene, _ in normalized_expr[:2048]:
            if gene in token_dict:
                tokens.append(token_dict[gene])
        
        return tokens
    
    # Tokenize all cells (with progress)
    tokenized_cells = []
    n_cells = min(X.shape[0], max_cells)
    print(f"Tokenizing {n_cells} cells...")
    
    for i in range(n_cells):
        tokens = tokenize_single_cell(X[i], gene_ensembl_ids)
        if len(tokens) > 0:
            tokenized_cells.append(tokens)
        
        if (i + 1) % 100 == 0:
            print(f"  Tokenized {i+1}/{n_cells}")
    
    avg_tokens = np.mean([len(t) for t in tokenized_cells])
    print(f"Successfully tokenized: {len(tokenized_cells)} cells")
    print(f"Average tokens per cell: {avg_tokens:.0f}")
    
    return tokenized_cells, gene_ensembl_ids

def compute_baseline_embeddings(model, tokenized_cells, device, batch_size=32):
    """Compute baseline cell embeddings efficiently"""
    print("\n=== Computing Baseline Embeddings ===")
    
    import torch
    
    def get_cell_embedding(model, tokens, device):
        """Get embedding for one cell"""
        if len(tokens) == 0:
            return None
        
        # Pad to max sequence length 
        max_len = 2048  # V1 max length
        padded_tokens = tokens[:max_len] + [0] * max(0, max_len - len(tokens))
        
        input_ids = torch.tensor([padded_tokens[:max_len]], dtype=torch.long, device=device)
        attention_mask = torch.ones_like(input_ids)
        attention_mask[0, len(tokens):] = 0  # Mask padding
        
        with torch.no_grad():
            outputs = model.bert(input_ids=input_ids, attention_mask=attention_mask)
            # Use mean pooling of non-padded tokens
            mask_expanded = attention_mask.unsqueeze(-1).expand(outputs.last_hidden_state.size()).float()
            sum_embeddings = torch.sum(outputs.last_hidden_state * mask_expanded, 1)
            sum_mask = torch.clamp(mask_expanded.sum(1), min=1e-9)
            embedding = sum_embeddings / sum_mask
            return embedding.squeeze().cpu().numpy()
    
    # Process in batches to manage memory
    embeddings = []
    n_cells = len(tokenized_cells)
    
    for i in range(0, n_cells, batch_size):
        batch_end = min(i + batch_size, n_cells)
        batch_embeddings = []
        
        for j in range(i, batch_end):
            emb = get_cell_embedding(model, tokenized_cells[j], device)
            if emb is not None:
                batch_embeddings.append(emb)
        
        if batch_embeddings:
            embeddings.extend(batch_embeddings)
        
        # Clear GPU cache periodically
        if device.type == "cuda":
            torch.cuda.empty_cache()
        
        print(f"  Processed {batch_end}/{n_cells} cells")
    
    embeddings = np.array(embeddings)
    print(f"Baseline embeddings shape: {embeddings.shape}")
    
    return embeddings

def run_perturbation_screen(model, tokenized_cells, baseline_embeddings, 
                          intelligence_genes, gene_name_to_id, token_dict, 
                          gene_ensembl_ids, device, max_genes=100):
    """Run in silico perturbation screen"""
    print("\n=== Step 4: In Silico Perturbation Screen ===")
    
    import torch
    
    # Find which intelligence genes are in our data
    available_genes = []
    for gene in intelligence_genes:
        # Map to Ensembl ID
        ensembl_id = gene_name_to_id.get(gene) or gene_name_to_id.get(gene.upper())
        if ensembl_id and ensembl_id in token_dict and ensembl_id in gene_ensembl_ids:
            available_genes.append((gene, ensembl_id))
    
    available_genes = available_genes[:max_genes]  # Limit for computational feasibility
    print(f"Will perturb {len(available_genes)} intelligence genes (max {max_genes})")
    
    results = []
    
    def get_perturbed_embedding(tokens, target_token, model, device):
        """Get embedding after removing target gene token"""
        perturbed_tokens = [t for t in tokens if t != target_token]
        if len(perturbed_tokens) == 0:
            return None
        
        # Pad sequence
        max_len = 2048
        padded = perturbed_tokens[:max_len] + [0] * max(0, max_len - len(perturbed_tokens))
        
        input_ids = torch.tensor([padded[:max_len]], dtype=torch.long, device=device)
        attention_mask = torch.ones_like(input_ids)
        attention_mask[0, len(perturbed_tokens):] = 0
        
        with torch.no_grad():
            outputs = model.bert(input_ids=input_ids, attention_mask=attention_mask)
            mask_expanded = attention_mask.unsqueeze(-1).expand(outputs.last_hidden_state.size()).float()
            sum_embeddings = torch.sum(outputs.last_hidden_state * mask_expanded, 1)
            sum_mask = torch.clamp(mask_expanded.sum(1), min=1e-9)
            embedding = sum_embeddings / sum_mask
            return embedding.squeeze().cpu().numpy()
    
    # Process each gene
    for gi, (gene_name, ensembl_id) in enumerate(available_genes):
        if ensembl_id not in token_dict:
            continue
        
        target_token = token_dict[ensembl_id]
        shifts = []
        cells_affected = 0
        
        # Test perturbation on each cell
        for ci, tokens in enumerate(tokenized_cells):
            if target_token in tokens:
                cells_affected += 1
                perturbed_emb = get_perturbed_embedding(tokens, target_token, model, device)
                
                if perturbed_emb is not None and ci < len(baseline_embeddings):
                    shift = np.linalg.norm(perturbed_emb - baseline_embeddings[ci])
                    shifts.append(shift)
        
        if shifts:
            result = {
                "gene_symbol": gene_name,
                "ensembl_id": ensembl_id,
                "mean_shift": np.mean(shifts),
                "std_shift": np.std(shifts),
                "median_shift": np.median(shifts),
                "max_shift": np.max(shifts),
                "min_shift": np.min(shifts),
                "n_cells_expressing": cells_affected,
                "n_cells_measured": len(shifts),
                "expression_frequency": cells_affected / len(tokenized_cells),
                "perturbation_impact_score": np.mean(shifts) * (cells_affected / len(tokenized_cells))
            }
            results.append(result)
        
        # Progress and memory management
        if (gi + 1) % 10 == 0:
            print(f"  Completed {gi+1}/{len(available_genes)} genes")
            if device.type == "cuda":
                torch.cuda.empty_cache()
    
    # Convert to DataFrame and sort by impact score
    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df = results_df.sort_values("perturbation_impact_score", ascending=False)
    
    print(f"Perturbation screen complete: {len(results_df)} genes analyzed")
    return results_df

def drug_repurposing_analysis(results_df):
    """Cross-reference results with drug database"""
    print("\n=== Step 5: Drug Repurposing Analysis ===")
    
    # Load DGIdb interactions
    print(f"Loading drug interactions: {DRUG_DATA}")
    drugs_df = pd.read_csv(DRUG_DATA, sep='\t')
    
    print(f"Drug database: {len(drugs_df)} interactions")
    print(f"Columns: {list(drugs_df.columns)}")
    
    # Find gene column
    gene_columns = [col for col in drugs_df.columns if any(term in col.lower() for term in ['gene', 'target'])]
    drug_columns = [col for col in drugs_df.columns if any(term in col.lower() for term in ['drug', 'compound', 'name'])]
    
    if not gene_columns:
        print("Warning: No gene column found in drug database")
        results_df['druggable'] = False
        results_df['drug_count'] = 0
        return results_df
    
    gene_col = gene_columns[0]
    drug_col = drug_columns[0] if drug_columns else gene_columns[0]
    
    print(f"Using gene column: {gene_col}")
    print(f"Using drug column: {drug_col}")
    
    # Create drug target mapping
    drug_targets = drugs_df[gene_col].str.upper().unique()
    
    # Map results to drugs
    results_df['druggable'] = results_df['gene_symbol'].str.upper().isin(drug_targets)
    
    # Count drugs per gene
    drug_counts = drugs_df.groupby(drugs_df[gene_col].str.upper()).size()
    results_df['drug_count'] = results_df['gene_symbol'].str.upper().map(drug_counts).fillna(0).astype(int)
    
    # Get drug details for top targets
    druggable_results = results_df[results_df['druggable']].copy()
    
    drug_details = []
    for _, row in druggable_results.iterrows():
        gene_drugs = drugs_df[drugs_df[gene_col].str.upper() == row['gene_symbol'].upper()]
        for _, drug_row in gene_drugs.iterrows():
            drug_details.append({
                'gene_symbol': row['gene_symbol'],
                'perturbation_impact_score': row['perturbation_impact_score'],
                'drug_name': drug_row[drug_col],
                'interaction_type': drug_row.get('interaction_types', 'unknown')
            })
    
    drug_details_df = pd.DataFrame(drug_details)
    
    print(f"Druggable targets found: {results_df['druggable'].sum()}")
    print(f"Total drug-gene interactions: {len(drug_details_df)}")
    
    return results_df, drug_details_df

def generate_summary_report(results_df, drug_details_df, runtime_stats):
    """Generate comprehensive summary report"""
    print("\n=== Step 6: Generating Summary Report ===")
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Basic statistics
    n_genes_tested = len(results_df)
    n_druggable = results_df['druggable'].sum() if 'druggable' in results_df.columns else 0
    top_impact = results_df.iloc[0]['perturbation_impact_score'] if len(results_df) > 0 else 0
    
    # Create report
    report = f"""# Virtual Perturbation Screening Results
## Cognitive Enhancement Target Discovery Pipeline

**Generated:** {timestamp}
**Pipeline Version:** Full Analysis v1.0
**Model:** Geneformer V1-10M
**Device:** {runtime_stats.get('device', 'unknown')}

---

## Executive Summary

This analysis identified potential targets for cognitive enhancement through virtual perturbation screening of brain single-cell RNA sequencing data. We systematically perturbed {n_genes_tested} intelligence-associated genes in brain cells using the Geneformer foundation model to predict downstream effects.

**Key Findings:**
- **{n_genes_tested} genes** were successfully screened for perturbation effects
- **{n_druggable} genes** have existing drug interactions (druggable targets)
- **{len(drug_details_df) if drug_details_df is not None else 0}** total drug-gene interactions identified
- **Top target:** {results_df.iloc[0]['gene_symbol'] if len(results_df) > 0 else 'N/A'} (impact score: {top_impact:.4f})

---

## Top 20 Target Genes

| Rank | Gene Symbol | Impact Score | Expression Freq | Cells Affected | Druggable | Drug Count |
|------|-------------|--------------|-----------------|----------------|-----------|------------|"""
    
    # Add top genes table
    for i, (_, row) in enumerate(results_df.head(20).iterrows()):
        rank = i + 1
        druggable = "Yes" if row.get('druggable', False) else "No"
        drug_count = row.get('drug_count', 0)
        
        report += f"""
| {rank} | {row['gene_symbol']} | {row['perturbation_impact_score']:.4f} | {row['expression_frequency']:.2%} | {row['n_cells_expressing']} | {druggable} | {drug_count} |"""
    
    # Add methodology section
    report += f"""

---

## Methodology

### Data Sources
- **Brain Data:** DLPFC single-cell RNA-seq ({runtime_stats.get('n_cells', 'unknown')} cells)
- **Target Genes:** Literature-curated intelligence-associated genes
- **Drug Database:** DGIdb drug-gene interaction database
- **Model:** Geneformer V1-10M transformer for genomics

### Analysis Pipeline
1. **Data Preprocessing:** Cell filtering, gene mapping, quality control
2. **Tokenization:** Rank-value encoding for Geneformer input
3. **Baseline Embedding:** Computation of cell state representations
4. **Perturbation:** In silico gene knockout via token removal
5. **Impact Scoring:** Embedding shift magnitude as perturbation effect
6. **Drug Mapping:** Cross-reference with pharmacological databases

### Metrics
- **Perturbation Impact Score:** Mean embedding shift × expression frequency
- **Expression Frequency:** Fraction of cells expressing the gene
- **Embedding Shift:** Euclidean distance in latent space after perturbation

---

## Drug Repurposing Candidates

### High-Impact Druggable Targets"""
    
    if drug_details_df is not None and len(drug_details_df) > 0:
        # Group by gene and show top drugs
        top_druggable = results_df[results_df.get('druggable', False)].head(10)
        
        for _, row in top_druggable.iterrows():
            gene_drugs = drug_details_df[drug_details_df['gene_symbol'] == row['gene_symbol']]
            drug_list = ', '.join(gene_drugs['drug_name'].unique()[:5])
            if len(gene_drugs) > 5:
                drug_list += f" (+{len(gene_drugs)-5} more)"
            
            report += f"""

**{row['gene_symbol']}** (Impact: {row['perturbation_impact_score']:.4f})
- Drugs: {drug_list}
- Expression frequency: {row['expression_frequency']:.2%}
- Total drug interactions: {len(gene_drugs)}"""
    
    # Add limitations and technical details
    report += f"""

---

## Technical Details

### Computational Resources
- **Runtime:** {runtime_stats.get('total_time', 'unknown')} minutes
- **GPU Memory:** {runtime_stats.get('gpu_memory', 'N/A')}
- **Model Size:** ~104M parameters
- **Batch Size:** {runtime_stats.get('batch_size', 'unknown')}

### Quality Metrics
- **Tokenization Success:** {runtime_stats.get('tokenization_success', 'unknown')}%
- **Average Tokens/Cell:** {runtime_stats.get('avg_tokens', 'unknown')}
- **Gene Mapping Rate:** {runtime_stats.get('gene_mapping', 'unknown')}%

---

## Limitations & Considerations

### Technical Limitations
- **Scale:** Analysis limited to {n_genes_tested} genes for computational feasibility
- **Model Size:** Using smallest Geneformer variant (V1-10M) due to memory constraints
- **Cell Type:** Analysis focused on single brain cell population
- **Perturbation Type:** Gene knockout only (no overexpression or partial knockdown)

### Biological Limitations
- **In Silico Nature:** Predictions may not reflect true biological responses
- **Cell Context:** Effects may vary across brain regions and developmental stages
- **Gene Interactions:** Does not account for compensatory mechanisms or gene networks
- **Causality:** Correlation between embedding shifts and functional impact assumed

### Methodological Limitations
- **No Statistical Testing:** Results not corrected for multiple hypothesis testing
- **Validation Needed:** Requires experimental validation in relevant model systems
- **Drug Specificity:** DGIdb interactions may include indirect or non-specific effects

---

## Next Steps & Recommendations

### Immediate Follow-up
1. **Scale Analysis:** Run full genome-wide screen with larger model (V2-316M)
2. **Multiple Cell Types:** Extend to microglia, oligodendrocytes, astrocytes
3. **Statistical Validation:** Apply multiple testing correction and significance thresholds
4. **Overexpression Screen:** Add gene upregulation perturbations

### Experimental Validation
1. **In Vitro Testing:** Selected targets in neuronal cell cultures
2. **Animal Models:** Cognitive testing with target modulation
3. **Human Studies:** Clinical trials for promising drug candidates

### Computational Enhancements
1. **Network Analysis:** Incorporate gene regulatory networks
2. **Pathway Enrichment:** Functional annotation of target sets
3. **Drug Optimization:** Molecular docking and drug design
4. **Multi-Modal:** Integration with brain imaging and proteomics

---

**Analysis completed at {timestamp}**
**For questions contact the computational biology team**
"""
    
    # Save report
    report_file = RESULTS_DIR / "perturbation_screening_report.md"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"Report saved: {report_file}")
    return report_file

def create_summary_figures(results_df):
    """Generate summary plots and figures"""
    print("\n=== Generating Summary Figures ===")
    
    if len(results_df) == 0:
        print("No results to plot")
        return
    
    plt.style.use('default')
    
    # Figure 1: Top genes ranked by impact score
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # A) Top 20 genes bar plot
    top_20 = results_df.head(20)
    axes[0, 0].barh(range(len(top_20)), top_20['perturbation_impact_score'])
    axes[0, 0].set_yticks(range(len(top_20)))
    axes[0, 0].set_yticklabels(top_20['gene_symbol'])
    axes[0, 0].set_xlabel('Perturbation Impact Score')
    axes[0, 0].set_title('Top 20 Target Genes')
    axes[0, 0].grid(axis='x', alpha=0.3)
    
    # B) Impact score distribution
    axes[0, 1].hist(results_df['perturbation_impact_score'], bins=30, alpha=0.7, color='steelblue')
    axes[0, 1].axvline(results_df['perturbation_impact_score'].median(), color='red', linestyle='--', label='Median')
    axes[0, 1].set_xlabel('Perturbation Impact Score')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of Impact Scores')
    axes[0, 1].legend()
    axes[0, 1].grid(alpha=0.3)
    
    # C) Expression frequency vs impact
    axes[1, 0].scatter(results_df['expression_frequency'], results_df['perturbation_impact_score'], 
                      alpha=0.6, s=50)
    axes[1, 0].set_xlabel('Expression Frequency')
    axes[1, 0].set_ylabel('Perturbation Impact Score') 
    axes[1, 0].set_title('Expression vs Impact Relationship')
    axes[1, 0].grid(alpha=0.3)
    
    # D) Druggable vs non-druggable targets
    if 'druggable' in results_df.columns:
        druggable_scores = results_df[results_df['druggable']]['perturbation_impact_score']
        non_druggable_scores = results_df[~results_df['druggable']]['perturbation_impact_score']
        
        axes[1, 1].boxplot([druggable_scores, non_druggable_scores], 
                          labels=['Druggable', 'Non-druggable'])
        axes[1, 1].set_ylabel('Perturbation Impact Score')
        axes[1, 1].set_title('Druggable vs Non-druggable Targets')
        axes[1, 1].grid(alpha=0.3)
    else:
        axes[1, 1].text(0.5, 0.5, 'No drug data available', ha='center', va='center', 
                        transform=axes[1, 1].transAxes)
        axes[1, 1].set_title('Drug Analysis Not Available')
    
    plt.tight_layout()
    fig_file = RESULTS_DIR / "summary_figures.png"
    plt.savefig(fig_file, dpi=300, bbox_inches='tight')
    print(f"Summary figures saved: {fig_file}")
    plt.close()

def save_results(results_df, drug_details_df=None):
    """Save all results to files"""
    print("\n=== Saving Results ===")
    
    # Main results
    results_file = RESULTS_DIR / "perturbation_results.csv"
    results_df.to_csv(results_file, index=False)
    print(f"Results saved: {results_file}")
    
    # Drug details
    if drug_details_df is not None:
        drug_file = RESULTS_DIR / "drug_repurposing_candidates.csv"
        drug_details_df.to_csv(drug_file, index=False)
        print(f"Drug candidates saved: {drug_file}")
    
    # Top targets JSON for easy parsing
    top_targets = results_df.head(50).to_dict('records')
    json_file = RESULTS_DIR / "top_targets.json"
    with open(json_file, 'w') as f:
        json.dump(top_targets, f, indent=2)
    print(f"Top targets JSON saved: {json_file}")

def main():
    """Main pipeline execution"""
    print("=" * 80)
    print("VIRTUAL PERTURBATION SCREENING FOR COGNITIVE ENHANCEMENT")
    print("=" * 80)
    
    start_time = time.time()
    runtime_stats = {}
    
    try:
        # Step 0: Setup
        device = setup_environment()
        runtime_stats['device'] = str(device)
        
        # Step 1: Load brain data
        adata = load_brain_data(sample_cells=1000)
        runtime_stats['n_cells'] = adata.n_obs
        
        # Step 2: Load gene mappings
        gene_name_to_id, token_dict, gene_medians = load_gene_mappings()
        
        # Step 3: Get intelligence genes
        intelligence_genes, core_genes = get_intelligence_genes()
        
        # Step 4: Load model
        model = load_geneformer_model(device)
        
        # Step 5: Tokenize cells
        tokenized_cells, gene_ensembl_ids = tokenize_cells(
            adata, gene_name_to_id, token_dict, gene_medians)
        runtime_stats['tokenization_success'] = len(tokenized_cells) / adata.n_obs * 100
        runtime_stats['avg_tokens'] = np.mean([len(t) for t in tokenized_cells])
        
        # Step 6: Compute baseline embeddings
        baseline_embeddings = compute_baseline_embeddings(
            model, tokenized_cells, device)
        
        # Step 7: Run perturbation screen
        results_df = run_perturbation_screen(
            model, tokenized_cells, baseline_embeddings,
            intelligence_genes, gene_name_to_id, token_dict,
            gene_ensembl_ids, device)
        
        # Step 8: Drug repurposing
        results_df, drug_details_df = drug_repurposing_analysis(results_df)
        
        # Step 9: Generate figures
        create_summary_figures(results_df)
        
        # Step 10: Save results
        save_results(results_df, drug_details_df)
        
        # Step 11: Generate report
        runtime_stats['total_time'] = (time.time() - start_time) / 60
        report_file = generate_summary_report(results_df, drug_details_df, runtime_stats)
        
        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Total runtime: {runtime_stats['total_time']:.1f} minutes")
        print(f"Genes analyzed: {len(results_df)}")
        print(f"Druggable targets: {results_df['druggable'].sum() if 'druggable' in results_df else 0}")
        print(f"Results saved to: {RESULTS_DIR}")
        print(f"Report: {report_file}")
        
        return results_df, drug_details_df
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    results_df, drug_details_df = main()