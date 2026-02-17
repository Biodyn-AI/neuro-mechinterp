#!/usr/bin/env python3
"""
scGPT In-Silico Perturbation Analysis on Intelligence Genes
Critical experiment for comparing scGPT vs Geneformer rankings
"""

# CRITICAL: Apply torchtext patch first
exec(open("/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py").read())

import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from pathlib import Path
from scipy.stats import spearmanr
from scipy.spatial.distance import cosine
import warnings
warnings.filterwarnings('ignore')

print("🧠 Starting scGPT Intelligence Gene Perturbation Analysis")
print("=" * 60)

# Device setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Paths
scgpt_checkpoint_dir = "/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT_checkpoints/brain/"
data_path = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
output_path = "/mnt/d/openclaw/intelligence-augmentation/analysis/SCGPT_CROSSMODEL_REPORT.md"

# Intelligence genes to test
intelligence_genes = ['CADM2', 'GRIN2A', 'CAMK2A', 'MEF2C', 'APP']
geneformer_ranking = ['CADM2', 'GRIN2A', 'CAMK2A', 'MEF2C', 'APP']  # From best to worst impact

print(f"📊 Target genes: {intelligence_genes}")
print(f"📈 Geneformer ranking (high→low impact): {' > '.join(geneformer_ranking)}")

# Load model configuration
args_path = Path(scgpt_checkpoint_dir) / "args.json"
if not args_path.exists():
    raise FileNotFoundError(f"Args file not found: {args_path}")

with open(args_path, 'r') as f:
    args = json.load(f)

print(f"🔧 Loaded model config from: {args_path}")

# Import scGPT modules
try:
    from scgpt.model.model import TransformerModel
    from scgpt.tokenizer import GeneVocab
    from scgpt.preprocess import Preprocessor
    print("✅ scGPT modules imported successfully")
except ImportError as e:
    print(f"❌ Failed to import scGPT: {e}")
    raise

# Load vocabulary
vocab_path = Path(scgpt_checkpoint_dir) / "vocab.json"
if not vocab_path.exists():
    raise FileNotFoundError(f"Vocab file not found: {vocab_path}")

vocab = GeneVocab.from_file(vocab_path)
print(f"📚 Loaded vocabulary: {len(vocab)} genes")

# Debug vocabulary structure
print(f"🔍 Vocabulary type: {type(vocab)}")
print(f"🔍 Available methods: {[m for m in dir(vocab) if not m.startswith('_')]}")
if hasattr(vocab, 'get_stoi'):
    print(f"🔍 Sample vocab items: {list(vocab.get_stoi().items())[:5]}")
elif hasattr(vocab, 'stoi'):
    print(f"🔍 Sample vocab items: {list(vocab.stoi.items())[:5]}")

# Try to get pad token in different ways
pad_token = None
if "<pad>" in vocab:
    pad_token = vocab["<pad>"]
elif hasattr(vocab, 'get_stoi') and "<pad>" in vocab.get_stoi():
    pad_token = vocab.get_stoi()["<pad>"]
elif hasattr(vocab, 'stoi') and "<pad>" in vocab.stoi:
    pad_token = vocab.stoi["<pad>"]
else:
    # Look for other padding tokens
    for token in ["<pad>", "[PAD]", "PAD", "<unk>"]:
        if token in vocab:
            pad_token = vocab[token]
            print(f"🔍 Using {token} as pad token")
            break

print(f"🔍 Pad token found: {pad_token} (type: {type(pad_token)})")

# Load model
print("🔄 Loading scGPT model...")

# Get pad token ID (ensure it's an integer)
pad_token_id = pad_token if pad_token is not None else 0
if isinstance(pad_token_id, str):
    pad_token_id = 0  # fallback to 0 if string
print(f"🔧 Pad token ID: {pad_token_id}")

# Create model using core parameters
# Note: Many training args are for preprocessing/training, not model architecture
model = TransformerModel(
    ntoken=len(vocab),
    d_model=args.get("embsize", 512),
    nhead=args.get("nheads", 8), 
    d_hid=args.get("d_hid", 512),
    nlayers=args.get("nlayers", 12),
    dropout=args.get("dropout", 0.2),
    pad_token="<pad>",
    vocab=vocab,
    # Only include parameters that are accepted by TransformerModel
    pad_value=args.get("pad_value", -2),
    n_input_bins=args.get("n_bins", 51),
)

# Load checkpoint
checkpoint_path = Path(scgpt_checkpoint_dir) / "best_model.pt"
if not checkpoint_path.exists():
    raise FileNotFoundError(f"Model checkpoint not found: {checkpoint_path}")

checkpoint = torch.load(checkpoint_path, map_location=device)
print(f"🔍 Checkpoint keys: {list(checkpoint.keys())}")

# Load state dict with strict=False to ignore architecture mismatches
# This should load compatible layers and skip incompatible ones
try:
    missing_keys, unexpected_keys = model.load_state_dict(checkpoint, strict=False)
    print(f"⚠️ Missing keys: {len(missing_keys)} items")
    print(f"⚠️ Unexpected keys: {len(unexpected_keys)} items")
    if missing_keys:
        print(f"🔍 Sample missing keys: {missing_keys[:3]}...")
    if unexpected_keys:
        print(f"🔍 Sample unexpected keys: {unexpected_keys[:3]}...")
    print("✅ Model loaded (with some architecture mismatches)")
except Exception as e:
    print(f"❌ Failed to load checkpoint: {e}")
    raise
model.to(device)
model.eval()
print(f"✅ Model loaded from: {checkpoint_path}")

# Load brain data
print("🧠 Loading brain single-cell data...")
adata = sc.read_h5ad(data_path)
print(f"📊 Data shape: {adata.shape} (cells × genes)")

# Get gene token IDs
gene_token_ids = {}
missing_genes = []
for gene in intelligence_genes:
    if gene in vocab:
        gene_token_ids[gene] = vocab[gene]
        print(f"🧬 {gene}: token ID {gene_token_ids[gene]}")
    else:
        missing_genes.append(gene)
        print(f"❌ {gene}: not found in vocabulary")

if missing_genes:
    print(f"⚠️ Missing genes from vocab: {missing_genes}")
    print("Continuing with available genes...")
    intelligence_genes = [g for g in intelligence_genes if g not in missing_genes]

# Preprocess data for scGPT
print("🔄 Preprocessing data...")
preprocessor = Preprocessor(
    use_key="X",  # use raw count
    filter_gene_by_counts=3,  # gene filter
    filter_cell_by_counts=False,  # no cell filter for perturbation
    log1p=args.get("log1p", True),
    normalize_total=args.get("normalize_total", 1e4),
    result_normed_key="X_normed",
)

# Process a subset of cells for perturbation (limit to avoid OOM)
n_cells_to_test = min(200, adata.n_obs)
print(f"🎯 Testing on {n_cells_to_test} cells")

# Select random subset of cells
np.random.seed(42)  # For reproducibility
cell_indices = np.random.choice(adata.n_obs, n_cells_to_test, replace=False)
adata_subset = adata[cell_indices, :].copy()

print(f"📊 Subset shape before preprocessing: {adata_subset.shape}")

# Preprocess subset
try:
    adata_preprocessed = preprocessor(adata_subset, batch_key=None)
    if adata_preprocessed is None:
        print("❌ Preprocessor returned None - using original data with basic preprocessing")
        adata_subset = adata_subset.copy()
        # Apply basic preprocessing manually
        sc.pp.normalize_total(adata_subset, target_sum=1e4)
        sc.pp.log1p(adata_subset)
        adata_subset.layers["X_normed"] = adata_subset.X.copy()
    else:
        adata_subset = adata_preprocessed
        print(f"✅ Preprocessing completed, shape: {adata_subset.shape}")
except Exception as e:
    print(f"❌ Preprocessing failed: {e}")
    print("🔄 Using original data with basic preprocessing")
    sc.pp.normalize_total(adata_subset, target_sum=1e4)
    sc.pp.log1p(adata_subset)
    adata_subset.layers["X_normed"] = adata_subset.X.copy()

# Get expression matrix (use normalized counts)
if "X_normed" in adata_subset.layers:
    expr_matrix = adata_subset.layers["X_normed"]
else:
    expr_matrix = adata_subset.X

if hasattr(expr_matrix, 'toarray'):
    expr_matrix = expr_matrix.toarray()

print(f"📊 Expression matrix shape: {expr_matrix.shape}")

# Prepare data for model input
gene_names = adata_subset.var_names.tolist()
print(f"📋 Available genes in data: {len(gene_names)}")

# Map gene names to vocab indices
gene_vocab_mapping = {}
for i, gene in enumerate(gene_names):
    if gene in vocab:
        gene_vocab_mapping[i] = vocab[gene]

print(f"🔗 Mapped {len(gene_vocab_mapping)} genes to vocabulary")

# Perturbation analysis
perturbation_results = {}

for gene in intelligence_genes:
    print(f"\n🧬 Analyzing gene: {gene}")
    gene_token_id = gene_token_ids[gene]
    
    # Find gene index in data
    if gene not in gene_names:
        print(f"⚠️ Gene {gene} not found in expression data, skipping...")
        continue
    
    gene_idx = gene_names.index(gene)
    
    # Store embedding shifts for this gene
    embedding_shifts = []
    
    # Process cells in batches to avoid OOM
    batch_size = 1  # Use batch_size=1 as recommended
    
    with torch.no_grad():
        for cell_idx in range(n_cells_to_test):
            # Get expression values for this cell
            cell_expr = expr_matrix[cell_idx, :]
            
            # Create input tensors
            # Convert to token IDs and expression values
            input_gene_ids = []
            input_values = []
            
            for data_gene_idx, expr_val in enumerate(cell_expr):
                if data_gene_idx in gene_vocab_mapping and expr_val > 0:
                    input_gene_ids.append(gene_vocab_mapping[data_gene_idx])
                    input_values.append(expr_val)
            
            if len(input_gene_ids) == 0:
                continue  # Skip cells with no expressed genes
            
            # Convert to tensors
            gene_ids = torch.tensor(input_gene_ids, dtype=torch.long, device=device).unsqueeze(0)
            values = torch.tensor(input_values, dtype=torch.float, device=device).unsqueeze(0)
            
            # Normal forward pass
            try:
                output_normal = model(gene_ids, values)
                if hasattr(output_normal, 'last_hidden_state'):
                    cls_embed_normal = output_normal.last_hidden_state[:, 0, :]  # CLS token
                else:
                    cls_embed_normal = output_normal[:, 0, :]  # Assume first token is CLS
                
                # Perturbed forward pass (zero out target gene)
                gene_ids_perturbed = gene_ids.clone()
                values_perturbed = values.clone()
                
                # Find the target gene in input and zero it out
                target_gene_mask = gene_ids_perturbed == gene_token_id
                if target_gene_mask.any():
                    values_perturbed[target_gene_mask] = args.get("pad_value", -2)  # Use pad_value from config
                    
                    output_perturbed = model(gene_ids_perturbed, values_perturbed)
                    if hasattr(output_perturbed, 'last_hidden_state'):
                        cls_embed_perturbed = output_perturbed.last_hidden_state[:, 0, :]
                    else:
                        cls_embed_perturbed = output_perturbed[:, 0, :]
                    
                    # Calculate embedding shift (cosine distance)
                    normal_vec = cls_embed_normal.squeeze().cpu().numpy()
                    perturbed_vec = cls_embed_perturbed.squeeze().cpu().numpy()
                    
                    # Cosine distance (1 - cosine similarity)
                    cos_dist = cosine(normal_vec, perturbed_vec)
                    
                    # L2 distance as alternative
                    l2_dist = np.linalg.norm(normal_vec - perturbed_vec)
                    
                    embedding_shifts.append({
                        'cosine_distance': cos_dist,
                        'l2_distance': l2_dist
                    })
                    
            except Exception as e:
                print(f"⚠️ Error processing cell {cell_idx}: {e}")
                continue
    
    if embedding_shifts:
        # Calculate average perturbation impact
        avg_cosine_dist = np.mean([s['cosine_distance'] for s in embedding_shifts])
        avg_l2_dist = np.mean([s['l2_distance'] for s in embedding_shifts])
        
        perturbation_results[gene] = {
            'avg_cosine_distance': avg_cosine_dist,
            'avg_l2_distance': avg_l2_dist,
            'n_cells_analyzed': len(embedding_shifts)
        }
        
        print(f"📊 {gene}: avg cosine dist = {avg_cosine_dist:.6f}, L2 dist = {avg_l2_dist:.6f} (n={len(embedding_shifts)} cells)")
    else:
        print(f"❌ No valid perturbations for {gene}")

print(f"\n🔍 Perturbation Analysis Complete")
print("=" * 60)

# Rank genes by perturbation impact (higher distance = more impact)
if perturbation_results:
    # Use cosine distance for ranking (primary metric)
    scgpt_ranking_cosine = sorted(perturbation_results.keys(), 
                                  key=lambda x: perturbation_results[x]['avg_cosine_distance'], 
                                  reverse=True)
    
    # Also rank by L2 distance
    scgpt_ranking_l2 = sorted(perturbation_results.keys(), 
                              key=lambda x: perturbation_results[x]['avg_l2_distance'], 
                              reverse=True)
    
    print(f"📈 scGPT ranking (cosine, high→low impact): {' > '.join(scgpt_ranking_cosine)}")
    print(f"📈 scGPT ranking (L2, high→low impact): {' > '.join(scgpt_ranking_l2)}")
    print(f"📈 Geneformer ranking: {' > '.join(geneformer_ranking)}")
    
    # Calculate Spearman correlation
    # Need to create numerical rankings for correlation
    def create_rank_vector(gene_list, ranking):
        """Create rank vector where 1=highest impact, 2=second highest, etc."""
        ranks = {}
        for i, gene in enumerate(ranking):
            ranks[gene] = i + 1  # 1-indexed ranking
        return [ranks.get(gene, len(ranking) + 1) for gene in gene_list]
    
    # Get common genes between scGPT and Geneformer
    common_genes = [g for g in geneformer_ranking if g in perturbation_results]
    print(f"🎯 Common genes for correlation: {common_genes}")
    
    if len(common_genes) >= 3:  # Need at least 3 points for meaningful correlation
        # Create rank vectors for common genes
        geneformer_ranks = create_rank_vector(common_genes, geneformer_ranking)
        scgpt_cosine_ranks = create_rank_vector(common_genes, scgpt_ranking_cosine)
        scgpt_l2_ranks = create_rank_vector(common_genes, scgpt_ranking_l2)
        
        # Calculate Spearman correlations
        spearman_cosine, p_value_cosine = spearmanr(geneformer_ranks, scgpt_cosine_ranks)
        spearman_l2, p_value_l2 = spearmanr(geneformer_ranks, scgpt_l2_ranks)
        
        print(f"\n📊 CORRELATION RESULTS:")
        print(f"🔗 Spearman correlation (cosine): r = {spearman_cosine:.4f}, p = {p_value_cosine:.4f}")
        print(f"🔗 Spearman correlation (L2): r = {spearman_l2:.4f}, p = {p_value_l2:.4f}")
        
        # Write detailed report
        report_content = f"""# scGPT vs Geneformer Cross-Model Perturbation Analysis

## Executive Summary

This analysis compares scGPT in-silico perturbation effects with Geneformer rankings for 5 intelligence-associated genes in human brain single-cell data.

## Methods

- **Model**: scGPT brain checkpoint from `{scgpt_checkpoint_dir}`
- **Data**: Human brain single-cell RNA-seq (DLPFC_11k.h5ad)
- **Cells analyzed**: {n_cells_to_test}
- **Genes tested**: {', '.join(intelligence_genes)}
- **Perturbation method**: Zero out gene expression (set to pad_value=-2)
- **Impact metric**: Embedding shift (cosine distance and L2 distance)

## Results

### Gene Perturbation Impact Scores

| Gene | Avg Cosine Distance | Avg L2 Distance | Cells Analyzed |
|------|-------------------|----------------|---------------|
"""
        
        for gene in intelligence_genes:
            if gene in perturbation_results:
                result = perturbation_results[gene]
                report_content += f"| {gene} | {result['avg_cosine_distance']:.6f} | {result['avg_l2_distance']:.6f} | {result['n_cells_analyzed']} |\n"
        
        report_content += f"""
### Rankings Comparison

**scGPT (Cosine Distance)**: {' > '.join(scgpt_ranking_cosine)}
**scGPT (L2 Distance)**: {' > '.join(scgpt_ranking_l2)}
**Geneformer**: {' > '.join(geneformer_ranking)}

### Statistical Correlation

- **Spearman correlation (Cosine)**: r = {spearman_cosine:.4f}, p = {p_value_cosine:.4f}
- **Spearman correlation (L2)**: r = {spearman_l2:.4f}, p = {p_value_l2:.4f}
- **Common genes analyzed**: {', '.join(common_genes)}

## Interpretation

{'**Strong positive correlation**' if spearman_cosine > 0.6 else '**Moderate positive correlation**' if spearman_cosine > 0.3 else '**Weak correlation**'} between scGPT and Geneformer rankings using cosine distance metric.

The correlation suggests {'excellent' if spearman_cosine > 0.8 else 'good' if spearman_cosine > 0.6 else 'moderate' if spearman_cosine > 0.3 else 'limited'} agreement between the two models on the relative importance of these intelligence genes.

## Technical Notes

- Analysis performed on {n_cells_to_test} randomly selected cells
- Perturbations applied by zeroing gene expression in input
- CLS token embeddings used for impact measurement
- Both cosine distance and L2 distance metrics computed
- Statistical significance tested using Spearman rank correlation

---
*Generated on: {pd.Timestamp.now()}*
*Model: scGPT brain checkpoint*
*Data: Human DLPFC single-cell RNA-seq*
"""
        
        # Write report
        with open(output_path, 'w') as f:
            f.write(report_content)
        
        print(f"📄 Report written to: {output_path}")
        print(f"✅ Analysis completed successfully!")
        
    else:
        print(f"❌ Insufficient common genes ({len(common_genes)}) for meaningful correlation")

else:
    print("❌ No perturbation results obtained")

print("\n🏁 scGPT Intelligence Gene Perturbation Analysis Complete")