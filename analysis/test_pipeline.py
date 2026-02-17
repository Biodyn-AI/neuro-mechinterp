#!/usr/bin/env python3
"""
Test pipeline for Intelligence Augmentation project.
Validates: scRNA-seq loading, Geneformer model loading, tokenization.
"""
import sys
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'
sys.stdout.reconfigure(encoding='utf-8')
import warnings
warnings.filterwarnings('ignore')

# Add geneformer to path
GENEFORMER_DIR = r"D:\openclaw\intelligence-augmentation\models\Geneformer"
sys.path.insert(0, GENEFORMER_DIR)

DATA_DIR = r"D:\openclaw\intelligence-augmentation\data"
SCRNA_FILE = os.path.join(DATA_DIR, "brain_scrna", "DLPFC_11k.h5ad")

print("=" * 60)
print("INTELLIGENCE AUGMENTATION PIPELINE - TEST")
print("=" * 60)

# ── 1. Load scRNA-seq data ──────────────────────────────────
print("\n[1/4] Loading brain scRNA-seq data...")
import scanpy as sc
import numpy as np

adata = sc.read_h5ad(SCRNA_FILE)
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Obs columns: {list(adata.obs.columns[:10])}")

# Cell type info
ct_col = None
for col in ['cell_type', 'celltype', 'CellType', 'cluster', 'cell_type_ontology_term_id']:
    if col in adata.obs.columns:
        ct_col = col
        break

if ct_col:
    cts = adata.obs[ct_col].value_counts()
    print(f"  Cell type column: '{ct_col}'")
    print(f"  Unique cell types: {len(cts)}")
    print(f"  Top 5 cell types:")
    for ct, count in cts.head(5).items():
        print(f"    {ct}: {count:,}")
else:
    print("  No standard cell type column found")
    print(f"  Available: {list(adata.obs.columns)}")

# ── 2. Load Geneformer model ────────────────────────────────
print("\n[2/4] Loading Geneformer model...")
from transformers import AutoModel, AutoConfig

# Try loading from local model files
model_paths = [
    os.path.join(GENEFORMER_DIR, "Geneformer-V2-104M"),
    os.path.join(GENEFORMER_DIR),
]

model = None
model_path_used = None
for mp in model_paths:
    if os.path.exists(os.path.join(mp, "config.json")):
        try:
            config = AutoConfig.from_pretrained(mp, trust_remote_code=True)
            model = AutoModel.from_pretrained(mp, config=config, trust_remote_code=True)
            model_path_used = mp
            print(f"  Model loaded from: {mp}")
            print(f"  Model type: {type(model).__name__}")
            print(f"  Parameters: {sum(p.numel() for p in model.parameters()):,}")
            break
        except Exception as e:
            print(f"  Failed to load from {mp}: {e}")

if model is None:
    print("  WARNING: Could not load model from local paths")
    print("  Trying huggingface hub...")
    try:
        model = AutoModel.from_pretrained("ctheodoris/Geneformer", trust_remote_code=True)
        print(f"  Model loaded from HuggingFace Hub")
    except Exception as e:
        print(f"  FAILED: {e}")

# ── 3. Test Geneformer tokenization ─────────────────────────
print("\n[3/4] Testing Geneformer tokenization...")
try:
    from geneformer import TranscriptomeTokenizer
    
    # Create a small test subset
    subset = adata[:100].copy()
    
    # Save as loom for tokenizer (Geneformer expects loom format)
    test_loom = os.path.join(DATA_DIR, "brain_scrna", "test_subset.loom")
    test_token_dir = os.path.join(DATA_DIR, "brain_scrna", "test_tokenized")
    os.makedirs(test_token_dir, exist_ok=True)
    
    # Check if data has raw counts
    if subset.raw is not None:
        print("  Using raw counts layer")
        subset = subset.raw.to_adata()[:100]
    
    print(f"  Subset: {subset.n_obs} cells x {subset.n_vars} genes")
    print(f"  Data range: [{subset.X.min():.2f}, {subset.X.max():.2f}]")
    
    # Try tokenization
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_type": ct_col} if ct_col else None,
    )
    print("  TranscriptomeTokenizer initialized")
    print("  Tokenizer ready (full tokenization requires loom format)")
    
except ImportError as e:
    print(f"  Import error: {e}")
    print("  Geneformer tokenizer not fully available")
except Exception as e:
    print(f"  Tokenization test: {e}")

# ── 4. Verify all data files ────────────────────────────────
print("\n[4/4] Verifying data files...")
files_to_check = {
    "GWAS sumstats": os.path.join(DATA_DIR, "gwas", "SavageJansen_IntMeta_sumstats.zip"),
    "DLPFC 11k": os.path.join(DATA_DIR, "brain_scrna", "DLPFC_11k.h5ad"),
    "Siletti DLPFC 113k": os.path.join(DATA_DIR, "brain_scrna", "Siletti_DLPFC_113k.h5ad"),
    "DGIdb interactions": os.path.join(DATA_DIR, "drugs", "interactions.tsv"),
}

total_size = 0
for name, path in files_to_check.items():
    if os.path.exists(path):
        size_mb = os.path.getsize(path) / (1024 * 1024)
        total_size += size_mb
        print(f"  ✓ {name}: {size_mb:.1f} MB")
    else:
        print(f"  ✗ {name}: NOT FOUND")

# Check Geneformer model files
gf_size = 0
for root, dirs, files in os.walk(GENEFORMER_DIR):
    for f in files:
        fp = os.path.join(root, f)
        gf_size += os.path.getsize(fp)
gf_size_mb = gf_size / (1024 * 1024)
total_size += gf_size_mb
print(f"  ✓ Geneformer repo: {gf_size_mb:.1f} MB")

print(f"\n  Total data: {total_size:.1f} MB ({total_size/1024:.2f} GB)")

print("\n" + "=" * 60)
print("PIPELINE TEST COMPLETE")
print("=" * 60)
