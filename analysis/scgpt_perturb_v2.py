"""scGPT in-silico perturbation of 5 intelligence genes — cross-model validation."""
exec(open("/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py").read())

import json, sys, torch, numpy as np, scanpy as sc
from pathlib import Path
from scgpt.tokenizer import GeneVocab

CKPT = Path("/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT_checkpoints/brain")
DATA = Path("/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")

# Intelligence genes: symbol -> ENSEMBL
GENES = {
    "CADM2": "ENSG00000175161",
    "GRIN2A": "ENSG00000183454",
    "CAMK2A": "ENSG00000070808",
    "MEF2C": "ENSG00000081189",
    "APP": "ENSG00000142192",
}

# Geneformer rankings (from existing results)
GENEFORMER_RANK = {"CADM2": 1, "GRIN2A": 2, "CAMK2A": 3, "MEF2C": 4, "APP": 5}

print("Loading vocab...")
vocab = GeneVocab.from_file(str(CKPT / "vocab.json"))
stoi = vocab.get_stoi()
print(f"Vocab: {len(vocab)} tokens")

# Check which genes are in vocab (by symbol)
for g in GENES:
    if g in stoi:
        print(f"  {g}: vocab_id={stoi[g]}")
    else:
        print(f"  {g}: NOT in vocab")

print("\nLoading data...")
adata = sc.read_h5ad(str(DATA))
print(f"Data: {adata.shape}")

# Map ENSEMBL -> symbol using feature_name if available
if "feature_name" in adata.var.columns:
    ens_to_sym = dict(zip(adata.var.index, adata.var["feature_name"]))
    print(f"Using feature_name column for mapping")
else:
    ens_to_sym = {}
    print("No feature_name column — trying var_names directly")

# Find genes in both dataset AND vocab
gene_map = {}  # symbol -> (adata_idx, vocab_id)
for sym, ens in GENES.items():
    # Check vocab
    if sym not in stoi:
        print(f"SKIP {sym}: not in scGPT vocab")
        continue
    
    # Check dataset — try ENSEMBL first, then feature_name
    if ens in adata.var_names:
        adata_idx = list(adata.var_names).index(ens)
        gene_map[sym] = (adata_idx, stoi[sym])
        print(f"OK {sym}: adata_idx={adata_idx}, vocab_id={stoi[sym]}")
    elif ens_to_sym:
        # Check if any var maps to this symbol
        matches = [k for k, v in ens_to_sym.items() if v == sym]
        if matches:
            adata_idx = list(adata.var_names).index(matches[0])
            gene_map[sym] = (adata_idx, stoi[sym])
            print(f"OK {sym}: adata_idx={adata_idx} (via feature_name), vocab_id={stoi[sym]}")
        else:
            print(f"SKIP {sym}: not found in dataset")
    else:
        print(f"SKIP {sym}: ENSEMBL {ens} not in var_names")

if not gene_map:
    print("ERROR: No genes found in both dataset and vocab!")
    sys.exit(1)

print(f"\n{len(gene_map)} genes available for perturbation")

# Load model
print("\nLoading scGPT model...")
with open(CKPT / "args.json") as f:
    args = json.load(f)

from scgpt.model.model import TransformerModel

pad_id = stoi.get("<pad>", len(vocab) - 3)
model = TransformerModel(
    ntoken=len(vocab),
    d_model=args["embsize"],
    nhead=args["nheads"],
    d_hid=args["d_hid"],
    nlayers=args["nlayers"],
    vocab=vocab,
    dropout=args.get("dropout", 0.2),
    pad_token="<pad>",
    pad_value=-2,
    do_mvc=False,
    do_dab=False,
    use_batch_labels=False,
    n_input_bins=args.get("n_bins", 51),
    input_emb_style=args.get("input_emb_style", "continuous"),
    cell_emb_style="cls",
    mvc_decoder_style="inner product",
    ecs_threshold=args.get("ecs_thres", 0.0),
    explicit_zero_prob=not args.get("no_explicit_zero_prob", False),
    use_fast_transformer=False,
)

state = torch.load(str(CKPT / "best_model.pt"), map_location="cpu")
# Filter mismatched keys
model_state = model.state_dict()
filtered = {}
for k, v in state.items():
    if k in model_state and v.shape == model_state[k].shape:
        filtered[k] = v
print(f"Loaded {len(filtered)}/{len(state)} parameters (filtered shape mismatches)")
model.load_state_dict(filtered, strict=False)
model.eval()

device = "cuda" if torch.cuda.is_available() else "cpu"
model = model.to(device)
print(f"Device: {device}")

# Sample cells
N_CELLS = 200
np.random.seed(42)
cell_idx = np.random.choice(adata.n_obs, min(N_CELLS, adata.n_obs), replace=False)
print(f"\nSampling {len(cell_idx)} cells for perturbation")

# Prepare expression matrix
from scipy.sparse import issparse
X = adata.X[cell_idx]
if issparse(X):
    X = X.toarray()
X = X.astype(np.float32)

# Get top expressed genes per cell (scGPT uses ranked genes)
# We'll use top 512 genes per cell
MAX_GENES = 512

def get_cell_input(cell_expr, gene_names_list):
    """Convert a cell's expression to scGPT input format."""
    # Get top expressed genes
    nonzero = np.where(cell_expr > 0)[0]
    if len(nonzero) == 0:
        return None, None, None
    
    # Sort by expression, take top MAX_GENES
    expr_vals = cell_expr[nonzero]
    top_idx = nonzero[np.argsort(-expr_vals)[:MAX_GENES]]
    
    # Map to vocab IDs
    gene_ids = []
    values = []
    for idx in top_idx:
        gene_name = gene_names_list[idx]
        # Try symbol first (via feature_name), then ENSEMBL
        sym = ens_to_sym.get(gene_name, gene_name) if ens_to_sym else gene_name
        if sym in stoi:
            gene_ids.append(stoi[sym])
            values.append(float(cell_expr[idx]))
    
    if not gene_ids:
        return None, None, None
    
    return gene_ids, values, len(gene_ids)

gene_names_list = list(adata.var_names)

# Run perturbation
print("\nRunning perturbations...")
results = {}

for gene_sym, (adata_gene_idx, vocab_id) in gene_map.items():
    print(f"\n--- Perturbing {gene_sym} ---")
    shifts = []
    
    for ci, cell_i in enumerate(cell_idx):
        cell_expr = X[ci]
        gene_ids, values, n = get_cell_input(cell_expr, gene_names_list)
        if gene_ids is None:
            continue
        
        # Check if this gene is in this cell's top genes
        if vocab_id not in gene_ids:
            continue
        
        # Normal forward pass
        gene_t = torch.tensor([gene_ids], dtype=torch.long, device=device)
        val_t = torch.tensor([values], dtype=torch.float, device=device)
        
        with torch.no_grad():
            # scGPT forward: gene_ids, values -> cell embedding
            src_key_padding_mask = torch.zeros(1, len(gene_ids), dtype=torch.bool, device=device)
            output_normal = model._encode(gene_t, val_t, src_key_padding_mask)
            emb_normal = output_normal[:, 0, :]  # CLS token
            
            # Perturbed: zero out the target gene's expression
            perturbed_values = values.copy()
            gene_pos = gene_ids.index(vocab_id)
            perturbed_values[gene_pos] = 0.0
            val_perturbed = torch.tensor([perturbed_values], dtype=torch.float, device=device)
            
            output_perturbed = model._encode(gene_t, val_perturbed, src_key_padding_mask)
            emb_perturbed = output_perturbed[:, 0, :]
            
            # Cosine distance
            cos_sim = torch.nn.functional.cosine_similarity(emb_normal, emb_perturbed).item()
            shift = 1.0 - cos_sim
            shifts.append(shift)
    
    mean_shift = np.mean(shifts) if shifts else 0.0
    std_shift = np.std(shifts) if shifts else 0.0
    results[gene_sym] = {
        "mean_shift": float(mean_shift),
        "std_shift": float(std_shift),
        "n_cells": len(shifts),
    }
    print(f"  {gene_sym}: mean_shift={mean_shift:.6f}, std={std_shift:.6f}, n={len(shifts)}")

# Rank by scGPT
scgpt_ranked = sorted(results.items(), key=lambda x: -x[1]["mean_shift"])
print("\n=== scGPT RANKINGS ===")
for i, (g, r) in enumerate(scgpt_ranked, 1):
    print(f"  {i}. {g}: shift={r['mean_shift']:.6f} (n={r['n_cells']})")

# Spearman correlation with Geneformer
from scipy.stats import spearmanr
scgpt_ranks = {g: i+1 for i, (g, _) in enumerate(scgpt_ranked)}
common = sorted(set(scgpt_ranks) & set(GENEFORMER_RANK))
if len(common) >= 3:
    x = [GENEFORMER_RANK[g] for g in common]
    y = [scgpt_ranks[g] for g in common]
    rho, pval = spearmanr(x, y)
    print(f"\nSpearman correlation (Geneformer vs scGPT): rho={rho:.4f}, p={pval:.4f}")
else:
    rho, pval = None, None
    print(f"\nNot enough common genes for correlation ({len(common)})")

# Save results
out = {
    "experiment": "scgpt_perturbation_crossmodel",
    "n_cells_sampled": len(cell_idx),
    "genes": results,
    "scgpt_ranking": [g for g, _ in scgpt_ranked],
    "geneformer_ranking": ["CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP"],
    "spearman_rho": rho,
    "spearman_pval": pval,
}
with open("/mnt/d/openclaw/intelligence-augmentation/analysis/scgpt_crossmodel_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nResults saved to scgpt_crossmodel_results.json")
