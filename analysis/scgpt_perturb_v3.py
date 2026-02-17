"""scGPT perturbation v3 — REMOVE gene from input sequence entirely."""
exec(open("/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py").read())

import json, sys, torch, numpy as np, scanpy as sc
from pathlib import Path
from scipy.sparse import issparse
from scgpt.tokenizer import GeneVocab

CKPT = Path("/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT_checkpoints/brain")
DATA = Path("/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")

GENES = {"CADM2": "ENSG00000175161", "GRIN2A": "ENSG00000183454", "CAMK2A": "ENSG00000070808", "MEF2C": "ENSG00000081189", "APP": "ENSG00000142192"}
GENEFORMER_RANK = {"CADM2": 1, "GRIN2A": 2, "CAMK2A": 3, "MEF2C": 4, "APP": 5}

print("Loading vocab...")
vocab = GeneVocab.from_file(str(CKPT / "vocab.json"))
stoi = vocab.get_stoi()

print("Loading data...")
adata = sc.read_h5ad(str(DATA))
ens_to_sym = dict(zip(adata.var.index, adata.var["feature_name"])) if "feature_name" in adata.var.columns else {}

# Map genes
gene_map = {}
for sym, ens in GENES.items():
    if sym not in stoi:
        continue
    matches = [k for k, v in ens_to_sym.items() if v == sym]
    if matches:
        gene_map[sym] = (list(adata.var_names).index(matches[0]), stoi[sym])
print(f"{len(gene_map)} genes mapped")

# Load model
with open(CKPT / "args.json") as f:
    args = json.load(f)

from scgpt.model.model import TransformerModel
model = TransformerModel(
    ntoken=len(vocab), d_model=args["embsize"], nhead=args["nheads"],
    d_hid=args["d_hid"], nlayers=args["nlayers"], vocab=vocab,
    dropout=0.0, pad_token="<pad>", pad_value=-2,
    do_mvc=False, do_dab=False, use_batch_labels=False,
    n_input_bins=args.get("n_bins", 51),
    input_emb_style=args.get("input_emb_style", "continuous"),
    cell_emb_style="cls", mvc_decoder_style="inner product",
    ecs_threshold=0.0, explicit_zero_prob=False, use_fast_transformer=False,
)
state = torch.load(str(CKPT / "best_model.pt"), map_location="cpu")
model_sd = model.state_dict()
filtered = {k: v for k, v in state.items() if k in model_sd and v.shape == model_sd[k].shape}
model.load_state_dict(filtered, strict=False)
model.eval()

device = "cuda" if torch.cuda.is_available() else "cpu"
model = model.to(device)
print(f"Device: {device}")

# Sample cells
N_CELLS = 200
MAX_GENES = 512
np.random.seed(42)
cell_idx = np.random.choice(adata.n_obs, N_CELLS, replace=False)
X = adata.X[cell_idx]
if issparse(X):
    X = X.toarray()
X = X.astype(np.float32)
gene_names = list(adata.var_names)

def get_cell_tokens(cell_expr):
    """Get top MAX_GENES expressed genes as (gene_ids, values)."""
    nonzero = np.where(cell_expr > 0)[0]
    if len(nonzero) == 0:
        return [], []
    expr_vals = cell_expr[nonzero]
    top_idx = nonzero[np.argsort(-expr_vals)[:MAX_GENES]]
    gene_ids, values = [], []
    for idx in top_idx:
        sym = ens_to_sym.get(gene_names[idx], gene_names[idx])
        if sym in stoi:
            gene_ids.append(stoi[sym])
            values.append(float(cell_expr[idx]))
    return gene_ids, values

def get_embedding(gene_ids, values):
    """Forward pass, return CLS embedding."""
    g = torch.tensor([gene_ids], dtype=torch.long, device=device)
    v = torch.tensor([values], dtype=torch.float, device=device)
    mask = torch.zeros(1, len(gene_ids), dtype=torch.bool, device=device)
    with torch.no_grad():
        out = model._encode(g, v, mask)
    return out[:, 0, :].squeeze()  # CLS

print(f"\nRunning perturbations on {N_CELLS} cells...")
results = {}

for gene_sym, (adata_gene_idx, vocab_id) in gene_map.items():
    print(f"\n--- {gene_sym} (vocab_id={vocab_id}) ---")
    shifts = []
    
    for ci in range(len(cell_idx)):
        gene_ids, values = get_cell_tokens(X[ci])
        if not gene_ids or vocab_id not in gene_ids:
            continue
        
        # Normal embedding
        emb_normal = get_embedding(gene_ids, values)
        
        # Perturbed: REMOVE the gene entirely from the sequence
        pos = gene_ids.index(vocab_id)
        perturbed_ids = gene_ids[:pos] + gene_ids[pos+1:]
        perturbed_vals = values[:pos] + values[pos+1:]
        
        if not perturbed_ids:
            continue
            
        emb_perturbed = get_embedding(perturbed_ids, perturbed_vals)
        
        # L2 distance (more sensitive than cosine for small changes)
        shift = torch.norm(emb_normal - emb_perturbed).item()
        shifts.append(shift)
    
    mean_shift = np.mean(shifts) if shifts else 0.0
    std_shift = np.std(shifts) if shifts else 0.0
    results[gene_sym] = {"mean_shift": float(mean_shift), "std_shift": float(std_shift), "n_cells": len(shifts)}
    print(f"  shift={mean_shift:.6f} +/- {std_shift:.6f} (n={len(shifts)})")

# Rankings
scgpt_ranked = sorted(results.items(), key=lambda x: -x[1]["mean_shift"])
print("\n=== scGPT RANKINGS (by L2 embedding shift) ===")
for i, (g, r) in enumerate(scgpt_ranked, 1):
    print(f"  {i}. {g}: shift={r['mean_shift']:.6f} (n={r['n_cells']})")

# Geneformer comparison
from scipy.stats import spearmanr
scgpt_ranks = {g: i+1 for i, (g, _) in enumerate(scgpt_ranked)}
common = sorted(set(scgpt_ranks) & set(GENEFORMER_RANK))
x = [GENEFORMER_RANK[g] for g in common]
y = [scgpt_ranks[g] for g in common]
rho, pval = spearmanr(x, y)
print(f"\nSpearman (Geneformer vs scGPT): rho={rho:.4f}, p={pval:.4f}")

# Save
out = {"experiment": "scgpt_perturbation_v3_removal", "n_cells": N_CELLS,
       "genes": results, "scgpt_ranking": [g for g, _ in scgpt_ranked],
       "geneformer_ranking": list(GENEFORMER_RANK.keys()),
       "spearman_rho": float(rho), "spearman_pval": float(pval)}
with open("/mnt/d/openclaw/intelligence-augmentation/analysis/scgpt_crossmodel_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("Saved to scgpt_crossmodel_results.json")
