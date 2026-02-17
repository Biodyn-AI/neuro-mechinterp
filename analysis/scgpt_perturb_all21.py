"""scGPT perturbation for ALL 21 intelligence genes — extends v3 approach."""
exec(open("/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py").read())

import json, sys, torch, numpy as np, scanpy as sc
from pathlib import Path
from scipy.sparse import issparse
from scgpt.tokenizer import GeneVocab

CKPT = Path("/mnt/d/openclaw/mechinterp-bio/biodyn-work/single_cell_mechinterp/external/scGPT_checkpoints/brain")
DATA = Path("/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad")

ALL_GENES = [
    "CADM2", "GRIN2A", "CAMK2A", "MEF2C", "APP", "SCN1A", "NRXN1", "GRIN2B",
    "HOMER1", "NEGR1", "TCF4", "NLGN1", "FOXO3", "MAPT", "SHANK3", "FMR1",
    "BDNF", "CREB1", "SYN1", "SLC6A4", "COMT"
]

print("Loading vocab...")
vocab = GeneVocab.from_file(str(CKPT / "vocab.json"))
stoi = vocab.get_stoi()

print("Loading data...")
adata = sc.read_h5ad(str(DATA))
ens_to_sym = dict(zip(adata.var.index, adata.var["feature_name"])) if "feature_name" in adata.var.columns else {}

# Map genes
gene_map = {}
for sym in ALL_GENES:
    if sym not in stoi:
        print(f"  {sym}: NOT in scGPT vocab, skipping")
        continue
    matches = [k for k, v in ens_to_sym.items() if v == sym]
    if matches:
        gene_map[sym] = (list(adata.var_names).index(matches[0]), stoi[sym])
        print(f"  {sym}: found, data_idx={gene_map[sym][0]}, vocab_idx={gene_map[sym][1]}")
    else:
        print(f"  {sym}: not in adata, skipping")

print(f"\nMapped {len(gene_map)}/{len(ALL_GENES)} genes")

# Load model
from scgpt.model import TransformerModel
import scgpt

model_file = CKPT / "best_model.pt"
model_config = json.load(open(CKPT / "args.json"))

ntoken = len(vocab)
d_model = model_config.get("embsize", 512)
nhead = model_config.get("nheads", 8)
d_hid = model_config.get("d_hid", 512)
nlayers = model_config.get("nlayers", 12)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")

model = TransformerModel(
    ntoken=ntoken,
    d_model=d_model,
    nhead=nhead,
    d_hid=d_hid,
    nlayers=nlayers,
    vocab=vocab,
    dropout=0.0,
    pad_token="<pad>",
    pad_value=-2,
    do_mvc=False, do_dab=False, use_batch_labels=False,
    n_input_bins=model_config.get("n_bins", 51),
    input_emb_style=model_config.get("input_emb_style", "continuous"),
    cell_emb_style="cls", mvc_decoder_style="inner product",
    ecs_threshold=0.0, explicit_zero_prob=False, use_fast_transformer=False,
)

state = torch.load(model_file, map_location="cpu")
model_sd = model.state_dict()
filtered_state = {k: v for k, v in state.items() if k in model_sd and v.shape == model_sd[k].shape}
model.load_state_dict(filtered_state, strict=False)
model.to(device)
model.eval()
print("Model loaded")

# Tokenize cells
N_CELLS = 500  # more cells for better stats
X = adata.X[:N_CELLS]
if issparse(X):
    X = X.toarray()
X = X.astype(np.float32)

gene_names = list(adata.var_names)
sym_names = [ens_to_sym.get(g, g) for g in gene_names]

# Build per-cell token sequences
MAX_SEQ = 512

def tokenize_cell(expr_vec):
    nonzero = np.nonzero(expr_vec)[0]
    if len(nonzero) == 0:
        return [], []
    vals = expr_vec[nonzero]
    order = np.argsort(-vals)
    tokens = []
    values = []
    for idx in order[:MAX_SEQ]:
        sym = sym_names[idx]
        if sym in stoi:
            tokens.append(stoi[sym])
            values.append(float(vals[idx]))
    return tokens, values

print("Tokenizing cells...")
all_tokens = []
all_values = []
for i in range(N_CELLS):
    t, v = tokenize_cell(X[i])
    if len(t) >= 10:
        all_tokens.append(t)
        all_values.append(v)
print(f"Tokenized {len(all_tokens)} cells")

def get_embedding(tokens, values):
    t = torch.tensor([tokens], dtype=torch.long, device=device)
    v = torch.tensor([values], dtype=torch.float, device=device)
    mask = torch.zeros(1, len(tokens), dtype=torch.bool, device=device)
    with torch.no_grad():
        out = model._encode(t, v, mask)
    return out[:, 0, :].squeeze().cpu().numpy()  # CLS token

# Perturbation: remove gene token and measure L2 shift
results = {}
for gene_sym, (data_idx, vocab_idx) in gene_map.items():
    shifts = []
    n_expressing = 0
    for i in range(len(all_tokens)):
        if vocab_idx not in all_tokens[i]:
            continue
        n_expressing += 1
        
        # Original embedding
        orig_emb = get_embedding(all_tokens[i], all_values[i])
        
        # Perturbed: remove gene token
        idx = all_tokens[i].index(vocab_idx)
        perturbed_tokens = all_tokens[i][:idx] + all_tokens[i][idx+1:]
        perturbed_values = all_values[i][:idx] + all_values[i][idx+1:]
        
        if len(perturbed_tokens) < 5:
            continue
        
        pert_emb = get_embedding(perturbed_tokens, perturbed_values)
        shift = np.linalg.norm(orig_emb - pert_emb)
        shifts.append(shift)
    
    if shifts:
        results[gene_sym] = {
            "mean_shift": float(np.mean(shifts)),
            "std_shift": float(np.std(shifts)),
            "n_cells": len(shifts),
            "n_expressing": n_expressing,
        }
        print(f"  {gene_sym}: shift={np.mean(shifts):.6f} ± {np.std(shifts):.6f}, n={len(shifts)}")
    else:
        print(f"  {gene_sym}: 0 expressing cells in sample")

# Save
output = {
    "experiment": "scgpt_all21_perturbation",
    "n_cells_sampled": N_CELLS,
    "n_cells_tokenized": len(all_tokens),
    "genes": results,
}
out_path = "/mnt/d/openclaw/intelligence-augmentation/analysis/scgpt_all21_results.json"
with open(out_path, "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")

# Cross-model comparison with Geneformer
geneformer_shifts = {
    'CADM2': 0.0196, 'GRIN2A': 0.0190, 'CAMK2A': 0.0189, 'MEF2C': 0.0184,
    'APP': 0.0183, 'SCN1A': 0.0179, 'NRXN1': 0.0178, 'GRIN2B': 0.0176,
    'HOMER1': 0.0175, 'NEGR1': 0.0166, 'TCF4': 0.0164, 'NLGN1': 0.0161,
    'FOXO3': 0.0159, 'MAPT': 0.0155, 'SHANK3': 0.0148, 'FMR1': 0.0144,
    'BDNF': 0.0140, 'CREB1': 0.0136, 'SYN1': 0.0130, 'SLC6A4': 0.0119,
    'COMT': 0.0105,
}

common = sorted(set(results.keys()) & set(geneformer_shifts.keys()))
if len(common) >= 5:
    from scipy.stats import spearmanr
    gf = [geneformer_shifts[g] for g in common]
    sc_vals = [results[g]["mean_shift"] for g in common]
    rho, p = spearmanr(gf, sc_vals)
    print(f"\nCross-model (n={len(common)}): Spearman rho={rho:.3f}, p={p:.4f}")
