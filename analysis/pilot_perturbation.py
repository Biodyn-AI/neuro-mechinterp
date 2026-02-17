"""
Pilot Perturbation Screening — Intelligence Enhancement Targets
Lean version: 200 cells, 50 genes, CPU, Geneformer V1-10M
"""
import sys
import os
import gc
import time
import pickle
import numpy as np
import pandas as pd

# Set memory-conservative options
os.environ["TOKENIZERS_PARALLELISM"] = "false"

PROJ = r"D:\openclaw\intelligence-augmentation"
MODEL_DIR = os.path.join(PROJ, "models", "Geneformer", "Geneformer-V1-10M")
DATA_FILE = os.path.join(PROJ, "data", "brain_scrna", "DLPFC_11k.h5ad")
GWAS_FILE = os.path.join(PROJ, "data", "gwas", "sumstats", "SavageJansen_2018_intelligence_metaanalysis.txt")
DRUG_FILE = os.path.join(PROJ, "data", "drugs", "interactions.tsv")
RESULTS_DIR = os.path.join(PROJ, "analysis", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Gene dictionaries path (V1 uses 30M dicts)
GENE_DICT_DIR = os.path.join(PROJ, "models", "Geneformer", "geneformer", "gene_dictionaries_30m")

print("=" * 60)
print("PILOT PERTURBATION SCREENING")
print("=" * 60)

# ---- STEP 1: Load and filter brain data ----
print("\n[1/6] Loading brain scRNA-seq data...")
t0 = time.time()
import scanpy as sc
adata = sc.read_h5ad(DATA_FILE)
print(f"  Full dataset: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Filter to glutamatergic neurons
glut = adata[adata.obs["cell_type"] == "glutamatergic neuron"].copy()
print(f"  Glutamatergic neurons: {glut.shape[0]} cells")

# Subsample to 200 cells for memory safety
if glut.shape[0] > 200:
    sc.pp.subsample(glut, n_obs=200, random_state=42)
    print(f"  Subsampled to: {glut.shape[0]} cells")

del adata
gc.collect()
print(f"  Time: {time.time()-t0:.1f}s")

# ---- STEP 2: Map gene names to Ensembl IDs ----
print("\n[2/6] Mapping genes to Ensembl IDs...")
t0 = time.time()

# Load Geneformer gene dictionaries
with open(os.path.join(GENE_DICT_DIR, "gene_name_id_dict_gc30M.pkl"), "rb") as f:
    gene_name_to_id = pickle.load(f)

with open(os.path.join(GENE_DICT_DIR, "token_dictionary_gc30M.pkl"), "rb") as f:
    token_dict = pickle.load(f)

# Check if gene names are already Ensembl or symbols
sample_genes = list(glut.var_names[:5])
print(f"  Sample gene names: {sample_genes}")
is_ensembl = sample_genes[0].startswith("ENSG")

if not is_ensembl:
    # Map gene symbols to Ensembl IDs
    mapped = {}
    for g in glut.var_names:
        g_upper = g.upper()
        if g_upper in gene_name_to_id:
            mapped[g] = gene_name_to_id[g_upper]
        elif g in gene_name_to_id:
            mapped[g] = gene_name_to_id[g]
    
    # Filter to mapped genes only
    keep_genes = [g for g in glut.var_names if g in mapped]
    glut = glut[:, keep_genes].copy()
    glut.var["ensembl_id"] = [mapped[g] for g in glut.var_names]
    print(f"  Mapped {len(keep_genes)} / {len(mapped)} genes to Ensembl IDs")
else:
    print(f"  Genes already in Ensembl format")

print(f"  Time: {time.time()-t0:.1f}s")

# ---- STEP 3: Define intelligence gene signature from GWAS ----
print("\n[3/6] Defining intelligence gene signature from GWAS...")
t0 = time.time()

# Load GWAS - only need significant hits
gwas = pd.read_csv(GWAS_FILE, sep="\t", usecols=["SNP", "CHR", "POS", "P", "Zscore", "stdBeta"])
print(f"  Total GWAS SNPs: {len(gwas)}")

# Genome-wide significant
sig = gwas[gwas["P"] < 5e-8].copy()
print(f"  Genome-wide significant (P < 5e-8): {len(sig)} SNPs")

# We need to map SNPs to genes. Since we don't have a full gene position file,
# use the known intelligence genes from literature instead (Savage 2018 top genes)
# These are the most replicated intelligence-associated genes
INTELLIGENCE_GENES = [
    "FOXO3", "MEF2C", "TCF4", "NEGR1", "KIBRA", "COMT", "BDNF", "NRXN1",
    "CADM2", "LRRC4C", "ADAM12", "CSE1L", "MAPT", "WNT4", "FGF5",
    "ATXN1", "DCC", "GRM8", "DLGAP1", "NCALD", "SHANK3", "GRIN2A",
    "SLC6A3", "DRD2", "DRD1", "CREB1", "PRKACA", "CAMK2A", "SYN1",
    "SH2B1", "RBFOX1", "MEF2A", "NRXN3", "CNTNAP2", "GRIK3",
    "ATP2B2", "KCNJ6", "CACNA1C", "SCN1A", "HCN1", "GABRA5",
    "SLC17A7", "GRM5", "HOMER1", "NLGN1", "SYNGAP1", "FMR1",
    "DYRK1A", "DISC1", "NRG1"
]

# Find which are in our data
if not is_ensembl:
    intel_in_data = [g for g in INTELLIGENCE_GENES if g in glut.var_names]
else:
    # Genes are Ensembl IDs in var_names — map intelligence gene names to Ensembl
    intel_ensembl = {g: gene_name_to_id.get(g) for g in INTELLIGENCE_GENES if g in gene_name_to_id}
    var_set = set(glut.var_names)
    intel_in_data = [g for g, eid in intel_ensembl.items() if eid in var_set]

print(f"  Intelligence genes in dataset: {len(intel_in_data)} / {len(INTELLIGENCE_GENES)}")
print(f"  Genes: {intel_in_data[:10]}...")

del gwas, sig
gc.collect()
print(f"  Time: {time.time()-t0:.1f}s")

# ---- STEP 4: Extract embeddings with Geneformer ----
print("\n[4/6] Computing Geneformer embeddings...")
t0 = time.time()

import torch
from transformers import BertForMaskedLM, BertConfig

# Load model
print(f"  Loading Geneformer V1-10M from {MODEL_DIR}...")
model = BertForMaskedLM.from_pretrained(MODEL_DIR)
model.eval()
print(f"  Model loaded. Parameters: {sum(p.numel() for p in model.parameters()):,}")

# Tokenize cells manually (rank-value encoding)
with open(os.path.join(GENE_DICT_DIR, "gene_median_dictionary_gc30M.pkl"), "rb") as f:
    gene_medians = pickle.load(f)

def tokenize_cell(expression_vector, gene_ensembl_ids, token_dict, gene_medians):
    """Rank-value encode a single cell for Geneformer."""
    # Get nonzero genes
    nonzero_mask = expression_vector > 0
    nonzero_genes = gene_ensembl_ids[nonzero_mask]
    nonzero_expr = expression_vector[nonzero_mask]
    
    # Normalize by median
    normalized = []
    for gene, expr in zip(nonzero_genes, nonzero_expr):
        median = gene_medians.get(gene, None)
        if median is not None and median > 0:
            normalized.append((gene, expr / median))
    
    # Sort by normalized expression (descending) — rank value encoding
    normalized.sort(key=lambda x: x[1], reverse=True)
    
    # Convert to token IDs (max 2048 for V1)
    tokens = []
    for gene, _ in normalized[:2048]:
        if gene in token_dict:
            tokens.append(token_dict[gene])
    
    return tokens

# Get expression matrix and gene IDs
if hasattr(glut.X, "toarray"):
    X = glut.X.toarray()
else:
    X = np.array(glut.X)

if "ensembl_id" in glut.var.columns:
    gene_ids = np.array(glut.var["ensembl_id"])
else:
    gene_ids = np.array(glut.var_names)

print(f"  Tokenizing {X.shape[0]} cells...")
tokenized_cells = []
for i in range(X.shape[0]):
    tokens = tokenize_cell(X[i], gene_ids, token_dict, gene_medians)
    if len(tokens) > 0:
        tokenized_cells.append(tokens)

print(f"  Successfully tokenized: {len(tokenized_cells)} cells")
print(f"  Avg tokens per cell: {np.mean([len(t) for t in tokenized_cells]):.0f}")

# Compute baseline embeddings
print(f"  Computing baseline embeddings...")

def get_embedding(model, token_ids, max_len=2048):
    """Get CLS embedding from Geneformer."""
    input_ids = torch.tensor([token_ids[:max_len]], dtype=torch.long)
    attention_mask = torch.ones_like(input_ids)
    with torch.no_grad():
        outputs = model.bert(input_ids=input_ids, attention_mask=attention_mask)
        # Use mean of last hidden state as embedding
        embedding = outputs.last_hidden_state.mean(dim=1).squeeze().numpy()
    return embedding

# Get baseline embeddings for all cells
baseline_embeddings = []
for i, tokens in enumerate(tokenized_cells):
    emb = get_embedding(model, tokens)
    baseline_embeddings.append(emb)
    if (i+1) % 50 == 0:
        print(f"    Embedded {i+1}/{len(tokenized_cells)} cells")

baseline_embeddings = np.array(baseline_embeddings)
mean_baseline = baseline_embeddings.mean(axis=0)
print(f"  Baseline embedding shape: {baseline_embeddings.shape}")
print(f"  Time: {time.time()-t0:.1f}s")

# ---- STEP 5: In-silico perturbation ----
print("\n[5/6] Running in-silico perturbation (gene deletion)...")
t0 = time.time()

# Perturb the intelligence genes that are in our data
genes_to_perturb = intel_in_data[:50]  # Max 50 for pilot
print(f"  Perturbing {len(genes_to_perturb)} intelligence-associated genes")

perturbation_results = []

for gi, gene_name in enumerate(genes_to_perturb):
    # Find the ensembl ID for this gene
    if is_ensembl:
        gene_ens = intel_ensembl.get(gene_name)
    else:
        gene_ens = gene_name_to_id.get(gene_name)
    
    if gene_ens is None or gene_ens not in token_dict:
        continue
    
    gene_token = token_dict[gene_ens]
    
    # For each cell, remove this gene's token and recompute embedding
    shifts = []
    for ci, tokens in enumerate(tokenized_cells):
        if gene_token in tokens:
            # Remove the gene token (simulate knockout)
            perturbed_tokens = [t for t in tokens if t != gene_token]
            if len(perturbed_tokens) > 0:
                perturbed_emb = get_embedding(model, perturbed_tokens)
                shift = np.linalg.norm(perturbed_emb - baseline_embeddings[ci])
                shifts.append(shift)
    
    if shifts:
        result = {
            "gene": gene_name,
            "ensembl_id": gene_ens,
            "mean_shift": np.mean(shifts),
            "std_shift": np.std(shifts),
            "max_shift": np.max(shifts),
            "n_cells_affected": len(shifts),
            "pct_cells_affected": len(shifts) / len(tokenized_cells) * 100
        }
        perturbation_results.append(result)
        
    if (gi+1) % 10 == 0:
        print(f"    Perturbed {gi+1}/{len(genes_to_perturb)} genes")

results_df = pd.DataFrame(perturbation_results)
results_df = results_df.sort_values("mean_shift", ascending=False)
print(f"\n  Results: {len(results_df)} genes successfully perturbed")
print(f"  Time: {time.time()-t0:.1f}s")

# ---- STEP 6: Cross-reference with drugs and generate report ----
print("\n[6/6] Cross-referencing with drug databases and generating report...")
t0 = time.time()

# Load drug interactions
drugs = pd.read_csv(DRUG_FILE, sep="\t")
drug_cols = list(drugs.columns)
print(f"  Drug interactions: {len(drugs)} entries, columns: {drug_cols[:5]}")

# Find drug interactions for our top genes
gene_col = [c for c in drug_cols if "gene" in c.lower()]
drug_col = [c for c in drug_cols if "drug" in c.lower() or "compound" in c.lower()]

if gene_col:
    gene_col = gene_col[0]
    drug_genes = set(drugs[gene_col].str.upper())
    results_df["has_drug"] = results_df["gene"].str.upper().isin(drug_genes)
    druggable = results_df[results_df["has_drug"]]
    print(f"  Druggable targets in results: {len(druggable)}")
else:
    results_df["has_drug"] = False
    print("  Warning: Could not identify gene column in drug database")

# Save results
results_df.to_csv(os.path.join(RESULTS_DIR, "perturbation_results.csv"), index=False)

# Generate report
report = f"""# Pilot Perturbation Screening Results
## Virtual Perturbation Screening for Cognitive Enhancement Targets

**Date:** {time.strftime('%Y-%m-%d %H:%M')}
**Model:** Geneformer V1-10M
**Data:** DLPFC glutamatergic neurons ({len(tokenized_cells)} cells)
**Genes perturbed:** {len(genes_to_perturb)} intelligence-associated genes
**Method:** Gene token deletion (knockout simulation)

---

## Top 20 Genes by Perturbation Impact

| Rank | Gene | Mean Embedding Shift | Cells Affected (%) | Druggable |
|------|------|---------------------|--------------------|-----------| 
"""

for i, row in results_df.head(20).iterrows():
    rank = results_df.index.tolist().index(i) + 1
    drug_flag = "Yes" if row.get("has_drug", False) else "No"
    report += f"| {rank} | {row['gene']} | {row['mean_shift']:.4f} | {row['pct_cells_affected']:.1f}% | {drug_flag} |\n"

report += f"""
## Summary Statistics
- Total genes perturbed: {len(results_df)}
- Mean embedding shift (all genes): {results_df['mean_shift'].mean():.4f}
- Max embedding shift: {results_df['mean_shift'].max():.4f} ({results_df.iloc[0]['gene']})
- Druggable targets in top 20: {len(results_df.head(20)[results_df.head(20).get('has_drug', False)])}

## Interpretation
Genes with larger embedding shifts have greater influence on the cell's transcriptomic state
when knocked out, suggesting they are key regulators of the gene programs active in DLPFC
pyramidal neurons. The top-ranked genes represent the highest-leverage targets for modulating
cognitive-relevant gene expression programs.

## Limitations
- Pilot scale only (200 cells, {len(genes_to_perturb)} genes)
- CPU-only, V1-10M model (smallest)
- Gene knockout only (no overexpression)
- Intelligence gene list is curated, not systematic GWAS mapping
- No statistical correction for multiple testing
- Embedding shift ≠ functional impact

## Next Steps
1. Scale to full dataset (2400+ cells) with GPU
2. Perturb all ~20K genes (genome-wide screen)
3. Use V2-316M model for better representations
4. Add overexpression perturbation mode
5. Systematic GWAS-to-gene mapping with MAGMA
6. Drug repurposing analysis with ChEMBL binding data
"""

with open(os.path.join(RESULTS_DIR, "pilot_report.md"), "w") as f:
    f.write(report)

print(f"  Report saved to {os.path.join(RESULTS_DIR, 'pilot_report.md')}")
print(f"  Time: {time.time()-t0:.1f}s")

print("\n" + "=" * 60)
print("PILOT COMPLETE")
print("=" * 60)
