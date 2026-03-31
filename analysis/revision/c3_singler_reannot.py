#!/usr/bin/env python3
"""
C3: Cell-type annotation sensitivity analysis using SingleR via rpy2.

Two modes:
  A) If raw h5ad is available: load it, run SingleR on real counts.
  B) If not: reconstruct a pseudo-expression matrix from Geneformer tokens.
     Geneformer ranks genes by expression (highest first), so token
     position encodes relative expression rank.  We invert ranks to
     create a pseudo-count matrix (rank 1 = highest pseudo-count).
     This is imperfect but lets us run SingleR as a sensitivity check.

Requirements:
  - R 4.x  with BiocManager, SingleR, celldex, SingleCellExperiment
  - Python: rpy2 (== 3.5.x for R 4.3), scanpy, anndata, datasets, numpy, pandas

Usage:
  python c3_singler_reannot.py                         # auto-detect mode
  python c3_singler_reannot.py --h5ad /path/to.h5ad    # force h5ad mode
"""

import argparse
import json
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "c3_results"

TOKENIZED_DATASET = (
    "/Volumes/Crucial X6/MacBook/Code/biodyn/"
    "neuro-mechinterp/analysis/results/insilico_wsl/brain.dataset"
)
GENE_DICT_DIR = "/tmp/geneformer_repo/geneformer/gene_dictionaries_30m"

# Candidate paths for raw h5ad (checked in order)
H5AD_CANDIDATES = [
    SCRIPT_DIR / "DLPFC_11k.h5ad",
    Path.home() / "data" / "DLPFC_11k.h5ad",
    Path("/tmp/DLPFC_11k.h5ad"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def find_h5ad(override: str | None = None) -> Path | None:
    """Return path to h5ad if it exists, else None."""
    if override:
        p = Path(override)
        return p if p.is_file() else None
    for p in H5AD_CANDIDATES:
        if p.is_file():
            print(f"[INFO] Found raw h5ad at {p}")
            return p
    return None


def load_geneformer_dicts():
    """Load Geneformer gene dictionaries.

    Returns
    -------
    token_to_ensembl : dict[int, str]
        Maps token ID -> Ensembl gene ID.
    ensembl_to_symbol : dict[str, str]
        Maps Ensembl ID -> gene symbol.
    """
    with open(os.path.join(GENE_DICT_DIR, "token_dictionary_gc30M.pkl"), "rb") as f:
        gene_to_token = pickle.load(f)  # ensembl -> token_id
    token_to_ensembl = {int(v): k for k, v in gene_to_token.items()
                        if k not in ("<pad>", "<mask>")}

    with open(os.path.join(GENE_DICT_DIR, "gene_name_id_dict_gc30M.pkl"), "rb") as f:
        name_to_id = pickle.load(f)  # symbol -> ensembl
    ensembl_to_symbol = {v: k for k, v in name_to_id.items()}

    return token_to_ensembl, ensembl_to_symbol


def build_pseudo_expression_from_tokens():
    """Build a pseudo-expression AnnData from the tokenized HF dataset.

    Geneformer tokenisation sorts genes by descending expression rank.
    Token at position 0 = most highly expressed gene in that cell.
    We assign pseudo-counts as: max_rank - position (so top gene gets
    the highest value).

    Returns
    -------
    adata : anndata.AnnData
        cells x genes, with .obs['cell_type_original'] from the dataset.
    """
    from datasets import load_from_disk
    import anndata as ad
    import scipy.sparse as sp

    print("[INFO] Loading tokenized dataset ...")
    ds = load_from_disk(TOKENIZED_DATASET)

    token_to_ensembl, ensembl_to_symbol = load_geneformer_dicts()

    # Collect all genes across cells to define the var axis
    all_ensembl = set()
    for row in ds:
        for tid in row["input_ids"][:row["length"]]:
            ens = token_to_ensembl.get(int(tid))
            if ens:
                all_ensembl.add(ens)
    all_ensembl = sorted(all_ensembl)
    gene_idx = {g: i for i, g in enumerate(all_ensembl)}
    n_genes = len(all_ensembl)
    n_cells = len(ds)

    print(f"[INFO] Constructing pseudo-expression matrix: "
          f"{n_cells} cells x {n_genes} genes")

    # Build sparse matrix (lil for row-wise construction)
    mat = sp.lil_matrix((n_cells, n_genes), dtype=np.float32)

    for ci, row in enumerate(ds):
        length = row["length"]
        tokens = row["input_ids"][:length]
        for pos, tid in enumerate(tokens):
            ens = token_to_ensembl.get(int(tid))
            if ens and ens in gene_idx:
                # Invert rank: position 0 (highest expressed) gets
                # the largest pseudo-count
                mat[ci, gene_idx[ens]] = float(length - pos)

    mat = mat.tocsr()

    # Map Ensembl IDs to gene symbols where possible
    symbols = [ensembl_to_symbol.get(e, e) for e in all_ensembl]

    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame({
            "cell_type_original": ds["cell_type"],
        }),
        var=pd.DataFrame({
            "ensembl_id": all_ensembl,
            "gene_symbol": symbols,
        }, index=symbols),
    )

    # Handle duplicate gene symbols (some Ensembl IDs map to same symbol)
    adata.var_names_make_unique()

    print(f"[INFO] Pseudo-expression AnnData: {adata.shape}")
    return adata


def load_raw_h5ad(path: Path):
    """Load raw h5ad and return AnnData."""
    import scanpy as sc

    print(f"[INFO] Loading raw h5ad: {path}")
    adata = sc.read_h5ad(str(path))
    print(f"[INFO] Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


# ---------------------------------------------------------------------------
# SingleR via rpy2
# ---------------------------------------------------------------------------
def run_singler(adata, use_ensembl: bool = False):
    """Run SingleR annotation via rpy2.

    Parameters
    ----------
    adata : AnnData
        Expression matrix (cells x genes).  Gene names should be symbols
        if use_ensembl=False, or Ensembl IDs if use_ensembl=True.
    use_ensembl : bool
        If True, use Ensembl IDs for gene matching (requires the reference
        to be queried with Ensembl IDs -- celldex uses symbols, so normally
        False).

    Returns
    -------
    pd.DataFrame with columns: labels, pruned_labels, max_score, delta
    """
    # Ensure ExperimentHub cache dir exists (avoids interactive prompt)
    import os as _os
    _eh_cache = _os.path.join(
        _os.path.expanduser("~"),
        "Library", "Caches", "org.R-project.R", "R", "ExperimentHub",
    )
    _os.makedirs(_eh_cache, exist_ok=True)

    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    # Set non-interactive so BiocFileCache / ExperimentHub don't prompt
    ro.r('options(BiocFileCache_interactive = FALSE)')
    ro.r('Sys.setenv(EXPERIMENT_HUB_CACHE = "{}")'.format(_eh_cache))

    importr("SingleR")
    importr("celldex")
    importr("SingleCellExperiment")

    # Extract dense matrix (SingleR needs genes x cells)
    if hasattr(adata.X, "toarray"):
        expr = adata.X.toarray()
    else:
        expr = np.asarray(adata.X)

    # genes x cells (R convention)
    expr_t = expr.T.astype(np.float64)
    gene_names = list(adata.var_names)

    n_genes, n_cells = expr_t.shape
    print(f"[INFO] Running SingleR on {n_cells} cells x {n_genes} genes ...")

    # Transfer to R
    ro.globalenv["expr_mat"] = ro.r["matrix"](
        numpy2ri.numpy2rpy(expr_t),
        nrow=n_genes,
        ncol=n_cells,
    )
    ro.globalenv["gene_names"] = ro.StrVector(gene_names)

    # Run SingleR in R
    ro.r("""
    rownames(expr_mat) <- gene_names

    # Use HumanPrimaryCellAtlasData as reference -- the broadest human
    # reference in celldex, includes neurons, glia, etc.
    message("[R] Fetching HumanPrimaryCellAtlasData reference ...")
    ref <- celldex::HumanPrimaryCellAtlasData()

    message("[R] Running SingleR ...")
    results <- SingleR::SingleR(
        test   = expr_mat,
        ref    = ref,
        labels = ref$label.main,
        de.method = "wilcox"
    )

    singler_labels  <- results$labels
    singler_pruned  <- results$pruned.labels
    singler_scores  <- apply(results$scores, 1, max)
    singler_delta   <- results$delta.next
    """)

    labels = list(ro.globalenv["singler_labels"])
    pruned = list(ro.globalenv["singler_pruned"])
    scores = list(ro.globalenv["singler_scores"])
    delta = list(ro.globalenv["singler_delta"])

    numpy2ri.deactivate()
    pandas2ri.deactivate()

    df = pd.DataFrame({
        "singler_label": labels,
        "singler_pruned": pruned,
        "singler_max_score": scores,
        "singler_delta": delta,
    })
    print(f"[INFO] SingleR finished. Unique labels: {df['singler_label'].nunique()}")
    return df


# ---------------------------------------------------------------------------
# Comparison & reporting
# ---------------------------------------------------------------------------
BROAD_TYPE_MAP = {
    # Map common SingleR HumanPrimaryCellAtlas labels -> our 6 broad types
    "Neurons": "neuron",
    "Neuroepithelial_cell": "neuron",
    "Astrocyte": "astrocyte",
    "Astrocytes": "astrocyte",
    "Oligodendrocytes": "oligodendrocyte",
    "Microglia": "microglial cell",
    "Macrophage": "macrophage",
    "Monocyte": "macrophage",
    "Smooth_muscle_cells": "mural cell",
    "Endothelial_cells": "endothelial cell",
    "Fibroblasts": "mural cell",
    "Tissue_stem_cells": "other",
    "T_cells": "immune",
    "B_cell": "immune",
    "NK_cell": "immune",
    "DC": "immune",
    "HSC_-G-CSF": "immune",
    "Platelets": "other",
    "Erythroblast": "other",
}

ORIGINAL_BROAD_MAP = {
    "oligodendrocyte": "oligodendrocyte",
    "glutamatergic neuron": "neuron",
    "GABAergic interneuron": "neuron",
    "pvalb GABAergic cortical interneuron": "neuron",
    "astrocyte": "astrocyte",
    "oligodendrocyte precursor cell": "oligodendrocyte",
    "microglial cell": "microglial cell",
    "mural cell": "mural cell",
    "macrophage": "macrophage",
}


def broad_label(label: str, mapping: dict) -> str:
    """Map a label to its broad category, falling back to 'other'."""
    return mapping.get(label, "other")


def compare_annotations(original_labels, singler_df):
    """Compare original vs SingleR annotations at broad-type level."""
    df = singler_df.copy()
    df["original"] = original_labels
    df["original_broad"] = df["original"].map(
        lambda x: broad_label(x, ORIGINAL_BROAD_MAP)
    )
    df["singler_broad"] = df["singler_label"].map(
        lambda x: broad_label(x, BROAD_TYPE_MAP)
    )

    # Exact broad-type agreement
    df["agree_broad"] = df["original_broad"] == df["singler_broad"]
    agreement_rate = df["agree_broad"].mean()

    # Confusion matrix (broad)
    confusion_broad = pd.crosstab(
        df["original_broad"], df["singler_broad"],
        rownames=["marker_based"], colnames=["singler"],
    )

    # Confusion matrix (fine-grained)
    confusion_fine = pd.crosstab(
        df["original"], df["singler_label"],
        rownames=["marker_based"], colnames=["singler"],
    )

    # Per-type agreement
    per_type = df.groupby("original_broad")["agree_broad"].mean()

    return df, confusion_broad, confusion_fine, agreement_rate, per_type


def write_report(df, confusion_broad, confusion_fine, agreement_rate,
                 per_type, mode, output_dir):
    """Write results and markdown report."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # CSVs
    df.to_csv(output_dir / "annotation_comparison.csv", index=False)
    confusion_broad.to_csv(output_dir / "confusion_matrix_broad.csv")
    confusion_fine.to_csv(output_dir / "confusion_matrix_fine.csv")

    # Summary JSON
    summary = {
        "mode": mode,
        "n_cells": len(df),
        "overall_agreement_broad": round(agreement_rate, 4),
        "per_type_agreement": {k: round(v, 4) for k, v in per_type.items()},
        "singler_label_distribution": df["singler_label"].value_counts().to_dict(),
        "original_label_distribution": df["original"].value_counts().to_dict(),
    }
    with open(output_dir / "c3_singler_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Markdown report
    lines = [
        "# C3: SingleR Re-annotation Sensitivity Analysis",
        "",
        f"**Mode:** {'Raw h5ad counts' if mode == 'h5ad' else 'Pseudo-expression from Geneformer tokens'}",
        f"**Cells:** {len(df)}",
        f"**Overall broad-type agreement:** {agreement_rate:.1%}",
        "",
        "## Per-type agreement (broad categories)",
        "",
        "| Broad type | Agreement |",
        "|---|---|",
    ]
    for ct in sorted(per_type.index):
        lines.append(f"| {ct} | {per_type[ct]:.1%} |")

    lines += [
        "",
        "## Confusion matrix (broad categories)",
        "",
        confusion_broad.to_markdown(),
        "",
        "## SingleR label distribution",
        "",
        "| SingleR label | Count |",
        "|---|---|",
    ]
    for lab, cnt in df["singler_label"].value_counts().items():
        lines.append(f"| {lab} | {cnt} |")

    if mode == "pseudo":
        lines += [
            "",
            "## Caveats",
            "",
            "This analysis used **pseudo-expression values** reconstructed from",
            "Geneformer token ranks, NOT raw UMI counts.  Geneformer tokenisation",
            "rank-orders genes by expression within each cell, so we inverted the",
            "ranks to produce ordinal pseudo-counts.  This preserves the relative",
            "ordering of genes but not the magnitude or distributional shape of",
            "real scRNA-seq counts.  Results should be interpreted as a directional",
            "sensitivity check, not a definitive re-annotation.",
        ]

    lines += [
        "",
        "## Interpretation for manuscript",
        "",
        f"The broad-type agreement between marker-based annotation and SingleR",
        f"(HumanPrimaryCellAtlasData reference) is **{agreement_rate:.1%}**.",
    ]
    if agreement_rate >= 0.80:
        lines.append(
            "This confirms that marker-based annotations are largely consistent "
            "with reference-based classification at the broad-type level, "
            "supporting the robustness of our cell-type-stratified analyses."
        )
    else:
        lines.append(
            "Agreement is moderate, suggesting that some cell-type-specific "
            "findings may be sensitive to annotation method. Detailed "
            "per-type concordance should be examined."
        )

    report_path = output_dir / "c3_singler_report.md"
    report_path.write_text("\n".join(lines) + "\n")
    print(f"[INFO] Report written to {report_path}")
    print(f"[INFO] CSVs and JSON written to {output_dir}/")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="C3: SingleR cell-type re-annotation sensitivity analysis"
    )
    parser.add_argument(
        "--h5ad", type=str, default=None,
        help="Path to raw DLPFC h5ad file (if available)"
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(OUTPUT_DIR),
        help="Output directory for results"
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    print("=" * 65)
    print("  C3: Cell-Type Annotation Sensitivity — SingleR Re-annotation")
    print("=" * 65)

    # Decide mode
    h5ad_path = find_h5ad(args.h5ad)

    if h5ad_path is not None:
        mode = "h5ad"
        adata = load_raw_h5ad(h5ad_path)
        original_labels = list(adata.obs.get(
            "cell_type",
            adata.obs.get("celltype", adata.obs.iloc[:, 0])
        ))
    else:
        mode = "pseudo"
        print("[INFO] No raw h5ad found. Using pseudo-expression from tokens.")
        adata = build_pseudo_expression_from_tokens()
        original_labels = list(adata.obs["cell_type_original"])

    # Log-normalise (SingleR expects log-normalised data)
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print("[INFO] Data log-normalised.")

    # Run SingleR
    singler_df = run_singler(adata, use_ensembl=False)

    # Compare
    (comparison_df, confusion_broad, confusion_fine,
     agreement_rate, per_type) = compare_annotations(original_labels, singler_df)

    # Report
    write_report(
        comparison_df, confusion_broad, confusion_fine,
        agreement_rate, per_type, mode, output_dir,
    )

    print(f"\n{'=' * 65}")
    print(f"  Overall broad-type agreement: {agreement_rate:.1%}")
    print(f"  Mode: {mode}")
    print(f"  Results in: {output_dir}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
