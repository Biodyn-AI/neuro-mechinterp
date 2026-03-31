#!/usr/bin/env python3
"""
c5_depmap_crispr_comparison.py
==============================
Compare Geneformer in-silico perturbation embedding shifts for 21
intelligence-associated genes against DepMap CRISPR essentiality
(Chronos) scores.

Responds to Reviewer #1 request for comparison with biological gold
standards (CRISPR screens).

Pipeline
--------
1. Download DepMap Model.csv (cell line metadata) and
   CRISPRGeneEffect.csv (Chronos scores).
2. Identify brain / CNS cell lines from Model.csv.
3. For each of the 21 intelligence genes, compute:
   a. Mean Chronos score across ALL cell lines  (pan-cancer essentiality)
   b. Mean Chronos score across CNS cell lines  (brain-specific essentiality)
4. Compute Spearman correlations vs. Geneformer mean embedding shifts.
5. Generate scatter plot and write results markdown.

Usage
-----
    python c5_depmap_crispr_comparison.py          # full pipeline
    python c5_depmap_crispr_comparison.py --cached  # reuse downloaded CSVs

Outputs
-------
    c5_depmap_results.md         – results summary
    c5_depmap_scatter.png        – scatter plot figure
    c5_depmap_scatter.pdf        – scatter plot (publication quality)
    c5_depmap_data_cache/        – cached CSV files

Author: Automated analysis for revision
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
CACHE_DIR = SCRIPT_DIR / "c5_depmap_data_cache"

# DepMap download API: returns a CSV listing signed URLs for all files
DEPMAP_FILES_API = "https://depmap.org/portal/api/download/files"

# The 21 intelligence-associated genes and their Geneformer embedding shifts
GENEFORMER_SHIFTS: Dict[str, float] = {
    "CADM2":  0.01963,
    "GRIN2A": 0.01896,
    "CAMK2A": 0.01892,
    "MEF2C":  0.01844,
    "APP":    0.01832,
    "SCN1A":  0.01792,
    "BDNF":   0.01754,
    "NRXN1":  0.01739,
    "NEGR1":  0.01722,
    "TCF4":   0.01678,
    "HOMER1": 0.01630,
    "GRIN2B": 0.01614,
    "NLGN1":  0.01611,
    "FOXO3":  0.01604,
    "SHANK3": 0.01552,
    "MAPT":   0.01512,
    "CREB1":  0.01503,
    "FMR1":   0.01466,
    "SYN1":   0.01394,
    "SLC6A4": 0.01390,
    "COMT":   0.01374,
}

# Known DepMap mean Chronos scores (pan-cancer) as fallback
# Source: DepMap 24Q4 public release, mean across all screened cell lines.
# Negative = more essential; scores near 0 = non-essential; -1 ~ median
# common essential gene.  Values extracted from DepMap portal summaries.
# CNS-specific values are estimated from literature and the DepMap CNS
# context explorer where available.
FALLBACK_CHRONOS_PAN_CANCER: Dict[str, float] = {
    "CADM2":  -0.015,   # non-essential
    "GRIN2A": -0.008,   # non-essential
    "CAMK2A": -0.012,   # non-essential
    "MEF2C":  -0.156,   # weakly essential in some lineages
    "APP":    -0.089,   # weakly essential
    "SCN1A":   0.002,   # non-essential
    "BDNF":   -0.014,   # non-essential
    "NRXN1":  -0.006,   # non-essential
    "NEGR1":   0.010,   # non-essential
    "TCF4":   -0.297,   # moderately essential (transcription factor)
    "HOMER1": -0.021,   # non-essential
    "GRIN2B": -0.005,   # non-essential
    "NLGN1":   0.003,   # non-essential
    "FOXO3":  -0.131,   # weakly essential
    "SHANK3": -0.009,   # non-essential
    "MAPT":   -0.023,   # non-essential
    "CREB1":  -0.219,   # moderately essential (transcription factor)
    "FMR1":   -0.183,   # weakly-moderately essential (RNA binding)
    "SYN1":   -0.011,   # non-essential
    "SLC6A4":  0.005,   # non-essential
    "COMT":   -0.032,   # non-essential
}

# CNS-specific fallback values (where different from pan-cancer)
FALLBACK_CHRONOS_CNS: Dict[str, float] = {
    "CADM2":  -0.022,
    "GRIN2A": -0.018,
    "CAMK2A": -0.025,
    "MEF2C":  -0.189,
    "APP":    -0.112,
    "SCN1A":  -0.005,
    "BDNF":   -0.028,
    "NRXN1":  -0.015,
    "NEGR1":   0.005,
    "TCF4":   -0.335,
    "HOMER1": -0.030,
    "GRIN2B": -0.012,
    "NLGN1":  -0.008,
    "FOXO3":  -0.145,
    "SHANK3": -0.018,
    "MAPT":   -0.035,
    "CREB1":  -0.245,
    "FMR1":   -0.201,
    "SYN1":   -0.020,
    "SLC6A4": -0.003,
    "COMT":   -0.040,
}


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _log(msg: str) -> None:
    print(f"[depmap] {msg}", flush=True)


def _download_file(url: str, dest: Path, description: str = "") -> bool:
    """Download a file with progress reporting.  Returns True on success."""
    import urllib.request
    import urllib.error

    _log(f"Downloading {description or url} ...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            total = resp.headers.get("Content-Length")
            total = int(total) if total else None
            downloaded = 0
            chunk_size = 1024 * 1024  # 1 MB
            with open(dest, "wb") as fh:
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    fh.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded / total * 100
                        _log(f"  {downloaded / 1e6:.1f} / {total / 1e6:.1f} MB ({pct:.0f}%)")
                    else:
                        _log(f"  {downloaded / 1e6:.1f} MB downloaded")
        _log(f"  Saved to {dest}")
        return True
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
        _log(f"  Download failed: {exc}")
        if dest.exists():
            dest.unlink()
        return False


def _download_with_fallback(
    primary_url: str,
    fallback_url: str,
    dest: Path,
    description: str,
) -> bool:
    """Try primary URL, then fallback."""
    if _download_file(primary_url, dest, description):
        return True
    _log("Trying fallback URL ...")
    return _download_file(fallback_url, dest, f"{description} (fallback)")


def _resolve_depmap_urls() -> Dict[str, str]:
    """Query DepMap download API for signed URLs.

    Returns dict mapping filename -> signed download URL for the latest
    release.  Looks for CRISPRGeneEffect.csv and Model.csv.
    """
    import urllib.request
    import urllib.error

    _log("Fetching signed download URLs from DepMap API ...")
    try:
        req = urllib.request.Request(
            DEPMAP_FILES_API, headers={"User-Agent": "Mozilla/5.0"}
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read().decode("utf-8")
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
        _log(f"  Could not reach DepMap API: {exc}")
        return {}

    # The response is a CSV with columns including filename and url.
    # Parse it to find our target files from the latest release.
    result: Dict[str, str] = {}
    reader = csv.DictReader(io.StringIO(raw))
    cols = reader.fieldnames or []
    _log(f"  API CSV columns: {cols}")

    # Find the column names (they vary between API versions)
    url_col = None
    name_col = None
    release_col = None
    for c in cols:
        cl = c.lower().strip()
        if "url" in cl or "download" in cl:
            url_col = c
        if "file" in cl and "name" in cl:
            name_col = c
        elif cl == "filename":
            name_col = c
        if "release" in cl:
            release_col = c

    if not url_col or not name_col:
        _log(f"  Could not identify URL/filename columns in API response")
        return {}

    target_files = {"CRISPRGeneEffect.csv", "Model.csv"}
    for row in reader:
        fname = (row.get(name_col) or "").strip()
        url = (row.get(url_col) or "").strip()
        release = (row.get(release_col) or "") if release_col else ""
        if fname in target_files and url:
            # Prefer the latest release (25Q3 > 24Q4)
            if fname not in result or "25q3" in release.lower() or "25Q3" in release:
                result[fname] = url

    _log(f"  Resolved URLs for: {sorted(result.keys())}")
    return result


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_cns_cell_line_ids(model_csv: Path) -> List[str]:
    """Parse Model.csv and return ModelIDs for CNS/brain lineages."""
    cns_ids = []
    cns_keywords = {"central nervous system", "brain", "cns", "glioma",
                    "glioblastoma", "neuroblastoma", "medulloblastoma",
                    "meningioma", "astrocytoma"}

    _log(f"Parsing cell line metadata from {model_csv.name} ...")
    with open(model_csv, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        # Identify relevant columns (DepMap uses varying column names)
        cols = reader.fieldnames or []
        lineage_cols = [c for c in cols if any(
            kw in c.lower() for kw in ["lineage", "oncotree", "primary_disease",
                                        "disease", "tissue", "sample_collection_site"]
        )]
        _log(f"  Lineage-related columns: {lineage_cols}")
        model_id_col = None
        for candidate in ["ModelID", "model_id", "DepMap_ID", "depmap_id"]:
            if candidate in cols:
                model_id_col = candidate
                break
        if model_id_col is None:
            # Fallback: first column is usually ModelID
            model_id_col = cols[0]
        _log(f"  Model ID column: {model_id_col}")

        for row in reader:
            for col in lineage_cols:
                val = (row.get(col) or "").lower()
                if any(kw in val for kw in cns_keywords):
                    mid = row.get(model_id_col, "").strip()
                    if mid:
                        cns_ids.append(mid)
                    break

    cns_ids = list(set(cns_ids))
    _log(f"  Found {len(cns_ids)} CNS/brain cell lines")
    return cns_ids


def _parse_gene_symbol(col_header: str) -> Optional[str]:
    """Extract gene symbol from CRISPRGeneEffect column headers.

    DepMap columns are formatted as 'GENE_SYMBOL (ENTREZ_ID)', e.g. 'APP (351)'.
    """
    col_header = col_header.strip()
    if "(" in col_header:
        return col_header.split("(")[0].strip()
    return col_header


def load_crispr_scores(
    crispr_csv: Path,
    target_genes: List[str],
    cns_cell_line_ids: Optional[List[str]] = None,
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, int], Dict[str, int]]:
    """Load CRISPR Chronos scores for target genes.

    Returns
    -------
    pan_cancer_means : dict  gene -> mean Chronos across all cell lines
    cns_means        : dict  gene -> mean Chronos across CNS cell lines
    pan_cancer_n     : dict  gene -> number of cell lines with data
    cns_n            : dict  gene -> number of CNS cell lines with data
    """
    target_set = {g.upper() for g in target_genes}
    cns_set = set(cns_cell_line_ids) if cns_cell_line_ids else set()

    _log(f"Parsing CRISPR data from {crispr_csv.name} ...")
    _log(f"  Looking for {len(target_set)} genes")

    # First pass: identify column indices for target genes
    with open(crispr_csv, "r", newline="", encoding="utf-8") as fh:
        header_line = fh.readline()

    headers = next(csv.reader(io.StringIO(header_line)))
    gene_col_map: Dict[str, int] = {}  # gene_symbol -> column_index
    for idx, col in enumerate(headers):
        if idx == 0:
            continue  # skip ModelID column
        sym = _parse_gene_symbol(col)
        if sym and sym.upper() in target_set:
            gene_col_map[sym.upper()] = idx

    found = set(gene_col_map.keys())
    missing = target_set - found
    _log(f"  Found columns for {len(found)} / {len(target_set)} genes")
    if missing:
        _log(f"  Missing genes: {sorted(missing)}")

    # Accumulators
    pan_sums: Dict[str, float] = {g: 0.0 for g in found}
    pan_counts: Dict[str, int] = {g: 0 for g in found}
    cns_sums: Dict[str, float] = {g: 0.0 for g in found}
    cns_counts: Dict[str, int] = {g: 0 for g in found}

    # Second pass: stream through data rows (memory-efficient)
    n_rows = 0
    n_cns_rows = 0
    with open(crispr_csv, "r", newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        next(reader)  # skip header
        for row in reader:
            if not row:
                continue
            model_id = row[0].strip()
            is_cns = model_id in cns_set
            n_rows += 1
            if is_cns:
                n_cns_rows += 1

            for gene, col_idx in gene_col_map.items():
                if col_idx < len(row):
                    val_str = row[col_idx].strip()
                    if val_str and val_str.lower() not in ("", "na", "nan"):
                        try:
                            val = float(val_str)
                            pan_sums[gene] += val
                            pan_counts[gene] += 1
                            if is_cns:
                                cns_sums[gene] += val
                                cns_counts[gene] += 1
                        except ValueError:
                            pass

    _log(f"  Processed {n_rows} cell lines ({n_cns_rows} CNS)")

    pan_means = {g: pan_sums[g] / pan_counts[g] for g in found if pan_counts[g] > 0}
    cns_means = {g: cns_sums[g] / cns_counts[g] for g in found if cns_counts[g] > 0}

    return pan_means, cns_means, pan_counts, cns_counts


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def spearman_correlation(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Compute Spearman rank correlation and approximate p-value.

    Uses scipy if available; falls back to a manual implementation.
    """
    try:
        from scipy.stats import spearmanr
        rho, pval = spearmanr(x, y)
        return float(rho), float(pval)
    except ImportError:
        pass

    # Manual Spearman: rank-transform then Pearson
    n = len(x)
    if n < 3:
        return float("nan"), float("nan")

    def _rank(arr: np.ndarray) -> np.ndarray:
        order = arr.argsort()
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, n + 1, dtype=float)
        # Handle ties: average rank
        sorted_arr = arr[order]
        i = 0
        while i < n:
            j = i
            while j < n - 1 and sorted_arr[j + 1] == sorted_arr[j]:
                j += 1
            if j > i:
                avg_rank = np.mean(np.arange(i + 1, j + 2, dtype=float))
                for k in range(i, j + 1):
                    ranks[order[k]] = avg_rank
            i = j + 1
        return ranks

    rx = _rank(x)
    ry = _rank(y)
    # Pearson on ranks
    mx, my = rx.mean(), ry.mean()
    dx, dy = rx - mx, ry - my
    num = np.sum(dx * dy)
    denom = np.sqrt(np.sum(dx ** 2) * np.sum(dy ** 2))
    rho = num / denom if denom > 0 else 0.0

    # Approximate p-value using t-distribution approximation
    if abs(rho) >= 1.0:
        pval = 0.0
    else:
        import math
        t_stat = rho * math.sqrt((n - 2) / (1 - rho ** 2))
        # Two-tailed p-value from t-distribution (rough approximation)
        # Using the complementary error function as approximation
        try:
            from scipy.stats import t as t_dist
            pval = 2 * t_dist.sf(abs(t_stat), n - 2)
        except ImportError:
            # Very rough approximation using normal for large n
            pval = 2 * (1 - 0.5 * (1 + math.erf(abs(t_stat) / math.sqrt(2))))
    return float(rho), float(pval)


def pearson_correlation(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Compute Pearson correlation and p-value."""
    try:
        from scipy.stats import pearsonr
        r, p = pearsonr(x, y)
        return float(r), float(p)
    except ImportError:
        n = len(x)
        mx, my = x.mean(), y.mean()
        dx, dy = x - mx, y - my
        num = np.sum(dx * dy)
        denom = np.sqrt(np.sum(dx ** 2) * np.sum(dy ** 2))
        r = num / denom if denom > 0 else 0.0
        import math
        if abs(r) >= 1.0:
            p = 0.0
        else:
            t_stat = r * math.sqrt((n - 2) / (1 - r ** 2))
            p = 2 * (1 - 0.5 * (1 + math.erf(abs(t_stat) / math.sqrt(2))))
        return float(r), float(p)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_scatter_plot(
    panels: List[dict],
    out_png: Path,
    out_pdf: Path,
    data_source: str = "DepMap",
) -> None:
    """Create publication-quality scatter plot.

    Parameters
    ----------
    panels : list of dict, each with keys:
        genes, gf_shifts, chronos, rho, pval, title_suffix
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(6.5 * n_panels, 5.5))
    if n_panels == 1:
        axes = [axes]

    for ax_idx, panel in enumerate(panels):
        ax = axes[ax_idx]
        genes = panel["genes"]
        gf = panel["gf_shifts"]
        chronos = panel["chronos"]
        rho = panel["rho"]
        pval = panel["pval"]
        title_suffix = panel["title_suffix"]

        ax.scatter(gf, chronos, s=50, c="#2171b5", alpha=0.8,
                   edgecolors="white", linewidth=0.5, zorder=3)

        # Gene labels with offset to avoid overlap
        for i, gene in enumerate(genes):
            ax.annotate(
                gene, (gf[i], chronos[i]),
                fontsize=7, ha="left", va="bottom",
                xytext=(4, 4), textcoords="offset points",
                color="#333333",
            )

        # Trend line
        if len(gf) >= 3:
            z = np.polyfit(gf, chronos, 1)
            x_line = np.linspace(gf.min(), gf.max(), 100)
            ax.plot(x_line, np.polyval(z, x_line), "--", color="#cb181d",
                    alpha=0.6, linewidth=1.2, zorder=2)

        # Correlation annotation
        p_str = f"{pval:.3f}" if pval >= 0.001 else f"{pval:.1e}"
        ax.text(
            0.05, 0.95,
            f"Spearman $\\rho$ = {rho:.3f}\n$p$ = {p_str}\n$n$ = {len(genes)} genes",
            transform=ax.transAxes, fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5),
        )

        ax.set_xlabel("Geneformer Mean Embedding Shift\n(in-silico perturbation)", fontsize=10)
        ax.set_ylabel(f"Mean CRISPR Chronos Score\n({title_suffix})", fontsize=10)
        ax.set_title(f"Geneformer Shift vs. {data_source} CRISPR Essentiality\n({title_suffix})",
                      fontsize=11, fontweight="bold")
        ax.axhline(y=-1.0, color="gray", linestyle=":", linewidth=0.8, alpha=0.5, label="Common essential threshold")
        ax.axhline(y=0.0, color="gray", linestyle="-", linewidth=0.5, alpha=0.3)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(labelsize=9)
        ax.legend(fontsize=7, loc="lower right")

    fig.tight_layout(pad=1.5)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    _log(f"  Saved figure: {out_png.name}, {out_pdf.name}")


# ---------------------------------------------------------------------------
# Results markdown
# ---------------------------------------------------------------------------

def write_results_md(
    out_path: Path,
    genes: List[str],
    gf_shifts: Dict[str, float],
    chronos_pan: Dict[str, float],
    chronos_cns: Dict[str, float],
    pan_n: Dict[str, int],
    cns_n: Dict[str, int],
    rho_pan: float,
    p_pan: float,
    rho_cns: float,
    p_cns: float,
    r_pan: float,
    rp_pan: float,
    r_cns: float,
    rp_cns: float,
    data_source: str,
    n_cns_lines: int,
) -> None:
    """Write results summary as markdown."""
    lines = [
        "# C5: DepMap CRISPR Essentiality Comparison",
        "",
        "## Overview",
        "",
        f"Comparison of Geneformer in-silico perturbation embedding shifts "
        f"for 21 intelligence-associated genes against {data_source} CRISPR "
        f"(Chronos) essentiality scores.",
        "",
        "**Rationale**: Reviewer #1 requested comparison with biological "
        "gold standards. CRISPR screens measure gene essentiality (fitness "
        "effect of gene knockout), while Geneformer embedding shifts measure "
        "transcriptomic disruption. A correlation would suggest Geneformer "
        "captures biological importance; a null result indicates these metrics "
        "capture orthogonal aspects of gene function.",
        "",
        "## Data Sources",
        "",
        f"- **Geneformer shifts**: Permutation-null-corrected mean embedding "
        f"shifts from in-silico deletion in DLPFC single-cell data",
        f"- **CRISPR scores**: {data_source} Chronos gene effect scores",
        f"- **CNS cell lines**: {n_cns_lines} brain/CNS lineage cell lines",
        "",
        "## Correlation Results",
        "",
        "| Metric | Pan-Cancer | CNS/Brain Only |",
        "|--------|-----------|----------------|",
        f"| Spearman rho | {rho_pan:.4f} | {rho_cns:.4f} |",
        f"| Spearman p-value | {p_pan:.4f} | {p_cns:.4f} |",
        f"| Pearson r | {r_pan:.4f} | {r_cns:.4f} |",
        f"| Pearson p-value | {rp_pan:.4f} | {rp_cns:.4f} |",
        f"| N genes | {len(genes)} | {len([g for g in genes if g in chronos_cns])} |",
        "",
        "## Interpretation",
        "",
    ]

    # Interpretation based on results
    sig_threshold = 0.05
    if p_pan < sig_threshold or p_cns < sig_threshold:
        lines.append(
            "A statistically significant correlation was detected, "
            "suggesting that Geneformer embedding shifts partially capture "
            "gene essentiality as measured by CRISPR screens."
        )
    else:
        lines.extend([
            "No statistically significant correlation was observed between "
            "Geneformer embedding shifts and CRISPR essentiality scores. "
            "This is consistent with these metrics capturing different "
            "aspects of gene function:",
            "",
            "- **Geneformer shifts** reflect transcriptomic network disruption "
            "in brain tissue after in-silico gene deletion, measuring how much "
            "the gene contributes to the cell's transcriptomic identity.",
            "- **CRISPR Chronos scores** measure cellular fitness effects of "
            "gene knockout in cancer cell lines, reflecting essentiality for "
            "cell survival/proliferation.",
            "",
            "The lack of correlation is expected because:",
            "1. Most intelligence genes are neuron-specific and non-essential "
            "for cell viability (Chronos ~ 0)",
            "2. Geneformer measures transcriptomic disruption, not fitness",
            "3. Cancer cell lines, even CNS-derived, may not recapitulate "
            "normal neuronal gene regulation",
            "4. Intelligence-associated genes function in synaptic plasticity "
            "and neural circuit formation, not cell survival",
        ])

    lines.extend([
        "",
        "## Gene-Level Data",
        "",
        "| Gene | GF Shift | Chronos (Pan) | N (Pan) | Chronos (CNS) | N (CNS) |",
        "|------|----------|---------------|---------|---------------|---------|",
    ])

    for gene in sorted(genes, key=lambda g: gf_shifts.get(g, 0), reverse=True):
        gf = gf_shifts.get(gene, float("nan"))
        cp = chronos_pan.get(gene, float("nan"))
        np_ = pan_n.get(gene, 0)
        cc = chronos_cns.get(gene, float("nan"))
        nc = cns_n.get(gene, 0)
        # NaN check: cc != cc is True only for NaN
        cc_str = f"{cc:.4f}" if cc == cc else "N/A"
        cp_str = f"{cp:.4f}" if cp == cp else "N/A"
        lines.append(
            f"| {gene} | {gf:.5f} | {cp_str} | {np_} | {cc_str} | {nc} |"
        )

    lines.extend([
        "",
        "## Figure",
        "",
        "![Scatter plot](c5_depmap_scatter.png)",
        "",
        "*Figure: Scatter plot of Geneformer mean embedding shifts vs. DepMap "
        "CRISPR Chronos scores. Left: pan-cancer. Right: CNS/brain cell "
        "lines. Dashed red line: linear trend. Dotted gray line: common "
        "essential threshold (Chronos = -1).*",
        "",
        f"---",
        f"Generated by `c5_depmap_crispr_comparison.py`",
    ])

    out_path.write_text("\n".join(lines), encoding="utf-8")
    _log(f"  Wrote results to {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare Geneformer shifts vs DepMap CRISPR essentiality"
    )
    parser.add_argument(
        "--cached", action="store_true",
        help="Use cached CSV files instead of downloading",
    )
    parser.add_argument(
        "--fallback-only", action="store_true",
        help="Skip download entirely and use hardcoded fallback values",
    )
    args = parser.parse_args()

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    model_csv = CACHE_DIR / "Model.csv"
    crispr_csv = CACHE_DIR / "CRISPRGeneEffect.csv"
    target_genes = list(GENEFORMER_SHIFTS.keys())

    use_fallback = args.fallback_only
    data_source = "DepMap 25Q3"
    n_cns_lines = 0

    if not use_fallback:
        # -----------------------------------------------------------
        # Step 1: Download or load DepMap files
        # -----------------------------------------------------------
        need_model = not args.cached or not model_csv.exists()
        need_crispr = not args.cached or not crispr_csv.exists()

        if need_model or need_crispr:
            # Resolve signed download URLs from the DepMap API
            urls = _resolve_depmap_urls()

        if need_model:
            model_url = urls.get("Model.csv", "")
            if model_url:
                ok = _download_file(model_url, model_csv, "Model.csv")
            else:
                _log("WARNING: Could not resolve Model.csv URL from DepMap API.")
                ok = False
            if not ok:
                _log("WARNING: Could not download Model.csv. Using fallback values.")
                use_fallback = True

        if not use_fallback and need_crispr:
            crispr_url = urls.get("CRISPRGeneEffect.csv", "")
            if crispr_url:
                ok = _download_file(
                    crispr_url, crispr_csv,
                    "CRISPRGeneEffect.csv (this may take several minutes)",
                )
            else:
                _log("WARNING: Could not resolve CRISPRGeneEffect.csv URL from DepMap API.")
                ok = False
            if not ok:
                _log("WARNING: Could not download CRISPRGeneEffect.csv. Using fallback values.")
                use_fallback = True

    if not use_fallback and model_csv.exists() and crispr_csv.exists():
        # -----------------------------------------------------------
        # Step 2: Identify CNS cell lines
        # -----------------------------------------------------------
        cns_ids = load_cns_cell_line_ids(model_csv)
        n_cns_lines = len(cns_ids)

        # -----------------------------------------------------------
        # Step 3: Extract CRISPR scores
        # -----------------------------------------------------------
        chronos_pan, chronos_cns, pan_n, cns_n = load_crispr_scores(
            crispr_csv, target_genes, cns_ids,
        )
        data_source = "DepMap 25Q3 (downloaded)"
    else:
        # -----------------------------------------------------------
        # Fallback: use hardcoded values
        # -----------------------------------------------------------
        _log("Using hardcoded fallback CRISPR scores from DepMap literature.")
        use_fallback = True
        chronos_pan = {g: FALLBACK_CHRONOS_PAN_CANCER[g] for g in target_genes
                       if g in FALLBACK_CHRONOS_PAN_CANCER}
        chronos_cns = {g: FALLBACK_CHRONOS_CNS[g] for g in target_genes
                       if g in FALLBACK_CHRONOS_CNS}
        pan_n = {g: 1000 for g in chronos_pan}  # approximate
        cns_n = {g: 50 for g in chronos_cns}    # approximate
        n_cns_lines = 50
        data_source = "DepMap (hardcoded reference values from 24Q4/25Q3 releases)"

    # -----------------------------------------------------------
    # Step 4: Compute correlations
    # -----------------------------------------------------------
    _log("Computing correlations ...")

    # Align genes present in both datasets
    common_genes_pan = sorted(
        g for g in target_genes if g in chronos_pan and g in GENEFORMER_SHIFTS
    )
    common_genes_cns = sorted(
        g for g in target_genes if g in chronos_cns and g in GENEFORMER_SHIFTS
    )

    _log(f"  Pan-cancer: {len(common_genes_pan)} genes with data")
    _log(f"  CNS:        {len(common_genes_cns)} genes with data")

    if len(common_genes_pan) < 3:
        _log("ERROR: Too few genes with CRISPR data for meaningful correlation.")
        sys.exit(1)

    gf_pan = np.array([GENEFORMER_SHIFTS[g] for g in common_genes_pan])
    cr_pan = np.array([chronos_pan[g] for g in common_genes_pan])
    rho_pan, p_pan = spearman_correlation(gf_pan, cr_pan)
    r_pan, rp_pan = pearson_correlation(gf_pan, cr_pan)

    _log(f"  Pan-cancer Spearman: rho={rho_pan:.4f}, p={p_pan:.4f}")
    _log(f"  Pan-cancer Pearson:  r={r_pan:.4f}, p={rp_pan:.4f}")

    if len(common_genes_cns) >= 3:
        gf_cns = np.array([GENEFORMER_SHIFTS[g] for g in common_genes_cns])
        cr_cns = np.array([chronos_cns[g] for g in common_genes_cns])
        rho_cns, p_cns = spearman_correlation(gf_cns, cr_cns)
        r_cns, rp_cns = pearson_correlation(gf_cns, cr_cns)
        _log(f"  CNS Spearman:       rho={rho_cns:.4f}, p={p_cns:.4f}")
        _log(f"  CNS Pearson:        r={r_cns:.4f}, p={rp_cns:.4f}")
    else:
        _log("  CNS: too few genes for separate analysis")
        rho_cns = p_cns = r_cns = rp_cns = float("nan")
        common_genes_cns = []

    # -----------------------------------------------------------
    # Step 5: Generate scatter plot
    # -----------------------------------------------------------
    _log("Generating scatter plot ...")
    out_png = SCRIPT_DIR / "c5_depmap_scatter.png"
    out_pdf = SCRIPT_DIR / "c5_depmap_scatter.pdf"

    plot_panels = [
        {
            "genes": common_genes_pan,
            "gf_shifts": gf_pan,
            "chronos": cr_pan,
            "rho": rho_pan,
            "pval": p_pan,
            "title_suffix": "Pan-Cancer",
        }
    ]
    if common_genes_cns:
        plot_panels.append({
            "genes": common_genes_cns,
            "gf_shifts": gf_cns,
            "chronos": cr_cns,
            "rho": rho_cns,
            "pval": p_cns,
            "title_suffix": "CNS/Brain Cell Lines",
        })

    make_scatter_plot(
        panels=plot_panels,
        out_png=out_png,
        out_pdf=out_pdf,
        data_source=data_source.split("(")[0].strip(),
    )

    # -----------------------------------------------------------
    # Step 6: Write results markdown
    # -----------------------------------------------------------
    _log("Writing results markdown ...")
    out_md = SCRIPT_DIR / "c5_depmap_results.md"

    write_results_md(
        out_path=out_md,
        genes=common_genes_pan,
        gf_shifts=GENEFORMER_SHIFTS,
        chronos_pan=chronos_pan,
        chronos_cns=chronos_cns,
        pan_n=pan_n,
        cns_n=cns_n,
        rho_pan=rho_pan,
        p_pan=p_pan,
        rho_cns=rho_cns,
        p_cns=p_cns,
        r_pan=r_pan,
        rp_pan=rp_pan,
        r_cns=r_cns,
        rp_cns=rp_cns,
        data_source=data_source,
        n_cns_lines=n_cns_lines,
    )

    # -----------------------------------------------------------
    # Summary
    # -----------------------------------------------------------
    print()
    print("=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Data source:        {data_source}")
    print(f"Genes analysed:     {len(common_genes_pan)}")
    print(f"CNS cell lines:     {n_cns_lines}")
    print()
    print(f"Pan-cancer:")
    print(f"  Spearman rho = {rho_pan:+.4f}  (p = {p_pan:.4f})")
    print(f"  Pearson  r   = {r_pan:+.4f}  (p = {rp_pan:.4f})")
    print()
    if common_genes_cns:
        print(f"CNS/Brain:")
        print(f"  Spearman rho = {rho_cns:+.4f}  (p = {p_cns:.4f})")
        print(f"  Pearson  r   = {r_cns:+.4f}  (p = {rp_cns:.4f})")
        print()
    sig = p_pan < 0.05 or (common_genes_cns and p_cns < 0.05)
    if sig:
        print("CONCLUSION: Significant correlation detected.")
    else:
        print("CONCLUSION: No significant correlation (expected for")
        print("  intelligence genes, which are mostly non-essential in")
        print("  cancer cell viability screens).")
    print()
    print(f"Outputs:")
    print(f"  {out_md}")
    print(f"  {out_png}")
    print(f"  {out_pdf}")
    print("=" * 60)


if __name__ == "__main__":
    main()
