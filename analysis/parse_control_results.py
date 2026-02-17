exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
#!/usr/bin/env python3
"""
Parse Geneformer perturbation results for control genes and compare
with intelligence gene results. Compute Mann-Whitney U test and Cohen's d.
"""
import os, pickle, json
import numpy as np
from scipy import stats

BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"

# Intelligence gene results from final_rankings.json
INTEL_RESULTS_FILE = os.path.join(BASE, "final_rankings.json")

# Matching info
MATCHED_TO = {
    "EXOC6B": "CADM2", "DOCK9": "APP", "PDE7B": "FOXO3", "ZDHHC9": "NLGN1",
    "HDAC11": "MAPT", "XRCC6": "TCF4", "WDR46": "CREB1", "CUEDC2": "DISC1",
    "IGSF10": "NEGR1", "SALL3": "FMR1", "MFSD6L": "COMT", "TMEM88": "MEF2C",
    "ACTB": "housekeeping", "GAPDH": "housekeeping", "B2M": "housekeeping",
    "RPL13A": "housekeeping", "YWHAZ": "housekeeping",
}

def parse_control_results():
    """Parse pickle files from control perturbation runs."""
    results = {}
    raw_shifts = {}  # Store per-cell shifts for statistical tests

    for d in sorted(os.listdir(BASE)):
        if not d.startswith("control_perturb_"):
            continue
        gene = d.replace("control_perturb_", "")
        path = os.path.join(BASE, d)
        files = [f for f in os.listdir(path) if f.endswith('.pickle')]
        if not files:
            print(f"  SKIP {gene} - no pickle files")
            continue

        pkl_file = os.path.join(path, files[0])
        with open(pkl_file, "rb") as f:
            data = pickle.load(f)

        for key, val in data.items():
            token_id = key[0] if isinstance(key, tuple) else key
            emb_type = key[1] if isinstance(key, tuple) else "unknown"

            if isinstance(val, list):
                arr = np.array(val)
                # These are cosine similarity values (1 - shift)
                cos_sim_mean = float(np.mean(arr))
                cos_sim_std = float(np.std(arr))
                shift_mean = 1.0 - cos_sim_mean
                n_cells = len(arr)
                results[gene] = {
                    "cosine_sim_mean": cos_sim_mean,
                    "cosine_sim_std": cos_sim_std,
                    "shift_mean": shift_mean,
                    "shift_std": cos_sim_std,
                    "n_cells": n_cells,
                    "matched_to": MATCHED_TO.get(gene, "unknown"),
                    "type": "housekeeping" if MATCHED_TO.get(gene) == "housekeeping" else "expression_matched",
                }
                raw_shifts[gene] = 1.0 - arr  # Per-cell shift values
                print(f"  {gene}: shift={shift_mean:.6f} +/- {cos_sim_std:.6f}, n_cells={n_cells}, matched_to={MATCHED_TO.get(gene)}")

    return results, raw_shifts


def parse_intelligence_results():
    """Parse existing intelligence gene results."""
    with open(INTEL_RESULTS_FILE) as f:
        intel_data = json.load(f)

    results = {}
    for gene, info in intel_data.items():
        results[gene] = {
            "cosine_sim_mean": info["cosine_sim_mean"],
            "cosine_sim_std": info["cosine_sim_std"],
            "shift_mean": info["shift_mean"],
            "shift_std": info["shift_std"],
            "n_cells": info["n_cells"],
        }

    # Also load per-cell raw data from pickle files for statistical tests
    raw_shifts = {}
    for d in sorted(os.listdir(BASE)):
        if d.startswith("perturb_") and not d.startswith("perturb_combo_"):
            gene = d.replace("perturb_", "")
            if gene in ("combo", ):
                continue
            path = os.path.join(BASE, d)
            files = [f for f in os.listdir(path) if f.endswith('.pickle')]
            if not files:
                continue
            pkl_file = os.path.join(path, files[0])
            with open(pkl_file, "rb") as f:
                data = pickle.load(f)
            for key, val in data.items():
                if isinstance(val, list):
                    arr = np.array(val)
                    raw_shifts[gene] = 1.0 - arr  # cosine sim -> shift

    return results, raw_shifts


def cohens_d(group1, group2):
    """Compute Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return float('nan')
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return float('nan')
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def run_comparison(intel_results, intel_raw, control_results, control_raw):
    """Run statistical comparisons between intelligence and control genes."""
    print("\n" + "=" * 70)
    print("STATISTICAL COMPARISON: INTELLIGENCE vs CONTROL GENES")
    print("=" * 70)

    # Separate expression-matched and housekeeping controls
    expr_matched = {g: r for g, r in control_results.items() if r["type"] == "expression_matched"}
    housekeeping = {g: r for g, r in control_results.items() if r["type"] == "housekeeping"}

    # Get mean shifts for each group
    intel_shifts = [r["shift_mean"] for r in intel_results.values()]
    expr_ctrl_shifts = [r["shift_mean"] for r in expr_matched.values()]
    house_shifts = [r["shift_mean"] for r in housekeeping.values()]
    all_ctrl_shifts = [r["shift_mean"] for r in control_results.values()]

    # ---- Group-level statistics (mean shifts) ----
    print("\n--- Group-Level Statistics (Gene-Mean Shifts) ---")
    print(f"Intelligence genes (n={len(intel_shifts)}): mean={np.mean(intel_shifts):.6f}, median={np.median(intel_shifts):.6f}, std={np.std(intel_shifts):.6f}")
    print(f"Expression-matched controls (n={len(expr_ctrl_shifts)}): mean={np.mean(expr_ctrl_shifts):.6f}, median={np.median(expr_ctrl_shifts):.6f}, std={np.std(expr_ctrl_shifts):.6f}")
    print(f"Housekeeping controls (n={len(house_shifts)}): mean={np.mean(house_shifts):.6f}, median={np.median(house_shifts):.6f}, std={np.std(house_shifts):.6f}")
    print(f"All controls (n={len(all_ctrl_shifts)}): mean={np.mean(all_ctrl_shifts):.6f}, median={np.median(all_ctrl_shifts):.6f}, std={np.std(all_ctrl_shifts):.6f}")

    comparisons = {}

    # ---- Mann-Whitney U tests (gene-level) ----
    print("\n--- Mann-Whitney U Tests (Gene-Level) ---")

    # Intelligence vs Expression-matched
    if len(expr_ctrl_shifts) >= 2:
        u_stat, p_val = stats.mannwhitneyu(intel_shifts, expr_ctrl_shifts, alternative='greater')
        d = cohens_d(intel_shifts, expr_ctrl_shifts)
        print(f"\nIntelligence vs Expression-Matched Controls:")
        print(f"  U = {u_stat:.1f}, p = {p_val:.6f} (one-sided, greater)")
        print(f"  Cohen's d = {d:.4f}")
        print(f"  Significant (p < 0.05): {'YES' if p_val < 0.05 else 'NO'}")
        comparisons["intel_vs_expr_matched"] = {
            "test": "Mann-Whitney U (one-sided, greater)",
            "U": float(u_stat), "p_value": float(p_val),
            "cohens_d": float(d),
            "intel_mean": float(np.mean(intel_shifts)),
            "control_mean": float(np.mean(expr_ctrl_shifts)),
            "significant": bool(p_val < 0.05),
        }

    # Intelligence vs Housekeeping
    if len(house_shifts) >= 2:
        u_stat, p_val = stats.mannwhitneyu(intel_shifts, house_shifts, alternative='greater')
        d = cohens_d(intel_shifts, house_shifts)
        print(f"\nIntelligence vs Housekeeping Controls:")
        print(f"  U = {u_stat:.1f}, p = {p_val:.6f} (one-sided, greater)")
        print(f"  Cohen's d = {d:.4f}")
        print(f"  Significant (p < 0.05): {'YES' if p_val < 0.05 else 'NO'}")
        comparisons["intel_vs_housekeeping"] = {
            "test": "Mann-Whitney U (one-sided, greater)",
            "U": float(u_stat), "p_value": float(p_val),
            "cohens_d": float(d),
            "intel_mean": float(np.mean(intel_shifts)),
            "control_mean": float(np.mean(house_shifts)),
            "significant": bool(p_val < 0.05),
        }

    # Intelligence vs All Controls
    if len(all_ctrl_shifts) >= 2:
        u_stat, p_val = stats.mannwhitneyu(intel_shifts, all_ctrl_shifts, alternative='greater')
        d = cohens_d(intel_shifts, all_ctrl_shifts)
        print(f"\nIntelligence vs All Controls:")
        print(f"  U = {u_stat:.1f}, p = {p_val:.6f} (one-sided, greater)")
        print(f"  Cohen's d = {d:.4f}")
        print(f"  Significant (p < 0.05): {'YES' if p_val < 0.05 else 'NO'}")
        comparisons["intel_vs_all_controls"] = {
            "test": "Mann-Whitney U (one-sided, greater)",
            "U": float(u_stat), "p_value": float(p_val),
            "cohens_d": float(d),
            "intel_mean": float(np.mean(intel_shifts)),
            "control_mean": float(np.mean(all_ctrl_shifts)),
            "significant": bool(p_val < 0.05),
        }

    # ---- Paired comparison (matched controls) ----
    print("\n--- Paired Comparison (Expression-Matched) ---")
    paired_intel = []
    paired_ctrl = []
    pair_names = []
    for ctrl_gene, ctrl_info in expr_matched.items():
        intel_gene = ctrl_info["matched_to"]
        if intel_gene in intel_results:
            paired_intel.append(intel_results[intel_gene]["shift_mean"])
            paired_ctrl.append(ctrl_info["shift_mean"])
            pair_names.append(f"{intel_gene} vs {ctrl_gene}")
            print(f"  {intel_gene}: {intel_results[intel_gene]['shift_mean']:.6f}  vs  {ctrl_gene}: {ctrl_info['shift_mean']:.6f}  (diff: {intel_results[intel_gene]['shift_mean'] - ctrl_info['shift_mean']:.6f})")

    if len(paired_intel) >= 2:
        w_stat, p_val = stats.wilcoxon(paired_intel, paired_ctrl, alternative='greater')
        differences = [a - b for a, b in zip(paired_intel, paired_ctrl)]
        mean_diff = np.mean(differences)
        n_intel_higher = sum(1 for d in differences if d > 0)
        print(f"\n  Wilcoxon signed-rank test:")
        print(f"  W = {w_stat:.1f}, p = {p_val:.6f} (one-sided, greater)")
        print(f"  Mean difference: {mean_diff:.6f}")
        print(f"  Intelligence higher in {n_intel_higher}/{len(differences)} pairs")
        comparisons["paired_wilcoxon"] = {
            "test": "Wilcoxon signed-rank (one-sided, greater)",
            "W": float(w_stat), "p_value": float(p_val),
            "mean_difference": float(mean_diff),
            "n_pairs": len(paired_intel),
            "n_intel_higher": int(n_intel_higher),
            "significant": bool(p_val < 0.05),
        }

    # ---- Cell-level pooled comparison ----
    print("\n--- Cell-Level Pooled Comparison ---")
    all_intel_cells = np.concatenate([v for v in intel_raw.values()])
    all_ctrl_cells = np.concatenate([v for k, v in control_raw.items() if MATCHED_TO.get(k) != "housekeeping"])
    all_house_cells = np.concatenate([v for k, v in control_raw.items() if MATCHED_TO.get(k) == "housekeeping"]) if any(MATCHED_TO.get(k) == "housekeeping" for k in control_raw) else np.array([])

    print(f"  Intelligence cells: n={len(all_intel_cells)}, mean={np.mean(all_intel_cells):.6f}")
    print(f"  Expr-matched control cells: n={len(all_ctrl_cells)}, mean={np.mean(all_ctrl_cells):.6f}")
    if len(all_house_cells) > 0:
        print(f"  Housekeeping control cells: n={len(all_house_cells)}, mean={np.mean(all_house_cells):.6f}")

    u_stat, p_val = stats.mannwhitneyu(all_intel_cells, all_ctrl_cells, alternative='greater')
    d = cohens_d(all_intel_cells, all_ctrl_cells)
    print(f"\n  Intelligence vs Expr-Matched (cell-level):")
    print(f"  U = {u_stat:.0f}, p = {p_val:.2e} (one-sided, greater)")
    print(f"  Cohen's d = {d:.4f}")
    comparisons["cell_level_intel_vs_expr_matched"] = {
        "test": "Mann-Whitney U cell-level (one-sided, greater)",
        "U": float(u_stat), "p_value": float(p_val),
        "cohens_d": float(d),
        "n_intel_cells": int(len(all_intel_cells)),
        "n_control_cells": int(len(all_ctrl_cells)),
        "significant": bool(p_val < 0.05),
    }

    return comparisons


def main():
    print("=" * 70)
    print("CONTROL GENE PERTURBATION RESULTS ANALYSIS")
    print("=" * 70)

    # Parse control results
    print("\n--- Parsing Control Gene Results ---")
    control_results, control_raw = parse_control_results()
    print(f"\nTotal control genes analyzed: {len(control_results)}")

    if len(control_results) == 0:
        print("ERROR: No control results found!")
        return

    # Parse intelligence gene results
    print("\n--- Parsing Intelligence Gene Results ---")
    intel_results, intel_raw = parse_intelligence_results()
    print(f"Total intelligence genes: {len(intel_results)}")

    # Run comparisons
    comparisons = run_comparison(intel_results, intel_raw, control_results, control_raw)

    # Save all results
    output = {
        "control_results": control_results,
        "comparisons": comparisons,
        "metadata": {
            "n_intelligence_genes": len(intel_results),
            "n_control_genes": len(control_results),
            "n_expression_matched": sum(1 for r in control_results.values() if r["type"] == "expression_matched"),
            "n_housekeeping": sum(1 for r in control_results.values() if r["type"] == "housekeeping"),
        }
    }

    out_file = "/mnt/d/openclaw/intelligence-augmentation/analysis/control_perturbation_results.json"
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
