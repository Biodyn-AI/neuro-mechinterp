exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
#!/usr/bin/env python3
"""
Compare intelligence gene perturbation effects vs control gene effects.
Uses Mann-Whitney U test and Cohen's d for statistical comparison.
"""
import os, sys, json, pickle, glob
import numpy as np
from scipy import stats
from collections import defaultdict

OUTPUT_BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"

# Intelligence gene map
INTELLIGENCE_GENES = {
    "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
    "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
    "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
    "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
    "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
    "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
    "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
}

# Control gene categories
EXPRESSION_MATCHED = {
    "EXOC6B": ("ENSG00000144036", "CADM2"),
    "DOCK9": ("ENSG00000088387", "APP"),
    "PDE7B": ("ENSG00000171408", "FOXO3"),
    "ZDHHC9": ("ENSG00000188706", "NLGN1"),
    "HDAC11": ("ENSG00000163517", "MAPT"),
    "XRCC6": ("ENSG00000196419", "TCF4"),
    "WDR46": ("ENSG00000227057", "CREB1"),
    "CUEDC2": ("ENSG00000107874", "DISC1"),
    "IGSF10": ("ENSG00000152580", "NEGR1"),
    "SALL3": ("ENSG00000256463", "FMR1"),
    "MFSD6L": ("ENSG00000185156", "COMT"),
    "TMEM88": ("ENSG00000167874", "MEF2C"),
}

HOUSEKEEPING = {
    "ACTB": "ENSG00000075624",
    "GAPDH": "ENSG00000111640",
    "B2M": "ENSG00000166710",
    "RPL13A": "ENSG00000142541",
    "YWHAZ": "ENSG00000164924",
}

def load_perturbation_result(result_dir):
    """Load perturbation result from a directory. Returns cosine similarities list."""
    if not os.path.exists(result_dir):
        return None
    pickles = [f for f in os.listdir(result_dir) if f.endswith('.pickle')]
    if not pickles:
        return None
    fpath = os.path.join(result_dir, pickles[0])
    with open(fpath, 'rb') as f:
        data = pickle.load(f)
    if len(data) == 0:
        return None
    key = list(data.keys())[0]
    cos_sims = data[key]
    return cos_sims

def compute_shift(cos_sims):
    """Convert cosine similarities to embedding shift (1 - cos_sim)."""
    return [1.0 - c for c in cos_sims]

def cohens_d(group1, group2):
    """Compute Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (mean1 - mean2) / pooled_std

def main():
    print("=" * 70)
    print("INTELLIGENCE vs CONTROL GENE PERTURBATION ANALYSIS")
    print("=" * 70)

    # Load intelligence gene results
    print("\n--- Loading Intelligence Gene Results ---")
    intel_results = {}
    for symbol in INTELLIGENCE_GENES:
        result_dir = os.path.join(OUTPUT_BASE, f"perturb_{symbol}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = compute_shift(cos_sims)
            intel_results[symbol] = {
                "cos_sims": cos_sims,
                "shifts": shifts,
                "n_cells": len(cos_sims),
                "mean_shift": np.mean(shifts),
                "std_shift": np.std(shifts),
                "mean_cos_sim": np.mean(cos_sims),
            }
            print(f"  {symbol}: n={len(cos_sims)}, mean_shift={np.mean(shifts):.6f}")
        else:
            print(f"  {symbol}: NO RESULT")

    # Load expression-matched control results
    print("\n--- Loading Expression-Matched Control Results ---")
    expr_ctrl_results = {}
    for symbol, (ensembl, matched_to) in EXPRESSION_MATCHED.items():
        result_dir = os.path.join(OUTPUT_BASE, f"control_perturb_{symbol}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = compute_shift(cos_sims)
            expr_ctrl_results[symbol] = {
                "cos_sims": cos_sims,
                "shifts": shifts,
                "n_cells": len(cos_sims),
                "mean_shift": np.mean(shifts),
                "std_shift": np.std(shifts),
                "mean_cos_sim": np.mean(cos_sims),
                "matched_to": matched_to,
            }
            print(f"  {symbol} (ctrl for {matched_to}): n={len(cos_sims)}, mean_shift={np.mean(shifts):.6f}")
        else:
            print(f"  {symbol} (ctrl for {matched_to}): NO RESULT")

    # Load housekeeping control results
    print("\n--- Loading Housekeeping Control Results ---")
    hk_results = {}
    for symbol in HOUSEKEEPING:
        result_dir = os.path.join(OUTPUT_BASE, f"control_perturb_{symbol}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = compute_shift(cos_sims)
            hk_results[symbol] = {
                "cos_sims": cos_sims,
                "shifts": shifts,
                "n_cells": len(cos_sims),
                "mean_shift": np.mean(shifts),
                "std_shift": np.std(shifts),
                "mean_cos_sim": np.mean(cos_sims),
            }
            print(f"  {symbol}: n={len(cos_sims)}, mean_shift={np.mean(shifts):.6f}")
        else:
            print(f"  {symbol}: NO RESULT")

    # Combine all controls
    all_ctrl_results = {**expr_ctrl_results, **hk_results}

    # ============================================================
    # ANALYSIS 1: Gene-level mean shift comparison
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 1: Gene-level Mean Embedding Shift Comparison")
    print("=" * 70)

    intel_mean_shifts = [r["mean_shift"] for r in intel_results.values()]
    expr_ctrl_mean_shifts = [r["mean_shift"] for r in expr_ctrl_results.values()]
    hk_mean_shifts = [r["mean_shift"] for r in hk_results.values()]
    all_ctrl_mean_shifts = [r["mean_shift"] for r in all_ctrl_results.values()]

    print(f"\nIntelligence genes (n={len(intel_mean_shifts)}):")
    print(f"  Mean of mean shifts: {np.mean(intel_mean_shifts):.6f}")
    print(f"  Median of mean shifts: {np.median(intel_mean_shifts):.6f}")
    print(f"  Std of mean shifts: {np.std(intel_mean_shifts):.6f}")
    print(f"  Range: [{np.min(intel_mean_shifts):.6f}, {np.max(intel_mean_shifts):.6f}]")

    print(f"\nExpression-matched controls (n={len(expr_ctrl_mean_shifts)}):")
    print(f"  Mean of mean shifts: {np.mean(expr_ctrl_mean_shifts):.6f}")
    print(f"  Median of mean shifts: {np.median(expr_ctrl_mean_shifts):.6f}")
    print(f"  Std of mean shifts: {np.std(expr_ctrl_mean_shifts):.6f}")
    print(f"  Range: [{np.min(expr_ctrl_mean_shifts):.6f}, {np.max(expr_ctrl_mean_shifts):.6f}]")

    print(f"\nHousekeeping controls (n={len(hk_mean_shifts)}):")
    print(f"  Mean of mean shifts: {np.mean(hk_mean_shifts):.6f}")
    print(f"  Median of mean shifts: {np.median(hk_mean_shifts):.6f}")
    print(f"  Std of mean shifts: {np.std(hk_mean_shifts):.6f}")
    print(f"  Range: [{np.min(hk_mean_shifts):.6f}, {np.max(hk_mean_shifts):.6f}]")

    print(f"\nAll controls (n={len(all_ctrl_mean_shifts)}):")
    print(f"  Mean of mean shifts: {np.mean(all_ctrl_mean_shifts):.6f}")
    print(f"  Median of mean shifts: {np.median(all_ctrl_mean_shifts):.6f}")

    # Mann-Whitney U tests
    print("\n--- Mann-Whitney U Tests (gene-level mean shifts) ---")

    # Intelligence vs expression-matched controls
    u_stat, p_val = stats.mannwhitneyu(intel_mean_shifts, expr_ctrl_mean_shifts, alternative='greater')
    d_val = cohens_d(intel_mean_shifts, expr_ctrl_mean_shifts)
    print(f"\nIntelligence vs Expression-matched controls:")
    print(f"  U statistic: {u_stat:.1f}")
    print(f"  p-value (one-tailed, greater): {p_val:.6f}")
    print(f"  Cohen's d: {d_val:.4f}")
    print(f"  Interpretation: {'Significant' if p_val < 0.05 else 'Not significant'} at p < 0.05")

    intel_vs_expr = {"u_stat": float(u_stat), "p_value": float(p_val), "cohens_d": float(d_val)}

    # Intelligence vs housekeeping controls
    u_stat2, p_val2 = stats.mannwhitneyu(intel_mean_shifts, hk_mean_shifts, alternative='greater')
    d_val2 = cohens_d(intel_mean_shifts, hk_mean_shifts)
    print(f"\nIntelligence vs Housekeeping controls:")
    print(f"  U statistic: {u_stat2:.1f}")
    print(f"  p-value (one-tailed, greater): {p_val2:.6f}")
    print(f"  Cohen's d: {d_val2:.4f}")
    print(f"  Interpretation: {'Significant' if p_val2 < 0.05 else 'Not significant'} at p < 0.05")

    intel_vs_hk = {"u_stat": float(u_stat2), "p_value": float(p_val2), "cohens_d": float(d_val2)}

    # Intelligence vs all controls
    u_stat3, p_val3 = stats.mannwhitneyu(intel_mean_shifts, all_ctrl_mean_shifts, alternative='greater')
    d_val3 = cohens_d(intel_mean_shifts, all_ctrl_mean_shifts)
    print(f"\nIntelligence vs All controls:")
    print(f"  U statistic: {u_stat3:.1f}")
    print(f"  p-value (one-tailed, greater): {p_val3:.6f}")
    print(f"  Cohen's d: {d_val3:.4f}")
    print(f"  Interpretation: {'Significant' if p_val3 < 0.05 else 'Not significant'} at p < 0.05")

    intel_vs_all = {"u_stat": float(u_stat3), "p_value": float(p_val3), "cohens_d": float(d_val3)}

    # ============================================================
    # ANALYSIS 2: Cell-level pooled comparison
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 2: Cell-level Pooled Embedding Shift Comparison")
    print("=" * 70)

    # Pool all cell-level shifts
    intel_all_shifts = []
    for r in intel_results.values():
        intel_all_shifts.extend(r["shifts"])

    ctrl_all_shifts = []
    for r in all_ctrl_results.values():
        ctrl_all_shifts.extend(r["shifts"])

    expr_ctrl_all_shifts = []
    for r in expr_ctrl_results.values():
        expr_ctrl_all_shifts.extend(r["shifts"])

    hk_all_shifts = []
    for r in hk_results.values():
        hk_all_shifts.extend(r["shifts"])

    print(f"\nIntelligence: {len(intel_all_shifts)} total cell observations")
    print(f"  Mean shift: {np.mean(intel_all_shifts):.6f}")
    print(f"  Median shift: {np.median(intel_all_shifts):.6f}")

    print(f"\nExpression-matched controls: {len(expr_ctrl_all_shifts)} total cell observations")
    print(f"  Mean shift: {np.mean(expr_ctrl_all_shifts):.6f}")
    print(f"  Median shift: {np.median(expr_ctrl_all_shifts):.6f}")

    print(f"\nHousekeeping controls: {len(hk_all_shifts)} total cell observations")
    print(f"  Mean shift: {np.mean(hk_all_shifts):.6f}")
    print(f"  Median shift: {np.median(hk_all_shifts):.6f}")

    # Cell-level Mann-Whitney
    u4, p4 = stats.mannwhitneyu(intel_all_shifts, expr_ctrl_all_shifts, alternative='greater')
    d4 = cohens_d(intel_all_shifts, expr_ctrl_all_shifts)
    print(f"\nCell-level: Intelligence vs Expression-matched controls:")
    print(f"  U: {u4:.0f}, p: {p4:.2e}, Cohen's d: {d4:.4f}")

    u5, p5 = stats.mannwhitneyu(intel_all_shifts, ctrl_all_shifts, alternative='greater')
    d5 = cohens_d(intel_all_shifts, ctrl_all_shifts)
    print(f"Cell-level: Intelligence vs All controls:")
    print(f"  U: {u5:.0f}, p: {p5:.2e}, Cohen's d: {d5:.4f}")

    # ============================================================
    # ANALYSIS 3: Paired comparison (expression-matched only)
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 3: Paired Comparison (Expression-Matched Controls)")
    print("=" * 70)

    paired_data = []
    for ctrl_sym, (_, matched_intel) in EXPRESSION_MATCHED.items():
        if ctrl_sym in expr_ctrl_results and matched_intel in intel_results:
            intel_shift = intel_results[matched_intel]["mean_shift"]
            ctrl_shift = expr_ctrl_results[ctrl_sym]["mean_shift"]
            diff = intel_shift - ctrl_shift
            paired_data.append({
                "intelligence_gene": matched_intel,
                "control_gene": ctrl_sym,
                "intel_shift": intel_shift,
                "ctrl_shift": ctrl_shift,
                "difference": diff,
                "intel_n": intel_results[matched_intel]["n_cells"],
                "ctrl_n": expr_ctrl_results[ctrl_sym]["n_cells"],
            })
            print(f"  {matched_intel} ({intel_shift:.6f}) vs {ctrl_sym} ({ctrl_shift:.6f}): diff={diff:+.6f}")

    if paired_data:
        diffs = [p["difference"] for p in paired_data]
        intel_wins = sum(1 for d in diffs if d > 0)
        print(f"\nIntelligence gene > control: {intel_wins}/{len(diffs)} pairs")
        print(f"Mean difference: {np.mean(diffs):.6f}")
        print(f"Median difference: {np.median(diffs):.6f}")

        # Wilcoxon signed-rank test on paired differences
        if len(diffs) >= 6:
            w_stat, w_pval = stats.wilcoxon(diffs, alternative='greater')
            print(f"\nWilcoxon signed-rank test (paired, one-tailed):")
            print(f"  W statistic: {w_stat:.1f}")
            print(f"  p-value: {w_pval:.6f}")
            print(f"  Interpretation: {'Significant' if w_pval < 0.05 else 'Not significant'} at p < 0.05")
        else:
            w_stat, w_pval = None, None
            print(f"\nToo few pairs ({len(diffs)}) for Wilcoxon test")

        # Sign test
        p_sign = stats.binomtest(intel_wins, len(diffs), 0.5, alternative='greater').pvalue
        print(f"\nSign test: {intel_wins}/{len(diffs)} intelligence > control")
        print(f"  p-value: {p_sign:.6f}")

    # ============================================================
    # ANALYSIS 4: Effect size interpretation
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYSIS 4: Effect Size Interpretation")
    print("=" * 70)

    def interpret_d(d):
        d = abs(d)
        if d < 0.2: return "negligible"
        elif d < 0.5: return "small"
        elif d < 0.8: return "medium"
        else: return "large"

    print(f"\nCohen's d interpretations:")
    print(f"  Intelligence vs Expr-matched: d={intel_vs_expr['cohens_d']:.4f} ({interpret_d(intel_vs_expr['cohens_d'])})")
    print(f"  Intelligence vs Housekeeping: d={intel_vs_hk['cohens_d']:.4f} ({interpret_d(intel_vs_hk['cohens_d'])})")
    print(f"  Intelligence vs All controls: d={intel_vs_all['cohens_d']:.4f} ({interpret_d(intel_vs_all['cohens_d'])})")

    # ============================================================
    # FULL RESULTS TABLE
    # ============================================================
    print("\n" + "=" * 70)
    print("COMPLETE RESULTS TABLE")
    print("=" * 70)

    print(f"\n{'Gene':>12} {'Type':>15} {'N cells':>8} {'Mean Shift':>12} {'Std':>10}")
    print("-" * 60)
    for sym in sorted(intel_results.keys(), key=lambda s: intel_results[s]["mean_shift"], reverse=True):
        r = intel_results[sym]
        print(f"{sym:>12} {'intelligence':>15} {r['n_cells']:>8} {r['mean_shift']:>12.6f} {r['std_shift']:>10.6f}")
    print("-" * 60)
    for sym in sorted(expr_ctrl_results.keys(), key=lambda s: expr_ctrl_results[s]["mean_shift"], reverse=True):
        r = expr_ctrl_results[sym]
        print(f"{sym:>12} {'expr-ctrl':>15} {r['n_cells']:>8} {r['mean_shift']:>12.6f} {r['std_shift']:>10.6f}")
    print("-" * 60)
    for sym in sorted(hk_results.keys(), key=lambda s: hk_results[s]["mean_shift"], reverse=True):
        r = hk_results[sym]
        print(f"{sym:>12} {'housekeeping':>15} {r['n_cells']:>8} {r['mean_shift']:>12.6f} {r['std_shift']:>10.6f}")

    # ============================================================
    # Save all results as JSON
    # ============================================================
    output_json = {
        "intelligence_genes": {k: {"mean_shift": v["mean_shift"], "std_shift": v["std_shift"],
                                    "n_cells": v["n_cells"], "mean_cos_sim": v["mean_cos_sim"]}
                               for k, v in intel_results.items()},
        "expression_matched_controls": {k: {"mean_shift": v["mean_shift"], "std_shift": v["std_shift"],
                                             "n_cells": v["n_cells"], "mean_cos_sim": v["mean_cos_sim"],
                                             "matched_to": v["matched_to"]}
                                        for k, v in expr_ctrl_results.items()},
        "housekeeping_controls": {k: {"mean_shift": v["mean_shift"], "std_shift": v["std_shift"],
                                       "n_cells": v["n_cells"], "mean_cos_sim": v["mean_cos_sim"]}
                                  for k, v in hk_results.items()},
        "statistical_tests": {
            "gene_level": {
                "intel_vs_expr_matched": intel_vs_expr,
                "intel_vs_housekeeping": intel_vs_hk,
                "intel_vs_all_controls": intel_vs_all,
            },
            "cell_level": {
                "intel_vs_expr_matched": {"u_stat": float(u4), "p_value": float(p4), "cohens_d": float(d4)},
                "intel_vs_all_controls": {"u_stat": float(u5), "p_value": float(p5), "cohens_d": float(d5)},
            },
            "paired_comparison": {
                "pairs": paired_data,
                "intel_wins": intel_wins if paired_data else 0,
                "total_pairs": len(paired_data),
                "mean_difference": float(np.mean(diffs)) if paired_data else 0,
                "wilcoxon_p": float(w_pval) if (paired_data and w_pval is not None) else None,
                "sign_test_p": float(p_sign) if paired_data else None,
            }
        },
        "summary_stats": {
            "intelligence": {
                "n_genes": len(intel_results),
                "mean_of_means": float(np.mean(intel_mean_shifts)),
                "median_of_means": float(np.median(intel_mean_shifts)),
                "total_cells": len(intel_all_shifts),
            },
            "expr_matched_controls": {
                "n_genes": len(expr_ctrl_results),
                "mean_of_means": float(np.mean(expr_ctrl_mean_shifts)),
                "median_of_means": float(np.median(expr_ctrl_mean_shifts)),
                "total_cells": len(expr_ctrl_all_shifts),
            },
            "housekeeping_controls": {
                "n_genes": len(hk_results),
                "mean_of_means": float(np.mean(hk_mean_shifts)),
                "median_of_means": float(np.median(hk_mean_shifts)),
                "total_cells": len(hk_all_shifts),
            },
        }
    }

    json_path = os.path.join(OUTPUT_BASE, "control_vs_intelligence_analysis.json")
    with open(json_path, "w") as f:
        json.dump(output_json, f, indent=2)
    print(f"\nResults saved to {json_path}")

    # Return for use by report generator
    return output_json

if __name__ == "__main__":
    results = main()
