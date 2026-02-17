exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
"""
Step 3: Analyze permutation null results.
Compute empirical p-values for each intelligence gene vs its expression-bin null.
Apply FDR correction. Generate PERMUTATION_NULL_FINAL.md.
"""
import os, sys, json, pickle, time
import numpy as np
from collections import defaultdict

OUTPUT_BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
ANALYSIS_DIR = "/mnt/d/openclaw/intelligence-augmentation/analysis"
SELECTION_JSON = os.path.join(ANALYSIS_DIR, "permutation_null_v2_selection.json")
INCREMENTAL_JSON = os.path.join(ANALYSIS_DIR, "permutation_null_v2_incremental.json")

INTELLIGENCE_GENES = {
    "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
    "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
    "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
    "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
    "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
    "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
    "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
}

EXISTING_CONTROLS = {
    "EXOC6B": "ENSG00000144036", "DOCK9": "ENSG00000088387", "PDE7B": "ENSG00000171408",
    "ZDHHC9": "ENSG00000188706", "HDAC11": "ENSG00000163517", "XRCC6": "ENSG00000196419",
    "WDR46": "ENSG00000227057", "CUEDC2": "ENSG00000107874", "IGSF10": "ENSG00000152580",
    "SALL3": "ENSG00000256463", "TMEM88": "ENSG00000167874",
}


def load_perturbation_result(result_dir):
    """Load perturbation result pickle. Returns list of cosine similarities."""
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
    return data[key]


def main():
    print("=== PERMUTATION NULL V2: ANALYSIS ===\n")

    # Load selection info
    with open(SELECTION_JSON) as f:
        selection = json.load(f)
    bin_info = selection['bin_info']
    intel_info = selection['intelligence_genes']
    gene_freq = selection.get('gene_frequencies', {})

    # Load incremental control results
    with open(INCREMENTAL_JSON) as f:
        control_results = json.load(f)

    # Count successes
    successful = {k: v for k, v in control_results.items()
                  if v.get('status') in ('success', 'reused')}
    print(f"Control genes: {len(control_results)} total, {len(successful)} successful\n")

    # Load intelligence gene results from pickle files
    print("Loading intelligence gene perturbation results...")
    intel_results = {}
    for symbol, ensembl in INTELLIGENCE_GENES.items():
        result_dir = os.path.join(OUTPUT_BASE, f"perturb_{symbol}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = [1.0 - c for c in cos_sims]
            freq_val = intel_info.get(symbol, {}).get('frequency')
            intel_results[symbol] = {
                'mean_shift': float(np.mean(shifts)),
                'std_shift': float(np.std(shifts)),
                'n_cells': len(shifts),
                'frequency': freq_val,
            }
            print(f"  {symbol}: shift={np.mean(shifts):.6f}, n={len(shifts)}")
        else:
            print(f"  WARNING: {symbol} - no perturbation result found")

    # Also load existing control results from previous experiment (for extra null power)
    print("\nLoading existing expression-matched control results...")
    existing_ctrl_shifts = {}
    for ctrl_sym, ctrl_ens in EXISTING_CONTROLS.items():
        result_dir = os.path.join(OUTPUT_BASE, f"control_perturb_{ctrl_sym}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = [1.0 - c for c in cos_sims]
            freq_val = gene_freq.get(ctrl_ens, {}).get('frequency')
            existing_ctrl_shifts[ctrl_sym] = {
                'mean_shift': float(np.mean(shifts)),
                'n_cells': len(shifts),
                'frequency': freq_val,
                'ensembl': ctrl_ens,
            }
            print(f"  {ctrl_sym}: shift={np.mean(shifts):.6f}, n={len(shifts)}")

    # Build null distribution per expression bin
    print("\nBuilding null distributions per expression bin...")
    bin_edges = bin_info['bin_edges']
    n_bins = len(bin_edges) - 1

    null_per_bin = defaultdict(list)  # bin -> list of mean_shifts

    # Add new permutation controls
    for sym, res in successful.items():
        b = res.get('bin')
        if b is not None and 'mean_shift' in res:
            null_per_bin[b].append(res['mean_shift'])

    # Add existing controls, assigned to bins by frequency
    for ctrl_sym, info in existing_ctrl_shifts.items():
        freq = info.get('frequency')
        if freq is not None:
            for b in range(n_bins):
                if bin_edges[b] <= freq <= bin_edges[b + 1]:
                    null_per_bin[b].append(info['mean_shift'])
                    break

    for b in sorted(null_per_bin.keys()):
        shifts = null_per_bin[b]
        print(f"  Bin {b}: {len(shifts)} controls, mean={np.mean(shifts):.6f}, std={np.std(shifts):.6f}")

    # Compute empirical p-values
    print("\nComputing empirical p-values...")
    pvalues = {}
    intel_bins = bin_info.get('intel_bins', {})

    for symbol, ires in intel_results.items():
        intel_shift = ires['mean_shift']
        intel_freq = ires.get('frequency')
        intel_bin = intel_bins.get(symbol, {}).get('bin')

        # Strategy: use same-bin null. If < 5 controls, expand to adjacent bins.
        use_null = []
        null_type = "bin_based"

        if intel_bin is not None:
            use_null = list(null_per_bin.get(intel_bin, []))

        # If insufficient, expand to adjacent bins
        if len(use_null) < 5 and intel_bin is not None:
            for delta in [1, -1, 2, -2]:
                adj_bin = intel_bin + delta
                if 0 <= adj_bin < n_bins:
                    use_null.extend(null_per_bin.get(adj_bin, []))
                if len(use_null) >= 5:
                    break
            if len(use_null) > len(null_per_bin.get(intel_bin, [])):
                null_type = "adjacent_bins"

        # Last resort: use all controls
        if len(use_null) < 3:
            use_null = []
            for b in range(n_bins):
                use_null.extend(null_per_bin.get(b, []))
            null_type = "all_controls"

        n_null = len(use_null)
        n_greater = sum(1 for s in use_null if s >= intel_shift)
        empirical_p = (n_greater + 1) / (n_null + 1)

        # Z-score
        if n_null > 1:
            null_mean = float(np.mean(use_null))
            null_std = float(np.std(use_null))
            z_score = (intel_shift - null_mean) / null_std if null_std > 0 else None
        else:
            null_mean = None
            null_std = None
            z_score = None

        pvalues[symbol] = {
            'mean_shift': intel_shift,
            'n_cells': ires['n_cells'],
            'expression_frequency': intel_freq,
            'null_type': null_type,
            'n_null_genes': n_null,
            'n_exceeding': n_greater,
            'empirical_p': float(empirical_p),
            'z_score': float(z_score) if z_score is not None else None,
            'null_mean': null_mean,
            'null_std': null_std,
            'bin': intel_bin,
        }

    # FDR correction (Benjamini-Hochberg)
    print("\nApplying Benjamini-Hochberg FDR correction...")
    genes = sorted(pvalues.keys())
    raw_ps = [pvalues[g]['empirical_p'] for g in genes]
    n = len(raw_ps)
    sorted_indices = np.argsort(raw_ps)
    sorted_ps = np.array(raw_ps)[sorted_indices]

    fdr_values = np.zeros(n)
    for i in range(n):
        rank = i + 1
        fdr_values[i] = sorted_ps[i] * n / rank
    for i in range(n - 2, -1, -1):
        fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])
    fdr_values = np.minimum(fdr_values, 1.0)

    fdr_adjusted = np.zeros(n)
    for i in range(n):
        fdr_adjusted[sorted_indices[i]] = fdr_values[i]

    for i, g in enumerate(genes):
        pvalues[g]['fdr_q'] = float(fdr_adjusted[i])
        sig = "***" if fdr_adjusted[i] < 0.001 else "**" if fdr_adjusted[i] < 0.01 else "*" if fdr_adjusted[i] < 0.05 else "+" if fdr_adjusted[i] < 0.10 else "ns"
        pvalues[g]['significance'] = sig

    # Sort by shift
    ranked = sorted(pvalues.keys(), key=lambda g: pvalues[g]['mean_shift'], reverse=True)

    # Print results table
    print("\n" + "=" * 100)
    print(f"{'Rank':>4} {'Gene':>10} {'Shift':>10} {'Null Mean':>10} {'Null N':>6} {'z':>7} {'p':>8} {'FDR q':>8} {'Sig':>4}")
    print("-" * 100)
    for rank, gene in enumerate(ranked, 1):
        p = pvalues[gene]
        z_str = f"{p['z_score']:.2f}" if p['z_score'] is not None else "N/A"
        nm_str = f"{p['null_mean']:.6f}" if p['null_mean'] is not None else "N/A"
        print(f"{rank:>4} {gene:>10} {p['mean_shift']:>10.6f} {nm_str:>10} {p['n_null_genes']:>6} {z_str:>7} {p['empirical_p']:>8.4f} {p['fdr_q']:>8.4f} {p['significance']:>4}")

    # Count significant
    n_sig_005 = sum(1 for g in ranked if pvalues[g]['fdr_q'] < 0.05)
    n_sig_01 = sum(1 for g in ranked if pvalues[g]['fdr_q'] < 0.01)
    n_sig_010 = sum(1 for g in ranked if pvalues[g]['fdr_q'] < 0.10)
    n_nom_005 = sum(1 for g in ranked if pvalues[g]['empirical_p'] < 0.05)

    print(f"\nNominally significant (p < 0.05): {n_nom_005}/{len(ranked)}")
    print(f"FDR < 0.10: {n_sig_010}/{len(ranked)}")
    print(f"FDR < 0.05: {n_sig_005}/{len(ranked)}")
    print(f"FDR < 0.01: {n_sig_01}/{len(ranked)}")

    # ===== Compute global statistics =====
    all_ctrl_shifts = []
    for b in range(n_bins):
        all_ctrl_shifts.extend(null_per_bin.get(b, []))
    all_intel_shifts = [pvalues[g]['mean_shift'] for g in ranked]

    # ===== Load pathway info =====
    rankings_path = os.path.join(OUTPUT_BASE, "final_rankings.json")
    pathways = {}
    if os.path.exists(rankings_path):
        with open(rankings_path) as f:
            rankings_data = json.load(f)
        for g, info in rankings_data.items():
            pathways[g] = info.get('pathway', '')

    # ===== Save JSON results =====
    json_output = {
        'metadata': {
            'analysis': 'Permutation null v2: 50 expression-matched control genes',
            'date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'n_intelligence_genes': len(pvalues),
            'n_control_genes_successful': len(successful),
            'n_existing_controls_included': len(existing_ctrl_shifts),
            'max_ncells': 500,
            'random_seed': 42,
        },
        'bin_info': bin_info,
        'intelligence_gene_results': {g: pvalues[g] for g in ranked},
        'control_gene_results': {k: v for k, v in control_results.items() if v.get('status') in ('success', 'reused')},
        'existing_control_results': existing_ctrl_shifts,
        'summary': {
            'n_significant_fdr005': n_sig_005,
            'n_significant_fdr01': n_sig_01,
            'n_significant_fdr010': n_sig_010,
            'n_nominally_significant': n_nom_005,
            'intelligence_mean_shift': float(np.mean(all_intel_shifts)),
            'intelligence_median_shift': float(np.median(all_intel_shifts)),
            'null_mean_shift': float(np.mean(all_ctrl_shifts)) if all_ctrl_shifts else None,
            'null_median_shift': float(np.median(all_ctrl_shifts)) if all_ctrl_shifts else None,
            'null_std_shift': float(np.std(all_ctrl_shifts)) if all_ctrl_shifts else None,
        }
    }

    json_path = os.path.join(ANALYSIS_DIR, "permutation_null_v2_results.json")
    with open(json_path, 'w') as f:
        json.dump(json_output, f, indent=2)
    print(f"\nSaved JSON: {json_path}")

    # ===== Generate Markdown report =====
    md = []
    md.append("# Permutation Null Experiment: Intelligence Gene Perturbation Effects")
    md.append("")
    md.append(f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}")
    md.append(f"**Dataset:** DLPFC_11k.h5ad (human dorsolateral prefrontal cortex)")
    md.append(f"**Model:** Geneformer V1 (10M parameters, pretrained)")
    md.append(f"**Method:** In silico gene deletion, 500 cells per gene, cell-level embedding shift (1 - cosine similarity)")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## Summary")
    md.append("")
    md.append("This analysis builds an empirical null distribution by running Geneformer")
    md.append("delete-perturbation on 50 random expression-matched control genes (10 per")
    md.append("expression quintile bin), plus 11 previously-run expression-matched controls.")
    md.append("Each control gene was run in a fresh subprocess to avoid OOM and broken-pipe")
    md.append("errors. Empirical p-values are computed per expression bin, with Benjamini-")
    md.append("Hochberg FDR correction across all 21 intelligence genes.")
    md.append("")
    md.append(f"- **Intelligence genes tested:** {len(pvalues)}")
    md.append(f"- **New control genes (permutation):** {len(successful)}")
    md.append(f"- **Existing controls reused:** {len(existing_ctrl_shifts)}")
    md.append(f"- **Total controls in null:** {len(all_ctrl_shifts)}")
    md.append("")
    md.append(f"- **Nominally significant (p < 0.05):** {n_nom_005}/{len(ranked)} genes")
    md.append(f"- **FDR < 0.10:** {n_sig_010}/{len(ranked)} genes")
    md.append(f"- **FDR < 0.05:** {n_sig_005}/{len(ranked)} genes")
    md.append(f"- **FDR < 0.01:** {n_sig_01}/{len(ranked)} genes")
    md.append("")
    md.append(f"- **Mean intelligence gene shift:** {np.mean(all_intel_shifts):.6f}")
    if all_ctrl_shifts:
        md.append(f"- **Mean null (control) shift:** {np.mean(all_ctrl_shifts):.6f}")
        md.append(f"- **Null std:** {np.std(all_ctrl_shifts):.6f}")
        md.append(f"- **Null median:** {np.median(all_ctrl_shifts):.6f}")
        md.append(f"- **Null 95th percentile:** {np.percentile(all_ctrl_shifts, 95):.6f}")
    md.append("")
    md.append("---")
    md.append("")

    md.append("## Method")
    md.append("")
    md.append("1. **Expression frequency computation:** For each gene in DLPFC_11k.h5ad,")
    md.append("   compute the fraction of cells expressing it (count > 0).")
    md.append("2. **Expression matching:** For each intelligence gene, identify all genes")
    md.append("   within +/- 20% expression frequency that are in Geneformer's token dictionary")
    md.append("   (excluding intelligence genes, existing controls, housekeeping genes).")
    md.append("3. **Stratified sampling:** Divide expression range into 5 equal-width bins.")
    md.append("   Sample 50 control genes proportionally across bins (approx 10 per bin).")
    md.append("4. **Perturbation:** Run Geneformer delete-perturbation on each control gene")
    md.append("   in a **fresh subprocess** (500 cells max, batch size 32). Each gene's result")
    md.append("   is saved incrementally to JSON before the next gene starts.")
    md.append("5. **Null augmentation:** Include 11 previously-run expression-matched controls")
    md.append("   (from initial paired control analysis) in the null distribution.")
    md.append("6. **Empirical p-value:** For each intelligence gene:")
    md.append("   p = (n_controls_in_same_bin_with_shift >= intel_shift + 1) / (n_controls_in_bin + 1)")
    md.append("7. **FDR correction:** Benjamini-Hochberg across all 21 genes.")
    md.append("")
    md.append("---")
    md.append("")

    md.append("## Results")
    md.append("")
    md.append("| Rank | Gene | Pathway | Shift | Null Mean | Null N | z-score | p (emp) | FDR q | Sig |")
    md.append("|------|------|---------|-------|-----------|--------|---------|---------|-------|-----|")
    for rank, gene in enumerate(ranked, 1):
        p = pvalues[gene]
        pathway = pathways.get(gene, '')
        z_str = f"{p['z_score']:.2f}" if p['z_score'] is not None else "N/A"
        nm_str = f"{p['null_mean']:.6f}" if p['null_mean'] is not None else "N/A"
        md.append(f"| {rank} | **{gene}** | {pathway} | {p['mean_shift']:.6f} | "
                  f"{nm_str} | {p['n_null_genes']} | {z_str} | "
                  f"{p['empirical_p']:.4f} | {p['fdr_q']:.4f} | {p['significance']} |")
    md.append("")
    md.append("*Significance: \\*\\*\\* FDR<0.001, \\*\\* FDR<0.01, \\* FDR<0.05, + FDR<0.10, ns = not significant*")
    md.append("")
    md.append("---")
    md.append("")

    md.append("## Expression Bin Details")
    md.append("")
    md.append("| Bin | Freq Range | Intel Genes | Controls in Null | Intel Gene Names |")
    md.append("|-----|-----------|-------------|-----------------|------------------|")
    for b in range(n_bins):
        intel_in_bin = bin_info['genes_per_bin'].get(str(b), [])
        n_ctrl = len(null_per_bin.get(b, []))
        gene_list = ", ".join(intel_in_bin) if intel_in_bin else "-"
        md.append(f"| {b} | {bin_edges[b]:.4f}-{bin_edges[b+1]:.4f} | "
                  f"{len(intel_in_bin)} | {n_ctrl} | {gene_list} |")
    md.append("")
    md.append("---")
    md.append("")

    md.append("## Null Distribution Statistics")
    md.append("")
    if all_ctrl_shifts:
        pcts = np.percentile(all_ctrl_shifts, [5, 25, 50, 75, 95])
        md.append(f"- **N total controls:** {len(all_ctrl_shifts)}")
        md.append(f"- **Mean:** {np.mean(all_ctrl_shifts):.6f}")
        md.append(f"- **Std:** {np.std(all_ctrl_shifts):.6f}")
        md.append(f"- **5th percentile:** {pcts[0]:.6f}")
        md.append(f"- **25th percentile:** {pcts[1]:.6f}")
        md.append(f"- **Median:** {pcts[2]:.6f}")
        md.append(f"- **75th percentile:** {pcts[3]:.6f}")
        md.append(f"- **95th percentile:** {pcts[4]:.6f}")
    md.append("")

    # Per-bin null stats
    md.append("### Per-Bin Null Statistics")
    md.append("")
    md.append("| Bin | N Controls | Mean Shift | Std | Min | Max |")
    md.append("|-----|-----------|-----------|-----|-----|-----|")
    for b in range(n_bins):
        shifts = null_per_bin.get(b, [])
        if shifts:
            md.append(f"| {b} | {len(shifts)} | {np.mean(shifts):.6f} | {np.std(shifts):.6f} | "
                      f"{min(shifts):.6f} | {max(shifts):.6f} |")
        else:
            md.append(f"| {b} | 0 | - | - | - | - |")
    md.append("")
    md.append("---")
    md.append("")

    md.append("## Interpretation")
    md.append("")
    if n_sig_005 > 0:
        md.append(f"Out of {len(pvalues)} intelligence-associated genes tested, **{n_sig_005}** show")
        md.append("significantly larger perturbation effects than expression-matched random genes")
        md.append("(FDR < 0.05). This provides evidence that intelligence-associated genes have")
        md.append("disproportionately large effects on brain cell transcriptomic states as modeled")
        md.append("by Geneformer, beyond what is expected from their expression level alone.")
    elif n_sig_010 > 0:
        md.append(f"Out of {len(pvalues)} intelligence-associated genes tested, **{n_sig_010}** show")
        md.append("marginally significant perturbation effects (FDR < 0.10) compared to expression-")
        md.append("matched random genes. While not reaching strict significance at FDR < 0.05,")
        md.append("this trend suggests intelligence genes may have larger effects than random genes")
        md.append("of similar expression level.")
    elif n_nom_005 > 0:
        md.append(f"While {n_nom_005} intelligence genes show nominally significant effects (p < 0.05),")
        md.append("none survive FDR correction. The intelligence genes as a group show a trend toward")
        md.append("larger perturbation effects, but individual gene significance is limited by the")
        md.append("size of the null distribution.")
    else:
        md.append("No individual intelligence genes reach significance against the expression-matched")
        md.append("null distribution. This does not invalidate the group-level finding (Mann-Whitney")
        md.append("p = 0.031 from the paired analysis) but indicates that the effect is distributed")
        md.append("across the gene set rather than driven by individual outliers.")
    md.append("")

    # Add comparison to previous analysis
    md.append("### Comparison with Previous Paired Analysis")
    md.append("")
    md.append("The original paired analysis (11 expression-matched controls) found:")
    md.append("- Mann-Whitney U: p = 0.031 (intelligence > controls)")
    md.append("- Cohen's d = 0.69 (medium effect)")
    md.append("- 9/10 intelligence genes > matched control (Wilcoxon p = 0.019)")
    md.append("")
    md.append(f"This expanded permutation null ({len(all_ctrl_shifts)} controls) provides a more")
    md.append("robust empirical null for gene-level inference, complementing the group-level")
    md.append("paired analysis.")
    md.append("")
    md.append("---")
    md.append("")

    md.append("## Technical Details")
    md.append("")
    md.append("- Each control gene was perturbed in a **fresh Python subprocess** to avoid")
    md.append("  memory leaks and broken-pipe errors from Geneformer's internal multiprocessing.")
    md.append("- Results were saved incrementally after each gene completed.")
    md.append("- Failed genes were logged and skipped (not retried).")
    md.append(f"- Total controls attempted: {len(control_results)}")
    md.append(f"- Successful: {len(successful)}")
    n_failed = len(control_results) - len(successful)
    md.append(f"- Failed: {n_failed}")
    md.append("")

    # List failed genes
    failed_genes = {k: v for k, v in control_results.items() if v.get('status') not in ('success', 'reused')}
    if failed_genes:
        md.append("### Failed Control Genes")
        md.append("")
        md.append("| Gene | Error |")
        md.append("|------|-------|")
        for sym, info in failed_genes.items():
            err = info.get('error', 'unknown')[:80]
            md.append(f"| {sym} | {err} |")
        md.append("")

    md.append("---")
    md.append("")
    md.append(f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*")

    md_path = os.path.join(ANALYSIS_DIR, "PERMUTATION_NULL_FINAL.md")
    with open(md_path, 'w') as f:
        f.write('\n'.join(md))
    print(f"Saved report: {md_path}")

    return json_path, md_path


if __name__ == "__main__":
    main()
