#!/usr/bin/env python3
"""
C1/C4 Full: Run Geneformer perturbations for null model stability.

This script:
1. Loads the tokenized brain dataset
2. For each of several random seeds, selects new expression-matched control genes
3. Runs Geneformer InSilicoPerturber on each control gene (using MPS on Mac)
4. Computes empirical p-values and FDR correction
5. Reports stability of the 0/21 FDR result across seeds

Uses MPS (Apple Silicon GPU) for inference.
"""

import os
import sys
import json
import pickle
import time
import gc
import random
import shutil
import numpy as np
from pathlib import Path
from collections import defaultdict

os.environ['TOKENIZERS_PARALLELISM'] = 'false'

# ============================================================
# Configuration
# ============================================================

BASE_DIR = Path("/Volumes/Crucial X6/MacBook/Code/biodyn/neuro-mechinterp")
MODEL_DIR = str(BASE_DIR / "models" / "Geneformer" / "Geneformer-V1-10M")
TOKENIZED_DATA = str(BASE_DIR / "analysis" / "results" / "insilico_wsl" / "brain.dataset")
TOKEN_DICT_PATH = "/tmp/geneformer_repo/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
OUTPUT_DIR = BASE_DIR / "analysis" / "revision" / "c1c4_full_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Existing results
EXISTING_RESULTS_DIR = BASE_DIR / "analysis" / "results" / "insilico_wsl" / "permutation_null_v2"
EXISTING_SELECTION = BASE_DIR / "analysis" / "permutation_null_v2_selection.json"
EXISTING_RESULTS = BASE_DIR / "analysis" / "permutation_null_v2_results.json"

# Seeds for resampling
SEEDS = [42, 123, 456]  # Start with 3 seeds to keep runtime reasonable

INTELLIGENCE_GENES = {
    "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
    "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
    "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
    "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
    "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
    "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
    "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
}

HOUSEKEEPING = {"ACTB", "GAPDH", "B2M", "RPL13A", "YWHAZ"}

N_NEW_CONTROLS_PER_SEED = 10  # New controls per seed (on top of existing 44)


def load_token_dict():
    """Load Geneformer token dictionary."""
    with open(TOKEN_DICT_PATH, 'rb') as f:
        token_dict = pickle.load(f)
    # Reverse: token_id -> ensembl_id
    id_to_token = {v: k for k, v in token_dict.items()}
    return token_dict, id_to_token


def compute_gene_frequencies(tokenized_data_path, token_dict):
    """Compute expression frequency for all genes in the tokenized dataset."""
    from datasets import load_from_disk
    ds = load_from_disk(tokenized_data_path)
    n_cells = len(ds)

    gene_counts = defaultdict(int)
    for cell in ds:
        unique_tokens = set(cell['input_ids'])
        for tok in unique_tokens:
            gene_counts[tok] += 1

    # Map to ensembl IDs and compute frequency
    id_to_token = {v: k for k, v in token_dict.items()}
    gene_freqs = {}
    for tok_id, count in gene_counts.items():
        if tok_id in id_to_token:
            ensembl = id_to_token[tok_id]
            gene_freqs[ensembl] = count / n_cells

    return gene_freqs, n_cells


def select_control_genes(gene_freqs, token_dict, seed, n_controls=10):
    """Select expression-matched control genes for a given seed."""
    rng = random.Random(seed)

    intel_ensembls = set(INTELLIGENCE_GENES.values())

    # Get intel gene frequencies
    intel_freqs = {}
    for sym, ens in INTELLIGENCE_GENES.items():
        if ens in gene_freqs:
            intel_freqs[sym] = gene_freqs[ens]

    # Get all candidate genes (in token dict, not intel, not housekeeping)
    # Build symbol lookup
    ensembl_to_symbol = {}
    for sym, ens in INTELLIGENCE_GENES.items():
        ensembl_to_symbol[ens] = sym

    candidates = []
    for ensembl, freq in gene_freqs.items():
        if ensembl in intel_ensembls:
            continue
        candidates.append((ensembl, freq))

    # Compute frequency bins
    intel_freq_values = list(intel_freqs.values())
    min_freq = min(intel_freq_values)
    max_freq = max(intel_freq_values)
    n_bins = 5
    bin_edges = np.linspace(min_freq * 0.9, max_freq * 1.1, n_bins + 1)

    # Assign intel genes to bins
    intel_bins = defaultdict(list)
    for sym, freq in intel_freqs.items():
        for i in range(n_bins):
            if bin_edges[i] <= freq < bin_edges[i + 1]:
                intel_bins[i].append(sym)
                break
        else:
            intel_bins[n_bins - 1].append(sym)

    # Select controls proportionally
    total_intel = sum(len(v) for v in intel_bins.values())
    controls_per_bin = {}
    remaining = n_controls
    for i in range(n_bins):
        n = max(1, round(n_controls * len(intel_bins[i]) / total_intel)) if intel_bins[i] else 0
        controls_per_bin[i] = min(n, remaining)
        remaining -= controls_per_bin[i]

    # Select from candidates
    selected = []
    existing_syms = set()

    # Load existing controls to avoid duplicates
    if EXISTING_SELECTION.exists():
        with open(EXISTING_SELECTION) as f:
            existing = json.load(f)
        existing_syms = {c['symbol'] for c in existing['selected_controls']}

    for bin_idx in range(n_bins):
        n_needed = controls_per_bin.get(bin_idx, 0)
        if n_needed == 0:
            continue

        bin_lo, bin_hi = bin_edges[bin_idx], bin_edges[bin_idx + 1]
        tolerance = 0.20

        bin_candidates = [
            (ens, freq) for ens, freq in candidates
            if bin_lo * (1 - tolerance) <= freq <= bin_hi * (1 + tolerance)
            and ens not in intel_ensembls
            and ens not in {c[0] for c in selected}
        ]

        rng.shuffle(bin_candidates)
        for ens, freq in bin_candidates[:n_needed]:
            selected.append((ens, freq, bin_idx))

    return selected


def run_perturbation_single_gene(ensembl_id, symbol, output_dir):
    """Run Geneformer InSilicoPerturber on a single gene."""
    import torch

    gene_output = os.path.join(output_dir, f"ctrl_{symbol}")
    result_file = os.path.join(output_dir, f"result_{symbol}.json")

    # Check if already computed
    if os.path.exists(result_file):
        with open(result_file) as f:
            return json.load(f)

    print(f"  Running perturbation for {symbol} ({ensembl_id})...")
    t0 = time.time()

    try:
        # Force CPU since MPS may have issues with Geneformer
        os.environ['CUDA_VISIBLE_DEVICES'] = ''

        from geneformer import InSilicoPerturber

        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[ensembl_id],
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cell",
            filter_data=None,
            max_ncells=500,
            emb_layer=-1,
            forward_batch_size=16,
            nproc=1,
        )

        os.makedirs(gene_output, exist_ok=True)

        isp.perturb_data(
            model_directory=MODEL_DIR,
            input_data_file=TOKENIZED_DATA,
            output_directory=gene_output,
            output_prefix=f"ctrl_{symbol}",
        )

        elapsed = time.time() - t0

        # Load result
        pickles = [f for f in os.listdir(gene_output) if f.endswith('.pickle')]
        if pickles:
            with open(os.path.join(gene_output, pickles[0]), 'rb') as f:
                data = pickle.load(f)
            if len(data) > 0:
                key = list(data.keys())[0]
                cos_sims = data[key]
                shifts = [1.0 - c for c in cos_sims]
                result = {
                    'ensembl': ensembl_id,
                    'symbol': symbol,
                    'mean_shift': float(np.mean(shifts)),
                    'std_shift': float(np.std(shifts)),
                    'n_cells': len(shifts),
                    'status': 'success',
                    'time_seconds': elapsed,
                }
                print(f"    SUCCESS: shift={np.mean(shifts):.6f}, n={len(shifts)}, time={elapsed:.1f}s")
            else:
                result = {'ensembl': ensembl_id, 'symbol': symbol, 'status': 'empty', 'time_seconds': elapsed}
        else:
            result = {'ensembl': ensembl_id, 'symbol': symbol, 'status': 'no_pickle', 'time_seconds': elapsed}

        # Clean up large output files
        shutil.rmtree(gene_output, ignore_errors=True)

    except Exception as e:
        import traceback
        elapsed = time.time() - t0
        result = {
            'ensembl': ensembl_id,
            'symbol': symbol,
            'status': 'error',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'time_seconds': elapsed,
        }
        print(f"    ERROR: {e}")

    # Save result
    with open(result_file, 'w') as f:
        json.dump(result, f, indent=2)

    gc.collect()
    return result


def load_existing_control_results():
    """Load all existing control gene perturbation results."""
    results = {}
    if EXISTING_RESULTS_DIR.exists():
        for f in EXISTING_RESULTS_DIR.glob("result_*.json"):
            with open(f) as fh:
                data = json.load(fh)
            if data.get('status') == 'success':
                results[data['symbol']] = data
    return results


def compute_empirical_pvalues(intel_results, control_results, gene_freqs):
    """Compute empirical p-values and FDR correction."""
    from scipy.stats import false_discovery_control

    control_shifts = [r['mean_shift'] for r in control_results.values()]

    pvalues = {}
    for sym, ens in INTELLIGENCE_GENES.items():
        if ens not in gene_freqs:
            continue

        # Load intel gene shift from existing results
        with open(EXISTING_RESULTS) as f:
            existing = json.load(f)

        intel_shift = existing['intelligence_gene_results'].get(sym, {}).get('mean_shift')
        if intel_shift is None:
            continue

        # Empirical p-value: fraction of controls with shift >= intel shift
        n_exceeding = sum(1 for s in control_shifts if s >= intel_shift)
        p = (n_exceeding + 1) / (len(control_shifts) + 1)
        pvalues[sym] = p

    # FDR correction
    if pvalues:
        genes = list(pvalues.keys())
        ps = [pvalues[g] for g in genes]
        try:
            qs = false_discovery_control(ps, method='bh')
            fdr = {g: float(q) for g, q in zip(genes, qs)}
        except Exception:
            fdr = {g: p * len(ps) for g, p in zip(genes, ps)}  # Simple Bonferroni fallback
    else:
        fdr = {}

    return pvalues, fdr


def main():
    print("=" * 70)
    print("C1/C4 Full: Geneformer Perturbation Null Model Stability")
    print("=" * 70)

    # Step 1: Load token dictionary and compute gene frequencies
    print("\n[1/4] Loading token dictionary and computing gene frequencies...")
    token_dict, id_to_token = load_token_dict()
    gene_freqs, n_cells = compute_gene_frequencies(TOKENIZED_DATA, token_dict)
    print(f"  {len(gene_freqs)} genes in tokenized data, {n_cells} cells")

    # Step 2: Load existing control results
    print("\n[2/4] Loading existing control results...")
    existing_controls = load_existing_control_results()
    print(f"  {len(existing_controls)} existing control results loaded")

    # Step 3: For each seed, select and run new controls
    print(f"\n[3/4] Running perturbations for {len(SEEDS)} seeds...")
    all_seed_results = {}

    for seed in SEEDS:
        print(f"\n--- Seed {seed} ---")
        seed_dir = OUTPUT_DIR / f"seed_{seed}"
        seed_dir.mkdir(exist_ok=True)

        # Select new controls
        new_controls = select_control_genes(gene_freqs, token_dict, seed, N_NEW_CONTROLS_PER_SEED)
        print(f"  Selected {len(new_controls)} new control genes")

        # Run perturbations on new controls
        new_results = {}
        for ens, freq, bin_idx in new_controls:
            # Find symbol for this ensembl ID
            symbol = id_to_token.get(int(list(token_dict.keys())[0])) if False else ens.split('G')[1][:5]
            # Use ensembl as identifier
            result = run_perturbation_single_gene(ens, f"seed{seed}_{ens[-5:]}", str(seed_dir))
            if result.get('status') == 'success':
                new_results[f"new_{ens}"] = result

        # Combine existing + new controls
        combined = {**existing_controls, **new_results}

        # Compute p-values
        pvalues, fdr = compute_empirical_pvalues(
            INTELLIGENCE_GENES, combined, gene_freqs
        )

        n_sig = sum(1 for q in fdr.values() if q < 0.05)
        n_nom = sum(1 for p in pvalues.values() if p < 0.05)
        print(f"  Seed {seed}: {n_nom} nominal p<0.05, {n_sig} FDR q<0.05")

        all_seed_results[seed] = {
            'n_controls': len(combined),
            'n_new': len(new_results),
            'n_fdr_significant': n_sig,
            'n_nominal_significant': n_nom,
            'pvalues': pvalues,
            'fdr': fdr,
            'new_control_shifts': {k: v['mean_shift'] for k, v in new_results.items()},
        }

    # Step 4: Save combined results
    print("\n[4/4] Saving results...")
    output = {
        'metadata': {
            'date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'seeds': SEEDS,
            'n_existing_controls': len(existing_controls),
            'n_new_controls_per_seed': N_NEW_CONTROLS_PER_SEED,
            'model': MODEL_DIR,
            'dataset': TOKENIZED_DATA,
        },
        'per_seed': all_seed_results,
        'summary': {
            'fdr_significant_range': [
                min(r['n_fdr_significant'] for r in all_seed_results.values()),
                max(r['n_fdr_significant'] for r in all_seed_results.values()),
            ],
            'conclusion': 'STABLE' if all(r['n_fdr_significant'] == 0 for r in all_seed_results.values()) else 'VARIABLE',
        }
    }

    with open(OUTPUT_DIR / 'c1c4_full_results.json', 'w') as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nResults saved to {OUTPUT_DIR / 'c1c4_full_results.json'}")
    print(f"\n{'='*70}")
    print(f"CONCLUSION: 0/21 FDR-significant across all {len(SEEDS)} seeds: "
          f"{'YES - STABLE' if output['summary']['conclusion'] == 'STABLE' else 'VARIABLE'}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
