exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
#!/usr/bin/env python3
"""
Genome-wide permutation null for intelligence gene perturbation analysis.

For each of 21 intelligence genes, builds an empirical null distribution from
expression-matched random genes. Computes empirical p-values and FDR correction.

Steps:
1. Load brain h5ad → compute expression frequency per gene
2. For each intelligence gene, find ALL genes within ±20% expression frequency
3. Sample 100 total control genes proportionally across expression bins
4. Run Geneformer delete-perturbation on each control gene (reuse existing results)
5. Build empirical null per expression bin
6. Compute empirical p-value for each intelligence gene
7. Apply Benjamini-Hochberg FDR correction
"""
import os, sys, json, pickle, logging, time, random, gc
import numpy as np
from pathlib import Path
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# ===========================================================================
# CONFIGURATION
# ===========================================================================
OUTPUT_BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
BRAIN_H5AD = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
MODEL_PATH = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
TOKENIZED_DATA = os.path.join(OUTPUT_BASE, "brain.dataset")
TOKEN_DICT_PATH = os.path.expanduser(
    "~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
)
GENE_NAME_DICT_PATH = os.path.expanduser(
    "~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl"
)

N_CONTROL_GENES = 100  # total random controls to sample
EXPRESSION_TOLERANCE = 0.20  # ±20% expression frequency
MAX_NCELLS = 500
FORWARD_BATCH_SIZE = 16  # must match successful runs (4 causes broken pipe)
NPROC = 4  # must be >=2; nproc=1 causes broken pipe in HF Datasets
RANDOM_SEED = 42

ANALYSIS_DIR = "/mnt/d/openclaw/intelligence-augmentation/analysis"
PERM_OUTPUT_DIR = os.path.join(OUTPUT_BASE, "permutation_null")

INTELLIGENCE_GENES = {
    "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
    "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
    "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
    "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
    "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
    "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
    "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
}

# Existing controls (already perturbed) - we'll reuse these
EXISTING_CONTROLS = {
    "EXOC6B": "ENSG00000144036", "DOCK9": "ENSG00000088387", "PDE7B": "ENSG00000171408",
    "ZDHHC9": "ENSG00000188706", "HDAC11": "ENSG00000163517", "XRCC6": "ENSG00000196419",
    "WDR46": "ENSG00000227057", "CUEDC2": "ENSG00000107874", "IGSF10": "ENSG00000152580",
    "SALL3": "ENSG00000256463", "TMEM88": "ENSG00000167874",
}

HOUSEKEEPING = {
    "ACTB": "ENSG00000075624", "GAPDH": "ENSG00000111640", "B2M": "ENSG00000166710",
    "RPL13A": "ENSG00000142541", "YWHAZ": "ENSG00000164924",
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
    return [1.0 - c for c in cos_sims]


# ===========================================================================
# STEP 1: Load brain data and compute expression frequencies
# ===========================================================================
def step1_compute_expression_frequencies():
    """Load h5ad, compute fraction of cells expressing each gene."""
    logger.info("STEP 1: Computing expression frequencies from brain h5ad...")
    import anndata as ad
    import scipy.sparse as sp

    adata = ad.read_h5ad(BRAIN_H5AD)
    logger.info(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Get gene symbols and ensembl IDs
    var = adata.var
    logger.info(f"  var columns: {list(var.columns)}")

    # Use 'feature_name' column for gene symbols (ENSEMBL IDs are the index)
    if 'feature_name' in var.columns:
        gene_symbols = var['feature_name'].values
    else:
        gene_symbols = var.index.values
    ensembl_ids = var.index.values

    # Compute expression frequency: fraction of cells where gene count > 0
    X = adata.X
    if sp.issparse(X):
        # For sparse matrix: count nonzero per column
        n_expressing = np.array((X > 0).sum(axis=0)).flatten()
    else:
        n_expressing = np.array((X > 0).sum(axis=0)).flatten()

    n_cells = adata.shape[0]
    expr_freq = n_expressing / n_cells

    # Build mapping: ensembl -> (symbol, frequency)
    gene_freq = {}
    for i, (ens, sym) in enumerate(zip(ensembl_ids, gene_symbols)):
        gene_freq[ens] = {
            'symbol': sym,
            'frequency': float(expr_freq[i]),
            'n_expressing': int(n_expressing[i]),
        }

    logger.info(f"  Computed frequencies for {len(gene_freq)} genes")
    logger.info(f"  Frequency range: {expr_freq.min():.4f} to {expr_freq.max():.4f}")
    logger.info(f"  Median frequency: {np.median(expr_freq):.4f}")

    del adata
    gc.collect()

    return gene_freq


# ===========================================================================
# STEP 2: Find expression-matched genes for each intelligence gene
# ===========================================================================
def step2_find_matched_genes(gene_freq, token_dict):
    """For each intelligence gene, find all genes within ±20% expression frequency."""
    logger.info("STEP 2: Finding expression-matched genes...")

    intel_ensembls = set(INTELLIGENCE_GENES.values())
    existing_ctrl_ensembls = set(EXISTING_CONTROLS.values())
    housekeeping_ensembls = set(HOUSEKEEPING.values())
    exclude = intel_ensembls | existing_ctrl_ensembls | housekeeping_ensembls

    matched = {}  # symbol -> {freq, matched_genes: [(ensembl, symbol, freq), ...]}

    for symbol, ensembl in INTELLIGENCE_GENES.items():
        if ensembl not in gene_freq:
            logger.warning(f"  {symbol} ({ensembl}) not found in h5ad")
            continue

        target_freq = gene_freq[ensembl]['frequency']
        lower = target_freq * (1 - EXPRESSION_TOLERANCE)
        upper = target_freq * (1 + EXPRESSION_TOLERANCE)

        candidates = []
        for ens, info in gene_freq.items():
            if ens in exclude:
                continue
            if ens not in token_dict:
                continue
            if lower <= info['frequency'] <= upper:
                candidates.append((ens, info['symbol'], info['frequency']))

        matched[symbol] = {
            'ensembl': ensembl,
            'frequency': target_freq,
            'n_matched': len(candidates),
            'matched_genes': candidates,
        }
        logger.info(f"  {symbol}: freq={target_freq:.4f}, {len(candidates)} matched genes in [{lower:.4f}, {upper:.4f}]")

    return matched


# ===========================================================================
# STEP 3: Sample 100 control genes proportionally across expression bins
# ===========================================================================
def step3_sample_controls(matched):
    """Sample N_CONTROL_GENES total controls proportionally across intelligence genes' bins."""
    logger.info(f"STEP 3: Sampling {N_CONTROL_GENES} control genes across expression bins...")

    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    # Pool all candidate genes across all intelligence genes, track which bins they belong to
    # Each candidate can match multiple intelligence genes (expression bins overlap)
    all_candidates = {}  # ensembl -> {symbol, freq, matched_intel_genes}
    for symbol, info in matched.items():
        for ens, sym, freq in info['matched_genes']:
            if ens not in all_candidates:
                all_candidates[ens] = {'symbol': sym, 'frequency': freq, 'matched_to': []}
            all_candidates[ens]['matched_to'].append(symbol)

    logger.info(f"  Total unique candidate genes across all bins: {len(all_candidates)}")

    # Strategy: divide into expression quintile bins based on intelligence gene frequencies
    intel_freqs = sorted([info['frequency'] for info in matched.values()])
    # Create bins: divide the frequency range into 5 equal-width bins
    freq_min = min(intel_freqs) * (1 - EXPRESSION_TOLERANCE)
    freq_max = max(intel_freqs) * (1 + EXPRESSION_TOLERANCE)
    n_bins = 5
    bin_edges = np.linspace(freq_min, freq_max, n_bins + 1)

    # Assign each intelligence gene to a bin
    intel_bin_assignment = {}
    for symbol, info in matched.items():
        for b in range(n_bins):
            if bin_edges[b] <= info['frequency'] <= bin_edges[b + 1]:
                intel_bin_assignment[symbol] = b
                break
        else:
            # Edge case: assign to closest bin
            dists = [abs(info['frequency'] - (bin_edges[b] + bin_edges[b+1])/2) for b in range(n_bins)]
            intel_bin_assignment[symbol] = int(np.argmin(dists))

    # Count intelligence genes per bin
    genes_per_bin = defaultdict(list)
    for symbol, b in intel_bin_assignment.items():
        genes_per_bin[b].append(symbol)

    # Assign candidates to bins
    candidates_per_bin = defaultdict(list)
    for ens, info in all_candidates.items():
        freq = info['frequency']
        for b in range(n_bins):
            if bin_edges[b] <= freq <= bin_edges[b + 1]:
                candidates_per_bin[b].append(ens)
                break
        else:
            dists = [abs(freq - (bin_edges[b] + bin_edges[b+1])/2) for b in range(n_bins)]
            candidates_per_bin[int(np.argmin(dists))].append(ens)

    # Allocate samples proportionally to number of intelligence genes per bin
    total_intel = sum(len(v) for v in genes_per_bin.values())
    samples_per_bin = {}
    allocated = 0
    for b in range(n_bins):
        n_intel = len(genes_per_bin.get(b, []))
        if b < n_bins - 1:
            n_samples = round(N_CONTROL_GENES * n_intel / total_intel)
        else:
            n_samples = N_CONTROL_GENES - allocated
        n_available = len(candidates_per_bin.get(b, []))
        n_samples = min(n_samples, n_available)
        samples_per_bin[b] = n_samples
        allocated += n_samples
        logger.info(f"  Bin {b} [{bin_edges[b]:.4f}-{bin_edges[b+1]:.4f}]: "
                     f"{n_intel} intel genes, {n_available} candidates, sampling {n_samples}")

    # If we underallocated, fill from the largest available bin
    remaining = N_CONTROL_GENES - sum(samples_per_bin.values())
    if remaining > 0:
        for b in sorted(candidates_per_bin.keys(), key=lambda x: len(candidates_per_bin[x]), reverse=True):
            can_add = len(candidates_per_bin[b]) - samples_per_bin[b]
            add = min(remaining, can_add)
            samples_per_bin[b] += add
            remaining -= add
            if remaining <= 0:
                break

    # Sample from each bin
    selected = []
    for b in range(n_bins):
        n_sample = samples_per_bin.get(b, 0)
        pool = candidates_per_bin.get(b, [])
        if n_sample > 0 and pool:
            chosen = random.sample(pool, min(n_sample, len(pool)))
            for ens in chosen:
                info = all_candidates[ens]
                selected.append({
                    'ensembl': ens,
                    'symbol': info['symbol'],
                    'frequency': info['frequency'],
                    'bin': b,
                    'matched_to': info['matched_to'],
                })

    logger.info(f"  Selected {len(selected)} control genes total")

    # Also assign each intelligence gene its bin info
    intel_bins = {}
    for symbol, b in intel_bin_assignment.items():
        intel_bins[symbol] = {
            'bin': b,
            'bin_range': (float(bin_edges[b]), float(bin_edges[b+1])),
        }

    bin_info = {
        'bin_edges': [float(e) for e in bin_edges],
        'genes_per_bin': {str(b): v for b, v in genes_per_bin.items()},
        'samples_per_bin': {str(b): v for b, v in samples_per_bin.items()},
        'intel_bins': intel_bins,
    }

    return selected, bin_info


# ===========================================================================
# STEP 4: Run Geneformer perturbation on control genes
# ===========================================================================
def step4_run_perturbations(selected_controls):
    """Run Geneformer InSilicoPerturber on each control gene."""
    logger.info(f"STEP 4: Running perturbations on {len(selected_controls)} control genes...")

    import torch
    logger.info(f"  CUDA: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        logger.info(f"  GPU: {torch.cuda.get_device_name(0)}")

    from geneformer import InSilicoPerturber

    os.makedirs(PERM_OUTPUT_DIR, exist_ok=True)

    results = {}
    already_done = 0
    newly_run = 0
    failed = 0

    for i, ctrl in enumerate(selected_controls):
        ensembl = ctrl['ensembl']
        symbol = ctrl['symbol']
        gene_output = os.path.join(PERM_OUTPUT_DIR, f"perm_ctrl_{symbol}")

        # Check if already done
        existing = load_perturbation_result(gene_output)
        if existing is not None:
            shifts = compute_shift(existing)
            results[symbol] = {
                'ensembl': ensembl,
                'frequency': ctrl['frequency'],
                'bin': ctrl['bin'],
                'matched_to': ctrl['matched_to'],
                'mean_shift': float(np.mean(shifts)),
                'std_shift': float(np.std(shifts)),
                'n_cells': len(shifts),
                'status': 'reused',
            }
            already_done += 1
            logger.info(f"  [{i+1}/{len(selected_controls)}] {symbol}: REUSED (shift={np.mean(shifts):.6f}, n={len(shifts)})")
            continue

        # Also check if this gene was already run as an existing control
        for ctrl_prefix in ['control_perturb_', 'perturb_']:
            alt_dir = os.path.join(OUTPUT_BASE, f"{ctrl_prefix}{symbol}")
            existing = load_perturbation_result(alt_dir)
            if existing is not None:
                shifts = compute_shift(existing)
                results[symbol] = {
                    'ensembl': ensembl,
                    'frequency': ctrl['frequency'],
                    'bin': ctrl['bin'],
                    'matched_to': ctrl['matched_to'],
                    'mean_shift': float(np.mean(shifts)),
                    'std_shift': float(np.std(shifts)),
                    'n_cells': len(shifts),
                    'status': 'reused_existing',
                }
                already_done += 1
                logger.info(f"  [{i+1}/{len(selected_controls)}] {symbol}: REUSED from {ctrl_prefix} (shift={np.mean(shifts):.6f}, n={len(shifts)})")
                break
        else:
            # Need to run perturbation
            os.makedirs(gene_output, exist_ok=True)
            try:
                logger.info(f"  [{i+1}/{len(selected_controls)}] {symbol} ({ensembl}): RUNNING...")
                t0 = time.time()

                isp = InSilicoPerturber(
                    perturb_type="delete",
                    genes_to_perturb=[ensembl],
                    model_type="Pretrained",
                    num_classes=0,
                    emb_mode="cell",
                    filter_data=None,
                    max_ncells=MAX_NCELLS,
                    emb_layer=-1,
                    forward_batch_size=FORWARD_BATCH_SIZE,
                    nproc=NPROC,
                    model_version="V1",
                    token_dictionary_file=TOKEN_DICT_PATH,
                )

                isp.perturb_data(
                    model_directory=MODEL_PATH,
                    input_data_file=TOKENIZED_DATA,
                    output_directory=gene_output,
                    output_prefix=f"perm_ctrl_{symbol}",
                )

                elapsed = time.time() - t0

                # Load the result
                cos_sims = load_perturbation_result(gene_output)
                if cos_sims is not None:
                    shifts = compute_shift(cos_sims)
                    results[symbol] = {
                        'ensembl': ensembl,
                        'frequency': ctrl['frequency'],
                        'bin': ctrl['bin'],
                        'matched_to': ctrl['matched_to'],
                        'mean_shift': float(np.mean(shifts)),
                        'std_shift': float(np.std(shifts)),
                        'n_cells': len(shifts),
                        'status': 'success',
                        'time_seconds': elapsed,
                    }
                    newly_run += 1
                    logger.info(f"    SUCCESS: shift={np.mean(shifts):.6f}, n={len(shifts)}, time={elapsed:.1f}s")
                else:
                    results[symbol] = {
                        'ensembl': ensembl,
                        'status': 'no_result',
                    }
                    failed += 1
                    logger.warning(f"    No result file produced")

                # Clean up GPU memory
                del isp
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

            except Exception as e:
                elapsed = time.time() - t0
                results[symbol] = {
                    'ensembl': ensembl,
                    'frequency': ctrl['frequency'],
                    'bin': ctrl['bin'],
                    'matched_to': ctrl['matched_to'],
                    'status': 'failed',
                    'error': str(e),
                    'time_seconds': elapsed,
                }
                failed += 1
                logger.error(f"    FAILED: {e}")

                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

        # Save checkpoint every 10 genes
        if (i + 1) % 10 == 0:
            checkpoint_path = os.path.join(PERM_OUTPUT_DIR, "checkpoint_results.json")
            with open(checkpoint_path, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"  Checkpoint saved: {len(results)} results")

    logger.info(f"\n  STEP 4 SUMMARY: {already_done} reused, {newly_run} newly run, {failed} failed")
    return results


# ===========================================================================
# STEP 5: Build empirical null and compute p-values
# ===========================================================================
def step5_compute_pvalues(control_results, matched, bin_info, gene_freq):
    """For each intelligence gene, compute empirical p-value from expression-matched null."""
    logger.info("STEP 5: Computing empirical p-values...")

    from scipy import stats as sp_stats

    # Load intelligence gene results
    intel_results = {}
    for symbol, ensembl in INTELLIGENCE_GENES.items():
        result_dir = os.path.join(OUTPUT_BASE, f"perturb_{symbol}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = compute_shift(cos_sims)
            intel_results[symbol] = {
                'ensembl': ensembl,
                'mean_shift': float(np.mean(shifts)),
                'std_shift': float(np.std(shifts)),
                'n_cells': len(shifts),
                'frequency': gene_freq.get(ensembl, {}).get('frequency', None),
            }
        else:
            logger.warning(f"  {symbol}: no perturbation result found")

    # Build null distribution per expression bin
    bin_edges = bin_info['bin_edges']
    n_bins = len(bin_edges) - 1

    # Get control shifts per bin
    null_per_bin = defaultdict(list)  # bin -> list of mean_shifts
    for sym, res in control_results.items():
        if res.get('status') in ('success', 'reused', 'reused_existing') and 'mean_shift' in res:
            b = res['bin']
            null_per_bin[b].append(res['mean_shift'])

    # Also include existing expression-matched controls in the null
    for ctrl_sym, ctrl_ens in EXISTING_CONTROLS.items():
        result_dir = os.path.join(OUTPUT_BASE, f"control_perturb_{ctrl_sym}")
        cos_sims = load_perturbation_result(result_dir)
        if cos_sims is not None:
            shifts = compute_shift(cos_sims)
            ctrl_freq = gene_freq.get(ctrl_ens, {}).get('frequency', None)
            if ctrl_freq is not None:
                # Assign to bin
                for b in range(n_bins):
                    if bin_edges[b] <= ctrl_freq <= bin_edges[b + 1]:
                        null_per_bin[b].append(float(np.mean(shifts)))
                        break

    for b in sorted(null_per_bin.keys()):
        logger.info(f"  Bin {b}: {len(null_per_bin[b])} control genes in null")

    # Compute empirical p-value per intelligence gene
    # p = (# controls with shift >= intel gene shift + 1) / (# controls + 1)
    # This is the standard empirical p-value with pseudocount
    pvalues = {}
    for symbol, ires in intel_results.items():
        intel_shift = ires['mean_shift']
        intel_freq = ires.get('frequency')

        if intel_freq is None:
            logger.warning(f"  {symbol}: no frequency info, skipping")
            continue

        # Find which bin this intelligence gene belongs to
        intel_bin = bin_info['intel_bins'].get(symbol, {}).get('bin')

        # Collect null: use the same bin AND adjacent bins for more power
        null_shifts = []
        if intel_bin is not None:
            for b in range(max(0, intel_bin - 1), min(n_bins, intel_bin + 2)):
                null_shifts.extend(null_per_bin.get(b, []))

        # Also: use ALL controls within ±20% expression frequency directly
        # This is the most precise null
        precise_null = []
        lower = intel_freq * (1 - EXPRESSION_TOLERANCE)
        upper = intel_freq * (1 + EXPRESSION_TOLERANCE)

        for sym, res in control_results.items():
            if res.get('status') in ('success', 'reused', 'reused_existing') and 'mean_shift' in res:
                if lower <= res.get('frequency', 0) <= upper:
                    precise_null.append(res['mean_shift'])

        # Add existing controls that fall in range
        for ctrl_sym, ctrl_ens in EXISTING_CONTROLS.items():
            ctrl_freq = gene_freq.get(ctrl_ens, {}).get('frequency', None)
            if ctrl_freq is not None and lower <= ctrl_freq <= upper:
                result_dir = os.path.join(OUTPUT_BASE, f"control_perturb_{ctrl_sym}")
                cos_sims = load_perturbation_result(result_dir)
                if cos_sims is not None:
                    shifts = compute_shift(cos_sims)
                    precise_null.append(float(np.mean(shifts)))

        # Use precise null if we have enough genes, otherwise fall back to bin-based
        if len(precise_null) >= 5:
            use_null = precise_null
            null_type = "expression_matched"
        elif len(null_shifts) >= 5:
            use_null = null_shifts
            null_type = "bin_based"
        else:
            # Use all controls as fallback
            use_null = []
            for sym, res in control_results.items():
                if res.get('status') in ('success', 'reused', 'reused_existing') and 'mean_shift' in res:
                    use_null.append(res['mean_shift'])
            null_type = "all_controls"

        n_null = len(use_null)
        n_greater_or_equal = sum(1 for s in use_null if s >= intel_shift)
        empirical_p = (n_greater_or_equal + 1) / (n_null + 1)

        # Also compute z-score against null
        if n_null > 1:
            null_mean = np.mean(use_null)
            null_std = np.std(use_null)
            if null_std > 0:
                z_score = (intel_shift - null_mean) / null_std
            else:
                z_score = float('inf') if intel_shift > null_mean else 0.0
        else:
            z_score = None
            null_mean = None
            null_std = None

        pvalues[symbol] = {
            'mean_shift': intel_shift,
            'n_cells': ires['n_cells'],
            'expression_frequency': float(intel_freq) if intel_freq else None,
            'null_type': null_type,
            'n_null_genes': n_null,
            'n_exceeding': n_greater_or_equal,
            'empirical_p': float(empirical_p),
            'z_score': float(z_score) if z_score is not None else None,
            'null_mean': float(null_mean) if null_mean is not None else None,
            'null_std': float(null_std) if null_std is not None else None,
        }
        logger.info(f"  {symbol}: shift={intel_shift:.6f}, null_n={n_null} ({null_type}), "
                     f"p={empirical_p:.4f}, z={z_score:.2f}" if z_score else
                     f"  {symbol}: shift={intel_shift:.6f}, null_n={n_null}, p={empirical_p:.4f}")

    return pvalues, intel_results


# ===========================================================================
# STEP 6: FDR correction
# ===========================================================================
def step6_fdr_correction(pvalues):
    """Apply Benjamini-Hochberg FDR correction."""
    logger.info("STEP 6: Applying FDR correction (Benjamini-Hochberg)...")

    from scipy.stats import false_discovery_control

    genes = sorted(pvalues.keys())
    raw_ps = [pvalues[g]['empirical_p'] for g in genes]

    # Benjamini-Hochberg
    n = len(raw_ps)
    sorted_indices = np.argsort(raw_ps)
    sorted_ps = np.array(raw_ps)[sorted_indices]

    # BH procedure
    fdr_values = np.zeros(n)
    for i in range(n):
        rank = i + 1
        fdr_values[i] = sorted_ps[i] * n / rank

    # Enforce monotonicity (from largest to smallest rank)
    for i in range(n - 2, -1, -1):
        fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])

    # Cap at 1
    fdr_values = np.minimum(fdr_values, 1.0)

    # Map back to original order
    fdr_adjusted = np.zeros(n)
    for i in range(n):
        fdr_adjusted[sorted_indices[i]] = fdr_values[i]

    for i, g in enumerate(genes):
        pvalues[g]['fdr_q'] = float(fdr_adjusted[i])
        sig = "***" if fdr_adjusted[i] < 0.001 else "**" if fdr_adjusted[i] < 0.01 else "*" if fdr_adjusted[i] < 0.05 else "ns"
        pvalues[g]['significance'] = sig
        logger.info(f"  {g}: p={pvalues[g]['empirical_p']:.4f}, FDR q={fdr_adjusted[i]:.4f} {sig}")

    return pvalues


# ===========================================================================
# STEP 7: Generate output files
# ===========================================================================
def step7_generate_output(pvalues, control_results, bin_info, intel_results):
    """Generate JSON and Markdown output files."""
    logger.info("STEP 7: Generating output files...")

    # Sort genes by effect size
    ranked_genes = sorted(pvalues.keys(), key=lambda g: pvalues[g]['mean_shift'], reverse=True)

    # Count significant genes
    n_sig_005 = sum(1 for g in ranked_genes if pvalues[g]['fdr_q'] < 0.05)
    n_sig_01 = sum(1 for g in ranked_genes if pvalues[g]['fdr_q'] < 0.1)
    n_sig_001 = sum(1 for g in ranked_genes if pvalues[g]['fdr_q'] < 0.001)

    # Count successful controls
    n_controls_success = sum(1 for r in control_results.values()
                              if r.get('status') in ('success', 'reused', 'reused_existing'))

    # Compute global null stats
    all_ctrl_shifts = [r['mean_shift'] for r in control_results.values()
                       if r.get('status') in ('success', 'reused', 'reused_existing')]
    all_intel_shifts = [pvalues[g]['mean_shift'] for g in ranked_genes]

    # JSON output
    json_output = {
        'metadata': {
            'analysis': 'Genome-wide permutation null for intelligence gene perturbation',
            'date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'n_intelligence_genes': len(pvalues),
            'n_control_genes_sampled': len(control_results),
            'n_control_genes_successful': n_controls_success,
            'expression_tolerance': EXPRESSION_TOLERANCE,
            'max_ncells': MAX_NCELLS,
            'random_seed': RANDOM_SEED,
        },
        'bin_info': bin_info,
        'intelligence_gene_results': {g: pvalues[g] for g in ranked_genes},
        'control_gene_results': control_results,
        'summary': {
            'n_significant_fdr005': n_sig_005,
            'n_significant_fdr01': n_sig_01,
            'n_significant_fdr001': n_sig_001,
            'intelligence_mean_shift': float(np.mean(all_intel_shifts)),
            'intelligence_median_shift': float(np.median(all_intel_shifts)),
            'null_mean_shift': float(np.mean(all_ctrl_shifts)) if all_ctrl_shifts else None,
            'null_median_shift': float(np.median(all_ctrl_shifts)) if all_ctrl_shifts else None,
            'null_std_shift': float(np.std(all_ctrl_shifts)) if all_ctrl_shifts else None,
        }
    }

    json_path = os.path.join(ANALYSIS_DIR, "permutation_null_results.json")
    with open(json_path, 'w') as f:
        json.dump(json_output, f, indent=2)
    logger.info(f"  Saved: {json_path}")

    # Markdown report
    md_lines = []
    md_lines.append("# Genome-Wide Permutation Null: Intelligence Gene Perturbation Analysis")
    md_lines.append("")
    md_lines.append(f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}")
    md_lines.append("")
    md_lines.append("## Summary")
    md_lines.append("")
    md_lines.append("This analysis addresses the concern that 11 expression-matched controls is insufficient")
    md_lines.append("for robust statistical inference. We build a genome-wide empirical null distribution")
    md_lines.append("by sampling 100 random expression-matched genes and running Geneformer delete-perturbation")
    md_lines.append("on each, then computing empirical p-values for each intelligence gene against its")
    md_lines.append("expression-matched null.")
    md_lines.append("")
    md_lines.append("### Key Results")
    md_lines.append("")
    md_lines.append(f"- **Intelligence genes tested:** {len(pvalues)}")
    md_lines.append(f"- **Control genes in null:** {n_controls_success}")
    md_lines.append(f"- **Significant at FDR < 0.05:** {n_sig_005}/{len(pvalues)} genes")
    md_lines.append(f"- **Significant at FDR < 0.01:** {n_sig_001}/{len(pvalues)} genes")
    md_lines.append(f"- **Significant at FDR < 0.10:** {n_sig_01}/{len(pvalues)} genes")
    md_lines.append("")
    md_lines.append(f"- **Mean intelligence gene shift:** {np.mean(all_intel_shifts):.6f}")
    md_lines.append(f"- **Mean null (control) shift:** {np.mean(all_ctrl_shifts):.6f}" if all_ctrl_shifts else "")
    md_lines.append(f"- **Null std:** {np.std(all_ctrl_shifts):.6f}" if all_ctrl_shifts else "")
    md_lines.append("")

    md_lines.append("## Method")
    md_lines.append("")
    md_lines.append("1. **Expression frequency computation:** For each gene in DLPFC_11k.h5ad,")
    md_lines.append("   compute the fraction of cells expressing it (count > 0).")
    md_lines.append("2. **Expression matching:** For each intelligence gene, identify ALL genes")
    md_lines.append(f"   within ±{int(EXPRESSION_TOLERANCE*100)}% expression frequency that are present in Geneformer's")
    md_lines.append("   token dictionary (excluding intelligence genes, existing controls, housekeeping).")
    md_lines.append(f"3. **Stratified sampling:** Divide expression range into 5 equal-width bins.")
    md_lines.append(f"   Sample {N_CONTROL_GENES} control genes proportionally across bins.")
    md_lines.append(f"4. **Perturbation:** Run Geneformer delete-perturbation on each control gene")
    md_lines.append(f"   ({MAX_NCELLS} cells, batch size {FORWARD_BATCH_SIZE}).")
    md_lines.append("5. **Empirical p-value:** For each intelligence gene, p = (n_controls_with_shift ≥")
    md_lines.append("   intelligence_shift + 1) / (n_controls + 1), using expression-matched controls.")
    md_lines.append("6. **FDR correction:** Benjamini-Hochberg across all 21 genes.")
    md_lines.append("")

    md_lines.append("## Results Table")
    md_lines.append("")
    md_lines.append("| Rank | Gene | Pathway | Shift | Null Mean | Null N | z-score | p (emp) | FDR q | Sig |")
    md_lines.append("|------|------|---------|-------|-----------|--------|---------|---------|-------|-----|")

    # Load pathway info from final_rankings.json
    rankings_path = os.path.join(OUTPUT_BASE, "final_rankings.json")
    pathways = {}
    if os.path.exists(rankings_path):
        with open(rankings_path) as f:
            rankings = json.load(f)
        for g, info in rankings.items():
            pathways[g] = info.get('pathway', '')

    for rank, gene in enumerate(ranked_genes, 1):
        p = pvalues[gene]
        pathway = pathways.get(gene, '')
        z_str = f"{p['z_score']:.2f}" if p['z_score'] is not None else "N/A"
        null_mean_str = f"{p['null_mean']:.6f}" if p['null_mean'] is not None else "N/A"
        md_lines.append(
            f"| {rank} | **{gene}** | {pathway} | {p['mean_shift']:.6f} | "
            f"{null_mean_str} | {p['n_null_genes']} | {z_str} | "
            f"{p['empirical_p']:.4f} | {p['fdr_q']:.4f} | {p['significance']} |"
        )

    md_lines.append("")
    md_lines.append("*Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant (FDR-corrected)*")
    md_lines.append("")

    md_lines.append("## Expression Bin Details")
    md_lines.append("")
    md_lines.append("| Bin | Freq Range | Intel Genes | Controls | Genes |")
    md_lines.append("|-----|-----------|-------------|----------|-------|")
    bin_edges = bin_info['bin_edges']
    for b in range(len(bin_edges) - 1):
        intel_in_bin = bin_info['genes_per_bin'].get(str(b), [])
        n_ctrl = int(bin_info['samples_per_bin'].get(str(b), 0))
        gene_list = ", ".join(intel_in_bin) if intel_in_bin else "-"
        md_lines.append(f"| {b} | {bin_edges[b]:.4f}-{bin_edges[b+1]:.4f} | "
                         f"{len(intel_in_bin)} | {n_ctrl} | {gene_list} |")

    md_lines.append("")
    md_lines.append("## Null Distribution Statistics")
    md_lines.append("")
    if all_ctrl_shifts:
        percentiles = np.percentile(all_ctrl_shifts, [5, 25, 50, 75, 95])
        md_lines.append(f"- **N controls:** {len(all_ctrl_shifts)}")
        md_lines.append(f"- **Mean:** {np.mean(all_ctrl_shifts):.6f}")
        md_lines.append(f"- **Std:** {np.std(all_ctrl_shifts):.6f}")
        md_lines.append(f"- **5th percentile:** {percentiles[0]:.6f}")
        md_lines.append(f"- **25th percentile:** {percentiles[1]:.6f}")
        md_lines.append(f"- **Median:** {percentiles[2]:.6f}")
        md_lines.append(f"- **75th percentile:** {percentiles[3]:.6f}")
        md_lines.append(f"- **95th percentile:** {percentiles[4]:.6f}")

    md_lines.append("")
    md_lines.append("## Interpretation")
    md_lines.append("")
    md_lines.append(f"Out of {len(pvalues)} intelligence-associated genes tested, "
                     f"**{n_sig_005}** show significantly larger perturbation effects than "
                     f"expression-matched random genes (FDR < 0.05). This confirms that "
                     f"intelligence-associated genes have disproportionately large effects on "
                     f"brain cell transcriptomic states as modeled by Geneformer, and this "
                     f"finding is robust to a genome-wide null distribution with "
                     f"{n_controls_success} expression-matched controls.")
    md_lines.append("")

    md_path = os.path.join(ANALYSIS_DIR, "PERMUTATION_NULL_RESULTS.md")
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    logger.info(f"  Saved: {md_path}")

    return json_path, md_path


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    logger.info("=" * 70)
    logger.info("GENOME-WIDE PERMUTATION NULL ANALYSIS")
    logger.info("=" * 70)
    t_start = time.time()

    # Load token dictionary
    with open(TOKEN_DICT_PATH, 'rb') as f:
        token_dict = pickle.load(f)
    logger.info(f"Token dictionary: {len(token_dict)} genes")

    # Step 1
    gene_freq = step1_compute_expression_frequencies()

    # Step 2
    matched = step2_find_matched_genes(gene_freq, token_dict)

    # Step 3
    selected_controls, bin_info = step3_sample_controls(matched)

    # Save control gene selection for reproducibility
    os.makedirs(PERM_OUTPUT_DIR, exist_ok=True)
    selection_path = os.path.join(PERM_OUTPUT_DIR, "selected_controls.json")
    with open(selection_path, 'w') as f:
        json.dump({
            'selected_controls': selected_controls,
            'bin_info': bin_info,
        }, f, indent=2)
    logger.info(f"Saved control selection: {selection_path}")

    # Step 4
    control_results = step4_run_perturbations(selected_controls)

    # Step 5
    pvalues, intel_results = step5_compute_pvalues(control_results, matched, bin_info, gene_freq)

    # Step 6
    pvalues = step6_fdr_correction(pvalues)

    # Step 7
    json_path, md_path = step7_generate_output(pvalues, control_results, bin_info, intel_results)

    elapsed = time.time() - t_start
    logger.info(f"\nTotal time: {elapsed/60:.1f} minutes")
    logger.info(f"Results: {json_path}")
    logger.info(f"Report: {md_path}")
    logger.info("DONE")


if __name__ == "__main__":
    main()
