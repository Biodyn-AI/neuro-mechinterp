exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
"""
Step 1: Select 50 expression-matched control genes for the permutation null experiment.
Saves the selection to a JSON file that the per-gene runner will read.
"""
import os, sys, json, pickle, gc, random
import numpy as np
from collections import defaultdict

# ===== CONFIGURATION =====
OUTPUT_BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
BRAIN_H5AD = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
TOKEN_DICT_PATH = os.path.expanduser(
    "~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl"
)
ANALYSIS_DIR = "/mnt/d/openclaw/intelligence-augmentation/analysis"
N_CONTROL_GENES = 50
EXPRESSION_TOLERANCE = 0.20
RANDOM_SEED = 42
N_BINS = 5

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

HOUSEKEEPING = {
    "ACTB": "ENSG00000075624", "GAPDH": "ENSG00000111640", "B2M": "ENSG00000166710",
    "RPL13A": "ENSG00000142541", "YWHAZ": "ENSG00000164924",
}


def main():
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    # Load token dictionary
    print("Loading token dictionary...")
    with open(TOKEN_DICT_PATH, 'rb') as f:
        token_dict = pickle.load(f)
    print(f"  Token dictionary: {len(token_dict)} genes")

    # Load brain h5ad and compute expression frequencies
    print("Loading brain h5ad...")
    import anndata as ad
    import scipy.sparse as sp

    adata = ad.read_h5ad(BRAIN_H5AD)
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    var = adata.var
    if 'feature_name' in var.columns:
        gene_symbols = var['feature_name'].values
    else:
        gene_symbols = var.index.values
    ensembl_ids = var.index.values

    X = adata.X
    if sp.issparse(X):
        n_expressing = np.array((X > 0).sum(axis=0)).flatten()
    else:
        n_expressing = np.array((X > 0).sum(axis=0)).flatten()

    n_cells = adata.shape[0]
    expr_freq = n_expressing / n_cells

    gene_freq = {}
    for i, (ens, sym) in enumerate(zip(ensembl_ids, gene_symbols)):
        gene_freq[ens] = {'symbol': sym, 'frequency': float(expr_freq[i])}

    del adata, X
    gc.collect()

    # Find expression-matched genes for each intelligence gene
    print("\nFinding expression-matched candidate genes...")
    intel_ensembls = set(INTELLIGENCE_GENES.values())
    existing_ctrl_ensembls = set(EXISTING_CONTROLS.values())
    housekeeping_ensembls = set(HOUSEKEEPING.values())
    exclude = intel_ensembls | existing_ctrl_ensembls | housekeeping_ensembls

    # Pool all candidates
    all_candidates = {}
    for symbol, ensembl in INTELLIGENCE_GENES.items():
        if ensembl not in gene_freq:
            print(f"  WARNING: {symbol} ({ensembl}) not in h5ad")
            continue
        target_freq = gene_freq[ensembl]['frequency']
        lower = target_freq * (1 - EXPRESSION_TOLERANCE)
        upper = target_freq * (1 + EXPRESSION_TOLERANCE)

        for ens, info in gene_freq.items():
            if ens in exclude:
                continue
            if ens not in token_dict:
                continue
            if lower <= info['frequency'] <= upper:
                if ens not in all_candidates:
                    all_candidates[ens] = {'symbol': info['symbol'], 'frequency': info['frequency'], 'matched_to': []}
                all_candidates[ens]['matched_to'].append(symbol)

    print(f"  Total unique candidate genes: {len(all_candidates)}")

    # Create expression bins
    intel_freqs = sorted([gene_freq[ens]['frequency'] for ens in INTELLIGENCE_GENES.values() if ens in gene_freq])
    freq_min = min(intel_freqs) * (1 - EXPRESSION_TOLERANCE)
    freq_max = max(intel_freqs) * (1 + EXPRESSION_TOLERANCE)
    bin_edges = np.linspace(freq_min, freq_max, N_BINS + 1)

    # Assign intelligence genes to bins
    intel_bin_assignment = {}
    for symbol, ensembl in INTELLIGENCE_GENES.items():
        if ensembl not in gene_freq:
            continue
        freq = gene_freq[ensembl]['frequency']
        assigned = False
        for b in range(N_BINS):
            if bin_edges[b] <= freq <= bin_edges[b + 1]:
                intel_bin_assignment[symbol] = b
                assigned = True
                break
        if not assigned:
            dists = [abs(freq - (bin_edges[b] + bin_edges[b+1])/2) for b in range(N_BINS)]
            intel_bin_assignment[symbol] = int(np.argmin(dists))

    genes_per_bin = defaultdict(list)
    for symbol, b in intel_bin_assignment.items():
        genes_per_bin[b].append(symbol)

    # Assign candidates to bins
    candidates_per_bin = defaultdict(list)
    for ens, info in all_candidates.items():
        freq = info['frequency']
        assigned = False
        for b in range(N_BINS):
            if bin_edges[b] <= freq <= bin_edges[b + 1]:
                candidates_per_bin[b].append(ens)
                assigned = True
                break
        if not assigned:
            dists = [abs(freq - (bin_edges[b] + bin_edges[b+1])/2) for b in range(N_BINS)]
            candidates_per_bin[int(np.argmin(dists))].append(ens)

    # Allocate 10 per bin (50 total / 5 bins = 10)
    # Or proportionally if some bins have fewer candidates
    samples_per_bin = {}
    total_intel = sum(len(v) for v in genes_per_bin.values())
    allocated = 0
    for b in range(N_BINS):
        n_intel = len(genes_per_bin.get(b, []))
        if b < N_BINS - 1:
            n_samples = max(1, round(N_CONTROL_GENES * n_intel / total_intel))
        else:
            n_samples = N_CONTROL_GENES - allocated
        n_available = len(candidates_per_bin.get(b, []))
        n_samples = min(n_samples, n_available)
        samples_per_bin[b] = n_samples
        allocated += n_samples

    # If underallocated, fill from largest bins
    remaining = N_CONTROL_GENES - sum(samples_per_bin.values())
    if remaining > 0:
        for b in sorted(candidates_per_bin.keys(), key=lambda x: len(candidates_per_bin[x]), reverse=True):
            can_add = len(candidates_per_bin[b]) - samples_per_bin.get(b, 0)
            add = min(remaining, can_add)
            samples_per_bin[b] = samples_per_bin.get(b, 0) + add
            remaining -= add
            if remaining <= 0:
                break

    # Sample
    selected = []
    for b in range(N_BINS):
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
        print(f"  Bin {b} [{bin_edges[b]:.4f}-{bin_edges[b+1]:.4f}]: "
              f"{len(genes_per_bin.get(b, []))} intel genes, "
              f"{len(pool)} candidates, sampled {min(n_sample, len(pool))}")

    print(f"\nTotal selected control genes: {len(selected)}")

    # Save intelligence gene info
    intel_info = {}
    for symbol, ensembl in INTELLIGENCE_GENES.items():
        freq = gene_freq.get(ensembl, {}).get('frequency', None)
        intel_info[symbol] = {
            'ensembl': ensembl,
            'frequency': freq,
            'bin': intel_bin_assignment.get(symbol),
        }

    # Also save gene_freq for the selected + intel genes (for p-value computation)
    freq_subset = {}
    for ens in [c['ensembl'] for c in selected]:
        if ens in gene_freq:
            freq_subset[ens] = gene_freq[ens]
    for ens in INTELLIGENCE_GENES.values():
        if ens in gene_freq:
            freq_subset[ens] = gene_freq[ens]
    for ens in EXISTING_CONTROLS.values():
        if ens in gene_freq:
            freq_subset[ens] = gene_freq[ens]

    output = {
        'selected_controls': selected,
        'bin_info': {
            'bin_edges': [float(e) for e in bin_edges],
            'genes_per_bin': {str(b): v for b, v in genes_per_bin.items()},
            'samples_per_bin': {str(b): v for b, v in samples_per_bin.items()},
            'intel_bins': {s: {'bin': b, 'bin_range': [float(bin_edges[b]), float(bin_edges[b+1])]}
                          for s, b in intel_bin_assignment.items()},
        },
        'intelligence_genes': intel_info,
        'gene_frequencies': freq_subset,
    }

    out_path = os.path.join(ANALYSIS_DIR, "permutation_null_v2_selection.json")
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved selection to: {out_path}")

    # Print the selected genes
    for i, c in enumerate(selected):
        print(f"  {i+1}. {c['symbol']} ({c['ensembl']}) freq={c['frequency']:.4f} bin={c['bin']}")

if __name__ == "__main__":
    main()
