#!/usr/bin/env python3
"""Phase 2: Combinatorial perturbations - pairwise deletion of top 5 genes."""
import os, sys, pickle, json, logging, torch
import numpy as np
from itertools import combinations

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
tokenized_data = os.path.join(output_dir, "brain.dataset")
token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")

top5 = {
    "CADM2": "ENSG00000175161",
    "GRIN2A": "ENSG00000183454",
    "CAMK2A": "ENSG00000070808",
    "MEF2C": "ENSG00000081189",
    "APP": "ENSG00000142192",
}

logger.info(f"CUDA: {torch.cuda.is_available()}, GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A'}")

from geneformer import InSilicoPerturber

pairs = list(combinations(top5.keys(), 2))
logger.info(f"Running {len(pairs)} pairwise perturbations")

results = {}
for i, (g1, g2) in enumerate(pairs):
    pair_name = f"{g1}_{g2}"
    logger.info(f"\n=== [{i+1}/{len(pairs)}] Perturbing {g1} + {g2} ===")
    
    pair_output = os.path.join(output_dir, f"perturb_combo_{pair_name}")
    os.makedirs(pair_output, exist_ok=True)
    
    try:
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[top5[g1], top5[g2]],
            combos=1,  # combinatorial mode
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cell",
            filter_data=None,
            max_ncells=500,
            emb_layer=-1,
            forward_batch_size=32,
            nproc=1,
            model_version="V1",
            token_dictionary_file=token_dict_path,
        )
        
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data,
            output_directory=pair_output,
            output_prefix=f"combo_{pair_name}",
        )
        
        # Load results
        pkl_files = [f for f in os.listdir(pair_output) if f.endswith('.pickle')]
        logger.info(f"  Output files: {pkl_files}")
        
        if pkl_files:
            with open(os.path.join(pair_output, pkl_files[0]), 'rb') as f:
                data = pickle.load(f)
            # Get cosine similarities
            for key in data:
                if 'cell_emb' in str(key):
                    cos_sims = data[key]
                    shift = 1.0 - np.mean(cos_sims)
                    results[pair_name] = {
                        'genes': [g1, g2],
                        'mean_shift': float(shift),
                        'std': float(np.std([1.0 - s for s in cos_sims])),
                        'n_cells': len(cos_sims),
                        'key': str(key),
                    }
                    logger.info(f"  {pair_name}: shift={shift:.4f}, n={len(cos_sims)}")
                    break
        
        results[pair_name]['status'] = 'success'
    except Exception as e:
        logger.error(f"  {pair_name}: FAILED - {e}")
        results[pair_name] = {'status': 'failed', 'error': str(e), 'genes': [g1, g2]}
    
    # Clear CUDA cache
    torch.cuda.empty_cache()

# Save results
with open(os.path.join(output_dir, 'combinatorial_results.json'), 'w') as f:
    json.dump(results, f, indent=2)

# Load individual results for synergy analysis
individual_shifts = {}
token_map = {"CADM2": 13697, "GRIN2A": 15103, "CAMK2A": 1135, "MEF2C": 1532, "APP": 8172}
for symbol, token_id in token_map.items():
    pdir = os.path.join(output_dir, f'perturb_{symbol}')
    files = os.listdir(pdir)
    with open(os.path.join(pdir, files[0]), 'rb') as f:
        data = pickle.load(f)
    key = (np.int16(token_id), 'cell_emb')
    individual_shifts[symbol] = 1.0 - np.mean(data[key])

# Synergy analysis
logger.info("\n" + "="*80)
logger.info("SYNERGY ANALYSIS")
logger.info("="*80)
for pair_name, r in results.items():
    if r.get('status') != 'success':
        continue
    g1, g2 = r['genes']
    expected = individual_shifts[g1] + individual_shifts[g2]
    observed = r['mean_shift']
    synergy = observed - expected
    synergy_pct = (synergy / expected) * 100
    label = "SYNERGISTIC" if synergy > 0 else "ANTAGONISTIC"
    logger.info(f"  {pair_name}: observed={observed:.4f}, expected={expected:.4f}, "
                f"synergy={synergy:+.4f} ({synergy_pct:+.1f}%) [{label}]")
    r['expected_additive'] = float(expected)
    r['synergy'] = float(synergy)
    r['synergy_pct'] = float(synergy_pct)

# Re-save with synergy info
with open(os.path.join(output_dir, 'combinatorial_results.json'), 'w') as f:
    json.dump(results, f, indent=2)

logger.info("\nPhase 2 complete!")
