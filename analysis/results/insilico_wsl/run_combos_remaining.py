#!/usr/bin/env python3
"""Run remaining combinatorial perturbations for top 5 genes pairwise."""
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

# Check which combos already have results
all_pairs = list(combinations(top5.keys(), 2))
remaining = []
for g1, g2 in all_pairs:
    combo_dir = os.path.join(output_dir, f"perturb_combo_{g1}_{g2}")
    pickles = [f for f in os.listdir(combo_dir) if f.endswith('.pickle')] if os.path.isdir(combo_dir) else []
    if not pickles:
        # Also check reverse order
        combo_dir2 = os.path.join(output_dir, f"perturb_combo_{g2}_{g1}")
        pickles2 = [f for f in os.listdir(combo_dir2) if f.endswith('.pickle')] if os.path.isdir(combo_dir2) else []
        if not pickles2:
            remaining.append((g1, g2))

logger.info(f"Total pairs: {len(all_pairs)}, remaining: {len(remaining)}")
logger.info(f"Remaining: {remaining}")
logger.info(f"CUDA: {torch.cuda.is_available()}")

from geneformer import InSilicoPerturber

for g1, g2 in remaining:
    combo_name = f"{g1}_{g2}"
    logger.info(f"\n=== Combo: {combo_name} ===")
    combo_output = os.path.join(output_dir, f"perturb_combo_{combo_name}")
    os.makedirs(combo_output, exist_ok=True)
    
    try:
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[top5[g1], top5[g2]],
            combos=1,  # Enable combinatorial mode
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
            output_directory=combo_output,
            output_prefix=f"combo_{combo_name}",
        )
        logger.info(f"  {combo_name}: SUCCESS")
    except Exception as e:
        logger.error(f"  {combo_name}: FAILED - {e}")
        import traceback; traceback.print_exc()

logger.info("=== ALL COMBOS DONE ===")
