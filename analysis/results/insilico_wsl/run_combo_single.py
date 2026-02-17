#!/usr/bin/env python3
"""Run a single combo perturbation. Usage: python run_combo_single.py GENE1 GENE2"""
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

import os, sys, pickle, json, logging, torch
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

def main():
    g1, g2 = sys.argv[1], sys.argv[2]
    
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    top5 = {
        "CADM2": "ENSG00000175161", "GRIN2A": "ENSG00000183454",
        "CAMK2A": "ENSG00000070808", "MEF2C": "ENSG00000081189",
        "APP": "ENSG00000142192",
    }
    
    combo_name = f"{g1}_{g2}"
    combo_output = os.path.join(output_dir, f"perturb_combo_{combo_name}")
    os.makedirs(combo_output, exist_ok=True)
    
    logger.info(f"CUDA: {torch.cuda.is_available()}")
    logger.info(f"Running combo: {combo_name}")
    
    from geneformer import InSilicoPerturber
    
    isp = InSilicoPerturber(
        perturb_type="delete",
        genes_to_perturb=[top5[g1], top5[g2]],
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

if __name__ == '__main__':
    main()
