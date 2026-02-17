#!/usr/bin/env python3
"""Run remaining genes that weren't completed"""
import multiprocessing
multiprocessing.set_start_method('fork', force=True)
import os, sys
os.environ["TOKENIZERS_PARALLELISM"] = "false"
import pickle, json, logging, torch
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

base = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
tokenized_data = os.path.join(base, "brain.dataset")
token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")

gene_map = {
    "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760", "TCF4": "ENSG00000196628",
    "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689", "CREB1": "ENSG00000118260",
    "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056", "SCN1A": "ENSG00000144285",
    "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
}

# Skip already completed
completed = set()
for d in os.listdir(base):
    if d.startswith("perturb_") and os.listdir(os.path.join(base, d)):
        completed.add(d.replace("perturb_", ""))

remaining = {k: v for k, v in gene_map.items() if k not in completed}
logger.info(f"Remaining: {len(remaining)} genes: {list(remaining.keys())}")

from geneformer import InSilicoPerturber

for i, (symbol, ensembl_id) in enumerate(remaining.items()):
    logger.info(f"=== {symbol} ({ensembl_id}) [{i+1}/{len(remaining)}] ===")
    gene_output = os.path.join(base, f"perturb_{symbol}")
    os.makedirs(gene_output, exist_ok=True)
    try:
        isp = InSilicoPerturber(
            perturb_type="delete", genes_to_perturb=[ensembl_id],
            model_type="Pretrained", num_classes=0, emb_mode="cell",
            filter_data=None, max_ncells=500, emb_layer=-1,
            forward_batch_size=32, nproc=1, model_version="V1",
            token_dictionary_file=token_dict_path,
        )
        isp.perturb_data(model_directory=model_path, input_data_file=tokenized_data,
                         output_directory=gene_output, output_prefix=f"perturb_{symbol}")
        logger.info(f"  {symbol}: SUCCESS")
    except Exception as e:
        logger.error(f"  {symbol}: FAILED - {e}")

logger.info("=== ALL DONE ===")
