#!/usr/bin/env python3
"""
Downstream gene analysis: For each top gene perturbation, identify which
other genes' embeddings shift the most (gene-level embedding analysis).
Uses Geneformer InSilicoPerturber with emb_mode="gene" to get gene-level shifts.
"""
import os, sys, pickle, json, logging, torch
import numpy as np
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    # Load token-to-gene mapping for interpreting results
    with open(token_dict_path, 'rb') as f:
        token_dict = pickle.load(f)
    # Reverse: token_id -> ensembl
    id_to_ensembl = {v: k for k, v in token_dict.items()}
    
    # Load gene name mapping
    gene_name_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl")
    if os.path.exists(gene_name_path):
        with open(gene_name_path, 'rb') as f:
            gene_name_dict = pickle.load(f)
        ensembl_to_name = {v: k for k, v in gene_name_dict.items()}
    else:
        ensembl_to_name = {}
    
    top5 = {
        "CADM2": "ENSG00000175161", "GRIN2A": "ENSG00000183454",
        "CAMK2A": "ENSG00000070808", "MEF2C": "ENSG00000081189",
        "APP": "ENSG00000142192",
    }
    
    logger.info(f"CUDA: {torch.cuda.is_available()}")
    from geneformer import InSilicoPerturber
    
    all_downstream = {}
    
    for symbol, ensembl_id in top5.items():
        logger.info(f"\n=== Downstream analysis for {symbol} ===")
        gene_output = os.path.join(output_dir, f"downstream_{symbol}")
        os.makedirs(gene_output, exist_ok=True)
        
        try:
            isp = InSilicoPerturber(
                perturb_type="delete",
                genes_to_perturb=[ensembl_id],
                model_type="Pretrained",
                num_classes=0,
                emb_mode="gene",
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
                output_directory=gene_output,
                output_prefix=f"downstream_{symbol}",
            )
            
            # Parse results
            pkl_files = [f for f in os.listdir(gene_output) if f.endswith('.pickle')]
            gene_shifts = {}
            for pf in pkl_files:
                with open(os.path.join(gene_output, pf), 'rb') as f:
                    data = pickle.load(f)
                for key, val in data.items():
                    token_id = key[0] if isinstance(key, tuple) else key
                    if isinstance(val, list) and len(val) > 0:
                        arr = np.array(val)
                        ensembl = id_to_ensembl.get(int(token_id), f"token_{token_id}")
                        name = ensembl_to_name.get(ensembl, ensembl)
                        gene_shifts[name] = {
                            "mean_shift": float(1 - np.mean(arr)),
                            "std": float(np.std(arr)),
                            "n_cells": len(arr),
                            "ensembl": ensembl
                        }
            
            # Sort by shift magnitude
            ranked = sorted(gene_shifts.items(), key=lambda x: abs(x[1]["mean_shift"]), reverse=True)
            top20 = ranked[:20]
            
            all_downstream[symbol] = {g: d for g, d in top20}
            logger.info(f"  Top 5 downstream of {symbol}:")
            for g, d in top20[:5]:
                logger.info(f"    {g}: shift={d['mean_shift']:.6f}")
            
        except Exception as e:
            logger.error(f"  {symbol}: FAILED - {e}")
            import traceback
            traceback.print_exc()
            all_downstream[symbol] = {"error": str(e)}
    
    with open(os.path.join(output_dir, "downstream_results.json"), "w") as f:
        json.dump(all_downstream, f, indent=2)
    logger.info("\nSaved downstream_results.json")
