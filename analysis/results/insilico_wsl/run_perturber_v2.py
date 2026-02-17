#!/usr/bin/env python3
"""
In silico perturbation: perturb each intelligence gene individually
"""
import os, sys, pickle, json, logging, torch
import numpy as np
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

def main():
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    gene_map = {
        "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
        "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
        "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
        "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
        "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
        "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
        "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
    }
    
    logger.info(f"CUDA: {torch.cuda.is_available()}, GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A'}")
    
    from geneformer import InSilicoPerturber
    
    results = {}
    for i, (symbol, ensembl_id) in enumerate(gene_map.items()):
        logger.info(f"\n=== Perturbing {symbol} ({ensembl_id}) [{i+1}/{len(gene_map)}] ===")
        gene_output = os.path.join(output_dir, f"perturb_{symbol}")
        os.makedirs(gene_output, exist_ok=True)
        
        try:
            isp = InSilicoPerturber(
                perturb_type="delete",
                genes_to_perturb=[ensembl_id],
                model_type="Pretrained",
                num_classes=0,
                emb_mode="cell",
                filter_data=None,
                max_ncells=500,
                emb_layer=-1,
                forward_batch_size=32,
                nproc=4,
                model_version="V1",
                token_dictionary_file=token_dict_path,
            )
            
            isp.perturb_data(
                model_directory=model_path,
                input_data_file=tokenized_data,
                output_directory=gene_output,
                output_prefix=f"perturb_{symbol}",
            )
            
            # Collect results from output files
            result_files = [f for f in os.listdir(gene_output) if f.endswith('.csv') or f.endswith('.pkl')]
            results[symbol] = {"status": "success", "ensembl": ensembl_id, "files": result_files}
            logger.info(f"  {symbol}: SUCCESS, files: {result_files}")
            
        except Exception as e:
            logger.error(f"  {symbol}: FAILED - {e}")
            results[symbol] = {"status": "failed", "ensembl": ensembl_id, "error": str(e)}
    
    # Save summary
    with open(os.path.join(output_dir, "perturbation_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"\n=== DONE ===")
    successes = sum(1 for r in results.values() if r["status"] == "success")
    logger.info(f"Results: {successes}/{len(results)} genes perturbed successfully")

if __name__ == "__main__":
    main()
