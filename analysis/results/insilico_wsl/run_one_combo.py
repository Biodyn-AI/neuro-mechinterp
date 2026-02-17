#!/usr/bin/env python3
"""Run a single combo perturbation with verbose error handling"""
import os, sys, pickle, json, logging, torch, traceback
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    g1, g2 = "CADM2", "APP"
    e1, e2 = "ENSG00000175161", "ENSG00000142192"
    
    logger.info(f"CUDA: {torch.cuda.is_available()}")
    logger.info("Importing InSilicoPerturber...")
    
    try:
        from geneformer import InSilicoPerturber
        logger.info("Import OK")
        
        gene_output = os.path.join(output_dir, f"perturb_combo_{g1}_{g2}")
        os.makedirs(gene_output, exist_ok=True)
        
        # Clear old files
        for f in os.listdir(gene_output):
            if f.endswith('.pickle'):
                os.remove(os.path.join(gene_output, f))
        
        logger.info(f"Creating ISP for {g1}+{g2}...")
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[e1, e2],
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
        logger.info("ISP created, running perturb_data...")
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data,
            output_directory=gene_output,
            output_prefix=f"combo_{g1}_{g2}",
        )
        logger.info("SUCCESS!")
        # Check output
        files = os.listdir(gene_output)
        logger.info(f"Output files: {files}")
    except Exception as e:
        logger.error(f"FAILED: {e}")
        traceback.print_exc()
        sys.exit(1)
