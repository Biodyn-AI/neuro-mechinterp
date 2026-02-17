#!/usr/bin/env python3
"""
In silico perturbation experiment with Geneformer
Testing impact of 22 intelligence genes on brain cell embeddings
"""

import os
import sys
import time
import logging
import pickle
import torch
import numpy as np
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('/tmp/perturbation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def main():
    logger.info("Starting in silico perturbation experiment...")
    
    # Check CUDA availability
    cuda_available = torch.cuda.is_available()
    device_info = f"CUDA available: {cuda_available}"
    if cuda_available:
        device_info += f", GPU: {torch.cuda.get_device_name(0)}"
    logger.info(device_info)
    
    # Define paths
    brain_data_path = "/mnt/d/openclaw/intelligence-augmentation/data/brain_scrna/DLPFC_11k.h5ad"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    tokenized_data_path = os.path.join(output_dir, "brain.dataset")
    
    # Intelligence genes - symbol to Ensembl mapping
    gene_symbol_to_ensembl = {
        "MEF2C": "ENSG00000081189", "BDNF": "ENSG00000049759", "GRIN2B": "ENSG00000273079",
        "CADM2": "ENSG00000175161", "NRXN1": "ENSG00000179915", "CAMK2A": "ENSG00000070808",
        "GRIN2A": "ENSG00000183454", "SHANK3": "ENSG00000251322", "HOMER1": "ENSG00000152413",
        "APP": "ENSG00000142192", "NEGR1": "ENSG00000172260", "NLGN1": "ENSG00000169760",
        "TCF4": "ENSG00000196628", "MAPT": "ENSG00000186868", "FOXO3": "ENSG00000118689",
        "CREB1": "ENSG00000118260", "FMR1": "ENSG00000102081", "SYN1": "ENSG00000008056",
        "SCN1A": "ENSG00000144285", "SLC6A4": "ENSG00000108576", "COMT": "ENSG00000093010",
    }
    
    # Load token dictionary to filter to genes that exist in it
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    with open(token_dict_path, "rb") as f:
        token_dict = pickle.load(f)
    
    # Filter to Ensembl IDs present in token dictionary
    intelligence_genes = []
    for sym, ens in gene_symbol_to_ensembl.items():
        if ens in token_dict:
            intelligence_genes.append(ens)
            logger.info(f"  {sym} -> {ens} (token {token_dict[ens]})")
        else:
            logger.warning(f"  {sym} -> {ens} NOT in token dictionary, skipping")
    
    logger.info(f"Intelligence genes to perturb: {len(intelligence_genes)} of {len(gene_symbol_to_ensembl)} genes")
    logger.info(f"Genes: {', '.join(intelligence_genes)}")
    
    try:
        # Step 1: Check paths exist
        if not os.path.exists(brain_data_path):
            raise FileNotFoundError(f"Brain data not found: {brain_data_path}")
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model not found: {model_path}")
        
        logger.info("Paths verified successfully")
        
        # Step 2: Tokenize the brain data if not already done
        if not os.path.exists(tokenized_data_path):
            logger.info("Tokenizing brain data...")
            from geneformer import TranscriptomeTokenizer
            
            tk = TranscriptomeTokenizer(
                custom_attr_name_dict={"cell_type": "cell_type"},
                nproc=4,
                model_version="V1",
            )
            
            # Simplified tokenization call - removed unsupported parameters
            tk.tokenize_data(
                data_directory="/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/prepared_data",
                output_directory=output_dir,
                output_prefix="brain",
                file_format="h5ad",
                use_generator=False
            )
            logger.info("Data tokenization completed!")
        else:
            logger.info("Using existing tokenized data...")
        
        # Step 3: Set up InSilicoPerturber
        logger.info("Setting up InSilicoPerturber...")
        from geneformer import InSilicoPerturber
        
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=intelligence_genes,
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cls",
            filter_data=None,
            max_ncells=500,  # Start with 500 cells for faster testing
            emb_layer=-1,
            forward_batch_size=16,  # Smaller batch for memory safety
            nproc=4,
            model_version="V1",
            token_dictionary_file=os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
        )
        logger.info("InSilicoPerturber initialized successfully!")
        
        # Step 4: Run perturbation analysis
        logger.info("Starting perturbation analysis...")
        logger.info(f"Model path: {model_path}")
        logger.info(f"Input data: {tokenized_data_path}")
        logger.info(f"Output directory: {output_dir}")
        
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data_path,
            output_directory=output_dir,
            output_prefix="intelligence_genes"
        )
        
        logger.info("Perturbation analysis completed successfully!")
        
        # Step 5: Analyze results
        logger.info("Analyzing perturbation results...")
        result_files = list(Path(output_dir).glob("intelligence_genes*.pickle"))
        logger.info(f"Found {len(result_files)} result files")
        
        for file_path in result_files:
            try:
                size_mb = file_path.stat().st_size / (1024*1024)
                logger.info(f"  - {file_path.name} ({size_mb:.2f} MB)")
                
                # Try to examine the structure
                with open(file_path, 'rb') as f:
                    data = pickle.load(f)
                if isinstance(data, dict):
                    logger.info(f"    Keys: {list(data.keys())}")
                else:
                    logger.info(f"    Type: {type(data)}")
            except Exception as e:
                logger.warning(f"    Could not examine file {file_path.name}: {e}")
        
        # Step 6: Save summary
        summary_file = os.path.join(output_dir, "perturbation_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("In Silico Perturbation Experiment Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Genes tested: {len(intelligence_genes)}\n")
            f.write(f"Gene list: {', '.join(intelligence_genes)}\n")
            f.write(f"Result files: {len(result_files)}\n")
            f.write(f"Model: {model_path}\n")
            f.write(f"Data: {brain_data_path}\n")
            f.write(f"CUDA available: {cuda_available}\n")
            f.write(f"Max cells processed: 500\n")
        
        logger.info(f"Summary saved to: {summary_file}")
        logger.info("Experiment completed successfully!")
        
    except Exception as e:
        logger.error(f"Error during perturbation experiment: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()