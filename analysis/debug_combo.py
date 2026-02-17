#!/usr/bin/env python3
import os, sys, traceback, torch
from multiprocessing import freeze_support

def main():
    freeze_support()
    
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    print(f"CUDA: {torch.cuda.is_available()}", flush=True)
    
    try:
        from geneformer import InSilicoPerturber
        print("Creating ISP...", flush=True)
        
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=["ENSG00000175161", "ENSG00000142192"],  # CADM2, APP
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cell",
            cell_emb_style="mean_pool",
            filter_data=None,
            max_ncells=500,
            emb_layer=-1,
            forward_batch_size=16,
            nproc=1,
            model_version="V1",
            token_dictionary_file=token_dict_path,
        )
        print(f"ISP created. combos={isp.combos}, perturb_type={isp.perturb_type}", flush=True)
        print(f"genes_to_perturb: {isp.genes_to_perturb}", flush=True)
        
        print("Running perturb_data...", flush=True)
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data,
            output_directory=os.path.join(output_dir, "perturb_combo_CADM2_APP"),
            output_prefix="combo_CADM2_APP",
        )
        print("SUCCESS!", flush=True)
    except Exception as e:
        print(f"ERROR: {e}", flush=True)
        traceback.print_exc()

if __name__ == '__main__':
    main()
