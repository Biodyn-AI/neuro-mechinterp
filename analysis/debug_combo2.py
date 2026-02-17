#!/usr/bin/env python3
import os, sys, traceback, torch, gc
from multiprocessing import freeze_support

def main():
    freeze_support()
    
    g1, g2 = sys.argv[1], sys.argv[2]
    gene_ensembl = {
        "CADM2": "ENSG00000175161", "GRIN2A": "ENSG00000183454",
        "CAMK2A": "ENSG00000070808", "MEF2C": "ENSG00000081189",
        "APP": "ENSG00000142192",
    }
    
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    pair_name = f"{g1}_{g2}"
    pair_output = os.path.join(output_dir, f"perturb_combo_{pair_name}")
    os.makedirs(pair_output, exist_ok=True)
    
    print(f"Perturbing {g1} + {g2}...", flush=True)
    torch.cuda.empty_cache()
    gc.collect()
    
    try:
        from geneformer import InSilicoPerturber
        import numpy as np, pickle
        
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[gene_ensembl[g1], gene_ensembl[g2]],
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cell",
            cell_emb_style="mean_pool",
            filter_data=None,
            max_ncells=500,
            emb_layer=-1,
            forward_batch_size=8,
            nproc=1,
            model_version="V1",
            token_dictionary_file=token_dict_path,
            clear_mem_ncells=50,
        )
        
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data,
            output_directory=pair_output,
            output_prefix=f"combo_{pair_name}",
        )
        
        pkl_files = [f for f in os.listdir(pair_output) if f.endswith('.pickle')]
        print(f"Output files: {pkl_files}", flush=True)
        for pf in pkl_files:
            with open(os.path.join(pair_output, pf), 'rb') as f:
                data = pickle.load(f)
            for key in data:
                if 'cell_emb' in str(key):
                    cos_sims = data[key]
                    shift = 1.0 - np.mean(cos_sims)
                    print(f"RESULT: {pair_name} shift={shift:.6f} n={len(cos_sims)}", flush=True)
        print("DONE", flush=True)
    except Exception as e:
        print(f"ERROR: {e}", flush=True)
        traceback.print_exc()

if __name__ == '__main__':
    main()
