#!/usr/bin/env python3
"""Run combo perturbation with reduced memory. Usage: python run_combo_v3.py GENE1 GENE2"""
import os, sys, pickle, torch, gc
import numpy as np
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
    
    from geneformer import InSilicoPerturber
    
    # Run in two batches of 250 cells to avoid OOM
    all_cos_sims = []
    for batch_idx, (start, end) in enumerate([(0, 100), (100, 200), (200, 300), (300, 400), (400, 500)]):
        print(f"  Batch {batch_idx+1}/5 (cells {start}-{end})...", flush=True)
        torch.cuda.empty_cache()
        gc.collect()
        
        batch_output = os.path.join(pair_output, f"batch_{batch_idx}")
        os.makedirs(batch_output, exist_ok=True)
        
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[gene_ensembl[g1], gene_ensembl[g2]],
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cell",
            cell_emb_style="mean_pool",
            filter_data=None,
            max_ncells=None,
            cell_inds_to_perturb={"start": start, "end": end},
            emb_layer=-1,
            forward_batch_size=16,
            nproc=1,
            model_version="V1",
            token_dictionary_file=token_dict_path,
        )
        
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data,
            output_directory=batch_output,
            output_prefix=f"combo_{pair_name}_b{batch_idx}",
        )
        
        # Load and collect results
        pkl_files = [f for f in os.listdir(batch_output) if f.endswith('.pickle')]
        for pf in pkl_files:
            with open(os.path.join(batch_output, pf), 'rb') as f:
                data = pickle.load(f)
            for key in data:
                if 'cell_emb' in str(key):
                    all_cos_sims.extend(data[key])
        
        del isp
        torch.cuda.empty_cache()
        gc.collect()
    
    if all_cos_sims:
        shift = 1.0 - np.mean(all_cos_sims)
        print(f"RESULT: {pair_name} shift={shift:.6f} n={len(all_cos_sims)}", flush=True)
        
        # Save combined result in same format
        combined_data = {((0, 0), 'cell_emb'): all_cos_sims}
        with open(os.path.join(pair_output, f"combined_{pair_name}.pickle"), 'wb') as f:
            pickle.dump(combined_data, f)
    
    print("DONE", flush=True)

if __name__ == '__main__':
    main()
