#!/usr/bin/env python3
"""Run perturbation with emb_mode='cell_and_gene' for downstream gene analysis.
Usage: python run_gene_level.py GENE_SYMBOL"""
import os, sys, pickle, torch
import numpy as np
from multiprocessing import freeze_support

def main():
    freeze_support()
    
    symbol = sys.argv[1]
    gene_ensembl = {
        "CADM2": "ENSG00000175161", "GRIN2A": "ENSG00000183454",
        "CAMK2A": "ENSG00000070808", "MEF2C": "ENSG00000081189",
        "APP": "ENSG00000142192",
    }
    
    output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
    model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
    tokenized_data = os.path.join(output_dir, "brain.dataset")
    token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
    
    gene_output = os.path.join(output_dir, f"perturb_genelevel_{symbol}")
    os.makedirs(gene_output, exist_ok=True)
    
    print(f"Gene-level perturbation for {symbol}...")
    print(f"CUDA: {torch.cuda.is_available()}")
    
    from geneformer import InSilicoPerturber
    
    isp = InSilicoPerturber(
        perturb_type="delete",
        genes_to_perturb=[gene_ensembl[symbol]],
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell_and_gene",
        cell_emb_style="mean_pool",
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
        output_prefix=f"genelevel_{symbol}",
    )
    
    pkl_files = [f for f in os.listdir(gene_output) if f.endswith('.pickle')]
    print(f"Output files: {pkl_files}")
    
    for pf in pkl_files:
        with open(os.path.join(gene_output, pf), 'rb') as f:
            data = pickle.load(f)
        print(f"\n{pf}:")
        print(f"  Type: {type(data)}")
        if isinstance(data, dict):
            print(f"  N keys: {len(data)}")
            for k in list(data.keys())[:5]:
                v = data[k]
                if isinstance(v, list):
                    print(f"  {k}: list[{len(v)}]")
                elif isinstance(v, np.ndarray):
                    print(f"  {k}: array{v.shape}")
                else:
                    print(f"  {k}: {type(v)}")
    
    print("DONE")

if __name__ == '__main__':
    main()
