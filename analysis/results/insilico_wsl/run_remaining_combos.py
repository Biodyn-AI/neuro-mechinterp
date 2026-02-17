#!/usr/bin/env python3
"""Run remaining combo perturbations sequentially, no debug logging"""
import os, sys, logging, torch, traceback

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

# Suppress all noisy loggers
for name in ['numba', 'transformers', 'datasets', 'torch', 'filelock', 'urllib3']:
    logging.getLogger(name).setLevel(logging.ERROR)

output_dir = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
model_path = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
tokenized_data = os.path.join(output_dir, "brain.dataset")
token_dict_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")

genes = {
    "CADM2": "ENSG00000175161",
    "GRIN2A": "ENSG00000183454",
    "CAMK2A": "ENSG00000070808",
    "MEF2C": "ENSG00000081189",
    "APP": "ENSG00000142192",
}

combos = [
    ("CADM2", "APP"),
    ("GRIN2A", "CAMK2A"),
    ("GRIN2A", "MEF2C"),
    ("GRIN2A", "APP"),
    ("CAMK2A", "MEF2C"),
    ("CAMK2A", "APP"),
    ("MEF2C", "APP"),
]

print(f"CUDA: {torch.cuda.is_available()}", flush=True)

from geneformer import InSilicoPerturber

for g1, g2 in combos:
    gene_output = os.path.join(output_dir, f"perturb_combo_{g1}_{g2}")
    # Check if already done
    if any(f.endswith('.pickle') for f in os.listdir(gene_output) if not os.path.isdir(os.path.join(gene_output, f))):
        print(f"SKIP {g1}+{g2}: already has pickle", flush=True)
        continue
    
    # Clean batch dirs
    for f in os.listdir(gene_output):
        fp = os.path.join(gene_output, f)
        if os.path.isdir(fp):
            import shutil
            shutil.rmtree(fp)
    
    e1, e2 = genes[g1], genes[g2]
    print(f"START {g1}+{g2}...", flush=True)
    try:
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
        isp.perturb_data(
            model_directory=model_path,
            input_data_file=tokenized_data,
            output_directory=gene_output,
            output_prefix=f"combo_{g1}_{g2}",
        )
        print(f"DONE {g1}+{g2}", flush=True)
    except Exception as e:
        print(f"FAIL {g1}+{g2}: {e}", flush=True)
        traceback.print_exc()

print("ALL COMBOS COMPLETE", flush=True)
