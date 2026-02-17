exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
"""
Run Geneformer delete-perturbation on a SINGLE control gene.
Called as: python permutation_null_v2_run_one.py <ensembl_id> <symbol> <output_json_path>

This script is meant to be called in a fresh subprocess for each gene,
so that memory is fully reclaimed between genes and no broken-pipe issues cascade.
"""
import os, sys, json, pickle, time, gc
import numpy as np

os.environ['TOKENIZERS_PARALLELISM'] = 'false'
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'

OUTPUT_BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
MODEL_PATH = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
TOKENIZED_DATA = os.path.join(OUTPUT_BASE, "brain.dataset")
PERM_OUTPUT_DIR = os.path.join(OUTPUT_BASE, "permutation_null_v2")


def load_perturbation_result(result_dir):
    """Load perturbation result from a directory. Returns cosine similarities list."""
    if not os.path.exists(result_dir):
        return None
    pickles = [f for f in os.listdir(result_dir) if f.endswith('.pickle')]
    if not pickles:
        return None
    fpath = os.path.join(result_dir, pickles[0])
    with open(fpath, 'rb') as f:
        data = pickle.load(f)
    if len(data) == 0:
        return None
    key = list(data.keys())[0]
    return data[key]


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <ensembl_id> <symbol> <output_json_path>")
        sys.exit(1)

    ensembl_id = sys.argv[1]
    symbol = sys.argv[2]
    output_json = sys.argv[3]

    print(f"=== Perturbing {symbol} ({ensembl_id}) ===")
    t0 = time.time()

    gene_output = os.path.join(PERM_OUTPUT_DIR, f"perm_ctrl_{symbol}")
    os.makedirs(gene_output, exist_ok=True)

    # Check if already done
    existing = load_perturbation_result(gene_output)
    if existing is not None:
        shifts = [1.0 - c for c in existing]
        result = {
            'ensembl': ensembl_id,
            'symbol': symbol,
            'mean_shift': float(np.mean(shifts)),
            'std_shift': float(np.std(shifts)),
            'n_cells': len(shifts),
            'status': 'reused',
            'time_seconds': 0.0,
        }
        with open(output_json, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"  REUSED: shift={np.mean(shifts):.6f}, n={len(shifts)}")
        return

    try:
        import torch
        print(f"  CUDA: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"  GPU: {torch.cuda.get_device_name(0)}")
            torch.cuda.empty_cache()

        from geneformer import InSilicoPerturber

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
            nproc=1,
            model_version="V1",
        )

        isp.perturb_data(
            model_directory=MODEL_PATH,
            input_data_file=TOKENIZED_DATA,
            output_directory=gene_output,
            output_prefix=f"perm_ctrl_{symbol}",
        )

        elapsed = time.time() - t0

        # Load result
        cos_sims = load_perturbation_result(gene_output)
        if cos_sims is not None:
            shifts = [1.0 - c for c in cos_sims]
            result = {
                'ensembl': ensembl_id,
                'symbol': symbol,
                'mean_shift': float(np.mean(shifts)),
                'std_shift': float(np.std(shifts)),
                'n_cells': len(shifts),
                'status': 'success',
                'time_seconds': elapsed,
            }
            print(f"  SUCCESS: shift={np.mean(shifts):.6f}, n={len(shifts)}, time={elapsed:.1f}s")
        else:
            result = {
                'ensembl': ensembl_id,
                'symbol': symbol,
                'status': 'no_result',
                'time_seconds': elapsed,
            }
            print(f"  WARNING: No result file produced")

    except Exception as e:
        import traceback
        elapsed = time.time() - t0
        result = {
            'ensembl': ensembl_id,
            'symbol': symbol,
            'status': 'failed',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'time_seconds': elapsed,
        }
        print(f"  FAILED: {e}")
        traceback.print_exc()

    with open(output_json, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Result saved to {output_json}")


if __name__ == "__main__":
    main()
