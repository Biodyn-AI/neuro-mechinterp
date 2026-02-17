exec(open('/mnt/c/Users/Agent/.openclaw/workspace/patch_torchtext.py').read())
#!/usr/bin/env python3
"""
Run Geneformer in silico perturbation on expression-matched control genes.
Uses identical pipeline to the intelligence gene analysis.
V2: Fixed version that works with current Geneformer install.
"""
import os, sys, json, time, logging, pickle
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

os.environ['TOKENIZERS_PARALLELISM'] = 'false'
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

# Paths - same as original intelligence gene analysis
OUTPUT_BASE = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
MODEL_PATH = "/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M"
TOKENIZED_DATA = os.path.join(OUTPUT_BASE, "brain.dataset")

# Control gene map: symbol -> ensembl ID
# 12 expression-matched controls + 5 housekeeping controls = 17 total
CONTROL_GENE_MAP = {
    # Expression-matched controls
    "EXOC6B": "ENSG00000144036",    # matched to CADM2
    "DOCK9": "ENSG00000088387",     # matched to APP
    "PDE7B": "ENSG00000171408",     # matched to FOXO3
    "ZDHHC9": "ENSG00000188706",    # matched to NLGN1
    "HDAC11": "ENSG00000163517",    # matched to MAPT
    "XRCC6": "ENSG00000196419",     # matched to TCF4
    "WDR46": "ENSG00000227057",     # matched to CREB1
    "CUEDC2": "ENSG00000107874",    # matched to DISC1
    "IGSF10": "ENSG00000152580",    # matched to NEGR1
    "SALL3": "ENSG00000256463",     # matched to FMR1
    "MFSD6L": "ENSG00000185156",    # matched to COMT
    "TMEM88": "ENSG00000167874",    # matched to MEF2C
    # Housekeeping controls
    "ACTB": "ENSG00000075624",
    "GAPDH": "ENSG00000111640",
    "B2M": "ENSG00000166710",
    "RPL13A": "ENSG00000142541",
    "YWHAZ": "ENSG00000164924",
}

# Which intelligence gene each control is matched to
MATCHED_TO = {
    "EXOC6B": "CADM2", "DOCK9": "APP", "PDE7B": "FOXO3", "ZDHHC9": "NLGN1",
    "HDAC11": "MAPT", "XRCC6": "TCF4", "WDR46": "CREB1", "CUEDC2": "DISC1",
    "IGSF10": "NEGR1", "SALL3": "FMR1", "MFSD6L": "COMT", "TMEM88": "MEF2C",
    "ACTB": "housekeeping", "GAPDH": "housekeeping", "B2M": "housekeeping",
    "RPL13A": "housekeeping", "YWHAZ": "housekeeping",
}

def check_existing_result(gene_output):
    """Check if a valid result pickle already exists."""
    if not os.path.exists(gene_output):
        return False
    pickles = [f for f in os.listdir(gene_output) if f.endswith('.pickle')]
    if not pickles:
        return False
    # Verify the pickle is valid by loading it
    try:
        fpath = os.path.join(gene_output, pickles[0])
        with open(fpath, 'rb') as f:
            data = pickle.load(f)
        # Should be a defaultdict with (token_id, 'cell_emb') keys
        if len(data) > 0:
            key = list(data.keys())[0]
            if isinstance(key, tuple) and 'cell_emb' in str(key):
                vals = list(data.values())[0]
                if isinstance(vals, list) and len(vals) > 0:
                    return True
    except Exception:
        pass
    return False

def main():
    import torch
    logger.info(f"CUDA: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")

    from geneformer import InSilicoPerturber

    results = {}
    total = len(CONTROL_GENE_MAP)
    start_time = time.time()

    for i, (symbol, ensembl_id) in enumerate(CONTROL_GENE_MAP.items()):
        logger.info(f"\n=== Perturbing CONTROL {symbol} ({ensembl_id}) [{i+1}/{total}] ===")
        logger.info(f"    Matched to: {MATCHED_TO[symbol]}")

        gene_output = os.path.join(OUTPUT_BASE, f"control_perturb_{symbol}")
        os.makedirs(gene_output, exist_ok=True)

        # Skip if already has valid results
        if check_existing_result(gene_output):
            existing = [f for f in os.listdir(gene_output) if f.endswith('.pickle')]
            logger.info(f"  SKIP {symbol} - already has valid results: {existing}")
            results[symbol] = {"status": "success", "ensembl": ensembl_id,
                               "matched_to": MATCHED_TO[symbol], "files": existing,
                               "skipped": True}
            continue

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
                nproc=1,
                model_version="V1",
            )

            isp.perturb_data(
                model_directory=MODEL_PATH,
                input_data_file=TOKENIZED_DATA,
                output_directory=gene_output,
                output_prefix=f"control_perturb_{symbol}",
            )

            result_files = [f for f in os.listdir(gene_output) if f.endswith('.pickle')]
            results[symbol] = {"status": "success", "ensembl": ensembl_id,
                               "matched_to": MATCHED_TO[symbol], "files": result_files}
            logger.info(f"  {symbol}: SUCCESS, files: {result_files}")

        except Exception as e:
            import traceback
            logger.error(f"  {symbol}: FAILED - {e}")
            logger.error(traceback.format_exc())
            results[symbol] = {"status": "failed", "ensembl": ensembl_id,
                               "matched_to": MATCHED_TO[symbol], "error": str(e)}

        # Clear GPU cache
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    elapsed = (time.time() - start_time) / 60
    logger.info(f"\n=== CONTROL PERTURBATION COMPLETE ===")
    successes = sum(1 for r in results.values() if r["status"] == "success")
    logger.info(f"Results: {successes}/{total} control genes perturbed successfully")
    logger.info(f"Elapsed: {elapsed:.1f} minutes")

    # Save summary
    summary_file = os.path.join(OUTPUT_BASE, "control_perturbation_summary_v2.json")
    with open(summary_file, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Summary saved to {summary_file}")

if __name__ == "__main__":
    main()
