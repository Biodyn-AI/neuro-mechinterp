#!/bin/bash
export TOKENIZERS_PARALLELISM=false
cd /mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl

for gene_line in "CREB1 ENSG00000118260" "FMR1 ENSG00000102081" "SYN1 ENSG00000008056" "SCN1A ENSG00000144285" "SLC6A4 ENSG00000108576" "COMT ENSG00000093010"; do
    set -- $gene_line
    symbol=$1
    ensembl=$2
    if [ -d "perturb_${symbol}" ] && [ "$(ls perturb_${symbol}/*.pickle 2>/dev/null | wc -l)" -gt 0 ]; then
        echo "SKIP $symbol"
        continue
    fi
    echo "=== $symbol ==="
    python3 -c "
import multiprocessing
multiprocessing.set_start_method('fork', force=True)
import os
os.environ['TOKENIZERS_PARALLELISM'] = 'false'
from geneformer import InSilicoPerturber
isp = InSilicoPerturber(
    perturb_type='delete', genes_to_perturb=['$ensembl'],
    model_type='Pretrained', num_classes=0, emb_mode='cell',
    filter_data=None, max_ncells=500, emb_layer=-1,
    forward_batch_size=32, nproc=1, model_version='V1',
    token_dictionary_file=os.path.expanduser('~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl'),
)
os.makedirs('perturb_$symbol', exist_ok=True)
isp.perturb_data(
    model_directory='/mnt/d/openclaw/intelligence-augmentation/models/Geneformer/Geneformer-V1-10M',
    input_data_file='brain.dataset',
    output_directory='perturb_$symbol',
    output_prefix='perturb_$symbol',
)
print('SUCCESS: $symbol')
" 2>&1
done
echo "=== ALL 6 COMPLETE ==="
