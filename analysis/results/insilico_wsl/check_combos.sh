#!/bin/bash
for d in /mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/perturb_combo_*/; do
    name=$(basename "$d")
    pkl=$(ls "$d"/*.pickle 2>/dev/null | wc -l)
    echo "$name: $pkl pickles"
done
