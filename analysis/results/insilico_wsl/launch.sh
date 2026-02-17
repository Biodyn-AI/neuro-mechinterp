#!/bin/bash
cd /mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl
# brain.dataset already tokenized
python3 run_perturber.py > /mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/stdout.log 2>&1
echo "EXIT CODE: $?" >> /mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/stdout.log
