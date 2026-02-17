#!/bin/bash
cd /mnt/d/openclaw/intelligence-augmentation/analysis

pairs=(
    "CADM2 CAMK2A"
    "CADM2 MEF2C"
    "CADM2 APP"
    "GRIN2A CAMK2A"
    "GRIN2A MEF2C"
    "GRIN2A APP"
    "CAMK2A MEF2C"
    "CAMK2A APP"
    "MEF2C APP"
)

for pair in "${pairs[@]}"; do
    echo "=== Running: $pair ==="
    python3 run_combo_single.py $pair 2>&1
    echo "=== Done: $pair ==="
    echo ""
done

echo "ALL COMBINATIONS COMPLETE"
