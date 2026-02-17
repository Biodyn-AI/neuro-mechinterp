#!/bin/bash
cd /mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl

pairs=(
  "GRIN2A CAMK2A"
  "GRIN2A MEF2C"
  "GRIN2A APP"
  "CAMK2A MEF2C"
  "CAMK2A APP"
  "MEF2C APP"
  "CADM2 APP"
)

for pair in "${pairs[@]}"; do
  g1=$(echo $pair | cut -d' ' -f1)
  g2=$(echo $pair | cut -d' ' -f2)
  dir="perturb_combo_${g1}_${g2}"
  # Check if already has pickle files
  if ls ${dir}/*.pickle 2>/dev/null | head -1 | grep -q pickle; then
    echo "SKIP: ${g1}_${g2} (already done)"
    continue
  fi
  echo "RUNNING: ${g1}_${g2}"
  python3 run_combo_single.py $g1 $g2 2>&1
  echo "DONE: ${g1}_${g2}"
done
echo "=== ALL COMPLETE ==="
