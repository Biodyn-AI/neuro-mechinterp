@echo off
echo =================================================================
echo Intelligence Augmentation - Scaled Perturbation Analysis
echo Target: 2000+ cells using Siletti_DLPFC_113k dataset
echo =================================================================

cd /d "D:\openclaw\intelligence-augmentation\analysis\results\scaled"

echo Activating conda environment 'bioinfo'...
call conda activate bioinfo

echo Starting analysis...
python scale_up_perturbation_2000plus.py

echo Analysis complete. Check logs in ./logs/ directory
pause