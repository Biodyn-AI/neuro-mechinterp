#!/usr/bin/env python3
"""
Phase 3 Completion Verification Script
Checks all outputs and provides final summary
"""

import os
import pandas as pd
from pathlib import Path

def verify_phase3_completion():
    """Verify all Phase 3 outputs are complete and valid"""
    
    phase3_dir = Path("D:/openclaw/intelligence-augmentation/analysis/results/phase3")
    
    print("=== PHASE 3 COMPLETION VERIFICATION ===\n")
    
    # Check required files
    required_files = {
        'Figures (PNG)': [
            'Figure1_Heatmap.png',
            'Figure2_RankedBars.png', 
            'Figure3_Volcano.png',
            'Figure4_DrugNetwork.png',
            'SupplementaryFigure1_ExpressionScatter.png',
            'SupplementaryFigure2_CellTypeProfiles.png'
        ],
        'Figures (PDF)': [
            'Figure1_Heatmap.pdf',
            'Figure2_RankedBars.pdf',
            'Figure3_Volcano.pdf', 
            'Figure4_DrugNetwork.pdf',
            'SupplementaryFigure1_ExpressionScatter.pdf',
            'SupplementaryFigure2_CellTypeProfiles.pdf'
        ],
        'Statistical Data': [
            'gene_statistics.csv',
            'gene_celltype_statistics.csv',
            'correlation_statistics.csv'
        ],
        'Publication Materials': [
            'results_table.tex',
            'PHASE3_SUMMARY_REPORT.md'
        ],
        'Source Code': [
            'phase3_analysis.py'
        ]
    }
    
    all_complete = True
    
    for category, files in required_files.items():
        print(f"[+] {category}:")
        for file in files:
            file_path = phase3_dir / file
            if file_path.exists():
                size = file_path.stat().st_size
                print(f"  [OK] {file} ({size:,} bytes)")
            else:
                print(f"  [MISSING] {file}")
                all_complete = False
        print()
    
    # Verify data file contents
    print("[DATA] Data Verification:")
    
    try:
        gene_stats = pd.read_csv(phase3_dir / 'gene_statistics.csv')
        print(f"  [OK] gene_statistics.csv: {len(gene_stats)} genes")
        
        gene_celltype_stats = pd.read_csv(phase3_dir / 'gene_celltype_statistics.csv') 
        print(f"  [OK] gene_celltype_statistics.csv: {len(gene_celltype_stats)} gene×celltype pairs")
        
        corr_stats = pd.read_csv(phase3_dir / 'correlation_statistics.csv')
        print(f"  [OK] correlation_statistics.csv: {len(corr_stats)} correlations")
        
        # Key statistics
        print(f"\n[RESULTS] Key Results:")
        top_gene = gene_stats.loc[gene_stats['mean_score'].idxmax()]
        print(f"  Top gene: {top_gene['gene']} (score: {top_gene['mean_score']:.6f})")
        
        significant_genes = gene_stats[gene_stats['p_value_fdr'] < 0.05]
        print(f"  Significant genes (FDR < 0.05): {len(significant_genes)}")
        
        max_corr = corr_stats.loc[corr_stats['correlation'].idxmax()]
        print(f"  Strongest correlation: {max_corr['gene']} (r = {max_corr['correlation']:.3f})")
        
    except Exception as e:
        print(f"  [ERROR] Data verification failed: {e}")
        all_complete = False
    
    print(f"\n{'='*50}")
    if all_complete:
        print("[SUCCESS] PHASE 3 ANALYSIS COMPLETE - ALL FILES GENERATED SUCCESSFULLY")
        print("[READY] FOR PUBLICATION SUBMISSION")
    else:
        print("[WARNING] SOME FILES ARE MISSING - CHECK ERRORS ABOVE")
    print(f"{'='*50}")
    
    return all_complete

if __name__ == "__main__":
    verify_phase3_completion()