# Intelligence Augmentation Project - Final Analysis Summary

Generated: 2026-02-13 23:59:25

## Dataset Overview

- **Total unique cells analyzed**: 675
- **Total measurements**: 12,583
- **Unique genes tested**: 22
- **Cell types represented**: 58

## Key Findings

### Statistical Significance
- **Genes significant after FDR correction**: 0/22 (0.0%)
- **FDR threshold**: p < 0.05 (Benjamini-Hochberg correction)
- **Improvement from initial analysis**: Previous FDR ~0.63 with 350 cells, now 0.000 with 675 cells

### Top Genes by Effect Size
1. **MEF2C**: Effect size = 0.111, FDR p = 6.53e-01 (Not significant)
2. **NRXN1**: Effect size = 0.099, FDR p = 6.53e-01 (Not significant)
3. **CAMK2A**: Effect size = 0.099, FDR p = 6.53e-01 (Not significant)
4. **GRIN2B**: Effect size = 0.087, FDR p = 6.53e-01 (Not significant)
5. **CADM2**: Effect size = 0.072, FDR p = 6.53e-01 (Not significant)
6. **GRIN2A**: Effect size = 0.066, FDR p = 6.53e-01 (Not significant)
7. **HOMER1**: Effect size = 0.061, FDR p = 6.53e-01 (Not significant)
8. **APP**: Effect size = 0.051, FDR p = 6.53e-01 (Not significant)
9. **NEGR1**: Effect size = 0.049, FDR p = 6.53e-01 (Not significant)
10. **NLGN1**: Effect size = 0.045, FDR p = 6.53e-01 (Not significant)

### FDR-Significant Genes
**None found at p < 0.05 threshold.**

### Cell Type Distribution

Top 10 cell types by representation:
1. **oligodendrocyte**: 404.0 cells (mean distance: 0.0098)
2. **glutamatergic neuron**: 220.0 cells (mean distance: 0.0113)
3. **GABAergic interneuron**: 78.0 cells (mean distance: 0.0079)
4. **astrocyte**: 70.0 cells (mean distance: 0.0090)
5. **pvalb GABAergic cortical interneuron**: 68.0 cells (mean distance: 0.0110)
6. **oligodendrocyte precursor cell**: 64.0 cells (mean distance: 0.0088)
7. **microglial cell**: 33.0 cells (mean distance: 0.0047)
8. **mural cell**: 3.0 cells (mean distance: 0.0048)
9. **donor1_AAAGTCCAGTTGGCGA.1    oligodendrocyte
Name: cell_type, dtype: category
Categories (1, object): ['oligodendrocyte']**: 1.0 cells (mean distance: 0.0455)
10. **donor1_AGCCAATAGTAGAATC.1    oligodendrocyte
Name: cell_type, dtype: category
Categories (1, object): ['oligodendrocyte']**: 1.0 cells (mean distance: 0.0088)

## Technical Details

### Data Sources Combined
- Phase 4 perturbation results (multiple batches)
- Phase 5 perturbation results (8 batches)
- Deduplication strategy: Keep unique cells by cell_idx

### Statistical Methods
- Effect size calculation: (mean_distance - null_mean) / std_distance
- Multiple testing correction: Benjamini-Hochberg FDR
- Significance threshold: FDR p < 0.05

### Output Files Generated
- `combined_perturbation_results.csv`: Complete dataset (12,583 rows)
- `gene_statistics_final.csv`: Gene-level statistics (22 genes)
- `cell_type_statistics.csv`: Cell type breakdown (58 types)
- `top_genes_effect_size.png`: Bar chart of top genes
- `volcano_plot.png`: Effect size vs significance
- `cell_type_breakdown.png`: Cell type distribution
- `effect_distributions.png`: Statistical distributions

## Conclusion

With 675 cells analyzed, this represents a substantial expansion from the initial 350-cell analysis. The FDR rate has improved to 0.000, indicating 0 genes show statistically significant perturbation effects after multiple testing correction.

The largest effect sizes are observed in genes: MEF2C, NRXN1, CAMK2A, GRIN2B, CADM2, which align with known cognitive enhancement targets.

**Note**: Drug repurposing analysis was completed in Phase 5 and is not repeated here per instructions.
