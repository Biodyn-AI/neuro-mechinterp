# Round 7 Robustness Checks

## 1) Residual confounding checks
- Spearman(shift, expression frequency) = 0.353, p = 0.116
- Spearman(shift, expressing-cell count) = 0.487, p = 0.0252
- Spearman(z-score vs matched null, expression frequency) = 0.116, p = 0.618
- Spearman(z-score vs matched null, expressing-cell count) = 0.256, p = 0.263

Interpretation: raw shifts remain partly linked to expression prevalence, but null-standardized z-scores show attenuated/non-significant coupling, supporting partial de-confounding by matched nulls.

## 2) Rank stability check (leave-one-gene-out)
- Top gene by null-standardized z-score: **CAMK2A**
- Top-gene identity preserved in 20/21 leave-one-out runs (95.2%)

## 3) Set-level intelligence-vs-control contrast
- Observed mean shift difference (intelligence - controls): 0.001400
- Bootstrap 95% CI: [-0.000427, 0.003167]
- Permutation p-value (two-sided): 0.1279

Interpretation: set-level contrast is modest and should be reported as exploratory unless replicated in independent cohorts.