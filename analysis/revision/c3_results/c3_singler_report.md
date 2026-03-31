# C3: SingleR Re-annotation Sensitivity Analysis

**Mode:** Raw h5ad counts
**Cells:** 11202
**Overall broad-type agreement:** 22.8%

## Per-type agreement (broad categories)

| Broad type | Agreement |
|---|---|
| astrocyte | 91.7% |
| macrophage | 70.0% |
| microglial cell | 0.0% |
| mural cell | 0.0% |
| neuron | 46.2% |
| oligodendrocyte | 0.0% |
| other | 0.0% |

## Confusion matrix (broad categories)

| marker_based    |   astrocyte |   endothelial cell |   immune |   macrophage |   mural cell |   neuron |   other |
|:----------------|------------:|-------------------:|---------:|-------------:|-------------:|---------:|--------:|
| astrocyte       |         717 |                  0 |        0 |            0 |            1 |       64 |       0 |
| macrophage      |           0 |                  0 |        1 |            7 |            0 |        0 |       2 |
| microglial cell |          11 |                  1 |       39 |          287 |            0 |        1 |      49 |
| mural cell      |           9 |                  3 |        0 |            0 |            0 |        5 |       1 |
| neuron          |        2130 |                  0 |        0 |            0 |            0 |     1833 |       5 |
| oligodendrocyte |        3046 |                 83 |        1 |            4 |            0 |     2800 |      93 |
| other           |           1 |                  0 |        8 |            0 |            0 |        0 |       0 |

## SingleR label distribution

| SingleR label | Count |
|---|---|
| Astrocyte | 5914 |
| Neurons | 4439 |
| Neuroepithelial_cell | 264 |
| Macrophage | 259 |
| Endothelial_cells | 87 |
| Monocyte | 39 |
| Gametocytes | 32 |
| Neutrophils | 30 |
| Platelets | 28 |
| DC | 24 |
| Embryonic_stem_cells | 19 |
| Erythroblast | 18 |
| HSC_CD34+ | 11 |
| HSC_-G-CSF | 10 |
| B_cell | 8 |
| NK_cell | 4 |
| T_cells | 3 |
| Pro-Myelocyte | 3 |
| MSC | 2 |
| Keratinocytes | 2 |
| Osteoblasts | 1 |
| iPS_cells | 1 |
| GMP | 1 |
| CMP | 1 |
| Smooth_muscle_cells | 1 |
| Pro-B_cell_CD34+ | 1 |

## Interpretation for manuscript

The broad-type agreement between marker-based annotation and SingleR
(HumanPrimaryCellAtlasData reference) is **22.8%**.
Agreement is moderate, suggesting that some cell-type-specific findings may be sensitive to annotation method. Detailed per-type concordance should be examined.
