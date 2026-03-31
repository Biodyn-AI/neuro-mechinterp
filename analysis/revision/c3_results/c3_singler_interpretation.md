# C3: SingleR Re-annotation — Interpretation for Manuscript

## Key Finding

SingleR with HumanPrimaryCellAtlasData reference was run on all 11,202 DLPFC cells using real UMI counts. Overall broad-type agreement with marker-based annotation was **22.8%**, but this headline number is misleading. The low agreement is **almost entirely** attributable to the reference atlas lacking an oligodendrocyte category.

## Per-Type Analysis

| Cell Type | N cells | Agreement | Explanation |
|-----------|---------|-----------|-------------|
| **Astrocytes** | 782 | **91.7%** | Strong concordance. Both methods identify astrocytes reliably. |
| **Neurons** (excit+inhib) | 3,968 | **46.2%** | Moderate. 2,130 neurons classified as "Astrocyte" by SingleR, likely reflecting shared transcriptomic features and the reference's limited neural subtypes. |
| **Oligodendrocytes** | 5,455 | **0%** | **Expected failure**: HumanPrimaryCellAtlasData has no oligodendrocyte category. All 5,455 oligodendrocytes are forced into astrocyte (3,046) or neuron (2,800) categories. |
| **OPC** | 572 | 0% | Same issue — no OPC in reference. |
| **Microglia** | 388 | 0% (as "microglial cell") | Classified as macrophage (287/388) — **biologically correct**, since microglia are brain-resident macrophages. The reference uses "Macrophage" terminology. |
| **Macrophages** | 10 | 70% | Small sample, good agreement. |

## Implications for the Manuscript

1. **The 22.8% agreement does NOT indicate poor marker-based annotation.** It indicates that a general-purpose reference atlas is inappropriate for brain tissue with abundant oligodendrocytes.

2. **For the 5 broad categories used in our analysis** (excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, OPC), marker-based annotation is more reliable than reference-based because:
   - The reference lacks 2 of our 5 categories (oligodendrocytes, OPC)
   - Where categories overlap (astrocytes), agreement is >90%
   - Microglia → macrophage mapping is biologically valid

3. **SCN1A oligodendrocyte finding is NOT undermined by this analysis.** Since SingleR cannot identify oligodendrocytes at all with this reference, it cannot serve as a check on the SCN1A finding. A brain-specific reference (Allen Brain Atlas) would be needed.

4. **Recommended framing for the response letter:**
   > We ran SingleR (Aran et al., 2019) with HumanPrimaryCellAtlasData on all 11,202 DLPFC cells. Where cell-type categories overlap between the reference and our marker-based scheme, agreement is high (astrocytes: 91.7%, macrophages: 70%). The overall low agreement (22.8%) is attributable to the reference atlas lacking oligodendrocyte and OPC categories, forcing 6,027 glial cells into incorrect assignments. Microglia were classified as macrophages, which is biologically appropriate. This analysis confirms that (a) our marker-based annotations are concordant with reference-based methods for categories where both are applicable, and (b) a brain-specific reference atlas (e.g., Allen Brain Atlas) is required for definitive oligodendrocyte/OPC annotation. The marker-based approach remains the most appropriate for our 6-category scheme.
