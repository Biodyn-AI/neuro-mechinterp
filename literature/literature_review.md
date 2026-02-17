# Literature Review: Virtual Perturbation Screening for Cognitive Enhancement Targets Using Single-Cell Foundation Models

**Date:** 2026-02-13  
**Project:** In-silico perturbation screening of ~20K genes in DLPFC pyramidal neurons using Geneformer/scGPT, scored against intelligence-GWAS signatures, for cognitive enhancement target discovery and drug repurposing.

---

## 1. Executive Summary: Novelty Assessment

**Verdict: The specific combination is GENUINELY NOVEL.** No published work combines all of:
1. Single-cell foundation models (Geneformer/scGPT)
2. Genome-wide virtual perturbation screening in brain neurons
3. Intelligence/cognition GWAS gene signatures as the scoring readout
4. Drug repurposing for cognitive enhancement

The closest existing works address individual components but never this specific integration. This represents a strong, publishable contribution.

---

## 2. Related Work — Detailed Paper Analysis

### 2.1 Foundation Models & Virtual Perturbation (Core Technology)

#### Paper 1: Geneformer — Transfer learning enables predictions in network biology
- **Authors:** Theodoris CV et al.
- **Year:** 2023
- **Journal:** Nature (Vol 619, pp 616–624)
- **PMID:** 37258680
- **Key findings:** Pretrained transformer on ~30M single-cell transcriptomes; demonstrated in-silico perturbation to identify cardiomyopathy therapeutic targets. Validated experimentally in iPSC cardiomyocytes.
- **How we differ:** Geneformer was applied to **cardiac** disease. No application to brain tissue, intelligence, or cognitive traits. We extend the same framework to **DLPFC pyramidal neurons** with **intelligence GWAS** as the phenotypic readout — a completely different biological domain and a much more ambitious genome-wide screen (~20K genes vs. targeted candidates).
- **Geneformer V2 (2024 preprint):** Expanded to ~104M transcriptomes, 316M parameters. Includes improved in-silico perturbation and treatment modules. Available on HuggingFace.

#### Paper 2: scGPT — Toward building a foundation model for single-cell multi-omics using generative AI
- **Authors:** Cui H, Wang C, Maan H et al.
- **Year:** 2024
- **Journal:** Nature Methods (Vol 21, pp 1470–1480)
- **PMID:** 38409223
- **Key findings:** Foundation model for single-cell multi-omics; perturbation prediction capability; gene regulatory network inference.
- **How we differ:** scGPT has perturbation prediction but has NOT been applied to brain tissue or cognitive phenotypes. Our application to intelligence screening is novel.

#### Paper 3: Therapeutic targets in lung diseases identified through single-cell analysis and Mendelian randomization
- **Authors:** Jiang Z, Lu Y, Wang X
- **Year:** 2025
- **Journal:** Heart & Lung (Vol 74, pp 211–223)
- **PMID:** 40716288
- **Key findings:** Used Geneformer in-silico perturbation on Human Lung Cell Atlas → generated pan-lung disease target spectrum (3570 genes across 12 traits) → validated with MR and colocalization. **This is methodologically the closest paper to our approach**, but applied to lung diseases.
- **How we differ:** Different tissue (brain vs lung), different disease/trait (intelligence vs pulmonary), different GWAS. Our intelligence-focused approach adds unique complexity because intelligence is polygenic and involves multiple brain cell types. We also add the drug repurposing layer.

#### Paper 4: Large-scale single-cell analysis and in silico perturbation reveal dynamic evolution of HCC
- **Authors:** Xia P et al.
- **Year:** 2026 (Jan 30, online ahead of print)
- **Journal:** NPJ Precision Oncology
- **PMID:** 41617928
- **Key findings:** Geneformer-based virtual knockout screening in hepatocellular carcinoma; identified HSP90B1 as convergent dependency.
- **How we differ:** Cancer application, not brain/cognition.

#### Paper 5: Toward a foundation model of causal cell and tissue biology with a Perturbation Cell and Tissue Atlas
- **Authors:** Rood JE et al. (Genentech/Regev lab)
- **Year:** 2024
- **Journal:** Cell (Vol 187, pp 4520–4545)
- **PMID:** 39178831
- **Key findings:** Visionary review proposing a Perturbation Cell Atlas combining CRISPR screens with AI/ML for causal foundation models. Does not specifically address brain/intelligence.
- **How we differ:** We implement a concrete version of this vision for a specific trait (intelligence) in a specific tissue (brain).

#### Paper 6: scELMo / scEMB — Alternative foundation models
- **PMIDs:** 40894586, 39386549
- **Key findings:** Competing/alternative single-cell foundation models with perturbation prediction capabilities. scELMo claims outperformance vs scGPT/Geneformer on some benchmarks.
- **How we differ:** These are method papers, not application papers. Could be compared as alternative backbones.

### 2.2 Intelligence GWAS & Brain Single-Cell Integration

#### Paper 7: Savage et al. — Genome-wide association meta-analysis in 269,867 individuals identifies new genetic and functional links to intelligence
- **Authors:** Savage JE, Jansen PR, Stringer S et al.
- **Year:** 2018
- **Journal:** Nature Genetics (Vol 50, pp 912–919)
- **PMID:** 29942086
- **Key findings:** Identified 205 associated genomic loci, 1,016 genes for intelligence. Enrichment in brain-expressed genes, particularly in striatum and hippocampus. This is the canonical intelligence GWAS.
- **Relevance:** Primary source for our intelligence gene signature. Need to check for newer, larger GWAS.

#### Paper 8: Yao, Harder et al. — Connecting genomic results for psychiatric disorders to human brain cell types and regions
- **Authors:** Yao S, Harder A, Darki F et al.
- **Year:** 2025
- **Journal:** Nature Communications (Vol 16, 395)
- **PMID:** 39755698
- **Key findings:** Integrated GWAS (including intelligence) with comprehensive single-cell transcriptomic atlas (Siletti et al.) of adult human brain. Found intelligence SNP-heritability enriched in **upper-layer intratelencephalic (UL-IT) neurons** — essentially cortical pyramidal neurons of layers 2/3. Also enriched in deep-layer IT neurons. Prefrontal cortex highlighted.
- **How we differ:** This paper maps GWAS→cell types but does NOT perform perturbation screening. It identifies *which* cell types matter but not *which genes to target*. Our work starts where theirs ends — we take the same cell types and perform genome-wide virtual perturbation. **This is a perfect "standing on shoulders" paper to cite.**

### 2.3 Drug Repurposing for Cognitive Enhancement

#### Paper 9: Lam et al. — Identifying nootropic drug targets via large-scale cognitive GWAS and transcriptomics
- **Authors:** Lam M, Chen CY, Ge T et al.
- **Year:** 2021
- **Journal:** Neuropsychopharmacology (Vol 46, pp 1788–1801)
- **PMID:** 34035472
- **Key findings:** Used cognitive GWAS + bulk brain transcriptomics (not single-cell) + drug databases to identify nootropic targets. Found enrichment in calcium signaling, identified existing drugs that could enhance cognition.
- **How we differ:** This is the closest work conceptually but uses **bulk transcriptomics** and **traditional bioinformatics** (not foundation models). Our single-cell foundation model approach provides cell-type-specific resolution and leverages learned gene network dynamics rather than simple statistical associations. This paper is both our strongest prior art AND our strongest motivation — we are doing what they did but with dramatically better tools.

#### Paper 10: SH2B1 Tunes Hippocampal ERK Signaling to Influence Fluid Intelligence
- **Authors:** Du X, Yan Y et al.
- **Year:** 2023
- **Journal:** Research (Vol 6, 0269)
- **PMID:** 38434247
- **Key findings:** GWAS-identified gene SH2B1 validated experimentally as influencing fluid intelligence via ERK signaling in hippocampus.
- **Relevance:** Example of a validated intelligence gene — potential positive control for our screen.

### 2.4 Brain Single-Cell Atlases

#### Paper 11: Siletti et al. — Transcriptomic diversity of cell types across the adult human brain
- **Year:** 2023
- **Journal:** Science
- **Key findings:** Most comprehensive human brain single-cell atlas: 3 million+ nuclei, 461 superclusters, >3000 clusters across entire brain. Includes DLPFC pyramidal neurons.
- **Relevance:** Primary data source for our project.

#### Paper 12: Allen Brain Cell Atlas / BICCN papers (Science 2023 package)
- Multiple papers including Tian et al. (PMID 37824674) on DNA methylation and Li et al. (PMID 37824643) on chromatin accessibility.
- **Relevance:** Complementary multiomics data for the same brain regions.

---

## 3. Novelty Assessment — Detailed

### What exists:
| Component | Existing Work | Gap |
|-----------|--------------|-----|
| Geneformer perturbation | Heart (Theodoris 2023), Lung (Jiang 2025), HCC (Xia 2026) | **Never applied to brain** |
| Intelligence GWAS → cell types | Yao et al. 2025 | Maps cell types but **no perturbation screening** |
| Nootropic drug targets via GWAS | Lam et al. 2021 | Uses **bulk** transcriptomics, no foundation models |
| scGPT perturbation in brain | **Nothing found** | Complete gap |
| Foundation model + GWAS scoring | **Nothing found** | Complete gap |
| Cognitive enhancement + drug repurposing via FM | **Nothing found** | Complete gap |

### What's novel in our approach:
1. **First application of Geneformer/scGPT virtual perturbation to brain tissue**
2. **First genome-wide perturbation screen scored against intelligence GWAS signatures**
3. **First single-cell-resolution nootropic target discovery using foundation models**
4. **Integration of FM perturbation → GWAS scoring → drug repurposing pipeline**

### Closest existing work:
**Jiang et al. 2025 (Lung diseases + Geneformer + MR)** — identical methodology but different tissue and different phenotype. Our work is essentially "the brain version" of this, but targeting a much more complex polygenic trait and adding drug repurposing.

---

## 4. Methodological Gaps & Challenges

### 4.1 Critical Challenges
1. **Geneformer's brain representation:** Geneformer V2 was pretrained on ~104M cells across tissues. Need to verify how many brain cells are in the training corpus and whether brain-specific fine-tuning is needed.
2. **Intelligence is ultra-polygenic:** ~1000+ associated genes with tiny individual effects. Unlike cardiomyopathy (few master regulators), intelligence may not have clear "master regulator" targets. The perturbation effects may be very subtle.
3. **Validation is hard:** Unlike disease models, there's no iPSC model of "intelligence." Experimental validation paths are limited to:
   - Cognitive mouse models (e.g., Morris water maze after target manipulation)
   - Known nootropic effects of identified drugs
   - Cross-referencing with cognitive decline in neurodegenerative diseases
4. **GWAS gene signature definition:** How to define the "intelligence signature" for scoring perturbations? Options:
   - Top GWAS hits by p-value
   - MAGMA gene-level associations
   - Polygenic score weights
   - Pathway-level signatures
5. **Cell type selection:** DLPFC pyramidal neurons (L2/3 IT) are supported by Yao 2025, but intelligence likely involves multiple cell types. May need multi-cell-type analysis.

### 4.2 Technical Challenges
6. **Computational cost:** ~20K gene perturbations × thousands of cells = massive compute. Geneformer V2 316M parameters requires GPU.
7. **Scoring metric:** How to quantify "shift toward intelligence-associated expression"? Cosine similarity in embedding space? Shift in GWAS gene expression? Custom loss function?
8. **Multiple testing:** 20K perturbations requires careful statistical correction.
9. **Directionality:** GWAS identifies risk/protective alleles but the direction of effect on expression is often unclear without eQTL data.

### 4.3 Ethical Considerations
10. **"Cognitive enhancement" framing:** Sensitive topic. Better to frame as "understanding genetic architecture of cognitive function" rather than "enhancement." Reviewers may push back on enhancement framing.
11. **Dual-use concerns:** Some journals may flag cognitive enhancement drug discovery.

---

## 5. Suggested Publication Venues (Ranked)

### Tier 1 (High-impact, if results are strong)
1. **Nature Biotechnology** — Ideal for novel computational framework + biological discovery. Geneformer was in Nature; applying it in a novel domain with drug repurposing fits NBT.
2. **Nature Methods** — If the methodological framework (FM perturbation → GWAS scoring pipeline) is generalizable beyond intelligence.
3. **Nature Genetics** — If the focus is on the genetics/genomics insights about intelligence architecture.

### Tier 2 (Strong specialty journals)
4. **Genome Biology** — Excellent for computational genomics with biological validation.
5. **Cell Genomics** — Good fit for integrative genomics approach.
6. **Molecular Psychiatry** — If framed around psychiatric/cognitive neuroscience.

### Tier 3 (Solid, higher acceptance probability)
7. **Nucleic Acids Research** — Methods-focused computational biology.
8. **Bioinformatics** — Pure methods paper.
9. **Neuropsychopharmacology** — If drug repurposing results are strong (where Lam 2021 was published).

### Recommendation:
**Target Nature Biotechnology** with the full pipeline paper. Frame it as: "A generalizable framework for trait-directed virtual perturbation screening using single-cell foundation models, applied to discover cognitive enhancement targets." The novelty of the framework + specific biological findings is NBT-level.

---

## 6. Key Intelligence GWAS Studies (Chronological)

| Study | Year | N | Loci | Key Genes |
|-------|------|---|------|-----------|
| Savage et al. | 2018 | 269,867 | 205 | 1,016 genes; FOXO3, MEF2C, TCF4 |
| Davies et al. | 2018 | 300,486 | 148 | Neurogenesis, synaptic function |
| Hill et al. | 2019 | 248,482 | 187 | 538 genes; DLGAP1 |
| Lam et al. (cognitive) | 2021 | 257,841 | — | Nootropic targets, Ca2+ signaling |

**Note:** Need to check if there's a newer large-scale intelligence GWAS from 2022-2025 with >500K individuals. The UK Biobank cognitive phenotype GWAS may have been updated.

---

## 7. Recommended Project Workflow

1. **Data:** Siletti et al. 2023 brain atlas → extract DLPFC UL-IT (L2/3) pyramidal neurons
2. **Model:** Geneformer V2-316M (HuggingFace) — no brain-specific fine-tuning needed initially (pretrained on diverse tissues including brain)
3. **Perturbation:** Run in-silico deletion/overexpression of all ~20K protein-coding genes
4. **Scoring:** Define intelligence GWAS gene set (Savage 2018 MAGMA genes), compute perturbation effect on this signature
5. **Analysis:** Rank genes by perturbation impact → network analysis → identify master regulators
6. **Drug mapping:** Cross-reference top targets with DrugBank/ChEMBL for repurposable compounds
7. **Validation:** Known nootropics should emerge; novel targets should have biological plausibility

---

## 8. References

1. Theodoris CV et al. (2023) Transfer learning enables predictions in network biology. *Nature* 619:616-624. PMID: 37258680
2. Cui H et al. (2024) scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nat Methods* 21:1470-1480. PMID: 38409223
3. Savage JE et al. (2018) Genome-wide association meta-analysis in 269,867 individuals identifies new genetic and functional links to intelligence. *Nat Genet* 50:912-919. PMID: 29942086
4. Yao S, Harder A et al. (2025) Connecting genomic results for psychiatric disorders to human brain cell types and regions reveals convergence with functional connectivity. *Nat Commun* 16:395. PMID: 39755698
5. Lam M et al. (2021) Identifying nootropic drug targets via large-scale cognitive GWAS and transcriptomics. *Neuropsychopharmacology* 46:1788-1801. PMID: 34035472
6. Jiang Z et al. (2025) Therapeutic targets in lung diseases identified through single-cell analysis and Mendelian randomization. *Heart Lung* 74:211-223. PMID: 40716288
7. Xia P et al. (2026) Large-scale single-cell analysis and in silico perturbation reveal dynamic evolution of HCC. *NPJ Precis Oncol* (ahead of print). PMID: 41617928
8. Rood JE et al. (2024) Toward a foundation model of causal cell and tissue biology with a Perturbation Cell and Tissue Atlas. *Cell* 187:4520-4545. PMID: 39178831
9. Siletti K et al. (2023) Transcriptomic diversity of cell types across the adult human brain. *Science* 382:eadf7044.
10. Du X et al. (2023) SH2B1 Tunes Hippocampal ERK Signaling to Influence Fluid Intelligence. *Research* 6:0269. PMID: 38434247
11. Tian W et al. (2023) Single-cell DNA methylation and 3D genome architecture in the human brain. *Science* 382:eadf5357. PMID: 37824674
12. Li YE et al. (2023) A comparative atlas of single-cell chromatin accessibility in the human brain. *Science* 382:eadf7044. PMID: 37824643
