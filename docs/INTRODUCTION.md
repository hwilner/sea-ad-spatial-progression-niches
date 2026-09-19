# Introduction — Paper 1: Spatial Niche Remodeling along Alzheimer's Progression (SEA-AD)

**Series note:** This is **Paper 1 of 4** in the SEA-AD spatial Alzheimer's progression series. It is the first paper and does not build on any former paper. Papers 2 (snRNA→MERFISH state mapping), 3 (glial–neuronal signaling rewiring), and 4 (cross-region extension) all build on the dataset freezing, donor-level QC, and pseudo-progression modeling established here.

## Background

Alzheimer's disease (AD) is the most common cause of dementia, and its defining brain pathology is well known: extracellular plaques of amyloid-β (Aβ) peptide and intracellular neurofibrillary tangles made of hyperphosphorylated tau protein. Pathology does not appear all at once — it accumulates over decades, following a stereotyped anatomical sequence first described by Braak and Braak, in which tau pathology begins in medial temporal regions and later engulfs the neocortex [1]. Yet pathology alone does not explain disease: what ultimately determines cognition is how the *cells* of the brain respond. The human cortex contains dozens of neuronal subtypes — excitatory glutamatergic neurons and inhibitory GABAergic interneurons of the Sst, Pvalb, and Vip families — embedded in a support tissue of glia: astrocytes (metabolic and synaptic support), microglia (resident immune cells), oligodendrocytes and their precursor cells, and vascular cells. Every one of these cell classes changes in AD.

Single-cell and single-nucleus RNA sequencing (snRNA-seq) transformed this picture. Early landmark studies of postmortem AD brain showed that transcriptional responses to pathology are strikingly cell-type specific: excitatory and inhibitory neurons, oligodendrocytes, microglia, and astrocytes each mount distinct programs [2, 3]. In mouse models, microglia were found to adopt a "disease-associated microglia" (DAM) state linked to Aβ plaques and gated by the risk gene TREM2 [4], and analogous "disease-associated astrocyte" (DAA) states were described in both aging and AD [5]. Large human cohorts — the ROSMAP studies (Religious Orders Study / Memory and Aging Project) and their follow-ups — scaled this to hundreds of donors and millions of nuclei, revealing genes and cell types associated with cognitive decline and resilience [6, 7].

These studies share one fundamental limitation: dissociation destroys geography. A nucleus tells you what it is, not where it was. Yet AD pathology is intrinsically spatial — plaques are focal objects around which microglia cluster, and vulnerable neurons die in specific layers and microenvironments. Spatial transcriptomics closes this gap. MERFISH (multiplexed error-robust fluorescence in situ hybridization) measures hundreds to thousands of RNA species in intact tissue at single-cell resolution, preserving each cell's coordinates and neighbors [8].

The Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD), led by the Allen Institute, is the largest resource combining these modalities for human AD [9]. Its flagship middle temporal gyrus (MTG) dataset includes snRNA-seq (~1.2 million nuclei) and MERFISH (~500,000 cells) from ~84 donors spanning the spectrum of AD neuropathological change, together with quantitative neuropathology (Aβ and phospho-tau burden), cognitive trajectories, and a continuous "pseudo-progression" score that orders donors along the disease continuum. SEA-AD showed that AD unfolds in at least two phases — an early phase marked by loss of specific Sst+ interneurons and inflammatory glial responses, and a later phase of accelerated pathology with loss of excitatory and Pvalb+/Vip+ neurons [9]. All processed data are openly available via CELLxGENE and the AWS Open Data Registry.

## Prior work and gap

Published analyses of SEA-AD and comparable atlases have concentrated on differential expression, cell-type abundance, and trajectory analyses of the dissociated data [6, 9]. The matched MERFISH component — the only large-scale human AD dataset in which disease stage and spatial position coexist — remains comparatively underexploited. Meanwhile, spatial-statistical methods for defining cellular neighborhoods ("niches") and quantifying cell–cell spatial co-occurrence have matured rapidly, including Bayesian spatial clustering (BayesSpace [10]) and graph-based frameworks for neighborhood enrichment and co-occurrence analysis (Squidpy [11]). What is missing is a systematic, donor-level description of how the *spatial organization* of MTG cell types — not merely their proportions — remodels along the AD continuum.

## Research questions

1. How do cell-type neighborhoods (spatial niches) in human MTG change along the AD pseudo-progression continuum?
2. Which niches are glia-dominated vs. neuron-dominated, and in what order do they remodel?
3. Can spatial graph features per donor predict neuropathological stage better than cell-type proportions alone?

## Data

| Resource | Scale | Access |
|---|---|---|
| SEA-AD snRNA-seq (MTG) | ~1.2M nuclei, ~84 donors | Open (CELLxGENE / AWS `allen-sea-ad-atlas`) [9] |
| SEA-AD MERFISH (MTG) | ~500K cells, matched subset | Open (AWS) [9] |
| SEA-AD quantitative neuropathology & cognitive metadata | per donor | Open [9] |
| Allen Brain Cell Atlas reference | 3M+ human nuclei | Open [9] |

## Methods

We first freeze a quality-controlled, donor-level analytical dataset (harmonized cell-type labels, covariates, pseudo-progression scores) that all later papers in this series inherit. On the MERFISH sections, we build k-nearest-neighbor spatial graphs, compute neighborhood composition vectors per cell, and cluster these into niches using graph-aware spatial clustering (SpaGCN-style / BayesSpace-style approaches [10, 11]). Neighbor-enrichment statistics are computed against label-shuffle null models. At the donor level, niche composition and spatial graph metrics are modeled against pseudo-progression using generalized additive models with covariate adjustment (age, sex, postmortem interval) and donor-level permutation tests to control false positives.

## Expected contributions

- A quantitative map of which spatial niches remodel, and in what order, along AD progression in MTG.
- A test of whether spatial organization carries disease-stage information beyond cell-type proportions.
- A frozen dataset and open, reproducible pipeline that serve as the shared foundation for Papers 2–4.

## Scope and boundary

This repository contains planning, software, and synthetic tests. Empirical outputs, result tables, and figures are produced under the project's data-use compliance rules and released only upon owner decision. Niche–progression associations are correlational; no causal or mechanistic claims are made in this paper.

## References

1. Braak H, Braak E. Neuropathological stageing of Alzheimer-related changes. *Acta Neuropathologica* 82, 239–259 (1991). doi:10.1007/BF00308809.
2. Mathys H, Davila-Velderrain J, Peng Z, et al. Single-cell transcriptomic analysis of Alzheimer's disease. *Nature* 570, 332–337 (2019). doi:10.1038/s41586-019-1195-2. PMID: 31042697.
3. Grubman A, Chew G, Ouyang JF, et al. A single-cell atlas of entorhinal cortex from individuals with Alzheimer's disease reveals cell-type-specific gene expression regulation. *Nature Neuroscience* 22, 2087–2097 (2019). doi:10.1038/s41593-019-0539-4.
4. Keren-Shaul H, Spinrad A, Weiner A, et al. A unique microglia type associated with restricting development of Alzheimer's disease. *Cell* 169, 1276–1290.e17 (2017). doi:10.1016/j.cell.2017.05.018.
5. Habib N, McCabe C, Medina S, et al. Disease-associated astrocytes in Alzheimer's disease and aging. *Nature Neuroscience* 23, 701–706 (2020). doi:10.1038/s41593-020-0624-8.
6. Mathys H, Boix CA, Akay LA, et al. Single-cell multiregion dissection of Alzheimer's disease. *Nature* 632, 858–868 (2024). doi:10.1038/s41586-024-07606-7. PMID: 39048816.
7. De Jager PL, Ma Y, McCabe C, et al. A multi-omic atlas of the human frontal cortex for aging and Alzheimer's disease research. *Scientific Data* 5, 180142 (2018). doi:10.1038/sdata.2018.142.
8. Chen KH, Boettiger AN, Moffitt JR, Wang S, Zhuang X. Spatially resolved, highly multiplexed RNA profiling in single cells. *Science* 348, aaa6090 (2015). doi:10.1126/science.aaa6090.
9. Gabitto MI, Travaglini KJ, Rachleff VM, et al. Integrated multimodal cell atlas of Alzheimer's disease. *Nature Neuroscience* 27, 2366–2383 (2024). doi:10.1038/s41593-024-01774-5.
10. Zhao E, Stone MR, Ren X, et al. Spatial transcriptomics at subspot resolution with BayesSpace. *Nature Biotechnology* 39, 1375–1384 (2021). doi:10.1038/s41587-021-00935-2.
11. Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nature Methods* 19, 171–178 (2022). doi:10.1038/s41592-021-01358-2.
