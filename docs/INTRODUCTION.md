# Introduction — Paper 1: Spatial Niche Remodeling along Alzheimer's Progression (SEA-AD)

**Series note:** This is **Paper 1 of 4** in the SEA-AD spatial Alzheimer's progression series. It is the first paper and does not build on any former paper. Papers 2 (snRNA→MERFISH state mapping), 3 (glial–neuronal signaling rewiring), and 4 (cross-region extension) all build on the dataset freezing, donor-level QC, and pseudo-progression modeling established here.

## Background

The Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD; Gabitto et al., *Nature Neuroscience* 2024) provides snRNA-seq, snATAC-seq, and MERFISH spatial transcriptomics from middle temporal gyrus of ~84 donors spanning Alzheimer's disease neuropathological change, together with quantitative neuropathology and a continuous "pseudo-progression" score. Processed data are openly available via CELLxGENE and the AWS Open Data Registry. Most published reuse of SEA-AD repeats differential-expression analyses; the matched spatial component is underexploited.

## Research questions

1. How do cell-type neighborhoods (spatial niches) in human MTG change along the AD pseudo-progression continuum?
2. Which niches are glia-dominated vs. neuron-dominated, and in what order do they remodel?
3. Can spatial graph features per donor predict neuropathological stage better than cell-type proportions alone?

## Data

| Resource | Scale | Access |
|---|---|---|
| SEA-AD snRNA-seq (MTG) | ~84 donors | Open (CELLxGENE / AWS `allen-sea-ad-atlas`) |
| SEA-AD MERFISH (MTG) | matched subset | Open (AWS) |
| Allen Brain Cell Atlas reference | 3M+ human nuclei | Open |

## Methods

Spatial neighborhood detection (SpaGCN-style graph methods), donor-level niche composition statistics, generalized additive models of niche metrics vs. pseudo-progression with covariate adjustment (age, sex, PMI), label-shuffle and donor-shuffle nulls.

## Expected contributions

- A quantitative map of niche remodeling order along AD progression.
- An open, reproducible analysis pipeline reusable by Papers 2–4.

## Scope and boundary

This repository contains planning, software, and synthetic tests. Empirical outputs, result tables, and figures are produced under the project's data-use compliance rules and released only upon owner decision. No causal or mechanistic claims are made from correlational niche–progression associations.
