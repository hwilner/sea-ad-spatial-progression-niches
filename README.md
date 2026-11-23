# SEA-AD Spatial Progression Niches

This independent research repository plans and tracks a spatial analysis of cell-type neighborhood remodeling along the Alzheimer's disease pseudo-progression continuum in the SEA-AD atlas. It provides analysis utilities plus a one-command pipeline over the **open** SEA-AD MTG MERFISH dataset (no DUA required).

## Introduction

Alzheimer's disease remodels not just cell-type abundances but the *spatial neighborhoods* cells live in. This project asks whether microanatomical niches in the middle temporal gyrus change systematically along the SEA-AD pseudo-progression continuum, using the open MERFISH atlas alone. See [docs/INTRODUCTION.md](docs/INTRODUCTION.md) for the full introduction.

## Extended introduction

A zero-background guide to the dataset, the spatial-niche concept, and the analysis strategy — written for readers new to spatial transcriptomics — is in [docs/EXTENDED_INTRODUCTION.md](docs/EXTENDED_INTRODUCTION.md).

## Materials and methods

Open SEA-AD MTG MERFISH cells (1,888,729 cells, 27 donors, 69 sections; scripted, checksummed download) → technical dataset freeze (1,887,728 QC-passing cells, all 27 donors and 69 sections retained) → per-section k-NN neighbor graphs → subclass neighborhood-enrichment z-scores against a label-permutation null → niche detection by clustering neighbor-composition vectors (Leiden if installed, deterministic k-means fallback) → donor-level niche-abundance vs. pseudo-progression association (pre-registered covariate-adjusted GAM + permutation nulls; associations only, no causal claims). Full scope, design decisions, and hygiene rules: [docs/METHODS.md](docs/METHODS.md); data access, provenance, and freeze criteria: [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md), [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md), [docs/DATASET_FREEZE.md](docs/DATASET_FREEZE.md).

## Results

Completed and pre-registered results are collected, in paper style with captions, in [docs/RESULTS.md](docs/RESULTS.md). Headline real-data findings from the committed aggregates in `reports/` (run over the open MERFISH dataset):

| Result | Value | Source |
|---|---|---|
| Analytical dataset | 1,888,729 cells, 27 donors, 69 sections, 24 subclasses | `reports/seaad_merfish_run_summary.json` (RESULTS.md Table 1) |
| Top enriched neighbor pair | L4 IT–L5 IT, mean z = 31.20 ± 1.74 | `reports/seaad_merfish_top_pairs.csv` (RESULTS.md Table 2) |
| Top depleted neighbor pair | L2/3 IT–Oligodendrocyte, mean z = −55.32 ± 3.44 | `reports/seaad_merfish_top_pairs.csv` (RESULTS.md Table 2) |
| Niche detection | 8 stable niches per section in all 69 sections; largest-niche fraction 0.202 ± 0.033 | `reports/seaad_merfish_niche_summary.csv` (RESULTS.md Table 3) |
| Validation | Suite: 32 passed, 0 failed, 0 skipped (with the optional h5ad reader installed) | RESULTS.md Table 6 |

Niche–progression association modeling is pre-registered (`docs/ANALYSIS_PLAN.md`), implemented, and synthetic-validated; the donor-level score/covariate freeze is committed (`data/manifest.json`, `reports/paper1_freeze/`). Real-data fitting now awaits QC-aware per-donor niche abundance aggregation and the pre-registered null families (RESULTS.md Tables 4–5, Figures 1–2).

### Planned work (issue-mapped)

| Planned work | Expected outcome | Issues |
|---|---|---|
| SEA-AD MERFISH ingestion, QC, and donor freeze | A fixed, versioned analytical dataset spanning the pseudo-progression score. | #5 (done), #15 (done), #16 (done), #17 (done), #2 |
| Spatial niche detection per donor | Quantitative niche maps per donor and disease stage. | #3 (done), #14 |
| Niche–progression association modeling | Identification of neighborhoods that remodel along pseudo-progression. | #7 (plan), #10 (model), #6, #8, #12, #22, #23 |
| Sensitivity and null analyses | Negative-control niches and label-shuffle nulls; no causal claims. | #9, #1, #19 (harness), #20, #21, #18 |
| Results package and manuscript | Frozen results package, reproducibility re-run, Paper 1 outline. | #4, #11, #13 |

**Current status:** the open SEA-AD MTG MERFISH dataset (Gabitto et al. 2024; AWS Open Data Registry) is scripted, checksummed, and cached by `scripts/download_seaad.py`; the Paper 1 donor/covariate freeze is committed at `data/manifest.json` with QC artifacts in `reports/paper1_freeze/`; neighbor enrichment and niche detection run end-to-end on real cells via `scripts/analyze_seaad.py` (aggregates in `reports/`). The pre-registered donor-level GAM and one-choice sensitivity harness are implemented and synthetic-smoke-tested; real-data progression modeling awaits QC-aware donor-level niche abundance aggregation and the pre-registered null families.

## What is included

| Path | Contents |
|---|---|
| `src/seaad_niches/spatial.py` | k-NN / radius neighbor-graph construction and neighborhood enrichment z-scores (label-permutation null). |
| `src/seaad_niches/niches.py` | Niche assignment by clustering neighbor-composition vectors (Leiden if installed, deterministic k-means fallback). |
| `src/seaad_niches/progression.py` | Donor-level niche abundance vs. pseudo-progression association (Spearman + score-permutation test). |
| `src/seaad_niches/niche_progression.py` | Pre-registered covariate-adjusted GAM (cubic B-spline + age/sex/PMI), returning effect sizes with uncertainty. |
| `src/seaad_niches/dataset_freeze.py` | Paper 1 donor inclusion, cell/section QC, freeze summary, and hashed manifest machinery. |
| `src/seaad_niches/simulate.py` | Synthetic tissue generator with planted niches and a planted progression gradient. |
| `src/seaad_niches/io.py` | SEA-AD loading interface; raises a clear "requires SEA-AD access" error until `scripts/download_seaad.py` has staged the open dataset (writes `MANIFEST.json`). |
| `src/seaad_niches/download.py` | Public SEA-AD MTG MERFISH download/cache from the AWS Open Data Registry bucket with size + SHA-256 integrity checks. |
| `scripts/download_seaad.py` | One-command download + verification of the open MERFISH h5ad (~500 MB). |
| `scripts/analyze_seaad.py` | Real-data pipeline: per-section neighbor enrichment + niche detection; writes small aggregates to `reports/`. |
| `scripts/freeze_paper1_dataset.py` | Paper 1 donor covariate freeze, QC report, and hashed `data/manifest.json` from the staged h5ad. |
| `Makefile` | `make real-pipeline` = download + analyze + dataset freeze; `make test` runs the suite. |
| `reports/` | Small derived CSV/JSON aggregates from the real-data run (enrichment z-scores, niche summaries). |
| `analysis/sensitivity.py` | One-choice-at-a-time sensitivity harness with an explicit varied-choice log. |
| `tests/` | Synthetic/fixture tests: planted neighbor enrichment, planted niche recovery, planted gradient detection, GAM recovery, sensitivity-harness smoke tests, dataset-freeze criteria/QC/manifest tests, io gating, download/cache machinery (mocked network). |
| `docs/` | Research status, methods scope, pre-registered analysis plan, results, data-access status, data sources/versions, and contribution guidance. |
| `pyproject.toml` / `requirements.txt` | Src-layout package `seaad_niches` (Python >= 3.10); core deps numpy, pandas, scipy, scikit-learn. |
| `CONTRIBUTING.md` | Setup, testing, backlog workflow, code style (ruff), and PR process. |
| `LICENSE` | MIT (copyright Harel Wilner 2026). |

## Use and validation

Install the Python requirements in an isolated environment, then run:

```bash
pip install -e ".[dev]"
python -m pytest -q
```

Real-data pipeline (open data, no DUA):

```bash
make real-pipeline
```

## Keywords

Alzheimer's disease, SEA-AD, MERFISH, spatial transcriptomics, cell niches, computational neuroscience, reproducible research.

## Documentation

- [Introduction for new readers](docs/INTRODUCTION.md)
- [Extended introduction — zero-background guide](docs/EXTENDED_INTRODUCTION.md)
- [Extended methods & materials](docs/METHODS.md)
- [Pre-registered niche–progression analysis plan](docs/ANALYSIS_PLAN.md)
- [Results — captioned paper-style tables](docs/RESULTS.md)
- [Paper 1 analytical dataset freeze](docs/DATASET_FREEZE.md)
- [SEA-AD data access status](docs/DATA_ACCESS.md)
- [Data sources and versions](docs/DATA_SOURCES.md)
