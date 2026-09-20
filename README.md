# SEA-AD Spatial Progression Niches (Paper 1)

This independent research repository plans and tracks a spatial analysis of cell-type neighborhood remodeling along the Alzheimer's disease pseudo-progression continuum in the SEA-AD atlas. It provides analysis utilities plus a one-command pipeline over the **open** SEA-AD MTG MERFISH dataset (no DUA required).

## Series position

This is **Paper 1** of the SEA-AD spatial Alzheimer's progression series (4 papers). It is the foundation of the series; Papers 2–4 build on its outputs.

## Research plan

| Planned work | Expected outcome |
|---|---|
| SEA-AD MERFISH + snRNA-seq ingestion and QC | A fixed, versioned analytical dataset spanning the pseudo-progression score. |
| Spatial niche detection per donor | Quantitative niche maps per donor and disease stage. |
| Niche–progression association modeling | Identification of neighborhoods that remodel along pseudo-progression. |
| Sensitivity and null analyses | Negative-control niches and label-shuffle nulls; no causal claims. |

**Current status:** the open SEA-AD MTG MERFISH dataset (Gabitto et al. 2024; AWS Open Data Registry) is scripted, checksummed, and cached by `scripts/download_seaad.py`; neighbor enrichment and niche detection run end-to-end on real cells via `scripts/analyze_seaad.py` (1,888,729 cells, 27 donors, 69 sections; aggregates in `reports/`). Pseudo-progression association modeling awaits donor-level scores.

## What is included

| Path | Contents |
|---|---|
| `src/seaad_niches/spatial.py` | k-NN / radius neighbor-graph construction and neighborhood enrichment z-scores (label-permutation null). |
| `src/seaad_niches/niches.py` | Niche assignment by clustering neighbor-composition vectors (Leiden if installed, deterministic k-means fallback). |
| `src/seaad_niches/progression.py` | Donor-level niche abundance vs. pseudo-progression association (Spearman + score-permutation test). |
| `src/seaad_niches/simulate.py` | Synthetic tissue generator with planted niches and a planted progression gradient. |
| `src/seaad_niches/io.py` | SEA-AD loading interface; raises a clear "requires SEA-AD access" error until `scripts/download_seaad.py` has staged the open dataset (writes `MANIFEST.json`). |
| `src/seaad_niches/download.py` | Public SEA-AD MTG MERFISH download/cache from the AWS Open Data Registry bucket with size + SHA-256 integrity checks. |
| `scripts/download_seaad.py` | One-command download + verification of the open MERFISH h5ad (~500 MB). |
| `scripts/analyze_seaad.py` | Real-data pipeline: per-section neighbor enrichment + niche detection; writes small aggregates to `reports/`. |
| `Makefile` | `make real-pipeline` = download + analyze; `make test` runs the suite. |
| `reports/` | Small derived CSV/JSON aggregates from the real-data run (enrichment z-scores, niche summaries). |
| `tests/` | Synthetic tests: planted neighbor enrichment, planted niche recovery, planted gradient detection, io gating, download/cache machinery (mocked network). |
| `docs/` | Research status, methods scope, data-access status, data sources/versions, and contribution guidance. |
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
- [SEA-AD data access status](docs/DATA_ACCESS.md)
- [Data sources and versions](docs/DATA_SOURCES.md)
