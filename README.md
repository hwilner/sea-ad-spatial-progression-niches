# SEA-AD Spatial Progression Niches (Paper 1)

This independent research repository plans and tracks a spatial analysis of cell-type neighborhood remodeling along the Alzheimer's disease pseudo-progression continuum in the SEA-AD atlas. It provides data-free analysis utilities for transparent review and extension.

## Series position

This is **Paper 1** of the SEA-AD spatial Alzheimer's progression series (4 papers). It is the foundation of the series; Papers 2–4 build on its outputs.

## Research plan

| Planned work | Expected outcome |
|---|---|
| SEA-AD MERFISH + snRNA-seq ingestion and QC | A fixed, versioned analytical dataset spanning the pseudo-progression score. |
| Spatial niche detection per donor | Quantitative niche maps per donor and disease stage. |
| Niche–progression association modeling | Identification of neighborhoods that remodel along pseudo-progression. |
| Sensitivity and null analyses | Negative-control niches and label-shuffle nulls; no causal claims. |

**Current status:** planning stage; no empirical calculation has been run.

## What is included

| Path | Contents |
|---|---|
| `src/` | In-memory spatial-graph, niche-detection, and statistics utilities (added as issues are completed). |
| `tests/` | Synthetic tests for spatial-graph and statistics invariants. |
| `docs/` | Research status, methods scope, deferred directions, and contribution guidance. |

## Use and validation

Install the Python requirements in an isolated environment, then run:

```bash
python -m pytest -q
```

## Keywords

Alzheimer's disease, SEA-AD, MERFISH, spatial transcriptomics, cell niches, computational neuroscience, reproducible research.

## Documentation

- [Introduction for new readers](docs/INTRODUCTION.md)
