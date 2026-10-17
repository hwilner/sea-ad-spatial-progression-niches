# Results

Paper-style results for the SEA-AD spatial progression niches project.
Tables 1–3 report **completed** real-data results, computed from the
committed aggregates in `reports/` (produced by
`scripts/analyze_seaad.py` on the open SEA-AD MTG MERFISH dataset; see
`docs/DATA_SOURCES.md`). Tables 4–5 and Figures 1–2 are
**pre-registered placeholders**: their structure is fixed in advance of
the analyses that will fill them, and each is linked to the issue that
authorizes the corresponding analysis. Table 6 records computational
validation of the codebase.

## Completed results

**Table 1. Real-data run summary: SEA-AD MTG MERFISH spatial pipeline.**
One end-to-end run of `scripts/analyze_seaad.py` over the open
all-donors MERFISH h5ad (Gabitto et al. 2024; AWS Open Data Registry;
no DUA). Per section: k-NN neighbor graph on (x, y) coordinates,
subclass neighborhood-enrichment z-scores against a label-permutation
null, and niche detection by k-means clustering of neighbor-composition
vectors. Values verbatim from `reports/seaad_merfish_run_summary.json`.

| Quantity | Value |
|---|---|
| Cells loaded | 1,888,729 |
| Donors | 27 |
| Sections analyzed | 69 |
| Cell-type labels (subclasses) | 24 |
| Neighbors per cell (k) | 6 |
| Label permutations per section | 100 |
| Niches per section (k-means) | 8 |
| Random seed | 0 |

**Table 2. Most enriched and most depleted subclass neighbor pairs.**
Mean ± SEM of the per-section enrichment z-score across the 69 sections
(label-permutation null, 100 permutations per section, k = 6 neighbor
graph). Top five enriched and five most depleted of the 276 unique
subclass pairs; the full ranked table is
`reports/seaad_merfish_top_pairs.csv`. Enrichment concentrates within
cortical excitatory lineages (intratelencephalic projections across
layers), within inhibitory lineages, and between endothelial cells and
VLMCs (vascular niche); the strongest depletion is between
oligodendrocytes and neuronal subclasses, consistent with gray/white
compartment segregation.

| Direction | Subclass A | Subclass B | Mean z | SEM z |
|---|---|---|---|---|
| Enriched | L4 IT | L5 IT | 31.20 | 1.74 |
| Enriched | L2/3 IT | Vip | 30.23 | 1.35 |
| Enriched | Endothelial | VLMC | 23.07 | 1.53 |
| Enriched | L5 IT | L5/6 NP | 18.38 | 0.89 |
| Enriched | L6 CT | L6 IT | 16.89 | 0.77 |
| Depleted | Oligodendrocyte | Vip | −19.72 | 1.12 |
| Depleted | L2/3 IT | L5 IT | −27.19 | 1.80 |
| Depleted | L5 IT | Oligodendrocyte | −28.63 | 2.10 |
| Depleted | L4 IT | Oligodendrocyte | −35.62 | 2.38 |
| Depleted | L2/3 IT | Oligodendrocyte | −55.32 | 3.44 |

**Table 3. Per-section niche-size summary (k-means, 8 niches per
section).** Niche assignments are stable in size: every one of the 69
sections yielded exactly 8 niches, and the largest/smallest niche
fractions vary modestly across sections and donors, with no degenerate
(single-cell or section-swallowing) clusters. Cells per section range
from 5,732 to 128,300 (27 donors). Summary statistics computed over the
69 rows of `reports/seaad_merfish_niche_summary.csv`; totals per the
bit-verified run summary (Table 1).

| Statistic | Largest niche fraction | Smallest niche fraction |
|---|---|---|
| Mean across sections | 0.202 | 0.069 |
| SD across sections | 0.033 | 0.015 |
| Minimum | 0.151 | 0.035 |
| Maximum | 0.330 | 0.098 |

## Pre-registered result structures (pending analyses)

**Table 4 (pre-registered; pending). Niche abundance vs.
pseudo-progression associations.** One row per niche: Spearman r between
donor-level niche abundance and the donor pseudo-progression score,
two-sided score-permutation p-value, Benjamini–Hochberg q-value across
niches, and n_donors. Analysis authorized by the pre-registration in
issue #7 (GAM formula and null design), modeled under issue #6, table
produced under issue #8. Blocked by the donor-level score/covariate
table (issues #2, #15) — no pseudo-progression scores are committed
yet, so this table is intentionally empty.

| Niche | Spearman r | Permutation p | BH q | n_donors |
|---|---|---|---|---|
| — (pending issues #2/#15 → #7 → #6 → #8) | | | | |

**Table 5 (pre-registered; pending). Null-model and sensitivity family.**
Label-shuffle and donor-shuffle null families (issue #9), the
sensitivity harness over neighbor-graph k/radius, covariate,
leave-one-out, and region-within-MTG variants (issues #19–#21), and the
effect-direction stability summary (issue #18). One row per analysis:
statistic, null mean ± SD, p-value. No causal claims are made from any
row of this table.

| Analysis | Statistic | Null mean ± SD | p |
|---|---|---|---|
| — (pending issues #9, #18–#21) | | | |

**Figure 1 (planned). Per-donor niche maps along the pseudo-progression
continuum.** Script under issue #14.

**Figure 2 (planned). Niche–progression summary figure and progression
curves.** Scripts under issues #12 and #22; docs summary under issue #23.

## Validation record

**Table 6. Computational validation of the analysis codebase
(completed).** Synthetic test suite (planted neighbor enrichment,
planted niche recovery, planted progression-gradient detection, io
gating, download/cache machinery with mocked network). Run this run on a
clean-room copy of the committed files (git blob-SHA parity verified for
all 13 source/test files against the `main` branch). The result
reproduces the independent verification recorded in `results.md` on
branch `mavis/verify-and-reproduce` (2026-09-29, Python 3.11). The
`ConstantInputWarning` emitted during the constant-score null test is
expected by design (undefined correlation on constant input, handled as
NaN-safe zero).

| Branch | Tests passed | Failed | Skipped | Wall clock | Environment |
|---|---|---|---|---|---|
| main | 18 | 0 | 1 | ~3 s | Python 3.12, CPU (numpy/pandas/scipy/scikit-learn) |
| mavis/verify-and-reproduce | 18 | 0 | 1 | 5.52 s | Python 3.11, CPU (recorded 2026-09-29 in branch `results.md`; branch = main + that file, so code is identical) |

The single skipped test (`test_load_merfish_cells_from_staged_cache`)
requires the optional `anndata` dependency and is skipped when it is not
installed; it does not count as a pass. Heavy-compute steps (full
real-data re-run, issue #11) are **deferred (needs environment/data)**:
the committed `reports/` aggregates above stand as the released result.
