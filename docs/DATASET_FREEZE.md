# Paper 1 analytical dataset freeze

**Status:** frozen on 2026-11-23 from the open SEA-AD MTG MERFISH h5ad.
**Issue scope:** #15 (donor inclusion/exclusion criteria and covariate table),
#16 (QC report), and #17 (versioned hashed manifest).

## Source dataset

- File: `middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad`
- Size: 497,737,587 bytes
- SHA-256: `3e3dac22446a8ce66afd07c209cd74df260b2f3d7c4085ffe4732390d2054a24`
- Access: open AWS Open Data Registry dataset; no DUA required
- Citation: Gabitto, Travaglini et al. (2024), *Nature Neuroscience* 27,
  2366–2383, doi:10.1038/s41593-024-01774-5

The raw h5ad remains local-only under `data/seaad/` and is not committed.

## Frozen inclusion and QC criteria

Criteria version `paper1-freeze-v1` is implemented by
`src/seaad_niches/dataset_freeze.py` and applied by
`scripts/freeze_paper1_dataset.py`.

1. A donor is eligible only if sex, age at death, PMI, and Continuous
   Pseudo-progression Score are complete and constant across that donor's
   cells.
2. A donor must retain at least 1,000 QC-passing analyzed cells.
3. A section must retain at least 1,000 QC-passing cells.
4. A cell passes technical QC only when it has finite x/y coordinates, a
   non-missing subclass, `Genes detected >= 1`, `Number of spots >= 1`, and
   `Supertype confidence >= 0`.
5. SEA-AD ages recorded as `90+` are right-censored in the source metadata.
   The model-ready `age` column uses the documented lower-bound encoding
   `90`, and `age_censored_90plus` records which donors carry that encoding.

These are technical freeze criteria, not biological exclusions. No downstream
niche-progression model was fit during the freeze.

## Freeze outcome

| Metric | Value |
|---|---:|
| Donors included | 27 / 27 |
| Donors excluded | 0 |
| Sections passing QC | 69 / 69 |
| Sections failing QC | 0 |
| Cells passing cell-level QC | 1,887,728 / 1,887,729 |
| QC-passing cells in passing sections | 1,887,728 |

One cell in section `H21.33.025.CX26.MTG.02.007.4.01.02` was excluded because
its recorded Supertype confidence was below zero. The section itself retained
31,446 QC-passing cells and passed the section threshold.

The existing spatial aggregates in RESULTS Tables 1–3 were generated before
this freeze and therefore include all 1,888,729 loaded cells. They remain a
valid historical pipeline result, but the donor-level progression analysis and
the final results package must use the frozen QC-passing set recorded in
`data/manifest.json` (1,887,728 cells).

The frozen donor cohort spans pseudo-progression scores 0.149616–0.889563,
includes 16 male and 11 female donors, and contains 14 donors with the
source-censored `90+` age encoding.

## Committed artifacts

| Artifact | Contents | SHA-256 |
|---|---|---|
| `reports/paper1_freeze/donor_covariates.csv` | One row per donor: model-ready `age`, `sex`, `pmi`, `pseudo_progression_score`, source age string, censoring flag, section/cell counts, inclusion status | `fb2d573aacb77592e28c9ad15e2af1fd63f29052d6985fa14bb3371354d7e6f2` |
| `reports/paper1_freeze/section_qc.csv` | One row per section: donor, cell counts, cell-level QC exclusions, technical failure counts, pass/fail status | `e3bf7b52487bf249aeabdb652a3e2e7e06c1ba8b2d135d5e94a41524156836e3` |
| `reports/paper1_freeze/qc_summary.json` | Machine-readable freeze outcome and threshold record | `b94e219d6886ae595fa5d5d90775cf14ec1d3db9ad963468310e38f5bc4e601c` |
| `reports/paper1_freeze/qc_report.md` | Generated human-readable QC report | `7151225d4edafdc88f508e7202cac3ef0d83d4e665410be9a2068deed7c0e7ae` |
| `data/manifest.json` | Versioned freeze manifest with donor list, artifact hashes, source h5ad hash, and schema/criteria versions | committed despite `data/` ignore rule; contains no raw data |

## Reproduction

```bash
python scripts/download_seaad.py --cache-dir data/seaad
python scripts/freeze_paper1_dataset.py \
  --cache-dir data/seaad \
  --output-dir reports/paper1_freeze \
  --manifest-path data/manifest.json
```

The freeze script requires the optional `anndata` dependency and reads the
h5ad in backed mode. A same-day clean-output re-run on 2026-11-23 regenerated
all four `reports/paper1_freeze/` artifacts byte-identically (SHA-256 match);
`data/manifest.json` itself differs only in its creation timestamp. The
committed `data/manifest.json` is the downstream Paper 1 freeze reference; it
contains hashes and donor/section lists only, not raw cell-level data.

## Next analysis boundary

The next modeling step is to regenerate per-donor niche abundances from the
QC-passing cells referenced by this manifest, then run the pre-registered
null families in `docs/ANALYSIS_PLAN.md`. Until those measured outputs exist,
RESULTS Tables 4–5 remain intentionally pending.
