# Paper 1 niche–progression analysis package (2026-11-28)

Measured outputs of `scripts/run_paper1_analysis.py`, the end-to-end
execution of `docs/ANALYSIS_PLAN.md` (incl. Amendment A1) on the frozen
SEA-AD analytical dataset (`data/manifest.json`, criteria version
`paper1-freeze-v1`). All values are associations only; no causal claims.

## Contents

| File | Contents |
|---|---|
| `niche_abundance_by_donor.csv` | Donor-by-niche abundance fractions (27 donors x 8 aligned niches; rows sum to 1) |
| `niche_progression_results.csv` | Table 4 source: per-niche GAM contrast, SE, CI, nominal p, donor-shuffle p + BH q, label-shuffle p + BH q (degenerate — see below) |
| `donor_shuffle_null_summary.csv` | Donor-shuffle family (B = 10,000, seed 0): observed contrast, null mean/SD, p, BH q per niche |
| `label_shuffle_null_summary.csv` | Label-shuffle family (B = 10,000, seed 1) as pre-registered — numerically invariant (degenerate) |
| `niche_alignment_map.csv` | Amendment A1 mapping: section, local niche, aligned niche, cell count (69 sections x 8). **Not committed** (bulky audit intermediate); regenerated deterministically on re-run under the pinned seeds/versions in `run_provenance.json` |
| `niche_alignment_quality.csv` | Alignment quality per niche (data description, not a result) |
| `run_provenance.json` | Seeds, permutation counts, package versions, manifest/h5ad/covariate SHA-256, timestamps |

## Headline

No niche shows a significant abundance–progression association under the
pre-registered donor-shuffle permutation null after BH FDR (smallest
q = 0.467). Niche 6's nominal OLS p (0.0445) does not survive the
pre-registered claim rule (permutation p = 0.0584, q = 0.467).

## Label-shuffle family: degenerate as pre-registered

The pre-registered label shuffle permutes niche labels across cells
within each donor/section, preserving marginal niche-label counts.
Because SEA-AD sections are nested within donors, this leaves every
donor-by-niche count — and therefore every fitted contrast — exactly
unchanged (measured null SD <= 1.3e-15; all p = 1). The family tests
nothing as specified and is **deferred** per the plan's failure-handling
rule; no post hoc null has been substituted. Repair requires an owner
plan amendment (issue #9): the natural candidate is a per-section
niche-identity permutation (re-draw each section's local-to-aligned
label mapping independently), which preserves per-section marginals
while allowing donor-level fractions to vary.

## Reproduce

```bash
python scripts/download_seaad.py --cache-dir data/seaad
python scripts/run_paper1_analysis.py --cache-dir data/seaad
```

The run verifies the h5ad size + SHA-256 against `data/manifest.json`
and fails closed on any mismatch.
