# Pre-registered analysis plan: niche–progression GAM

**Registered:** 2026-11-18  
**Issue:** #7  
**Scope:** planning and synthetic validation only. This plan freezes the model
and null design **before** any niche–progression model is fit to the frozen
SEA-AD donor dataset.

## Scientific question

Do donor-level microanatomical niche abundances in the SEA-AD middle temporal
gyrus MERFISH atlas change monotonically along the donor pseudo-progression
continuum, after adjustment for donor age, sex, and postmortem interval (PMI)?

All outputs are **associations only**. No causal or mechanistic claim is
licensed by this design.

## Primary model

The unit of replication is the **donor**, not the cell. Let:

- `d = 1, ..., n_donors` index donors;
- `k = 1, ..., K` index niches detected in the frozen per-section pipeline;
- `y[d, k]` be donor `d`'s abundance fraction for niche `k`, computed as the
  fraction of that donor's analyzed cells assigned to niche `k`;
- `s[d]` be the frozen donor pseudo-progression score;
- `z_age[d]` be age z-scored across the analyzed donors. The source metadata
  right-censor ages at `90+`; the frozen covariate table uses the documented
  lower-bound encoding `age = 90` and retains `age_censored_90plus` for audit;
- `female_c[d]` be sex coded female = 1, male = 0, then centered at the
  analyzed-cohort female proportion;
- `z_pmi[d]` be postmortem interval z-scored across the analyzed donors.

For each niche `k`, fit one Gaussian generalized additive model with identity
link:

```text
E[y[d, k]] = f_k(s[d])
           + beta_age[k]  * z_age[d]
           + beta_sex[k]  * female_c[d]
           + beta_pmi[k]  * z_pmi[d]
```

`f_k` is a cubic B-spline with **4 basis functions** on min–max scaled
pseudo-progression scores. The primary basis has no interior knots. If a
sensitivity run changes the basis size, that changed basis size is logged as
the sensitivity choice.

Estimation uses ordinary least squares on the fixed spline-plus-covariate
design matrix. Classical OLS standard errors are used for uncertainty. If the
design matrix is rank-deficient, the run fails closed rather than silently
dropping a term.

## Primary effect size and uncertainty

The primary effect for niche `k` is the fitted abundance-fraction contrast:

```text
Delta_k = fitted_mean_k(s = Q75) - fitted_mean_k(s = Q25)
```

where `Q25` and `Q75` are the observed donor-score quartiles and centered
covariates are held at their cohort means. The implementation reports:

- `effect_size` (`Delta_k`, in abundance-fraction units);
- `std_error`;
- Wald 95% confidence interval (`ci95_low`, `ci95_high`);
- nominal two-sided normal `p_value`;
- `n_donors`, spline degrees of freedom, and the exact score values used in
  the contrast.

The nominal p-value is descriptive. The inferential p-values come from the
permutation null families below.

## Null families

Both nulls preserve donor-level covariates and use the same GAM, effect-size
contrast, and donor set as the observed fit. Unless compute forces a documented
reduction before unblinding, each final null family uses **10,000 permutations**
and a fixed seed recorded with the output.

### Donor-shuffle null

1. Keep the donor-by-niche abundance matrix and covariate table fixed.
2. Permute pseudo-progression scores across donors.
3. Refit every niche GAM and recompute `Delta_k`.
4. The two-sided permutation p-value for niche `k` is
   `(1 + count(|Delta_null,k| >= |Delta_observed,k|)) / (1 + B)`.

This tests exchangeability of progression scores at the correct unit of
replication: the donor.

### Label-shuffle null

1. Within each donor/section, permute assigned niche labels across cells while
   preserving cell count and marginal niche-label counts.
2. Recompute donor-by-niche abundance fractions.
3. Keep donor scores and covariates fixed.
4. Refit every niche GAM and recompute `Delta_k` and the same permutation
   p-value formula as above.

This tests whether the result depends on meaningful niche identity rather than
the compositional abundance machinery.

## Multiple testing

For each null family separately, adjust the per-niche permutation p-values
across the `K` tested niches with the **Benjamini–Hochberg** procedure at FDR
0.05. Report the number of niches tested, the raw p-value, and the BH q-value
for every niche. Do not report a nominally significant nominal p-value as a
confirmed result unless its pre-registered permutation q-value is also below
0.05.

## Required controls and failure handling

- The implementation was fit and smoke-tested only on synthetic donor-level
  data before the analytical freeze; real-data fits must use the donor list,
  covariate table, and source hash recorded in `data/manifest.json`.
- Fail closed if donors are missing, covariates contain missing/non-finite
  values, abundance values are negative, scores are constant, or the model has
  no residual degrees of freedom.
- Record package versions, random seed, permutation count, and the SHA-256 of
  the frozen manifest with every real-data run.
- If either null family cannot be run exactly as specified, mark the affected
  analysis deferred and state the blocker; do not substitute a post hoc null.

## Explicit non-goals

- No cell-level progression model.
- No causal language.
- No selection of niches, covariates, spline complexity, or null family based
  on real-data significance.
- No claims from synthetic validation beyond computational correctness.

## Amendment A1: cross-section niche alignment (registered 2026-11-28, pre-unblinding)

**Status:** registered before any niche–progression model was fit to the
frozen donor dataset. No real-data niche abundance, effect size, or p-value
had been computed when this amendment was committed. Issue: #7 (parent plan),
implemented under #6.

### The gap this closes

The primary plan indexes niches as "niches detected in the frozen per-section
pipeline", but niche detection runs k-means independently per section, so raw
per-section niche labels are arbitrary: label 3 in one section has no
relationship to label 3 in another. A donor-level abundance table requires
labels that mean the same niche across sections. This amendment freezes the
alignment method.

### Frozen alignment procedure

1. **Per-section niches.** On the QC-passing cells of each frozen section,
   build the k-NN graph (k = 6), compute neighbor-composition vectors over
   the fixed 24-subclass vocabulary (columns present in the frozen dataset,
   in sorted order; missing subclasses are explicit zero columns), and
   cluster with k-means (k = 8, `n_init = 10`, seed 0) — the same settings
   as the completed pipeline in `scripts/analyze_seaad.py`.
2. **Per-section niche centroids.** For section `s` and local niche `j`,
   the centroid is the cell-count-weighted mean of the neighbor-composition
   vectors of the cells assigned to `j`, in the fixed subclass vocabulary.
3. **Reference section.** The frozen section with the largest QC-passing
   cell count; ties broken by lexicographically smallest section name.
4. **First-pass alignment.** Each section's 8 centroids are matched to the
   reference section's 8 centroids by solving the linear assignment problem
   (Hungarian algorithm, `scipy.optimize.linear_sum_assignment`) on the
   pairwise Euclidean centroid-distance matrix. The matching is a bijection
   between local labels and reference labels.
5. **One consensus refinement.** The global consensus centroid for aligned
   niche `k` is the cell-count-weighted mean of the aligned per-section
   centroids across all 69 sections. Every section (including the reference)
   is then re-aligned to the consensus centroids with the same Hungarian
   rule. This second mapping is final; no further iterations.
6. **Abundance.** `y[d, k]` is the fraction of donor `d`'s QC-passing
   analyzed cells whose final aligned label is `k`. All other plan
   definitions (model, contrast, nulls, FDR) are unchanged.

### Determinism and failure handling

- Every step is deterministic given the frozen dataset: k-means uses a
  fixed seed, `linear_sum_assignment` is deterministic, and the reference
  choice rule is explicit.
- The run fails closed if any section yields fewer than 8 niches, if any
  niche is empty after alignment in any section (abundance still defined,
  but the alignment is flagged), if the assignment is not a bijection, or
  if any donor ends with zero analyzed cells.
- Alignment quality is reported as data description (mean within-niche
  centroid distance to the consensus per niche), not as a result.

### Null-family run parameters

Also fixed here, before unblinding: both null families use **B = 10,000
permutations** as specified in the primary plan. The donor-shuffle family
uses RNG seed 0; the label-shuffle family uses RNG seed 1. Seeds and B are
recorded in the run provenance together with the frozen-manifest SHA-256.
