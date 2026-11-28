"""Pre-registered permutation null families for the niche-progression GAM.

Implements the two null families frozen in ``docs/ANALYSIS_PLAN.md``
(donor shuffle and label shuffle), the plus-one permutation p-value, and
Benjamini-Hochberg FDR adjustment across niches. Both families refit the
exact pre-registered model and contrast; the donor-level batch solves are
algebraically identical to
:func:`seaad_niches.niche_progression.fit_niche_progression_gam` (the
equivalence is pinned by tests). No post hoc substitutions: if a family
cannot be run exactly as specified it raises instead of degrading.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .niche_progression import (
    _contrast_vector,
    _covariate_matrix,
    _make_spline_basis,
    _require_aligned_inputs,
)


def _prepare_fixed_design(
    scores: pd.Series,
    covariates: pd.DataFrame,
    n_spline_basis: int,
    contrast_quantiles: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the permutation-invariant design pieces.

    Args:
        scores: Donor pseudo-progression scores (validated).
        covariates: Donor covariates (validated).
        n_spline_basis: Number of spline basis functions (pre-registered: 4).
        contrast_quantiles: Effect-contrast quantiles (pre-registered:
            ``(0.25, 0.75)``).

    Returns:
        A tuple ``(spline_design, covariate_block, contrast)`` where
        ``spline_design`` is ``(n_donors, n_spline_basis)``,
        ``covariate_block`` is ``(n_donors, 3)``, and ``contrast`` is the
        length-``n_spline_basis + 3`` contrast vector. The spline basis,
        its min-max scaling, and the contrast quantile values all derive
        from the score *values*, so they are invariant under permuting
        scores across donors.
    """
    spline = _make_spline_basis(scores, n_basis=n_spline_basis)
    spline_design = spline.transform(scores)
    covariate_block = _covariate_matrix(covariates)
    contrast, _, _ = _contrast_vector(spline, scores, contrast_quantiles)
    return spline_design, covariate_block, contrast


def _batch_contrasts(
    design: np.ndarray, outcomes: np.ndarray, contrast: np.ndarray
) -> np.ndarray:
    """Fit every niche GAM at once and return the effect contrasts.

    Args:
        design: ``(n_donors, n_parameters)`` design matrix.
        outcomes: ``(n_donors, n_niches)`` abundance matrix.
        contrast: Length-``n_parameters`` contrast vector.

    Returns:
        Array of length ``n_niches`` with the fitted
        ``Q75 - Q25`` abundance contrast per niche.

    Raises:
        ValueError: If the design is rank-deficient (fail closed, per the
            analysis plan).
    """
    beta, _, rank, _ = np.linalg.lstsq(design, outcomes, rcond=None)
    if int(rank) < design.shape[1]:
        raise ValueError("model design matrix is rank-deficient")
    return contrast @ beta


def observed_contrasts(
    abundance: pd.DataFrame,
    scores: pd.Series,
    covariates: pd.DataFrame,
    n_spline_basis: int = 4,
    contrast_quantiles: tuple[float, float] = (0.25, 0.75),
) -> np.ndarray:
    """Return the observed per-niche contrasts via the batch path.

    Algebraically identical to the ``effect_size`` column of
    :func:`fit_niche_progression_gam`; the equivalence is pinned by tests.
    """
    abundance, scores, covariates = _require_aligned_inputs(
        abundance, scores, covariates
    )
    spline_design, covariate_block, contrast = _prepare_fixed_design(
        scores, covariates, n_spline_basis, contrast_quantiles
    )
    design = np.column_stack([spline_design, covariate_block])
    if design.shape[0] <= design.shape[1]:
        raise ValueError("need more donors than model parameters")
    return _batch_contrasts(design, abundance.to_numpy(dtype=float), contrast)


def donor_shuffle_null(
    abundance: pd.DataFrame,
    scores: pd.Series,
    covariates: pd.DataFrame,
    n_permutations: int = 10_000,
    seed: int = 0,
    n_spline_basis: int = 4,
    contrast_quantiles: tuple[float, float] = (0.25, 0.75),
) -> np.ndarray:
    """Donor-shuffle null: permute scores across donors, refit, recompute.

    The abundance matrix and covariate table stay fixed; only the
    pseudo-progression scores are permuted across donors. This tests
    exchangeability of progression scores at the correct unit of
    replication.

    Args:
        abundance: Donor-by-niche abundance fractions.
        scores: Pseudo-progression score per donor.
        covariates: Donor covariates (``age``, ``sex``, ``pmi``).
        n_permutations: Number of permutations (pre-registered: 10,000).
        seed: RNG seed (pre-registered: 0 for this family).
        n_spline_basis: Spline basis size (pre-registered: 4).
        contrast_quantiles: Contrast quantiles (pre-registered).

    Returns:
        ``(n_permutations, n_niches)`` array of null contrasts.
    """
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")
    abundance, scores, covariates = _require_aligned_inputs(
        abundance, scores, covariates
    )
    spline_design, covariate_block, contrast = _prepare_fixed_design(
        scores, covariates, n_spline_basis, contrast_quantiles
    )
    n_donors = len(scores)
    if n_donors <= spline_design.shape[1] + covariate_block.shape[1]:
        raise ValueError("need more donors than model parameters")
    outcomes = abundance.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    nulls = np.empty((n_permutations, outcomes.shape[1]))
    for b in range(n_permutations):
        permuted_design = np.column_stack([spline_design[rng.permutation(n_donors)], covariate_block])
        nulls[b] = _batch_contrasts(permuted_design, outcomes, contrast)
    return nulls


def label_shuffle_null(
    aligned_labels: np.ndarray,
    section: pd.Series,
    donor: pd.Series,
    scores: pd.Series,
    covariates: pd.DataFrame,
    n_niches: int,
    n_permutations: int = 10_000,
    seed: int = 1,
    n_spline_basis: int = 4,
    contrast_quantiles: tuple[float, float] = (0.25, 0.75),
) -> np.ndarray:
    """Label-shuffle null: permute niche labels within sections, re-aggregate.

    Within each donor/section, assigned niche labels are permuted across
    cells, preserving cell count and marginal niche-label counts. Donor
    scores and covariates stay fixed. This tests whether any association
    depends on meaningful niche identity rather than the compositional
    abundance machinery.

    Args:
        aligned_labels: Final aligned niche label per cell.
        section: Section name per cell.
        donor: Donor ID per cell.
        scores: Donor pseudo-progression scores (indexed by donor).
        covariates: Donor covariates (indexed by donor).
        n_niches: Number of aligned niches.
        n_permutations: Number of permutations (pre-registered: 10,000).
        seed: RNG seed (pre-registered: 1 for this family).
        n_spline_basis: Spline basis size (pre-registered: 4).
        contrast_quantiles: Contrast quantiles (pre-registered).

    Returns:
        ``(n_permutations, n_niches)`` array of null contrasts.

    Raises:
        ValueError: If cells reference donors missing from ``scores`` or
            labels are out of range.
    """
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")
    aligned_labels = np.asarray(aligned_labels, dtype=int)
    section = pd.Series(np.asarray(section).astype(str))
    donor = pd.Series(np.asarray(donor).astype(str))
    n_cells = len(aligned_labels)
    if not (len(section) == len(donor) == n_cells):
        raise ValueError("labels, section, and donor must have equal length")
    if n_cells and (aligned_labels.min() < 0 or aligned_labels.max() >= n_niches):
        raise ValueError("aligned labels must lie in range(n_niches)")
    donors = sorted(pd.unique(donor))
    if set(donors) != set(scores.index):
        raise ValueError("cells must reference exactly the scored donors")
    if set(donors) != set(covariates.index):
        raise ValueError("cells must reference exactly the covariate donors")

    # Fixed design: scores and covariates never permute in this family.
    scores = scores.loc[donors]
    covariates = covariates.loc[donors]
    spline_design, covariate_block, contrast = _prepare_fixed_design(
        scores, covariates, n_spline_basis, contrast_quantiles
    )
    design = np.column_stack([spline_design, covariate_block])
    n_donors = len(donors)
    if n_donors <= design.shape[1]:
        raise ValueError("need more donors than model parameters")

    donor_code = donor.map({d: i for i, d in enumerate(donors)}).to_numpy()
    donor_totals = np.bincount(donor_code, minlength=n_donors).astype(float)
    if (donor_totals == 0).any():
        raise ValueError("a donor has zero cells")

    # Pre-group cell positions by section for within-section shuffles.
    section_codes, _ = pd.factorize(section, sort=True)
    order = np.argsort(section_codes, kind="stable")
    sorted_codes = section_codes[order]
    boundaries = np.flatnonzero(np.r_[True, sorted_codes[1:] != sorted_codes[:-1]])
    groups = np.split(order, boundaries[1:])

    rng = np.random.default_rng(seed)
    nulls = np.empty((n_permutations, n_niches))
    shuffled = np.empty(n_cells, dtype=int)
    for b in range(n_permutations):
        for group in groups:
            shuffled[group] = aligned_labels[group][rng.permutation(len(group))]
        counts = np.bincount(
            donor_code * n_niches + shuffled, minlength=n_donors * n_niches
        ).reshape(n_donors, n_niches)
        fractions = counts / donor_totals[:, None]
        nulls[b] = _batch_contrasts(design, fractions, contrast)
    return nulls


def permutation_pvalues(observed: np.ndarray, nulls: np.ndarray) -> np.ndarray:
    """Two-sided plus-one permutation p-values per niche.

    ``p_k = (1 + count(|null_k| >= |observed_k|)) / (1 + B)``, exactly as
    pre-registered.

    Args:
        observed: Length-``n_niches`` observed contrasts.
        nulls: ``(B, n_niches)`` null contrasts.

    Returns:
        Length-``n_niches`` array of p-values in ``(0, 1]``.
    """
    observed = np.asarray(observed, dtype=float)
    nulls = np.asarray(nulls, dtype=float)
    if nulls.ndim != 2 or nulls.shape[1] != observed.shape[0]:
        raise ValueError("nulls must be (B, n_niches) matching observed")
    exceed = (np.abs(nulls) >= np.abs(observed)[None, :]).sum(axis=0)
    return (1.0 + exceed) / (1.0 + nulls.shape[0])


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values across the tested niches.

    Args:
        pvalues: Per-niche p-values.

    Returns:
        Array of BH-adjusted q-values (monotone, <= 1).
    """
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or p.size == 0:
        raise ValueError("pvalues must be a non-empty 1-D array")
    if ((p < 0) | (p > 1)).any():
        raise ValueError("pvalues must lie in [0, 1]")
    n = p.size
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def summarize_null_family(
    observed: np.ndarray, nulls: np.ndarray
) -> pd.DataFrame:
    """Per-niche summary of one null family (for the committed report).

    Args:
        observed: Length-``n_niches`` observed contrasts.
        nulls: ``(B, n_niches)`` null contrasts.

    Returns:
        A ``pandas.DataFrame`` with the observed contrast, null mean,
        null SD, permutation p-value, and BH q-value per niche.
    """
    p_values = permutation_pvalues(observed, nulls)
    return pd.DataFrame(
        {
            "observed_contrast": observed,
            "null_mean": nulls.mean(axis=0),
            "null_sd": nulls.std(axis=0),
            "permutation_p": p_values,
            "bh_q": benjamini_hochberg(p_values),
            "n_permutations": nulls.shape[0],
        }
    )
