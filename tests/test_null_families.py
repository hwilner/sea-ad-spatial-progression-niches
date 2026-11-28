"""Tests for seaad_niches.null_families (pre-registered null machinery)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from seaad_niches.niche_progression import fit_niche_progression_gam
from seaad_niches.null_families import (
    benjamini_hochberg,
    donor_shuffle_null,
    label_shuffle_null,
    observed_contrasts,
    permutation_pvalues,
    summarize_null_family,
)


def _donor_level_data(
    n_donors: int = 16, n_niches: int = 3, slope: float = 0.3, seed: int = 3
):
    """Synthetic donor-level cohort with a planted gradient in niche_0."""
    rng = np.random.default_rng(seed)
    donors = [f"d{i:02d}" for i in range(n_donors)]
    scores = pd.Series(
        np.linspace(0.1, 0.9, n_donors), index=pd.Index(donors, name="donor_id")
    )
    covariates = pd.DataFrame(
        {
            "age": rng.normal(80, 5, n_donors),
            "sex": ["Female" if i % 3 == 0 else "Male" for i in range(n_donors)],
            "pmi": rng.normal(8, 2, n_donors),
        },
        index=pd.Index(donors, name="donor_id"),
    )
    base = rng.uniform(0.05, 0.3, (n_donors, n_niches))
    base[:, 0] += slope * scores.to_numpy()
    abundance = pd.DataFrame(
        base / base.sum(axis=1, keepdims=True),
        index=pd.Index(donors, name="donor_id"),
        columns=[f"niche_{k}" for k in range(n_niches)],
    )
    return abundance, scores, covariates


def test_permutation_pvalues_known_answers() -> None:
    nulls = np.zeros((100, 2))
    assert permutation_pvalues(np.array([0.0, 0.0]), nulls).tolist() == [1.0, 1.0]
    observed = np.array([10.0, -10.0])
    p = permutation_pvalues(observed, nulls)
    assert p.tolist() == [1 / 101, 1 / 101]
    nulls[:5, 0] = 15.0  # 5 exceedances in niche 0
    p = permutation_pvalues(observed, nulls)
    assert p[0] == pytest.approx(6 / 101)
    assert p[1] == pytest.approx(1 / 101)


def test_benjamini_hochberg_known_answers() -> None:
    p = np.array([0.001, 0.01, 0.04, 0.2])
    q = benjamini_hochberg(p)
    assert q.tolist() == pytest.approx([0.004, 0.02, 0.04 * 4 / 3, 0.2])
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)
    with pytest.raises(ValueError):
        benjamini_hochberg(np.array([1.2]))


def test_observed_contrasts_match_reference_gam() -> None:
    abundance, scores, covariates = _donor_level_data()
    gam = fit_niche_progression_gam(abundance, scores, covariates)
    batch = observed_contrasts(abundance, scores, covariates)
    assert np.allclose(gam["effect_size"].to_numpy(), batch, rtol=1e-10, atol=1e-12)


def test_donor_shuffle_null_reproducible_and_detects_planted_effect() -> None:
    # Unsaturated compositional gradient: niche_0 rises steeply with the
    # score while all fractions stay positive and rows sum to 1.
    n_donors = 16
    rng = np.random.default_rng(8)
    donors = [f"d{i:02d}" for i in range(n_donors)]
    scores = pd.Series(
        np.linspace(0.1, 0.9, n_donors), index=pd.Index(donors, name="donor_id")
    )
    covariates = pd.DataFrame(
        {
            "age": rng.normal(80, 5, n_donors),
            "sex": ["Female" if i % 3 == 0 else "Male" for i in range(n_donors)],
            "pmi": rng.normal(8, 2, n_donors),
        },
        index=pd.Index(donors, name="donor_id"),
    )
    noise = rng.normal(0, 0.005, (n_donors, 3))
    s_values = scores.to_numpy()
    abundance = pd.DataFrame(
        {
            "niche_0": 0.1 + 0.6 * s_values + noise[:, 0],
            "niche_1": 0.5 - 0.3 * s_values + noise[:, 1],
            "niche_2": 0.4 - 0.3 * s_values + noise[:, 2],
        },
        index=pd.Index(donors, name="donor_id"),
    )
    obs = observed_contrasts(abundance, scores, covariates)
    nulls_a = donor_shuffle_null(abundance, scores, covariates, n_permutations=200, seed=0)
    nulls_b = donor_shuffle_null(abundance, scores, covariates, n_permutations=200, seed=0)
    assert np.array_equal(nulls_a, nulls_b)
    assert nulls_a.shape == (200, 3)
    p = permutation_pvalues(obs, nulls_a)
    # The planted niche_0 gradient should be extreme under the donor null.
    assert p[0] <= 0.05
    # Null contrasts are centered near zero relative to the planted effect.
    assert abs(nulls_a[:, 0].mean()) < abs(obs[0])


def test_donor_shuffle_null_no_effect_is_not_tiny() -> None:
    abundance, scores, covariates = _donor_level_data(slope=0.0)
    obs = observed_contrasts(abundance, scores, covariates)
    nulls = donor_shuffle_null(abundance, scores, covariates, n_permutations=200, seed=1)
    p = permutation_pvalues(obs, nulls)
    assert p.min() > 1 / 201  # no maximally extreme false positive here


def test_label_shuffle_null_preserves_marginals_and_shape() -> None:
    rng = np.random.default_rng(11)
    n_donors, n_cells_per, n_niches = 10, 200, 3
    donors = [f"d{i:02d}" for i in range(n_donors)]
    donor = np.repeat(donors, n_cells_per)
    section = np.array([f"{d}_s{j}" for d in donors for j in range(2) for _ in range(n_cells_per // 2)])
    # Unequal marginal label counts per section make the shuffle non-trivial.
    labels = np.concatenate(
        [rng.choice(n_niches, size=n_cells_per, p=[0.6, 0.3, 0.1]) for _ in donors]
    )
    scores = pd.Series(np.linspace(0.1, 0.9, n_donors), index=pd.Index(donors, name="donor_id"))
    covariates = pd.DataFrame(
        {
            "age": rng.normal(80, 5, n_donors),
            "sex": ["Female" if i % 2 else "Male" for i in range(n_donors)],
            "pmi": rng.normal(8, 2, n_donors),
        },
        index=pd.Index(donors, name="donor_id"),
    )
    nulls = label_shuffle_null(
        labels, pd.Series(section), pd.Series(donor), scores, covariates,
        n_niches=n_niches, n_permutations=50, seed=1,
    )
    assert nulls.shape == (50, n_niches)
    assert np.isfinite(nulls).all()
    repeat = label_shuffle_null(
        labels, pd.Series(section), pd.Series(donor), scores, covariates,
        n_niches=n_niches, n_permutations=50, seed=1,
    )
    assert np.array_equal(nulls, repeat)


def test_label_shuffle_null_rejects_unknown_donors() -> None:
    scores = pd.Series([0.5], index=pd.Index(["d0"], name="donor_id"))
    covariates = pd.DataFrame(
        {"age": [80.0], "sex": ["Male"], "pmi": [8.0]},
        index=pd.Index(["d0"], name="donor_id"),
    )
    with pytest.raises(ValueError, match="scored donors"):
        label_shuffle_null(
            np.array([0, 1]), pd.Series(["s", "s"]), pd.Series(["d0", "dX"]),
            scores, covariates, n_niches=2, n_permutations=5,
        )


def test_summarize_null_family_columns() -> None:
    rng = np.random.default_rng(5)
    nulls = rng.normal(size=(100, 2))
    observed = np.array([2.0, 0.1])
    summary = summarize_null_family(observed, nulls)
    assert list(summary.columns) == [
        "observed_contrast",
        "null_mean",
        "null_sd",
        "permutation_p",
        "bh_q",
        "n_permutations",
    ]
    assert (summary["permutation_p"] > 0).all()
    assert (summary["bh_q"] >= summary["permutation_p"] - 1e-12).all()
