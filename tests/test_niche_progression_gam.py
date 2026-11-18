"""Synthetic tests for the pre-registered niche-progression GAM."""

import numpy as np
import pandas as pd
import pytest

from seaad_niches.niche_progression import fit_niche_progression_gam


def _synthetic_donor_data(
    n_donors: int = 28,
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Create donor-level data with a planted positive niche association."""
    donors = [f"donor_{i:02d}" for i in range(n_donors)]
    score = np.linspace(0.0, 1.0, n_donors)
    rng = np.random.default_rng(41)
    age = 68 + 9 * np.linspace(-1.0, 1.0, n_donors) + rng.normal(0, 0.4, n_donors)
    pmi = 8 + 4 * np.linspace(1.0, -1.0, n_donors) + rng.normal(0, 0.2, n_donors)
    sex = np.where(np.arange(n_donors) % 2 == 0, "F", "M")
    age_z = (age - age.mean()) / age.std()

    planted = 1.2 + 1.5 * (score - 0.5) + 0.08 * age_z + rng.normal(
        0, 0.025, n_donors
    )
    negative = 0.8 - 0.5 * (score - 0.5) + rng.normal(0, 0.025, n_donors)
    reference = np.full(n_donors, 1.1) + rng.normal(0, 0.025, n_donors)
    raw = np.column_stack([planted, negative, reference])
    abundance = pd.DataFrame(
        raw / raw.sum(axis=1, keepdims=True),
        index=donors,
        columns=["planted", "negative", "reference"],
    )
    scores = pd.Series(score, index=donors, name="pseudo_progression")
    covariates = pd.DataFrame(
        {"age": age, "sex": sex, "pmi": pmi}, index=donors
    )
    return abundance, scores, covariates


def test_gam_recovers_planted_association_with_uncertainty():
    """GAM returns a positive, uncertainty-qualified planted effect."""
    abundance, scores, covariates = _synthetic_donor_data()
    result = fit_niche_progression_gam(abundance, scores, covariates)

    planted = result.loc["planted"]
    assert planted["effect_size"] > 0.10
    assert planted["ci95_low"] > 0.0
    assert planted["p_value"] < 0.01
    assert planted["n_donors"] == len(abundance)
    assert planted["spline_df"] == 4
    assert result.loc["planted", "effect_size"] == result["effect_size"].max()


def test_gam_validates_required_covariates_and_donor_alignment():
    """Missing pre-registered covariates or donors fail closed."""
    abundance, scores, covariates = _synthetic_donor_data()
    with pytest.raises(ValueError, match="missing required columns"):
        fit_niche_progression_gam(
            abundance, scores, covariates.drop(columns=["pmi"])
        )
    with pytest.raises(ValueError, match="exactly the abundance donors"):
        fit_niche_progression_gam(
            abundance, scores.iloc[:-1], covariates
        )
