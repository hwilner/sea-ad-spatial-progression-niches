"""Smoke tests for the one-choice-at-a-time sensitivity harness."""

import pytest

from analysis.sensitivity import SensitivityScenario, run_sensitivity_family
from tests.test_niche_progression_gam import _synthetic_donor_data


def test_sensitivity_harness_runs_and_logs_single_choice_variants():
    """Harness executes primary and one-choice variants with metadata."""
    abundance, scores, covariates = _synthetic_donor_data()
    base_config = {
        "n_spline_basis": 4,
        "contrast_quantiles": (0.25, 0.75),
    }
    scenarios = [
        SensitivityScenario("spline_df_5", {"n_spline_basis": 5}),
        SensitivityScenario(
            "wider_contrast", {"contrast_quantiles": (0.20, 0.80)}
        ),
    ]
    run = run_sensitivity_family(
        abundance, scores, covariates, base_config, scenarios
    )

    assert run.log["scenario"].tolist() == [
        "primary",
        "spline_df_5",
        "wider_contrast",
    ]
    assert run.log["varied_key"].tolist() == [
        "none",
        "n_spline_basis",
        "contrast_quantiles",
    ]
    assert set(run.results["scenario"]) == {
        "primary",
        "spline_df_5",
        "wider_contrast",
    }
    planted = run.results[run.results["niche"] == "planted"]
    assert len(planted) == 3
    assert (planted["effect_size"] > 0).all()
    assert (planted["varied_key"].iloc[1:] != "none").all()


def test_sensitivity_harness_rejects_multi_choice_or_noop_variants():
    """Ambiguous sensitivity scenarios fail before any model is fit."""
    abundance, scores, covariates = _synthetic_donor_data()
    base_config = {"n_spline_basis": 4, "contrast_quantiles": (0.25, 0.75)}
    with pytest.raises(ValueError, match="exactly one choice"):
        run_sensitivity_family(
            abundance,
            scores,
            covariates,
            base_config,
            [
                SensitivityScenario(
                    "two_changes",
                    {
                        "n_spline_basis": 5,
                        "contrast_quantiles": (0.20, 0.80),
                    },
                )
            ],
        )
    with pytest.raises(ValueError, match="does not change"):
        run_sensitivity_family(
            abundance,
            scores,
            covariates,
            base_config,
            [SensitivityScenario("noop", {"n_spline_basis": 4})],
        )
