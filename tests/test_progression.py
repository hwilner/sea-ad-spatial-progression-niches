"""Synthetic tests: progression association detects the planted gradient."""

import numpy as np

from seaad_niches.progression import niche_progression_association
from seaad_niches.simulate import simulate_tissue


def test_detects_planted_progression_gradient():
    tissue = simulate_tissue(
        n_donors=12,
        cells_per_donor=800,
        progression_niche=1,
        progression_slope=0.8,
        seed=13,
    )
    # Use ground-truth niches directly: abundance of niche 2 (index 1)
    # scales with the donor score; other niches do not.
    res = niche_progression_association(
        tissue.true_niche,
        tissue.donor,
        tissue.progression_score,
        n_permutations=500,
        seed=13,
    )
    planted = res.loc[2]
    assert planted["spearman_r"] > 0.8
    assert planted["p_value"] < 0.05
    # The planted niche is the strongest positive association; the
    # background niche (0) anti-correlates by compositionality.
    assert res.loc[2, "spearman_r"] == res["spearman_r"].max()
    assert res.loc[0, "spearman_r"] < 0


def test_null_when_scores_shuffled_by_design():
    # Constant-score cohort: nothing to associate; r should be ~0/NaN-safe.
    tissue = simulate_tissue(n_donors=6, cells_per_donor=300, seed=2)
    scores = tissue.progression_score.copy()
    scores[:] = 0.5
    res = niche_progression_association(
        tissue.true_niche, tissue.donor, scores, n_permutations=100, seed=2
    )
    assert np.all(np.nan_to_num(res["spearman_r"].to_numpy()) == 0.0)
