"""Associate niche abundance with pseudo-progression scores.

Donor-level analysis: for each donor, the abundance (fraction) of each
niche is correlated with the donor's pseudo-progression score. A
label-shuffle permutation test on the scores provides a null p-value.
Associations only -- no causal claims.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def niche_progression_association(
    niche_labels: np.ndarray,
    donor_ids: np.ndarray,
    scores: pd.Series,
    n_permutations: int = 1000,
    seed: int = 0,
) -> pd.DataFrame:
    """Correlate per-donor niche abundance with pseudo-progression scores.

    Parameters
    ----------
    niche_labels:
    Niche label per cell.
    donor_ids:
    Donor identifier per cell (same length as ``niche_labels``).
    scores:
    Pseudo-progression score per donor, indexed by donor id.
    n_permutations:
    Number of score permutations for the null distribution.
    seed:
    RNG seed.

    Returns:
    -------
    pandas.DataFrame
    One row per niche with columns: ``spearman_r``, ``p_value``
    (permutation-based, two-sided), and ``n_donors``. Multiple-testing
    correction is left to the caller (e.g. Benjamini-Hochberg).
    """
    niche_labels = np.asarray(niche_labels)
    donor_ids = pd.Series(np.asarray(donor_ids))
    df = pd.DataFrame({"donor": donor_ids, "niche": niche_labels})
    abundance = (
        df.groupby(["donor", "niche"]).size().unstack(fill_value=0)
    )
    abundance = abundance.div(abundance.sum(axis=1), axis=0)
    missing = set(abundance.index) - set(scores.index)
    if missing:
        raise ValueError(f"scores missing for donors: {sorted(missing)}")

    s = scores.loc[abundance.index].to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    rows = []
    for niche in abundance.columns:
        a = abundance[niche].to_numpy(dtype=float)
        r_obs = stats.spearmanr(a, s).statistic
        null = np.empty(n_permutations)
        for i in range(n_permutations):
            null[i] = stats.spearmanr(a, rng.permutation(s)).statistic
        p = (1 + np.sum(np.abs(null) >= abs(r_obs))) / (1 + n_permutations)
        rows.append(
            {
                "niche": niche,
                "spearman_r": r_obs,
                "p_value": p,
                "n_donors": len(s),
            }
        )
    return pd.DataFrame(rows).set_index("niche").sort_index()
