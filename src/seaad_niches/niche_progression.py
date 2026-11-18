"""Covariate-adjusted GAM for donor-level niche-progression associations.

This module implements the pre-registered Paper 1 model in
``docs/ANALYSIS_PLAN.md``. The unit of replication is the donor: inputs are
one niche-abundance row per donor, one pseudo-progression score per donor,
and donor-level age, sex, and postmortem interval (PMI) covariates. No
cell-level statistics are pooled here, which avoids pseudoreplication.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import BSpline

_REQUIRED_COVARIATES = ("age", "sex", "pmi")


@dataclass(frozen=True)
class _SplineBasis:
    """Fixed cubic B-spline basis on min-max scaled scores.

    Attributes:
        lower: Lower score bound used for scaling.
        upper: Upper score bound used for scaling.
        knots: Interior knots in scaled coordinates.
        n_basis: Number of spline basis functions.
        degree: Spline degree.
    """

    lower: float
    upper: float
    knots: tuple[float, ...]
    n_basis: int
    degree: int = 3

    def transform(self, values: np.ndarray | pd.Series) -> np.ndarray:
        """Evaluate the fixed spline basis at new score values.

        Args:
            values: Pseudo-progression scores in the original scale.

        Returns:
            A dense ``(len(values), n_basis)`` design matrix.
        """
        x = np.asarray(values, dtype=float)
        if self.upper == self.lower:
            raise ValueError("pseudo-progression scores must not be constant")
        scaled = np.clip((x - self.lower) / (self.upper - self.lower), 0.0, 1.0)
        knot_vector = np.concatenate(
            [
                np.repeat(0.0, self.degree + 1),
                np.asarray(self.knots, dtype=float),
                np.repeat(1.0, self.degree + 1),
            ]
        )
        matrix = BSpline.design_matrix(
            scaled, knot_vector, self.degree, extrapolate=False
        )
        return np.asarray(matrix.todense())


def _make_spline_basis(scores: pd.Series, n_basis: int, degree: int = 3) -> _SplineBasis:
    """Build the deterministic spline basis specified by the analysis plan."""
    if n_basis < degree + 1:
        raise ValueError("n_basis must be at least degree + 1")
    values = scores.to_numpy(dtype=float)
    lower = float(np.min(values))
    upper = float(np.max(values))
    n_interior = n_basis - degree - 1
    if n_interior:
        scaled = (values - lower) / (upper - lower)
        quantiles = np.linspace(0.0, 1.0, n_interior + 2)[1:-1]
        knots = tuple(float(q) for q in np.quantile(scaled, quantiles))
    else:
        knots = ()
    return _SplineBasis(
        lower=lower, upper=upper, knots=knots, n_basis=n_basis, degree=degree
    )


def _require_aligned_inputs(
    abundance: pd.DataFrame,
    scores: pd.Series,
    covariates: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Validate and align donor-level inputs without changing their order."""
    if not isinstance(abundance, pd.DataFrame):
        raise TypeError("abundance must be a pandas DataFrame indexed by donor")
    if not isinstance(scores, pd.Series):
        raise TypeError("scores must be a pandas Series indexed by donor")
    if not isinstance(covariates, pd.DataFrame):
        raise TypeError("covariates must be a pandas DataFrame indexed by donor")
    if abundance.index.has_duplicates or scores.index.has_duplicates:
        raise ValueError("donor indices must be unique")
    if covariates.index.has_duplicates:
        raise ValueError("donor indices must be unique")
    if abundance.empty or abundance.shape[1] == 0:
        raise ValueError("abundance must contain at least one donor and one niche")
    if set(scores.index) != set(abundance.index):
        raise ValueError("scores must contain exactly the abundance donors")
    scores = scores.loc[abundance.index]
    missing = [name for name in _REQUIRED_COVARIATES if name not in covariates.columns]
    if missing:
        raise ValueError(f"covariates missing required columns: {missing}")
    if set(covariates.index) != set(abundance.index):
        raise ValueError("covariates must contain exactly the abundance donors")
    covariates = covariates.loc[abundance.index, list(_REQUIRED_COVARIATES)]

    abundance = abundance.astype(float)
    scores = scores.astype(float)
    if not np.isfinite(abundance.to_numpy()).all():
        raise ValueError("abundance values must be finite")
    if (abundance.to_numpy() < 0).any():
        raise ValueError("abundance values must be non-negative fractions")
    if not np.isfinite(scores.to_numpy()).all():
        raise ValueError("pseudo-progression scores must be finite")
    return abundance, scores, covariates


def _zscore(series: pd.Series, name: str) -> np.ndarray:
    """Return z-scores for a numeric covariate, or zeros if it is constant."""
    values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{name} must contain finite numeric values")
    sd = float(np.std(values))
    if sd == 0.0:
        return np.zeros_like(values)
    return (values - float(np.mean(values))) / sd


def _female_indicator(series: pd.Series) -> np.ndarray:
    """Encode sex as the pre-registered centered female indicator."""
    if pd.api.types.is_numeric_dtype(series):
        values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float)
        if not set(np.unique(values)).issubset({0.0, 1.0}):
            raise ValueError("numeric sex coding must use female=1 and male=0")
    else:
        mapping = {"female": 1.0, "f": 1.0, "male": 0.0, "m": 0.0}
        normalized = series.map(lambda value: str(value).strip().lower())
        invalid = sorted(set(normalized) - set(mapping))
        if invalid:
            raise ValueError(f"sex values must be female/F or male/M: {invalid}")
        values = normalized.map(mapping).to_numpy(dtype=float)
    return values - float(np.mean(values))


def _covariate_matrix(covariates: pd.DataFrame) -> np.ndarray:
    """Build the pre-registered centered covariate block."""
    age = _zscore(covariates["age"], "age")
    sex = _female_indicator(covariates["sex"])
    pmi = _zscore(covariates["pmi"], "pmi")
    return np.column_stack([age, sex, pmi])


def _design_matrix(
    scores: pd.Series,
    covariates: pd.DataFrame,
    n_spline_basis: int,
) -> tuple[np.ndarray, _SplineBasis]:
    """Construct the fixed spline-plus-covariate design matrix."""
    spline = _make_spline_basis(scores, n_basis=n_spline_basis)
    return (
        np.column_stack([spline.transform(scores), _covariate_matrix(covariates)]),
        spline,
    )


def _contrast_vector(
    spline: _SplineBasis,
    scores: pd.Series,
    contrast_quantiles: tuple[float, float],
) -> tuple[np.ndarray, float, float]:
    """Build the high-minus-low fitted-abundance contrast."""
    low_q, high_q = contrast_quantiles
    if not 0.0 <= low_q < high_q <= 1.0:
        raise ValueError("contrast_quantiles must satisfy 0 <= low < high <= 1")
    low_score = float(np.quantile(scores.to_numpy(dtype=float), low_q))
    high_score = float(np.quantile(scores.to_numpy(dtype=float), high_q))
    reference_covariates = np.zeros((2, len(_REQUIRED_COVARIATES)))
    boundary_design = np.column_stack(
        [spline.transform([low_score, high_score]), reference_covariates]
    )
    return boundary_design[1] - boundary_design[0], low_score, high_score


def fit_niche_progression_gam(
    abundance: pd.DataFrame,
    scores: pd.Series,
    covariates: pd.DataFrame,
    n_spline_basis: int = 4,
    contrast_quantiles: tuple[float, float] = (0.25, 0.75),
) -> pd.DataFrame:
    """Fit the pre-registered donor-level GAM for every niche.

    For niche ``k`` and donor ``d``, the fitted model is::

        abundance[d, k] = f_k(score[d])
            + beta_age[k] * z_age[d]
            + beta_sex[k] * centered_female[d]
            + beta_pmi[k] * z_pmi[d]
            + error[d, k]

    where ``f_k`` is a cubic B-spline with ``n_spline_basis`` basis functions
    and quantile-spaced interior knots. The reported effect size is the fitted
    abundance change from the low to the high score quantile, holding centered
    covariates at their cohort means.

    Args:
        abundance: Donor-by-niche abundance fractions. Rows are donors and
            columns are niche labels.
        scores: Pseudo-progression score per donor, indexed by donor.
        covariates: Donor-level covariates with columns ``age``, ``sex``, and
            ``pmi``. Sex accepts female/F/male/M or 1/0 coding.
        n_spline_basis: Number of cubic B-spline basis functions. The
            pre-registered primary value is 4.
        contrast_quantiles: Low and high score quantiles for the primary
            effect-size contrast. The pre-registered value is ``(0.25, 0.75)``.

    Returns:
        A ``pandas.DataFrame`` with one row per niche and columns for the
        effect size, standard error, 95% confidence interval, normal z-score,
        nominal p-value, number of donors, spline degrees of freedom, and the
        score values used in the contrast.

    Raises:
        TypeError: If an input has the wrong container type.
        ValueError: If donors, covariates, values, or residual degrees of
            freedom are invalid.
    """
    abundance, scores, covariates = _require_aligned_inputs(
        abundance, scores, covariates
    )
    x, spline = _design_matrix(scores, covariates, n_spline_basis)
    n_donors, n_parameters = x.shape
    if n_donors <= n_parameters:
        raise ValueError(
            "need more donors than model parameters "
            f"({n_donors} donors vs {n_parameters} parameters)"
        )

    contrast, low_score, high_score = _contrast_vector(
        spline, scores, contrast_quantiles
    )
    xtx_pinv = np.linalg.pinv(x.T @ x)
    rows: list[dict[str, float | int | str]] = []
    for niche in abundance.columns:
        y = abundance[niche].to_numpy(dtype=float)
        beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
        if int(rank) < n_parameters:
            raise ValueError("model design matrix is rank-deficient")
        residual = y - x @ beta
        residual_df = n_donors - int(rank)
        if residual_df <= 0:
            raise ValueError("model has no residual degrees of freedom")
        sigma2 = float(residual @ residual) / residual_df
        beta_covariance = sigma2 * xtx_pinv
        effect = float(contrast @ beta)
        variance = float(contrast @ beta_covariance @ contrast.T)
        std_error = float(np.sqrt(max(variance, 0.0)))
        if std_error > 0:
            z_score = effect / std_error
            p_value = float(2.0 * stats.norm.sf(abs(z_score)))
            ci_margin = float(stats.norm.ppf(0.975) * std_error)
            ci_low = effect - ci_margin
            ci_high = effect + ci_margin
        else:
            z_score = np.nan
            p_value = np.nan
            ci_low = np.nan
            ci_high = np.nan
        rows.append(
            {
                "niche": str(niche),
                "effect_size": effect,
                "std_error": std_error,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "z_score": z_score,
                "p_value": p_value,
                "n_donors": n_donors,
                "spline_df": spline.n_basis,
                "model_edf": int(rank),
                "contrast_low_score": low_score,
                "contrast_high_score": high_score,
            }
        )
    return pd.DataFrame(rows).set_index("niche").sort_index()
