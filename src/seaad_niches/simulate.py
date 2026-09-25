"""Synthetic spatial tissue generator with planted niches and a progression gradient.

Used for data-free testing of the spatial, niche, and progression modules.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class SyntheticTissue:
    """Container for one synthetic tissue section.

    Attributes:
    ----------
    coords:
    (n_cells, 2) x, y coordinates.
    cell_type:
    Cell-type label per cell (pandas Series of category dtype).
    true_niche:
    Ground-truth planted niche label per cell (0 = background).
    donor:
    Donor id per cell.
    progression_score:
    Donor-level pseudo-progression score (pandas Series indexed by donor).
    """

    coords: np.ndarray
    cell_type: pd.Series
    true_niche: np.ndarray
    donor: np.ndarray
    progression_score: pd.Series


def simulate_tissue(
    n_donors: int = 8,
    cells_per_donor: int = 600,
    cell_types: tuple[str, ...] = ("Exc", "Inh", "Astro", "Micro"),
    planted_niche_types: tuple[str, str] = ("Astro", "Micro"),
    niche_type_pairs: list[tuple[str, str]] | None = None,
    n_niche_centers: int = 3,
    niche_radius: float = 0.12,
    niche_fraction: float = 0.35,
    progression_niche: int = 1,
    progression_slope: float = 0.6,
    seed: int = 0,
) -> SyntheticTissue:
    """Simulate tissues with planted niches and a donor progression gradient.

    Each donor's tissue is a unit square. A fraction of cells is drawn
    around random niche centers; inside a niche, the two
    ``planted_niche_types`` are over-represented. Niche
    ``progression_niche``'s planted area fraction scales linearly with the
    donor's pseudo-progression score, giving a planted abundance gradient.

    Parameters
    ----------
    n_donors, cells_per_donor:
    Cohort and tissue size.
    cell_types:
    Background cell-type vocabulary (sampled uniformly).
    planted_niche_types:
    The two cell types enriched inside planted niches.
    niche_type_pairs:
    Optional per-niche cell-type signature pairs (cycled over niche
    centers). If None, every niche uses ``planted_niche_types``.
    n_niche_centers, niche_radius:
    Number of niche centers per tissue and their radius (in unit-square
    coordinates).
    niche_fraction:
    Fraction of cells assigned to planted niches at score = 0.5.
    progression_niche:
    Which planted niche index scales with the progression score.
    progression_slope:
    Strength of the abundance-vs-score relationship.
    seed:
    RNG seed.
    """
    rng = np.random.default_rng(seed)
    cell_types = list(cell_types)
    pairs = niche_type_pairs or [planted_niche_types]
    for pair in pairs:
        if pair[0] not in cell_types or pair[1] not in cell_types:
            raise ValueError("niche cell types must be in cell_types")

    scores = np.linspace(0.0, 1.0, n_donors)
    coords_all, types_all, niche_all, donor_all = [], [], [], []

    for d in range(n_donors):
        score = scores[d]
        coords = rng.uniform(0, 1, size=(cells_per_donor, 2))
        labels = rng.choice(cell_types, size=cells_per_donor)
        niche = np.zeros(cells_per_donor, dtype=int)

        centers = rng.uniform(0.1, 0.9, size=(n_niche_centers, 2))
        for j, c in enumerate(centers, start=1):
            # Niche `progression_niche` grows with the donor's score.
            frac = niche_fraction
            if j - 1 == progression_niche:
                frac = np.clip(
                    niche_fraction + progression_slope * (score - 0.5), 0.02, 0.9
                )
            r = niche_radius * np.sqrt(frac / niche_fraction)
            inside = np.linalg.norm(coords - c, axis=1) < r
            niche[inside] = j
            # Over-represent the two planted types inside the niche.
            n_in = inside.sum()
            flip = rng.random(n_in) < 0.95
            pair = pairs[(j - 1) % len(pairs)]
            planted = rng.choice(pair, size=n_in)
            sub = labels[inside]
            sub[flip] = planted[flip]
            labels[inside] = sub

        coords_all.append(coords)
        types_all.append(labels)
        niche_all.append(niche)
        donor_all.append(np.full(cells_per_donor, f"donor_{d}"))

    return SyntheticTissue(
        coords=np.vstack(coords_all),
        cell_type=pd.Series(np.concatenate(types_all), dtype="category"),
        true_niche=np.concatenate(niche_all),
        donor=np.concatenate(donor_all),
        progression_score=pd.Series(
            scores, index=[f"donor_{d}" for d in range(n_donors)], name="score"
        ),
    )
