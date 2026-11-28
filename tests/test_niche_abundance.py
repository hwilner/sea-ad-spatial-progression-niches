"""Tests for seaad_niches.niche_abundance (QC-aware aggregation)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from seaad_niches.niche_abundance import (
    aggregate_frozen_niche_abundance,
    donor_niche_abundance,
    fixed_subclass_vocabulary,
    section_niche_labels,
)
from seaad_niches.simulate import simulate_tissue


def test_fixed_subclass_vocabulary_sorted_unique() -> None:
    labels = pd.Series(["b", "a", "b", "c", "a"])
    assert fixed_subclass_vocabulary(labels) == ["a", "b", "c"]
    with pytest.raises(ValueError):
        fixed_subclass_vocabulary(pd.Series([], dtype=str))


def test_donor_niche_abundance_known_answer() -> None:
    donor = pd.Series(["d1", "d1", "d1", "d1", "d2", "d2"])
    labels = np.array([0, 0, 1, 2, 1, 1])
    abundance = donor_niche_abundance(donor, labels, n_niches=3)
    assert list(abundance.columns) == ["niche_0", "niche_1", "niche_2"]
    assert abundance.loc["d1"].tolist() == pytest.approx([0.5, 0.25, 0.25])
    assert abundance.loc["d2"].tolist() == pytest.approx([0.0, 1.0, 0.0])
    assert np.allclose(abundance.sum(axis=1), 1.0)


def test_donor_niche_abundance_length_mismatch() -> None:
    with pytest.raises(ValueError, match="same length"):
        donor_niche_abundance(pd.Series(["d1"]), np.array([0, 1]), n_niches=2)


def _synthetic_cells(n_donors: int = 3, cells_per_donor: int = 300) -> pd.DataFrame:
    """Build a small two-sections-per-donor synthetic cell table."""
    tissue = simulate_tissue(
        n_donors=n_donors, cells_per_donor=cells_per_donor, seed=7
    )
    donor_ids = np.asarray(tissue.donor).astype(str)
    frames = []
    for donor in np.unique(donor_ids):
        sel = donor_ids == donor
        half = int(sel.sum()) // 2
        idx = np.flatnonzero(sel)
        for half_i, part in enumerate((idx[:half], idx[half:])):
            frames.append(
                pd.DataFrame(
                    {
                        "donor": donor,
                        "section": f"{donor}_sec{half_i}",
                        "x": tissue.coords[part, 0],
                        "y": tissue.coords[part, 1],
                        "cell_type": tissue.cell_type.iloc[part].astype(str).to_numpy(),
                    }
                )
            )
    return pd.concat(frames, ignore_index=True)


def test_section_niche_labels_shapes_and_consistency() -> None:
    cells = _synthetic_cells(n_donors=1, cells_per_donor=400)
    section = cells[cells["section"] == cells["section"].iloc[0]]
    vocab = fixed_subclass_vocabulary(cells["cell_type"])
    labels, centroids, sizes = section_niche_labels(
        section, vocab, k_neighbors=6, n_niches=4, seed=0
    )
    assert len(labels) == len(section)
    assert centroids.shape == (4, len(vocab))
    assert sizes.sum() == len(section)
    assert labels.min() >= 0 and labels.max() < 4


def test_aggregate_frozen_niche_abundance_end_to_end() -> None:
    cells = _synthetic_cells()
    # Non-positional index simulates a QC-filtered frame.
    cells = cells.iloc[::1].set_index(pd.Index(np.arange(len(cells)) * 3 + 1))
    abundance, details = aggregate_frozen_niche_abundance(
        cells, k_neighbors=6, n_niches=4, seed=0
    )
    assert abundance.shape[0] == cells["donor"].nunique()
    assert abundance.shape[1] == 4
    assert np.allclose(abundance.sum(axis=1), 1.0)
    assert (abundance.to_numpy() >= 0).all()
    assert len(details["aligned_labels"]) == len(cells)
    assert details["reference"] in set(cells["section"].unique())
    assert len(details["mappings"]) == cells["section"].nunique()


def test_aggregate_frozen_niche_abundance_missing_columns() -> None:
    with pytest.raises(ValueError, match="missing required columns"):
        aggregate_frozen_niche_abundance(pd.DataFrame({"donor": ["d1"]}))
