"""Tests for seaad_niches.niche_alignment (Amendment A1)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from seaad_niches.niche_alignment import (
    align_niche_labels,
    alignment_quality,
    choose_reference_section,
    compute_niche_centroids,
)


def _composition_frame(rows: list[list[float]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["A", "B", "C"])


def test_compute_niche_centroids_weighted_means() -> None:
    comp = _composition_frame(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        ]
    )
    labels = np.array([0, 1, 1, 0])
    centroids, sizes = compute_niche_centroids(comp, labels, n_niches=2)
    assert sizes.tolist() == [2, 2]
    assert centroids.loc[0].tolist() == [0.5, 0.5, 0.0]
    assert centroids.loc[1].tolist() == [0.0, 0.5, 0.5]


def test_compute_niche_centroids_rejects_empty_niche() -> None:
    comp = _composition_frame([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    with pytest.raises(ValueError, match="empty niche"):
        compute_niche_centroids(comp, np.array([0, 0]), n_niches=2)


def test_compute_niche_centroids_rejects_out_of_range_labels() -> None:
    comp = _composition_frame([[1.0, 0.0, 0.0]])
    with pytest.raises(ValueError, match="range"):
        compute_niche_centroids(comp, np.array([5]), n_niches=2)


def test_choose_reference_section_largest_then_lexicographic() -> None:
    assert choose_reference_section({"s2": 10, "s1": 20, "s3": 5}) == "s1"
    assert choose_reference_section({"b": 10, "a": 10}) == "a"
    with pytest.raises(ValueError):
        choose_reference_section({})


def test_alignment_recovers_label_permutation() -> None:
    """Sections whose labels are permuted must map back to one label set."""
    base = _composition_frame(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    permuted = _composition_frame(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]
    )
    centroids = {"s1": base, "s2": permuted}
    sizes = {"s1": np.array([10, 10, 10]), "s2": np.array([10, 10, 10])}
    mappings, consensus, reference = align_niche_labels(
        centroids, sizes, reference="s1"
    )
    assert reference == "s1"
    # s1 niche 0 is the [1,0,0] centroid; s2's [1,0,0] row is local label 2.
    assert mappings["s2"][2] == mappings["s1"][0]
    # The mapping is a bijection.
    assert sorted(mappings["s2"].tolist()) == [0, 1, 2]
    assert consensus.shape == (3, 3)


def test_alignment_is_deterministic() -> None:
    rng = np.random.default_rng(42)
    centroids = {
        f"s{i}": _composition_frame(rng.random((4, 3)).tolist()) for i in range(4)
    }
    sizes = {f"s{i}": np.array([5, 6, 7, 8]) for i in range(4)}
    first = align_niche_labels(centroids, sizes)
    second = align_niche_labels(centroids, sizes)
    for section in centroids:
        assert (first[0][section] == second[0][section]).all()
    assert first[2] == second[2]


def test_alignment_rejects_inconsistent_frames() -> None:
    good = _composition_frame([[1.0, 0.0, 0.0]])
    bad = pd.DataFrame([[1.0, 0.0, 0.0]], columns=["A", "B", "D"])
    with pytest.raises(ValueError, match="inconsistent"):
        align_niche_labels(
            {"s1": good, "s2": bad},
            {"s1": np.array([3]), "s2": np.array([3])},
            reference="s1",
        )


def test_alignment_quality_reports_distances() -> None:
    base = _composition_frame([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    centroids = {"s1": base, "s2": base.copy()}
    sizes = {"s1": np.array([4, 4]), "s2": np.array([4, 4])}
    mappings, consensus, _ = align_niche_labels(centroids, sizes, reference="s1")
    quality = alignment_quality(centroids, sizes, mappings, consensus)
    assert quality["weighted_mean_distance"].tolist() == [0.0, 0.0]
    assert quality["n_sections"].tolist() == [2, 2]
