"""Synthetic tests for neighbor-graph construction and enrichment."""

import numpy as np
import pandas as pd

from seaad_niches.simulate import simulate_tissue
from seaad_niches.spatial import (
    build_neighbor_graph,
    neighbor_composition,
    neighborhood_enrichment,
)


def test_knn_graph_is_symmetric_with_expected_degrees():
    rng = np.random.default_rng(1)
    coords = rng.uniform(0, 1, size=(200, 2))
    adj = build_neighbor_graph(coords, k=6)
    assert (adj != adj.T).nnz == 0
    assert adj.diagonal().sum() == 0
    degrees = np.asarray(adj.sum(axis=1)).ravel()
    assert degrees.min() >= 6  # symmetrized k-NN: degree >= k


def test_radius_graph_connects_close_pairs():
    coords = np.array([[0.0, 0.0], [0.1, 0.0], [5.0, 5.0]])
    adj = build_neighbor_graph(coords, radius=0.5)
    assert adj[0, 1] == 1 and adj[1, 0] == 1
    assert adj[0, 2] == 0 and adj[1, 2] == 0


def test_composition_rows_are_fractions():
    coords = np.array(
        [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float
    )
    labels = pd.Series(["A", "A", "B", "B", "A", "B"])
    adj = build_neighbor_graph(coords, k=2)
    comp = neighbor_composition(adj, labels)
    assert list(comp.columns) == ["A", "B"]
    rowsum = comp.sum(axis=1)
    assert np.allclose(rowsum, 1.0)


def test_enrichment_recovers_planted_neighbor_pair():
    # Planted niches over-represent Astro-Micro co-occurrence.
    tissue = simulate_tissue(n_donors=3, cells_per_donor=800, seed=7)
    adj = build_neighbor_graph(tissue.coords, k=6)
    z = neighborhood_enrichment(adj, tissue.cell_type, n_permutations=50, seed=7)
    assert z.loc["Astro", "Micro"] > 2.5
    assert z.loc["Astro", "Micro"] == z.loc["Micro", "Astro"]
    # The depleted cross pair (Exc-Inh background is uniform) should not be
    # as strongly enriched as the planted pair.
    assert z.loc["Astro", "Micro"] > z.loc["Exc", "Inh"]


def test_enrichment_is_null_on_random_labels():
    rng = np.random.default_rng(3)
    coords = rng.uniform(0, 1, size=(1000, 2))
    labels = pd.Series(rng.choice(["A", "B", "C"], size=1000))
    adj = build_neighbor_graph(coords, k=6)
    z = neighborhood_enrichment(adj, labels, n_permutations=50, seed=3)
    assert np.abs(z.to_numpy()).max() < 4.0
