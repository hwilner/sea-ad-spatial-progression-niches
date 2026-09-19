"""Synthetic point-cloud tests for niche detection (task: implement-niche-detection).

These tests exercise the same pipeline as tests/test_niches.py and are kept
under the file name referenced by the backlog task card.
"""

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from seaad_niches.niches import assign_niches
from seaad_niches.simulate import simulate_tissue
from seaad_niches.spatial import build_neighbor_graph, neighbor_composition


def test_deterministic_niche_assignment_on_synthetic_point_cloud():
    rng = np.random.default_rng(0)
    coords = rng.uniform(0, 1, size=(400, 2))
    labels = pd.Series(rng.choice(["A", "B"], size=400)).astype("category")
    adj = build_neighbor_graph(coords, k=6)
    comp = neighbor_composition(adj, labels)
    first = assign_niches(comp, n_clusters=2, method="kmeans", seed=0)
    second = assign_niches(comp, n_clusters=2, method="kmeans", seed=0)
    assert np.array_equal(first, second)


def test_graph_parameters_k_and_radius_are_exposed():
    coords = np.array([[0.0, 0.0], [0.05, 0.0], [3.0, 3.0], [3.05, 3.0]])
    adj_k = build_neighbor_graph(coords, k=1)
    adj_r = build_neighbor_graph(coords, radius=0.1)
    assert adj_k.shape == adj_r.shape == (4, 4)
    assert adj_r.nnz == 4  # two close pairs, symmetric
