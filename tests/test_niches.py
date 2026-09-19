"""Synthetic tests: niche clustering recovers planted niches."""

import numpy as np
from scipy.optimize import linear_sum_assignment

from seaad_niches.niches import assign_niches
from seaad_niches.simulate import simulate_tissue
from seaad_niches.spatial import build_neighbor_graph, neighbor_composition


def _best_overlap(true_labels, pred_labels):
    true_vals = np.unique(true_labels)
    pred_vals = np.unique(pred_labels)
    m = np.zeros((len(true_vals), len(pred_vals)))
    for i, t in enumerate(true_vals):
        for j, p in enumerate(pred_vals):
            m[i, j] = np.sum((true_labels == t) & (pred_labels == p))
    r, c = linear_sum_assignment(-m)
    return m[r, c].sum() / len(true_labels)


def test_kmeans_recovers_planted_niches():
    tissue = simulate_tissue(
        n_donors=1,
        cells_per_donor=2000,
        n_niche_centers=3,
        niche_radius=0.18,
        niche_type_pairs=[("Astro", "Micro"), ("Exc", "Astro"), ("Inh", "Micro")],
        progression_slope=0.0,
        seed=11,
    )
    adj = build_neighbor_graph(tissue.coords, k=15)
    comp = neighbor_composition(adj, tissue.cell_type)
    niches = assign_niches(comp, n_clusters=4, method="kmeans", seed=11)
    assert len(niches) == len(tissue.true_niche)
    assert len(np.unique(niches)) == 4
    overlap = _best_overlap(tissue.true_niche, niches)
    assert overlap > 0.7  # deterministic planted structure is recoverable


def test_assign_niches_is_deterministic():
    tissue = simulate_tissue(n_donors=1, cells_per_donor=500, seed=5)
    adj = build_neighbor_graph(tissue.coords, k=8)
    comp = neighbor_composition(adj, tissue.cell_type)
    a = assign_niches(comp, n_clusters=3, method="kmeans", seed=42)
    b = assign_niches(comp, n_clusters=3, method="kmeans", seed=42)
    assert np.array_equal(a, b)
