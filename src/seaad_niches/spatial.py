"""Spatial neighbor-graph construction and neighborhood enrichment.

All functions are pure in-memory and dataset-agnostic: they operate on
(x, y) coordinate arrays and categorical cell labels. No dataset-specific
constants are hard-coded.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def build_neighbor_graph(
    coords: np.ndarray,
    k: int = 6,
    radius: float | None = None,
) -> csr_matrix:
    """Build a symmetric neighbor graph over 2D coordinates.

    Parameters
    ----------
    coords:
    Array of shape (n_cells, 2) with x, y coordinates.
    k:
    Number of nearest neighbors per cell (mutual symmetrization).
    Ignored when ``radius`` is given.
    radius:
    If provided, connect all cell pairs within this distance
    (radius graph) instead of a k-NN graph.

    Returns:
    -------
    scipy.sparse.csr_matrix
    Symmetric binary adjacency matrix (n_cells, n_cells), zero diagonal.
    """
    from sklearn.neighbors import NearestNeighbors

    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must have shape (n_cells, 2)")
    n = coords.shape[0]
    if n == 0:
        return csr_matrix((0, 0))

    if radius is not None:
        if radius <= 0:
            raise ValueError("radius must be positive")
        nn = NearestNeighbors(radius=radius).fit(coords)
        adj = nn.radius_neighbors_graph(coords, mode="connectivity")
    else:
        if not 1 <= k < max(n, 2):
            raise ValueError("k must be >= 1 and < n_cells")
        nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
        adj = nn.kneighbors_graph(coords, mode="connectivity")

    adj = adj.maximum(adj.T).tocsr()
    adj.setdiag(0)
    adj.eliminate_zeros()
    adj.data[:] = 1.0
    return adj.astype(np.float64)


def neighbor_composition(
    adj: csr_matrix,
    labels: np.ndarray | pd.Series,
) -> pd.DataFrame:
    """Compute each cell's neighbor cell-type composition.

    Parameters
    ----------
    adj:
    Binary adjacency matrix from :func:`build_neighbor_graph`.
    labels:
    Categorical cell-type label per cell (length n_cells).

    Returns:
    -------
    pandas.DataFrame
    Rows = cells, columns = cell types, values = fraction of neighbors
    of each type (rows sum to 1 for cells with >= 1 neighbor; rows for
    isolated cells are all zero).
    """
    labels = pd.Series(np.asarray(labels)).astype("category")
    cats = labels.cat.categories
    codes = labels.cat.codes.to_numpy()
    n_types = len(cats)
    onehot = np.zeros((adj.shape[0], n_types))
    onehot[np.arange(adj.shape[0]), codes] = 1.0
    counts = adj @ onehot
    totals = counts.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        comp = np.where(totals > 0, counts / np.where(totals == 0, 1, totals), 0.0)
    return pd.DataFrame(comp, columns=list(cats))


def neighborhood_enrichment(
    adj: csr_matrix,
    labels: np.ndarray | pd.Series,
    n_permutations: int = 100,
    seed: int = 0,
) -> pd.DataFrame:
    """Z-score enrichment of cell-type co-occurrence vs label permutation.

    For each ordered pair of cell types (a, b), the observed count of
    a-b adjacencies is compared against the distribution obtained by
    permuting labels while keeping the graph fixed.

    Parameters
    ----------
    adj:
    Binary adjacency matrix.
    labels:
    Cell-type label per cell.
    n_permutations:
    Number of label permutations for the null distribution.
    seed:
    RNG seed for reproducibility.

    Returns:
    -------
    pandas.DataFrame
    Square DataFrame (types x types) of z-scores. Positive values mean
    the pair co-occurs as neighbors more often than expected by chance.
    """
    labels = pd.Series(np.asarray(labels)).astype("category")
    cats = list(labels.cat.categories)
    codes = labels.cat.codes.to_numpy().astype(np.int64)
    n_types = len(cats)
    rng = np.random.default_rng(seed)

    rows, cols = adj.nonzero()
    keep = rows < cols  # count each undirected edge once
    ei, ej = rows[keep], cols[keep]

    def pair_counts(c: np.ndarray) -> np.ndarray:
        m = np.zeros((n_types, n_types))
        a, b = c[ei], c[ej]
        np.add.at(m, (a, b), 1.0)
        np.add.at(m, (b, a), 1.0)
        return m

    obs = pair_counts(codes)
    null = np.empty((n_permutations, n_types, n_types))
    for i in range(n_permutations):
        null[i] = pair_counts(rng.permutation(codes))
    mu = null.mean(axis=0)
    sd = null.std(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > 0, (obs - mu) / np.where(sd == 0, 1, sd), 0.0)
    return pd.DataFrame(z, index=cats, columns=cats)
