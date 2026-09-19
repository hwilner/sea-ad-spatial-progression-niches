"""Niche detection: cluster cells by neighbor cell-type composition.

Primary backend is Leiden community detection on a composition similarity
graph when the optional ``leidenalg`` + ``igraph`` dependencies are
installed; otherwise a deterministic k-means fallback (scikit-learn) is
used. Both paths are pure in-memory and dataset-agnostic.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _leiden_clusters(x: np.ndarray, n_clusters: int, seed: int) -> np.ndarray | None:
    """Leiden clustering on a k-NN graph of composition vectors, if available."""
    try:
        import igraph as ig
        import leidenalg
        from sklearn.neighbors import NearestNeighbors
    except ImportError:
        return None
    k = min(15, max(2, x.shape[0] // 10))
    nn = NearestNeighbors(n_neighbors=k).fit(x)
    adj = nn.kneighbors_graph(x, mode="distance")
    adj = adj.minimum(adj.T)
    sources, targets = adj.nonzero()
    g = ig.Graph(
        n=x.shape[0],
        edges=list(zip(sources.tolist(), targets.tolist())),
        edge_attrs={"weight": (1.0 / (1.0 + adj.data)).tolist()},
        directed=False,
    )
    part = leidenalg.find_partition(
        g,
        leidenalg.CPMVertexPartition,
        weights="weight",
        resolution_parameter=0.5,
        seed=seed,
        n_iterations=-1,
    )
    labels = np.asarray(part.membership)
    # If Leiden's count differs from the requested n_clusters we still
    # return its partition (documented behavior).
    return labels


def assign_niches(
    composition: pd.DataFrame,
    n_clusters: int,
    method: str = "auto",
    seed: int = 0,
) -> np.ndarray:
    """Assign each cell to a niche by clustering composition vectors.

    Parameters
    ----------
    composition:
        Output of :func:`seaad_niches.spatial.neighbor_composition`
        (cells x cell types).
    n_clusters:
        Number of niches. With ``method="leiden"`` this is a target only
        if the k-means fallback ends up being used; Leiden chooses its own
        partition size.
    method:
        ``"auto"`` (Leiden if importable, else k-means), ``"leiden"``
        (raise if unavailable), or ``"kmeans"``.
    seed:
        RNG seed for reproducibility.

    Returns
    -------
    numpy.ndarray
        Integer niche label per cell, length n_cells.
    """
    x = composition.to_numpy(dtype=float)
    if x.shape[0] < n_clusters:
        raise ValueError("n_clusters exceeds number of cells")

    if method in ("auto", "leiden"):
        labels = _leiden_clusters(x, n_clusters, seed)
        if labels is not None:
            return labels
        if method == "leiden":
            raise ImportError(
                "method='leiden' requires the optional 'leidenalg' and "
                "'igraph' packages; install seaad_niches[leiden] or use "
                "method='kmeans'."
            )

    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=seed)
    return km.fit_predict(x)
