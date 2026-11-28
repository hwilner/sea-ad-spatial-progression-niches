"""Cross-section niche label alignment (Analysis Plan Amendment A1).

Niche detection runs k-means independently per tissue section, so raw
per-section niche labels are arbitrary across sections. This module
implements the deterministic alignment procedure frozen in
``docs/ANALYSIS_PLAN.md`` Amendment A1 (registered 2026-11-28,
pre-unblinding):

1. each section's niche centroids (cell-count-weighted mean
   neighbor-composition vectors) are matched to a reference section's
   centroids with the Hungarian algorithm;
2. a cell-count-weighted consensus centroid set is computed from the
   aligned centroids;
3. every section (including the reference) is re-aligned to the
   consensus, and that second mapping is final.

All steps are deterministic given the inputs. No real-data significance
information is used anywhere in this module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def compute_niche_centroids(
    composition: pd.DataFrame,
    labels: np.ndarray,
    n_niches: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute cell-count-weighted mean composition per niche.

    Args:
        composition: Neighbor-composition vectors (cells x subclasses),
            the output of :func:`seaad_niches.spatial.neighbor_composition`
            reindexed to the fixed subclass vocabulary.
        labels: Integer niche label per cell, values in
            ``range(n_niches)``.
        n_niches: Number of niches expected.

    Returns:
        A tuple ``(centroids, sizes)`` where ``centroids`` is a
        ``pandas.DataFrame`` of shape ``(n_niches, n_subclasses)`` whose
        row ``j`` is the mean composition of cells assigned to niche
        ``j``, and ``sizes`` is the integer cell count per niche.

    Raises:
        ValueError: If inputs are inconsistent or any niche is empty.
    """
    labels = np.asarray(labels)
    if labels.ndim != 1 or len(labels) != len(composition):
        raise ValueError("labels must be a 1-D array matching composition rows")
    if n_niches < 1:
        raise ValueError("n_niches must be >= 1")
    if labels.size and (labels.min() < 0 or labels.max() >= n_niches):
        raise ValueError("labels must lie in range(n_niches)")
    sizes = np.bincount(labels, minlength=n_niches)
    if (sizes == 0).any():
        empty = np.flatnonzero(sizes == 0).tolist()
        raise ValueError(f"empty niche(s) {empty}: cannot align degenerate labels")
    x = composition.to_numpy(dtype=float)
    if not np.isfinite(x).all():
        raise ValueError("composition values must be finite")
    centroids = np.zeros((n_niches, x.shape[1]))
    np.add.at(centroids, labels, x)
    centroids /= sizes[:, None]
    return (
        pd.DataFrame(centroids, columns=list(composition.columns)),
        sizes.astype(int),
    )


def _hungarian_mapping(
    local: pd.DataFrame, target: pd.DataFrame
) -> np.ndarray:
    """Match local centroids to target centroids by minimum-distance assignment.

    Args:
        local: Centroids (n_niches x n_subclasses) to be relabeled.
        target: Reference/consensus centroids with the same columns.

    Returns:
        ``mapping`` such that local niche ``j`` receives aligned label
        ``mapping[j]``. The mapping is always a permutation of
        ``range(n_niches)``.

    Raises:
        ValueError: If shapes differ or the assignment is not a bijection.
    """
    if local.shape != target.shape:
        raise ValueError("centroid frames must have identical shape")
    if list(local.columns) != list(target.columns):
        raise ValueError("centroid frames must share the subclass vocabulary")
    a = local.to_numpy(dtype=float)
    b = target.to_numpy(dtype=float)
    distance = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(distance)
    mapping = np.empty(local.shape[0], dtype=int)
    mapping[rows] = cols
    if sorted(mapping.tolist()) != list(range(local.shape[0])):
        raise ValueError("assignment is not a bijection")  # pragma: no cover
    return mapping


def choose_reference_section(section_cell_counts: dict[str, int]) -> str:
    """Choose the reference section per Amendment A1.

    The reference is the section with the largest QC-passing cell count;
    ties are broken by lexicographically smallest section name.

    Args:
        section_cell_counts: Mapping of section name to QC-passing cells.

    Returns:
        The reference section name.

    Raises:
        ValueError: If the mapping is empty or any count is negative.
    """
    if not section_cell_counts:
        raise ValueError("no sections provided")
    if any(int(n) < 0 for n in section_cell_counts.values()):
        raise ValueError("cell counts must be non-negative")
    return max(sorted(section_cell_counts), key=lambda s: int(section_cell_counts[s]))


def align_niche_labels(
    centroids: dict[str, pd.DataFrame],
    sizes: dict[str, np.ndarray],
    reference: str | None = None,
) -> tuple[dict[str, np.ndarray], pd.DataFrame, str]:
    """Align per-section niche labels to a common label set.

    Implements the two-pass procedure of Amendment A1: Hungarian matching
    to the reference section, then one consensus-refinement re-alignment.

    Args:
        centroids: Mapping of section name to its ``(n_niches, n_subclasses)``
            centroid frame from :func:`compute_niche_centroids`. All frames
            must share shape and columns.
        sizes: Mapping of section name to per-niche cell counts.
        reference: Optional explicit reference section; defaults to
            :func:`choose_reference_section` applied to the section cell
            totals.

    Returns:
        A tuple ``(mappings, consensus, reference)`` where ``mappings``
        maps each section name to its final label mapping (local label ->
        aligned label), ``consensus`` is the consensus centroid frame, and
        ``reference`` is the reference section actually used.

    Raises:
        ValueError: If inputs are inconsistent or the reference is unknown.
    """
    sections = sorted(centroids)
    if not sections:
        raise ValueError("no sections provided")
    if set(sizes) != set(sections):
        raise ValueError("centroids and sizes must cover the same sections")
    first = centroids[sections[0]]
    for name in sections:
        if centroids[name].shape != first.shape or list(centroids[name].columns) != list(first.columns):
            raise ValueError(f"centroid frame for section {name} is inconsistent")
        if len(sizes[name]) != first.shape[0]:
            raise ValueError(f"sizes for section {name} do not match n_niches")
        if (np.asarray(sizes[name]) <= 0).any():
            raise ValueError(f"section {name} has an empty niche")
    if reference is None:
        totals = {s: int(np.sum(sizes[s])) for s in sections}
        reference = choose_reference_section(totals)
    if reference not in centroids:
        raise ValueError(f"unknown reference section: {reference}")

    # Pass 1: align every section to the reference.
    pass1 = {
        name: _hungarian_mapping(centroids[name], centroids[reference])
        for name in sections
    }

    # Consensus: cell-count-weighted mean of pass-1-aligned centroids.
    n_niches, n_cols = first.shape
    weighted = np.zeros((n_niches, n_cols))
    weights = np.zeros(n_niches)
    for name in sections:
        c = centroids[name].to_numpy(dtype=float)
        w = np.asarray(sizes[name], dtype=float)
        for local_label in range(n_niches):
            aligned = int(pass1[name][local_label])
            weighted[aligned] += w[local_label] * c[local_label]
            weights[aligned] += w[local_label]
    if (weights <= 0).any():
        raise ValueError("consensus has an empty niche")  # pragma: no cover
    consensus = pd.DataFrame(weighted / weights[:, None], columns=list(first.columns))

    # Pass 2 (final): re-align every section, including the reference, to
    # the consensus centroids.
    final = {
        name: _hungarian_mapping(centroids[name], consensus) for name in sections
    }
    return final, consensus, reference


def alignment_quality(
    centroids: dict[str, pd.DataFrame],
    sizes: dict[str, np.ndarray],
    mappings: dict[str, np.ndarray],
    consensus: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize alignment quality as data description (not a result).

    Args:
        centroids: Per-section centroid frames.
        sizes: Per-section per-niche cell counts.
        mappings: Final mappings from :func:`align_niche_labels`.
        consensus: Consensus centroid frame.

    Returns:
        A ``pandas.DataFrame`` with one row per aligned niche: the
        cell-count-weighted mean Euclidean distance of the per-section
        centroids to the consensus, and the maximum such distance.
    """
    rows = []
    cons = consensus.to_numpy(dtype=float)
    for k in range(consensus.shape[0]):
        dists, wts = [], []
        for name, frame in centroids.items():
            c = frame.to_numpy(dtype=float)
            w = np.asarray(sizes[name], dtype=float)
            for local_label in range(frame.shape[0]):
                if int(mappings[name][local_label]) == k:
                    dists.append(float(np.linalg.norm(c[local_label] - cons[k])))
                    wts.append(float(w[local_label]))
        wts_arr = np.asarray(wts)
        rows.append(
            {
                "niche": k,
                "weighted_mean_distance": float(
                    np.average(np.asarray(dists), weights=wts_arr)
                ),
                "max_distance": float(np.max(dists)),
                "n_sections": len(dists),
            }
        )
    return pd.DataFrame(rows).set_index("niche")
