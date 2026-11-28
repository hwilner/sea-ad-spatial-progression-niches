"""QC-aware donor-level niche abundance aggregation (frozen dataset).

This module turns the frozen SEA-AD per-cell table into the donor-by-niche
abundance matrix consumed by the pre-registered niche-progression GAM
(``docs/ANALYSIS_PLAN.md``). Per section it reproduces the completed
pipeline settings (k-NN graph with k = 6, neighbor composition over the
fixed subclass vocabulary, k-means with 8 niches and seed 0), aligns the
per-section labels with Amendment A1
(:mod:`seaad_niches.niche_alignment`), and aggregates aligned labels to
donor-level abundance fractions. It fits no model and computes no
significance measure.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .niche_alignment import (
    align_niche_labels,
    compute_niche_centroids,
)
from .niches import assign_niches
from .spatial import build_neighbor_graph, neighbor_composition


def fixed_subclass_vocabulary(cell_type: pd.Series) -> list[str]:
    """Return the sorted subclass vocabulary used for all composition vectors.

    Args:
        cell_type: Subclass label per QC-passing cell (all sections).

    Returns:
        Sorted list of unique subclass names. Sorting makes the column
        order independent of section processing order.
    """
    vocab = sorted(pd.unique(cell_type.astype(str)))
    if not vocab:
        raise ValueError("empty subclass vocabulary")
    return vocab


def section_niche_labels(
    section_cells: pd.DataFrame,
    vocabulary: list[str],
    k_neighbors: int = 6,
    n_niches: int = 8,
    seed: int = 0,
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """Run the frozen per-section niche pipeline on one section.

    Args:
        section_cells: QC-passing cells of one section with columns
            ``x``, ``y``, and ``cell_type``.
        vocabulary: Fixed subclass vocabulary from
            :func:`fixed_subclass_vocabulary`.
        k_neighbors: Neighbor-graph degree (pre-registered: 6).
        n_niches: Number of k-means niches (pre-registered: 8).
        seed: k-means seed (pre-registered: 0).

    Returns:
        A tuple ``(labels, centroids, sizes)``: integer local niche label
        per cell (aligned to ``section_cells`` row order), the
        ``(n_niches, len(vocabulary))`` centroid frame, and per-niche cell
        counts.

    Raises:
        ValueError: If the section is too small or yields an empty niche.
    """
    if len(section_cells) < n_niches:
        raise ValueError(
            f"section has {len(section_cells)} cells, fewer than n_niches={n_niches}"
        )
    coords = section_cells[["x", "y"]].to_numpy(dtype=float)
    adj = build_neighbor_graph(coords, k=k_neighbors)
    composition = neighbor_composition(adj, section_cells["cell_type"].astype(str))
    composition = composition.reindex(columns=vocabulary, fill_value=0.0)
    labels = assign_niches(composition, n_clusters=n_niches, method="kmeans", seed=seed)
    centroids, sizes = compute_niche_centroids(composition, labels, n_niches)
    return labels, centroids, sizes


def donor_niche_abundance(
    donor: pd.Series,
    aligned_labels: np.ndarray,
    n_niches: int,
) -> pd.DataFrame:
    """Aggregate aligned per-cell niche labels to donor abundance fractions.

    Args:
        donor: Donor ID per cell.
        aligned_labels: Final aligned niche label per cell.
        n_niches: Number of aligned niches.

    Returns:
        A ``pandas.DataFrame`` indexed by donor ID with one column per
        aligned niche (``niche_0`` ... ``niche_{n_niches - 1}``); values
        are the fraction of that donor's analyzed cells assigned to the
        niche. Each donor row sums to 1.

    Raises:
        ValueError: If any donor ends with zero analyzed cells.
    """
    donor = pd.Series(np.asarray(donor).astype(str))
    aligned_labels = np.asarray(aligned_labels, dtype=int)
    if len(donor) != len(aligned_labels):
        raise ValueError("donor and labels must have the same length")
    donors = sorted(pd.unique(donor))
    donor_code = donor.map({d: i for i, d in enumerate(donors)}).to_numpy()
    counts = np.bincount(
        donor_code * n_niches + aligned_labels,
        minlength=len(donors) * n_niches,
    ).reshape(len(donors), n_niches)
    totals = counts.sum(axis=1)
    if (totals == 0).any():
        empty = [donors[i] for i in np.flatnonzero(totals == 0)]
        raise ValueError(f"donor(s) with zero analyzed cells: {empty}")
    fractions = counts / totals[:, None]
    return pd.DataFrame(
        fractions,
        index=pd.Index(donors, name="donor_id"),
        columns=[f"niche_{k}" for k in range(n_niches)],
    )


def aggregate_frozen_niche_abundance(
    cells: pd.DataFrame,
    k_neighbors: int = 6,
    n_niches: int = 8,
    seed: int = 0,
    reference: str | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Run the full QC-aware aggregation on the frozen cell table.

    Args:
        cells: QC-passing cells with columns ``donor``, ``section``,
            ``x``, ``y``, and ``cell_type``. Callers must apply the freeze
            QC mask (:func:`seaad_niches.dataset_freeze.cell_qc_mask`)
            before calling this function.
        k_neighbors: Neighbor-graph degree (pre-registered: 6).
        n_niches: Number of niches per section (pre-registered: 8).
        seed: k-means seed (pre-registered: 0).
        reference: Optional explicit alignment reference section;
            defaults to the Amendment A1 rule.

    Returns:
        A tuple ``(abundance, details)`` where ``abundance`` is the
        donor-by-niche fraction matrix and ``details`` is a dict with the
        per-section local labels/centroids/sizes, the final label
        mappings, the consensus centroids, the reference section, and the
        fixed vocabulary — everything needed to attach per-cell aligned
        labels or audit the aggregation.

    Raises:
        ValueError: On any inconsistency (empty niches, unknown donors,
            non-finite coordinates).
    """
    required = {"donor", "section", "x", "y", "cell_type"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"cells missing required columns: {sorted(missing)}")
    if not np.isfinite(cells[["x", "y"]].to_numpy(dtype=float)).all():
        raise ValueError("non-finite coordinates in QC-passing cells")
    cells = cells.reset_index(drop=True)

    vocabulary = fixed_subclass_vocabulary(cells["cell_type"])
    per_section: dict[str, dict] = {}
    for section, frame in cells.groupby("section", sort=True):
        labels, centroids, sizes = section_niche_labels(
            frame, vocabulary, k_neighbors=k_neighbors, n_niches=n_niches, seed=seed
        )
        per_section[str(section)] = {
            "labels": labels,
            "centroids": centroids,
            "sizes": sizes,
            "index": frame.index.to_numpy(),
        }

    mappings, consensus, used_reference = align_niche_labels(
        {s: d["centroids"] for s, d in per_section.items()},
        {s: d["sizes"] for s, d in per_section.items()},
        reference=reference,
    )

    aligned = np.empty(len(cells), dtype=int)
    donor = cells["donor"].astype(str).to_numpy()
    for section, data in per_section.items():
        aligned[data["index"]] = mappings[section][data["labels"]]

    abundance = donor_niche_abundance(
        pd.Series(donor), aligned, n_niches=n_niches
    )
    details = {
        "vocabulary": vocabulary,
        "per_section": per_section,
        "mappings": mappings,
        "consensus": consensus,
        "reference": used_reference,
        "aligned_labels": aligned,
    }
    return abundance, details
