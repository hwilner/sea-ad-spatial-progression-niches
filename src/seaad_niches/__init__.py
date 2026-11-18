"""Data-free spatial niche utilities for the SEA-AD progression project (Paper 1)."""

from seaad_niches.niche_progression import fit_niche_progression_gam
from seaad_niches.niches import assign_niches
from seaad_niches.progression import niche_progression_association
from seaad_niches.simulate import SyntheticTissue, simulate_tissue
from seaad_niches.spatial import (
    build_neighbor_graph,
    neighbor_composition,
    neighborhood_enrichment,
)

__all__ = [
    "SyntheticTissue",
    "assign_niches",
    "build_neighbor_graph",
    "fit_niche_progression_gam",
    "neighbor_composition",
    "neighborhood_enrichment",
    "niche_progression_association",
    "simulate_tissue",
]

__version__ = "0.1.0"
