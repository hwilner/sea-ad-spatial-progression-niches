"""Data-free spatial niche utilities for the SEA-AD progression project (Paper 1)."""

from seaad_niches.spatial import build_neighbor_graph, neighbor_composition, neighborhood_enrichment
from seaad_niches.niches import assign_niches
from seaad_niches.progression import niche_progression_association
from seaad_niches.simulate import SyntheticTissue, simulate_tissue

__all__ = [
    "build_neighbor_graph",
    "neighbor_composition",
    "neighborhood_enrichment",
    "assign_niches",
    "niche_progression_association",
    "SyntheticTissue",
    "simulate_tissue",
]

__version__ = "0.1.0"
