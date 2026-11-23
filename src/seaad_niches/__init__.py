"""Data-free spatial niche utilities for the SEA-AD progression project (Paper 1)."""

from seaad_niches.dataset_freeze import (
    FreezeCriteria,
    build_donor_covariates,
    build_freeze_manifest,
    build_section_qc,
    cell_qc_mask,
    summarize_freeze,
)
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
    "FreezeCriteria",
    "SyntheticTissue",
    "assign_niches",
    "build_donor_covariates",
    "build_freeze_manifest",
    "build_neighbor_graph",
    "build_section_qc",
    "cell_qc_mask",
    "fit_niche_progression_gam",
    "neighbor_composition",
    "neighborhood_enrichment",
    "niche_progression_association",
    "simulate_tissue",
    "summarize_freeze",
]

__version__ = "0.1.0"
