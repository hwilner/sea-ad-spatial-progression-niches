"""Data-loading interface for the SEA-AD atlas.

Real-data loading is blocked until the repository owner completes the
SEA-AD data-use agreement (DUA) / open-data access steps tracked in the
data-access backlog issue. These functions fail loudly and informatively
until then. See docs/DATA_ACCESS.md.
"""

from __future__ import annotations

from pathlib import Path

ACCESS_ERROR = (
    "SEA-AD data access is not configured for this repository. "
    "Loading real SEA-AD MERFISH/snRNA-seq data requires the project owner "
    "to complete the data-access steps documented in docs/DATA_ACCESS.md "
    "(tracked in the setup-data-access-env backlog issue). "
    "Until then, use seaad_niches.simulate.simulate_tissue for development "
    "and testing."
)


def load_merfish_cells(cache_dir: str | Path | None = None):
    """Load the SEA-AD MERFISH (MTG) cell table with x/y coordinates.

    Requires SEA-AD data access; raises RuntimeError otherwise.
    Returns an anndata.AnnData or pandas.DataFrame once implemented.
    """
    _require_access(cache_dir)


def load_snrnaseq_reference(cache_dir: str | Path | None = None):
    """Load the SEA-AD snRNA-seq reference used for cell-type mapping.

    Requires SEA-AD data access; raises RuntimeError otherwise.
    """
    _require_access(cache_dir)


def _require_access(cache_dir: str | Path | None) -> None:
    manifest = None if cache_dir is None else Path(cache_dir) / "MANIFEST.json"
    if manifest is None or not manifest.exists():
        raise RuntimeError(ACCESS_ERROR)
    raise NotImplementedError(
        "Manifest found but the real-data loader is not implemented yet; "
        "this is part of the frozen-dataset backlog task."
    )
