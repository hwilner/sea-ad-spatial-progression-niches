"""Data-loading interface for the SEA-AD atlas.

The SEA-AD MTG MERFISH dataset is *open access* (AWS Open Data Registry,
``sea-ad-spatial-transcriptomics`` bucket; no DUA required). Before the
data are staged locally these functions fail loudly and informatively.
Stage the data with::

python scripts/download_seaad.py --cache-dir data/seaad

which downloads the public h5ad, verifies its checksum, and writes a
``MANIFEST.json`` into the cache directory. See docs/DATA_ACCESS.md.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ACCESS_ERROR = (
    "SEA-AD data access is not configured for this repository. "
    "Loading real SEA-AD MERFISH/snRNA-seq data requires the project owner "
    "to complete the data-access steps documented in docs/DATA_ACCESS.md "
    "(tracked in the setup-data-access-env backlog issue). "
    "Until then, use seaad_niches.simulate.simulate_tissue for development "
    "and testing."
)

# Candidate column names in the SEA-AD MERFISH h5ad ``obs`` table. The
# exact spelling has varied between releases, so resolution is tolerant.
_X_CANDIDATES = ("x", "center_x", "x_coord", "x_um", "global_x")
_Y_CANDIDATES = ("y", "center_y", "y_coord", "y_um", "global_y")
# obsm keys holding (x, y) per-section coordinates when obs lacks them.
_OBSM_XY_CANDIDATES = ("X_spatial_raw", "spatial", "X_spatial")
_DONOR_CANDIDATES = ("Donor ID", "donor", "donor_label", "Donor", "donor_id",
                     "donor_alias")
_SECTION_CANDIDATES = (
    "Section",
    "section",
    "brain_section_label",
    "section_id",
    "sample",
    "roi",
)
_CELLTYPE_CANDIDATES = ("Subclass", "subclass", "Supertype", "supertype",
                        "cell_type", "cluster")


def _pick(columns, candidates, kind):
    for c in candidates:
        if c in columns:
            return c
    raise KeyError(
        f"Could not find a {kind} column; tried {candidates}, "
        f"available: {list(columns)[:40]}"
    )


def load_merfish_cells(
    cache_dir: str | Path | None = None,
    cell_type_level: str | None = None,
) -> pd.DataFrame:
    """Load the SEA-AD MERFISH (MTG) cell table with x/y coordinates.

    Requires the staged open dataset (``MANIFEST.json`` in ``cache_dir``,
    written by ``scripts/download_seaad.py``); raises RuntimeError
    otherwise.

    Returns:
    -------
    pandas.DataFrame
    One row per cell with columns ``cell_id``, ``donor``, ``section``,
    ``x``, ``y`` and ``cell_type`` (subclass-level annotation by
    default; pass ``cell_type_level`` to choose another obs column,
    e.g. ``"supertype"``). Coordinates are per-section microns.
    """
    cache_dir = _require_access(cache_dir)
    import json

    manifest = json.loads((cache_dir / "MANIFEST.json").read_text())
    h5ad_path = cache_dir / manifest["file"]
    obs, obsm_xy = _read_obs(h5ad_path)

    if obsm_xy is None:
        xcol = _pick(obs.columns, _X_CANDIDATES, "x coordinate")
        ycol = _pick(obs.columns, _Y_CANDIDATES, "y coordinate")
        x = pd.to_numeric(obs[xcol])
        y = pd.to_numeric(obs[ycol])
    else:
        x = pd.to_numeric(pd.Series(obsm_xy[:, 0], index=obs.index))
        y = pd.to_numeric(pd.Series(obsm_xy[:, 1], index=obs.index))
    dcol = _pick(obs.columns, _DONOR_CANDIDATES, "donor")
    scol = _pick(obs.columns, _SECTION_CANDIDATES, "section")
    tcol = cell_type_level or _pick(
        obs.columns, _CELLTYPE_CANDIDATES, "cell-type annotation"
    )

    out = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "donor": obs[dcol].astype(str),
            "section": obs[scol].astype(str),
            "x": x,
            "y": y,
            "cell_type": obs[tcol].astype(str),
        }
    )
    return out


def load_snrnaseq_reference(cache_dir: str | Path | None = None):
    """Load the SEA-AD snRNA-seq reference used for cell-type mapping.

    Requires SEA-AD data access; raises RuntimeError otherwise.
    """
    _require_access(cache_dir)
    raise NotImplementedError(
        "snRNA-seq reference loading is not implemented yet; the MERFISH "
        "cells loaded by load_merfish_cells already carry mapped "
        "subclass/supertype annotations, which is what the spatial "
        "pipeline consumes."
    )


def _require_access(cache_dir: str | Path | None) -> Path:
    manifest = None if cache_dir is None else Path(cache_dir) / "MANIFEST.json"
    if manifest is None or not manifest.exists():
        raise RuntimeError(ACCESS_ERROR)
    return Path(cache_dir)


def _read_obs(h5ad_path: Path):
    """Read ``obs`` plus (x, y) coordinates from an h5ad (backed mode).

    Returns ``(obs, obsm_xy)`` where ``obsm_xy`` is an (n_cells, 2) array
    or ``None`` when coordinates live in ``obs`` columns instead.
    """
    try:
        import anndata as ad
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "Reading SEA-AD h5ad files requires the optional 'anndata' "
            "package: pip install anndata (or seaad_niches[spatial])."
        ) from e
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        obs = adata.obs.copy()
        obsm_xy = None
        if not any(c in obs.columns for c in _X_CANDIDATES):
            for key in _OBSM_XY_CANDIDATES:
                if key in adata.obsm:
                    import numpy as np

                    obsm_xy = np.asarray(adata.obsm[key])
                    break
    finally:
        adata.file.close()
    return obs, obsm_xy
