"""Freeze the Paper 1 donor covariates, QC report, and hashed manifest.

Reads the staged open SEA-AD MTG MERFISH h5ad in backed mode, applies the
technical criteria in :mod:`seaad_niches.dataset_freeze`, and writes small
committable artifacts under ``reports/paper1_freeze/`` plus
``data/manifest.json``. The raw h5ad remains local-only and is never committed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from seaad_niches.dataset_freeze import (
    DONOR_COLUMNS,
    SECTION_COLUMNS,
    FreezeCriteria,
    build_donor_covariates,
    build_freeze_manifest,
    build_section_qc,
    cell_qc_mask,
    render_qc_report,
    sha256_file,
    summarize_freeze,
)

_OBSM_XY_CANDIDATES = ("X_spatial_raw", "spatial", "X_spatial")
_X_CANDIDATES = ("x", "center_x", "x_coord", "x_um", "global_x")
_Y_CANDIDATES = ("y", "center_y", "y_coord", "y_um", "global_y")


def _load_freeze_obs(cache_dir: Path) -> tuple[pd.DataFrame, dict]:
    """Load only the cell metadata needed for the freeze from a staged h5ad.

    Args:
        cache_dir: Directory containing the download ``MANIFEST.json``.

    Returns:
        A tuple of ``(cell_metadata, source_manifest)`` where cell metadata has
        the donor, section, QC, and x/y columns required by the freeze module.
    """
    manifest_path = cache_dir / "MANIFEST.json"
    if not manifest_path.exists():
        raise RuntimeError(
            f"staged SEA-AD manifest not found: {manifest_path}. "
            "Run scripts/download_seaad.py first."
        )
    source_manifest = json.loads(manifest_path.read_text())
    h5ad_path = cache_dir / source_manifest["file"]

    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "The dataset freeze requires the optional 'anndata' package: "
            "pip install anndata (or seaad_niches[spatial])."
        ) from exc

    wanted = set(DONOR_COLUMNS) | set(SECTION_COLUMNS)
    wanted.discard("x")
    wanted.discard("y")
    wanted.update(_X_CANDIDATES)
    wanted.update(_Y_CANDIDATES)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        available = [column for column in adata.obs.columns if column in wanted]
        obs = adata.obs[available].copy()
        if not any(column in obs.columns for column in _X_CANDIDATES):
            for key in _OBSM_XY_CANDIDATES:
                if key in adata.obsm:
                    xy = np.asarray(adata.obsm[key])
                    obs["x"] = xy[:, 0]
                    obs["y"] = xy[:, 1]
                    break
        else:
            xcol = next(column for column in _X_CANDIDATES if column in obs.columns)
            ycol = next(column for column in _Y_CANDIDATES if column in obs.columns)
            obs = obs.rename(columns={xcol: "x", ycol: "y"})
    finally:
        adata.file.close()

    if "x" not in obs.columns or "y" not in obs.columns:
        raise KeyError(
            "could not resolve x/y coordinates from obs columns or obsm keys "
            f"{_OBSM_XY_CANDIDATES}"
        )
    return obs, source_manifest


def _repo_relative(path: Path) -> str:
    """Return a stable repository-relative path when possible."""
    try:
        return path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def main() -> int:
    """Run the Paper 1 dataset freeze and write committable artifacts.

    Returns:
        Process exit code (0 on successful freeze, 2 if QC criteria fail).
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", default="data/seaad", type=Path)
    parser.add_argument(
        "--output-dir", default=Path("reports/paper1_freeze"), type=Path
    )
    parser.add_argument(
        "--manifest-path", default=Path("data/manifest.json"), type=Path
    )
    args = parser.parse_args()

    criteria = FreezeCriteria()
    obs, source_manifest = _load_freeze_obs(args.cache_dir)
    qc_obs = obs.loc[cell_qc_mask(obs, criteria)].copy()
    donor_table = build_donor_covariates(qc_obs, criteria)
    section_qc = build_section_qc(obs, criteria)
    summary = summarize_freeze(donor_table, section_qc, criteria)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    donor_path = args.output_dir / "donor_covariates.csv"
    section_path = args.output_dir / "section_qc.csv"
    summary_path = args.output_dir / "qc_summary.json"
    report_path = args.output_dir / "qc_report.md"

    donor_table.to_csv(donor_path, index=False, float_format="%.10g")
    section_qc.to_csv(section_path, index=False, float_format="%.10g")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    report_path.write_text(render_qc_report(summary))

    artifact_hashes = {
        _repo_relative(path): sha256_file(path)
        for path in (donor_path, section_path, summary_path, report_path)
    }
    manifest = build_freeze_manifest(
        donor_table=donor_table,
        section_qc=section_qc,
        source_manifest=source_manifest,
        artifact_hashes=artifact_hashes,
        criteria=criteria,
    )
    args.manifest_path.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(
        "Paper 1 freeze: "
        f"{summary['n_donors_included']}/{summary['n_donors_total']} donors, "
        f"{summary['n_sections_passed']}/{summary['n_sections_total']} sections, "
        f"{summary['n_cells_in_passing_sections']} cells in passing sections"
    )
    print(f"Wrote {donor_path}, {section_path}, {summary_path}, {report_path}")
    print(f"Wrote {args.manifest_path}")

    if summary["n_donors_excluded"] or summary["n_sections_failed"]:
        print("QC criteria flagged exclusions; inspect the report before analysis.")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
