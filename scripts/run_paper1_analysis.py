#!/usr/bin/env python
"""Run the pre-registered Paper 1 niche-progression analysis on real data.

Executes, end to end, the analysis frozen in ``docs/ANALYSIS_PLAN.md``
(including Amendment A1):

1. verify the staged h5ad against the committed freeze manifest
   (``data/manifest.json``) and re-apply the freeze QC mask;
2. per-section niche pipeline (k = 6 neighbor graph, 24-subclass neighbor
   composition, k-means k = 8, seed 0) on QC-passing cells only;
3. Amendment A1 cross-section niche alignment;
4. donor-by-niche abundance aggregation;
5. observed covariate-adjusted GAM fits plus the donor-shuffle (B = 10,000,
   seed 0) and label-shuffle (B = 10,000, seed 1) null families with BH FDR;
6. write the small committable artifacts under ``reports/paper1_analysis/``
   with full provenance.

This script writes measured numbers. It must only be run after the plan,
the amendment, and the implementing code are committed (pre-unblinding).

Usage
-----
    python scripts/run_paper1_analysis.py --cache-dir data/seaad
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import scipy
import sklearn

from seaad_niches.dataset_freeze import FreezeCriteria, cell_qc_mask
from seaad_niches.niche_abundance import aggregate_frozen_niche_abundance
from seaad_niches.niche_alignment import alignment_quality
from seaad_niches.niche_progression import fit_niche_progression_gam
from seaad_niches.null_families import (
    donor_shuffle_null,
    label_shuffle_null,
    observed_contrasts,
    summarize_null_family,
)

N_PERMUTATIONS = 10_000  # pre-registered (Amendment A1)
DONOR_SHUFFLE_SEED = 0  # pre-registered (Amendment A1)
LABEL_SHUFFLE_SEED = 1  # pre-registered (Amendment A1)
K_NEIGHBORS = 6  # pre-registered (primary pipeline)
N_NICHES = 8  # pre-registered (primary pipeline)
KMEANS_SEED = 0  # pre-registered (primary pipeline)

_CELL_QC_COLUMNS = (
    "Donor ID",
    "Section",
    "Subclass",
    "Genes detected",
    "Number of spots",
    "Supertype confidence",
)
# Coordinate resolution mirrors scripts/freeze_paper1_dataset.py exactly.
_OBSM_XY_CANDIDATES = ("X_spatial_raw", "spatial", "X_spatial")
_X_CANDIDATES = ("x", "center_x", "x_coord", "x_um", "global_x")
_Y_CANDIDATES = ("y", "center_y", "y_coord", "y_um", "global_y")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_frozen_cells(cache_dir: Path, manifest: dict) -> pd.DataFrame:
    """Load QC-passing cells of the frozen dataset from the staged h5ad."""
    source = manifest["source_dataset"]
    h5ad_path = cache_dir / source["file"]
    if not h5ad_path.exists():
        raise RuntimeError(f"staged h5ad not found: {h5ad_path}")
    size = h5ad_path.stat().st_size
    if size != source["size_bytes"]:
        raise RuntimeError(f"h5ad size mismatch: {size} != {source['size_bytes']}")
    digest = _sha256(h5ad_path)
    if digest != source["sha256"]:
        raise RuntimeError("h5ad SHA-256 does not match the frozen manifest")

    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover
        raise ImportError("run_paper1_analysis requires anndata") from exc

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        wanted = [c for c in _CELL_QC_COLUMNS if c in adata.obs.columns]
        wanted += [c for c in _X_CANDIDATES + _Y_CANDIDATES if c in adata.obs.columns]
        obs = adata.obs[wanted].copy()
        if any(c in obs.columns for c in _X_CANDIDATES):
            xcol = next(c for c in _X_CANDIDATES if c in obs.columns)
            ycol = next(c for c in _Y_CANDIDATES if c in obs.columns)
            obs = obs.rename(columns={xcol: "x", ycol: "y"})
        else:
            for key in _OBSM_XY_CANDIDATES:
                if key in adata.obsm:
                    xy = np.asarray(adata.obsm[key])
                    obs["x"] = xy[:, 0]
                    obs["y"] = xy[:, 1]
                    break
    finally:
        adata.file.close()
    if "x" not in obs.columns or "y" not in obs.columns:
        raise RuntimeError("could not resolve x/y coordinates")

    mask = cell_qc_mask(obs, FreezeCriteria())
    cells = pd.DataFrame(
        {
            "donor": obs.loc[mask, "Donor ID"].astype(str).to_numpy(),
            "section": obs.loc[mask, "Section"].astype(str).to_numpy(),
            "x": obs.loc[mask, "x"].to_numpy(dtype=float),
            "y": obs.loc[mask, "y"].to_numpy(dtype=float),
            "cell_type": obs.loc[mask, "Subclass"].astype(str).to_numpy(),
        }
    )
    expected_donors = set(manifest["donors"]["included"])
    if set(cells["donor"].unique()) != expected_donors:
        raise RuntimeError("QC-passing donors do not match the frozen manifest")
    if len(cells) != manifest["n_cells_in_passing_sections"]:
        raise RuntimeError(
            f"QC-passing cell count {len(cells)} does not match the frozen "
            f"manifest ({manifest['n_cells_in_passing_sections']})"
        )
    return cells


def main() -> int:
    """Run the pre-registered analysis and write measured artifacts.

    Returns:
        Process exit code (0 on success).
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", default="data/seaad", type=Path)
    parser.add_argument("--freeze-manifest", default="data/manifest.json", type=Path)
    parser.add_argument(
        "--covariates", default="reports/paper1_freeze/donor_covariates.csv", type=Path
    )
    parser.add_argument(
        "--output-dir", default="reports/paper1_analysis", type=Path
    )
    args = parser.parse_args()

    started = datetime.now(timezone.utc)
    manifest = json.loads(args.freeze_manifest.read_text())
    manifest_sha = hashlib.sha256(args.freeze_manifest.read_bytes()).hexdigest()
    covariate_sha = _sha256(args.covariates)
    recorded_hashes = manifest.get("artifact_sha256", {})
    if recorded_hashes.get("reports/paper1_freeze/donor_covariates.csv") not in (
        None,
        covariate_sha,
    ):
        raise RuntimeError("covariate table does not match the frozen manifest")

    cells = _load_frozen_cells(args.cache_dir, manifest)
    print(f"Loaded {len(cells)} QC-passing cells "
          f"({cells['donor'].nunique()} donors, {cells['section'].nunique()} sections)")

    covariate_table = pd.read_csv(args.covariates)
    covariate_table = covariate_table[covariate_table["included"]]
    covariates = covariate_table.set_index("donor_id")[["age", "sex", "pmi"]]
    scores = covariate_table.set_index("donor_id")["pseudo_progression_score"]

    print("Running per-section niche pipeline + Amendment A1 alignment ...")
    abundance, details = aggregate_frozen_niche_abundance(
        cells, k_neighbors=K_NEIGHBORS, n_niches=N_NICHES, seed=KMEANS_SEED
    )
    abundance = abundance.sort_index()
    scores = scores.loc[abundance.index]
    covariates = covariates.loc[abundance.index]
    print(f"Abundance matrix: {abundance.shape[0]} donors x "
          f"{abundance.shape[1]} niches; reference section {details['reference']}")

    print("Fitting observed GAM ...")
    gam = fit_niche_progression_gam(abundance, scores, covariates)
    obs = observed_contrasts(abundance, scores, covariates)
    if not np.allclose(gam["effect_size"].to_numpy(), obs, rtol=1e-10, atol=1e-12):
        raise RuntimeError("batch and reference GAM implementations disagree")

    print(f"Running donor-shuffle null (B={N_PERMUTATIONS}, seed={DONOR_SHUFFLE_SEED}) ...")
    donor_nulls = donor_shuffle_null(
        abundance, scores, covariates,
        n_permutations=N_PERMUTATIONS, seed=DONOR_SHUFFLE_SEED,
    )
    print(f"Running label-shuffle null (B={N_PERMUTATIONS}, seed={LABEL_SHUFFLE_SEED}) ...")
    label_nulls = label_shuffle_null(
        details["aligned_labels"], cells["section"], cells["donor"],
        scores, covariates, n_niches=abundance.shape[1],
        n_permutations=N_PERMUTATIONS, seed=LABEL_SHUFFLE_SEED,
    )

    donor_summary = summarize_null_family(obs, donor_nulls)
    label_summary = summarize_null_family(obs, label_nulls)

    results = gam.copy()
    results["donor_shuffle_p"] = donor_summary["permutation_p"].to_numpy()
    results["donor_shuffle_bh_q"] = donor_summary["bh_q"].to_numpy()
    results["label_shuffle_p"] = label_summary["permutation_p"].to_numpy()
    results["label_shuffle_bh_q"] = label_summary["bh_q"].to_numpy()

    quality = alignment_quality(
        {s: d["centroids"] for s, d in details["per_section"].items()},
        {s: d["sizes"] for s, d in details["per_section"].items()},
        details["mappings"],
        details["consensus"],
    )

    mapping_rows = []
    for section, mapping in details["mappings"].items():
        sizes = details["per_section"][section]["sizes"]
        for local_label, aligned_label in enumerate(mapping):
            mapping_rows.append(
                {
                    "section": section,
                    "local_niche": local_label,
                    "aligned_niche": int(aligned_label),
                    "n_cells": int(sizes[local_label]),
                }
            )
    mappings = pd.DataFrame(mapping_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    abundance.to_csv(args.output_dir / "niche_abundance_by_donor.csv", float_format="%.10g")
    results.to_csv(args.output_dir / "niche_progression_results.csv", float_format="%.10g")
    quality.to_csv(args.output_dir / "niche_alignment_quality.csv", float_format="%.10g")
    mappings.to_csv(args.output_dir / "niche_alignment_map.csv", index=False)
    donor_summary.to_csv(args.output_dir / "donor_shuffle_null_summary.csv", float_format="%.10g")
    label_summary.to_csv(args.output_dir / "label_shuffle_null_summary.csv", float_format="%.10g")

    provenance = {
        "analysis": "paper1 niche-progression GAM + null families",
        "analysis_plan": "docs/ANALYSIS_PLAN.md (incl. Amendment A1, 2026-11-28)",
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "freeze_manifest_sha256": manifest_sha,
        "freeze_criteria_version": manifest.get("criteria_version"),
        "covariate_table_sha256": covariate_sha,
        "source_h5ad_sha256": manifest["source_dataset"]["sha256"],
        "n_cells": len(cells),
        "n_donors": int(abundance.shape[0]),
        "n_sections": int(cells["section"].nunique()),
        "n_niches": int(abundance.shape[1]),
        "k_neighbors": K_NEIGHBORS,
        "kmeans_seed": KMEANS_SEED,
        "alignment_reference_section": details["reference"],
        "n_permutations": N_PERMUTATIONS,
        "donor_shuffle_seed": DONOR_SHUFFLE_SEED,
        "label_shuffle_seed": LABEL_SHUFFLE_SEED,
        "package_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "claims": "associations only; no causal language",
    }
    (args.output_dir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(f"Wrote artifacts to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
