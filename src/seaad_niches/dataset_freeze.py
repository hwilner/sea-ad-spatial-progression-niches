"""Freeze the Paper 1 donor-level analytical dataset.

This module turns the staged, open SEA-AD MTG MERFISH cell metadata into the
small, committable donor-level artifacts required by the pre-registered
niche-progression analysis. It performs technical QC only: donor metadata
completeness/consistency, per-section cell counts, coordinate finiteness, and
non-missing/positive cell metadata. It does not run hypothesis tests or fit
models.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

DONOR_COLUMNS = (
    "Donor ID",
    "Sex",
    "Age at Death",
    "PMI",
    "Continuous Pseudo-progression Score",
    "Section",
)
SECTION_COLUMNS = (
    "Donor ID",
    "Section",
    "Subclass",
    "Genes detected",
    "Number of spots",
    "Supertype confidence",
    "x",
    "y",
)


@dataclass(frozen=True)
class FreezeCriteria:
    """Technical inclusion and QC thresholds for the Paper 1 freeze.

    Attributes:
        min_cells_per_donor: Minimum analyzed cells required for a donor.
        min_cells_per_section: Minimum analyzed cells required for a section.
        min_genes_detected: Minimum allowed value of ``Genes detected``.
        min_spots: Minimum allowed value of ``Number of spots``.
        min_supertype_confidence: Minimum allowed mapping confidence.
        criteria_version: Version string for this frozen criterion set.
        covariate_version: Version string for the donor covariate schema.
        qc_version: Version string for the QC report schema.
    """

    min_cells_per_donor: int = 1_000
    min_cells_per_section: int = 1_000
    min_genes_detected: float = 1.0
    min_spots: float = 1.0
    min_supertype_confidence: float = 0.0
    criteria_version: str = "paper1-freeze-v1"
    covariate_version: str = "paper1-donor-covariates-v1"
    qc_version: str = "paper1-qc-v1"


def _require_columns(frame: pd.DataFrame, columns: tuple[str, ...], kind: str) -> None:
    """Raise a clear error when required metadata columns are absent."""
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{kind} metadata missing required columns: {missing}")


def parse_age_at_death(values: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Parse SEA-AD age-at-death values into numeric years and a censoring flag.

    SEA-AD records ages of 90 years and older as the right-censored string
    ``"90+"``. The Paper 1 covariate table uses the documented lower-bound
    encoding ``age = 90`` for those donors and sets
    ``age_censored_90plus = True`` so the convention is explicit.

    Args:
        values: Raw ``Age at Death`` values, one per cell or one per donor.

    Returns:
        A tuple ``(age_years, age_censored_90plus)``. ``age_years`` is numeric;
        the flag is boolean.
    """
    raw = values.astype(str).str.strip()
    censored = raw.str.contains("+", regex=False)
    numeric = pd.to_numeric(raw.str.replace("+", "", regex=False), errors="raise")
    return numeric.astype(float), censored.astype(bool)


def cell_qc_mask(
    obs: pd.DataFrame,
    criteria: FreezeCriteria | None = None,
) -> pd.Series:
    """Return the boolean mask of cells passing technical QC.

    Args:
        obs: Cell-level metadata containing :data:`SECTION_COLUMNS`.
        criteria: Freeze thresholds; defaults to :class:`FreezeCriteria`.

    Returns:
        A boolean ``pandas.Series`` aligned to ``obs``. Failing cells are
        excluded before donor/section abundance summaries are computed.
    """
    criteria = criteria or FreezeCriteria()
    _require_columns(obs, SECTION_COLUMNS, "cell QC")
    coords = obs[["x", "y"]].to_numpy(dtype=float)
    return (
        pd.Series(np.isfinite(coords).all(axis=1), index=obs.index)
        & obs["Subclass"].notna()
        & (obs["Genes detected"] >= criteria.min_genes_detected)
        & (obs["Number of spots"] >= criteria.min_spots)
        & (obs["Supertype confidence"] >= criteria.min_supertype_confidence)
    )


def _first_constant(group: pd.DataFrame, column: str, donor_id: str) -> Any:
    """Return a donor-constant value, failing closed on inconsistency."""
    values = group[column]
    if values.isna().any():
        raise ValueError(f"donor {donor_id!r} has missing values in {column!r}")
    if values.nunique(dropna=False) != 1:
        raise ValueError(
            f"donor {donor_id!r} has non-constant values in {column!r}"
        )
    return values.iloc[0]


def build_donor_covariates(
    obs: pd.DataFrame,
    criteria: FreezeCriteria | None = None,
) -> pd.DataFrame:
    """Build the donor-level Paper 1 covariate table from cell metadata.

    Args:
        obs: Cell-level metadata containing :data:`DONOR_COLUMNS`.
        criteria: Freeze thresholds; defaults to :class:`FreezeCriteria`.

    Returns:
        One row per donor, sorted by ``donor_id``. The table includes the
        model-ready columns ``age``, ``sex``, ``pmi``, and
        ``pseudo_progression_score`` plus donor QC counts and inclusion status.
    """
    criteria = criteria or FreezeCriteria()
    _require_columns(obs, DONOR_COLUMNS, "donor")

    rows: list[dict[str, Any]] = []
    grouped = obs.groupby("Donor ID", observed=True, sort=True)
    for donor_id, group in grouped:
        sex = _first_constant(group, "Sex", str(donor_id))
        age_raw = _first_constant(group, "Age at Death", str(donor_id))
        pmi = float(_first_constant(group, "PMI", str(donor_id)))
        score = float(
            _first_constant(
                group, "Continuous Pseudo-progression Score", str(donor_id)
            )
        )
        age_years, age_censored = parse_age_at_death(pd.Series([age_raw]))
        n_cells = len(group)
        reasons: list[str] = []
        if not np.isfinite(pmi):
            reasons.append("non-finite PMI")
        if not np.isfinite(score):
            reasons.append("non-finite pseudo-progression score")
        if n_cells < criteria.min_cells_per_donor:
            reasons.append(
                f"fewer than {criteria.min_cells_per_donor} analyzed cells"
            )
        rows.append(
            {
                "donor_id": str(donor_id),
                "age": float(age_years.iloc[0]),
                "age_at_death": str(age_raw),
                "age_censored_90plus": bool(age_censored.iloc[0]),
                "sex": str(sex),
                "pmi": pmi,
                "pseudo_progression_score": score,
                "n_sections": int(group["Section"].nunique()),
                "n_cells": n_cells,
                "included": not reasons,
                "exclusion_reason": "; ".join(reasons),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("no donors found in cell metadata")
    if not out["included"].any():
        raise ValueError("freeze criteria excluded every donor")
    return out


def build_section_qc(
    obs: pd.DataFrame,
    criteria: FreezeCriteria | None = None,
) -> pd.DataFrame:
    """Build per-section technical QC summaries for the Paper 1 freeze.

    Args:
        obs: Cell-level metadata containing :data:`SECTION_COLUMNS`.
        criteria: Freeze thresholds; defaults to :class:`FreezeCriteria`.

    Returns:
        One row per section with cell counts, metadata failure counts, and a
        boolean ``passes_qc`` flag. Sections are not dropped here; the report
        makes pass/fail status explicit.
    """
    criteria = criteria or FreezeCriteria()
    _require_columns(obs, SECTION_COLUMNS, "section QC")

    rows: list[dict[str, Any]] = []
    grouped = obs.groupby("Section", observed=True, sort=True)
    for section_id, group in grouped:
        donor_values = group["Donor ID"].dropna().unique()
        if len(donor_values) != 1:
            raise ValueError(f"section {section_id!r} maps to {len(donor_values)} donors")
        n_cells = len(group)
        passing_mask = cell_qc_mask(group, criteria)
        n_cells_passing = int(passing_mask.sum())
        failures = {
            "nonfinite_coordinates": int(
                (~np.isfinite(group[["x", "y"]].to_numpy(dtype=float))).any(axis=1).sum()
            ),
            "missing_subclass": int(group["Subclass"].isna().sum()),
            "low_genes_detected": int(
                (
                    group["Genes detected"].isna()
                    | (group["Genes detected"] < criteria.min_genes_detected)
                ).sum()
            ),
            "low_number_of_spots": int(
                (
                    group["Number of spots"].isna()
                    | (group["Number of spots"] < criteria.min_spots)
                ).sum()
            ),
            "low_supertype_confidence": int(
                (
                    group["Supertype confidence"].isna()
                    | (
                        group["Supertype confidence"]
                        < criteria.min_supertype_confidence
                    )
                ).sum()
            ),
        }
        reasons = []
        if n_cells_passing < criteria.min_cells_per_section:
            reasons.append(
                f"fewer than {criteria.min_cells_per_section} QC-passing cells"
            )
        rows.append(
            {
                "section_id": str(section_id),
                "donor_id": str(donor_values[0]),
                "n_cells": n_cells,
                "n_cells_passing_cell_qc": n_cells_passing,
                "n_cells_excluded": n_cells - n_cells_passing,
                "n_subclasses": int(group["Subclass"].nunique()),
                **failures,
                "passes_qc": not reasons,
                "qc_failure_reason": "; ".join(reasons),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("no sections found in cell metadata")
    return out


def summarize_freeze(
    donor_table: pd.DataFrame,
    section_qc: pd.DataFrame,
    criteria: FreezeCriteria | None = None,
) -> dict[str, Any]:
    """Summarize donor and section QC outcomes for reports and manifests.

    Args:
        donor_table: Output of :func:`build_donor_covariates`.
        section_qc: Output of :func:`build_section_qc`.
        criteria: Freeze thresholds; defaults to :class:`FreezeCriteria`.

    Returns:
        JSON-serializable summary with donor/section counts and explicit
        inclusion and pass/fail lists.
    """
    criteria = criteria or FreezeCriteria()
    included = donor_table[donor_table["included"]]
    excluded = donor_table[~donor_table["included"]]
    passed = section_qc[section_qc["passes_qc"]]
    failed = section_qc[~section_qc["passes_qc"]]
    return {
        "criteria_version": criteria.criteria_version,
        "covariate_version": criteria.covariate_version,
        "qc_version": criteria.qc_version,
        "thresholds": {
            "min_cells_per_donor": criteria.min_cells_per_donor,
            "min_cells_per_section": criteria.min_cells_per_section,
            "min_genes_detected": criteria.min_genes_detected,
            "min_spots": criteria.min_spots,
            "min_supertype_confidence": criteria.min_supertype_confidence,
        },
        "n_donors_total": len(donor_table),
        "n_donors_included": len(included),
        "n_donors_excluded": len(excluded),
        "included_donors": included["donor_id"].tolist(),
        "excluded_donors": excluded[["donor_id", "exclusion_reason"]].to_dict("records"),
        "n_sections_total": len(section_qc),
        "n_sections_passed": len(passed),
        "n_sections_failed": len(failed),
        "failed_sections": failed[["section_id", "qc_failure_reason"]].to_dict("records"),
        "n_cells_total": int(section_qc["n_cells"].sum()),
        "n_cells_passing_cell_qc": int(section_qc["n_cells_passing_cell_qc"].sum()),
        "n_cells_in_passing_sections": int(
            passed["n_cells_passing_cell_qc"].sum()
        ),
    }


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 hex digest of a file.

    Args:
        path: File path to hash.

    Returns:
        Lowercase SHA-256 hex digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def build_freeze_manifest(
    donor_table: pd.DataFrame,
    section_qc: pd.DataFrame,
    source_manifest: dict[str, Any],
    artifact_hashes: dict[str, str],
    criteria: FreezeCriteria | None = None,
    created_utc: str | None = None,
) -> dict[str, Any]:
    """Build the versioned, hashed Paper 1 dataset manifest.

    Args:
        donor_table: Output of :func:`build_donor_covariates`.
        section_qc: Output of :func:`build_section_qc`.
        source_manifest: Cache ``MANIFEST.json`` for the source h5ad.
        artifact_hashes: Mapping from committed artifact path to SHA-256.
        criteria: Freeze thresholds; defaults to :class:`FreezeCriteria`.
        created_utc: Optional fixed timestamp for deterministic tests.

    Returns:
        JSON-serializable manifest. It contains donor lists, artifact hashes,
        versions, and source-data provenance, but no raw cell-level data.
    """
    criteria = criteria or FreezeCriteria()
    summary = summarize_freeze(donor_table, section_qc, criteria)
    return {
        "manifest_type": "seaad_paper1_analytical_dataset",
        "manifest_version": "1.0.0",
        "created_utc": created_utc
        or datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "criteria_version": criteria.criteria_version,
        "covariate_table_version": criteria.covariate_version,
        "qc_report_version": criteria.qc_version,
        "source_dataset": {
            "dataset": source_manifest.get("dataset"),
            "key": source_manifest.get("key"),
            "source_url": source_manifest.get("source_url"),
            "registry": source_manifest.get("registry"),
            "citation": source_manifest.get("citation"),
            "file": source_manifest.get("file"),
            "size_bytes": source_manifest.get("size_bytes"),
            "sha256": source_manifest.get("sha256"),
            "access": source_manifest.get("access"),
        },
        "donors": {
            "included": summary["included_donors"],
            "excluded": summary["excluded_donors"],
            "n_included": summary["n_donors_included"],
        },
        "sections": {
            "n_total": summary["n_sections_total"],
            "n_passed": summary["n_sections_passed"],
            "failed": summary["failed_sections"],
        },
        "n_cells_in_passing_sections": summary["n_cells_in_passing_sections"],
        "artifact_sha256": artifact_hashes,
        "contains_raw_data": False,
    }


def render_qc_report(summary: dict[str, Any]) -> str:
    """Render the freeze QC summary as a short Markdown report.

    Args:
        summary: Output of :func:`summarize_freeze`.

    Returns:
        Markdown text suitable for committing under ``reports/``.
    """
    excluded = summary["excluded_donors"] or ["none"]
    failed = summary["failed_sections"] or ["none"]
    return "\n".join(
        [
            "# Paper 1 analytical dataset QC report",
            "",
            f"Criteria version: `{summary['criteria_version']}`",
            f"Covariate table version: `{summary['covariate_version']}`",
            f"QC report version: `{summary['qc_version']}`",
            "",
            "## Inclusion and technical QC thresholds",
            "",
            ("- Donors require complete, donor-constant sex, age-at-death, PMI, and "
            "Continuous Pseudo-progression Score metadata."),
            (f"- Donors require at least {summary['thresholds']['min_cells_per_donor']} "
            "analyzed cells."),
            (f"- Sections require at least {summary['thresholds']['min_cells_per_section']} "
            "cells."),
            ("- Cells must have finite x/y coordinates, a non-missing subclass, "
            f"Genes detected >= {summary['thresholds']['min_genes_detected']}, "
            f"Number of spots >= {summary['thresholds']['min_spots']}, and "
            "Supertype confidence >= "
            f"{summary['thresholds']['min_supertype_confidence']}."),
            ("- SEA-AD ages recorded as `90+` are encoded as age 90 with "
            "`age_censored_90plus = True` in the donor covariate table."),
            "",
            "## Outcome",
            "",
            (f"- Donors included: {summary['n_donors_included']} / "
            f"{summary['n_donors_total']}"),
            f"- Donors excluded: {excluded}",
            (f"- Sections passing QC: {summary['n_sections_passed']} / "
            f"{summary['n_sections_total']}"),
            f"- Sections failing QC: {failed}",
            (f"- Cells passing cell-level QC: {summary['n_cells_passing_cell_qc']} / "
            f"{summary['n_cells_total']}"),
            (f"- QC-passing cells in passing sections: "
            f"{summary['n_cells_in_passing_sections']} / "
            f"{summary['n_cells_total']}"),
            "",
            "No downstream niche-progression model is fit in this QC report.",
            "",
        ]
    )
