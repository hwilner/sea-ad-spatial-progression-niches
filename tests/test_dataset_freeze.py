"""Tests for the Paper 1 analytical-dataset freeze machinery."""

import pandas as pd
import pytest

from seaad_niches.dataset_freeze import (
    FreezeCriteria,
    build_donor_covariates,
    build_freeze_manifest,
    build_section_qc,
    cell_qc_mask,
    parse_age_at_death,
    render_qc_report,
    summarize_freeze,
)


def _freeze_obs() -> pd.DataFrame:
    """Return a tiny cell-metadata fixture with two donors and three sections."""
    return pd.DataFrame(
        {
            "Donor ID": ["D1", "D1", "D1", "D1", "D2", "D2"],
            "Sex": ["Female"] * 4 + ["Male"] * 2,
            "Age at Death": ["90+"] * 4 + ["82"] * 2,
            "PMI": [5.0] * 4 + [7.5] * 2,
            "Continuous Pseudo-progression Score": [0.2] * 4 + [0.8] * 2,
            "Section": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "Subclass": ["Astro", "Micro", "Astro", "Oligo", "Astro", "Micro"],
            "Genes detected": [10.0, 12.0, 8.0, 9.0, 10.0, 11.0],
            "Number of spots": [11.0, 13.0, 9.0, 10.0, 11.0, 12.0],
            "Supertype confidence": [0.8, 0.9, 0.7, 0.8, 0.9, 0.8],
            "x": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
            "y": [0.0, 0.0, 1.0, 1.0, 0.0, 1.0],
        }
    )


def _tiny_criteria() -> FreezeCriteria:
    """Return thresholds suitable for the six-cell test fixture."""
    return FreezeCriteria(min_cells_per_donor=1, min_cells_per_section=1)


def test_parse_age_at_death_marks_90plus_censoring():
    """Age parsing uses the documented 90+ lower-bound encoding."""
    age, censored = parse_age_at_death(pd.Series(["72", "90+", "88"]))
    assert age.tolist() == [72.0, 90.0, 88.0]
    assert censored.tolist() == [False, True, False]


def test_cell_qc_mask_excludes_only_failing_cells():
    """Cell-level QC masks failing cells without dropping whole sections."""
    obs = _freeze_obs()
    obs.loc[0, "Genes detected"] = 0.0
    mask = cell_qc_mask(obs, _tiny_criteria())
    assert mask.tolist() == [False, True, True, True, True, True]


def test_build_donor_covariates_outputs_model_ready_schema():
    """The donor table exposes age, sex, PMI, score, and inclusion metadata."""
    table = build_donor_covariates(_freeze_obs(), _tiny_criteria())
    assert {"age", "sex", "pmi", "pseudo_progression_score"}.issubset(table.columns)
    assert table["donor_id"].tolist() == ["D1", "D2"]
    assert table.loc[table["donor_id"] == "D1", "age"].item() == 90.0
    assert table.loc[table["donor_id"] == "D1", "age_censored_90plus"].item()
    assert table["included"].all()


def test_build_donor_covariates_applies_min_cell_threshold():
    """Donor inclusion is explicit rather than silently dropping rows."""
    table = build_donor_covariates(
        _freeze_obs(), FreezeCriteria(min_cells_per_donor=3)
    )
    assert table.loc[table["donor_id"] == "D1", "included"].item()
    assert not table.loc[table["donor_id"] == "D2", "included"].item()
    assert "fewer than 3" in table.loc[
        table["donor_id"] == "D2", "exclusion_reason"
    ].item()


def test_build_donor_covariates_fails_on_inconstant_score():
    """Donor-level covariates must be constant across cells from that donor."""
    obs = _freeze_obs()
    obs.loc[1, "Continuous Pseudo-progression Score"] = 0.3
    with pytest.raises(ValueError, match="non-constant"):
        build_donor_covariates(obs, _tiny_criteria())


def test_build_section_qc_counts_technical_failures():
    """Section QC reports low metadata counts without mutating the input."""
    obs = _freeze_obs()
    obs.loc[0, "Genes detected"] = 0.0
    obs.loc[1, "Supertype confidence"] = -0.1
    qc = build_section_qc(obs, _tiny_criteria())
    s1 = qc[qc["section_id"] == "S1"].iloc[0]
    assert s1["n_cells"] == 2
    assert s1["n_cells_passing_cell_qc"] == 0
    assert s1["n_cells_excluded"] == 2
    assert s1["low_genes_detected"] == 1
    assert s1["low_supertype_confidence"] == 1
    assert not s1["passes_qc"]
    assert "fewer than 1 QC-passing cells" in s1["qc_failure_reason"]


def test_build_section_qc_rejects_multi_donor_sections():
    """A section must map to exactly one donor for donor-level aggregation."""
    obs = _freeze_obs()
    obs.loc[3, "Donor ID"] = "D2"
    with pytest.raises(ValueError, match="maps to 2 donors"):
        build_section_qc(obs, _tiny_criteria())


def test_summarize_and_manifest_record_versions_donors_and_hashes():
    """The freeze summary and manifest carry the promised frozen references."""
    criteria = _tiny_criteria()
    donors = build_donor_covariates(_freeze_obs(), criteria)
    sections = build_section_qc(_freeze_obs(), criteria)
    summary = summarize_freeze(donors, sections, criteria)
    manifest = build_freeze_manifest(
        donors,
        sections,
        source_manifest={
            "dataset": "fixture",
            "key": "fixture.h5ad",
            "source_url": "https://example.test/fixture.h5ad",
            "registry": "fixture-registry",
            "citation": "fixture citation",
            "file": "fixture.h5ad",
            "size_bytes": 123,
            "sha256": "abc123",
            "access": "open fixture",
        },
        artifact_hashes={"reports/paper1_freeze/donor_covariates.csv": "def456"},
        criteria=criteria,
        created_utc="2026-11-23T00:00:00+00:00",
    )
    assert summary["n_donors_included"] == 2
    assert summary["n_sections_passed"] == 3
    assert manifest["created_utc"] == "2026-11-23T00:00:00+00:00"
    assert manifest["covariate_table_version"] == "paper1-donor-covariates-v1"
    assert manifest["qc_report_version"] == "paper1-qc-v1"
    assert manifest["donors"]["included"] == ["D1", "D2"]
    assert manifest["artifact_sha256"] == {
        "reports/paper1_freeze/donor_covariates.csv": "def456"
    }
    assert manifest["contains_raw_data"] is False


def test_render_qc_report_states_no_model_was_fit():
    """The QC report keeps data-freeze work separate from inference."""
    criteria = _tiny_criteria()
    summary = summarize_freeze(
        build_donor_covariates(_freeze_obs(), criteria),
        build_section_qc(_freeze_obs(), criteria),
        criteria,
    )
    report = render_qc_report(summary)
    assert "Donors included: 2 / 2" in report
    assert "No downstream niche-progression model is fit" in report
