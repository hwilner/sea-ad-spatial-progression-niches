"""One-choice-at-a-time sensitivity harness for niche-progression models.

The harness is intentionally metadata-first: every sensitivity run records
exactly which analysis choice changed relative to the primary configuration.
That makes it impossible to silently compare runs that differ in more than one
choice. The smoke tests exercise the driver on synthetic donor-level data; no
real SEA-AD result is produced here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd

from seaad_niches.niche_progression import fit_niche_progression_gam


@dataclass(frozen=True)
class SensitivityScenario:
    """A named sensitivity variant with exactly one changed model choice.

    Attributes:
        name: Human-readable scenario identifier.
        overrides: Mapping containing exactly one model keyword to change.
    """

    name: str
    overrides: Mapping[str, Any]


@dataclass(frozen=True)
class SensitivityRun:
    """Sensitivity results plus the one-choice change log.

    Attributes:
        results: Long-form model results with scenario metadata columns.
        log: One row per run recording the varied choice and its values.
    """

    results: pd.DataFrame
    log: pd.DataFrame


def _display_value(value: Any) -> str:
    """Return a stable, readable representation of a configuration value."""
    if isinstance(value, tuple):
        return "(" + ", ".join(str(item) for item in value) + ")"
    return str(value)


def _validate_scenarios(
    base_config: Mapping[str, Any],
    scenarios: Sequence[SensitivityScenario],
) -> pd.DataFrame:
    """Validate that each sensitivity run varies exactly one choice."""
    rows: list[dict[str, str]] = [
        {
            "scenario": "primary",
            "varied_key": "none",
            "base_value": "",
            "variant_value": "",
        }
    ]
    seen_names = {"primary"}
    for scenario in scenarios:
        if not scenario.name:
            raise ValueError("sensitivity scenarios must have non-empty names")
        if scenario.name in seen_names:
            raise ValueError(f"duplicate sensitivity scenario name: {scenario.name}")
        seen_names.add(scenario.name)
        overrides = dict(scenario.overrides)
        if len(overrides) != 1:
            raise ValueError(
                f"scenario {scenario.name!r} must vary exactly one choice; "
                f"got {len(overrides)} overrides"
            )
        key, variant = next(iter(overrides.items()))
        if key not in base_config:
            raise ValueError(f"scenario {scenario.name!r} varies unknown choice {key!r}")
        base_value = base_config[key]
        if variant == base_value:
            raise ValueError(
                f"scenario {scenario.name!r} does not change {key!r} from the primary run"
            )
        rows.append(
            {
                "scenario": scenario.name,
                "varied_key": key,
                "base_value": _display_value(base_value),
                "variant_value": _display_value(variant),
            }
        )
    return pd.DataFrame(rows)


def run_sensitivity_family(
    abundance: pd.DataFrame,
    scores: pd.Series,
    covariates: pd.DataFrame,
    base_config: Mapping[str, Any],
    scenarios: Sequence[SensitivityScenario],
) -> SensitivityRun:
    """Run the primary model and one-choice sensitivity variants.

    Args:
        abundance: Donor-by-niche abundance fractions.
        scores: Donor-level pseudo-progression scores.
        covariates: Donor-level age, sex, and PMI covariates.
        base_config: Keyword arguments for the primary
            :func:`fit_niche_progression_gam` run.
        scenarios: Sensitivity variants. Each variant must override exactly
            one key already present in ``base_config``.

    Returns:
        A :class:`SensitivityRun` containing long-form results and a run log.

    Raises:
        ValueError: If a scenario varies zero or multiple choices, varies an
            unknown choice, duplicates a name, or repeats the primary value.
    """
    base_config = dict(base_config)
    log = _validate_scenarios(base_config, scenarios)
    frames: list[pd.DataFrame] = []
    all_runs = [SensitivityScenario("primary", {})] + list(scenarios)
    for scenario in all_runs:
        config = {**base_config, **scenario.overrides}
        result = fit_niche_progression_gam(
            abundance=abundance,
            scores=scores,
            covariates=covariates,
            **config,
        ).reset_index()
        result.insert(0, "scenario", scenario.name)
        if scenario.name == "primary":
            result["varied_key"] = "none"
            result["base_value"] = ""
            result["variant_value"] = ""
        else:
            key, variant = next(iter(scenario.overrides.items()))
            result["varied_key"] = key
            result["base_value"] = _display_value(base_config[key])
            result["variant_value"] = _display_value(variant)
        frames.append(result)
    return SensitivityRun(results=pd.concat(frames, ignore_index=True), log=log)
