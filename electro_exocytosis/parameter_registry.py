"""Auditable registries for framework models and runtime parameters.

The simulator's authoritative mechanistic defaults remain in
``parameters/default_parameters.yaml``.  This module does not introduce a
second source of numerical defaults; it turns that mapping into stable tables
that can be written beside a calibration result.  Optional effective
parameters and fitted-parameter rows can be joined to the same snapshot.

The registries are descriptive.  They deliberately do not validate parameter
overrides or change simulation behaviour.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from copy import deepcopy
from dataclasses import dataclass
import json
import math
import re
from typing import Any

import pandas as pd

from electro_exocytosis.io.readers import load_default_parameters, merge_parameters


PARAMETER_SNAPSHOT_COLUMNS = (
    "layer",
    "module",
    "submodule",
    "parameter_path",
    "qualified_name",
    "parameter",
    "default_value",
    "effective_value",
    "unit",
    "parameter_class",
    "runtime_status",
    "fit_status",
    "fit_variant",
    "status",
    "fit_initial",
    "fit_lower_bound",
    "fit_upper_bound",
    "fit_final",
    "source",
    "notes",
)


_BUILTIN_PARAMETER_METADATA: dict[str, dict[str, str]] = {
    "pulse.geometry_reference": {
        "runtime_status": "declared_not_consumed",
        "notes": "Retained in the defaults file, but the current pulse model does not read it.",
    },
    "pulse.temperature_capacity_J_m3_K": {
        "runtime_status": "declared_not_consumed",
        "notes": (
            "The current pulse model uses the exposure configuration (or the legacy "
            "4.0e6 constant), not this YAML value."
        ),
    },
    "dosimetry.cuvette_factor": {
        "runtime_status": "declared_not_consumed",
        "notes": "The current dosimetry model uses the GEOMETRY_FIELD_FACTORS code constant.",
    },
    "dosimetry.dish_factor": {
        "runtime_status": "declared_not_consumed",
        "notes": "The current dosimetry model uses the GEOMETRY_FIELD_FACTORS code constant.",
    },
    "dosimetry.flow_factor": {
        "runtime_status": "declared_not_consumed",
        "notes": "The current dosimetry model uses the GEOMETRY_FIELD_FACTORS code constant.",
    },
    "remodeling_repair.K_PS_uM": {
        "runtime_status": "declared_not_consumed",
        "notes": "Declared in RemodelingParams but not referenced by the current reduced equations.",
    },
    "remodeling_repair.tau_repair_s": {
        "runtime_status": "declared_not_consumed",
        "notes": "Declared in RemodelingParams but not referenced by the current reduced equations.",
    },
}


@dataclass(frozen=True, slots=True)
class ModelRegistryEntry:
    layer: str
    module: str
    submodule: str
    model: str
    implementation_status: str
    runtime_role: str
    default_section: str
    notes: str = ""


_MODEL_REGISTRY = (
    ModelRegistryEntry(
        "1",
        "pulse",
        "pulse descriptors",
        "waveform, train, dose, and energy descriptors",
        "implemented_preliminary",
        "core",
        "pulse",
    ),
    ModelRegistryEntry(
        "1",
        "dosimetry",
        "field geometry",
        "geometry-adjusted electric-field exposure",
        "implemented_preliminary",
        "core",
        "dosimetry",
        "Geometry factors are currently code constants rather than YAML parameters.",
    ),
    ModelRegistryEntry(
        "1",
        "dosimetry",
        "thermal dose",
        "legacy, adiabatic Joule, or lumped-thermal dose",
        "implemented_preliminary",
        "core",
        "dosimetry",
    ),
    ModelRegistryEntry(
        "2",
        "electrodynamics",
        "membrane charging",
        "reduced Schwan charging response",
        "implemented_preliminary",
        "core",
        "electrodynamics",
    ),
    ModelRegistryEntry(
        "2",
        "electrodynamics",
        "organelle polarization",
        "fractional ER, mitochondrial, and MVB voltages",
        "implemented_preliminary",
        "core",
        "electrodynamics",
    ),
    ModelRegistryEntry(
        "2",
        "electrodynamics",
        "electroporation",
        "phenomenological permeability and pore-density gates",
        "implemented_preliminary",
        "core",
        "electrodynamics",
    ),
    ModelRegistryEntry(
        "3",
        "ion_transport",
        "calcium influx and stores",
        "pore influx, ER release, SERCA, PMCA, and NCX",
        "implemented_preliminary",
        "core",
        "ion_transport",
    ),
    ModelRegistryEntry(
        "3",
        "ion_transport",
        "mitochondrial calcium",
        "uptake, release, mPTP, and potential recovery",
        "implemented_preliminary",
        "core",
        "ion_transport",
    ),
    ModelRegistryEntry(
        "3",
        "ion_transport",
        "monovalent ions and osmosis",
        "Na/K/Cl pore flux and osmotic recovery",
        "implemented_preliminary",
        "core",
        "ion_transport",
    ),
    ModelRegistryEntry(
        "3",
        "ion_transport",
        "metabolic stress",
        "ROS and ATP production, depletion, and recovery",
        "implemented_preliminary",
        "core",
        "ion_transport",
    ),
    ModelRegistryEntry(
        "4",
        "remodeling_repair",
        "calcium microdomains and PS",
        "local calcium, scramblase/flippase, and PS exposure",
        "implemented_preliminary",
        "core",
        "remodeling_repair",
    ),
    ModelRegistryEntry(
        "4",
        "remodeling_repair",
        "cytoskeletal remodeling",
        "calpain, actomyosin, and actin disruption",
        "implemented_preliminary",
        "core",
        "remodeling_repair",
    ),
    ModelRegistryEntry(
        "4",
        "remodeling_repair",
        "membrane repair",
        "annexin/lysosomal resealing and repair shedding",
        "implemented_preliminary",
        "core",
        "remodeling_repair",
    ),
    ModelRegistryEntry(
        "5",
        "ev_release",
        "endosomal sEV pathway",
        "MVB maturation, ILV loading, docking, fusion, and routing",
        "implemented_preliminary",
        "core",
        "ev_release",
    ),
    ModelRegistryEntry(
        "5",
        "ev_release",
        "plasma-membrane mlEV pathway",
        "budding, scission, and turnover",
        "implemented_preliminary",
        "core",
        "ev_release",
    ),
    ModelRegistryEntry(
        "5",
        "ev_release",
        "apoptotic-body pathway",
        "damage-gated commitment and resolution",
        "implemented_preliminary",
        "core",
        "ev_release",
    ),
    ModelRegistryEntry(
        "5b",
        "extracellular_kinetics",
        "extracellular stocks",
        "well-mixed sEV, mlEV, and AB source-to-stock balance",
        "implemented_preliminary",
        "core",
        "extracellular_kinetics",
    ),
    ModelRegistryEntry(
        "5b",
        "extracellular_kinetics",
        "particle loss",
        "effective loss, uptake, degradation, and adsorption",
        "implemented_preliminary",
        "core_optional",
        "extracellular_kinetics",
        "Loss terms default to zero for backward compatibility.",
    ),
    ModelRegistryEntry(
        "5b",
        "extracellular_kinetics",
        "aggregation",
        "reduced sEV-to-mlEV transfer",
        "implemented_preliminary",
        "core_optional",
        "extracellular_kinetics",
    ),
    ModelRegistryEntry(
        "5b",
        "extracellular_kinetics",
        "sampling and assay",
        "sample/replacement events, dilution, recovery, and background",
        "implemented_preliminary",
        "core_optional",
        "extracellular_kinetics",
    ),
    ModelRegistryEntry(
        "observation",
        "ev_size_observation",
        "source conversion",
        "pathway model units to particle-equivalent source scaling",
        "implemented_preliminary",
        "observation_calibration",
        "",
        "Parameters currently live in the size-resolved calibration pipeline, outside the authoritative defaults YAML.",
    ),
    ModelRegistryEntry(
        "observation",
        "ev_size_observation",
        "pathway-to-size kernel",
        "state-conditioned lognormal pathway kernels",
        "implemented_preliminary",
        "observation_calibration",
        "",
        "Parameters currently live in the size-resolved calibration pipeline, outside the authoritative defaults YAML.",
    ),
    ModelRegistryEntry(
        "observation",
        "ev_size_observation",
        "size-dependent loss",
        "smooth diameter-dependent extracellular loss",
        "implemented_preliminary",
        "observation_calibration",
        "",
        "Parameters currently live in the size-resolved calibration pipeline.",
    ),
    ModelRegistryEntry(
        "observation",
        "ev_size_observation",
        "state adapter",
        "pathway-specific cell-state shifts of size-kernel medians",
        "implemented_preliminary",
        "observation_calibration_optional",
        "",
        "Enabled only in state-conditioned fit variants.",
    ),
    ModelRegistryEntry(
        "observation",
        "ev_size_observation",
        "condition adapter",
        "shared linear and quadratic pulse-dose response correction",
        "implemented_diagnostic",
        "observation_calibration_optional",
        "",
        "Empirical diagnostic adapter, not a validated constitutive mechanism.",
    ),
    ModelRegistryEntry(
        "observation",
        "ev_size_observation",
        "instrument observation",
        "instrument broadening and common-bin rebinning",
        "implemented_preliminary",
        "observation_calibration",
        "",
        "Parameters currently live in the size-resolved calibration pipeline.",
    ),
    ModelRegistryEntry(
        "bridge",
        "experimental_bridge",
        "experimental-to-model scaling",
        "volume, cell-count, concentration, and normalization transforms",
        "implemented_preliminary",
        "data_interface",
        "",
    ),
    ModelRegistryEntry(
        "cross-cutting",
        "cell_state",
        "cell-state modifiers",
        "membrane, calcium handling, basal EV release, and stress sensitivity modifiers",
        "implemented_preliminary",
        "scenario_input",
        "",
        "Values are scenario inputs, not entries in the authoritative defaults YAML.",
    ),
    ModelRegistryEntry(
        "6",
        "cargo_potency",
        "cargo composition",
        "protein, RNA, lipid, antigen, and direct-load enrichment",
        "implemented_preliminary",
        "downstream",
        "cargo_potency",
    ),
    ModelRegistryEntry(
        "6",
        "cargo_potency",
        "recipient potency",
        "subtype-weighted saturating potency proxy",
        "implemented_preliminary",
        "downstream",
        "cargo_potency",
    ),
    ModelRegistryEntry(
        "7",
        "injury_quality",
        "cell injury and fate",
        "damage, stress, apoptosis, necrosis, and viability",
        "implemented_preliminary",
        "downstream",
        "injury_quality",
    ),
    ModelRegistryEntry(
        "7",
        "injury_quality",
        "particle quality",
        "debris/aggregate mixture, marker score, purity gate",
        "implemented_preliminary",
        "downstream",
        "injury_quality",
    ),
    ModelRegistryEntry(
        "8",
        "manufacturing_qc",
        "recovery and yield",
        "isolation recovery, yield, purity, and batch scaling",
        "implemented_preliminary",
        "downstream",
        "manufacturing_qc",
    ),
    ModelRegistryEntry(
        "8",
        "manufacturing_qc",
        "manufacturing objective",
        "weighted potency/yield/purity/viability objective",
        "implemented_preliminary",
        "downstream",
        "manufacturing_qc",
    ),
    ModelRegistryEntry(
        "input",
        "scenario",
        "pulse and exposure",
        "pulse protocol, exposure geometry, medium, and thermal inputs",
        "implemented",
        "scenario_input",
        "",
        "Scenario inputs are validated separately and are not part of the 241 scalar runtime defaults.",
    ),
    ModelRegistryEntry(
        "input",
        "scenario",
        "extracellular medium",
        "culture volume, viability handling, and sampling events",
        "implemented",
        "scenario_input",
        "",
        "Scenario inputs are validated separately and are not part of the 241 scalar runtime defaults.",
    ),
    ModelRegistryEntry(
        "orchestration",
        "simulation",
        "numerical integration",
        "coupled solve_ivp integration",
        "implemented_fixed",
        "core",
        "",
        "SimulationConfig.numerical_method is currently not used to select a solver.",
    ),
)


def flatten_parameter_mapping(parameters: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Flatten a nested parameter mapping while retaining its top-level module.

    Non-mapping leaves, including strings and lists, are retained.  Numeric
    leaves therefore correspond one-to-one with scalar defaults in the YAML.
    """

    rows: list[dict[str, Any]] = []
    for module, value in parameters.items():
        _flatten_value(str(module), (), value, rows)
    return rows


def build_parameter_snapshot(
    defaults: Mapping[str, Any] | None = None,
    *,
    effective_parameters: Mapping[str, Any] | None = None,
    fit_parameters: pd.DataFrame | Iterable[Mapping[str, Any]] | None = None,
    fit_module: str | None = None,
    selected_fit_variant: str | None = None,
    metadata: Mapping[str, Mapping[str, Any]] | None = None,
) -> pd.DataFrame:
    """Build a comprehensive runtime-default and calibration snapshot.

    ``effective_parameters`` may be either a sparse override or a complete
    parameter mapping; it is deep-merged onto ``defaults``.  Fit rows may use
    ``parameter``/``name``, ``initial``/``initial_guess``,
    ``lower``/``lower_bound``, ``upper``/``upper_bound``, and
    ``fitted``/``fitted_value``/``final`` column names.  Fit-specific
    ``unit``/``units``, ``parameter_class``/``role``, and ``submodule`` values
    override inferred metadata and are carried into the snapshot.

    When ``fit_module`` is supplied, unqualified fit parameters are kept in
    that namespace.  This is useful for observation-layer parameters that are
    intentionally separate from the mechanistic defaults.  With no
    ``fit_module``, an unqualified fit name is joined to an authoritative
    default only when its leaf name is unique; otherwise it is placed in the
    ``fit_parameters`` namespace.

    ``selected_fit_variant`` filters a multi-variant fit table before joining;
    callers should use it when the same parameter appears in several candidate
    model variants.
    """

    defaults_mapping = deepcopy(dict(defaults or load_default_parameters()))
    effective_mapping = merge_parameters(
        defaults_mapping, dict(effective_parameters or {})
    )
    default_rows = flatten_parameter_mapping(defaults_mapping)
    effective_rows = flatten_parameter_mapping(effective_mapping)
    effective_lookup = {row["qualified_name"]: row["value"] for row in effective_rows}
    combined_metadata = deepcopy(_BUILTIN_PARAMETER_METADATA)
    for key, value in (metadata or {}).items():
        combined_metadata.setdefault(str(key), {}).update(dict(value))

    records: list[dict[str, Any]] = []
    for row in default_rows:
        qualified_name = row["qualified_name"]
        default_value = _normalise_scalar(row["value"])
        effective_value = _normalise_scalar(
            effective_lookup.get(qualified_name, default_value)
        )
        changed = not _values_equal(default_value, effective_value)
        record = _empty_parameter_record()
        record.update(
            {
                "module": row["module"],
                "layer": _module_layer(row["module"]),
                "submodule": _infer_submodule(row["module"], row["parameter"]),
                "parameter_path": row["parameter_path"],
                "qualified_name": qualified_name,
                "parameter": row["parameter"],
                "default_value": _table_value(default_value),
                "effective_value": _table_value(effective_value),
                "unit": infer_parameter_unit(row["parameter"]),
                "parameter_class": classify_parameter(row["parameter"]),
                "runtime_status": "active",
                "fit_status": "not_fitted",
                "status": "fixed_override" if changed else "fixed_default",
                "source": "runtime_default_yaml",
            }
        )
        record.update(combined_metadata.get(qualified_name, {}))
        records.append(record)

    default_qualified_names = {str(row["qualified_name"]) for row in default_rows}
    for row in effective_rows:
        qualified_name = str(row["qualified_name"])
        if qualified_name in default_qualified_names:
            continue
        effective_value = _normalise_scalar(row["value"])
        record = _empty_parameter_record()
        record.update(
            {
                "layer": _module_layer(row["module"]),
                "module": row["module"],
                "submodule": _infer_submodule(row["module"], row["parameter"]),
                "parameter_path": row["parameter_path"],
                "qualified_name": qualified_name,
                "parameter": row["parameter"],
                "effective_value": _table_value(effective_value),
                "unit": infer_parameter_unit(row["parameter"]),
                "parameter_class": classify_parameter(row["parameter"]),
                "runtime_status": "unknown_override",
                "fit_status": "not_fitted",
                "status": "fixed_override",
                "source": "run_override",
                "notes": (
                    "Override key is absent from the authoritative runtime defaults; "
                    "the current permissive merge may leave it unused."
                ),
            }
        )
        records.append(record)

    fit_rows = _normalise_fit_rows(fit_parameters)
    if selected_fit_variant is not None:
        fit_rows = [row for row in fit_rows if row["variant"] == selected_fit_variant]
    if fit_rows:
        _join_fit_rows(records, fit_rows, fit_module=fit_module)

    frame = pd.DataFrame(records, columns=PARAMETER_SNAPSHOT_COLUMNS)
    return frame.sort_values(
        ["module", "parameter_path"], kind="stable", ignore_index=True
    )


def build_model_registry(
    fit_engagement: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Return the framework model registry with optional run-specific fit roles.

    Engagement keys may be a module (for example ``ev_release``) or a qualified
    ``module.submodule`` key.  The mapped value should be a concise role such as
    ``fixed upstream trajectory`` or ``fitted observation layer``.
    """

    engagement = dict(fit_engagement or {})
    rows: list[dict[str, Any]] = []
    for entry in _MODEL_REGISTRY:
        qualified = f"{entry.module}.{entry.submodule}"
        fit_role = engagement.get(
            qualified, engagement.get(entry.module, "not_engaged")
        )
        rows.append(
            {
                "layer": entry.layer,
                "module": entry.module,
                "submodule": entry.submodule,
                "model": entry.model,
                "implementation_status": entry.implementation_status,
                "runtime_role": entry.runtime_role,
                "default_section": entry.default_section,
                "fit_engaged": fit_role != "not_engaged",
                "fit_role": fit_role,
                "notes": entry.notes,
            }
        )
    return pd.DataFrame(rows)


def infer_parameter_unit(parameter: str) -> str:
    """Infer only units that are explicit in, or strongly implied by, a name."""

    name = str(parameter)
    lower = name.lower()
    explicit_suffixes = (
        ("particles_per_model_unit", "particles model-unit^-1"),
        ("particles_per_ml", "particles mL^-1"),
        ("j_m3_k", "J m^-3 K^-1"),
        ("j_m3", "J m^-3"),
        ("kv_cm", "kV cm^-1"),
        ("v_m", "V m^-1"),
        ("s_m", "S m^-1"),
        ("um_s", "uM s^-1"),
        ("mm_s", "mM s^-1"),
        ("per_ml", "mL^-1"),
        ("_um", "uM"),
        ("_mm", "mM"),
        ("_nm", "nm"),
        ("_m2", "m^-2"),
        ("_v", "V"),
        ("_c", "degC"),
        ("_hz", "Hz"),
        ("_ns", "ns"),
        ("_h", "h"),
        ("_m", "m"),
    )
    for suffix, unit in explicit_suffixes:
        if lower.endswith(suffix):
            return unit
    if lower.startswith("tau_") and lower.endswith("_s"):
        return "s"
    if lower.endswith("_time_s") or lower in {"t_start_s", "t_end_s", "output_dt_s"}:
        return "s"
    if (
        lower.startswith("k_")
        or lower.endswith("_rate_s")
        or lower.endswith("_factor_s")
    ):
        return "s^-1"
    if lower in {
        "baseline_sev_rate",
        "baseline_mlev_rate",
        "baseline_ab_rate",
        "damage_rate",
        "repair_rate",
        "shedding_rate_scale",
    }:
        return "model-unit s^-1"
    if lower == "cell_count" or lower.endswith("_number"):
        return "count"
    if _is_dimensionless_name(lower):
        return "dimensionless"
    return "unspecified/model units"


def classify_parameter(parameter: str) -> str:
    """Assign a conservative constitutive/kinetic/coupling/other class."""

    name = str(parameter)
    lower = name.lower()
    if (
        name.startswith("k_")
        or lower.startswith("tau_")
        or lower.startswith("j_")
        or "vmax" in lower
        or lower.endswith("_rate_s")
        or lower.endswith("_factor_s")
        or lower
        in {
            "baseline_sev_rate",
            "baseline_mlev_rate",
            "baseline_ab_rate",
            "damage_rate",
            "repair_rate",
            "shedding_rate_scale",
        }
    ):
        return "kinetic"
    if any(
        token in lower
        for token in (
            "_weight",
            "_coupling",
            "_modifier",
            "_gain",
            "_penalty",
            "_multiplier",
            "_coefficient",
        )
    ):
        return "coupling"
    if (
        name.startswith("K_")
        or "baseline" in lower
        or "initial" in lower
        or "threshold" in lower
        or "reference" in lower
        or "capacity" in lower
        or lower.endswith("_max")
        or "half_max" in lower
        or "half_voltage" in lower
        or lower in {"cell_count", "harvest_time_h"}
    ):
        return "constitutive"
    return "other"


def _module_layer(module: str) -> str:
    return {
        "pulse": "1",
        "dosimetry": "1",
        "electrodynamics": "2",
        "ion_transport": "3",
        "remodeling_repair": "4",
        "ev_release": "5",
        "extracellular_kinetics": "5b",
        "ev_size_observation": "observation",
        "cargo_potency": "6",
        "injury_quality": "7",
        "manufacturing_qc": "8",
    }.get(str(module), "cross-cutting")


def _infer_submodule(module: str, parameter: str) -> str:
    """Provide a broad mechanistic grouping without inventing new defaults."""

    lower = str(parameter).lower()
    if module == "pulse":
        return "pulse descriptors"
    if module == "dosimetry":
        return "field geometry"
    if module == "electrodynamics":
        if "pore" in lower or "permeability" in lower:
            return "electroporation"
        if "delta_v" in lower:
            return "organelle polarization"
        return "membrane charging"
    if module == "ion_transport":
        if any(token in lower for token in ("mito", "mptp")):
            return "mitochondrial calcium"
        if any(token in lower for token in ("na_", "k_", "cl_", "osmotic")):
            return "monovalent ions and osmosis"
        if "ros" in lower or "atp" in lower:
            return "metabolic stress"
        return "calcium influx and stores"
    if module == "remodeling_repair":
        if any(
            token in lower
            for token in ("annex", "lysosomal", "reseal", "repair", "shedding")
        ):
            return "membrane repair"
        if any(token in lower for token in ("calpain", "actin", "actomyo")):
            return "cytoskeletal remodeling"
        return "calcium microdomains and PS"
    if module == "ev_release":
        if "apopt" in lower or lower.endswith("ab_rate"):
            return "apoptotic-body pathway"
        if any(
            token in lower
            for token in ("mlev", "budding", "scission", "arf6", "rhoa", "tension")
        ):
            return "plasma-membrane mlEV pathway"
        return "endosomal sEV pathway"
    if module == "extracellular_kinetics":
        if "aggregation" in lower:
            return "aggregation"
        if any(
            token in lower for token in ("loss", "uptake", "degradation", "adsorption")
        ):
            return "particle loss"
        if any(
            token in lower for token in ("assay", "dilution", "recovery", "background")
        ):
            return "sampling and assay"
        return "extracellular stocks"
    if module == "ev_size_observation":
        if "state_shift" in lower:
            return "state adapter"
        if "dose_response" in lower:
            return "condition adapter"
        if "loss" in lower or "half_life" in lower:
            return "size-dependent loss"
        if "instrument" in lower or "assay" in lower:
            return "instrument observation"
        if any(token in lower for token in ("median", "geometric_sd", "source_scale")):
            return "pathway-to-size kernel"
        return "observation adapter"
    if module == "cargo_potency":
        if "recipient" in lower or "potency" in lower or "subtype" in lower:
            return "recipient potency"
        return "cargo composition"
    if module == "injury_quality":
        if any(
            token in lower
            for token in ("debris", "aggregate", "contamination", "marker", "quality")
        ):
            return "particle quality"
        return "cell injury and fate"
    if module == "manufacturing_qc":
        if any(
            token in lower
            for token in ("weight", "objective", "potency", "batch", "scalability")
        ):
            return "manufacturing objective"
        return "recovery and yield"
    return "unclassified"


def _flatten_value(
    module: str,
    path: tuple[str, ...],
    value: Any,
    rows: list[dict[str, Any]],
) -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            _flatten_value(module, (*path, str(key)), child, rows)
        return
    parameter_path = ".".join(path)
    rows.append(
        {
            "module": module,
            "parameter_path": parameter_path,
            "qualified_name": f"{module}.{parameter_path}",
            "parameter": path[-1] if path else module,
            "value": value,
        }
    )


def _empty_parameter_record() -> dict[str, Any]:
    return {column: pd.NA for column in PARAMETER_SNAPSHOT_COLUMNS}


def _normalise_fit_rows(
    fit_parameters: pd.DataFrame | Iterable[Mapping[str, Any]] | None,
) -> list[dict[str, Any]]:
    if fit_parameters is None:
        return []
    if isinstance(fit_parameters, pd.DataFrame):
        raw_rows = fit_parameters.to_dict(orient="records")
    else:
        raw_rows = [dict(row) for row in fit_parameters]
    normalised: list[dict[str, Any]] = []
    for raw in raw_rows:
        parameter = _first_present(
            raw, "qualified_name", "parameter_path", "parameter", "name"
        )
        if parameter is None or str(parameter).strip() == "":
            raise ValueError("Each fit-parameter row must provide a parameter name")
        normalised.append(
            {
                "parameter": str(parameter),
                "module": _first_present(raw, "module"),
                "submodule": _first_present(raw, "submodule"),
                "unit": _first_present(raw, "unit", "units"),
                "parameter_class": _first_present(raw, "parameter_class", "role"),
                "initial": _first_present(raw, "initial", "initial_guess"),
                "lower": _first_present(raw, "lower", "lower_bound"),
                "upper": _first_present(raw, "upper", "upper_bound"),
                "fitted": _first_present(
                    raw, "fitted", "fitted_value", "final", "final_value"
                ),
                "notes": _first_present(raw, "notes"),
                "variant": _first_present(raw, "variant"),
            }
        )
    return normalised


def _join_fit_rows(
    records: list[dict[str, Any]],
    fit_rows: list[dict[str, Any]],
    *,
    fit_module: str | None,
) -> None:
    by_qualified = {str(row["qualified_name"]): row for row in records}
    leaf_matches: dict[str, list[dict[str, Any]]] = {}
    for row in records:
        leaf_matches.setdefault(str(row["parameter"]), []).append(row)

    for fit in fit_rows:
        supplied = str(fit["parameter"])
        explicit_module = None if _is_missing(fit["module"]) else str(fit["module"])
        target: dict[str, Any] | None = None
        if supplied in by_qualified:
            target = by_qualified[supplied]
        elif explicit_module and f"{explicit_module}.{supplied}" in by_qualified:
            target = by_qualified[f"{explicit_module}.{supplied}"]
        elif fit_module and f"{fit_module}.{supplied}" in by_qualified:
            target = by_qualified[f"{fit_module}.{supplied}"]
        elif (
            fit_module is None
            and not explicit_module
            and len(leaf_matches.get(supplied, [])) == 1
        ):
            target = leaf_matches[supplied][0]

        if target is None:
            module = explicit_module or fit_module or "fit_parameters"
            parameter_path = supplied
            if supplied.startswith(f"{module}."):
                parameter_path = supplied[len(module) + 1 :]
            elif "." in supplied and explicit_module is None and fit_module is None:
                module, parameter_path = supplied.split(".", 1)
            target = _empty_parameter_record()
            target.update(
                {
                    "module": module,
                    "layer": _module_layer(module),
                    "submodule": _infer_submodule(
                        module, parameter_path.rsplit(".", 1)[-1]
                    ),
                    "parameter_path": parameter_path,
                    "qualified_name": f"{module}.{parameter_path}",
                    "parameter": parameter_path.rsplit(".", 1)[-1],
                    "unit": infer_parameter_unit(parameter_path.rsplit(".", 1)[-1]),
                    "parameter_class": classify_parameter(
                        parameter_path.rsplit(".", 1)[-1]
                    ),
                    "runtime_status": "observation_calibration",
                    "status": "configured_for_fit",
                    "source": "fit_registry",
                    "notes": "Calibration parameter outside the authoritative runtime defaults YAML.",
                }
            )
            records.append(target)
            by_qualified[str(target["qualified_name"])] = target
            leaf_matches.setdefault(str(target["parameter"]), []).append(target)

        target["fit_initial"] = fit["initial"]
        target["fit_lower_bound"] = fit["lower"]
        target["fit_upper_bound"] = fit["upper"]
        target["fit_final"] = fit["fitted"]
        has_final = not _is_missing(fit["fitted"])
        target["fit_status"] = "fitted" if has_final else "configured_for_fit"
        target["fit_variant"] = fit["variant"]
        target["status"] = target["fit_status"]
        if not _is_missing(fit["unit"]):
            target["unit"] = fit["unit"]
        if not _is_missing(fit["parameter_class"]):
            target["parameter_class"] = fit["parameter_class"]
        if not _is_missing(fit["submodule"]):
            target["submodule"] = fit["submodule"]
        if fit["variant"] is not None and not _is_missing(fit["variant"]):
            variant_note = f"Fit variant: {fit['variant']}."
            target["notes"] = _append_note(target.get("notes"), variant_note)
        if fit["notes"] is not None and not _is_missing(fit["notes"]):
            target["notes"] = _append_note(target.get("notes"), str(fit["notes"]))


def _first_present(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping and not _is_missing(mapping[key]):
            return mapping[key]
    return None


def _is_missing(value: Any) -> bool:
    if value is None or value is pd.NA:
        return True
    try:
        result = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return bool(result) if isinstance(result, (bool, type(pd.NA))) else False


def _values_equal(left: Any, right: Any) -> bool:
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=0.0)
    return left == right


def _table_value(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True)
    return value


def _normalise_scalar(value: Any) -> Any:
    """Represent YAML numeric-looking strings numerically in audit tables.

    PyYAML follows YAML 1.1 number resolution, where exponent forms such as
    ``1.0e12`` may be loaded as strings unless the exponent has an explicit
    sign.  The runtime coercion layer accepts those values; the snapshot should
    nevertheless expose them as numbers for sorting and downstream analysis.
    """

    if isinstance(value, str) and re.fullmatch(
        r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)[eE][+-]?\d+", value.strip()
    ):
        return float(value)
    return value


def _append_note(current: Any, addition: str) -> str:
    if current is None or current is pd.NA or str(current).strip() == "":
        return addition
    return f"{str(current).rstrip()} {addition}"


def _is_dimensionless_name(lower: str) -> bool:
    return bool(
        re.search(
            r"(?:fraction|factor|weight|scale|slope|coefficient|modifier|gain|"
            r"penalty|multiplier|efficiency|activity|consistency|variability)$",
            lower,
        )
        or lower.startswith("n_")
        or "hill" in lower
        or lower in {"atp_baseline", "ros_baseline", "ps_max", "geometric_sd"}
    )
