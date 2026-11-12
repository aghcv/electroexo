from __future__ import annotations

import pandas as pd
import pytest

from electro_exocytosis.io.readers import load_default_parameters
from electro_exocytosis.parameter_registry import (
    build_model_registry,
    build_parameter_snapshot,
    classify_parameter,
    flatten_parameter_mapping,
    infer_parameter_unit,
)


def test_flatten_parameter_mapping_preserves_nested_paths() -> None:
    rows = flatten_parameter_mapping(load_default_parameters())
    lookup = {row["qualified_name"]: row for row in rows}

    assert len(rows) == 241
    assert lookup["cargo_potency.potency_weights.RNA"]["value"] == pytest.approx(0.25)
    assert (
        lookup["cargo_potency.potency_weights.RNA"]["parameter_path"]
        == "potency_weights.RNA"
    )


def test_parameter_snapshot_joins_overrides_and_separate_fit_layer() -> None:
    fit = pd.DataFrame(
        [
            {
                "variant": "state_conditioned",
                "parameter": "effective_half_life_h",
                "initial": 2.0,
                "lower_bound": 0.1,
                "upper_bound": 48.0,
                "fitted": 1.25,
                "units": "hours",
                "role": "calibrated kinetic",
                "submodule": "effective disappearance",
            }
        ]
    )
    snapshot = build_parameter_snapshot(
        effective_parameters={
            "extracellular_kinetics": {"effective_loss_rate_s": 0.001}
        },
        fit_parameters=fit,
        fit_module="ev_size_observation",
    ).set_index("qualified_name")

    loss = snapshot.loc["extracellular_kinetics.effective_loss_rate_s"]
    assert loss["default_value"] == pytest.approx(0.0)
    assert loss["effective_value"] == pytest.approx(0.001)
    assert loss["status"] == "fixed_override"
    assert loss["parameter_class"] == "kinetic"

    half_life = snapshot.loc["ev_size_observation.effective_half_life_h"]
    assert pd.isna(half_life["default_value"])
    assert half_life["fit_initial"] == pytest.approx(2.0)
    assert half_life["fit_final"] == pytest.approx(1.25)
    assert half_life["fit_status"] == "fitted"
    assert half_life["fit_variant"] == "state_conditioned"
    assert half_life["runtime_status"] == "observation_calibration"
    assert half_life["unit"] == "hours"
    assert half_life["parameter_class"] == "calibrated kinetic"
    assert half_life["submodule"] == "effective disappearance"
    assert "state_conditioned" in half_life["notes"]


def test_parameter_snapshot_can_select_one_fit_variant() -> None:
    snapshot = build_parameter_snapshot(
        fit_parameters=[
            {"variant": "first", "parameter": "scale", "fitted": 1.0},
            {"variant": "second", "parameter": "scale", "fitted": 2.0},
        ],
        fit_module="observation_fit",
        selected_fit_variant="second",
    )
    fitted = snapshot[snapshot["fit_status"] == "fitted"]

    assert len(fitted) == 1
    assert fitted.iloc[0]["qualified_name"] == "observation_fit.scale"
    assert fitted.iloc[0]["fit_final"] == pytest.approx(2.0)
    assert fitted.iloc[0]["fit_variant"] == "second"


def test_size_observation_parameters_receive_specific_submodules() -> None:
    snapshot = build_parameter_snapshot(
        fit_parameters=[
            {"variant": "fit", "parameter": "sEV_state_shift", "fitted": 0.2},
            {
                "variant": "fit",
                "parameter": "dose_response_linear",
                "fitted": -0.1,
            },
            {
                "variant": "fit",
                "parameter": "sEV_median_diameter_nm",
                "fitted": 110.0,
            },
            {
                "variant": "fit",
                "parameter": "sEV_geometric_sd",
                "fitted": 1.4,
            },
            {
                "variant": "fit",
                "parameter": "sEV_source_scale_particles_per_model_unit",
                "fitted": 1.0e4,
            },
        ],
        fit_module="ev_size_observation",
        selected_fit_variant="fit",
    ).set_index("parameter")

    assert snapshot.loc["sEV_state_shift", "submodule"] == "state adapter"
    assert snapshot.loc["dose_response_linear", "submodule"] == "condition adapter"
    assert (
        snapshot.loc["sEV_median_diameter_nm", "submodule"] == "pathway-to-size kernel"
    )
    assert snapshot.loc["sEV_geometric_sd", "submodule"] == "pathway-to-size kernel"
    assert (
        snapshot.loc["sEV_source_scale_particles_per_model_unit", "submodule"]
        == "pathway-to-size kernel"
    )


def test_parameter_snapshot_marks_known_declared_but_unused_defaults() -> None:
    snapshot = build_parameter_snapshot().set_index("qualified_name")

    assert (
        snapshot.loc["pulse.geometry_reference", "runtime_status"]
        == "declared_not_consumed"
    )
    assert (
        snapshot.loc["dosimetry.dish_factor", "runtime_status"]
        == "declared_not_consumed"
    )
    assert (
        snapshot.loc["remodeling_repair.tau_repair_s", "runtime_status"]
        == "declared_not_consumed"
    )
    assert snapshot.loc["ion_transport.Ca_baseline_uM", "runtime_status"] == "active"
    assert snapshot.loc[
        "electrodynamics.pore_density_max_m2", "default_value"
    ] == pytest.approx(1.0e12)
    assert snapshot.loc[
        "manufacturing_qc.cell_count", "default_value"
    ] == pytest.approx(1.0e6)


def test_parameter_snapshot_retains_unknown_sparse_override_leaves() -> None:
    snapshot = build_parameter_snapshot(
        effective_parameters={
            "extracellular_kinetics": {"future_sink_rate_s": 0.004},
            "future_module": {"new_coupling_weight": 0.3},
        }
    ).set_index("qualified_name")

    future_sink = snapshot.loc["extracellular_kinetics.future_sink_rate_s"]
    assert pd.isna(future_sink["default_value"])
    assert future_sink["effective_value"] == pytest.approx(0.004)
    assert future_sink["source"] == "run_override"
    assert future_sink["runtime_status"] == "unknown_override"
    assert "permissive merge" in future_sink["notes"]

    future_module = snapshot.loc["future_module.new_coupling_weight"]
    assert future_module["effective_value"] == pytest.approx(0.3)
    assert future_module["source"] == "run_override"


def test_unqualified_fit_name_matches_unique_authoritative_parameter() -> None:
    snapshot = build_parameter_snapshot(
        fit_parameters=[
            {
                "name": "effective_loss_rate_s",
                "initial_guess": 0.0,
                "lower": 0.0,
                "upper": 0.1,
                "final": 0.002,
            }
        ]
    ).set_index("qualified_name")

    fitted = snapshot.loc["extracellular_kinetics.effective_loss_rate_s"]
    assert fitted["fit_final"] == pytest.approx(0.002)
    assert fitted["source"] == "runtime_default_yaml"


def test_model_registry_supports_module_and_submodule_fit_roles() -> None:
    registry = build_model_registry(
        {
            "ev_release": "fixed upstream trajectory",
            "ev_size_observation.pathway-to-size kernel": "fitted observation layer",
        }
    )

    ev_release = registry[registry["module"] == "ev_release"]
    assert ev_release["fit_engaged"].all()
    assert set(ev_release["fit_role"]) == {"fixed upstream trajectory"}
    size_kernel = registry[
        (registry["module"] == "ev_size_observation")
        & (registry["submodule"] == "pathway-to-size kernel")
    ].iloc[0]
    assert bool(size_kernel["fit_engaged"])
    assert size_kernel["fit_role"] == "fitted observation layer"
    size_submodules = set(
        registry.loc[registry["module"] == "ev_size_observation", "submodule"]
    )
    assert {
        "source conversion",
        "pathway-to-size kernel",
        "size-dependent loss",
        "state adapter",
        "condition adapter",
        "instrument observation",
    }.issubset(size_submodules)
    solver = registry[registry["module"] == "simulation"].iloc[0]
    assert "numerical_method" in solver["notes"]


def test_parameter_unit_and_class_are_conservative() -> None:
    assert infer_parameter_unit("tau_Ca_homeostasis_s") == "s"
    assert infer_parameter_unit("effective_loss_rate_s") == "s^-1"
    assert infer_parameter_unit("Ca_baseline_uM") == "uM"
    assert infer_parameter_unit("mystery") == "unspecified/model units"
    assert classify_parameter("k_MVB_docking_s") == "kinetic"
    assert classify_parameter("K_fusion_Ca_uM") == "constitutive"
    assert classify_parameter("rab27_weight") == "coupling"
