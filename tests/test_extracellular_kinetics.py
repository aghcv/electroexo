from __future__ import annotations

import math

import numpy as np
import pytest

from electro_exocytosis.models.extracellular_kinetics import (
    ExtracellularKineticsParams,
    MediumSamplingEvent,
    coerce_extracellular_kinetics_params,
    simulate_extracellular_kinetics,
)


def _zero_release(time_s: np.ndarray) -> dict[str, np.ndarray]:
    zeros = np.zeros_like(time_s, dtype=float)
    return {"sEV": zeros, "mlEV": zeros, "AB": zeros}


def test_zero_loss_matches_population_scaled_cumulative_release() -> None:
    time_s = np.array([0.0, 1800.0, 3600.0])
    cumulative = {
        "sEV": np.array([0.0, 2.0, 5.0]),
        "mlEV": np.array([0.0, 1.0, 1.0]),
        "AB": np.zeros(3),
    }

    result = simulate_extracellular_kinetics(
        time_s,
        cumulative,
        ExtracellularKineticsParams(),
        initial_cell_density_per_ml=10.0,
        initial_volume_ml=2.0,
    )

    assert np.allclose(
        result.timeseries["sEV_extracellular_concentration_particles_per_ml"],
        10.0 * cumulative["sEV"],
    )
    assert np.allclose(
        result.timeseries["total_extracellular_concentration_particles_per_ml"],
        10.0 * (cumulative["sEV"] + cumulative["mlEV"]),
    )


def test_first_order_loss_reproduces_known_half_life() -> None:
    time_s = np.array([0.0, 3600.0])
    params = ExtracellularKineticsParams(
        initial_sEV_concentration_particles_per_ml=100.0,
        effective_loss_rate_s=math.log(2.0) / 3600.0,
    )

    result = simulate_extracellular_kinetics(
        time_s,
        _zero_release(time_s),
        params,
        initial_cell_density_per_ml=1.0,
        initial_volume_ml=1.0,
    )

    concentration = result.timeseries[
        "sEV_extracellular_concentration_particles_per_ml"
    ]
    assert concentration[-1] == pytest.approx(50.0)


def test_transient_source_followed_by_loss_can_produce_a_peak_and_decline() -> None:
    time_s = np.array([0.0, 3600.0, 7200.0])
    cumulative = {
        "sEV": np.array([0.0, 100.0, 100.0]),
        "mlEV": np.zeros(3),
        "AB": np.zeros(3),
    }
    params = ExtracellularKineticsParams(
        degradation_rate_s=math.log(2.0) / 1800.0,
    )

    result = simulate_extracellular_kinetics(
        time_s,
        cumulative,
        params,
        initial_cell_density_per_ml=1.0,
        initial_volume_ml=1.0,
    )

    concentration = result.timeseries[
        "total_extracellular_concentration_particles_per_ml"
    ]
    assert concentration[1] > concentration[0]
    assert concentration[2] == pytest.approx(0.25 * concentration[1])


def test_sample_then_replace_dilutes_subsequent_concentration() -> None:
    time_s = np.array([0.0, 3600.0, 7200.0])
    params = ExtracellularKineticsParams(
        initial_sEV_concentration_particles_per_ml=100.0,
    )

    result = simulate_extracellular_kinetics(
        time_s,
        _zero_release(time_s),
        params,
        initial_cell_density_per_ml=1.0,
        initial_volume_ml=2.0,
        sampling_events=[
            MediumSamplingEvent(
                time_s=3600.0,
                sampled_volume_ml=1.0,
                replacement_volume_ml=1.0,
            )
        ],
    )

    concentration = result.timeseries[
        "sEV_extracellular_concentration_particles_per_ml"
    ]
    assert np.allclose(concentration, [100.0, 100.0, 50.0])
    assert result.event_log[0]["particles_removed"] == pytest.approx(100.0)
    assert result.event_log[0]["volume_after_ml"] == pytest.approx(2.0)


def test_reduced_aggregation_moves_small_particles_and_reduces_total_count() -> None:
    time_s = np.array([0.0, 3600.0])
    params = ExtracellularKineticsParams(
        initial_sEV_concentration_particles_per_ml=100.0,
        sEV_to_mlEV_aggregation_rate_s=math.log(2.0) / 3600.0,
        aggregation_particle_yield=0.5,
    )

    result = simulate_extracellular_kinetics(
        time_s,
        _zero_release(time_s),
        params,
        initial_cell_density_per_ml=1.0,
        initial_volume_ml=1.0,
    )

    assert result.timeseries[
        "sEV_extracellular_concentration_particles_per_ml"
    ][-1] == pytest.approx(50.0)
    assert result.timeseries[
        "mlEV_extracellular_concentration_particles_per_ml"
    ][-1] == pytest.approx(25.0)
    assert result.timeseries[
        "total_extracellular_concentration_particles_per_ml"
    ][-1] == pytest.approx(75.0)


def test_assay_transform_is_separate_from_true_extracellular_stock() -> None:
    time_s = np.array([0.0, 1.0])
    params = ExtracellularKineticsParams(
        initial_sEV_concentration_particles_per_ml=100.0,
        assay_recovery_fraction=0.5,
        assay_dilution_factor=2.0,
        assay_background_concentration_particles_per_ml=10.0,
    )

    result = simulate_extracellular_kinetics(
        time_s,
        _zero_release(time_s),
        params,
        initial_cell_density_per_ml=1.0,
        initial_volume_ml=1.0,
    )

    assert np.allclose(
        result.timeseries["measured_particle_concentration_particles_per_ml"],
        35.0,
    )


def test_nested_parameters_are_coerced_and_validated() -> None:
    params = coerce_extracellular_kinetics_params(
        {
            "extracellular_kinetics": {
                "effective_loss_rate_s": "0.001",
                "uptake_rate_s": "0.002",
                "assay_recovery_fraction": "0.7",
            }
        }
    )
    assert params.effective_loss_rate_s == pytest.approx(0.001)
    assert params.uptake_rate_s == pytest.approx(0.002)
    assert params.assay_recovery_fraction == pytest.approx(0.7)

    with pytest.raises(ValueError, match="cannot exceed one"):
        coerce_extracellular_kinetics_params(
            {"extracellular_kinetics": {"assay_recovery_fraction": 1.2}}
        )


def test_sampling_cannot_remove_more_than_the_available_volume() -> None:
    time_s = np.array([0.0, 10.0])
    with pytest.raises(ValueError, match="removes"):
        simulate_extracellular_kinetics(
            time_s,
            _zero_release(time_s),
            ExtracellularKineticsParams(),
            initial_cell_density_per_ml=1.0,
            initial_volume_ml=1.0,
            sampling_events=[
                MediumSamplingEvent(time_s=0.0, sampled_volume_ml=2.0)
            ],
        )
