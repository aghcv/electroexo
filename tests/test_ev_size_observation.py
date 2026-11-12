from __future__ import annotations

import math

import numpy as np
import pytest

from electro_exocytosis.models.ev_size_observation import (
    SIZE_PATHWAYS,
    PathwaySizeKernelParams,
    SizeResolvedKineticsParams,
    build_size_observation_matrix,
    lognormal_bin_probabilities,
    simulate_size_resolved_extracellular_kinetics,
)


def _pathway_params(
    *, state_shift_coefficient: float = 0.0
) -> dict[str, PathwaySizeKernelParams]:
    return {
        "sEV": PathwaySizeKernelParams(
            90.0,
            1.45,
            source_scale_particles_per_model_unit=2.0,
            state_shift_coefficient=state_shift_coefficient,
        ),
        "mlEV": PathwaySizeKernelParams(180.0, 1.50),
        "AB": PathwaySizeKernelParams(
            350.0, 1.60, source_scale_particles_per_model_unit=0.5
        ),
    }


def _zero_release(time_s: np.ndarray) -> dict[str, np.ndarray]:
    return {pathway: np.zeros_like(time_s) for pathway in SIZE_PATHWAYS}


def test_lognormal_kernel_is_nonnegative_normalized_and_integrated() -> None:
    edges = np.array([20.0, 50.0, 100.0, 200.0, 500.0])
    probabilities = lognormal_bin_probabilities(edges, 100.0, 1.5)

    assert probabilities.shape == (4,)
    assert np.all(probabilities >= 0.0)
    assert probabilities.sum() == pytest.approx(1.0)
    assert probabilities[1] == pytest.approx(probabilities[2])

    with pytest.raises(ValueError, match="greater than one"):
        lognormal_bin_probabilities(edges, 100.0, 1.0)
    with pytest.raises(ValueError, match="strictly positive"):
        lognormal_bin_probabilities([0.0, 100.0], 50.0, 1.5)


def test_zero_width_observation_matrix_is_identity_and_conserves_mass() -> None:
    edges = np.array([50.0, 100.0, 200.0, 400.0])
    identity = build_size_observation_matrix(edges, edges, 0.0)
    assert np.allclose(identity, np.eye(3))

    coarse = build_size_observation_matrix(edges, [50.0, 200.0, 400.0], 0.0)
    latent_particles = np.array([2.0, 3.0, 5.0])
    observed_particles = coarse @ latent_particles
    assert np.allclose(coarse.sum(axis=0), 1.0)
    assert observed_particles.sum() == pytest.approx(latent_particles.sum())


def test_no_loss_conserves_population_scaled_pathway_source() -> None:
    time_s = np.array([0.0, 1.0, 2.0])
    cumulative = {
        "sEV": np.array([0.0, 2.0, 5.0]),
        "mlEV": np.array([0.0, 1.0, 2.0]),
        "AB": np.array([0.0, 0.0, 1.0]),
    }
    result = simulate_size_resolved_extracellular_kinetics(
        time_s,
        cumulative,
        [20.0, 50.0, 100.0, 200.0, 500.0, 1000.0],
        _pathway_params(),
        SizeResolvedKineticsParams(),
        initial_cell_density_per_ml=10.0,
    )

    expected_by_pathway = {"sEV": 100.0, "mlEV": 20.0, "AB": 5.0}
    for pathway, expected in expected_by_pathway.items():
        assert result.pathway_latent_concentration_particles_per_ml[pathway][
            -1
        ].sum() == pytest.approx(expected)
        assert np.allclose(result.kernel_probabilities[pathway].sum(axis=1), 1.0)
    assert result.total_latent_concentration_particles_per_ml[
        -1
    ].sum() == pytest.approx(sum(expected_by_pathway.values()))


def test_positive_loss_exponent_removes_larger_bins_faster() -> None:
    time_s = np.array([0.0, 1.0])
    result = simulate_size_resolved_extracellular_kinetics(
        time_s,
        _zero_release(time_s),
        [50.0, 100.0, 200.0],
        _pathway_params(),
        SizeResolvedKineticsParams(
            effective_loss_rate_s=math.log(2.0),
            loss_reference_diameter_nm=100.0,
            loss_size_exponent=1.0,
        ),
        initial_cell_density_per_ml=1.0,
        initial_concentration_particles_per_ml=[100.0, 100.0],
    )

    terminal = result.ambient_latent_concentration_particles_per_ml[-1]
    assert result.loss_rates_s[1] > result.loss_rates_s[0]
    assert terminal[1] < terminal[0]
    assert np.allclose(
        terminal,
        100.0 * np.exp(-result.loss_rates_s),
    )


def test_midpoint_state_shift_moves_kernel_toward_larger_bins() -> None:
    time_s = np.array([0.0, 1.0])
    cumulative = _zero_release(time_s)
    cumulative["sEV"] = np.array([0.0, 10.0])
    edges = [20.0, 50.0, 100.0, 200.0, 500.0]

    baseline = simulate_size_resolved_extracellular_kinetics(
        time_s,
        cumulative,
        edges,
        _pathway_params(state_shift_coefficient=0.0),
        initial_cell_density_per_ml=1.0,
        state_signals={"sEV": [0.0, 2.0]},
    )
    shifted = simulate_size_resolved_extracellular_kinetics(
        time_s,
        cumulative,
        edges,
        _pathway_params(state_shift_coefficient=math.log(2.0)),
        initial_cell_density_per_ml=1.0,
        state_signals={"sEV": [0.0, 2.0]},
    )

    baseline_kernel = baseline.kernel_probabilities["sEV"][0]
    shifted_kernel = shifted.kernel_probabilities["sEV"][0]
    centers = np.sqrt(np.asarray(edges[:-1]) * np.asarray(edges[1:]))
    assert shifted_kernel @ centers > baseline_kernel @ centers
    assert shifted_kernel[-1] > baseline_kernel[-1]


def test_assay_recovery_and_background_produce_finite_observations() -> None:
    time_s = np.array([0.0, 1.0])
    result = simulate_size_resolved_extracellular_kinetics(
        time_s,
        _zero_release(time_s),
        [50.0, 100.0, 200.0],
        _pathway_params(),
        SizeResolvedKineticsParams(assay_recovery_fraction=0.5),
        initial_cell_density_per_ml=1.0,
        initial_concentration_particles_per_ml=[10.0, 20.0],
        observed_background_concentration_particles_per_ml=[1.0, 2.0],
    )

    assert np.all(np.isfinite(result.observed_concentration_particles_per_ml))
    assert np.allclose(
        result.observed_concentration_particles_per_ml,
        [[6.0, 12.0], [6.0, 12.0]],
    )


def test_instrument_error_matrix_is_finite_and_nonnegative() -> None:
    matrix = build_size_observation_matrix(
        [25.0, 50.0, 100.0, 200.0, 400.0],
        [20.0, 40.0, 80.0, 160.0, 320.0, 640.0],
        instrument_log_diameter_sd=0.15,
    )
    assert matrix.shape == (5, 4)
    assert np.all(np.isfinite(matrix))
    assert np.all(matrix >= 0.0)
    assert np.all(matrix.sum(axis=0) <= 1.0 + 1.0e-12)
