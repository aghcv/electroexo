from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electro_exocytosis.experimental_bridge import ExperimentalObservationBridge
from electro_exocytosis.models.ev_size_observation import (
    PathwaySizeKernelParams,
    SizeResolvedKineticsParams,
    simulate_size_resolved_extracellular_kinetics,
)
from tools.fit_ffrci_size_resolved import (
    COMMON_SIZE_BIN_EDGES_NM,
    DEFAULT_BRIDGE_CONFIG,
    DEFAULT_DATA_DIR,
    EXOID_FILENAME,
    EXPOSURES,
    FIT_VARIANTS,
    OBSERVATION_TIMES_H,
    SizeResolvedFFRCIPredictor,
    build_upstream_trajectories,
    hellinger_distance,
    load_size_resolved_observations,
    plot_size_time_error_contours,
    plot_size_time_surface_overlay,
    prediction_metrics,
    signed_normalized_error_percent,
)


PRIVATE_EXOID_PATH = Path(DEFAULT_DATA_DIR) / EXOID_FILENAME
EXPECTED_CONDITION_TIMES = {
    ("1p20kV", 0.5),
    ("1p20kV", 1.0),
    ("1p20kV", 3.0),
    ("3p40kV", 0.5),
    ("3p40kV", 1.0),
    ("3p40kV", 3.0),
    ("5p40kV", 0.5),
    ("5p40kV", 1.0),
    ("5p40kV", 3.0),
}


def _require_private_exoid_data() -> None:
    if not PRIVATE_EXOID_PATH.exists():
        pytest.skip("gitignored FFRCI Exoid source CSV is not available")


@pytest.fixture(scope="module")
def real_observations():
    _require_private_exoid_data()
    return load_size_resolved_observations(Path(DEFAULT_DATA_DIR))


@pytest.fixture(scope="module")
def provisional_bridge() -> ExperimentalObservationBridge:
    return ExperimentalObservationBridge.from_yaml(DEFAULT_BRIDGE_CONFIG)


@pytest.fixture(scope="module")
def cached_upstream_trajectories(real_observations, provisional_bridge):
    return build_upstream_trajectories(
        provisional_bridge,
        pulse_width_ns=60.0,
        repetition_rate_hz=1.0,
    )


def test_loader_retains_all_available_histograms_and_common_bins(
    real_observations,
) -> None:
    observations = real_observations
    distribution_keys = observations[
        ["condition", "time_h", "replicate"]
    ].drop_duplicates()
    condition_time_keys = observations[["condition", "time_h"]].drop_duplicates()

    assert len(observations) == 390
    assert len(distribution_keys) == 26
    assert len(condition_time_keys) == 9
    assert set(condition_time_keys.itertuples(index=False, name=None)) == (
        EXPECTED_CONDITION_TIMES
    )

    missing_cell = distribution_keys[
        (distribution_keys["condition"] == "3p40kV")
        & np.isclose(distribution_keys["time_h"], 1.0)
    ]
    assert set(missing_cell["replicate"]) == {1, 2}
    assert not (
        (distribution_keys["condition"] == "3p40kV")
        & np.isclose(distribution_keys["time_h"], 1.0)
        & (distribution_keys["replicate"] == 3)
    ).any()

    expected_lower = np.arange(80.0, 380.0, 20.0)
    expected_upper = np.arange(100.0, 400.0, 20.0)
    expected_centers = np.arange(90.0, 380.0, 20.0)
    assert np.array_equal(COMMON_SIZE_BIN_EDGES_NM, np.arange(80.0, 400.0, 20.0))
    assert np.array_equal(
        np.sort(observations["size_bin_lower_nm"].unique()), expected_lower
    )
    assert np.array_equal(
        np.sort(observations["size_bin_upper_nm"].unique()), expected_upper
    )
    assert np.array_equal(
        np.sort(observations["size_bin_center_nm"].unique()), expected_centers
    )
    assert (
        observations.groupby(["condition", "time_h", "replicate"]).size()
        == len(expected_centers)
    ).all()


def test_hellinger_distance_handles_proportional_vectors_and_zero_bins() -> None:
    assert hellinger_distance(
        np.array([1.0, 2.0, 3.0]), np.array([7.0, 14.0, 21.0])
    ) == pytest.approx(0.0, abs=1.0e-15)

    distance = hellinger_distance(
        np.array([0.0, 2.0, 0.0, 3.0]),
        np.array([1.0, 2.0, 4.0, 0.0]),
    )
    assert np.isfinite(distance)
    assert 0.0 <= distance <= 1.0


def test_signed_normalized_error_percent_sign_and_zero_edges() -> None:
    observed = np.array([2.0, 1.0, 2.0, 0.0, 5.0, 0.0])
    predicted = np.array([2.0, 2.0, 1.0, 5.0, 0.0, 0.0])

    error = signed_normalized_error_percent(observed, predicted)
    expected = np.array([0.0, 100.0 / 3.0, -100.0 / 3.0, 100.0, -100.0])

    assert error[:5] == pytest.approx(expected)
    assert np.isnan(error[5])


def test_signed_normalized_error_percent_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="same shape"):
        signed_normalized_error_percent(np.array([1.0]), np.array([1.0, 2.0]))
    with pytest.raises(ValueError, match="nonnegative"):
        signed_normalized_error_percent(np.array([-1.0]), np.array([1.0]))
    with pytest.raises(ValueError, match="nonnegative"):
        signed_normalized_error_percent(np.array([1.0]), np.array([-1.0]))


def test_cached_upstream_trajectories_cover_all_conditions(
    cached_upstream_trajectories,
) -> None:
    assert set(cached_upstream_trajectories) == {
        "1p20kV",
        "3p40kV",
        "5p40kV",
    }
    for trajectory in cached_upstream_trajectories.values():
        assert np.all(np.diff(trajectory.time_s) > 0.0)
        assert np.all(np.isfinite(trajectory.viable_fraction))
        for cumulative in trajectory.cumulative_release.values():
            assert len(cumulative) == len(trajectory.time_s)
            assert np.all(np.isfinite(cumulative))
            assert np.all(np.diff(cumulative) >= -1.0e-9)


@pytest.mark.parametrize("variant_name", ["static_kernel", "state_conditioned"])
def test_initial_predictors_are_finite_nonnegative_and_scaled_once(
    variant_name,
    real_observations,
    provisional_bridge,
    cached_upstream_trajectories,
) -> None:
    variant = next(item for item in FIT_VARIANTS if item.name == variant_name)
    predictor = SizeResolvedFFRCIPredictor(
        real_observations,
        cached_upstream_trajectories,
        provisional_bridge,
        variant,
    )
    prediction = predictor.predict(predictor.initial_vector)
    metrics, total_frame, sample_frame = prediction_metrics(prediction)

    assert len(prediction) == 390
    assert len(total_frame) == 9
    assert len(sample_frame) == 26
    assert metrics["independent_total_targets"] == 9
    predicted_condition_times = set(
        total_frame[["condition", "time_h"]].itertuples(index=False, name=None)
    )
    assert predicted_condition_times == EXPECTED_CONDITION_TIMES

    concentration_columns = [
        "predicted_particles_per_ml",
        "predicted_true_particles_per_ml",
        "observed_particle_equivalents_per_initial_cell",
        "predicted_particle_equivalents_per_initial_cell",
    ]
    values = prediction[concentration_columns].to_numpy(dtype=float)
    assert np.all(np.isfinite(values))
    assert np.all(values >= 0.0)

    initial_cell_density_per_ml = (
        provisional_bridge.initial_cell_count / provisional_bridge.medium_volume_ml
    )
    assert np.allclose(
        prediction["observed_particle_equivalents_per_initial_cell"],
        prediction["treatment_particles_per_ml"] / initial_cell_density_per_ml,
    )
    assert np.allclose(
        prediction["predicted_particle_equivalents_per_initial_cell"],
        prediction["predicted_particles_per_ml"] / initial_cell_density_per_ml,
    )


def test_positive_loss_can_create_decline_from_nondecreasing_release() -> None:
    time_s = np.array([0.0, 1.0, 2.0])
    cumulative_release = {
        "sEV": np.array([0.0, 1.0, 1.0]),
        "mlEV": np.zeros(3),
        "AB": np.zeros(3),
    }
    assert np.all(np.diff(cumulative_release["sEV"]) >= 0.0)

    pathway_params = {
        pathway: PathwaySizeKernelParams(100.0, 1.5)
        for pathway in ("sEV", "mlEV", "AB")
    }
    result = simulate_size_resolved_extracellular_kinetics(
        time_s,
        cumulative_release,
        [50.0, 200.0],
        pathway_params,
        SizeResolvedKineticsParams(effective_loss_rate_s=math.log(2.0)),
        initial_cell_density_per_ml=1.0,
    )
    total = result.observed_concentration_particles_per_ml.sum(axis=1)

    assert total[1] > 0.0
    assert total[2] < total[1]


def _synthetic_size_time_prediction() -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    size_centers_nm = np.arange(90.0, 180.0, 20.0)
    for condition_index, exposure in enumerate(EXPOSURES):
        for time_index, time_h in enumerate(OBSERVATION_TIMES_H):
            for replicate in (1, 2, 3):
                is_missing_p3 = (
                    exposure.condition == "3p40kV" and time_h == 1.0 and replicate == 3
                )
                if is_missing_p3:
                    continue
                for size_index, size_nm in enumerate(size_centers_nm):
                    observed = (
                        1.0e8
                        * (1.0 + 0.20 * condition_index + 0.10 * time_index)
                        * (1.0 + 0.03 * (replicate - 2))
                        * np.exp(-0.20 * size_index)
                    )
                    predicted = observed * (
                        0.82 + 0.04 * size_index + 0.03 * time_index
                    )
                    rows.append(
                        {
                            "condition": exposure.condition,
                            "time_h": float(time_h),
                            "replicate": replicate,
                            "size_bin_center_nm": size_nm,
                            "treatment_particles_per_ml": observed,
                            "predicted_particles_per_ml": predicted,
                        }
                    )
    return pd.DataFrame(rows)


def test_size_time_visualizations_render_headlessly(tmp_path: Path) -> None:
    prediction = _synthetic_size_time_prediction()
    surface_path = tmp_path / "size_time_surface_overlay.png"
    error_path = tmp_path / "size_time_error_contours.png"

    plot_size_time_surface_overlay(prediction, surface_path, "synthetic")
    plot_size_time_error_contours(prediction, error_path, "synthetic")

    for output in (surface_path, error_path):
        assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert output.stat().st_size > 1_000
