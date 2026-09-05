from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electro_exocytosis.experimental_bridge import (
    ExperimentalObservationBridge,
    aggregate_repeated_observations,
)


def test_repeated_observations_are_aggregated_without_losing_raw_rows() -> None:
    frame = pd.DataFrame(
        {
            "condition": ["control", "control", "treated", "treated", "treated"],
            "time_h": [1.0, 1.0, 1.0, 1.0, 1.0],
            "concentration": [10.0, 14.0, 20.0, 24.0, 28.0],
        }
    )
    original = frame.copy(deep=True)

    result = aggregate_repeated_observations(
        frame,
        group_columns=["condition", "time_h"],
        value_columns=["concentration"],
    )

    pd.testing.assert_frame_equal(frame, original)
    pd.testing.assert_frame_equal(result.raw_observations, original)
    assert result.raw_observations is not frame
    assert result.group_columns == ("condition", "time_h")
    assert result.value_columns == ("concentration",)
    assert list(result.summary) == [
        "condition",
        "time_h",
        "concentration_mean",
        "concentration_sd",
        "concentration_se",
        "concentration_n",
    ]
    control, treated = result.summary.itertuples(index=False)
    assert control.concentration_mean == pytest.approx(12.0)
    assert control.concentration_sd == pytest.approx(np.sqrt(8.0))
    assert control.concentration_se == pytest.approx(2.0)
    assert control.concentration_n == 2
    assert treated.concentration_mean == pytest.approx(24.0)
    assert treated.concentration_sd == pytest.approx(4.0)
    assert treated.concentration_se == pytest.approx(4.0 / np.sqrt(3.0))
    assert treated.concentration_n == 3


def test_aggregation_handles_missing_values_per_measurement_and_singletons() -> None:
    frame = pd.DataFrame(
        {
            "condition": ["control", "control", "treated", None],
            "signal_a": [1.0, np.nan, 3.0, 5.0],
            "signal_b": [2.0, 4.0, np.nan, 8.0],
        }
    )

    summary = aggregate_repeated_observations(
        frame,
        group_columns=["condition"],
        value_columns=["signal_a", "signal_b"],
    ).summary

    control = summary.iloc[0]
    assert control["signal_a_mean"] == pytest.approx(1.0)
    assert control["signal_a_n"] == 1
    assert np.isnan(control["signal_a_sd"])
    assert np.isnan(control["signal_a_se"])
    assert control["signal_b_mean"] == pytest.approx(3.0)
    assert control["signal_b_sd"] == pytest.approx(np.sqrt(2.0))
    assert control["signal_b_se"] == pytest.approx(1.0)
    assert control["signal_b_n"] == 2

    missing_key = summary.loc[summary["condition"].isna()].iloc[0]
    assert missing_key["signal_a_mean"] == pytest.approx(5.0)
    assert missing_key["signal_b_mean"] == pytest.approx(8.0)


@pytest.mark.parametrize(
    ("group_columns", "value_columns", "error_type", "match"),
    [
        ([], ["value"], ValueError, "group_columns"),
        (["group"], [], ValueError, "value_columns"),
        (["missing"], ["value"], ValueError, "Missing aggregation columns"),
        (["group"], ["group"], ValueError, "must be distinct"),
        ("group", ["value"], TypeError, "not a string"),
    ],
)
def test_aggregation_validates_column_selections(
    group_columns: list[str] | str,
    value_columns: list[str],
    error_type: type[Exception],
    match: str,
) -> None:
    frame = pd.DataFrame({"group": ["a"], "value": [1.0]})

    with pytest.raises(error_type, match=match):
        aggregate_repeated_observations(
            frame,
            group_columns=group_columns,
            value_columns=value_columns,
        )


@pytest.mark.parametrize(
    ("values", "error_type", "match"),
    [
        (["not numeric"], TypeError, "must be numeric"),
        ([np.inf], ValueError, "infinite observations"),
        ([True], TypeError, "must be numeric"),
    ],
)
def test_aggregation_rejects_invalid_numeric_observations(
    values: list[object], error_type: type[Exception], match: str
) -> None:
    frame = pd.DataFrame({"group": ["a"], "value": values})

    with pytest.raises(error_type, match=match):
        aggregate_repeated_observations(
            frame,
            group_columns=["group"],
            value_columns=["value"],
        )


def test_concentration_to_initial_cell_output_and_round_trip() -> None:
    bridge = ExperimentalObservationBridge(
        initial_cell_count=5.0e6,
        medium_volume_ml=5.0,
    )

    per_cell = bridge.concentration_to_particles_per_cell([1.0e9, 2.0e9])

    assert np.allclose(per_cell, [1000.0, 2000.0])
    assert np.allclose(
        bridge.particles_per_cell_to_concentration(per_cell), [1.0e9, 2.0e9]
    )
    assert bridge.effective_cell_density_per_ml == pytest.approx(1.0e6)


def test_corrections_and_viable_cell_basis_are_explicit() -> None:
    bridge = ExperimentalObservationBridge(
        initial_cell_count=5.0e6,
        medium_volume_ml=5.0,
        viability_fraction=0.5,
        recovery_fraction=0.5,
        dilution_factor=2.0,
        background_concentration_particles_per_ml=1.0e8,
        cell_basis="viable",
    )

    # (1.1e9 - 0.1e9) * 2 / 0.5 = 4e9 corrected particles/mL;
    # 2.5e6 viable cells in 5 mL = 0.5e6 cells/mL.
    assert bridge.concentration_to_particles_per_cell(1.1e9) == pytest.approx(8000.0)
    assert bridge.particles_per_cell_to_concentration(8000.0) == pytest.approx(1.1e9)


def test_yaml_and_dataframe_interface() -> None:
    path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "ffrci_experimental_bridge.yml"
    )
    bridge = ExperimentalObservationBridge.from_yaml(path)
    frame = pd.DataFrame({"concentration": [2.0e9]})

    transformed = bridge.transform_frame(frame, concentration_column="concentration")

    assert transformed["observed_particles_per_cell"].iloc[0] == pytest.approx(2000.0)
    assert bridge.to_metadata()["effective_cell_density_per_ml"] == pytest.approx(1.0e6)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"initial_cell_count": 0, "medium_volume_ml": 1},
        {"initial_cell_count": 1, "medium_volume_ml": 0},
        {"initial_cell_count": 1, "medium_volume_ml": 1, "viability_fraction": 0},
        {"initial_cell_count": 1, "medium_volume_ml": 1, "recovery_fraction": 1.1},
    ],
)
def test_invalid_bridge_settings_fail(kwargs: dict[str, float]) -> None:
    with pytest.raises(ValueError):
        ExperimentalObservationBridge(**kwargs)
