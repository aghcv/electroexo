from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electro_exocytosis.experimental_bridge import ExperimentalObservationBridge


def test_concentration_to_initial_cell_output_and_round_trip() -> None:
    bridge = ExperimentalObservationBridge(
        initial_cell_count=5.0e6,
        medium_volume_ml=5.0,
    )

    per_cell = bridge.concentration_to_particles_per_cell([1.0e9, 2.0e9])

    assert np.allclose(per_cell, [1000.0, 2000.0])
    assert np.allclose(bridge.particles_per_cell_to_concentration(per_cell), [1.0e9, 2.0e9])
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
    path = Path(__file__).resolve().parents[1] / "examples" / "ffrci_experimental_bridge.yml"
    bridge = ExperimentalObservationBridge.from_yaml(path)
    frame = pd.DataFrame({"concentration": [2.0e9]})

    transformed = bridge.transform_frame(
        frame, concentration_column="concentration"
    )

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
