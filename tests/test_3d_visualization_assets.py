from __future__ import annotations

from pathlib import Path

from electro_exocytosis.config import (
    CellStateConfig,
    ExposureConfig,
    PulseConfig,
    ScenarioConfig,
    SimulationConfig,
    SimulationScenario,
)
from electro_exocytosis.simulation import Simulation
from electro_exocytosis.visualization.three_d import (
    build_default_event_schedule,
    load_human_cell_manifest,
    prepare_human_cell_vtp_assets,
    scene_values_from_result,
    write_event_schedule_metadata,
)


def test_prepare_human_cell_vtp_assets_extracts_referenced_blocks(tmp_path: Path) -> None:
    source_dir = Path(__file__).resolve().parents[1] / "docs"
    asset_dir = tmp_path / "human_cell_vtp"

    manifest = prepare_human_cell_vtp_assets(
        source_dir,
        asset_dir,
        components=("surface", "mitocondria"),
    )
    loaded = load_human_cell_manifest(asset_dir)

    assert manifest.file_count == 4
    assert loaded.file_count == 4
    assert {component.name for component in loaded.components} == {"surface", "mitocondria"}
    for component in loaded.components:
        for rel_path in component.files:
            assert (asset_dir / rel_path).exists()


def test_scene_values_from_result_returns_normalized_controls() -> None:
    scenario = SimulationScenario(
        scenario=ScenarioConfig(name="visualization_smoke"),
        pulse=PulseConfig(
            amplitude_kV_cm=15.0,
            pulse_width_ns=200.0,
            pulse_number=5,
            repetition_rate_Hz=5.0,
        ),
        exposure=ExposureConfig(dosimetry_model="joule_lumped_thermal"),
        cell_state=CellStateConfig(),
        simulation=SimulationConfig(t_start_s=0.0, t_end_s=30.0, output_dt_s=5.0),
    )
    result = Simulation(scenario).run()

    values = scene_values_from_result(result)

    assert values["scenario"] == "visualization_smoke"
    for key in (
        "membrane_activation",
        "calcium_signal",
        "ps_exposure",
        "repair_state",
        "ros_signal",
        "mitochondrial_stress",
        "ev_release_signal",
        "viability_fraction",
    ):
        assert 0.0 <= float(values[key]) <= 1.0


def test_default_event_schedule_is_ordered_and_serializable(tmp_path: Path) -> None:
    scenario = SimulationScenario(
        scenario=ScenarioConfig(name="visualization_schedule_smoke"),
        pulse=PulseConfig(
            amplitude_kV_cm=15.0,
            pulse_width_ns=200.0,
            pulse_number=5,
            repetition_rate_Hz=5.0,
        ),
        exposure=ExposureConfig(dosimetry_model="joule_lumped_thermal"),
        cell_state=CellStateConfig(),
        simulation=SimulationConfig(t_start_s=0.0, t_end_s=30.0, output_dt_s=5.0),
    )
    result = Simulation(scenario).run()

    schedule = build_default_event_schedule(result)
    written_schedule = write_event_schedule_metadata(result, tmp_path / "storyboard_schedule.json")

    assert len(schedule) == 6
    assert schedule == written_schedule
    assert schedule[0].event_id == "01_nsPEF_exposure"
    assert schedule[-1].event_id == "06_EV_release"
    assert {event.signal_key for event in schedule} == {
        "field_exposure",
        "membrane_electrodynamics",
        "calcium_transport",
        "mitochondrial_bioenergetics",
        "remodeling_repair",
        "ev_release",
    }
    assert (tmp_path / "storyboard_schedule.json").exists()
