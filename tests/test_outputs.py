from __future__ import annotations

from pathlib import Path

from electro_exocytosis.io.readers import load_scenario
from electro_exocytosis.io.writers import (
    save_ev_outputs_csv,
    save_parameters_yaml,
    save_run_metadata,
    save_summary_json,
    save_timeseries_csv,
)
from electro_exocytosis.simulation import Simulation



def test_writers_create_outputs(tmp_path: Path) -> None:
    scenario = load_scenario(Path("examples/scenario_baseline.yaml"))
    result = Simulation(scenario).run()
    save_summary_json(result, tmp_path)
    save_timeseries_csv(result, tmp_path)
    save_ev_outputs_csv(result, tmp_path)
    save_parameters_yaml(result, tmp_path)
    save_run_metadata(result, tmp_path)
    for name in [
        "summary.json",
        "state_timeseries.csv",
        "ev_outputs.csv",
        "parameters_used.yaml",
        "run_metadata.json",
    ]:
        assert (tmp_path / name).exists()
