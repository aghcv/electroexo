from __future__ import annotations

from pathlib import Path

from electro_exocytosis.io.readers import load_scenario
from electro_exocytosis.simulation import Simulation



def test_solver_smoke() -> None:
    scenario = load_scenario(Path("examples/scenario_baseline.yaml"))
    result = Simulation(scenario).run()
    assert len(result.t_array) > 0
    assert not result.state_timeseries.empty
    assert not result.ev_timeseries.empty
    for key in [
        "scenario_name",
        "mode",
        "dose_index",
        "peak_ca_i",
        "cumulative_small_EV",
        "purity_score",
        "warnings",
    ]:
        assert key in result.summary
