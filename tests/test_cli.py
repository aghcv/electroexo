import json
from pathlib import Path

from typer.testing import CliRunner

from electro_exocytosis.cli import app

DEMO = Path(__file__).resolve().parents[1] / "examples" / "scenario_demo.yml"


def test_headless_demo_writes_standard_artifacts(tmp_path: Path) -> None:
    output = tmp_path / "demo"
    result = CliRunner().invoke(app, [str(DEMO), "--out", str(output), "--no-plots"])

    assert result.exit_code == 0, result.output
    assert "release_demo" in result.output
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["scenario_name"] == "release_demo"
    assert metadata["n_timepoints"] == 7
    for name in ("summary.json", "state_timeseries.csv", "ev_outputs.csv", "parameters_used.yaml"):
        assert (output / name).is_file()
    assert not (output / "plots").exists()


def test_invalid_scenario_exits_before_creating_run_directory(tmp_path: Path) -> None:
    scenario = tmp_path / "invalid.yml"
    scenario.write_text(
        DEMO.read_text(encoding="utf-8").replace("output_dt_s: 5.0", "output_dt_s: 0"),
        encoding="utf-8",
    )
    output = tmp_path / "invalid"
    result = CliRunner().invoke(app, [str(scenario), "--out", str(output), "--no-plots"])

    assert result.exit_code == 2
    assert "output_dt_s" in result.output
    assert not output.exists()
