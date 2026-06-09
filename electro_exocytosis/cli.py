from __future__ import annotations

from pathlib import Path

import typer

from electro_exocytosis.evidence.evidence_loader import EvidenceLoader
from electro_exocytosis.io.readers import load_scenario
from electro_exocytosis.io.writers import (
    save_ev_outputs_csv,
    save_parameters_yaml,
    save_run_metadata,
    save_summary_json,
    save_timeseries_csv,
)
from electro_exocytosis.simulation import Simulation
from electro_exocytosis.visualization.plots import generate_all_plots

app = typer.Typer()


@app.command()
def run(
    scenario_file: str,
    out: str = "results/output",
    no_plots: bool = False,
    evidence_file: str | None = None,
):
    """Run an nsPEF EV release simulation from a YAML scenario file."""
    scenario = load_scenario(scenario_file)
    simulation = Simulation(scenario)
    result = simulation.run()

    outdir = Path(out)
    outdir.mkdir(parents=True, exist_ok=True)
    save_summary_json(result, outdir)
    save_timeseries_csv(result, outdir)
    save_ev_outputs_csv(result, outdir)
    save_parameters_yaml(result, outdir)
    save_run_metadata(result, outdir)

    if evidence_file:
        loader = EvidenceLoader(evidence_file)
        if loader.load():
            typer.echo(loader.summarize())
        else:
            typer.echo(f"Evidence workbook not found: {evidence_file}")

    if not no_plots:
        plots_dir = outdir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        generate_all_plots(result, plots_dir)

    typer.echo(f"Completed simulation '{result.scenario_name}' -> {outdir}")
