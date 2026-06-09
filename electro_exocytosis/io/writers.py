from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from electro_exocytosis.simulation import SimulationResult



def save_summary_json(result: SimulationResult, outdir: Path) -> None:
    """Save summary metrics as JSON."""
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "summary.json").write_text(json.dumps(result.summary, indent=2), encoding="utf-8")



def save_timeseries_csv(result: SimulationResult, outdir: Path) -> None:
    """Save state time series as CSV."""
    outdir.mkdir(parents=True, exist_ok=True)
    result.state_timeseries.to_csv(outdir / "state_timeseries.csv", index=False)



def save_ev_outputs_csv(result: SimulationResult, outdir: Path) -> None:
    """Save EV outputs as CSV."""
    outdir.mkdir(parents=True, exist_ok=True)
    ev_frame = result.ev_timeseries.copy()
    ev_frame["viability_fraction"] = result.summary["viability_fraction"]
    ev_frame["quality_pass"] = result.parameters_used["terminal_quality"]["quality_pass"]
    ev_frame.to_csv(outdir / "ev_outputs.csv", index=False)



def save_parameters_yaml(result: SimulationResult, outdir: Path) -> None:
    """Save merged parameters as YAML."""
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "parameters_used.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(result.parameters_used, handle, sort_keys=True)



def save_run_metadata(result: SimulationResult, outdir: Path) -> None:
    """Save run metadata as JSON."""
    outdir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "scenario_name": result.scenario_name,
        "mode": result.mode,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "warnings": result.warnings,
        "n_timepoints": int(len(result.t_array)),
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
