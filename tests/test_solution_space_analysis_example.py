from __future__ import annotations

from pathlib import Path

from examples.solution_space_analysis.generate_solution_space import (
    build_solution_space_dataset,
    build_solution_space_specs,
    write_outputs,
    write_scenario_files,
)


def test_solution_space_analysis_small_grid(tmp_path: Path) -> None:
    specs = build_solution_space_specs(
        amplitudes_kV_cm=(8.0, 20.0),
        pulse_widths_ns=(100.0,),
        pulse_numbers=(5, 30),
    )

    scenario_paths = write_scenario_files(specs, tmp_path / "scenarios")
    summary, results = build_solution_space_dataset(specs)

    assert len(scenario_paths) == len(specs)
    assert len(results) == len(specs)
    assert {
        "scenario_name",
        "scenario_family",
        "dose_index",
        "total_extracellular_vesicle_output",
        "viability_fraction",
        "optimization_objective",
    }.issubset(summary.columns)

    low = summary[
        (summary["scenario_family"] == "dose_grid")
        & (summary["amplitude_kV_cm"] == 8.0)
        & (summary["pulse_number"] == 5)
    ].iloc[0]
    high = summary[
        (summary["scenario_family"] == "dose_grid")
        & (summary["amplitude_kV_cm"] == 20.0)
        & (summary["pulse_number"] == 30)
    ].iloc[0]
    assert high["dose_index"] > low["dose_index"]
    assert high["total_extracellular_vesicle_output"] >= low["total_extracellular_vesicle_output"]

    write_outputs(summary, results, tmp_path / "results", make_plots=False)
    assert (tmp_path / "results" / "solution_space_summary.csv").exists()
    assert (tmp_path / "results" / "top_scenarios_by_objective.csv").exists()
    assert (tmp_path / "results" / "analysis_manifest.json").exists()
