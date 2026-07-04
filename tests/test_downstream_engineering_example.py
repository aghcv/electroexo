from __future__ import annotations

from pathlib import Path

from examples.compare_downstream_engineering import build_response_table, write_outputs


def test_downstream_engineering_example_outputs(tmp_path: Path) -> None:
    summary = build_response_table()
    assert {"potency_score", "purity_score", "optimization_objective"}.issubset(summary.columns)
    assert summary["optimization_objective"].max() > summary["optimization_objective"].min()

    write_outputs(summary, tmp_path, make_plots=False)
    assert (tmp_path / "downstream_engineering_summary.csv").exists()
    assert (tmp_path / "abbreviations.json").exists()
