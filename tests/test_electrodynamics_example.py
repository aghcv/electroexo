from __future__ import annotations

from pathlib import Path

from examples.compare_membrane_electrodynamics import build_response_table, write_outputs


def test_electrodynamics_example_builds_summary(tmp_path: Path) -> None:
    summary = build_response_table()

    assert not summary.empty
    assert {
        "amplitude_kV_cm",
        "pulse_width_ns",
        "cell_radius_um",
        "schwan_limit_V",
        "membrane_charging_factor",
        "delta_Vm_V",
        "membrane_permeability",
        "pore_density_um2",
    }.issubset(summary.columns)
    assert summary["delta_Vm_V"].max() > summary["delta_Vm_V"].min()

    write_outputs(summary, tmp_path, make_plots=False)
    assert (tmp_path / "electrodynamics_response_summary.csv").exists()
