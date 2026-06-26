from __future__ import annotations

from pathlib import Path

from examples.map_ev_release_regimes import build_regime_table, write_outputs


def test_ev_release_regime_map_small_grid(tmp_path: Path) -> None:
    summary = build_regime_table(amplitudes_kV_cm=(8.0, 20.0), pulse_numbers=(5, 30))

    assert len(summary) == 4
    assert {
        "amplitude_kV_cm",
        "pulse_number",
        "cumulative_small_EV",
        "cumulative_apoptotic_body",
        "viability_fraction",
        "dominant_subtype",
    }.issubset(summary.columns)

    low = summary[(summary["amplitude_kV_cm"] == 8.0) & (summary["pulse_number"] == 5)].iloc[0]
    high = summary[(summary["amplitude_kV_cm"] == 20.0) & (summary["pulse_number"] == 30)].iloc[0]

    assert high["dose_index"] > low["dose_index"]
    assert high["cumulative_apoptotic_body"] >= low["cumulative_apoptotic_body"]

    write_outputs(summary, tmp_path, make_plots=False)
    assert (tmp_path / "ev_release_regime_map.csv").exists()
