from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import pandas as pd

from electro_exocytosis.abbreviations import STANDARD_ABBREVIATIONS
from electro_exocytosis.config import (
    CellStateConfig,
    ExposureConfig,
    PulseConfig,
    ScenarioConfig,
    SimulationConfig,
    SimulationScenario,
)
from electro_exocytosis.simulation import Simulation
from electro_exocytosis.visualization.style import (
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    line_styles,
    save_manuscript_figure,
)


@dataclass(frozen=True)
class EVScenarioSpec:
    name: str
    parameter_overrides: dict = field(default_factory=dict)
    cell_state_overrides: dict = field(default_factory=dict)
    pulse_overrides: dict = field(default_factory=dict)


SCENARIOS = (
    EVScenarioSpec(name="baseline_ev_biogenesis"),
    EVScenarioSpec(
        name="rab_escrt_limited",
        parameter_overrides={
            "ev_release": {
                "rab_conversion_baseline": 0.15,
                "escrt_baseline": 0.20,
                "rab_docking_baseline": 0.15,
                "baseline_sEV_rate": 0.70,
            }
        },
    ),
    EVScenarioSpec(
        name="ceramide_secretory_bias",
        parameter_overrides={
            "ev_release": {
                "ceramide_baseline": 0.50,
                "ceramide_PS_weight": 0.45,
                "acidification_ceramide_relief": 0.75,
                "baseline_sEV_rate": 1.20,
            }
        },
    ),
    EVScenarioSpec(
        name="injury_apoptotic_shift",
        parameter_overrides={
            "ev_release": {
                "baseline_AB_rate": 0.03,
                "apoptosis_damage_weight": 0.65,
                "apoptosis_ATP_loss_weight": 0.30,
            }
        },
        cell_state_overrides={"stress_sensitivity_modifier": 1.6},
        pulse_overrides={"amplitude_kV_cm": 24.0, "pulse_number": 40},
    ),
)


def build_response_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run EV-biogenesis comparison scenarios."""
    summary_rows: list[dict[str, float | str]] = []
    timeseries_frames: list[pd.DataFrame] = []

    for spec in SCENARIOS:
        result = Simulation(_build_scenario(spec), params_override=spec.parameter_overrides).run()
        summary_rows.append(
            {
                "scenario": spec.name,
                "peak_secretory_bias": result.summary["peak_secretory_bias"],
                "peak_lysosomal_routing": result.summary["peak_lysosomal_routing"],
                "peak_docked_MVB_pool": result.summary["peak_docked_MVB_pool"],
                "peak_budding_pool": result.summary["peak_budding_pool"],
                "peak_apoptotic_commitment": result.summary["peak_apoptotic_commitment"],
                "cumulative_small_EV": result.summary["cumulative_small_EV"],
                "cumulative_medium_large_EV": result.summary["cumulative_medium_large_EV"],
                "cumulative_apoptotic_body": result.summary["cumulative_apoptotic_body"],
                "viability_fraction": result.summary["viability_fraction"],
            }
        )
        frame = result.ev_timeseries[
            [
                "t",
                "MVB_pool",
                "ILV_load",
                "docked_MVB_pool",
                "budding_pool",
                "apoptotic_commitment",
                "secretory_bias",
                "lysosomal_routing",
                "sEV_rate",
                "mlEV_rate",
                "AB_rate",
            ]
        ].copy()
        frame.insert(0, "scenario", spec.name)
        timeseries_frames.append(frame)

    return pd.DataFrame(summary_rows), pd.concat(timeseries_frames, ignore_index=True)


def write_outputs(summary: pd.DataFrame, timeseries: pd.DataFrame, outdir: Path, make_plots: bool = True) -> None:
    """Write EV-biogenesis comparison outputs."""
    outdir.mkdir(parents=True, exist_ok=True)
    STANDARD_ABBREVIATIONS.rename_columns(summary).to_csv(outdir / "ev_biogenesis_summary.csv", index=False)
    STANDARD_ABBREVIATIONS.rename_columns(timeseries).to_csv(
        outdir / "ev_biogenesis_timeseries.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.write_bundle(outdir, keys=("EV", "MVB", "ILV", "ESCRT", "sEV", "m/lEV", "AB"))
    if make_plots:
        _plot_summary(summary, outdir)
        _plot_timeseries(timeseries, outdir)


def _build_scenario(spec: EVScenarioSpec) -> SimulationScenario:
    pulse_kwargs = {
        "amplitude_kV_cm": 15.0,
        "pulse_width_ns": 200.0,
        "pulse_number": 20,
        "repetition_rate_Hz": 5.0,
    }
    pulse_kwargs.update(spec.pulse_overrides)
    cell_state_kwargs = {"cell_type": "generic"}
    cell_state_kwargs.update(spec.cell_state_overrides)

    return SimulationScenario(
        scenario=ScenarioConfig(name=spec.name),
        pulse=PulseConfig(**pulse_kwargs),
        exposure=ExposureConfig(dosimetry_model="joule_lumped_thermal"),
        cell_state=CellStateConfig(**cell_state_kwargs),
        simulation=SimulationConfig(t_start_s=0.0, t_end_s=1800.0, output_dt_s=5.0),
    )


def _plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = [
        "cumulative_small_EV",
        "cumulative_medium_large_EV",
        "cumulative_apoptotic_body",
        "peak_secretory_bias",
        "peak_lysosomal_routing",
    ]
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(len(metrics), include_markers=False)
    x = range(len(summary))
    for index, metric in enumerate(metrics):
        values = summary[metric] / max(float(summary[metric].max()), 1e-12)
        ax.plot(x, values, label=STANDARD_ABBREVIATIONS.plot_label(metric), **styles[index])
    ax.set_xticks(list(x))
    ax.set_xticklabels(summary["scenario"], rotation=20, ha="right")
    ax.set_ylabel("Normalized response")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(
        fig,
        outdir / "ev_biogenesis_summary.png",
        abbreviation_keys=("EV", "MVB", "ESCRT", "sEV", "m/lEV", "AB"),
    )
    plt.close(fig)


def _plot_timeseries(timeseries: pd.DataFrame, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(timeseries["scenario"].nunique())
    for index, (scenario, frame) in enumerate(timeseries.groupby("scenario", sort=False)):
        style = styles[index]
        axes[0, 0].plot(frame["t"], frame["MVB_pool"], label=scenario, **style)
        axes[0, 1].plot(frame["t"], frame["secretory_bias"], label=scenario, **style)
        axes[1, 0].plot(frame["t"], frame["sEV_rate"], label=scenario, **style)
        axes[1, 1].plot(frame["t"], frame["AB_rate"], label=scenario, **style)

    axes[0, 0].set_ylabel("Multivesicular-body pool")
    axes[0, 1].set_ylabel("Secretory bias")
    axes[1, 0].set_ylabel("Small-EV release rate")
    axes[1, 1].set_ylabel("Apoptotic-body release rate")
    for ax in axes[1, :]:
        ax.set_xlabel("Time (s)")
    axes[0, 0].legend()
    fig.tight_layout()
    save_manuscript_figure(
        fig,
        outdir / "ev_biogenesis_timeseries.png",
        abbreviation_keys=("EV", "MVB", "AB", "sEV"),
    )
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/ev_biogenesis_comparison"))
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary, timeseries = build_response_table()
    write_outputs(summary, timeseries, args.out, make_plots=not args.no_plots)
    print(f"Wrote EV-biogenesis comparison outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
