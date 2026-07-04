from __future__ import annotations

import argparse
from dataclasses import dataclass
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
class RepairScenarioSpec:
    name: str
    parameter_overrides: dict | None = None


SCENARIOS = (
    RepairScenarioSpec(name="baseline_repair"),
    RepairScenarioSpec(
        name="buffered_microdomain",
        parameter_overrides={"remodeling_repair": {"microdomain_gain": 0.5, "microdomain_pore_gain": 0.5}},
    ),
    RepairScenarioSpec(
        name="calpain_inhibited",
        parameter_overrides={"remodeling_repair": {"K_calpain_uM": 8.0, "actin_calpain_weight": 0.2}},
    ),
    RepairScenarioSpec(
        name="annexin_lysosome_strong",
        parameter_overrides={
            "remodeling_repair": {
                "K_annex_uM": 0.2,
                "K_lysosomal_repair_uM": 0.25,
                "resealing_annexin_weight": 0.6,
                "resealing_lysosome_weight": 0.45,
            }
        },
    ),
)


def build_response_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run remodeling/repair comparison scenarios."""
    summary_rows: list[dict[str, float | str]] = []
    timeseries_frames: list[pd.DataFrame] = []

    for spec in SCENARIOS:
        result = Simulation(_build_scenario(spec), params_override=spec.parameter_overrides).run()
        state = result.state_timeseries
        summary_rows.append(
            {
                "scenario": spec.name,
                "peak_ca_submembrane_uM": float(state["Ca_submembrane"].max()),
                "peak_ps_exposure": result.summary["peak_ps_exposure"],
                "peak_calpain_activity": result.summary["peak_calpain_activity"],
                "peak_annexin_activity": result.summary["peak_annexin_activity"],
                "peak_actin_disruption": result.summary["peak_actin_disruption"],
                "peak_repair_state": result.summary["peak_repair_state"],
                "cumulative_repair_shedding": result.summary["cumulative_repair_shedding"],
                "cumulative_medium_large_EV": result.summary["cumulative_medium_large_EV"],
                "viability_fraction": result.summary["viability_fraction"],
            }
        )
        frame = state[
            [
                "t",
                "Ca_i",
                "Ca_submembrane",
                "PS_exposure",
                "calpain_activity",
                "annexin_activity",
                "lysosomal_repair_activity",
                "actin_disruption",
                "actomyosin_tension",
                "repair_state",
                "repair_shedding_rate",
            ]
        ].copy()
        frame.insert(0, "scenario", spec.name)
        timeseries_frames.append(frame)

    return pd.DataFrame(summary_rows), pd.concat(timeseries_frames, ignore_index=True)


def write_outputs(summary: pd.DataFrame, timeseries: pd.DataFrame, outdir: Path, make_plots: bool = True) -> None:
    """Write remodeling/repair comparison outputs."""
    outdir.mkdir(parents=True, exist_ok=True)
    STANDARD_ABBREVIATIONS.rename_columns(summary).to_csv(
        outdir / "remodeling_repair_summary.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.rename_columns(timeseries).to_csv(
        outdir / "remodeling_repair_timeseries.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.write_bundle(outdir, keys=("PS", "EV"))
    if make_plots:
        _plot_peak_responses(summary, outdir)
        _plot_timeseries(timeseries, outdir)


def _build_scenario(spec: RepairScenarioSpec) -> SimulationScenario:
    return SimulationScenario(
        scenario=ScenarioConfig(name=spec.name),
        pulse=PulseConfig(
            amplitude_kV_cm=15.0,
            pulse_width_ns=200.0,
            pulse_number=20,
            repetition_rate_Hz=5.0,
        ),
        exposure=ExposureConfig(dosimetry_model="joule_lumped_thermal"),
        cell_state=CellStateConfig(),
        simulation=SimulationConfig(t_start_s=0.0, t_end_s=1800.0, output_dt_s=5.0),
    )


def _plot_peak_responses(summary: pd.DataFrame, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = [
        "peak_ps_exposure",
        "peak_calpain_activity",
        "peak_annexin_activity",
        "peak_actin_disruption",
        "peak_repair_state",
        "cumulative_repair_shedding",
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
    save_manuscript_figure(fig, outdir / "remodeling_repair_peak_response.png", abbreviation_keys=("PS",))
    plt.close(fig)


def _plot_timeseries(timeseries: pd.DataFrame, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(timeseries["scenario"].nunique())
    for index, (scenario, frame) in enumerate(timeseries.groupby("scenario", sort=False)):
        style = styles[index]
        axes[0, 0].plot(frame["t"], frame["Ca_submembrane"], label=scenario, **style)
        axes[0, 1].plot(frame["t"], frame["PS_exposure"], label=scenario, **style)
        axes[1, 0].plot(frame["t"], frame["actin_disruption"], label=scenario, **style)
        axes[1, 1].plot(frame["t"], frame["repair_state"], label=scenario, **style)

    axes[0, 0].set_ylabel("Submembrane calcium (uM)")
    axes[0, 1].set_ylabel("Phosphatidylserine exposure")
    axes[1, 0].set_ylabel("Actin disruption")
    axes[1, 1].set_ylabel("Resealing state")
    for ax in axes[1, :]:
        ax.set_xlabel("Time (s)")
    axes[0, 0].legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "remodeling_repair_timeseries.png", abbreviation_keys=("PS",))
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/remodeling_repair_comparison"))
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary, timeseries = build_response_table()
    write_outputs(summary, timeseries, args.out, make_plots=not args.no_plots)
    print(f"Wrote remodeling/repair comparison outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
