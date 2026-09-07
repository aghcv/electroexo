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
    MANUSCRIPT_SINGLE_COLUMN_WIDTH_IN,
    line_styles,
    manuscript_style_context,
    save_manuscript_figure,
)

SCENARIO_LABELS = {
    "baseline_100ns": "Reference",
    "wide_pulse_high_dose": "Higher dose",
    "calcium_limited_medium": r"Low-Ca$^{2+}$ medium",
    "osmotic_recovery_slow": "Slow ion recovery",
}


@dataclass(frozen=True)
class Layer3ScenarioSpec:
    name: str
    amplitude_kV_cm: float
    pulse_width_ns: float
    pulse_number: int
    repetition_rate_Hz: float
    parameter_overrides: dict | None = None


SCENARIOS = (
    Layer3ScenarioSpec(
        name="baseline_100ns",
        amplitude_kV_cm=10.0,
        pulse_width_ns=100.0,
        pulse_number=10,
        repetition_rate_Hz=1.0,
    ),
    Layer3ScenarioSpec(
        name="wide_pulse_high_dose",
        amplitude_kV_cm=20.0,
        pulse_width_ns=300.0,
        pulse_number=20,
        repetition_rate_Hz=5.0,
    ),
    Layer3ScenarioSpec(
        name="calcium_limited_medium",
        amplitude_kV_cm=10.0,
        pulse_width_ns=100.0,
        pulse_number=10,
        repetition_rate_Hz=1.0,
        parameter_overrides={"ion_transport": {"Ca_ext_uM": 5.0}},
    ),
    Layer3ScenarioSpec(
        name="osmotic_recovery_slow",
        amplitude_kV_cm=10.0,
        pulse_width_ns=100.0,
        pulse_number=10,
        repetition_rate_Hz=1.0,
        parameter_overrides={"ion_transport": {"tau_ion_recovery_s": 1200.0}},
    ),
)


def build_response_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run comparison scenarios and return summary and time-series tables."""
    summary_rows: list[dict[str, float | str]] = []
    timeseries_frames: list[pd.DataFrame] = []

    for spec in SCENARIOS:
        scenario = _build_scenario(spec)
        result = Simulation(scenario, params_override=spec.parameter_overrides).run()
        state = result.state_timeseries
        summary_rows.append(
            {
                "scenario": spec.name,
                "amplitude_kV_cm": spec.amplitude_kV_cm,
                "pulse_width_ns": spec.pulse_width_ns,
                "pulse_number": spec.pulse_number,
                "repetition_rate_Hz": spec.repetition_rate_Hz,
                "peak_ca_i_uM": result.summary["peak_ca_i"],
                "peak_ca_mito_uM": result.summary["peak_ca_mito"],
                "peak_ros": result.summary["peak_ros"],
                "min_atp": result.summary["min_atp"],
                "min_mitochondrial_potential": result.summary[
                    "min_mitochondrial_potential"
                ],
                "peak_osmotic_stress": result.summary["peak_osmotic_stress"],
                "peak_pore_activation": float(state["pore_activation"].max()),
                "peak_J_Ca_pore_uM_s": float(state["J_Ca_pore"].max()),
                "peak_J_ER_release_uM_s": float(state["J_ER_release"].max()),
                "terminal_damage": float(state["damage"].iloc[-1]),
                "cumulative_small_EV": result.summary["cumulative_small_EV"],
                "viability_fraction": result.summary["viability_fraction"],
            }
        )
        frame = state[
            [
                "t",
                "Ca_i",
                "Ca_ER",
                "Ca_mito",
                "mitochondrial_potential",
                "ROS",
                "ATP",
                "Na_i",
                "K_i",
                "Cl_i",
                "osmotic_stress",
                "damage",
            ]
        ].copy()
        frame.insert(0, "scenario", spec.name)
        timeseries_frames.append(frame)

    return pd.DataFrame(summary_rows), pd.concat(timeseries_frames, ignore_index=True)


def write_outputs(
    summary: pd.DataFrame,
    timeseries: pd.DataFrame,
    outdir: Path,
    make_plots: bool = True,
) -> None:
    """Write Layer 3 comparison outputs."""
    outdir.mkdir(parents=True, exist_ok=True)
    STANDARD_ABBREVIATIONS.rename_columns(summary).to_csv(
        outdir / "ion_transport_bioenergetics_summary.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.rename_columns(timeseries).to_csv(
        outdir / "ion_transport_bioenergetics_timeseries.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.write_bundle(outdir, keys=("Ca_i", "ROS", "ATP", "EV", "ER"))
    if make_plots:
        _plot_peak_responses(summary, outdir)
        _plot_timeseries(timeseries, outdir)


def _build_scenario(spec: Layer3ScenarioSpec) -> SimulationScenario:
    return SimulationScenario(
        scenario=ScenarioConfig(name=spec.name),
        pulse=PulseConfig(
            amplitude_kV_cm=spec.amplitude_kV_cm,
            pulse_width_ns=spec.pulse_width_ns,
            pulse_number=spec.pulse_number,
            repetition_rate_Hz=spec.repetition_rate_Hz,
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
        "peak_ca_i_uM",
        "peak_ca_mito_uM",
        "peak_ros",
        "peak_osmotic_stress",
        "min_atp",
        "min_mitochondrial_potential",
    ]
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(len(metrics), include_markers=False)
    x = range(len(summary))
    for index, metric in enumerate(metrics):
        values = summary[metric] / max(float(summary[metric].max()), 1e-12)
        ax.plot(
            x, values, label=STANDARD_ABBREVIATIONS.plot_label(metric), **styles[index]
        )
    ax.set_xticks(list(x))
    ax.set_xticklabels(summary["scenario"], rotation=20, ha="right")
    ax.set_ylabel("Normalized response")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(
        fig,
        outdir / "layer3_peak_response_comparison.png",
        abbreviation_keys=("Ca_i", "ROS", "ATP", "EV"),
    )
    plt.close(fig)


def _plot_timeseries(timeseries: pd.DataFrame, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with manuscript_style_context():
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(MANUSCRIPT_SINGLE_COLUMN_WIDTH_IN, 4.15),
            sharex=True,
        )
        styles = line_styles(timeseries["scenario"].nunique())
        for index, (scenario, frame) in enumerate(
            timeseries.groupby("scenario", sort=False)
        ):
            style = styles[index]
            label = SCENARIO_LABELS.get(scenario, scenario.replace("_", " ").title())
            axes[0, 0].plot(frame["t"], frame["Ca_i"], label=label, **style)
            axes[0, 1].plot(frame["t"], frame["osmotic_stress"], label=label, **style)
            axes[1, 0].plot(
                frame["t"], frame["mitochondrial_potential"], label=label, **style
            )
            axes[1, 1].plot(frame["t"], frame["ROS"], label=label, **style)

        panel_titles = (
            r"Cytosolic Ca$^{2+}$",
            "Osmotic stress",
            "Mitochondrial\npotential",
            "Reactive oxygen\nspecies",
        )
        y_labels = (r"$\mu$M", "Relative", "Relative", "Relative")
        for axis, title, y_label in zip(axes.flat, panel_titles, y_labels, strict=True):
            axis.set_title(title, fontsize=9, pad=3)
            axis.set_ylabel(y_label, fontsize=9)
        for ax in axes[1, :]:
            ax.set_xlabel("Time (s)", fontsize=9)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
            ncol=2,
            fontsize=8,
            frameon=True,
            facecolor="white",
            edgecolor="#555555",
            framealpha=0.96,
            borderaxespad=0.0,
            columnspacing=0.9,
            handlelength=1.8,
        )
        for axis in axes.flat:
            axis.tick_params(labelsize=8.5)
            axis.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.6)
        fig.subplots_adjust(
            left=0.16,
            right=0.98,
            bottom=0.11,
            top=0.78,
            wspace=0.40,
            hspace=0.52,
        )
        save_manuscript_figure(fig, outdir / "layer3_timeseries_comparison.png")
        plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/ion_transport_bioenergetics_comparison"),
    )
    parser.add_argument(
        "--no-plots", action="store_true", help="Skip PNG plot generation."
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary, timeseries = build_response_table()
    write_outputs(summary, timeseries, args.out, make_plots=not args.no_plots)
    print(f"Wrote Layer 3 ion/bioenergetics comparison outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
