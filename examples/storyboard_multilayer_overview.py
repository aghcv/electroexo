from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
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
from electro_exocytosis.simulation import Simulation, SimulationResult
from electro_exocytosis.visualization.style import (
    MANUSCRIPT_DOUBLE_COLUMN_WIDTH_IN,
    manuscript_style_context,
    save_manuscript_figure,
)


@dataclass(frozen=True)
class StoryScenarioSpec:
    name: str
    label: str
    pulse_overrides: dict = field(default_factory=dict)
    exposure_overrides: dict = field(default_factory=dict)
    cell_state_overrides: dict = field(default_factory=dict)
    parameter_overrides: dict = field(default_factory=dict)
    narrative: str = ""


SCENARIOS = (
    StoryScenarioSpec(
        name="mild_reversible_window",
        label="Mild reversible window",
        pulse_overrides={
            "amplitude_kV_cm": 10.0,
            "pulse_width_ns": 100.0,
            "pulse_number": 10,
            "repetition_rate_Hz": 1.0,
        },
        narrative="Near-threshold charging with limited downstream amplification.",
    ),
    StoryScenarioSpec(
        name="productive_secretory_window",
        label="Productive secretory window",
        pulse_overrides={
            "amplitude_kV_cm": 15.0,
            "pulse_width_ns": 200.0,
            "pulse_number": 20,
            "repetition_rate_Hz": 5.0,
        },
        parameter_overrides={
            "remodeling_repair": {
                "K_annex_uM": 0.2,
                "K_lysosomal_repair_uM": 0.25,
                "resealing_annexin_weight": 0.6,
                "resealing_lysosome_weight": 0.45,
            },
            "ev_release": {
                "ceramide_baseline": 0.50,
                "ceramide_PS_weight": 0.45,
                "acidification_ceramide_relief": 0.75,
                "baseline_sEV_rate": 1.20,
            },
        },
        narrative="Coordinated Ca2+, repair, and secretory routing bias favor small-EV release.",
    ),
    StoryScenarioSpec(
        name="injury_apoptotic_window",
        label="Injury-dominant window",
        pulse_overrides={
            "amplitude_kV_cm": 24.0,
            "pulse_width_ns": 200.0,
            "pulse_number": 40,
            "repetition_rate_Hz": 5.0,
        },
        cell_state_overrides={"stress_sensitivity_modifier": 1.6},
        parameter_overrides={
            "ev_release": {
                "baseline_AB_rate": 0.03,
                "apoptosis_damage_weight": 0.65,
                "apoptosis_ATP_loss_weight": 0.30,
            }
        },
        narrative="Stronger electroporation and stress coupling overwhelm repair and shift output toward apoptotic vesiculation.",
    ),
)


def build_storyboard_dataset(
    scenarios: Sequence[StoryScenarioSpec] = SCENARIOS,
) -> tuple[pd.DataFrame, list[tuple[StoryScenarioSpec, SimulationResult]]]:
    """Run three illustrative scenarios spanning mild, productive, and injurious regimes."""
    rows: list[dict[str, float | str]] = []
    outputs: list[tuple[StoryScenarioSpec, SimulationResult]] = []
    for spec in scenarios:
        result = Simulation(
            _build_scenario(spec), params_override=spec.parameter_overrides
        ).run()
        outputs.append((spec, result))
        rows.append(
            {
                "scenario": spec.name,
                "label": spec.label,
                "narrative": spec.narrative,
                "dose_index": result.summary["dose_index"],
                "peak_ca_i": result.summary["peak_ca_i"],
                "peak_ros": result.summary["peak_ros"],
                "peak_repair_state": result.summary["peak_repair_state"],
                "peak_secretory_bias": result.summary["peak_secretory_bias"],
                "cumulative_small_EV": result.summary["cumulative_small_EV"],
                "cumulative_medium_large_EV": result.summary[
                    "cumulative_medium_large_EV"
                ],
                "cumulative_apoptotic_body": result.summary[
                    "cumulative_apoptotic_body"
                ],
                "viability_fraction": result.summary["viability_fraction"],
            }
        )
    return pd.DataFrame(rows), outputs


def write_outputs(
    summary: pd.DataFrame,
    outputs: list[tuple[StoryScenarioSpec, SimulationResult]],
    outdir: Path,
    make_plots: bool = True,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    STANDARD_ABBREVIATIONS.rename_columns(summary).to_csv(
        outdir / "multilayer_storyboard_summary.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.write_bundle(
        outdir,
        keys=(
            "nsPEF",
            "EV",
            "PM",
            "ER",
            "MVB",
            "ROS",
            "ATP",
            "PS",
            "sEV",
            "m/lEV",
            "AB",
        ),
    )
    if make_plots:
        plot_storyboard(outputs, outdir / "multilayer_storyboard.png")


def plot_storyboard(
    outputs: list[tuple[StoryScenarioSpec, SimulationResult]],
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_rows = len(outputs)
    n_cols = 5
    with manuscript_style_context():
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(MANUSCRIPT_DOUBLE_COLUMN_WIDTH_IN, 5.2),
        )
        if n_rows == 1:
            axes = np.array([axes])

        layer3_ca_max = max(
            float(result.state_timeseries["Ca_i"].max()) for _, result in outputs
        )
        layer3_ros_max = max(
            float(result.state_timeseries["ROS"].max()) for _, result in outputs
        )

        col_titles = [
            "L1\nDosimetry",
            "L2\nElectrodynamics",
            "L3\nCa/ROS/ATP",
            "L4\nRepair",
            "L5\nEV output",
        ]
        for col, title in enumerate(col_titles):
            axes[0, col].set_title(title, fontsize=9, pad=5)

        concise_row_labels = {
            "mild_reversible_window": "Mild",
            "productive_secretory_window": "Productive",
            "injury_apoptotic_window": "Injury",
        }
        row_labels = tuple(
            concise_row_labels.get(spec.name, spec.label) for spec, _ in outputs
        )
        for row, (_, result) in enumerate(outputs):
            _plot_layer1_panel(axes[row, 0], result)
            _plot_layer2_panel(axes[row, 1], result)
            _plot_layer3_panel(
                axes[row, 2],
                result,
                layer3_ca_max=layer3_ca_max,
                layer3_ros_max=layer3_ros_max,
            )
            _plot_layer4_panel(axes[row, 3], result)
            _plot_layer5_panel(axes[row, 4], result)
            for col in range(n_cols):
                if row != n_rows - 1:
                    axes[row, col].set_xlabel("")
                    axes[row, col].tick_params(axis="x", labelbottom=False)
                if row != n_rows // 2:
                    axes[row, col].set_ylabel("")

        fig.subplots_adjust(
            left=0.11,
            right=0.995,
            bottom=0.16,
            top=0.89,
            wspace=0.46,
            hspace=0.24,
        )

        for row, label in enumerate(row_labels):
            position = axes[row, 0].get_position()
            fig.text(
                0.018,
                (position.y0 + position.y1) / 2.0,
                label,
                rotation=90,
                va="center",
                ha="center",
                fontsize=8.5,
                fontweight="bold",
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": "white",
                    "edgecolor": "#666666",
                    "linewidth": 0.7,
                },
            )

        layer3_handles, layer3_labels = axes[0, 2].get_legend_handles_labels()
        layer4_handles, layer4_labels = axes[0, 3].get_legend_handles_labels()
        fig.legend(
            layer3_handles + layer4_handles,
            layer3_labels + layer4_labels,
            loc="lower center",
            bbox_to_anchor=(0.55, 0.018),
            ncol=6,
            fontsize=8,
            frameon=True,
            facecolor="white",
            edgecolor="#555555",
            framealpha=1.0,
            borderpad=0.35,
            handlelength=1.8,
            handletextpad=0.4,
            columnspacing=0.9,
        )
        save_manuscript_figure(fig, output_path)
        plt.close(fig)


def _build_scenario(spec: StoryScenarioSpec) -> SimulationScenario:
    pulse_kwargs = {
        "amplitude_kV_cm": 15.0,
        "pulse_width_ns": 200.0,
        "pulse_number": 20,
        "repetition_rate_Hz": 5.0,
    }
    pulse_kwargs.update(spec.pulse_overrides)
    exposure_kwargs = {"dosimetry_model": "joule_lumped_thermal"}
    exposure_kwargs.update(spec.exposure_overrides)
    cell_state_kwargs = {"cell_type": "generic"}
    cell_state_kwargs.update(spec.cell_state_overrides)
    return SimulationScenario(
        scenario=ScenarioConfig(name=spec.name),
        pulse=PulseConfig(**pulse_kwargs),
        exposure=ExposureConfig(**exposure_kwargs),
        cell_state=CellStateConfig(**cell_state_kwargs),
        simulation=SimulationConfig(t_start_s=0.0, t_end_s=1800.0, output_dt_s=5.0),
    )


def _plot_layer1_panel(ax, result: SimulationResult) -> None:
    dosimetry = result.parameters_used["dosimetry"]
    scenario = result.parameters_used["scenario"]
    descriptors = result.parameters_used["pulse_descriptors"]
    temp0 = float(scenario["exposure"]["temperature_C"])
    t_train = max(
        float(descriptors["train_duration_s"]), float(descriptors["pulse_width_s"])
    )
    horizon = t_train + max(t_train, 10.0)
    times = np.linspace(0.0, horizon, 160)
    tau = float(scenario["exposure"]["thermal_relaxation_time_s"])
    eta = float(scenario["exposure"]["thermal_efficiency"])
    adiabatic = float(dosimetry["adiabatic_temperature_rise_K"])
    heating_rate = adiabatic / max(t_train, 1e-12)
    rise = eta * heating_rate * tau * (1.0 - np.exp(-np.minimum(times, t_train) / tau))
    cooling = times > t_train
    rise[cooling] = rise[cooling] * np.exp(-(times[cooling] - t_train) / tau)
    ax.plot(times, temp0 + rise, color="#4477AA", linewidth=1.6)
    ax.axvline(t_train, color="0.5", linestyle=":", linewidth=1.0)
    ax.set_ylabel("Temperature (°C)", fontsize=8, labelpad=2)
    ax.set_xlabel("Time (s)", fontsize=8, labelpad=2)
    ax.locator_params(axis="both", nbins=3)
    _style_storyboard_axis(ax)


def _plot_layer2_panel(ax, result: SimulationResult) -> None:
    electro = result.parameters_used["electrodynamics"]
    labels = ["PM", "ER", "Mito", "MVB"]
    values = [
        float(electro["delta_Vm"]),
        float(electro["delta_V_ER"]),
        float(electro["delta_V_mito"]),
        float(electro["delta_V_MVB"]),
    ]
    x = np.arange(len(labels))
    ax.bar(x, values, color=["#111111", "#555555", "#888888", "#BBBBBB"], width=0.68)
    ax.axhline(0.25, color="#CC6677", linestyle=":", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Voltage (V)", fontsize=8, labelpad=2)
    ax.locator_params(axis="y", nbins=3)
    _style_storyboard_axis(ax)


def _plot_layer3_panel(
    ax, result: SimulationResult, *, layer3_ca_max: float, layer3_ros_max: float
) -> None:
    state = result.state_timeseries
    t = state["t"]
    ax.plot(
        t,
        state["Ca_i"] / max(layer3_ca_max, 1e-12),
        color="#4477AA",
        linewidth=1.5,
        label=r"Ca$^{2+}$",
    )
    ax.plot(
        t,
        state["ROS"] / max(layer3_ros_max, 1e-12),
        color="#EE6677",
        linewidth=1.5,
        linestyle="--",
        label="ROS",
    )
    ax.plot(
        t,
        state["ATP"],
        color="#228833",
        linewidth=1.5,
        linestyle="-.",
        label="ATP",
    )
    ax.set_xlabel("Time (s)", fontsize=8, labelpad=2)
    ax.set_ylabel("Normalized state", fontsize=8, labelpad=2)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(0.0, 300.0)
    ax.set_xticks([0.0, 150.0, 300.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    _style_storyboard_axis(ax)


def _plot_layer4_panel(ax, result: SimulationResult) -> None:
    state = result.state_timeseries
    t = state["t"]
    ax.plot(
        t,
        state["PS_exposure"],
        color="#AA3377",
        linewidth=1.5,
        label="PS",
    )
    ax.plot(
        t,
        state["repair_state"],
        color="#228833",
        linewidth=1.5,
        linestyle="--",
        label="Repair",
    )
    ax.plot(
        t,
        state["actin_disruption"],
        color="#CCBB44",
        linewidth=1.5,
        linestyle="-.",
        label="Actin",
    )
    ax.set_xlabel("Time (s)", fontsize=8, labelpad=2)
    ax.set_ylabel("State", fontsize=8, labelpad=2)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(0.0, 300.0)
    ax.set_xticks([0.0, 150.0, 300.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    _style_storyboard_axis(ax)


def _plot_layer5_panel(ax, result: SimulationResult) -> None:
    sev = float(result.summary["cumulative_small_EV"])
    mlev = float(result.summary["cumulative_medium_large_EV"])
    ab = float(result.summary["cumulative_apoptotic_body"])
    labels = ["sEV", "m/lEV", "AB"]
    colors = ["#4477AA", "#228833", "#EE6677"]
    values = [sev, mlev, ab]
    ax.bar(labels, values, color=colors, width=0.72)
    ax.set_yscale("log")
    ax.set_ylabel("EV output", fontsize=8, labelpad=2)
    _style_storyboard_axis(ax)


def _style_storyboard_axis(ax) -> None:
    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        labelsize=7.5,
        length=2.5,
        width=0.7,
        pad=1.5,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=Path("results/multilayer_storyboard")
    )
    parser.add_argument(
        "--no-plots", action="store_true", help="Skip PNG plot generation."
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary, outputs = build_storyboard_dataset()
    write_outputs(summary, outputs, args.out, make_plots=not args.no_plots)
    print(f"Wrote multilayer storyboard outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
