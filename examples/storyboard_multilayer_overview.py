from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from electro_exocytosis.config import (
    CellStateConfig,
    ExposureConfig,
    PulseConfig,
    ScenarioConfig,
    SimulationConfig,
    SimulationScenario,
)
from electro_exocytosis.simulation import Simulation, SimulationResult
from electro_exocytosis.visualization.style import save_manuscript_figure


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
        result = Simulation(_build_scenario(spec), params_override=spec.parameter_overrides).run()
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
                "cumulative_medium_large_EV": result.summary["cumulative_medium_large_EV"],
                "cumulative_apoptotic_body": result.summary["cumulative_apoptotic_body"],
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
    summary.to_csv(outdir / "multilayer_storyboard_summary.csv", index=False)
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
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16.0, 8.8))
    if n_rows == 1:
        axes = np.array([axes])

    layer3_ca_max = max(float(result.state_timeseries["Ca_i"].max()) for _, result in outputs)
    layer3_ros_max = max(float(result.state_timeseries["ROS"].max()) for _, result in outputs)

    col_titles = [
        "Layer 1\nDosimetry",
        "Layer 2\nElectrodynamics",
        "Layer 3\nCa/ROS/ATP",
        "Layer 4\nRepair State",
        "Layer 5\nEV Output",
    ]
    for col, title in enumerate(col_titles):
        axes[0, col].set_title(title, fontsize=10, pad=10)

    for row, (spec, result) in enumerate(outputs):
        _plot_layer1_panel(axes[row, 0], result)
        _plot_layer2_panel(axes[row, 1], result)
        _plot_layer3_panel(axes[row, 2], result, layer3_ca_max=layer3_ca_max, layer3_ros_max=layer3_ros_max)
        _plot_layer4_panel(axes[row, 3], result)
        _plot_layer5_panel(axes[row, 4], result)
        axes[row, 0].text(
            -0.35,
            0.5,
            spec.label,
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=10,
            fontweight="bold",
        )

    fig.tight_layout(rect=(0.04, 0.02, 1.0, 1.0))
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
    t_train = max(float(descriptors["train_duration_s"]), float(descriptors["pulse_width_s"]))
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
    ax.set_ylabel("T (C)", fontsize=8)
    ax.set_xlabel("t (s)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.text(
        0.03,
        0.08,
        f"dose={result.summary['dose_index']:.2f}\nend DT={dosimetry['temperature_rise_K']:.2f} C",
        transform=ax.transAxes,
        fontsize=7,
        bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.8},
    )


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
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("DV (V)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.text(
        0.03,
        0.88,
        f"Pm={electro['membrane_permeability']:.2f}\nNp={electro['pore_density']/1e12:.2f} um^-2",
        transform=ax.transAxes,
        fontsize=7,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.8},
    )


def _plot_layer3_panel(ax, result: SimulationResult, *, layer3_ca_max: float, layer3_ros_max: float) -> None:
    state = result.state_timeseries
    t = state["t"]
    ax.plot(t, state["Ca_i"] / max(layer3_ca_max, 1e-12), color="#4477AA", linewidth=1.5, label="Ca_i")
    ax.plot(t, state["ROS"] / max(layer3_ros_max, 1e-12), color="#EE6677", linewidth=1.5, linestyle="--", label="ROS")
    ax.plot(t, state["ATP"], color="#228833", linewidth=1.5, linestyle="-.", label="ATP")
    ax.set_xlabel("t (s)", fontsize=8)
    ax.set_ylabel("Norm.", fontsize=8)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(0.0, 300.0)
    ax.tick_params(labelsize=7)
    if ax is ax.figure.axes[2]:
        ax.legend(fontsize=6, loc="upper right")


def _plot_layer4_panel(ax, result: SimulationResult) -> None:
    state = result.state_timeseries
    t = state["t"]
    ax.plot(t, state["PS_exposure"], color="#AA3377", linewidth=1.5, label="PS")
    ax.plot(t, state["repair_state"], color="#228833", linewidth=1.5, linestyle="--", label="repair")
    ax.plot(t, state["actin_disruption"], color="#CCBB44", linewidth=1.5, linestyle="-.", label="actin")
    ax.set_xlabel("t (s)", fontsize=8)
    ax.set_ylabel("State", fontsize=8)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(0.0, 300.0)
    ax.tick_params(labelsize=7)
    if ax is ax.figure.axes[3]:
        ax.legend(fontsize=6, loc="upper right")


def _plot_layer5_panel(ax, result: SimulationResult) -> None:
    sev = float(result.summary["cumulative_small_EV"])
    mlev = float(result.summary["cumulative_medium_large_EV"])
    ab = float(result.summary["cumulative_apoptotic_body"])
    total = sev + mlev + ab
    labels = ["sEV", "m/lEV", "AB"]
    colors = ["#4477AA", "#228833", "#EE6677"]
    values = [sev, mlev, ab]
    ax.bar(labels, values, color=colors, width=0.72)
    ax.set_yscale("log")
    ax.set_ylabel("Cumulative EV", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.tick_params(labelsize=7)
    ax.text(
        0.02,
        0.92,
        f"total={total:.2f}\nviab={result.summary['viability_fraction']:.2f}",
        transform=ax.transAxes,
        fontsize=7,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.8},
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/multilayer_storyboard"))
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary, outputs = build_storyboard_dataset()
    write_outputs(summary, outputs, args.out, make_plots=not args.no_plots)
    print(f"Wrote multilayer storyboard outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
