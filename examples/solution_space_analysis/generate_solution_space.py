from __future__ import annotations

import argparse
import json
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib"))

import numpy as np
import pandas as pd
import yaml

from electro_exocytosis.abbreviations import STANDARD_ABBREVIATIONS
from electro_exocytosis.config import SimulationScenario
from electro_exocytosis.simulation import Simulation, SimulationResult
from electro_exocytosis.visualization.style import save_manuscript_figure


DEFAULT_SCENARIO_DIR = Path(__file__).resolve().parent
DEFAULT_OUTDIR = Path("results/solution_space_analysis")

AMPLITUDES_KV_CM = (5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
PULSE_WIDTHS_NS = (50.0, 100.0, 200.0)
PULSE_NUMBERS = (5, 10, 20, 40, 60)
REPETITION_RATES_HZ = (1.0, 5.0, 20.0)
MEDIUM_CONDUCTIVITIES_S_M = (0.5, 1.0, 1.5, 2.0)
WAVEFORMS = ("square", "bipolar", "exponential")
DOSIMETRY_MODELS = ("legacy", "joule_adiabatic", "joule_lumped_thermal")
MODIFIER_LEVELS = (0.7, 1.0, 1.3)

BASE_PULSE = {
    "amplitude_kV_cm": 15.0,
    "pulse_width_ns": 100.0,
    "pulse_number": 20,
    "repetition_rate_Hz": 5.0,
    "waveform": "square",
}
BASE_EXPOSURE = {
    "geometry": "cuvette",
    "medium_conductivity_S_m": 1.5,
    "temperature_C": 37.0,
    "cell_density_per_ml": 1_000_000,
    "dosimetry_model": "joule_lumped_thermal",
}
BASE_CELL_STATE = {
    "cell_type": "generic",
    "membrane_modifier": 1.0,
    "calcium_handling_modifier": 1.0,
    "baseline_EV_release_modifier": 1.0,
    "stress_sensitivity_modifier": 1.0,
}
BASE_SIMULATION = {
    "t_start_s": 0.0,
    "t_end_s": 1800.0,
    "output_dt_s": 20.0,
    "numerical_method": "solve_ivp",
}


@dataclass(frozen=True)
class SolutionSpaceScenarioSpec:
    name: str
    label: str
    family: str
    pulse_overrides: dict[str, Any] = field(default_factory=dict)
    exposure_overrides: dict[str, Any] = field(default_factory=dict)
    cell_state_overrides: dict[str, Any] = field(default_factory=dict)
    sweep_axis: str = "dose_grid"
    sweep_value: str = "nominal"


def build_solution_space_specs(
    amplitudes_kV_cm: Sequence[float] = AMPLITUDES_KV_CM,
    pulse_widths_ns: Sequence[float] = PULSE_WIDTHS_NS,
    pulse_numbers: Sequence[int] = PULSE_NUMBERS,
) -> tuple[SolutionSpaceScenarioSpec, ...]:
    """Build a broad but bounded set of scenario definitions for the current model."""
    specs: list[SolutionSpaceScenarioSpec] = []
    specs.extend(_dose_grid_specs(amplitudes_kV_cm, pulse_widths_ns, pulse_numbers))
    specs.extend(_waveform_conductivity_specs())
    specs.extend(_repetition_rate_specs())
    specs.extend(_dosimetry_model_specs())
    specs.extend(_cell_modifier_specs())
    return tuple(specs)


def _dose_grid_specs(
    amplitudes_kV_cm: Sequence[float],
    pulse_widths_ns: Sequence[float],
    pulse_numbers: Sequence[int],
) -> list[SolutionSpaceScenarioSpec]:
    specs: list[SolutionSpaceScenarioSpec] = []
    for pulse_width_ns in pulse_widths_ns:
        for amplitude_kV_cm in amplitudes_kV_cm:
            for pulse_number in pulse_numbers:
                name = (
                    f"dose_grid_w{_int_token(pulse_width_ns, 3)}ns_"
                    f"a{_int_token(amplitude_kV_cm, 2)}kvcm_n{_int_token(pulse_number, 3)}"
                )
                specs.append(
                    SolutionSpaceScenarioSpec(
                        name=name,
                        label=f"{amplitude_kV_cm:g} kV/cm, {pulse_width_ns:g} ns, {pulse_number} pulses",
                        family="dose_grid",
                        pulse_overrides={
                            "amplitude_kV_cm": float(amplitude_kV_cm),
                            "pulse_width_ns": float(pulse_width_ns),
                            "pulse_number": int(pulse_number),
                            "repetition_rate_Hz": 5.0,
                        },
                        sweep_axis="amplitude_width_pulse_number",
                        sweep_value=f"{amplitude_kV_cm:g}_{pulse_width_ns:g}_{pulse_number}",
                    )
                )
    return specs


def _waveform_conductivity_specs() -> list[SolutionSpaceScenarioSpec]:
    specs: list[SolutionSpaceScenarioSpec] = []
    for waveform in WAVEFORMS:
        for conductivity in MEDIUM_CONDUCTIVITIES_S_M:
            name = f"waveform_{waveform}_sigma{_decimal_token(conductivity)}"
            specs.append(
                SolutionSpaceScenarioSpec(
                    name=name,
                    label=f"{waveform} waveform, conductivity {conductivity:g} S/m",
                    family="waveform_conductivity",
                    pulse_overrides={"waveform": waveform},
                    exposure_overrides={"medium_conductivity_S_m": float(conductivity)},
                    sweep_axis="waveform_conductivity",
                    sweep_value=f"{waveform}_{conductivity:g}",
                )
            )
    return specs


def _repetition_rate_specs() -> list[SolutionSpaceScenarioSpec]:
    specs: list[SolutionSpaceScenarioSpec] = []
    for repetition_rate in REPETITION_RATES_HZ:
        name = f"repetition_rate_{_decimal_token(repetition_rate)}hz"
        specs.append(
            SolutionSpaceScenarioSpec(
                name=name,
                label=f"{repetition_rate:g} Hz pulse repetition rate",
                family="pulse_timing",
                pulse_overrides={"repetition_rate_Hz": float(repetition_rate)},
                sweep_axis="repetition_rate_Hz",
                sweep_value=f"{repetition_rate:g}",
            )
        )
    return specs


def _dosimetry_model_specs() -> list[SolutionSpaceScenarioSpec]:
    specs: list[SolutionSpaceScenarioSpec] = []
    for dosimetry_model in DOSIMETRY_MODELS:
        specs.append(
            SolutionSpaceScenarioSpec(
                name=f"dosimetry_{dosimetry_model}",
                label=f"{dosimetry_model.replace('_', ' ')} dosimetry",
                family="dosimetry_model",
                exposure_overrides={"dosimetry_model": dosimetry_model},
                sweep_axis="dosimetry_model",
                sweep_value=dosimetry_model,
            )
        )
    return specs


def _cell_modifier_specs() -> list[SolutionSpaceScenarioSpec]:
    modifier_labels = {
        "membrane_modifier": "membrane sensitivity",
        "calcium_handling_modifier": "calcium handling",
        "baseline_EV_release_modifier": "baseline EV release",
        "stress_sensitivity_modifier": "stress sensitivity",
    }
    specs: list[SolutionSpaceScenarioSpec] = [
        SolutionSpaceScenarioSpec(
            name="cell_state_nominal",
            label="Nominal cell state",
            family="cell_state_modifier",
            sweep_axis="cell_state_baseline",
            sweep_value="1.0",
        )
    ]
    for modifier_name, modifier_label in modifier_labels.items():
        for value in MODIFIER_LEVELS:
            if value == 1.0:
                continue
            level = "low" if value < 1.0 else "high"
            specs.append(
                SolutionSpaceScenarioSpec(
                    name=f"cell_state_{modifier_name.replace('_modifier', '')}_{level}",
                    label=f"{modifier_label}: {value:g}x",
                    family="cell_state_modifier",
                    cell_state_overrides={modifier_name: float(value)},
                    sweep_axis=modifier_name,
                    sweep_value=f"{value:g}",
                )
            )
    combined_states = (
        (
            "cell_state_resilient_secretory",
            "Resilient secretory-biased cell state",
            {
                "membrane_modifier": 0.8,
                "calcium_handling_modifier": 1.2,
                "baseline_EV_release_modifier": 1.3,
                "stress_sensitivity_modifier": 0.7,
            },
        ),
        (
            "cell_state_fragile_stress_biased",
            "Fragile stress-biased cell state",
            {
                "membrane_modifier": 1.3,
                "calcium_handling_modifier": 0.8,
                "baseline_EV_release_modifier": 0.8,
                "stress_sensitivity_modifier": 1.3,
            },
        ),
    )
    for name, label, overrides in combined_states:
        specs.append(
            SolutionSpaceScenarioSpec(
                name=name,
                label=label,
                family="cell_state_modifier",
                cell_state_overrides=overrides,
                sweep_axis="combined_cell_state",
                sweep_value=name.replace("cell_state_", ""),
            )
        )
    return specs


def scenario_payload(spec: SolutionSpaceScenarioSpec) -> dict[str, Any]:
    pulse = {**BASE_PULSE, **spec.pulse_overrides}
    exposure = {**BASE_EXPOSURE, **spec.exposure_overrides}
    cell_state = {**BASE_CELL_STATE, **spec.cell_state_overrides}
    return {
        "scenario": {
            "name": spec.name,
            "mode": "cell_based_electro_exocytosis",
        },
        "pulse": pulse,
        "exposure": exposure,
        "cell_state": cell_state,
        "simulation": dict(BASE_SIMULATION),
    }


def write_scenario_files(
    specs: Sequence[SolutionSpaceScenarioSpec],
    scenario_dir: Path = DEFAULT_SCENARIO_DIR,
) -> list[Path]:
    scenario_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for spec in specs:
        path = scenario_dir / f"scenario_{spec.name}.yml"
        with path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(scenario_payload(spec), handle, sort_keys=False)
        written.append(path)
    return written


def build_solution_space_dataset(
    specs: Sequence[SolutionSpaceScenarioSpec],
) -> tuple[pd.DataFrame, dict[str, SimulationResult]]:
    rows: list[dict[str, Any]] = []
    results: dict[str, SimulationResult] = {}
    for spec in specs:
        scenario = SimulationScenario.model_validate(scenario_payload(spec))
        result = Simulation(scenario).run()
        results[spec.name] = result
        rows.append(_summary_row(spec, result))
    return pd.DataFrame(rows), results


def write_outputs(
    summary: pd.DataFrame,
    results: dict[str, SimulationResult],
    outdir: Path = DEFAULT_OUTDIR,
    make_plots: bool = True,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary = summary.sort_values(["scenario_family", "scenario_name"]).reset_index(drop=True)
    summary.to_csv(outdir / "solution_space_summary_raw.csv", index=False)
    STANDARD_ABBREVIATIONS.rename_columns(summary).to_csv(outdir / "solution_space_summary.csv", index=False)
    _write_rankings(summary, outdir)
    _write_manifest(summary, outdir)
    STANDARD_ABBREVIATIONS.write_bundle(
        outdir,
        keys=("nsPEF", "EV", "ER", "MVB", "ROS", "ATP", "PS", "sEV", "m/lEV", "AB"),
    )

    if make_plots:
        plot_dose_grid_heatmaps(summary, outdir / "solution_space_dose_heatmaps.png")
        plot_tradeoff_space(summary, outdir / "solution_space_tradeoffs.png")
        plot_modifier_sensitivity(summary, outdir / "solution_space_cell_state_sensitivity.png")
        plot_representative_timeseries(summary, results, outdir / "solution_space_representative_timeseries.png")


def _summary_row(spec: SolutionSpaceScenarioSpec, result: SimulationResult) -> dict[str, Any]:
    scenario = result.parameters_used["scenario"]
    pulse = scenario["pulse"]
    exposure = scenario["exposure"]
    cell_state = scenario["cell_state"]
    dosimetry = result.parameters_used["dosimetry"]
    electrodynamics = result.parameters_used["electrodynamics"]
    summary = result.summary
    cumulative_small_ev = float(summary["cumulative_small_EV"])
    cumulative_medium_large_ev = float(summary["cumulative_medium_large_EV"])
    cumulative_apoptotic_body = float(summary["cumulative_apoptotic_body"])
    total_ev = cumulative_small_ev + cumulative_medium_large_ev + cumulative_apoptotic_body
    total_ev_safe = max(total_ev, 1e-12)
    quality_pass = bool(result.parameters_used["terminal_quality"]["quality_pass"])
    return {
        "scenario_name": spec.name,
        "scenario_label": spec.label,
        "scenario_family": spec.family,
        "sweep_axis": spec.sweep_axis,
        "sweep_value": spec.sweep_value,
        "amplitude_kV_cm": float(pulse["amplitude_kV_cm"]),
        "pulse_width_ns": float(pulse["pulse_width_ns"]),
        "pulse_number": int(pulse["pulse_number"]),
        "repetition_rate_Hz": float(pulse["repetition_rate_Hz"]),
        "waveform": str(pulse["waveform"]),
        "medium_conductivity_S_m": float(exposure["medium_conductivity_S_m"]),
        "dosimetry_model": str(exposure["dosimetry_model"]),
        "temperature_C": float(exposure["temperature_C"]),
        "membrane_modifier": float(cell_state["membrane_modifier"]),
        "calcium_handling_modifier": float(cell_state["calcium_handling_modifier"]),
        "baseline_EV_release_modifier": float(cell_state["baseline_EV_release_modifier"]),
        "stress_sensitivity_modifier": float(cell_state["stress_sensitivity_modifier"]),
        "dose_index": float(summary["dose_index"]),
        "temperature_rise_C": float(dosimetry["temperature_rise_K"]),
        "adiabatic_temperature_rise_C": float(dosimetry["adiabatic_temperature_rise_K"]),
        "thermal_retention_factor": float(dosimetry["thermal_retention_factor"]),
        "plasma_membrane_voltage_V": float(electrodynamics["delta_Vm"]),
        "endoplasmic_reticulum_voltage_V": float(electrodynamics["delta_V_ER"]),
        "mitochondrial_voltage_V": float(electrodynamics["delta_V_mito"]),
        "multivesicular_body_voltage_V": float(electrodynamics["delta_V_MVB"]),
        "membrane_permeability": float(electrodynamics["membrane_permeability"]),
        "pore_density_per_m2": float(electrodynamics["pore_density"]),
        "peak_cytosolic_calcium_uM": float(summary["peak_ca_i"]),
        "peak_reactive_oxygen_species": float(summary["peak_ros"]),
        "minimum_adenosine_triphosphate_state": float(summary["min_atp"]),
        "peak_phosphatidylserine_exposure": float(summary["peak_ps_exposure"]),
        "peak_resealing_state": float(summary["peak_repair_state"]),
        "peak_secretory_routing_bias": float(summary["peak_secretory_bias"]),
        "cumulative_small_extracellular_vesicle_output": cumulative_small_ev,
        "cumulative_medium_large_extracellular_vesicle_output": cumulative_medium_large_ev,
        "cumulative_apoptotic_body_output": cumulative_apoptotic_body,
        "total_extracellular_vesicle_output": total_ev,
        "small_extracellular_vesicle_fraction": cumulative_small_ev / total_ev_safe,
        "medium_large_extracellular_vesicle_fraction": cumulative_medium_large_ev / total_ev_safe,
        "apoptotic_body_fraction": cumulative_apoptotic_body / total_ev_safe,
        "viability_fraction": float(summary["viability_fraction"]),
        "purity_score": float(summary["purity_score"]),
        "potency_score": float(summary["potency_score"]),
        "cell_normalized_yield": float(summary["cell_normalized_yield"]),
        "optimization_objective": float(summary["optimization_objective"]),
        "quality_pass": quality_pass,
    }


def plot_dose_grid_heatmaps(summary: pd.DataFrame, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dose_grid = summary[summary["scenario_family"] == "dose_grid"].copy()
    pulse_widths = sorted(dose_grid["pulse_width_ns"].unique())
    pulse_numbers = sorted(dose_grid["pulse_number"].unique())
    amplitudes = sorted(dose_grid["amplitude_kV_cm"].unique())
    metrics = (
        ("dose_index", "log10 dose index", True),
        ("plasma_membrane_voltage_V", "PM voltage (V)", False),
        ("peak_cytosolic_calcium_uM", "peak cytosolic calcium (uM)", False),
        ("total_extracellular_vesicle_output", "log10 total EV output", True),
        ("viability_fraction", "viability fraction", False),
        ("optimization_objective", "optimization objective", False),
    )
    fig, axes = plt.subplots(len(pulse_widths), len(metrics), figsize=(20.0, 9.6), squeeze=False)
    for row_index, pulse_width in enumerate(pulse_widths):
        width_frame = dose_grid[dose_grid["pulse_width_ns"] == pulse_width]
        for col_index, (metric, title, log_scale) in enumerate(metrics):
            ax = axes[row_index, col_index]
            grid = _pivot_grid(width_frame, metric, amplitudes, pulse_numbers)
            plotted_grid = np.log10(np.clip(grid, 1e-12, None)) if log_scale else grid
            im = ax.imshow(plotted_grid, origin="lower", aspect="auto", cmap="viridis")
            if row_index == 0:
                ax.set_title(title, fontsize=10)
            if col_index == 0:
                ax.set_ylabel(f"{pulse_width:g} ns\nAmplitude (kV/cm)", fontsize=9)
            else:
                ax.set_ylabel("")
            ax.set_xlabel("Pulse number", fontsize=8)
            ax.set_xticks(range(len(pulse_numbers)))
            ax.set_xticklabels([str(value) for value in pulse_numbers], fontsize=7)
            ax.set_yticks(range(len(amplitudes)))
            ax.set_yticklabels([f"{value:g}" for value in amplitudes], fontsize=7)
            cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
            cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    save_manuscript_figure(fig, output_path, abbreviation_keys=("PM", "EV"))
    plt.close(fig)


def plot_tradeoff_space(summary: pd.DataFrame, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.2))
    family_markers = {
        "dose_grid": "o",
        "waveform_conductivity": "s",
        "pulse_timing": "^",
        "dosimetry_model": "D",
        "cell_state_modifier": "P",
    }
    objective_norm = plt.Normalize(summary["optimization_objective"].min(), summary["optimization_objective"].max())
    for family, marker in family_markers.items():
        frame = summary[summary["scenario_family"] == family]
        if frame.empty:
            continue
        axes[0].scatter(
            frame["total_extracellular_vesicle_output"],
            frame["viability_fraction"],
            c=frame["optimization_objective"],
            cmap="viridis",
            norm=objective_norm,
            marker=marker,
            edgecolor="0.2",
            linewidth=0.3,
            alpha=0.85,
            label=family.replace("_", " "),
        )
    axes[0].set_xlabel("Total EV output")
    axes[0].set_ylabel("Viability fraction")
    axes[0].legend(fontsize=7)
    axes[0].set_title("Yield versus viability")
    objective_mappable = plt.cm.ScalarMappable(norm=objective_norm, cmap="viridis")
    objective_mappable.set_array([])
    cbar = fig.colorbar(objective_mappable, ax=axes[0], fraction=0.046, pad=0.04)
    cbar.set_label("Optimization objective", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    dose_grid = summary[summary["scenario_family"] == "dose_grid"]
    sc = axes[1].scatter(
        dose_grid["small_extracellular_vesicle_fraction"],
        dose_grid["apoptotic_body_fraction"],
        c=dose_grid["dose_index"],
        s=30 + 80 * np.clip(dose_grid["potency_score"], 0, 1),
        cmap="plasma",
        edgecolor="0.2",
        linewidth=0.25,
        alpha=0.82,
    )
    axes[1].set_xlabel("Small EV fraction")
    axes[1].set_ylabel("Apoptotic-body fraction")
    axes[1].set_title("Subtype balance across dose grid")
    cbar = fig.colorbar(sc, ax=axes[1], fraction=0.046, pad=0.04)
    cbar.set_label("Dose index", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    axes[2].scatter(
        summary["potency_score"],
        summary["purity_score"],
        c=summary["viability_fraction"],
        s=30 + 80 * np.clip(summary["apoptotic_body_fraction"], 0, 1),
        cmap="cividis",
        edgecolor="0.2",
        linewidth=0.25,
        alpha=0.82,
    )
    axes[2].set_xlabel("Potency score")
    axes[2].set_ylabel("Purity score")
    axes[2].set_title("Potency, purity, and injury pressure")
    cbar = fig.colorbar(axes[2].collections[0], ax=axes[2], fraction=0.046, pad=0.04)
    cbar.set_label("Viability fraction", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    for ax in axes:
        ax.grid(True, color="0.9", linewidth=0.7)
    fig.tight_layout()
    save_manuscript_figure(fig, output_path, abbreviation_keys=("EV", "sEV", "AB"))
    plt.close(fig)


def plot_modifier_sensitivity(summary: pd.DataFrame, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modifier_frame = summary[summary["scenario_family"] == "cell_state_modifier"].copy()
    baseline = modifier_frame[modifier_frame["scenario_name"] == "cell_state_nominal"].iloc[0]
    metrics = (
        ("cumulative_small_extracellular_vesicle_output", "small EV output"),
        ("viability_fraction", "viability"),
        ("apoptotic_body_fraction", "apoptotic-body fraction"),
        ("optimization_objective", "objective"),
    )
    axes_to_show = (
        "membrane_modifier",
        "calcium_handling_modifier",
        "baseline_EV_release_modifier",
        "stress_sensitivity_modifier",
    )
    sensitivity = np.zeros((len(axes_to_show), len(metrics)))
    annotations: list[list[str]] = []
    for row_index, axis_name in enumerate(axes_to_show):
        axis_frame = modifier_frame[modifier_frame["sweep_axis"] == axis_name]
        low = axis_frame[axis_frame["sweep_value"].astype(float) < 1.0].iloc[0]
        high = axis_frame[axis_frame["sweep_value"].astype(float) > 1.0].iloc[0]
        row_annotations: list[str] = []
        for col_index, (metric, _) in enumerate(metrics):
            base = max(float(baseline[metric]), 1e-12)
            relative_low = float(low[metric]) / base
            relative_high = float(high[metric]) / base
            value = relative_high - relative_low
            sensitivity[row_index, col_index] = value
            row_annotations.append(f"{value:+.2f}")
        annotations.append(row_annotations)

    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    im = ax.imshow(sensitivity, cmap="coolwarm", aspect="auto", vmin=-np.max(np.abs(sensitivity)), vmax=np.max(np.abs(sensitivity)))
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([title for _, title in metrics], rotation=25, ha="right", fontsize=8)
    ax.set_yticks(range(len(axes_to_show)))
    ax.set_yticklabels([axis_name.replace("_modifier", "").replace("_", " ") for axis_name in axes_to_show], fontsize=8)
    ax.set_title("One-at-a-time cell-state sensitivity: high modifier response minus low modifier response")
    for row_index, row in enumerate(annotations):
        for col_index, label in enumerate(row):
            ax.text(col_index, row_index, label, ha="center", va="center", fontsize=8, color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label("Change in metric relative to nominal", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    save_manuscript_figure(fig, output_path, abbreviation_keys=("EV", "sEV", "AB"))
    plt.close(fig)


def plot_representative_timeseries(
    summary: pd.DataFrame,
    results: dict[str, SimulationResult],
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = _select_representative_scenarios(summary)
    fig, axes = plt.subplots(len(selected), 4, figsize=(16.0, 8.6), squeeze=False)
    for row_index, (regime_label, scenario_name) in enumerate(selected):
        result = results[scenario_name]
        state = result.state_timeseries
        ev = result.ev_timeseries
        t_min = state["t"] / 60.0
        ev_t_min = ev["t"] / 60.0

        axes[row_index, 0].plot(t_min, state["Ca_i"], color="#4477AA", linewidth=1.4, label="cytosolic calcium")
        axes[row_index, 0].plot(t_min, state["Ca_submembrane"], color="#66CCEE", linewidth=1.4, linestyle="--", label="submembrane calcium")
        axes[row_index, 0].set_ylabel(f"{regime_label}\nCalcium (uM)", fontsize=8)

        axes[row_index, 1].plot(t_min, state["ROS"], color="#EE6677", linewidth=1.4, label="reactive oxygen species")
        axes[row_index, 1].plot(t_min, state["ATP"], color="#228833", linewidth=1.4, linestyle="--", label="adenosine triphosphate")
        axes[row_index, 1].plot(t_min, state["damage"], color="#111111", linewidth=1.2, linestyle=":", label="damage")
        axes[row_index, 1].set_ylabel("Stress state", fontsize=8)

        axes[row_index, 2].plot(t_min, state["PS_exposure"], color="#AA3377", linewidth=1.4, label="phosphatidylserine")
        axes[row_index, 2].plot(t_min, state["repair_state"], color="#228833", linewidth=1.4, linestyle="--", label="repair")
        axes[row_index, 2].plot(t_min, state["actin_disruption"], color="#CCBB44", linewidth=1.4, linestyle=":", label="actin disruption")
        axes[row_index, 2].set_ylabel("Remodeling state", fontsize=8)

        axes[row_index, 3].plot(ev_t_min, ev["sEV_cumulative"], color="#4477AA", linewidth=1.4, label="small EV")
        axes[row_index, 3].plot(ev_t_min, ev["mlEV_cumulative"], color="#228833", linewidth=1.4, linestyle="--", label="medium/large EV")
        axes[row_index, 3].plot(ev_t_min, ev["AB_cumulative"], color="#EE6677", linewidth=1.4, linestyle=":", label="apoptotic body")
        axes[row_index, 3].set_yscale("log")
        axes[row_index, 3].set_ylabel("Cumulative output", fontsize=8)

        for col_index in range(4):
            axes[row_index, col_index].grid(True, color="0.9", linewidth=0.7)
            axes[row_index, col_index].tick_params(labelsize=7)
            if row_index == len(selected) - 1:
                axes[row_index, col_index].set_xlabel("Time (min)", fontsize=8)
            if row_index == 0:
                axes[row_index, col_index].legend(fontsize=6)

    for col_index, title in enumerate(("Calcium", "Stress and energy", "Repair/remodeling", "EV subtype output")):
        axes[0, col_index].set_title(title, fontsize=10)

    fig.tight_layout()
    save_manuscript_figure(fig, output_path, abbreviation_keys=("EV", "ROS", "ATP", "PS", "sEV", "m/lEV", "AB"))
    plt.close(fig)


def _write_rankings(summary: pd.DataFrame, outdir: Path) -> None:
    ranked_objective = summary.sort_values("optimization_objective", ascending=False).head(12)
    viable = summary[summary["viability_fraction"] >= 0.75]
    viable_ranked = viable.sort_values(
        ["cumulative_small_extracellular_vesicle_output", "optimization_objective"],
        ascending=False,
    ).head(12)
    injury_ranked = summary.sort_values("apoptotic_body_fraction", ascending=False).head(12)
    for frame, name in (
        (ranked_objective, "top_scenarios_by_objective.csv"),
        (viable_ranked, "top_viable_small_ev_scenarios.csv"),
        (injury_ranked, "injury_skewed_scenarios.csv"),
    ):
        STANDARD_ABBREVIATIONS.rename_columns(frame).to_csv(outdir / name, index=False)


def _write_manifest(summary: pd.DataFrame, outdir: Path) -> None:
    payload = {
        "analysis": "solution_space_analysis",
        "scenario_count": int(len(summary)),
        "scenario_families": {
            str(family): int(count) for family, count in summary["scenario_family"].value_counts().sort_index().items()
        },
        "parameter_ranges": {
            "amplitude_kV_cm": [float(summary["amplitude_kV_cm"].min()), float(summary["amplitude_kV_cm"].max())],
            "pulse_width_ns": [float(summary["pulse_width_ns"].min()), float(summary["pulse_width_ns"].max())],
            "pulse_number": [int(summary["pulse_number"].min()), int(summary["pulse_number"].max())],
            "repetition_rate_Hz": [float(summary["repetition_rate_Hz"].min()), float(summary["repetition_rate_Hz"].max())],
            "medium_conductivity_S_m": [
                float(summary["medium_conductivity_S_m"].min()),
                float(summary["medium_conductivity_S_m"].max()),
            ],
        },
        "model_status": "Exploratory reduced equations; not experimentally validated.",
    }
    (outdir / "analysis_manifest.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _select_representative_scenarios(summary: pd.DataFrame) -> list[tuple[str, str]]:
    dose_grid = summary[summary["scenario_family"] == "dose_grid"].copy()
    mild = dose_grid.sort_values(["dose_index", "total_extracellular_vesicle_output"]).iloc[0]
    viable = dose_grid[(dose_grid["viability_fraction"] >= 0.75) & (dose_grid["apoptotic_body_fraction"] < 0.35)]
    if viable.empty:
        viable = dose_grid
    productive = viable.sort_values("optimization_objective", ascending=False).iloc[0]
    injury = dose_grid.sort_values(["apoptotic_body_fraction", "total_extracellular_vesicle_output"], ascending=False).iloc[0]
    return [
        ("Mild low-dose", str(mild["scenario_name"])),
        ("Productive viable", str(productive["scenario_name"])),
        ("Injury-skewed", str(injury["scenario_name"])),
    ]


def _pivot_grid(
    frame: pd.DataFrame,
    metric: str,
    amplitudes: Sequence[float],
    pulse_numbers: Sequence[int],
) -> np.ndarray:
    table = frame.pivot(index="amplitude_kV_cm", columns="pulse_number", values=metric)
    table = table.reindex(index=amplitudes, columns=pulse_numbers)
    return table.to_numpy(dtype=float)


def _int_token(value: float | int, width: int) -> str:
    return f"{int(round(float(value))):0{width}d}"


def _decimal_token(value: float | int) -> str:
    return f"{float(value):g}".replace(".", "p")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario-dir", type=Path, default=DEFAULT_SCENARIO_DIR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG figure generation.")
    parser.add_argument("--no-scenarios", action="store_true", help="Do not write scenario_*.yml files.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    specs = build_solution_space_specs()
    if not args.no_scenarios:
        written = write_scenario_files(specs, args.scenario_dir)
    else:
        written = []
    summary, results = build_solution_space_dataset(specs)
    write_outputs(summary, results, args.out, make_plots=not args.no_plots)
    if written:
        print(f"Wrote {len(written)} scenario files to {args.scenario_dir}")
    print(f"Wrote solution-space outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
