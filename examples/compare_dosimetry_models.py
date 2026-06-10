from __future__ import annotations

import argparse
import os
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from electro_exocytosis.config import ExposureConfig, PulseConfig
from electro_exocytosis.models.dosimetry import DosimetryResult, compute_dosimetry
from electro_exocytosis.models.pulse import PulseDescriptors, compute_pulse_descriptors
from electro_exocytosis.visualization.style import (
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE,
    bar_colors,
    bar_hatch,
    line_styles,
    save_manuscript_figure,
)

DOSIMETRY_MODELS = ("legacy", "joule_adiabatic", "joule_lumped_thermal")

MODEL_NOTES = {
    "legacy": "Original absorbed-energy equation; linear geometry scaling; no waveform or cooling correction.",
    "joule_adiabatic": "Waveform-aware Joule energy; geometry heat scales with mean E^2; no cooling correction.",
    "joule_lumped_thermal": "Waveform-aware Joule energy plus finite-train heat retention and post-train cooling.",
}


@dataclass(frozen=True)
class ScenarioSpec:
    name: str
    label: str
    pulse_kwargs: dict[str, float | int | str]
    exposure_kwargs: dict[str, float | str]


def scenario_specs() -> list[ScenarioSpec]:
    """Return nsPEF cases chosen to expose different model sensitivities."""
    empirical_in_vitro_thermal_params = {
        "thermal_relaxation_time_s": 5.0,
        "thermal_efficiency": 0.35,
    }
    return [
        ScenarioSpec(
            name="zap_like_low_conductivity_20Hz",
            label="Low conductivity buffer, 20 Hz",
            pulse_kwargs={
                "amplitude_kV_cm": 90,
                "pulse_width_ns": 10,
                "pulse_number": 500,
                "repetition_rate_Hz": 20,
                "waveform": "square",
            },
            exposure_kwargs={
                "geometry": "cuvette",
                "medium_conductivity_S_m": 0.2,
                "temperature_C": 22.0,
                **empirical_in_vitro_thermal_params,
            },
        ),
        ScenarioSpec(
            name="hbss_like_high_conductivity_20Hz",
            label="High conductivity buffer, 20 Hz",
            pulse_kwargs={
                "amplitude_kV_cm": 100,
                "pulse_width_ns": 10,
                "pulse_number": 500,
                "repetition_rate_Hz": 20,
                "waveform": "square",
            },
            exposure_kwargs={
                "geometry": "cuvette",
                "medium_conductivity_S_m": 1.4,
                "temperature_C": 22.0,
                **empirical_in_vitro_thermal_params,
            },
        ),
        ScenarioSpec(
            name="hbss_like_high_conductivity_200Hz",
            label="High conductivity buffer, 200 Hz",
            pulse_kwargs={
                "amplitude_kV_cm": 100,
                "pulse_width_ns": 10,
                "pulse_number": 500,
                "repetition_rate_Hz": 200,
                "waveform": "square",
            },
            exposure_kwargs={
                "geometry": "cuvette",
                "medium_conductivity_S_m": 1.4,
                "temperature_C": 22.0,
                **empirical_in_vitro_thermal_params,
            },
        ),
        ScenarioSpec(
            name="dish_exponential_high_prr",
            label="Dish geometry, exponential pulse, 200 Hz",
            pulse_kwargs={
                "amplitude_kV_cm": 100,
                "pulse_width_ns": 10,
                "pulse_number": 500,
                "repetition_rate_Hz": 200,
                "waveform": "exponential",
            },
            exposure_kwargs={
                "geometry": "dish",
                "medium_conductivity_S_m": 1.4,
                "temperature_C": 22.0,
                **empirical_in_vitro_thermal_params,
            },
        ),
    ]


def build_comparison_tables(
    scenarios: Sequence[ScenarioSpec] | None = None,
    profile_points: int = 240,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute model comparison summary and inferred temperature profiles."""
    summary_rows: list[dict[str, float | int | str]] = []
    profile_rows: list[dict[str, float | str]] = []
    for scenario in scenarios or scenario_specs():
        pulse = PulseConfig(**scenario.pulse_kwargs)
        for model in DOSIMETRY_MODELS:
            exposure = ExposureConfig(**scenario.exposure_kwargs, dosimetry_model=model)
            descriptors = compute_pulse_descriptors(pulse, exposure)
            dosimetry = compute_dosimetry(descriptors, exposure)
            summary_rows.append(_summary_row(scenario, pulse, exposure, descriptors, dosimetry))
            profile_rows.extend(
                _temperature_profile_rows(
                    scenario,
                    exposure,
                    descriptors,
                    dosimetry,
                    profile_points=profile_points,
                )
            )
    return pd.DataFrame(summary_rows), pd.DataFrame(profile_rows)


def _summary_row(
    scenario: ScenarioSpec,
    pulse: PulseConfig,
    exposure: ExposureConfig,
    descriptors: PulseDescriptors,
    dosimetry: DosimetryResult,
) -> dict[str, float | int | str]:
    return {
        "scenario": scenario.name,
        "scenario_label": scenario.label,
        "dosimetry_model": exposure.dosimetry_model,
        "model_note": MODEL_NOTES[exposure.dosimetry_model],
        "waveform": pulse.waveform,
        "geometry": exposure.geometry,
        "E_peak_kV_cm": pulse.amplitude_kV_cm,
        "pulse_width_ns": pulse.pulse_width_ns,
        "pulse_number": pulse.pulse_number,
        "repetition_rate_Hz": pulse.repetition_rate_Hz,
        "train_duration_s": descriptors.train_duration_s,
        "medium_conductivity_S_m": exposure.medium_conductivity_S_m,
        "waveform_energy_factor": descriptors.waveform_energy_factor,
        "field_uniformity_factor": dosimetry.field_uniformity_factor,
        "mean_E2_factor": dosimetry.mean_E2_factor,
        "absorbed_energy_density_mJ_mm3": descriptors.energy_density_J_m3 / 1.0e6,
        "geometry_heat_density_mJ_mm3": dosimetry.joule_heat_density_J_m3 / 1.0e6,
        "adiabatic_temp_rise_C": dosimetry.adiabatic_temperature_rise_K,
        "end_temp_rise_C": dosimetry.temperature_rise_K,
        "end_temperature_C": exposure.temperature_C + dosimetry.temperature_rise_K,
        "thermal_retention_factor": dosimetry.thermal_retention_factor,
    }


def _temperature_profile_rows(
    scenario: ScenarioSpec,
    exposure: ExposureConfig,
    descriptors: PulseDescriptors,
    dosimetry: DosimetryResult,
    profile_points: int,
) -> list[dict[str, float | str]]:
    train_duration_s = max(descriptors.train_duration_s, descriptors.pulse_width_s)
    horizon_s = train_duration_s + max(train_duration_s, 10.0)
    times = np.linspace(0.0, horizon_s, profile_points)
    temp_rise = _temperature_rise_profile(times, train_duration_s, exposure, dosimetry)
    return [
        {
            "scenario": scenario.name,
            "scenario_label": scenario.label,
            "dosimetry_model": exposure.dosimetry_model,
            "time_s": float(time_s),
            "temp_rise_C": float(rise_C),
            "temperature_C": float(exposure.temperature_C + rise_C),
            "train_duration_s": train_duration_s,
        }
        for time_s, rise_C in zip(times, temp_rise, strict=True)
    ]


def _temperature_rise_profile(
    times_s: np.ndarray,
    train_duration_s: float,
    exposure: ExposureConfig,
    dosimetry: DosimetryResult,
) -> np.ndarray:
    """Infer a simple temperature trace from each endpoint dosimetry model."""
    heating_times_s = np.minimum(times_s, train_duration_s)
    if exposure.dosimetry_model == "joule_lumped_thermal":
        tau_s = exposure.thermal_relaxation_time_s
        heating_rate_C_s = dosimetry.adiabatic_temperature_rise_K / train_duration_s
        rise_C = exposure.thermal_efficiency * heating_rate_C_s * tau_s * (
            1.0 - np.exp(-heating_times_s / tau_s)
        )
        cooling_mask = times_s > train_duration_s
        rise_C[cooling_mask] = rise_C[cooling_mask] * np.exp(
            -(times_s[cooling_mask] - train_duration_s) / tau_s
        )
        return rise_C

    return dosimetry.adiabatic_temperature_rise_K * (heating_times_s / train_duration_s)


def plot_temperature_profiles(profiles: pd.DataFrame, outdir: Path) -> None:
    """Plot inferred temperature trajectories for each scenario and model."""
    outdir.mkdir(parents=True, exist_ok=True)
    scenario_labels = profiles["scenario_label"].drop_duplicates().tolist()
    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE, sharey=False)
    for ax, scenario_label in zip(axes.flat, scenario_labels, strict=True):
        subset = profiles[profiles["scenario_label"] == scenario_label]
        train_duration_s = float(subset["train_duration_s"].iloc[0])
        styles = line_styles(len(DOSIMETRY_MODELS))
        for index, model in enumerate(DOSIMETRY_MODELS):
            model_subset = subset[subset["dosimetry_model"] == model]
            ax.plot(model_subset["time_s"], model_subset["temp_rise_C"], label=model, **styles[index])
        ax.axvline(train_duration_s, color="0.4", linestyle=":", linewidth=1)
        ax.set_title(scenario_label)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temperature rise (C)")
    axes.flat[0].legend(loc="best", fontsize=8)
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "temperature_profiles.png")
    plt.close(fig)


def plot_endpoint_summary(summary: pd.DataFrame, outdir: Path) -> None:
    """Plot end-of-train temperature rise by scenario and model."""
    outdir.mkdir(parents=True, exist_ok=True)
    pivot = summary.pivot(index="scenario_label", columns="dosimetry_model", values="end_temp_rise_C")
    pivot = pivot.loc[summary["scenario_label"].drop_duplicates()]
    pivot = pivot.loc[:, list(DOSIMETRY_MODELS)]
    ax = pivot.plot(kind="bar", figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE, color=bar_colors(len(DOSIMETRY_MODELS)))
    for index, container in enumerate(ax.containers):
        for patch in container:
            patch.set_hatch(bar_hatch(index))
            patch.set_edgecolor("#222222")
            patch.set_linewidth(0.6)
    ax.set_xlabel("")
    ax.set_ylabel("End-of-train temperature rise (C)")
    ax.legend(title="Dosimetry model")
    ax.figure.tight_layout()
    save_manuscript_figure(ax.figure, outdir / "end_temperature_rise.png")
    plt.close(ax.figure)


def write_outputs(summary: pd.DataFrame, profiles: pd.DataFrame, outdir: Path, make_plots: bool) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "dosimetry_model_summary.csv", index=False)
    profiles.to_csv(outdir / "dosimetry_temperature_profiles.csv", index=False)
    if make_plots:
        plot_temperature_profiles(profiles, outdir)
        plot_endpoint_summary(summary, outdir)


def print_console_summary(summary: pd.DataFrame) -> None:
    columns = [
        "scenario",
        "dosimetry_model",
        "absorbed_energy_density_mJ_mm3",
        "geometry_heat_density_mJ_mm3",
        "adiabatic_temp_rise_C",
        "end_temp_rise_C",
        "thermal_retention_factor",
    ]
    print(summary[columns].to_string(index=False, float_format=lambda value: f"{value:0.3f}"))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare nsPEF pulse/dosimetry models across heat-sensitive scenarios."
    )
    parser.add_argument(
        "--out",
        default="results/dosimetry_model_comparison",
        help="Directory for CSV and PNG outputs.",
    )
    parser.add_argument("--no-plots", action="store_true", help="Write CSV outputs only.")
    parser.add_argument(
        "--profile-points",
        type=int,
        default=240,
        help="Number of time points per scenario/model temperature profile.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    summary, profiles = build_comparison_tables(profile_points=args.profile_points)
    outdir = Path(args.out)
    write_outputs(summary, profiles, outdir, make_plots=not args.no_plots)
    print_console_summary(summary)
    print(f"\nWrote dosimetry comparison outputs to {outdir}")


if __name__ == "__main__":
    main()
