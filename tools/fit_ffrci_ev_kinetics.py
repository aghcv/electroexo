#!/usr/bin/env python3
"""Fit the existing EV-release model to the FFRi longitudinal EV totals.

This is deliberately a calibration driver, not a new model.  It imports the
current simulator unchanged, adjusts existing ``ev_release`` parameters, and
fits one shared parameter vector to the three measured 0.5/1/3-hour series.

Example
-------
MPLCONFIGDIR=/tmp/electroexo-mpl PYTHONPATH=. python \
    tools/fit_ffrci_ev_kinetics.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.optimize import least_squares


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from electro_exocytosis.config import SimulationScenario  # noqa: E402
from electro_exocytosis.experimental_bridge import (  # noqa: E402
    ExperimentalObservationBridge,
)
from electro_exocytosis.io.readers import load_default_parameters  # noqa: E402
from electro_exocytosis.simulation import Simulation  # noqa: E402
from tools.analyze_ffrci_data import parse_exoid_file, summarize_exoid  # noqa: E402


DEFAULT_DATA_DIR = REPO_ROOT / "data" / "experimental" / "ffrci_data_sharing"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "ffrci_ev_kinetics_fit"
EXOID_FILENAME = "Exoid_particle number and size in CD4 sham and different fields.csv"


@dataclass(frozen=True)
class ParameterSpec:
    """An existing EV-model parameter and its multiplicative fit bounds."""

    name: str
    module: str
    lower_multiplier: float
    upper_multiplier: float


# These six parameters control the three modeled EV pathways. Broad baseline
# rate bounds support both legacy concentration-scale fits and physically
# scaled single-cell observations. The experimental bridge should be preferred
# so baseline rates do not absorb cell-number and volume conversions.
FIT_PARAMETERS = (
    ParameterSpec("baseline_sEV_rate", "small_EV_MVB_release", 1.0e-6, 1.0e12),
    ParameterSpec("baseline_mlEV_rate", "medium_large_EV_budding", 1.0e-6, 1.0e12),
    ParameterSpec("baseline_AB_rate", "apoptotic_body_release", 1.0e-6, 1.0e12),
    ParameterSpec("k_ILV_release_s", "small_EV_MVB_release", 0.05, 20.0),
    ParameterSpec("k_budding_s", "medium_large_EV_budding", 0.05, 20.0),
    ParameterSpec("k_apoptotic_commitment_s", "apoptotic_body_release", 0.05, 20.0),
)


@dataclass(frozen=True)
class Exposure:
    condition: str
    pulse_count: int
    amplitude_kV_cm: float


EXPOSURES = (
    Exposure("1p20kV", 1, 20.0),
    Exposure("3p40kV", 3, 40.0),
    Exposure("5p40kV", 5, 40.0),
)


def load_replicate_observations(
    data_dir: Path,
    *,
    max_particle_diameter_nm: float | None = None,
    observation_bridge: ExperimentalObservationBridge | None = None,
) -> pd.DataFrame:
    """Return treated EV totals while retaining the file's p1/p2/p3 labels.

    When ``max_particle_diameter_nm`` is supplied, concentrations are summed
    only across bins whose reported particle diameter is strictly below the
    cutoff.  The strict comparison matches wording such as "smaller than
    200 nm" and avoids silently including a bin centered exactly at the
    boundary.
    """

    path = data_dir / EXOID_FILENAME
    if not path.exists():
        raise FileNotFoundError(f"Exoid CSV not found: {path}")

    particle_bins = parse_exoid_file(path)
    if max_particle_diameter_nm is not None:
        if max_particle_diameter_nm <= 0.0:
            raise ValueError("max_particle_diameter_nm must be positive")
        particle_bins = particle_bins.loc[
            particle_bins["particle_diameter_nm"] < max_particle_diameter_nm
        ].copy()
        if particle_bins.empty:
            raise ValueError(
                "Particle-diameter filter removed every reported Exoid bin"
            )
    sample_summary = summarize_exoid(particle_bins)
    replicate_totals = sample_summary.loc[
        sample_summary["sample_type"] == "treatment",
        [
            "condition",
            "harvest_time_h",
            "replicate",
            "summed_bin_concentration_particles_per_ml",
        ],
    ].rename(
        columns={
            "harvest_time_h": "time_h",
            "summed_bin_concentration_particles_per_ml": (
                "measured_concentration_particles_per_ml"
            ),
        }
    )
    exposure_table = pd.DataFrame(asdict(exposure) for exposure in EXPOSURES)
    observations = replicate_totals.merge(exposure_table, on="condition", how="inner")
    observations = observations.sort_values(
        ["pulse_count", "replicate", "time_h"]
    ).reset_index(drop=True)
    observations["particle_diameter_filter"] = (
        "all_reported_bins"
        if max_particle_diameter_nm is None
        else f"diameter_lt_{max_particle_diameter_nm:g}_nm"
    )
    measured = observations[
        "measured_concentration_particles_per_ml"
    ].to_numpy(dtype=float)
    if observation_bridge is None:
        observations["observed_concentration"] = measured
        observations["observation_unit"] = "particles_per_ml"
    else:
        observations["observed_concentration"] = (
            observation_bridge.concentration_to_particles_per_cell(measured)
        )
        observations["observation_unit"] = "particle_equivalents_per_cell"
    return observations


def load_longitudinal_observations(
    data_dir: Path,
    *,
    max_particle_diameter_nm: float | None = None,
    observation_bridge: ExperimentalObservationBridge | None = None,
) -> pd.DataFrame:
    """Return the nine condition/time means used for the joint fit."""

    replicate_observations = load_replicate_observations(
        data_dir,
        max_particle_diameter_nm=max_particle_diameter_nm,
        observation_bridge=observation_bridge,
    )
    observations = (
        replicate_observations.groupby(
            ["condition", "time_h", "pulse_count", "amplitude_kV_cm"], as_index=False
        )
        .agg(
            observed_concentration=("observed_concentration", "mean"),
            replicate_sd=("observed_concentration", "std"),
            measured_concentration_particles_per_ml=(
                "measured_concentration_particles_per_ml",
                "mean",
            ),
            measured_concentration_sd_particles_per_ml=(
                "measured_concentration_particles_per_ml",
                "std",
            ),
            n_replicates=("replicate", "nunique"),
        )
        .sort_values(["pulse_count", "time_h"])
        .reset_index(drop=True)
    )
    observations["particle_diameter_filter"] = (
        "all_reported_bins"
        if max_particle_diameter_nm is None
        else f"diameter_lt_{max_particle_diameter_nm:g}_nm"
    )
    observations["observation_unit"] = (
        "particles_per_ml"
        if observation_bridge is None
        else "particle_equivalents_per_cell"
    )

    expected = {
        (item.condition, time_h) for item in EXPOSURES for time_h in (0.5, 1.0, 3.0)
    }
    found = set(zip(observations["condition"], observations["time_h"], strict=False))
    if found != expected:
        raise ValueError(
            f"Expected nine condition/time observations; missing={expected - found}"
        )
    return observations


def build_scenario(
    exposure: Exposure, *, pulse_width_ns: float, repetition_rate_hz: float
) -> SimulationScenario:
    """Construct the current Ruben-experiment scenario for one exposure."""

    return SimulationScenario.model_validate(
        {
            "scenario": {
                "name": f"FFRCI fit: {exposure.condition}",
                "mode": "cell_based_electro_exocytosis",
            },
            "pulse": {
                "amplitude_kV_cm": exposure.amplitude_kV_cm,
                "pulse_width_ns": pulse_width_ns,
                "pulse_number": exposure.pulse_count,
                "repetition_rate_Hz": repetition_rate_hz,
                "waveform": "square",
            },
            "exposure": {
                "geometry": "cuvette",
                "medium_conductivity_S_m": 1.6,
                "temperature_C": 37.0,
                "cell_density_per_ml": 1000000.0,
            },
            "cell_state": {
                "cell_type": "generic",
                "membrane_modifier": 1.0,
                "calcium_handling_modifier": 1.0,
                "baseline_EV_release_modifier": 1.0,
                "stress_sensitivity_modifier": 1.0,
            },
            "simulation": {
                "t_start_s": 0.0,
                "t_end_s": 3.0 * 3600.0,
                "output_dt_s": 60.0,
                "numerical_method": "solve_ivp",
            },
        }
    )


class CurrentModelPredictor:
    """Run the unmodified current simulator for the three exposure conditions."""

    def __init__(
        self,
        observations: pd.DataFrame,
        *,
        pulse_width_ns: float,
        repetition_rate_hz: float,
    ) -> None:
        self.observations = observations.copy()
        self.parameters = load_default_parameters()
        self.defaults = {
            spec.name: float(self.parameters["ev_release"][spec.name])
            for spec in FIT_PARAMETERS
        }
        self.scenarios = {
            exposure.condition: build_scenario(
                exposure,
                pulse_width_ns=pulse_width_ns,
                repetition_rate_hz=repetition_rate_hz,
            )
            for exposure in EXPOSURES
        }
        self.evaluations = 0
        self.simulator_runs = 0
        self._basis_cache: dict[tuple[float, ...], np.ndarray] = {}

    def values_from_log_multipliers(
        self, log_multipliers: Iterable[float]
    ) -> dict[str, float]:
        return {
            spec.name: self.defaults[spec.name] * math.exp(float(value))
            for spec, value in zip(FIT_PARAMETERS, log_multipliers, strict=True)
        }

    def predict(self, ev_parameters: dict[str, float]) -> np.ndarray:
        dynamic_names = [spec.name for spec in FIT_PARAMETERS[3:]]
        cache_key = tuple(float(ev_parameters[name]) for name in dynamic_names)
        basis = self._basis_cache.get(cache_key)
        if basis is None:
            predictions: dict[tuple[str, float], tuple[float, float, float]] = {}
            basis_parameters = dict(ev_parameters)
            for spec in FIT_PARAMETERS[:3]:
                basis_parameters[spec.name] = 1.0
            override = {"ev_release": basis_parameters}

            condition_names = set(self.observations["condition"])
            for exposure in EXPOSURES:
                if exposure.condition not in condition_names:
                    continue
                result = Simulation(
                    self.scenarios[exposure.condition], params_override=override
                ).run()
                self.simulator_runs += 1
                frame = result.ev_timeseries
                time_s = frame["t"].to_numpy(dtype=float)
                subtype_cumulative = frame[
                    ["sEV_cumulative", "mlEV_cumulative", "AB_cumulative"]
                ].to_numpy(dtype=float)
                condition_times = self.observations.loc[
                    self.observations["condition"] == exposure.condition, "time_h"
                ]
                for time_h in condition_times:
                    predictions[(exposure.condition, float(time_h))] = tuple(
                        float(
                            np.interp(
                                float(time_h) * 3600.0,
                                time_s,
                                subtype_cumulative[:, column],
                            )
                        )
                        for column in range(3)
                    )

            basis = np.array(
                [
                    predictions[(row.condition, float(row.time_h))]
                    for row in self.observations.itertuples()
                ],
                dtype=float,
            )
            self._basis_cache[cache_key] = basis

        baseline_rates = np.array(
            [float(ev_parameters[spec.name]) for spec in FIT_PARAMETERS[:3]],
            dtype=float,
        )
        self.evaluations += 1
        return basis @ baseline_rates

    def residuals(self, log_multipliers: np.ndarray) -> np.ndarray:
        try:
            predicted = self.predict(self.values_from_log_multipliers(log_multipliers))
            observed = self.observations["observed_concentration"].to_numpy(dtype=float)
            if not np.all(np.isfinite(predicted)) or np.any(predicted <= 0.0):
                raise FloatingPointError("non-positive or non-finite prediction")
            return np.log10(predicted) - np.log10(observed)
        except Exception:
            # Keep a failed ODE solve away from the optimizer without hiding a
            # failure in the final selected parameter vector.
            return np.full(len(self.observations), 1.0e6, dtype=float)


def _metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    log_observed = np.log10(observed)
    log_predicted = np.log10(predicted)
    residual = log_predicted - log_observed
    denominator = float(np.sum((log_observed - np.mean(log_observed)) ** 2))
    return {
        "rmse_log10": float(np.sqrt(np.mean(residual**2))),
        "mae_log10": float(np.mean(np.abs(residual))),
        "median_absolute_percent_error": float(
            100.0 * np.median(np.abs(predicted - observed) / observed)
        ),
        "r2_log10": (
            float(1.0 - np.sum(residual**2) / denominator) if denominator else math.nan
        ),
    }


def fit_model(
    predictor: CurrentModelPredictor,
    *,
    starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[object, np.ndarray, np.ndarray, list[object]]:
    """Run deterministic multistart bounded least squares in log space."""

    defaults = predictor.defaults
    initial_prediction = predictor.predict(defaults)
    observed = predictor.observations["observed_concentration"].to_numpy(dtype=float)
    common_scale = float(np.exp(np.mean(np.log(observed) - np.log(initial_prediction))))

    lower = np.log([spec.lower_multiplier for spec in FIT_PARAMETERS])
    upper = np.log([spec.upper_multiplier for spec in FIT_PARAMETERS])
    primary = np.zeros(len(FIT_PARAMETERS), dtype=float)
    primary[:3] = math.log(common_scale)
    primary = np.clip(primary, lower + 1.0e-8, upper - 1.0e-8)

    rng = np.random.default_rng(seed)
    guesses = [primary]
    for _ in range(max(0, starts - 1)):
        guess = primary.copy()
        guess[:3] += rng.normal(0.0, 1.5, size=3)
        guess[3:] += rng.normal(0.0, 0.9, size=len(FIT_PARAMETERS) - 3)
        guesses.append(np.clip(guess, lower + 1.0e-8, upper - 1.0e-8))

    best = None
    candidates = []
    for guess in guesses:
        candidate = least_squares(
            predictor.residuals,
            guess,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            diff_step=0.025,
            max_nfev=max_nfev,
            ftol=1.0e-8,
            xtol=1.0e-8,
            gtol=1.0e-8,
        )
        candidates.append(candidate)
        if best is None or candidate.cost < best.cost:
            best = candidate

    if best is None:
        raise RuntimeError("Optimizer did not produce a candidate")
    fitted_prediction = predictor.predict(predictor.values_from_log_multipliers(best.x))
    return best, initial_prediction, fitted_prediction, candidates


def _write_outputs(
    output_dir: Path,
    predictor: CurrentModelPredictor,
    optimizer_result: object,
    initial_prediction: np.ndarray,
    fitted_prediction: np.ndarray,
    *,
    starts: int,
    max_nfev: int,
    pulse_width_ns: float,
    repetition_rate_hz: float,
    max_particle_diameter_nm: float | None,
    observation_bridge: ExperimentalObservationBridge | None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    observed = predictor.observations["observed_concentration"].to_numpy(dtype=float)
    fitted_values = predictor.values_from_log_multipliers(optimizer_result.x)

    parameter_rows = []
    for spec in FIT_PARAMETERS:
        default = predictor.defaults[spec.name]
        fitted = fitted_values[spec.name]
        parameter_rows.append(
            {
                "parameter": spec.name,
                "module": spec.module,
                "default": default,
                "fitted": fitted,
                "multiplier": fitted / default,
                "lower_multiplier": spec.lower_multiplier,
                "upper_multiplier": spec.upper_multiplier,
            }
        )
    pd.DataFrame(parameter_rows).to_csv(
        output_dir / "fitted_parameters.csv", index=False
    )
    with (output_dir / "fitted_parameters.yml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump({"ev_release": fitted_values}, handle, sort_keys=False)

    comparison = predictor.observations.copy()
    comparison["predicted_default"] = initial_prediction
    comparison["predicted_fitted"] = fitted_prediction
    comparison["fitted_over_observed"] = fitted_prediction / observed
    comparison["log10_residual"] = np.log10(fitted_prediction) - np.log10(observed)
    comparison.to_csv(output_dir / "observed_vs_predicted.csv", index=False)

    monotonic_conflicts = []
    for exposure in EXPOSURES:
        values = comparison.loc[
            comparison["condition"] == exposure.condition, "observed_concentration"
        ].to_numpy(dtype=float)
        if np.any(np.diff(values) < 0.0):
            monotonic_conflicts.append(exposure.condition)

    summary = {
        "data_points": int(len(observed)),
        "conditions": [item.condition for item in EXPOSURES],
        "shared_fitted_parameters": [spec.name for spec in FIT_PARAMETERS],
        "scenario_assumptions": {
            "pulse_width_ns": pulse_width_ns,
            "repetition_rate_Hz": repetition_rate_hz,
            "amplitude_labels_interpreted_as_kV_cm": True,
        },
        "experimental_observation_filter": {
            "particle_diameter_nm": (
                "all reported bins"
                if max_particle_diameter_nm is None
                else f"strictly less than {max_particle_diameter_nm:g} nm"
            )
        },
        "experimental_observation_bridge": (
            None
            if observation_bridge is None
            else observation_bridge.to_metadata()
        ),
        "optimizer": {
            "method": "bounded multistart scipy.optimize.least_squares",
            "objective": "unweighted residuals in log10 observation units",
            "starts": starts,
            "max_function_evaluations_per_start": max_nfev,
            "success": bool(optimizer_result.success),
            "status": int(optimizer_result.status),
            "message": str(optimizer_result.message),
            "cost": float(optimizer_result.cost),
            "selected_start_nfev": int(optimizer_result.nfev),
            "simulator_parameter_evaluations_total": predictor.evaluations,
            "simulator_runs_total": predictor.simulator_runs,
        },
        "default_metrics": _metrics(observed, initial_prediction),
        "fitted_metrics": _metrics(observed, fitted_prediction),
        "structural_note": {
            "observed_conditions_with_a_decline": monotonic_conflicts,
            "meaning": (
                "The current outputs are cumulative EV release and therefore non-decreasing. "
                "A decrease in a measured concentration trajectory cannot be reproduced exactly "
                "by parameter fitting alone."
            ),
        },
    }
    (output_dir / "fit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    def plot_observed_and_fitted(*, logarithmic: bool, filename: str) -> None:
        fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), sharey=False)
        for axis, exposure in zip(axes, EXPOSURES, strict=True):
            subset = comparison[comparison["condition"] == exposure.condition]
            axis.plot(
                subset["time_h"],
                subset["observed_concentration"],
                "o-",
                lw=2,
                label="Observed",
            )
            axis.plot(
                subset["time_h"],
                subset["predicted_fitted"],
                "s--",
                lw=2,
                label="Fitted model",
            )
            visible_values = np.concatenate(
                [
                    subset["observed_concentration"].to_numpy(dtype=float),
                    subset["predicted_fitted"].to_numpy(dtype=float),
                ]
            )
            if logarithmic:
                axis.set_yscale("log")
                axis.set_ylim(visible_values.min() / 1.18, visible_values.max() * 1.18)
            else:
                span = float(visible_values.max() - visible_values.min())
                padding = max(0.08 * span, 0.03 * float(visible_values.max()))
                axis.set_ylim(
                    max(0.0, float(visible_values.min()) - padding),
                    float(visible_values.max()) + padding,
                )
            axis.set_title(exposure.condition)
            axis.set_xlabel("Time after nsPEF (h)")
            axis.grid(alpha=0.25)
        if observation_bridge is not None:
            concentration_label = (
                "Particle equivalents per initial cell"
                if observation_bridge.cell_basis == "initial"
                else "Particle equivalents per viable cell"
            )
            if max_particle_diameter_nm is not None:
                concentration_label += f" (<{max_particle_diameter_nm:g} nm)"
        else:
            concentration_label = (
                "Total EV concentration (particles/mL)"
                if max_particle_diameter_nm is None
                else f"Particle concentration <{max_particle_diameter_nm:g} nm (particles/mL)"
            )
        axes[0].set_ylabel(concentration_label)
        axes[0].legend(frameon=False, fontsize=8)
        scale_label = "logarithmic" if logarithmic else "linear"
        size_label = (
            "all reported sizes"
            if max_particle_diameter_nm is None
            else f"particles <{max_particle_diameter_nm:g} nm"
        )
        observation_label = (
            "single-cell-equivalent EV output"
            if observation_bridge is not None
            else "EV concentration"
        )
        fig.suptitle(
            f"Observed and fitted {observation_label}: {size_label} ({scale_label} scale)"
        )
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=220, bbox_inches="tight")
        plt.close(fig)

    plot_observed_and_fitted(logarithmic=True, filename="longitudinal_fit.png")
    plot_observed_and_fitted(logarithmic=True, filename="longitudinal_fit_log.png")
    plot_observed_and_fitted(logarithmic=False, filename="longitudinal_fit_linear.png")

    readme = f"""# FFRi longitudinal EV kinetic fit

This directory was generated by `tools/fit_ffrci_ev_kinetics.py`. The existing
simulation and EV-release model source were not modified. One parameter vector
was fit jointly to all nine treated condition/time totals.

- Objective: unweighted log10 residuals in {"particle equivalents per cell" if observation_bridge is not None else "particles per mL"}
- Pulse assumptions: {pulse_width_ns:g} ns, {repetition_rate_hz:g} Hz
- Label interpretation: 20/40 kV labels were treated as 20/40 kV/cm
- Experimental size filter: {"all reported bins" if max_particle_diameter_nm is None else f"particle diameter strictly below {max_particle_diameter_nm:g} nm"}
- Experimental bridge: {"not applied" if observation_bridge is None else f"{observation_bridge.initial_cell_count:g} initial cells in {observation_bridge.medium_volume_ml:g} mL; {observation_bridge.cell_basis}-cell basis"}
- Fit output: `observed_vs_predicted.csv`
- Parameter override: `fitted_parameters.yml`
- Diagnostics: `fit_summary.json`, `longitudinal_fit_log.png`, and
  `longitudinal_fit_linear.png`

Important interpretation: the model output is cumulative released EVs and must
increase with time. The observed concentration decreases in at least one
interval for {', '.join(monotonic_conflicts)}. Those decreases cannot be matched
exactly by changing current parameter values; the fitted curves are therefore
the best compromise available to the present cumulative-output model.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")


def fit_sample_sets(
    replicate_observations: pd.DataFrame,
    output_dir: Path,
    *,
    starts: int,
    max_nfev: int,
    seed: int,
    pulse_width_ns: float,
    repetition_rate_hz: float,
    near_optimal_delta_log10: float,
    max_particle_diameter_nm: float | None,
    observation_bridge: ExperimentalObservationBridge | None,
) -> dict[str, object]:
    """Fit each condition/replicate trajectory and summarize parameter ranges."""

    fit_cache: dict[tuple[object, ...], dict[str, object]] = {}
    comparison_frames: list[pd.DataFrame] = []
    best_parameter_rows: list[dict[str, object]] = []
    ensemble_rows: list[dict[str, object]] = []
    fit_rows: list[dict[str, object]] = []
    total_simulator_runs = 0

    grouped = replicate_observations.groupby(["condition", "replicate"], sort=True)
    for group_index, ((condition, replicate), group) in enumerate(grouped):
        group = group.sort_values("time_h").reset_index(drop=True)
        signature = (
            condition,
            tuple(group["time_h"].astype(float)),
            # Collapse only numerical summation noise while retaining genuine
            # replicate differences in either concentration or per-cell units.
            tuple(
                float(f"{float(value):.9g}")
                for value in group["observed_concentration"]
            ),
        )
        reused = signature in fit_cache
        if not reused:
            predictor = CurrentModelPredictor(
                group,
                pulse_width_ns=pulse_width_ns,
                repetition_rate_hz=repetition_rate_hz,
            )
            best, initial, fitted, candidates = fit_model(
                predictor,
                starts=starts,
                max_nfev=max_nfev,
                seed=seed + group_index,
            )
            fit_cache[signature] = {
                "predictor": predictor,
                "best": best,
                "initial": initial,
                "fitted": fitted,
                "candidates": candidates,
            }
            total_simulator_runs += predictor.simulator_runs
        fit_data = fit_cache[signature]
        predictor = fit_data["predictor"]
        best = fit_data["best"]
        initial = np.asarray(fit_data["initial"], dtype=float)
        fitted = np.asarray(fit_data["fitted"], dtype=float)
        candidates = fit_data["candidates"]
        observed = group["observed_concentration"].to_numpy(dtype=float)
        best_metrics = _metrics(observed, fitted)

        comparison = group.copy()
        comparison["predicted_default"] = initial
        comparison["predicted_fitted"] = fitted
        comparison["fitted_over_observed"] = fitted / observed
        comparison["log10_residual"] = np.log10(fitted) - np.log10(observed)
        comparison["fit_reused_for_identical_trajectory"] = reused
        comparison_frames.append(comparison)

        fit_rows.append(
            {
                "condition": condition,
                "replicate": int(replicate),
                "n_time_points": int(len(group)),
                "times_h": ";".join(f"{value:g}" for value in group["time_h"]),
                "rmse_log10": best_metrics["rmse_log10"],
                "median_absolute_percent_error": best_metrics[
                    "median_absolute_percent_error"
                ],
                "optimizer_success": bool(best.success),
                "fit_reused_for_identical_trajectory": reused,
            }
        )

        best_values = predictor.values_from_log_multipliers(best.x)
        for spec in FIT_PARAMETERS:
            best_parameter_rows.append(
                {
                    "condition": condition,
                    "replicate": int(replicate),
                    "n_time_points": int(len(group)),
                    "module": spec.module,
                    "parameter": spec.name,
                    "default": predictor.defaults[spec.name],
                    "fitted": best_values[spec.name],
                    "multiplier": best_values[spec.name]
                    / predictor.defaults[spec.name],
                    "rmse_log10": best_metrics["rmse_log10"],
                    "fit_reused_for_identical_trajectory": reused,
                }
            )

        candidate_records = []
        for start_index, candidate in enumerate(candidates, start=1):
            values = predictor.values_from_log_multipliers(candidate.x)
            candidate_prediction = predictor.predict(values)
            candidate_rmse = _metrics(observed, candidate_prediction)["rmse_log10"]
            candidate_records.append((start_index, candidate, values, candidate_rmse))
        best_candidate_rmse = min(record[3] for record in candidate_records)
        for start_index, candidate, values, candidate_rmse in candidate_records:
            near_optimal = (
                candidate_rmse <= best_candidate_rmse + near_optimal_delta_log10
            )
            for spec in FIT_PARAMETERS:
                ensemble_rows.append(
                    {
                        "condition": condition,
                        "replicate": int(replicate),
                        "start": start_index,
                        "module": spec.module,
                        "parameter": spec.name,
                        "value": values[spec.name],
                        "multiplier": values[spec.name] / predictor.defaults[spec.name],
                        "rmse_log10": candidate_rmse,
                        "near_optimal": near_optimal,
                        "optimizer_success": bool(candidate.success),
                    }
                )

    comparisons = pd.concat(comparison_frames, ignore_index=True)
    best_parameters = pd.DataFrame(best_parameter_rows)
    ensemble = pd.DataFrame(ensemble_rows)
    fit_summary = pd.DataFrame(fit_rows)
    comparisons.to_csv(output_dir / "sample_set_observed_vs_predicted.csv", index=False)
    best_parameters.to_csv(output_dir / "sample_set_best_parameters.csv", index=False)
    ensemble.to_csv(output_dir / "sample_set_parameter_ensemble.csv", index=False)
    fit_summary.to_csv(output_dir / "sample_set_fit_metrics.csv", index=False)

    near_optimal = ensemble[ensemble["near_optimal"]].copy()
    near_optimal_ranges = near_optimal.groupby(
        ["condition", "replicate", "module", "parameter"], as_index=False
    ).agg(
        n_near_optimal_solutions=("value", "count"),
        minimum=("value", "min"),
        median=("value", "median"),
        maximum=("value", "max"),
        minimum_multiplier=("multiplier", "min"),
        median_multiplier=("multiplier", "median"),
        maximum_multiplier=("multiplier", "max"),
    )
    near_optimal_ranges.to_csv(
        output_dir / "sample_set_near_optimal_ranges.csv", index=False
    )

    distribution_frames = []
    for scope, scoped in [
        ("all_sample_sets", best_parameters),
        *[
            (condition, best_parameters[best_parameters["condition"] == condition])
            for condition in best_parameters["condition"].unique()
        ],
    ]:
        distribution = scoped.groupby(["module", "parameter"], as_index=False).agg(
            n_sample_sets=("fitted", "count"),
            minimum=("fitted", "min"),
            q25=("fitted", lambda values: values.quantile(0.25)),
            median=("fitted", "median"),
            q75=("fitted", lambda values: values.quantile(0.75)),
            maximum=("fitted", "max"),
            minimum_multiplier=("multiplier", "min"),
            median_multiplier=("multiplier", "median"),
            maximum_multiplier=("multiplier", "max"),
        )
        distribution.insert(0, "scope", scope)
        distribution_frames.append(distribution)
    parameter_distributions = pd.concat(distribution_frames, ignore_index=True)
    parameter_distributions.to_csv(
        output_dir / "constitutive_parameter_distributions.csv", index=False
    )

    fig, axes = plt.subplots(3, 3, figsize=(11.5, 9.0), sharex=True, sharey=True)
    for row_index, exposure in enumerate(EXPOSURES):
        for column_index, replicate in enumerate((1, 2, 3)):
            axis = axes[row_index, column_index]
            subset = comparisons[
                (comparisons["condition"] == exposure.condition)
                & (comparisons["replicate"] == replicate)
            ]
            if subset.empty:
                axis.set_axis_off()
                continue
            axis.plot(
                subset["time_h"],
                subset["observed_concentration"],
                "o-",
                label="Observed",
            )
            axis.plot(
                subset["time_h"], subset["predicted_fitted"], "s--", label="Fitted"
            )
            axis.set_yscale("log")
            axis.grid(alpha=0.25)
            axis.set_title(f"{exposure.condition} p{replicate}")
            if row_index == 2:
                axis.set_xlabel("Time after nsPEF (h)")
            if column_index == 0:
                if observation_bridge is not None:
                    axis.set_ylabel(
                        "Particles / initial cell"
                        if observation_bridge.cell_basis == "initial"
                        else "Particles / viable cell"
                    )
                else:
                    axis.set_ylabel(
                        "EV concentration"
                        if max_particle_diameter_nm is None
                        else f"Concentration <{max_particle_diameter_nm:g} nm"
                    )
    axes[0, 0].legend(frameon=False, fontsize=8)
    size_label = (
        "all reported particle sizes"
        if max_particle_diameter_nm is None
        else f"particles <{max_particle_diameter_nm:g} nm"
    )
    fig.suptitle(
        f"Per sample set fits using the current cumulative EV model: {size_label}"
    )
    fig.tight_layout()
    fig.savefig(
        output_dir / "sample_set_longitudinal_fits.png", dpi=220, bbox_inches="tight"
    )
    plt.close(fig)

    per_cell = replicate_observations.groupby(["condition", "time_h"])[
        "observed_concentration"
    ]
    per_cell_relative_spread = (per_cell.max() - per_cell.min()) / per_cell.mean()
    indistinguishable_cells = int((per_cell_relative_spread <= 1.0e-8).sum())
    all_cells_indistinguishable = indistinguishable_cells == len(
        per_cell_relative_spread
    )
    missing_points = [
        {"condition": exposure.condition, "replicate": replicate, "time_h": time_h}
        for exposure in EXPOSURES
        for replicate in (1, 2, 3)
        for time_h in (0.5, 1.0, 3.0)
        if not (
            (replicate_observations["condition"] == exposure.condition)
            & (replicate_observations["replicate"] == replicate)
            & (replicate_observations["time_h"] == time_h)
        ).any()
    ]
    summary = {
        "nominal_sample_sets": int(len(fit_summary)),
        "unique_observed_trajectories_optimized": int(len(fit_cache)),
        "complete_three_time_point_sets": int(
            (fit_summary["n_time_points"] == 3).sum()
        ),
        "missing_condition_replicate_time_points": missing_points,
        "condition_time_cells_with_numerically_indistinguishable_replicate_totals": (
            indistinguishable_cells
        ),
        "condition_time_cells_total": int(len(per_cell_relative_spread)),
        "maximum_relative_spread_across_replicate_totals": float(
            per_cell_relative_spread.max()
        ),
        "near_optimal_definition": (
            f"multistart solution RMSE within {near_optimal_delta_log10:g} log10 units of the "
            "best solution for that sample set"
        ),
        "experimental_observation_filter": {
            "particle_diameter_nm": (
                "all reported bins"
                if max_particle_diameter_nm is None
                else f"strictly less than {max_particle_diameter_nm:g} nm"
            )
        },
        "experimental_observation_bridge": (
            None
            if observation_bridge is None
            else observation_bridge.to_metadata()
        ),
        "simulator_runs_total_after_duplicate_trajectory_reuse": total_simulator_runs,
        "interpretation": (
            "The p labels can be aligned into nominal trajectories, but the file does not establish "
            "whether they are repeated measurements of the same culture, matched biological "
            "replicates, or independent endpoint wells. All available replicate totals are "
            "numerically indistinguishable within each condition/time cell (relative spread below "
            "1e-8), so between-replicate parameter ranges do not represent measured biological "
            "variability."
            if all_cells_indistinguishable
            else (
                "The p labels can be aligned into nominal trajectories, but the file does not "
                "establish whether they are repeated measurements of the same culture, matched "
                "biological replicates, or independent endpoint wells. The size-filtered totals "
                "differ across p labels because their reported size distributions differ; these "
                "parameter ranges should not be interpreted as biological population distributions "
                "until the replicate structure is confirmed."
            )
        ),
    }
    (output_dir / "sample_set_fit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--fit-mode", choices=("joint", "sample-sets", "both"), default="both"
    )
    parser.add_argument(
        "--experimental-bridge-config",
        type=Path,
        default=None,
        help=(
            "YAML configuration that converts experimental particles/mL to "
            "single-cell-equivalent model observations"
        ),
    )
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--max-nfev", type=int, default=60)
    parser.add_argument(
        "--near-optimal-delta-log10",
        type=float,
        default=0.02,
        help="RMSE tolerance used to retain multistart sample-set solutions",
    )
    parser.add_argument("--seed", type=int, default=20260903)
    parser.add_argument("--pulse-width-ns", type=float, default=60.0)
    parser.add_argument("--repetition-rate-hz", type=float, default=1.0)
    parser.add_argument(
        "--max-particle-diameter-nm",
        type=float,
        default=None,
        help=(
            "Sum only experimental Exoid bins with particle diameter strictly "
            "below this cutoff; the simulated EV model is otherwise unchanged"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.starts < 1 or args.max_nfev < 1:
        raise ValueError("--starts and --max-nfev must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    observation_bridge = (
        None
        if args.experimental_bridge_config is None
        else ExperimentalObservationBridge.from_yaml(args.experimental_bridge_config)
    )
    if args.fit_mode in {"joint", "both"}:
        observations = load_longitudinal_observations(
            args.data_dir,
            max_particle_diameter_nm=args.max_particle_diameter_nm,
            observation_bridge=observation_bridge,
        )
        predictor = CurrentModelPredictor(
            observations,
            pulse_width_ns=args.pulse_width_ns,
            repetition_rate_hz=args.repetition_rate_hz,
        )
        result, initial_prediction, fitted_prediction, _ = fit_model(
            predictor,
            starts=args.starts,
            max_nfev=args.max_nfev,
            seed=args.seed,
        )
        _write_outputs(
            args.output_dir,
            predictor,
            result,
            initial_prediction,
            fitted_prediction,
            starts=args.starts,
            max_nfev=args.max_nfev,
            pulse_width_ns=args.pulse_width_ns,
            repetition_rate_hz=args.repetition_rate_hz,
            max_particle_diameter_nm=args.max_particle_diameter_nm,
            observation_bridge=observation_bridge,
        )
        fitted_metrics = _metrics(
            observations["observed_concentration"].to_numpy(dtype=float),
            fitted_prediction,
        )
        print(f"Joint optimizer success: {result.success} ({result.message})")
        print(f"Joint fitted log10 RMSE: {fitted_metrics['rmse_log10']:.4f}")
        print(
            "Joint median absolute percent error: "
            f"{fitted_metrics['median_absolute_percent_error']:.1f}%"
        )

    if args.fit_mode in {"sample-sets", "both"}:
        sample_summary = fit_sample_sets(
            load_replicate_observations(
                args.data_dir,
                max_particle_diameter_nm=args.max_particle_diameter_nm,
                observation_bridge=observation_bridge,
            ),
            args.output_dir,
            starts=args.starts,
            max_nfev=args.max_nfev,
            seed=args.seed,
            pulse_width_ns=args.pulse_width_ns,
            repetition_rate_hz=args.repetition_rate_hz,
            near_optimal_delta_log10=args.near_optimal_delta_log10,
            max_particle_diameter_nm=args.max_particle_diameter_nm,
            observation_bridge=observation_bridge,
        )
        print(
            "Sample-set fits: "
            f"{sample_summary['nominal_sample_sets']} nominal sets, "
            f"{sample_summary['unique_observed_trajectories_optimized']} unique trajectories"
        )
    print(f"Fit written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
