#!/usr/bin/env python3
"""Fit stimulated-to-sham EV concentration ratios with the current EV model.

The Exoid file contains two control labels without harvest times.  This driver
therefore treats each control type as a separate, time-invariant concentration
reference and fits the stimulated/sham fold change.  Model predictions are
normalized to a near-zero-field sham proxy simulated at the same 0.5, 1, and
3 hour times.  No model source equations are changed.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from electro_exocytosis.io.readers import load_default_parameters  # noqa: E402
from electro_exocytosis.experimental_bridge import (  # noqa: E402
    ExperimentalObservationBridge,
)
from electro_exocytosis.simulation import Simulation  # noqa: E402
from tools.analyze_ffrci_data import parse_exoid_file, summarize_exoid  # noqa: E402
from tools.fit_ffrci_ev_kinetics import (  # noqa: E402
    DEFAULT_DATA_DIR,
    EXOID_FILENAME,
    EXPOSURES,
    FIT_PARAMETERS,
    Exposure,
    build_scenario,
)


DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "ffrci_ev_kinetics_normalized_lt200nm"
CONTROL_ORDER = ("sham2", "sham_media")
CONTROL_LABELS = {"sham2": "Sham2", "sham_media": "sham media"}

# A ratio cannot identify the common scale of all three baseline release rates.
# Fixing baseline_sEV_rate is a gauge choice; the remaining two baseline rates
# estimate pathway mixture relative to that reference.
FREE_PARAMETER_NAMES = (
    "baseline_mlEV_rate",
    "baseline_AB_rate",
    "k_ILV_release_s",
    "k_budding_s",
    "k_apoptotic_commitment_s",
)
FREE_PARAMETER_BOUNDS = {
    "baseline_mlEV_rate": (1.0e-6, 1.0e6),
    "baseline_AB_rate": (1.0e-6, 1.0e6),
    "k_ILV_release_s": (0.05, 20.0),
    "k_budding_s": (0.05, 20.0),
    "k_apoptotic_commitment_s": (0.05, 20.0),
}


def load_normalized_observations(
    data_dir: Path,
    *,
    max_particle_diameter_nm: float | None,
    observation_bridge: ExperimentalObservationBridge | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return condition means normalized separately to each control type."""

    bins = parse_exoid_file(data_dir / EXOID_FILENAME)
    if max_particle_diameter_nm is not None:
        if max_particle_diameter_nm <= 0.0:
            raise ValueError("max_particle_diameter_nm must be positive")
        bins = bins.loc[
            bins["particle_diameter_nm"] < max_particle_diameter_nm
        ].copy()
    samples = summarize_exoid(bins)
    concentration = "summed_bin_concentration_particles_per_ml"

    treated = (
        samples.loc[samples["sample_type"] == "treatment"]
        .groupby(
            ["condition", "harvest_time_h", "pulse_count", "amplitude_label_kV"],
            as_index=False,
        )
        .agg(
            stimulated_concentration=(concentration, "mean"),
            stimulated_sd=(concentration, "std"),
            stimulated_n=("dataset_number", "count"),
        )
        .rename(
            columns={
                "harvest_time_h": "time_h",
                "amplitude_label_kV": "amplitude_kV_cm",
            }
        )
    )
    controls = (
        samples.loc[samples["sample_type"] == "control"]
        .groupby("condition", as_index=False)
        .agg(
            sham_concentration=(concentration, "mean"),
            sham_sd=(concentration, "std"),
            sham_n=("dataset_number", "count"),
        )
        .rename(columns={"condition": "control"})
    )
    if observation_bridge is None:
        treated["stimulated_particles_per_cell"] = treated[
            "stimulated_concentration"
        ]
        controls["sham_particles_per_cell"] = controls["sham_concentration"]
        observation_unit = "particles_per_ml"
    else:
        treated["stimulated_particles_per_cell"] = (
            observation_bridge.concentration_to_particles_per_cell(
                treated["stimulated_concentration"].to_numpy(dtype=float)
            )
        )
        controls["sham_particles_per_cell"] = (
            observation_bridge.concentration_to_particles_per_cell(
                controls["sham_concentration"].to_numpy(dtype=float)
            )
        )
        observation_unit = "particle_equivalents_per_cell"
    missing = set(CONTROL_ORDER) - set(controls["control"])
    if missing:
        raise ValueError(f"Missing expected control labels: {sorted(missing)}")

    frames = []
    for control in CONTROL_ORDER:
        control_row = controls.loc[controls["control"] == control].iloc[0]
        normalized = treated.copy()
        normalized.insert(0, "control", control)
        normalized["sham_concentration"] = float(control_row["sham_concentration"])
        normalized["sham_particles_per_cell"] = float(
            control_row["sham_particles_per_cell"]
        )
        normalized["sham_sd"] = float(control_row["sham_sd"])
        normalized["sham_n"] = int(control_row["sham_n"])
        normalized["observed_fold_of_sham"] = (
            normalized["stimulated_particles_per_cell"]
            / normalized["sham_particles_per_cell"]
        )
        normalized["observation_unit_before_ratio"] = observation_unit
        frames.append(normalized)
    observations = pd.concat(frames, ignore_index=True).sort_values(
        ["control", "pulse_count", "time_h"]
    )
    return observations.reset_index(drop=True), controls


class NormalizedModelPredictor:
    """Predict stimulated/sham ratios using the unchanged current simulator."""

    def __init__(
        self,
        observations: pd.DataFrame,
        *,
        pulse_width_ns: float,
        repetition_rate_hz: float,
    ) -> None:
        self.observations = observations.copy().reset_index(drop=True)
        self.defaults = {
            name: float(load_default_parameters()["ev_release"][name])
            for name in ("baseline_sEV_rate", *FREE_PARAMETER_NAMES)
        }
        # The scenario schema requires a positive field and at least one pulse;
        # 1e-9 kV/cm is used as a numerically negligible unexposed proxy.
        self.sham_proxy = Exposure("sham_proxy", 1, 1.0e-9)
        self.scenarios = {
            exposure.condition: build_scenario(
                exposure,
                pulse_width_ns=pulse_width_ns,
                repetition_rate_hz=repetition_rate_hz,
            )
            for exposure in (*EXPOSURES, self.sham_proxy)
        }
        self._basis_cache: dict[tuple[float, ...], tuple[np.ndarray, np.ndarray]] = {}
        self.simulator_runs = 0
        self.evaluations = 0

    def values_from_log_multipliers(
        self, log_multipliers: np.ndarray
    ) -> dict[str, float]:
        values = dict(self.defaults)
        for name, log_multiplier in zip(
            FREE_PARAMETER_NAMES, log_multipliers, strict=True
        ):
            values[name] = self.defaults[name] * math.exp(float(log_multiplier))
        return values

    def predict(self, parameters: dict[str, float]) -> np.ndarray:
        dynamic_names = FREE_PARAMETER_NAMES[2:]
        cache_key = tuple(float(parameters[name]) for name in dynamic_names)
        cached = self._basis_cache.get(cache_key)
        if cached is None:
            basis_parameters = dict(parameters)
            for spec in FIT_PARAMETERS[:3]:
                basis_parameters[spec.name] = 1.0
            override = {"ev_release": basis_parameters}
            times_h = sorted(set(self.observations["time_h"].astype(float)))
            times_s = np.asarray(times_h) * 3600.0
            by_condition: dict[str, np.ndarray] = {}
            for exposure in (*EXPOSURES, self.sham_proxy):
                result = Simulation(
                    self.scenarios[exposure.condition], params_override=override
                ).run()
                self.simulator_runs += 1
                frame = result.ev_timeseries
                by_condition[exposure.condition] = np.column_stack(
                    [
                        np.interp(times_s, frame["t"], frame[column])
                        for column in (
                            "sEV_cumulative",
                            "mlEV_cumulative",
                            "AB_cumulative",
                        )
                    ]
                )
            time_index = {time_h: index for index, time_h in enumerate(times_h)}
            treated_basis = np.vstack(
                [
                    by_condition[row.condition][time_index[float(row.time_h)]]
                    for row in self.observations.itertuples()
                ]
            )
            sham_basis = np.vstack(
                [
                    by_condition[self.sham_proxy.condition][
                        time_index[float(row.time_h)]
                    ]
                    for row in self.observations.itertuples()
                ]
            )
            cached = (treated_basis, sham_basis)
            self._basis_cache[cache_key] = cached

        rates = np.asarray(
            [
                parameters["baseline_sEV_rate"],
                parameters["baseline_mlEV_rate"],
                parameters["baseline_AB_rate"],
            ],
            dtype=float,
        )
        treated_basis, sham_basis = cached
        self.evaluations += 1
        return (treated_basis @ rates) / (sham_basis @ rates)

    def residuals(self, log_multipliers: np.ndarray) -> np.ndarray:
        try:
            predicted = self.predict(
                self.values_from_log_multipliers(log_multipliers)
            )
            observed = self.observations["observed_fold_of_sham"].to_numpy(
                dtype=float
            )
            if np.any(predicted <= 0.0) or not np.all(np.isfinite(predicted)):
                raise FloatingPointError("Invalid normalized prediction")
            return np.log10(predicted) - np.log10(observed)
        except Exception:
            return np.full(len(self.observations), 1.0e6, dtype=float)


def fit_normalized(
    predictor: NormalizedModelPredictor, *, starts: int, max_nfev: int, seed: int
) -> tuple[object, list[object]]:
    lower = np.log(
        [FREE_PARAMETER_BOUNDS[name][0] for name in FREE_PARAMETER_NAMES]
    )
    upper = np.log(
        [FREE_PARAMETER_BOUNDS[name][1] for name in FREE_PARAMETER_NAMES]
    )
    rng = np.random.default_rng(seed)
    guesses = [np.zeros(len(FREE_PARAMETER_NAMES), dtype=float)]
    for _ in range(starts - 1):
        guesses.append(rng.uniform(lower, upper))
    candidates = [
        least_squares(
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
        for guess in guesses
    ]
    return min(candidates, key=lambda candidate: candidate.cost), candidates


def ratio_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = np.log10(predicted) - np.log10(observed)
    return {
        "rmse_log10_fold": float(np.sqrt(np.mean(residual**2))),
        "mae_log10_fold": float(np.mean(np.abs(residual))),
        "median_absolute_percent_error": float(
            100.0 * np.median(np.abs(predicted - observed) / observed)
        ),
        "maximum_absolute_fold_error": float(
            np.max(np.maximum(predicted / observed, observed / predicted))
        ),
    }


def plot_normalized_fits(
    comparison: pd.DataFrame,
    output_dir: Path,
    *,
    max_particle_diameter_nm: float | None,
) -> None:
    for logarithmic, filename in (
        (False, "normalized_observed_vs_fitted_linear.png"),
        (True, "normalized_observed_vs_fitted_log.png"),
    ):
        fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6), sharex=True)
        for row_index, control in enumerate(CONTROL_ORDER):
            for column_index, exposure in enumerate(EXPOSURES):
                axis = axes[row_index, column_index]
                subset = comparison.loc[
                    (comparison["control"] == control)
                    & (comparison["condition"] == exposure.condition)
                ].sort_values("time_h")
                axis.plot(
                    subset["time_h"],
                    subset["observed_fold_of_sham"],
                    "o-",
                    lw=2,
                    label="Observed stimulated/sham",
                )
                axis.plot(
                    subset["time_h"],
                    subset["predicted_fold_of_sham"],
                    "s--",
                    lw=2,
                    label="Fitted model stimulated/sham",
                )
                axis.axhline(1.0, color="0.35", ls=":", lw=1.2, label="Sham reference")
                if logarithmic:
                    axis.set_yscale("log")
                values = np.r_[
                    subset["observed_fold_of_sham"],
                    subset["predicted_fold_of_sham"],
                    1.0,
                ]
                if logarithmic:
                    axis.set_ylim(values.min() / 1.18, values.max() * 1.18)
                else:
                    span = max(float(values.max() - values.min()), 0.2)
                    axis.set_ylim(max(0.0, float(values.min()) - 0.12 * span), float(values.max()) + 0.12 * span)
                axis.set_title(exposure.condition)
                axis.grid(alpha=0.25)
                if row_index == 1:
                    axis.set_xlabel("Time after nsPEF (h)")
                if column_index == 0:
                    axis.set_ylabel(
                        f"Fold of {CONTROL_LABELS[control]}\nconcentration"
                    )
        axes[0, 0].legend(frameon=False, fontsize=7)
        scale = "logarithmic" if logarithmic else "linear"
        size_label = (
            "all reported sizes"
            if max_particle_diameter_nm is None
            else f"<{max_particle_diameter_nm:g} nm"
        )
        fig.suptitle(
            f"Normalized experimental and fitted EV behavior ({size_label}; {scale} scale)"
        )
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=220, bbox_inches="tight")
        plt.close(fig)


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    observation_bridge = (
        None
        if args.experimental_bridge_config is None
        else ExperimentalObservationBridge.from_yaml(args.experimental_bridge_config)
    )
    observations, controls = load_normalized_observations(
        args.data_dir,
        max_particle_diameter_nm=args.max_particle_diameter_nm,
        observation_bridge=observation_bridge,
    )
    observations.to_csv(args.output_dir / "normalized_observations.csv", index=False)
    controls.to_csv(args.output_dir / "sham_reference_summary.csv", index=False)

    comparison_frames = []
    parameter_rows = []
    summaries: dict[str, object] = {}
    for control_index, control in enumerate(CONTROL_ORDER):
        control_observations = observations.loc[
            observations["control"] == control
        ].reset_index(drop=True)
        predictor = NormalizedModelPredictor(
            control_observations,
            pulse_width_ns=args.pulse_width_ns,
            repetition_rate_hz=args.repetition_rate_hz,
        )
        default_prediction = predictor.predict(predictor.defaults)
        best, candidates = fit_normalized(
            predictor,
            starts=args.starts,
            max_nfev=args.max_nfev,
            seed=args.seed + control_index,
        )
        fitted_parameters = predictor.values_from_log_multipliers(best.x)
        fitted_prediction = predictor.predict(fitted_parameters)
        observed = control_observations["observed_fold_of_sham"].to_numpy(
            dtype=float
        )
        comparison = control_observations.copy()
        comparison["predicted_default_fold_of_sham"] = default_prediction
        comparison["predicted_fold_of_sham"] = fitted_prediction
        comparison["log10_residual"] = np.log10(fitted_prediction) - np.log10(
            observed
        )
        comparison_frames.append(comparison)

        for name in ("baseline_sEV_rate", *FREE_PARAMETER_NAMES):
            default = predictor.defaults[name]
            parameter_rows.append(
                {
                    "control_reference": control,
                    "parameter": name,
                    "default": default,
                    "fitted": fitted_parameters[name],
                    "multiplier": fitted_parameters[name] / default,
                    "fit_status": (
                        "fixed as ratio-scale reference"
                        if name == "baseline_sEV_rate"
                        else "fitted"
                    ),
                }
            )
        summaries[control] = {
            "control_label": CONTROL_LABELS[control],
            "control_has_harvest_time": False,
            "control_reference_assumption": (
                "The measured control mean is treated as a time-invariant reference at "
                "0.5, 1, and 3 hours."
            ),
            "optimizer_success": bool(best.success),
            "optimizer_message": str(best.message),
            "starts": args.starts,
            "max_nfev_per_start": args.max_nfev,
            "candidate_costs": [float(candidate.cost) for candidate in candidates],
            "metrics": ratio_metrics(observed, fitted_prediction),
            "simulator_runs": predictor.simulator_runs,
        }

    comparisons = pd.concat(comparison_frames, ignore_index=True)
    comparisons.to_csv(
        args.output_dir / "normalized_observed_vs_predicted.csv", index=False
    )
    pd.DataFrame(parameter_rows).to_csv(
        args.output_dir / "normalized_fitted_parameters.csv", index=False
    )
    plot_normalized_fits(
        comparisons,
        args.output_dir,
        max_particle_diameter_nm=args.max_particle_diameter_nm,
    )

    summary = {
        "experimental_particle_filter": (
            "all reported sizes"
            if args.max_particle_diameter_nm is None
            else f"diameter strictly less than {args.max_particle_diameter_nm:g} nm"
        ),
        "normalization": "stimulated concentration divided by one control-type mean",
        "experimental_observation_bridge": (
            None
            if observation_bridge is None
            else observation_bridge.to_metadata()
        ),
        "controls_fitted_separately": list(CONTROL_ORDER),
        "model_sham_proxy": {
            "amplitude_kV_cm": 1.0e-9,
            "pulse_number": 1,
            "reason": (
                "The current scenario schema requires positive amplitude and at least one pulse; "
                "this field is numerically negligible."
            ),
        },
        "identifiability_note": (
            "Stimulated/sham ratios cannot identify the common absolute scale of baseline EV "
            "release. baseline_sEV_rate was fixed as the reference while the relative mlEV and "
            "AB baseline rates and three kinetic constants were fitted."
        ),
        "fits": summaries,
    }
    (args.output_dir / "normalized_fit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    size_description = (
        "all reported particle sizes"
        if args.max_particle_diameter_nm is None
        else f"particles strictly below {args.max_particle_diameter_nm:g} nm"
    )
    bridge_description = (
        "No experimental-to-cell bridge was applied."
        if observation_bridge is None
        else (
            f"Concentrations were converted using {observation_bridge.initial_cell_count:g} "
            f"initial cells in {observation_bridge.medium_volume_ml:g} mL on an "
            f"{observation_bridge.cell_basis}-cell basis."
        )
    )
    readme = f"""# FFRi normalized EV kinetic fits

This experiment fits stimulated-to-sham ratios using Exoid {size_description}.
{bridge_description} `Sham2` and
`sham media` are analyzed separately because the file does not define their
protocol relationship.

Neither control has a harvest-time label. Each measured control mean is
therefore used as a constant reference at 0.5, 1, and 3 hours. This is a
sensitivity analysis, not a validated time-matched sham comparison.

The existing model source is unchanged. A near-zero-field simulation is used
as its unpulsed sham proxy. Because normalization removes absolute scale,
`baseline_sEV_rate` is fixed as a reference and only relative pathway mixture
and kinetic parameters are fitted.

- Main figure: `normalized_observed_vs_fitted_linear.png`
- Log figure: `normalized_observed_vs_fitted_log.png`
- Comparison table: `normalized_observed_vs_predicted.csv`
- Fit metrics and assumptions: `normalized_fit_summary.json`
"""
    (args.output_dir / "README.md").write_text(readme, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--max-particle-diameter-nm", type=float, default=200.0)
    parser.add_argument(
        "--experimental-bridge-config",
        type=Path,
        default=None,
        help="YAML configuration for particles/mL to particles/cell conversion",
    )
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--max-nfev", type=int, default=60)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--pulse-width-ns", type=float, default=60.0)
    parser.add_argument("--repetition-rate-hz", type=float, default=1.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.starts < 1 or args.max_nfev < 1:
        raise ValueError("--starts and --max-nfev must be positive")
    run(args)
    print(f"Normalized fits written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
