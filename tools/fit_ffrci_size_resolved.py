#!/usr/bin/env python3
"""Fit the v1.1 size-resolved extracellular bridge to FFRCI Exoid data.

The fit keeps the intracellular ODE model unchanged.  Its three cumulative
release pathways are mapped through state-conditioned lognormal size kernels,
propagated as extracellular size-bin stocks with smooth size-dependent loss,
and passed through an explicit observation matrix. Repeated measurements are
rebinned, summarized as batch means with sample SD, and fitted with total
concentration and conditional size composition as distinct information.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

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

from electro_exocytosis.experimental_bridge import (  # noqa: E402
    ExperimentalObservationBridge,
    aggregate_repeated_observations,
)
from electro_exocytosis.models.ev_size_observation import (  # noqa: E402
    PathwaySizeKernelParams,
    SizeResolvedKineticsParams,
    simulate_size_resolved_extracellular_kinetics,
)
from electro_exocytosis.simulation import Simulation  # noqa: E402
from electro_exocytosis.visualization.style import (  # noqa: E402
    FITTED_COLOR,
    FITTED_MODEL_LABEL,
    MANUSCRIPT_COLOR_DPI,
    OBSERVED_COLOR,
    OBSERVED_MEAN_LABEL,
    OBSERVED_MEAN_SD_LABEL,
    add_figure_note,
    manuscript_style_context,
    place_manuscript_legend,
    save_manuscript_figure,
    style_manuscript_axis,
)
from tools.analyze_ffrci_data import parse_exoid_file, parse_sample_label  # noqa: E402
from tools.fit_ffrci_ev_kinetics import (  # noqa: E402
    DEFAULT_DATA_DIR,
    EXOID_FILENAME,
    EXPOSURES,
    Exposure,
    build_scenario,
)


DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "ffrci_size_resolved_fit_v1_1"
DEFAULT_BRIDGE_CONFIG = REPO_ROOT / "examples" / "ffrci_experimental_bridge.yml"
COMMON_SIZE_BIN_EDGES_NM = np.arange(80.0, 380.0 + 20.0, 20.0)
LATENT_SIZE_BIN_EDGES_NM = np.r_[
    40.0,
    60.0,
    COMMON_SIZE_BIN_EDGES_NM,
    500.0,
    650.0,
    900.0,
    1300.0,
    2000.0,
]
PATHWAYS = ("sEV", "mlEV", "AB")
OBSERVATION_TIMES_H = (0.5, 1.0, 3.0)
CONDITION_LABELS = {
    "1p20kV": "1 pulse, 20 kV",
    "3p40kV": "3 pulses, 40 kV",
    "5p40kV": "5 pulses, 40 kV",
}


@dataclass(frozen=True)
class UpstreamTrajectory:
    condition: str
    dose_index: float
    time_s: np.ndarray
    cumulative_release: dict[str, np.ndarray]
    viable_fraction: np.ndarray
    state_signals: dict[str, np.ndarray]


@dataclass(frozen=True)
class FitVariant:
    name: str
    fit_state_shifts: bool
    fit_dose_response: bool


FIT_VARIANTS = (
    FitVariant("static_kernel", False, False),
    FitVariant("state_conditioned", True, False),
    FitVariant("state_conditioned_dose_response", True, True),
)


@dataclass(frozen=True)
class ParameterDefinition:
    name: str
    initial: float
    lower: float
    upper: float
    prior_center: float | None = None
    prior_scale: float | None = None


BASE_PARAMETER_DEFINITIONS = (
    ParameterDefinition(
        "log_sEV_source_scale", math.log(1.0e4), math.log(1.0e-2), math.log(1.0e10)
    ),
    ParameterDefinition(
        "log_mlEV_source_scale", math.log(1.0e2), math.log(1.0e-2), math.log(1.0e8)
    ),
    ParameterDefinition(
        "log_AB_source_scale", math.log(1.0e4), math.log(1.0e-2), math.log(1.0e10)
    ),
    ParameterDefinition(
        "log_effective_half_life_h", math.log(2.0), math.log(0.10), math.log(48.0)
    ),
    ParameterDefinition("loss_size_exponent", 0.0, -3.0, 3.0, 0.0, 1.0),
    ParameterDefinition(
        "log_sEV_median_nm",
        math.log(105.0),
        math.log(50.0),
        math.log(190.0),
        math.log(105.0),
        math.log(1.35),
    ),
    ParameterDefinition(
        "log_mlEV_median_nm",
        math.log(210.0),
        math.log(100.0),
        math.log(420.0),
        math.log(210.0),
        math.log(1.45),
    ),
    ParameterDefinition(
        "log_AB_median_nm",
        math.log(420.0),
        math.log(180.0),
        math.log(900.0),
        math.log(420.0),
        math.log(1.60),
    ),
    ParameterDefinition(
        "log_sEV_gsd_minus_one",
        math.log(0.35),
        math.log(0.08),
        math.log(1.5),
        math.log(0.35),
        0.65,
    ),
    ParameterDefinition(
        "log_mlEV_gsd_minus_one",
        math.log(0.45),
        math.log(0.08),
        math.log(1.5),
        math.log(0.45),
        0.65,
    ),
    ParameterDefinition(
        "log_AB_gsd_minus_one",
        math.log(0.60),
        math.log(0.08),
        math.log(1.5),
        math.log(0.60),
        0.65,
    ),
)

STATE_PARAMETER_DEFINITIONS = tuple(
    ParameterDefinition(f"{pathway}_state_shift", 0.0, -4.0, 4.0, 0.0, 1.0)
    for pathway in PATHWAYS
)

DOSE_PARAMETER_DEFINITIONS = (
    ParameterDefinition("dose_response_linear", 0.0, -8.0, 8.0, 0.0, 2.0),
    ParameterDefinition("dose_response_quadratic", 0.0, -8.0, 8.0, 0.0, 2.0),
)


def parameter_definitions(variant: FitVariant) -> tuple[ParameterDefinition, ...]:
    definitions = BASE_PARAMETER_DEFINITIONS
    if variant.fit_state_shifts:
        definitions += STATE_PARAMETER_DEFINITIONS
    if variant.fit_dose_response:
        definitions += DOSE_PARAMETER_DEFINITIONS
    return definitions


def load_size_resolved_observations(data_dir: Path) -> pd.DataFrame:
    """Rebin treatment and cell-containing controls without pairing them."""

    raw = parse_exoid_file(data_dir / EXOID_FILENAME)
    concentration = raw["concentration_particles_per_ml"].to_numpy(dtype=float)
    if np.any(~np.isfinite(concentration)) or np.any(concentration < 0.0):
        raise ValueError("Particle concentrations must be finite and nonnegative")

    metadata_rows: list[dict[str, object]] = []
    for measurement_id, group in raw.groupby("dataset_number", sort=True):
        parsed = parse_sample_label(str(group["label"].iloc[0]))
        sample_type = str(parsed["sample_type"])
        condition = str(parsed["condition"])
        if sample_type == "treatment":
            sample_role = "treatment"
        elif condition == "sham2":
            sample_role = "initial_control"
        else:
            continue
        metadata_rows.append(
            {
                "measurement_id": int(measurement_id),
                "sample_role": sample_role,
                "condition": condition,
                "time_h": parsed["harvest_time_h"],
            }
        )
    metadata = pd.DataFrame(metadata_rows)
    if metadata.empty:
        raise ValueError("No treatment or cell-containing control measurements found")

    retained = raw[raw["dataset_number"].isin(metadata["measurement_id"])].copy()
    retained["size_bin_center_nm"] = pd.cut(
        retained["particle_diameter_nm"],
        bins=COMMON_SIZE_BIN_EDGES_NM,
        labels=(COMMON_SIZE_BIN_EDGES_NM[:-1] + COMMON_SIZE_BIN_EDGES_NM[1:]) / 2.0,
        right=False,
    ).astype(float)
    binned = (
        retained.dropna(subset=["size_bin_center_nm"])
        .groupby(["dataset_number", "size_bin_center_nm"], observed=True, sort=True)[
            "concentration_particles_per_ml"
        ]
        .sum()
    )
    complete_index = pd.MultiIndex.from_product(
        [
            metadata["measurement_id"].to_numpy(dtype=int),
            (COMMON_SIZE_BIN_EDGES_NM[:-1] + COMMON_SIZE_BIN_EDGES_NM[1:]) / 2.0,
        ],
        names=["dataset_number", "size_bin_center_nm"],
    )
    observations = (
        binned.reindex(complete_index, fill_value=0.0)
        .rename("concentration_particles_per_ml")
        .reset_index()
        .rename(columns={"dataset_number": "measurement_id"})
        .merge(metadata, on="measurement_id", how="left", validate="many_to_one")
    )
    observations["size_bin_lower_nm"] = observations["size_bin_center_nm"] - 10.0
    observations["size_bin_upper_nm"] = observations["size_bin_center_nm"] + 10.0
    observations = observations[
        [
            "sample_role",
            "condition",
            "time_h",
            "measurement_id",
            "size_bin_lower_nm",
            "size_bin_upper_nm",
            "size_bin_center_nm",
            "concentration_particles_per_ml",
        ]
    ].sort_values(
        ["sample_role", "condition", "time_h", "measurement_id", "size_bin_center_nm"],
        na_position="first",
    )
    observations = observations.reset_index(drop=True)
    expected_bins = len(COMMON_SIZE_BIN_EDGES_NM) - 1
    group_sizes = observations.groupby("measurement_id").size()
    if not group_sizes.eq(expected_bins).all():
        raise ValueError(
            "Every retained measurement must contain every common size bin"
        )
    if not np.allclose(
        np.sort(observations["size_bin_lower_nm"].unique()),
        COMMON_SIZE_BIN_EDGES_NM[:-1],
    ):
        raise ValueError("Unexpected lower edges in the common size rebinning")
    treatment = observations[observations["sample_role"] == "treatment"]
    if set(treatment["condition"].unique()) != {
        exposure.condition for exposure in EXPOSURES
    }:
        raise ValueError("Treatment conditions do not match the configured exposures")
    if not np.all(treatment.groupby("measurement_id")["time_h"].nunique().eq(1)):
        raise ValueError("Each treatment measurement must have one harvest time")
    if not np.all(
        observations.groupby("measurement_id")["concentration_particles_per_ml"].sum()
        > 0.0
    ):
        raise ValueError("Every retained measurement must have positive total mass")
    return observations


def aggregate_size_resolved_observations(
    raw_observations: pd.DataFrame,
) -> pd.DataFrame:
    """Build one model-facing mean distribution per condition and harvest time.

    Raw repeated histograms remain available to the caller. The returned table
    contains the mean, sample SD, SE, and count for every 20-nm concentration
    band, an all-control mean initial distribution, and independently computed
    total-concentration statistics that retain covariance across size bands.
    """

    required = {
        "sample_role",
        "condition",
        "time_h",
        "measurement_id",
        "size_bin_lower_nm",
        "size_bin_upper_nm",
        "size_bin_center_nm",
        "concentration_particles_per_ml",
    }
    missing = sorted(required - set(raw_observations.columns))
    if missing:
        raise ValueError(f"Missing size-resolved observation columns: {missing}")
    raw_concentration = raw_observations["concentration_particles_per_ml"].to_numpy(
        dtype=float
    )
    if np.any(~np.isfinite(raw_concentration)) or np.any(raw_concentration < 0.0):
        raise ValueError("Raw concentrations must be finite and nonnegative")

    treatment_rows = raw_observations[
        raw_observations["sample_role"] == "treatment"
    ].copy()
    control_rows = raw_observations[
        raw_observations["sample_role"] == "initial_control"
    ].copy()
    if treatment_rows.empty or control_rows.empty:
        raise ValueError("Both treatment and initial-control measurements are required")

    bin_summary = aggregate_repeated_observations(
        treatment_rows,
        group_columns=(
            "condition",
            "time_h",
            "size_bin_lower_nm",
            "size_bin_upper_nm",
            "size_bin_center_nm",
        ),
        value_columns=("concentration_particles_per_ml",),
        sort_groups=True,
    ).summary.rename(
        columns={
            "concentration_particles_per_ml_mean": "observed_particles_per_ml_mean",
            "concentration_particles_per_ml_sd": "observed_particles_per_ml_sd",
            "concentration_particles_per_ml_se": "observed_particles_per_ml_se",
            "concentration_particles_per_ml_n": "observed_measurement_count",
        }
    )

    control_summary = aggregate_repeated_observations(
        control_rows,
        group_columns=(
            "size_bin_lower_nm",
            "size_bin_upper_nm",
            "size_bin_center_nm",
        ),
        value_columns=("concentration_particles_per_ml",),
        sort_groups=True,
    ).summary.rename(
        columns={
            "concentration_particles_per_ml_mean": "initial_control_particles_per_ml_mean",
            "concentration_particles_per_ml_sd": "initial_control_particles_per_ml_sd",
            "concentration_particles_per_ml_se": "initial_control_particles_per_ml_se",
            "concentration_particles_per_ml_n": "initial_control_measurement_count",
        }
    )

    total_rows = (
        treatment_rows.groupby(
            ["condition", "time_h", "measurement_id"],
            observed=True,
            sort=True,
        )["concentration_particles_per_ml"]
        .sum()
        .rename("observed_total_particles_per_ml")
        .reset_index()
    )
    total_summary = aggregate_repeated_observations(
        total_rows,
        group_columns=("condition", "time_h"),
        value_columns=("observed_total_particles_per_ml",),
        sort_groups=True,
    ).summary.rename(
        columns={
            "observed_total_particles_per_ml_n": "total_measurement_count",
        }
    )

    summary = bin_summary.merge(
        control_summary,
        on=[
            "size_bin_lower_nm",
            "size_bin_upper_nm",
            "size_bin_center_nm",
        ],
        how="left",
        validate="many_to_one",
    ).merge(
        total_summary,
        on=["condition", "time_h"],
        how="left",
        validate="many_to_one",
    )
    summary = summary.sort_values(
        ["condition", "time_h", "size_bin_center_nm"]
    ).reset_index(drop=True)
    expected_rows = (
        len(EXPOSURES) * len(OBSERVATION_TIMES_H) * (len(COMMON_SIZE_BIN_EDGES_NM) - 1)
    )
    if len(summary) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} batch-mean size-bin observations, got {len(summary)}"
        )
    mean_columns = (
        "observed_particles_per_ml_mean",
        "initial_control_particles_per_ml_mean",
        "observed_total_particles_per_ml_mean",
    )
    count_columns = (
        "observed_measurement_count",
        "initial_control_measurement_count",
        "total_measurement_count",
    )
    means = summary.loc[:, mean_columns].to_numpy(dtype=float)
    counts = summary.loc[:, count_columns].to_numpy(dtype=float)
    if np.any(~np.isfinite(means)) or np.any(means < 0.0):
        raise ValueError(
            "Aggregated concentration means must be finite and nonnegative"
        )
    if np.any(~np.isfinite(counts)) or np.any(counts < 1.0):
        raise ValueError("Aggregated measurement counts must be finite and positive")
    uncertainty_groups = (
        (
            "observed_particles_per_ml_sd",
            "observed_particles_per_ml_se",
            "observed_measurement_count",
        ),
        (
            "initial_control_particles_per_ml_sd",
            "initial_control_particles_per_ml_se",
            "initial_control_measurement_count",
        ),
        (
            "observed_total_particles_per_ml_sd",
            "observed_total_particles_per_ml_se",
            "total_measurement_count",
        ),
    )
    for sd_column, se_column, count_column in uncertainty_groups:
        sd = summary[sd_column].to_numpy(dtype=float)
        se = summary[se_column].to_numpy(dtype=float)
        n = summary[count_column].to_numpy(dtype=float)
        invalid = (n >= 2.0) & (
            ~np.isfinite(sd) | ~np.isfinite(se) | (sd < 0.0) | (se < 0.0)
        )
        if np.any(invalid):
            raise ValueError(
                f"{sd_column} and {se_column} must be finite for groups with n >= 2"
            )
    if np.any(
        summary.groupby(["condition", "time_h"])["observed_particles_per_ml_mean"]
        .sum()
        .to_numpy(dtype=float)
        <= 0.0
    ):
        raise ValueError("Every mean treatment distribution must have positive mass")
    if float(control_summary["initial_control_particles_per_ml_mean"].sum()) <= 0.0:
        raise ValueError(
            "The mean initial-control distribution must have positive mass"
        )
    return summary


def _scenario_for_bridge(
    exposure: Exposure,
    bridge: ExperimentalObservationBridge,
    *,
    pulse_width_ns: float,
    repetition_rate_hz: float,
) -> object:
    scenario = build_scenario(
        exposure,
        pulse_width_ns=pulse_width_ns,
        repetition_rate_hz=repetition_rate_hz,
    )
    scenario.exposure.cell_density_per_ml = (
        bridge.initial_cell_count / bridge.medium_volume_ml
    )
    scenario.extracellular_medium.initial_volume_ml = bridge.medium_volume_ml
    scenario.extracellular_medium.use_time_varying_viability = True
    return scenario


def build_upstream_trajectories(
    bridge: ExperimentalObservationBridge,
    *,
    pulse_width_ns: float,
    repetition_rate_hz: float,
) -> dict[str, UpstreamTrajectory]:
    """Run the unchanged intracellular model once per exposure and one sham proxy."""

    sham_exposure = Exposure("sham_proxy", 1, 1.0e-9)
    sham_result = Simulation(
        _scenario_for_bridge(
            sham_exposure,
            bridge,
            pulse_width_ns=pulse_width_ns,
            repetition_rate_hz=repetition_rate_hz,
        )
    ).run()
    sham_frame = sham_result.ev_timeseries

    trajectories: dict[str, UpstreamTrajectory] = {}
    for exposure in EXPOSURES:
        result = Simulation(
            _scenario_for_bridge(
                exposure,
                bridge,
                pulse_width_ns=pulse_width_ns,
                repetition_rate_hz=repetition_rate_hz,
            )
        ).run()
        frame = result.ev_timeseries
        if not np.allclose(frame["t"], sham_frame["t"]):
            raise ValueError(
                "Treatment and sham simulations produced different time grids"
            )
        state_signals = {
            "sEV": 0.5
            * (
                frame["escrt_dependent_signal"].to_numpy(dtype=float)
                + frame["ceramide_signal"].to_numpy(dtype=float)
                - sham_frame["escrt_dependent_signal"].to_numpy(dtype=float)
                - sham_frame["ceramide_signal"].to_numpy(dtype=float)
            ),
            "mlEV": 0.5
            * (
                frame["budding_signal"].to_numpy(dtype=float)
                + frame["scission_signal"].to_numpy(dtype=float)
                - sham_frame["budding_signal"].to_numpy(dtype=float)
                - sham_frame["scission_signal"].to_numpy(dtype=float)
            ),
            "AB": (
                frame["apoptotic_blebbing_signal"].to_numpy(dtype=float)
                - sham_frame["apoptotic_blebbing_signal"].to_numpy(dtype=float)
            ),
        }
        trajectories[exposure.condition] = UpstreamTrajectory(
            condition=exposure.condition,
            dose_index=float(result.summary["dose_index"]),
            time_s=frame["t"].to_numpy(dtype=float),
            cumulative_release={
                "sEV": frame["sEV_cumulative"].to_numpy(dtype=float),
                "mlEV": frame["mlEV_cumulative"].to_numpy(dtype=float),
                "AB": frame["AB_cumulative"].to_numpy(dtype=float),
            },
            viable_fraction=frame["viable_producer_fraction"].to_numpy(dtype=float),
            state_signals=state_signals,
        )
    return trajectories


def decode_parameters(
    vector: np.ndarray,
    variant: FitVariant,
) -> dict[str, float]:
    definitions = parameter_definitions(variant)
    if len(vector) != len(definitions):
        raise ValueError("Parameter vector length does not match the fit variant")
    raw = {
        definition.name: float(value)
        for definition, value in zip(definitions, vector, strict=True)
    }
    decoded = {
        "sEV_source_scale_particles_per_model_unit": math.exp(
            raw["log_sEV_source_scale"]
        ),
        "mlEV_source_scale_particles_per_model_unit": math.exp(
            raw["log_mlEV_source_scale"]
        ),
        "AB_source_scale_particles_per_model_unit": math.exp(
            raw["log_AB_source_scale"]
        ),
        "effective_half_life_h": math.exp(raw["log_effective_half_life_h"]),
        "loss_size_exponent": raw["loss_size_exponent"],
        "sEV_median_diameter_nm": math.exp(raw["log_sEV_median_nm"]),
        "mlEV_median_diameter_nm": math.exp(raw["log_mlEV_median_nm"]),
        "AB_median_diameter_nm": math.exp(raw["log_AB_median_nm"]),
        "sEV_geometric_sd": 1.0 + math.exp(raw["log_sEV_gsd_minus_one"]),
        "mlEV_geometric_sd": 1.0 + math.exp(raw["log_mlEV_gsd_minus_one"]),
        "AB_geometric_sd": 1.0 + math.exp(raw["log_AB_gsd_minus_one"]),
        "sEV_state_shift": raw.get("sEV_state_shift", 0.0),
        "mlEV_state_shift": raw.get("mlEV_state_shift", 0.0),
        "AB_state_shift": raw.get("AB_state_shift", 0.0),
        "dose_response_linear": raw.get("dose_response_linear", 0.0),
        "dose_response_quadratic": raw.get("dose_response_quadratic", 0.0),
    }
    return decoded


def _pathway_parameters(values: dict[str, float]) -> dict[str, PathwaySizeKernelParams]:
    return {
        pathway: PathwaySizeKernelParams(
            median_diameter_nm=values[f"{pathway}_median_diameter_nm"],
            geometric_sd=values[f"{pathway}_geometric_sd"],
            source_scale_particles_per_model_unit=values[
                f"{pathway}_source_scale_particles_per_model_unit"
            ],
            state_shift_coefficient=values[f"{pathway}_state_shift"],
        )
        for pathway in PATHWAYS
    }


def _kinetics_parameters(values: dict[str, float]) -> SizeResolvedKineticsParams:
    return SizeResolvedKineticsParams(
        effective_loss_rate_s=math.log(2.0)
        / (values["effective_half_life_h"] * 3600.0),
        loss_reference_diameter_nm=150.0,
        loss_size_exponent=values["loss_size_exponent"],
        instrument_log_diameter_sd=0.0,
        assay_recovery_fraction=1.0,
    )


def _common_to_latent_initial(common_concentration: np.ndarray) -> np.ndarray:
    common = np.asarray(common_concentration, dtype=float)
    if common.shape != (len(COMMON_SIZE_BIN_EDGES_NM) - 1,):
        raise ValueError("Common initial concentration has the wrong number of bins")
    latent = np.zeros(len(LATENT_SIZE_BIN_EDGES_NM) - 1, dtype=float)
    for common_index, (lower, upper) in enumerate(
        zip(
            COMMON_SIZE_BIN_EDGES_NM[:-1],
            COMMON_SIZE_BIN_EDGES_NM[1:],
            strict=True,
        )
    ):
        latent_index = np.flatnonzero(
            np.isclose(LATENT_SIZE_BIN_EDGES_NM[:-1], lower)
            & np.isclose(LATENT_SIZE_BIN_EDGES_NM[1:], upper)
        )
        if len(latent_index) != 1:
            raise ValueError("Common observation bins must be exact latent-bin subsets")
        latent[int(latent_index[0])] = common[common_index]
    return latent


def _dose_gain(values: dict[str, float], dose_index: float) -> float:
    centered_dose = float(dose_index) - 0.72
    exponent = (
        values["dose_response_linear"] * centered_dose
        + values["dose_response_quadratic"] * centered_dose**2
    )
    return float(math.exp(np.clip(exponent, -20.0, 20.0)))


def hellinger_distance(observed: np.ndarray, predicted: np.ndarray) -> float:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if observed.shape != predicted.shape or observed.ndim != 1:
        raise ValueError(
            "Hellinger inputs must be one-dimensional arrays of equal length"
        )
    if np.any(observed < 0.0) or np.any(predicted < 0.0):
        raise ValueError("Hellinger inputs must be nonnegative")
    observed_total = float(observed.sum())
    predicted_total = float(predicted.sum())
    if observed_total <= 0.0 or predicted_total <= 0.0:
        raise ValueError("Hellinger inputs must each have positive mass")
    observed_fraction = observed / observed_total
    predicted_fraction = predicted / predicted_total
    return float(
        np.sqrt(
            0.5
            * np.sum((np.sqrt(observed_fraction) - np.sqrt(predicted_fraction)) ** 2)
        )
    )


class SizeResolvedFFRCIPredictor:
    """Size/loss/observation fit to batch-mean extracellular measurements."""

    def __init__(
        self,
        observations: pd.DataFrame,
        trajectories: dict[str, UpstreamTrajectory],
        bridge: ExperimentalObservationBridge,
        variant: FitVariant,
    ) -> None:
        self.observations = observations.copy().reset_index(drop=True)
        self.trajectories = trajectories
        self.bridge = bridge
        self.variant = variant
        self.definitions = parameter_definitions(variant)
        required = {
            "condition",
            "time_h",
            "size_bin_center_nm",
            "observed_particles_per_ml_mean",
            "observed_particles_per_ml_sd",
            "initial_control_particles_per_ml_mean",
            "observed_total_particles_per_ml_mean",
        }
        missing = sorted(required - set(self.observations.columns))
        if missing:
            raise ValueError(f"Missing batch-mean observation columns: {missing}")
        self.cell_keys = list(
            self.observations[["condition", "time_h"]]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )
        self._row_index = {
            key: self.observations.index[
                (self.observations["condition"] == key[0])
                & (self.observations["time_h"] == key[1])
            ].to_numpy()
            for key in self.cell_keys
        }
        self.evaluations = 0

    @property
    def initial_vector(self) -> np.ndarray:
        return np.asarray(
            [definition.initial for definition in self.definitions], dtype=float
        )

    @property
    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        return (
            np.asarray(
                [definition.lower for definition in self.definitions], dtype=float
            ),
            np.asarray(
                [definition.upper for definition in self.definitions], dtype=float
            ),
        )

    def predict(self, vector: np.ndarray) -> pd.DataFrame:
        values = decode_parameters(vector, self.variant)
        pathway_params = _pathway_parameters(values)
        kinetics_params = _kinetics_parameters(values)
        prediction = self.observations.copy()
        prediction["predicted_particles_per_ml"] = np.nan
        prediction["predicted_true_particles_per_ml"] = np.nan
        prediction["dose_response_gain"] = np.nan
        initial_by_bin = (
            prediction.drop_duplicates("size_bin_center_nm")
            .sort_values("size_bin_center_nm")["initial_control_particles_per_ml_mean"]
            .to_numpy(dtype=float)
        )
        initial_latent = _common_to_latent_initial(initial_by_bin)

        for condition in prediction["condition"].unique():
            trajectory = self.trajectories[str(condition)]
            gain = _dose_gain(values, trajectory.dose_index)
            scaled_cumulative = {
                pathway: values_array * gain
                for pathway, values_array in trajectory.cumulative_release.items()
            }
            source_result = simulate_size_resolved_extracellular_kinetics(
                trajectory.time_s,
                scaled_cumulative,
                LATENT_SIZE_BIN_EDGES_NM,
                pathway_params,
                kinetics_params,
                initial_cell_density_per_ml=(
                    self.bridge.initial_cell_count / self.bridge.medium_volume_ml
                ),
                viable_fraction=trajectory.viable_fraction,
                state_signals=trajectory.state_signals,
                observed_bin_edges_nm=COMMON_SIZE_BIN_EDGES_NM,
            )
            source_true_common = (
                source_result.total_latent_concentration_particles_per_ml
                @ source_result.observation_matrix.T
            )
            ambient_latent = (
                np.exp(
                    -trajectory.time_s[:, None] * source_result.loss_rates_s[None, :]
                )
                * initial_latent[None, :]
            )
            ambient_common = ambient_latent @ source_result.observation_matrix.T
            observed_common = (
                source_result.observed_concentration_particles_per_ml
                + ambient_common * kinetics_params.assay_recovery_fraction
            )
            true_common = source_true_common + ambient_common
            group = prediction[prediction["condition"] == condition]
            for time_h in sorted(group["time_h"].unique()):
                row_mask = (prediction["condition"] == condition) & (
                    prediction["time_h"] == time_h
                )
                row_indices = prediction.index[row_mask].to_numpy()
                row_indices = row_indices[
                    np.argsort(
                        prediction.loc[row_indices, "size_bin_center_nm"].to_numpy()
                    )
                ]
                time_index = int(
                    np.argmin(np.abs(trajectory.time_s - float(time_h) * 3600.0))
                )
                prediction.loc[row_indices, "predicted_particles_per_ml"] = (
                    observed_common[time_index]
                )
                prediction.loc[row_indices, "predicted_true_particles_per_ml"] = (
                    true_common[time_index]
                )
                prediction.loc[row_indices, "dose_response_gain"] = gain

        if prediction["predicted_particles_per_ml"].isna().any():
            raise RuntimeError("Some size-resolved observations were not predicted")
        density = self.bridge.initial_cell_count / self.bridge.medium_volume_ml
        prediction["observed_particle_equivalents_per_initial_cell"] = (
            prediction["observed_particles_per_ml_mean"] / density
        )
        prediction["observed_particle_equivalents_per_initial_cell_sd"] = (
            prediction["observed_particles_per_ml_sd"] / density
        )
        prediction["observed_particle_equivalents_per_initial_cell_se"] = (
            prediction["observed_particles_per_ml_se"] / density
        )
        prediction["predicted_particle_equivalents_per_initial_cell"] = (
            prediction["predicted_particles_per_ml"] / density
        )
        self.evaluations += 1
        return prediction

    def data_residuals(self, prediction: pd.DataFrame) -> np.ndarray:
        residuals: list[float] = []
        for condition, time_h in self.cell_keys:
            cell = prediction[
                (prediction["condition"] == condition)
                & (prediction["time_h"] == time_h)
            ]
            observed_total = float(cell["observed_total_particles_per_ml_mean"].iloc[0])
            predicted_total = float(cell["predicted_particles_per_ml"].sum())
            residuals.append(float(np.log(predicted_total / observed_total)))

        for key in self.cell_keys:
            group = prediction.loc[self._row_index[key]].sort_values(
                "size_bin_center_nm"
            )
            observed = group["observed_particles_per_ml_mean"].to_numpy(dtype=float)
            predicted = group["predicted_particles_per_ml"].to_numpy(dtype=float)
            observed_fraction = observed / observed.sum()
            predicted_fraction = predicted / predicted.sum()
            residuals.extend(
                (np.sqrt(predicted_fraction) - np.sqrt(observed_fraction)).tolist()
            )
        return np.asarray(residuals, dtype=float)

    def prior_residuals(self, vector: np.ndarray) -> np.ndarray:
        return np.asarray(
            [
                (float(value) - float(definition.prior_center))
                / float(definition.prior_scale)
                for definition, value in zip(self.definitions, vector, strict=True)
                if definition.prior_center is not None
                and definition.prior_scale is not None
            ],
            dtype=float,
        )

    def residuals(self, vector: np.ndarray) -> np.ndarray:
        try:
            prediction = self.predict(vector)
            data = self.data_residuals(prediction)
            priors = 0.25 * self.prior_residuals(vector)
            if not np.all(np.isfinite(data)):
                raise FloatingPointError("Non-finite size-resolved residual")
            return np.r_[data, priors]
        except Exception:
            residual_count = len(self.cell_keys) * (
                1 + len(COMMON_SIZE_BIN_EDGES_NM) - 1
            )
            prior_count = sum(
                definition.prior_center is not None for definition in self.definitions
            )
            return np.full(residual_count + prior_count, 1.0e6, dtype=float)


def fit_variant(
    predictor: SizeResolvedFFRCIPredictor,
    *,
    starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[object, list[object]]:
    lower, upper = predictor.bounds
    rng = np.random.default_rng(seed)
    guesses = [predictor.initial_vector]
    for _ in range(starts - 1):
        guess = predictor.initial_vector + rng.normal(0.0, 0.55, len(lower))
        guesses.append(np.clip(guess, lower + 1.0e-8, upper - 1.0e-8))
    candidates = [
        least_squares(
            predictor.residuals,
            guess,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            diff_step=0.02,
            max_nfev=max_nfev,
            ftol=1.0e-8,
            xtol=1.0e-8,
            gtol=1.0e-8,
        )
        for guess in guesses
    ]
    return min(candidates, key=lambda candidate: candidate.cost), candidates


def prediction_metrics(
    prediction: pd.DataFrame,
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame]:
    total_rows: list[dict[str, float | str]] = []
    for (condition, time_h), group in prediction.groupby(["condition", "time_h"]):
        observed = float(group["observed_total_particles_per_ml_mean"].iloc[0])
        observed_sd = float(group["observed_total_particles_per_ml_sd"].iloc[0])
        observed_se = float(group["observed_total_particles_per_ml_se"].iloc[0])
        measurement_count = int(group["total_measurement_count"].iloc[0])
        predicted = float(group["predicted_particles_per_ml"].sum())
        total_rows.append(
            {
                "condition": str(condition),
                "time_h": float(time_h),
                "measurement_count": measurement_count,
                "observed_total_particles_per_ml": observed,
                "observed_total_particles_per_ml_sd": observed_sd,
                "observed_total_particles_per_ml_se": observed_se,
                "predicted_total_particles_per_ml": predicted,
                "predicted_over_observed": predicted / observed,
                "log10_residual": math.log10(predicted / observed),
            }
        )
    total_frame = pd.DataFrame(total_rows)
    total_difference = (
        total_frame["predicted_total_particles_per_ml"]
        - total_frame["observed_total_particles_per_ml"]
    )
    total_sd = total_frame["observed_total_particles_per_ml_sd"]
    total_frame["model_within_observed_sd"] = total_difference.abs() <= total_sd
    total_frame["sd_scaled_difference"] = np.divide(
        total_difference,
        total_sd,
        out=np.full(len(total_frame), np.nan, dtype=float),
        where=total_sd.to_numpy(dtype=float) > 0.0,
    )

    centers = np.asarray(sorted(prediction["size_bin_center_nm"].unique()), dtype=float)
    distribution_rows: list[dict[str, float | str]] = []
    for (condition, time_h), group in prediction.groupby(["condition", "time_h"]):
        group = group.sort_values("size_bin_center_nm")
        observed = group["observed_particles_per_ml_mean"].to_numpy(dtype=float)
        predicted = group["predicted_particles_per_ml"].to_numpy(dtype=float)
        observed_fraction = observed / observed.sum()
        predicted_fraction = predicted / predicted.sum()
        counts = group["observed_measurement_count"].to_numpy(dtype=int)
        if np.any(counts != counts[0]):
            raise ValueError(
                "Measurement counts vary across bins within a distribution"
            )
        distribution_rows.append(
            {
                "condition": str(condition),
                "time_h": float(time_h),
                "measurement_count": int(counts[0]),
                "hellinger_distance": hellinger_distance(observed, predicted),
                "observed_mean_diameter_nm": float(np.sum(centers * observed_fraction)),
                "predicted_mean_diameter_nm": float(
                    np.sum(centers * predicted_fraction)
                ),
            }
        )
    distribution_frame = pd.DataFrame(distribution_rows)
    distribution_frame["mean_diameter_error_nm"] = (
        distribution_frame["predicted_mean_diameter_nm"]
        - distribution_frame["observed_mean_diameter_nm"]
    )

    log_total_residual = total_frame["log10_residual"].to_numpy(dtype=float)
    total_observed = total_frame["observed_total_particles_per_ml"].to_numpy(
        dtype=float
    )
    total_predicted = total_frame["predicted_total_particles_per_ml"].to_numpy(
        dtype=float
    )
    metrics = {
        "model_facing_total_targets": int(len(total_frame)),
        "batch_mean_size_distributions": int(len(distribution_frame)),
        "rmse_log10_total": float(np.sqrt(np.mean(log_total_residual**2))),
        "mae_log10_total": float(np.mean(np.abs(log_total_residual))),
        "median_absolute_percent_total_error": float(
            100.0 * np.median(np.abs(total_predicted - total_observed) / total_observed)
        ),
        "mean_hellinger_size_distance": float(
            distribution_frame["hellinger_distance"].mean()
        ),
        "median_hellinger_size_distance": float(
            distribution_frame["hellinger_distance"].median()
        ),
        "mae_mean_diameter_nm": float(
            distribution_frame["mean_diameter_error_nm"].abs().mean()
        ),
        "descriptive_fraction_totals_within_observed_sd": float(
            total_frame["model_within_observed_sd"].mean()
        ),
        "descriptive_fraction_size_bins_within_observed_sd": float(
            np.mean(
                np.abs(
                    prediction["predicted_particles_per_ml"].to_numpy(dtype=float)
                    - prediction["observed_particles_per_ml_mean"].to_numpy(dtype=float)
                )
                <= prediction["observed_particles_per_ml_sd"].to_numpy(dtype=float)
            )
        ),
    }
    return metrics, total_frame, distribution_frame


def _parameter_rows(
    predictor: SizeResolvedFFRCIPredictor,
    vector: np.ndarray,
) -> list[dict[str, float | str | bool]]:
    variant = predictor.variant
    lower, upper = predictor.bounds
    fitted = decode_parameters(vector, variant)
    initial = decode_parameters(predictor.initial_vector, variant)
    lower_values = decode_parameters(lower, variant)
    upper_values = decode_parameters(upper, variant)
    active_names = [
        "sEV_source_scale_particles_per_model_unit",
        "mlEV_source_scale_particles_per_model_unit",
        "AB_source_scale_particles_per_model_unit",
        "effective_half_life_h",
        "loss_size_exponent",
        "sEV_median_diameter_nm",
        "mlEV_median_diameter_nm",
        "AB_median_diameter_nm",
        "sEV_geometric_sd",
        "mlEV_geometric_sd",
        "AB_geometric_sd",
    ]
    if variant.fit_state_shifts:
        active_names += [f"{pathway}_state_shift" for pathway in PATHWAYS]
    if variant.fit_dose_response:
        active_names += ["dose_response_linear", "dose_response_quadratic"]
    return [
        {
            "variant": variant.name,
            "parameter": name,
            "initial": initial[name],
            "fitted": fitted[name],
            "lower_bound": lower_values[name],
            "upper_bound": upper_values[name],
            "at_lower_bound": math.isclose(
                fitted[name], lower_values[name], rel_tol=1e-4, abs_tol=1e-10
            ),
            "at_upper_bound": math.isclose(
                fitted[name], upper_values[name], rel_tol=1e-4, abs_tol=1e-10
            ),
        }
        for name in active_names
    ]


def build_kernel_frame(
    predictor: SizeResolvedFFRCIPredictor,
    vector: np.ndarray,
) -> pd.DataFrame:
    values = decode_parameters(vector, predictor.variant)
    pathway_params = _pathway_parameters(values)
    kinetics_params = _kinetics_parameters(values)
    rows: list[dict[str, float | str]] = []
    first_initial = (
        predictor.observations.drop_duplicates("size_bin_center_nm")
        .sort_values("size_bin_center_nm")["initial_control_particles_per_ml_mean"]
        .to_numpy(dtype=float)
    )
    for condition, trajectory in predictor.trajectories.items():
        gain = _dose_gain(values, trajectory.dose_index)
        result = simulate_size_resolved_extracellular_kinetics(
            trajectory.time_s,
            {
                pathway: cumulative * gain
                for pathway, cumulative in trajectory.cumulative_release.items()
            },
            LATENT_SIZE_BIN_EDGES_NM,
            pathway_params,
            kinetics_params,
            initial_cell_density_per_ml=(
                predictor.bridge.initial_cell_count / predictor.bridge.medium_volume_ml
            ),
            initial_concentration_particles_per_ml=_common_to_latent_initial(
                first_initial
            ),
            viable_fraction=trajectory.viable_fraction,
            state_signals=trajectory.state_signals,
            observed_bin_edges_nm=COMMON_SIZE_BIN_EDGES_NM,
        )
        for time_h in OBSERVATION_TIMES_H:
            time_index = int(
                np.argmin(np.abs(trajectory.time_s - float(time_h) * 3600.0))
            )
            kernel_index = max(time_index - 1, 0)
            for pathway in PATHWAYS:
                probabilities = result.kernel_probabilities[pathway][kernel_index]
                observed_window_mass = float(
                    np.sum(result.observation_matrix @ probabilities)
                )
                for bin_index, probability in enumerate(probabilities):
                    rows.append(
                        {
                            "variant": predictor.variant.name,
                            "condition": condition,
                            "time_h": float(time_h),
                            "pathway": pathway,
                            "size_bin_lower_nm": float(
                                LATENT_SIZE_BIN_EDGES_NM[bin_index]
                            ),
                            "size_bin_upper_nm": float(
                                LATENT_SIZE_BIN_EDGES_NM[bin_index + 1]
                            ),
                            "size_bin_center_nm": float(
                                0.5
                                * (
                                    LATENT_SIZE_BIN_EDGES_NM[bin_index]
                                    + LATENT_SIZE_BIN_EDGES_NM[bin_index + 1]
                                )
                            ),
                            "kernel_probability": float(probability),
                            "kernel_mass_in_latent_domain": float(probabilities.sum()),
                            "kernel_mass_in_common_window": observed_window_mass,
                            "loss_rate_s": float(result.loss_rates_s[bin_index]),
                            "dose_response_gain": gain,
                        }
                    )
    return pd.DataFrame(rows)


def _measurement_count_note(frame: pd.DataFrame, column: str) -> str:
    counts = sorted({int(value) for value in frame[column].dropna().unique()})
    if not counts:
        return "sample count unavailable"
    if len(counts) == 1:
        return f"n = {counts[0]}"
    if counts == list(range(counts[0], counts[-1] + 1)):
        return f"n = {counts[0]}–{counts[-1]}"
    return "n = " + ", ".join(str(value) for value in counts)


def plot_total_fit(total_frame: pd.DataFrame, output: Path) -> None:
    with manuscript_style_context():
        figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), sharey=False)
        for axis, exposure in zip(axes, EXPOSURES, strict=True):
            subset = total_frame[
                total_frame["condition"] == exposure.condition
            ].sort_values("time_h")
            observed = (
                subset["observed_total_particles_per_ml"].to_numpy(dtype=float) / 1.0e9
            )
            observed_sd = (
                subset["observed_total_particles_per_ml_sd"].to_numpy(dtype=float)
                / 1.0e9
            )
            predicted = (
                subset["predicted_total_particles_per_ml"].to_numpy(dtype=float) / 1.0e9
            )
            axis.errorbar(
                subset["time_h"],
                observed,
                yerr=observed_sd,
                color=OBSERVED_COLOR,
                marker="o",
                linestyle="-",
                linewidth=1.8,
                capsize=3.0,
                label=OBSERVED_MEAN_SD_LABEL,
            )
            axis.plot(
                subset["time_h"],
                predicted,
                color=FITTED_COLOR,
                marker="s",
                linestyle="--",
                linewidth=1.8,
                label=FITTED_MODEL_LABEL,
            )
            values = np.r_[observed - observed_sd, observed + observed_sd, predicted]
            span = max(float(np.nanmax(values) - np.nanmin(values)), 0.2)
            axis.set_ylim(
                max(0.0, float(np.nanmin(values)) - 0.12 * span),
                float(np.nanmax(values)) + 0.12 * span,
            )
            style_manuscript_axis(
                axis,
                x_label="Time after pulse exposure (h)",
                title=CONDITION_LABELS[exposure.condition],
            )
            axis.set_xticks(OBSERVATION_TIMES_H)
        axes[0].set_ylabel(
            r"Particle concentration, 80--380 nm ($10^9$ particles mL$^{-1}$)"
        )
        place_manuscript_legend(figure, axes, multi_panel=True, location="right")
        figure.suptitle("Extracellular particle concentration", fontsize=11)
        figure.subplots_adjust(
            left=0.08, right=0.97, bottom=0.20, top=0.78, wspace=0.28
        )
        add_figure_note(
            figure,
            "Points and error bars show mean ± SD across available measurements "
            f"({_measurement_count_note(total_frame, 'measurement_count')}; independence unconfirmed).",
            x=0.08,
            y=0.02,
            reserve_bottom=0.18,
        )
        save_manuscript_figure(figure, output, dpi=MANUSCRIPT_COLOR_DPI)
        plt.close(figure)


def plot_size_profiles(prediction: pd.DataFrame, output: Path) -> None:
    with manuscript_style_context():
        figure, axes = plt.subplots(
            3, 3, figsize=(12.0, 8.8), sharex=True, sharey=False
        )
        for row_index, exposure in enumerate(EXPOSURES):
            for column_index, time_h in enumerate(OBSERVATION_TIMES_H):
                axis = axes[row_index, column_index]
                subset = prediction[
                    (prediction["condition"] == exposure.condition)
                    & (prediction["time_h"] == time_h)
                ].sort_values("size_bin_center_nm")
                diameter = subset["size_bin_center_nm"].to_numpy(dtype=float)
                observed = (
                    subset["observed_particles_per_ml_mean"].to_numpy(dtype=float)
                    / 1.0e9
                )
                observed_sd = (
                    subset["observed_particles_per_ml_sd"].to_numpy(dtype=float) / 1.0e9
                )
                predicted = (
                    subset["predicted_particles_per_ml"].to_numpy(dtype=float) / 1.0e9
                )
                axis.errorbar(
                    diameter,
                    observed,
                    yerr=observed_sd,
                    color=OBSERVED_COLOR,
                    marker="o",
                    linestyle="-",
                    linewidth=1.5,
                    markersize=3.5,
                    capsize=2.5,
                    elinewidth=0.9,
                    label=OBSERVED_MEAN_SD_LABEL,
                )
                axis.plot(
                    diameter,
                    predicted,
                    color=FITTED_COLOR,
                    marker="s",
                    linestyle="--",
                    linewidth=1.5,
                    markersize=3.5,
                    label=FITTED_MODEL_LABEL,
                )
                style_manuscript_axis(
                    axis,
                    title=(f"{CONDITION_LABELS[exposure.condition]}, {time_h:g} h"),
                )
                axis.set_ylim(bottom=0.0)
                if row_index == 2:
                    axis.set_xlabel("Particle diameter (nm)")
                if column_index == 0:
                    axis.set_ylabel(r"Bin concentration ($10^9$ particles mL$^{-1}$)")
        place_manuscript_legend(figure, axes, multi_panel=True, location="right")
        figure.suptitle("Particle size distributions", fontsize=11)
        figure.subplots_adjust(
            left=0.08, right=0.97, bottom=0.12, top=0.92, wspace=0.24, hspace=0.34
        )
        add_figure_note(
            figure,
            "Observed curves show mean ± SD across available measurements "
            f"({_measurement_count_note(prediction, 'observed_measurement_count')}; independence unconfirmed). "
            "Concentrations are integrated within 20-nm bands.",
            x=0.08,
            y=0.015,
            reserve_bottom=0.10,
        )
        save_manuscript_figure(figure, output, dpi=MANUSCRIPT_COLOR_DPI)
        plt.close(figure)


def plot_size_time_heatmaps(prediction: pd.DataFrame, output: Path) -> None:
    grids: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    finite_values: list[np.ndarray] = []
    for exposure in EXPOSURES:
        for column in (
            "observed_particles_per_ml_mean",
            "predicted_particles_per_ml",
        ):
            size_nm, time_h, values = _mean_size_time_grid(
                prediction, exposure.condition, column
            )
            log_values = np.log10(np.clip(values, 1.0, None))
            grids[(exposure.condition, column)] = (size_nm, time_h, log_values)
            finite_values.append(log_values[np.isfinite(log_values)])
    color_min = float(np.min(np.concatenate(finite_values)))
    color_max = float(np.max(np.concatenate(finite_values)))

    with manuscript_style_context():
        figure = plt.figure(figsize=(12.8, 6.5))
        layout = figure.add_gridspec(
            2,
            4,
            width_ratios=(1.0, 1.0, 1.0, 0.055),
            left=0.07,
            right=0.94,
            bottom=0.14,
            top=0.84,
            wspace=0.20,
            hspace=0.30,
        )
        axes = np.empty((2, 3), dtype=object)
        color_axis = figure.add_subplot(layout[:, 3])
        image = None
        for column_index, exposure in enumerate(EXPOSURES):
            for row_index, column in enumerate(
                (
                    "observed_particles_per_ml_mean",
                    "predicted_particles_per_ml",
                )
            ):
                axis = figure.add_subplot(layout[row_index, column_index])
                axes[row_index, column_index] = axis
                size_nm, time_h, log_values = grids[(exposure.condition, column)]
                image = axis.pcolormesh(
                    size_nm,
                    time_h,
                    log_values,
                    shading="nearest",
                    cmap="viridis",
                    vmin=color_min,
                    vmax=color_max,
                )
                if row_index == 0:
                    axis.set_title(CONDITION_LABELS[exposure.condition])
                if row_index == 0 and column_index == 0:
                    axis.set_ylabel(f"{OBSERVED_MEAN_LABEL}\nTime after exposure (h)")
                elif row_index == 1 and column_index == 0:
                    axis.set_ylabel(f"{FITTED_MODEL_LABEL}\nTime after exposure (h)")
                elif column_index == 0:
                    axis.set_ylabel("Time after exposure (h)")
                if row_index == 1:
                    axis.set_xlabel("Particle diameter (nm)")
                axis.set_xticks([90, 150, 210, 270, 330, 370])
                axis.set_yticks(OBSERVATION_TIMES_H)
                axis.tick_params(direction="out", width=0.8, length=3.0)
        if image is not None:
            colorbar = figure.colorbar(image, cax=color_axis)
            colorbar.set_label(r"log$_{10}$ bin concentration (particles mL$^{-1}$)")
        figure.suptitle("Particle concentration by size and time", fontsize=11)
        add_figure_note(
            figure,
            "Observed cells show means at 0.5, 1, and 3 h; colored cells do not add intermediate observations. SD is shown in the profile figure and exported table.",
            x=0.07,
            y=0.02,
            reserve_bottom=0.12,
        )
        save_manuscript_figure(figure, output, dpi=MANUSCRIPT_COLOR_DPI)
        plt.close(figure)


def _mean_size_time_grid(
    prediction: pd.DataFrame,
    condition: str,
    value_column: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a model-facing batch mean on the measured time-by-size grid."""

    subset = prediction[prediction["condition"] == condition]
    pivot = subset.pivot(
        index="time_h",
        columns="size_bin_center_nm",
        values=value_column,
    ).reindex(index=OBSERVATION_TIMES_H)
    if pivot.isna().any().any():
        raise ValueError(
            f"Missing values in {condition} size-time grid for {value_column}"
        )
    sizes_nm = pivot.columns.to_numpy(dtype=float)
    times_h = pivot.index.to_numpy(dtype=float)
    return sizes_nm, times_h, pivot.to_numpy(dtype=float)


def signed_normalized_error_percent(
    observed: np.ndarray,
    predicted: np.ndarray,
) -> np.ndarray:
    """Return bounded signed error, 100*(predicted-observed)/(predicted+observed)."""

    observed_array = np.asarray(observed, dtype=float)
    predicted_array = np.asarray(predicted, dtype=float)
    if observed_array.shape != predicted_array.shape:
        raise ValueError("Observed and predicted arrays must have the same shape")
    if np.any(observed_array < 0.0) or np.any(predicted_array < 0.0):
        raise ValueError("Observed and predicted concentrations must be nonnegative")
    denominator = predicted_array + observed_array
    return np.divide(
        100.0 * (predicted_array - observed_array),
        denominator,
        out=np.full_like(denominator, np.nan, dtype=float),
        where=denominator > 0.0,
    )


def plot_size_time_surface_overlay(
    prediction: pd.DataFrame,
    output: Path,
) -> None:
    """Overlay measured and fitted size-time concentration surfaces."""

    grids: list[tuple[Exposure, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    global_max_billions = 0.0
    for exposure in EXPOSURES:
        sizes_nm, times_h, observed = _mean_size_time_grid(
            prediction,
            exposure.condition,
            "observed_particles_per_ml_mean",
        )
        _, _, fitted = _mean_size_time_grid(
            prediction,
            exposure.condition,
            "predicted_particles_per_ml",
        )
        observed_billions = observed / 1.0e9
        fitted_billions = fitted / 1.0e9
        global_max_billions = max(
            global_max_billions,
            float(np.nanmax(observed_billions)),
            float(np.nanmax(fitted_billions)),
        )
        grids.append((exposure, sizes_nm, times_h, observed_billions, fitted_billions))

    with manuscript_style_context():
        figure = plt.figure(figsize=(14.8, 4.8))
        axes = []
        for panel_index, grid_values in enumerate(grids, start=1):
            exposure, sizes_nm, times_h, observed_billions, fitted_billions = (
                grid_values
            )
            axis = figure.add_subplot(1, 3, panel_index, projection="3d")
            axes.append(axis)
            size_grid, time_grid = np.meshgrid(sizes_nm, times_h)
            axis.plot_surface(
                size_grid,
                time_grid,
                observed_billions,
                color=OBSERVED_COLOR,
                alpha=0.58,
                linewidth=0.35,
                edgecolor=OBSERVED_COLOR,
                antialiased=True,
            )
            axis.plot_surface(
                size_grid,
                time_grid,
                fitted_billions,
                color=FITTED_COLOR,
                alpha=0.32,
                linewidth=0.45,
                edgecolor=FITTED_COLOR,
                antialiased=True,
            )
            axis.plot_wireframe(
                size_grid,
                time_grid,
                fitted_billions,
                color=FITTED_COLOR,
                linewidth=0.7,
                alpha=0.85,
            )
            axis.scatter(
                size_grid,
                time_grid,
                observed_billions,
                color=OBSERVED_COLOR,
                marker="o",
                s=8,
                depthshade=False,
                label=OBSERVED_MEAN_LABEL if panel_index == 1 else None,
            )
            axis.scatter(
                size_grid,
                time_grid,
                fitted_billions,
                color=FITTED_COLOR,
                marker="^",
                s=8,
                alpha=0.85,
                depthshade=False,
                label=FITTED_MODEL_LABEL if panel_index == 1 else None,
            )
            axis.set_title(CONDITION_LABELS[exposure.condition])
            axis.set_xlabel("Particle diameter (nm)", labelpad=7)
            axis.set_ylabel("Time after exposure (h)", labelpad=7)
            axis.set_yticks(OBSERVATION_TIMES_H)
            axis.set_zlim(0.0, global_max_billions * 1.06)
            axis.view_init(elev=27, azim=-132)
        place_manuscript_legend(figure, axes, multi_panel=True, location="right")
        figure.text(
            0.012,
            0.50,
            r"Bin concentration ($10^9$ particles mL$^{-1}$)",
            ha="center",
            va="center",
            rotation="vertical",
        )
        figure.suptitle("Particle concentration by size and time", fontsize=11)
        figure.subplots_adjust(
            left=0.01, right=0.96, bottom=0.18, top=0.82, wspace=0.08
        )
        add_figure_note(
            figure,
            "Surfaces connect observed means at 0.5, 1, and 3 h; SD is shown in the profile figure and exported table. Values are concentrations within 20-nm bands.",
            x=0.06,
            y=0.015,
            reserve_bottom=0.16,
        )
        save_manuscript_figure(figure, output, dpi=MANUSCRIPT_COLOR_DPI)
        plt.close(figure)


def plot_size_time_error_contours(
    prediction: pd.DataFrame,
    output: Path,
) -> None:
    """Plot signed and absolute bounded error on the measured size-time grid."""

    with manuscript_style_context():
        figure = plt.figure(figsize=(13.2, 6.8))
        layout = figure.add_gridspec(
            2,
            4,
            width_ratios=(1.0, 1.0, 1.0, 0.055),
            left=0.07,
            right=0.94,
            bottom=0.15,
            top=0.84,
            wspace=0.18,
            hspace=0.34,
        )
        axes = np.empty((2, 3), dtype=object)
        for row_index in range(2):
            for column_index in range(3):
                axes[row_index, column_index] = figure.add_subplot(
                    layout[row_index, column_index]
                )
        signed_color_axis = figure.add_subplot(layout[0, 3])
        absolute_color_axis = figure.add_subplot(layout[1, 3])
        signed_levels = np.linspace(-100.0, 100.0, 9)
        absolute_levels = np.linspace(0.0, 100.0, 6)
        signed_image = None
        absolute_image = None
        for column_index, exposure in enumerate(EXPOSURES):
            sizes_nm, times_h, observed = _mean_size_time_grid(
                prediction,
                exposure.condition,
                "observed_particles_per_ml_mean",
            )
            _, _, fitted = _mean_size_time_grid(
                prediction,
                exposure.condition,
                "predicted_particles_per_ml",
            )
            size_grid, time_grid = np.meshgrid(sizes_nm, times_h)
            signed_error = signed_normalized_error_percent(observed, fitted)
            signed_image = axes[0, column_index].contourf(
                size_grid,
                time_grid,
                signed_error,
                levels=signed_levels,
                cmap="RdBu_r",
            )
            absolute_image = axes[1, column_index].contourf(
                size_grid,
                time_grid,
                np.abs(signed_error),
                levels=absolute_levels,
                cmap="magma",
            )
            axes[0, column_index].set_title(CONDITION_LABELS[exposure.condition])
            axes[1, column_index].set_xlabel("Particle diameter (nm)")
            for row_index in range(2):
                axes[row_index, column_index].set_yticks(OBSERVATION_TIMES_H)
                axes[row_index, column_index].scatter(
                    size_grid,
                    time_grid,
                    color="black",
                    s=4,
                    alpha=0.28,
                )
                axes[row_index, column_index].tick_params(
                    direction="out", width=0.8, length=3.0
                )
            if column_index == 0:
                axes[0, column_index].set_ylabel("Time after exposure (h)")
                axes[1, column_index].set_ylabel("Time after exposure (h)")
        if signed_image is not None:
            signed_bar = figure.colorbar(signed_image, cax=signed_color_axis)
            signed_bar.set_label("Signed normalized difference (%)")
        if absolute_image is not None:
            absolute_bar = figure.colorbar(absolute_image, cax=absolute_color_axis)
            absolute_bar.set_label("Absolute normalized difference (%)")
        figure.suptitle("Model error by particle size and time", fontsize=11)
        add_figure_note(
            figure,
            "Normalized difference = 100(model − experiment)/(model + experiment). Blue indicates underprediction; red indicates overprediction. Contours connect only the measured 0.5, 1, and 3 h slices. Low-concentration tail bands can show large relative differences; no detection-limit mask was available.",
            x=0.07,
            y=0.015,
            reserve_bottom=0.14,
        )
        save_manuscript_figure(figure, output, dpi=MANUSCRIPT_COLOR_DPI)
        plt.close(figure)


def write_outputs(
    output_dir: Path,
    fitted: dict[
        str, tuple[SizeResolvedFFRCIPredictor, object, list[object], pd.DataFrame]
    ],
    bridge: ExperimentalObservationBridge,
    *,
    data_source: Path,
    pulse_width_ns: float,
    repetition_rate_hz: float,
    starts: int,
    max_nfev: int,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    comparison_rows: list[pd.DataFrame] = []
    total_rows: list[pd.DataFrame] = []
    distribution_rows: list[pd.DataFrame] = []
    parameter_rows: list[dict[str, float | str | bool]] = []
    kernel_rows: list[pd.DataFrame] = []
    model_rows: list[dict[str, float | str | int | bool]] = []
    summaries: dict[str, object] = {}

    for variant_name, (predictor, best, candidates, prediction) in fitted.items():
        metrics, totals, distributions = prediction_metrics(prediction)
        data_residuals = predictor.data_residuals(prediction)
        prior_residuals = predictor.prior_residuals(best.x)
        prediction = prediction.copy()
        prediction.insert(0, "variant", variant_name)
        prediction["predicted_over_observed"] = prediction[
            "predicted_particles_per_ml"
        ] / prediction["observed_particles_per_ml_mean"].replace(0.0, np.nan)
        prediction["signed_normalized_error_percent"] = signed_normalized_error_percent(
            prediction["observed_particles_per_ml_mean"].to_numpy(dtype=float),
            prediction["predicted_particles_per_ml"].to_numpy(dtype=float),
        )
        prediction["model_within_observed_sd"] = (
            np.abs(
                prediction["predicted_particles_per_ml"]
                - prediction["observed_particles_per_ml_mean"]
            )
            <= prediction["observed_particles_per_ml_sd"]
        )
        prediction["sd_scaled_difference"] = np.divide(
            prediction["predicted_particles_per_ml"]
            - prediction["observed_particles_per_ml_mean"],
            prediction["observed_particles_per_ml_sd"],
            out=np.full(len(prediction), np.nan, dtype=float),
            where=prediction["observed_particles_per_ml_sd"].to_numpy(dtype=float)
            > 0.0,
        )
        prediction["sqrt_fraction_residual"] = np.nan
        for _, group in prediction.groupby(["condition", "time_h"]):
            observed = group["observed_particles_per_ml_mean"].to_numpy(dtype=float)
            predicted_values = group["predicted_particles_per_ml"].to_numpy(dtype=float)
            prediction.loc[group.index, "sqrt_fraction_residual"] = np.sqrt(
                predicted_values / predicted_values.sum()
            ) - np.sqrt(observed / observed.sum())
        totals.insert(0, "variant", variant_name)
        distributions.insert(0, "variant", variant_name)
        comparison_rows.append(prediction)
        total_rows.append(totals)
        distribution_rows.append(distributions)
        parameter_rows.extend(_parameter_rows(predictor, best.x))
        kernel_rows.append(build_kernel_frame(predictor, best.x))
        n_data = int(len(data_residuals))
        n_parameters = int(len(best.x))
        data_sse = float(np.sum(data_residuals**2))
        descriptive_bic = float(
            n_data * math.log(max(data_sse / n_data, 1.0e-300))
            + n_parameters * math.log(n_data)
        )
        model_rows.append(
            {
                "variant": variant_name,
                "n_parameters": n_parameters,
                "n_data_residual_components": n_data,
                "data_residual_sse": data_sse,
                "prior_residual_sse_unweighted": float(np.sum(prior_residuals**2)),
                "optimizer_cost_with_priors": float(best.cost),
                "descriptive_bic_on_composite_residual": descriptive_bic,
                "optimizer_success": bool(best.success),
                **metrics,
            }
        )
        summaries[variant_name] = {
            "metrics": metrics,
            "optimizer": {
                "success": bool(best.success),
                "message": str(best.message),
                "requested_starts": int(starts),
                "max_nfev_per_start": int(max_nfev),
                "nfev_selected_start": int(best.nfev),
                "candidate_costs": [float(candidate.cost) for candidate in candidates],
                "predictor_evaluations": predictor.evaluations,
            },
            "parameters": decode_parameters(best.x, predictor.variant),
        }

    model_comparison = pd.DataFrame(model_rows).sort_values(
        "descriptive_bic_on_composite_residual"
    )
    selected_variant = str(model_comparison.iloc[0]["variant"])
    selected_predictor, selected_best, _, selected_prediction = fitted[selected_variant]
    model_observations = selected_predictor.observations
    measurement_counts = (
        model_observations[["condition", "time_h", "observed_measurement_count"]]
        .drop_duplicates()
        .sort_values(["condition", "time_h"])
    )
    source_histogram_count = int(measurement_counts["observed_measurement_count"].sum())
    control_histogram_count = int(
        model_observations["initial_control_measurement_count"].max()
    )
    model_distribution_count = int(len(measurement_counts))
    measurement_count_records = [
        {
            "condition": str(row.condition),
            "condition_label": CONDITION_LABELS[str(row.condition)],
            "time_h": float(row.time_h),
            "n": int(row.observed_measurement_count),
        }
        for row in measurement_counts.itertuples(index=False)
    ]

    comparisons = pd.concat(comparison_rows, ignore_index=True)
    totals = pd.concat(total_rows, ignore_index=True)
    distributions = pd.concat(distribution_rows, ignore_index=True)
    parameters = pd.DataFrame(parameter_rows)
    kernels = pd.concat(kernel_rows, ignore_index=True)
    model_observations.to_csv(
        output_dir / "experimental_batch_summary.csv", index=False
    )
    comparisons.to_csv(output_dir / "observed_vs_predicted_size_bins.csv", index=False)
    totals.to_csv(output_dir / "total_concentration_fit.csv", index=False)
    distributions.to_csv(output_dir / "size_distribution_fit_metrics.csv", index=False)
    parameters.to_csv(output_dir / "fitted_parameters.csv", index=False)
    kernels.to_csv(output_dir / "pathway_size_kernels.csv", index=False)
    model_comparison.to_csv(output_dir / "model_comparison.csv", index=False)

    selected_values = decode_parameters(selected_best.x, selected_predictor.variant)
    with (output_dir / "fitted_parameters.yml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            {
                "selected_variant": selected_variant,
                "size_resolved_kinetics": selected_values,
            },
            handle,
            sort_keys=False,
        )

    selected_totals = totals[totals["variant"] == selected_variant]
    plot_total_fit(
        selected_totals,
        output_dir / "longitudinal_total_fit.png",
    )
    plot_size_profiles(
        selected_prediction,
        output_dir / "size_profile_fit.png",
    )
    plot_size_time_heatmaps(
        selected_prediction,
        output_dir / "size_time_fit.png",
    )
    plot_size_time_surface_overlay(
        selected_prediction,
        output_dir / "size_time_surface_overlay.png",
    )
    plot_size_time_error_contours(
        selected_prediction,
        output_dir / "size_time_fit_error_contours.png",
    )

    summary = {
        "selected_variant_by_descriptive_composite_bic": selected_variant,
        "fit_variants": summaries,
        "data": {
            "source": str(data_source),
            "common_size_window_nm": [80.0, 380.0],
            "common_bin_width_nm": 20.0,
            "source_treatment_histograms": source_histogram_count,
            "source_initial_control_histograms": control_histogram_count,
            "batch_mean_size_distributions": model_distribution_count,
            "model_facing_total_concentration_targets": model_distribution_count,
            "measurement_counts_by_condition_time": measurement_count_records,
            "aggregation": "arithmetic mean with sample SD, SE, and count after common-bin rebinning",
            "fit_target": "one batch-mean distribution per condition/time",
            "control_initialization": "all-batch mean of the undated cell-containing control size distribution",
        },
        "experimental_bridge": bridge.to_metadata(),
        "scenario_assumptions": {
            "pulse_width_ns": pulse_width_ns,
            "repetition_rate_Hz": repetition_rate_hz,
            "source_voltage_label_used_as_model_field_kV_per_cm": True,
            "cell_density_per_ml": bridge.initial_cell_count / bridge.medium_volume_ml,
        },
        "objective": {
            "total_component": "one natural-log total-concentration residual per condition/time",
            "size_component": "square-root composition residuals comparing the model with one batch-mean size distribution per condition/time",
            "observed_variability": "sample SD is exported and plotted descriptively but does not inverse-variance weight the objective",
            "sd_scaled_differences": "descriptive differences divided by sample SD; they are not z-scores or confidence coverage",
            "regularization": "weak literature-centered penalties on pathway medians/widths, size-loss slope, state shifts, and optional dose-response correction",
            "zero_bins": "handled in composition space without logarithms or pseudocounts",
        },
        "visualizations": {
            "size_time_surface_overlay": "experimental mean surface with a translucent fitted surface, wireframe, and markers on a common concentration scale",
            "size_time_fit_error_contours": "signed and absolute normalized difference between the model and experimental mean; both use 100*(model-experiment)/(model+experiment) and are diagnostic rather than optimizer residuals",
        },
        "interpretation_limits": [
            "The cell-containing controls have no harvest time and are used as a provisional initial distribution.",
            "The repeated records are summarized as batch measurements, but their biological independence is not established by the source file.",
            "The fitted pathway decomposition is not identifiable from diameter alone and must be treated as a reduced-order explanation.",
            "The dose-response correction is a diagnostic adapter for weak pulse-condition separation in the current intracellular model, not a validated constitutive mechanism.",
            "The observation matrix is explicit but fixed to identity on the common 20-nm bands because Exoid calibration and pore metadata were not supplied.",
            "Relative error can appear large in low-concentration tail bins because no instrument detection-limit mask was available.",
        ],
    }
    (output_dir / "fit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    selected_metrics = summaries[selected_variant]["metrics"]
    readme = f"""# Size-resolved extracellular fit (v1.1)

This analysis maps the unchanged intracellular sEV, mlEV, and AB release
trajectories through state-conditioned lognormal size kernels, a smooth
size-dependent extracellular loss law, and an explicit observation operator.
Repeated measurements are rebinned first and summarized as one experimental
mean ± SD distribution per condition and time. The simulation is initialized
from the all-batch mean of the undated cell-containing control distribution.

Selected descriptive variant: `{selected_variant}`

- Model-facing total targets: {model_distribution_count}
- Batch-mean size distributions: {model_distribution_count}
- Source treatment histograms: {source_histogram_count}
- Source initial-control histograms: {control_histogram_count}
- Total log10 RMSE (80–380 nm): {selected_metrics["rmse_log10_total"]:.4f}
- Median absolute total error: {selected_metrics["median_absolute_percent_total_error"]:.1f}%
- Mean Hellinger size distance: {selected_metrics["mean_hellinger_size_distance"]:.4f}
- Mean absolute error in distribution mean diameter (80–380 nm): {selected_metrics["mae_mean_diameter_nm"]:.1f} nm
- Total targets within observed mean ± SD: {100.0 * selected_metrics["descriptive_fraction_totals_within_observed_sd"]:.1f}% (descriptive)
- Size bins within observed mean ± SD: {100.0 * selected_metrics["descriptive_fraction_size_bins_within_observed_sd"]:.1f}% (descriptive)

The selection score is descriptive because it combines total and composition
residuals rather than a fully specified sampling likelihood. Sample SD is
reported and plotted but does not weight the objective because the repeated
measurement count is small and size bins are correlated. The source gives no
time-zero control, confirmed medium volume, viability trajectory, instrument
response, dilution, or recovery. The fitted parameters are therefore
provisional reduced-order values, not constitutive biological estimates.

Files:

- `experimental_batch_summary.csv`: model-facing mean, SD, SE, and count after common-bin aggregation.
- `observed_vs_predicted_size_bins.csv`: every fitted 20-nm bin for all variants.
- `total_concentration_fit.csv`: the nine model-facing longitudinal total targets.
- `size_distribution_fit_metrics.csv`: condition/time mean-distribution errors.
- `pathway_size_kernels.csv`: pathway kernel mass and smooth loss by condition/time/bin.
- `fitted_parameters.csv` and `fitted_parameters.yml`: fitted values and bounds.
- `model_comparison.csv` and `fit_summary.json`: variant metrics and assumptions.
- `longitudinal_total_fit.png`, `size_profile_fit.png`, and `size_time_fit.png`: conventional fit figures.
- `size_time_surface_overlay.png`: experimental mean size-time surface with a translucent fitted layer.
- `size_time_fit_error_contours.png`: signed and absolute normalized-difference contours. These bounded diagnostics are distinct from the optimized residuals.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--experimental-bridge-config",
        type=Path,
        default=DEFAULT_BRIDGE_CONFIG,
    )
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--max-nfev", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--pulse-width-ns", type=float, default=60.0)
    parser.add_argument("--repetition-rate-hz", type=float, default=1.0)
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=[variant.name for variant in FIT_VARIANTS],
        default=[variant.name for variant in FIT_VARIANTS],
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.starts < 1 or args.max_nfev < 1:
        raise ValueError("--starts and --max-nfev must be positive")
    bridge = ExperimentalObservationBridge.from_yaml(args.experimental_bridge_config)
    raw_observations = load_size_resolved_observations(args.data_dir)
    observations = aggregate_size_resolved_observations(raw_observations)
    trajectories = build_upstream_trajectories(
        bridge,
        pulse_width_ns=args.pulse_width_ns,
        repetition_rate_hz=args.repetition_rate_hz,
    )
    fitted: dict[
        str,
        tuple[SizeResolvedFFRCIPredictor, object, list[object], pd.DataFrame],
    ] = {}
    variants_by_name = {variant.name: variant for variant in FIT_VARIANTS}
    for variant_index, variant_name in enumerate(args.variants):
        variant = variants_by_name[variant_name]
        predictor = SizeResolvedFFRCIPredictor(
            observations,
            trajectories,
            bridge,
            variant,
        )
        best, candidates = fit_variant(
            predictor,
            starts=args.starts,
            max_nfev=args.max_nfev,
            seed=args.seed + variant_index,
        )
        prediction = predictor.predict(best.x)
        metrics, _, _ = prediction_metrics(prediction)
        fitted[variant.name] = (predictor, best, candidates, prediction)
        print(
            f"{variant.name}: total log10 RMSE={metrics['rmse_log10_total']:.4f}, "
            f"mean size Hellinger={metrics['mean_hellinger_size_distance']:.4f}, "
            f"success={best.success}"
        )
    summary = write_outputs(
        args.output_dir,
        fitted,
        bridge,
        data_source=args.data_dir / EXOID_FILENAME,
        pulse_width_ns=args.pulse_width_ns,
        repetition_rate_hz=args.repetition_rate_hz,
        starts=args.starts,
        max_nfev=args.max_nfev,
    )
    print(
        "Selected descriptive variant: "
        f"{summary['selected_variant_by_descriptive_composite_bic']}"
    )
    print(f"Results written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
