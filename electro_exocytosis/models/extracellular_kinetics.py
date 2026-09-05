from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from typing import Any

import numpy as np
from scipy.linalg import expm


EXTRACELLULAR_CLASSES = ("sEV", "mlEV", "AB")


@dataclass(slots=True)
class ExtracellularKineticsParams:
    """Parameters for the extracellular EV stock and observation layer.

    The three stocks are particle equivalents in the complete culture volume.
    Release is supplied by the intracellular model on a per-cell basis. Uptake,
    degradation, and adsorption are deliberately separate so experiments can
    later distinguish them. ``effective_loss_rate_s`` is the recommended first
    calibration target when the data identify disappearance but cannot
    distinguish its mechanism.
    """

    source_scale_particles_per_model_unit: float = 1.0
    initial_sEV_concentration_particles_per_ml: float = 0.0
    initial_mlEV_concentration_particles_per_ml: float = 0.0
    initial_AB_concentration_particles_per_ml: float = 0.0
    effective_loss_rate_s: float = 0.0
    uptake_rate_s: float = 0.0
    degradation_rate_s: float = 0.0
    adsorption_rate_s: float = 0.0
    sEV_loss_multiplier: float = 1.0
    mlEV_loss_multiplier: float = 1.0
    AB_loss_multiplier: float = 1.0
    sEV_to_mlEV_aggregation_rate_s: float = 0.0
    aggregation_particle_yield: float = 0.5
    assay_recovery_fraction: float = 1.0
    assay_dilution_factor: float = 1.0
    assay_background_concentration_particles_per_ml: float = 0.0


@dataclass(frozen=True, slots=True)
class MediumSamplingEvent:
    """Well-mixed removal followed by optional particle-free replacement."""

    time_s: float
    sampled_volume_ml: float = 0.0
    replacement_volume_ml: float = 0.0


@dataclass(slots=True)
class ExtracellularKineticsResult:
    """Extracellular time series and an auditable event log."""

    timeseries: dict[str, np.ndarray]
    event_log: list[dict[str, float]]


def simulate_extracellular_kinetics(
    time_s: Sequence[float] | np.ndarray,
    cumulative_release_per_cell: Mapping[str, Sequence[float] | np.ndarray],
    params: ExtracellularKineticsParams,
    *,
    initial_cell_density_per_ml: float,
    initial_volume_ml: float,
    viable_fraction: Sequence[float] | np.ndarray | None = None,
    sampling_events: Sequence[MediumSamplingEvent] = (),
) -> ExtracellularKineticsResult:
    """Propagate extracellular particle stocks from intracellular release.

    Release increments from the upstream ODE are used as the source so the
    zero-loss result exactly retains the upstream cumulative release at its
    output times. Within each interval, that source is treated as constant and
    the linear loss/aggregation system is integrated analytically.

    A sampling event is applied *after* the state at its time is recorded. This
    matches a sample-then-replace protocol: the reported concentration is the
    withdrawn sample, while removal/replacement affects subsequent states.
    """

    t = np.asarray(time_s, dtype=float)
    if t.ndim != 1 or len(t) < 2 or np.any(~np.isfinite(t)):
        raise ValueError("time_s must be a finite one-dimensional array with at least two points")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("time_s must be strictly increasing")
    if not np.isfinite(initial_cell_density_per_ml) or initial_cell_density_per_ml <= 0.0:
        raise ValueError("initial_cell_density_per_ml must be finite and positive")
    if not np.isfinite(initial_volume_ml) or initial_volume_ml <= 0.0:
        raise ValueError("initial_volume_ml must be finite and positive")

    _validate_params(params)
    cumulative = _coerce_cumulative_release(cumulative_release_per_cell, len(t))
    viability = _coerce_viability(viable_fraction, len(t))
    events = _coerce_events(sampling_events, float(t[0]), float(t[-1]))

    event_times = np.array([event.time_s for event in events], dtype=float)
    timeline = np.unique(np.concatenate([t, event_times])) if len(events) else t.copy()
    cumulative_timeline = np.column_stack(
        [np.interp(timeline, t, cumulative[:, index]) for index in range(3)]
    )
    viability_timeline = np.interp(timeline, t, viability)

    initial_cells = float(initial_cell_density_per_ml) * float(initial_volume_ml)
    volume_ml = float(initial_volume_ml)
    stocks = np.array(
        [
            params.initial_sEV_concentration_particles_per_ml * volume_ml,
            params.initial_mlEV_concentration_particles_per_ml * volume_ml,
            params.initial_AB_concentration_particles_per_ml * volume_ml,
        ],
        dtype=float,
    )
    dynamics = _stock_dynamics_matrix(params)
    events_by_time: dict[float, list[MediumSamplingEvent]] = {}
    for event in events:
        events_by_time.setdefault(float(event.time_s), []).append(event)

    records: dict[float, tuple[np.ndarray, float, float]] = {}
    event_log: list[dict[str, float]] = []
    previous_source_rate = 0.0

    for index, current_time in enumerate(timeline):
        if index > 0:
            previous_time = float(timeline[index - 1])
            dt = float(current_time - previous_time)
            release_increment = cumulative_timeline[index] - cumulative_timeline[index - 1]
            if np.any(release_increment < -1.0e-8):
                raise ValueError("cumulative release must be non-decreasing for every EV class")
            release_increment = np.clip(release_increment, 0.0, None)
            mean_viability = 0.5 * (
                float(viability_timeline[index - 1]) + float(viability_timeline[index])
            )
            source_rate = (
                release_increment
                * initial_cells
                * params.source_scale_particles_per_model_unit
                * mean_viability
                / dt
            )
            stocks = _advance_linear_stock(stocks, source_rate, dynamics, dt)
            previous_source_rate = float(np.sum(source_rate) / volume_ml)

        current_time_float = float(current_time)
        if np.any(t == current_time_float):
            records[current_time_float] = (
                np.clip(stocks.copy(), 0.0, None),
                volume_ml,
                previous_source_rate,
            )

        for event in events_by_time.get(current_time_float, []):
            if event.sampled_volume_ml > volume_ml + 1.0e-12:
                raise ValueError(
                    f"sampling event at {event.time_s:g} s removes "
                    f"{event.sampled_volume_ml:g} mL from only {volume_ml:g} mL"
                )
            volume_before = volume_ml
            stocks_before = stocks.copy()
            retained_fraction = max(
                0.0, 1.0 - float(event.sampled_volume_ml) / volume_before
            )
            stocks *= retained_fraction
            volume_ml = volume_before - float(event.sampled_volume_ml)
            volume_ml += float(event.replacement_volume_ml)
            if volume_ml <= 0.0:
                raise ValueError(
                    f"sampling event at {event.time_s:g} s leaves no culture medium"
                )
            event_log.append(
                {
                    "time_s": current_time_float,
                    "sampled_volume_ml": float(event.sampled_volume_ml),
                    "replacement_volume_ml": float(event.replacement_volume_ml),
                    "volume_before_ml": float(volume_before),
                    "volume_after_ml": float(volume_ml),
                    "particles_removed": float(np.sum(stocks_before - stocks)),
                }
            )

    stock_rows = np.vstack([records[float(value)][0] for value in t])
    volumes = np.array([records[float(value)][1] for value in t], dtype=float)
    source_rates = np.array([records[float(value)][2] for value in t], dtype=float)
    concentrations = stock_rows / volumes[:, None]

    loss_components = _loss_component_rates(params)
    effective_loss_flux = concentrations * loss_components["effective"][None, :]
    uptake_flux = concentrations * loss_components["uptake"][None, :]
    degradation_flux = concentrations * loss_components["degradation"][None, :]
    adsorption_flux = concentrations * loss_components["adsorption"][None, :]
    aggregation_loss_flux = (
        concentrations[:, 0] * params.sEV_to_mlEV_aggregation_rate_s
    )
    aggregation_gain_flux = aggregation_loss_flux * params.aggregation_particle_yield
    total_concentration = np.sum(concentrations, axis=1)
    measured_concentration = (
        params.assay_recovery_fraction
        * total_concentration
        / params.assay_dilution_factor
        + params.assay_background_concentration_particles_per_ml
    )

    timeseries = {
        "sEV_extracellular_particles": stock_rows[:, 0],
        "mlEV_extracellular_particles": stock_rows[:, 1],
        "AB_extracellular_particles": stock_rows[:, 2],
        "sEV_extracellular_concentration_particles_per_ml": concentrations[:, 0],
        "mlEV_extracellular_concentration_particles_per_ml": concentrations[:, 1],
        "AB_extracellular_concentration_particles_per_ml": concentrations[:, 2],
        "total_extracellular_concentration_particles_per_ml": total_concentration,
        "measured_particle_concentration_particles_per_ml": measured_concentration,
        "extracellular_medium_volume_ml": volumes,
        "viable_producer_fraction": viability,
        "extracellular_source_rate_particles_per_ml_s": source_rates,
        "sEV_effective_loss_flux_particles_per_ml_s": effective_loss_flux[:, 0],
        "mlEV_effective_loss_flux_particles_per_ml_s": effective_loss_flux[:, 1],
        "AB_effective_loss_flux_particles_per_ml_s": effective_loss_flux[:, 2],
        "sEV_uptake_flux_particles_per_ml_s": uptake_flux[:, 0],
        "mlEV_uptake_flux_particles_per_ml_s": uptake_flux[:, 1],
        "AB_uptake_flux_particles_per_ml_s": uptake_flux[:, 2],
        "sEV_degradation_flux_particles_per_ml_s": degradation_flux[:, 0],
        "mlEV_degradation_flux_particles_per_ml_s": degradation_flux[:, 1],
        "AB_degradation_flux_particles_per_ml_s": degradation_flux[:, 2],
        "sEV_adsorption_flux_particles_per_ml_s": adsorption_flux[:, 0],
        "mlEV_adsorption_flux_particles_per_ml_s": adsorption_flux[:, 1],
        "AB_adsorption_flux_particles_per_ml_s": adsorption_flux[:, 2],
        "sEV_aggregation_loss_flux_particles_per_ml_s": aggregation_loss_flux,
        "mlEV_aggregation_gain_flux_particles_per_ml_s": aggregation_gain_flux,
    }
    return ExtracellularKineticsResult(timeseries=timeseries, event_log=event_log)


def coerce_extracellular_kinetics_params(
    params: Mapping[str, Any] | ExtracellularKineticsParams | None,
) -> ExtracellularKineticsParams:
    """Build extracellular parameters from flat or nested configuration."""

    if params is None:
        result = ExtracellularKineticsParams()
    elif isinstance(params, ExtracellularKineticsParams):
        result = params
    else:
        nested_params = params.get("extracellular_kinetics")
        if isinstance(nested_params, Mapping):
            params = nested_params
        defaults = ExtracellularKineticsParams()
        values = {
            field.name: float(params.get(field.name, getattr(defaults, field.name)))
            for field in fields(ExtracellularKineticsParams)
        }
        result = ExtracellularKineticsParams(**values)
    _validate_params(result)
    return result


def extracellular_kinetics_defaults() -> dict[str, float]:
    """Return extracellular kinetics defaults as a plain dictionary."""

    return asdict(ExtracellularKineticsParams())


def _advance_linear_stock(
    stocks: np.ndarray,
    source_rate: np.ndarray,
    dynamics: np.ndarray,
    dt: float,
) -> np.ndarray:
    augmented = np.zeros((4, 4), dtype=float)
    augmented[:3, :3] = dynamics
    augmented[:3, 3] = source_rate
    transition = expm(augmented * float(dt))
    next_stocks = transition[:3, :3] @ stocks + transition[:3, 3]
    return np.clip(next_stocks, 0.0, None)


def _stock_dynamics_matrix(params: ExtracellularKineticsParams) -> np.ndarray:
    component_rates = _loss_component_rates(params)
    total_loss = (
        component_rates["effective"]
        + component_rates["uptake"]
        + component_rates["degradation"]
        + component_rates["adsorption"]
    )
    aggregation = params.sEV_to_mlEV_aggregation_rate_s
    dynamics = np.diag(-total_loss)
    dynamics[0, 0] -= aggregation
    dynamics[1, 0] += params.aggregation_particle_yield * aggregation
    return dynamics


def _loss_component_rates(params: ExtracellularKineticsParams) -> dict[str, np.ndarray]:
    multipliers = np.array(
        [params.sEV_loss_multiplier, params.mlEV_loss_multiplier, params.AB_loss_multiplier],
        dtype=float,
    )
    return {
        "effective": params.effective_loss_rate_s * multipliers,
        "uptake": params.uptake_rate_s * multipliers,
        "degradation": params.degradation_rate_s * multipliers,
        "adsorption": params.adsorption_rate_s * multipliers,
    }


def _coerce_cumulative_release(
    values: Mapping[str, Sequence[float] | np.ndarray],
    expected_length: int,
) -> np.ndarray:
    missing = [name for name in EXTRACELLULAR_CLASSES if name not in values]
    if missing:
        raise ValueError(f"missing cumulative release classes: {missing}")
    columns = []
    for name in EXTRACELLULAR_CLASSES:
        column = np.asarray(values[name], dtype=float)
        if column.shape != (expected_length,) or np.any(~np.isfinite(column)):
            raise ValueError(f"cumulative release for {name} has invalid shape or values")
        if np.any(column < -1.0e-12):
            raise ValueError(f"cumulative release for {name} must be nonnegative")
        columns.append(column)
    return np.column_stack(columns)


def _coerce_viability(
    viable_fraction: Sequence[float] | np.ndarray | None,
    expected_length: int,
) -> np.ndarray:
    if viable_fraction is None:
        return np.ones(expected_length, dtype=float)
    viability = np.asarray(viable_fraction, dtype=float)
    if viability.shape != (expected_length,) or np.any(~np.isfinite(viability)):
        raise ValueError("viable_fraction has invalid shape or values")
    if np.any((viability < 0.0) | (viability > 1.0)):
        raise ValueError("viable_fraction must remain between zero and one")
    return viability


def _coerce_events(
    events: Sequence[MediumSamplingEvent],
    start_time_s: float,
    end_time_s: float,
) -> list[MediumSamplingEvent]:
    coerced = []
    for event in events:
        if isinstance(event, MediumSamplingEvent):
            current = event
        elif isinstance(event, Mapping):
            current = MediumSamplingEvent(**dict(event))
        elif hasattr(event, "model_dump"):
            current = MediumSamplingEvent(**event.model_dump())
        else:
            raise TypeError("sampling events must be mappings or event objects")
        if not np.isfinite(current.time_s) or not start_time_s <= current.time_s <= end_time_s:
            raise ValueError("sampling event time must lie inside the simulation interval")
        if (
            not np.isfinite(current.sampled_volume_ml)
            or not np.isfinite(current.replacement_volume_ml)
            or current.sampled_volume_ml < 0.0
            or current.replacement_volume_ml < 0.0
        ):
            raise ValueError("sampling and replacement volumes must be finite and nonnegative")
        coerced.append(current)
    return sorted(coerced, key=lambda event: event.time_s)


def _validate_params(params: ExtracellularKineticsParams) -> None:
    values = asdict(params)
    for name, value in values.items():
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")
        if value < 0.0:
            raise ValueError(f"{name} must be nonnegative")
    if params.assay_recovery_fraction > 1.0:
        raise ValueError("assay_recovery_fraction cannot exceed one")
    if params.assay_dilution_factor <= 0.0:
        raise ValueError("assay_dilution_factor must be positive")
    if params.aggregation_particle_yield > 1.0:
        raise ValueError("aggregation_particle_yield cannot exceed one")
