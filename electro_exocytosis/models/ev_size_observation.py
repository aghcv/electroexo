"""Reduced pathway-to-size bridge for extracellular particle observations.

The intracellular EV model reports cumulative release for three mechanistic
pathways.  This module maps increments in those pathways to a common latent
diameter grid, propagates their extracellular concentrations through a smooth
size-dependent loss, and applies a separate instrument-size observation
operator.  It is intentionally standalone so existing cumulative and
three-class extracellular outputs remain backward compatible.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr


SIZE_PATHWAYS = ("sEV", "mlEV", "AB")


@dataclass(frozen=True, slots=True)
class PathwaySizeKernelParams:
    """Lognormal size kernel and release scale for one EV pathway.

    ``geometric_sd`` is the multiplicative geometric standard deviation and
    must exceed one.  At an interval midpoint the kernel median is

    ``median_diameter_nm * exp(state_shift_coefficient * state_signal)``.
    """

    median_diameter_nm: float
    geometric_sd: float
    source_scale_particles_per_model_unit: float = 1.0
    state_shift_coefficient: float = 0.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.median_diameter_nm) or self.median_diameter_nm <= 0.0:
            raise ValueError("median_diameter_nm must be finite and positive")
        if not np.isfinite(self.geometric_sd) or self.geometric_sd <= 1.0:
            raise ValueError("geometric_sd must be finite and greater than one")
        if (
            not np.isfinite(self.source_scale_particles_per_model_unit)
            or self.source_scale_particles_per_model_unit < 0.0
        ):
            raise ValueError(
                "source_scale_particles_per_model_unit must be finite and nonnegative"
            )
        if not np.isfinite(self.state_shift_coefficient):
            raise ValueError("state_shift_coefficient must be finite")


@dataclass(frozen=True, slots=True)
class SizeResolvedKineticsParams:
    """Shared extracellular-loss and instrument-response parameters."""

    effective_loss_rate_s: float = 0.0
    loss_reference_diameter_nm: float = 150.0
    loss_size_exponent: float = 0.0
    instrument_log_diameter_sd: float = 0.0
    assay_recovery_fraction: float = 1.0

    def __post_init__(self) -> None:
        if (
            not np.isfinite(self.effective_loss_rate_s)
            or self.effective_loss_rate_s < 0.0
        ):
            raise ValueError("effective_loss_rate_s must be finite and nonnegative")
        if (
            not np.isfinite(self.loss_reference_diameter_nm)
            or self.loss_reference_diameter_nm <= 0.0
        ):
            raise ValueError("loss_reference_diameter_nm must be finite and positive")
        if not np.isfinite(self.loss_size_exponent):
            raise ValueError("loss_size_exponent must be finite")
        if (
            not np.isfinite(self.instrument_log_diameter_sd)
            or self.instrument_log_diameter_sd < 0.0
        ):
            raise ValueError(
                "instrument_log_diameter_sd must be finite and nonnegative"
            )
        if (
            not np.isfinite(self.assay_recovery_fraction)
            or not 0.0 <= self.assay_recovery_fraction <= 1.0
        ):
            raise ValueError("assay_recovery_fraction must lie between zero and one")


@dataclass(frozen=True, slots=True)
class SizeResolvedKineticsResult:
    """Latent size stocks, pathway attribution, and assay-facing output.

    Concentration arrays use time along axis 0 and size bin along axis 1.
    Each value in ``pathway_latent_concentration_particles_per_ml`` has shape
    ``(n_time, n_latent_bins)``.  Each kernel array has shape
    ``(n_time - 1, n_latent_bins)`` and records the midpoint kernel used for
    the corresponding integration interval.  The observation matrix has shape
    ``(n_observed_bins, n_latent_bins)``.
    """

    time_s: np.ndarray
    latent_bin_edges_nm: np.ndarray
    observed_bin_edges_nm: np.ndarray
    ambient_latent_concentration_particles_per_ml: np.ndarray
    pathway_latent_concentration_particles_per_ml: dict[str, np.ndarray]
    total_latent_concentration_particles_per_ml: np.ndarray
    observed_concentration_particles_per_ml: np.ndarray
    kernel_probabilities: dict[str, np.ndarray]
    loss_rates_s: np.ndarray
    observation_matrix: np.ndarray

    @property
    def time(self) -> np.ndarray:
        """Alias for consumers that use an unqualified time field."""

        return self.time_s

    @property
    def bin_edges_nm(self) -> np.ndarray:
        """Alias for the latent-bin edges."""

        return self.latent_bin_edges_nm

    @property
    def ambient_latent_concentrations(self) -> np.ndarray:
        """Compact alias for the ambient latent concentration history."""

        return self.ambient_latent_concentration_particles_per_ml

    @property
    def pathway_latent_concentrations(self) -> dict[str, np.ndarray]:
        """Compact alias for pathway-attributed concentration histories."""

        return self.pathway_latent_concentration_particles_per_ml

    @property
    def total_latent_concentrations(self) -> np.ndarray:
        """Compact alias for the total latent concentration history."""

        return self.total_latent_concentration_particles_per_ml

    @property
    def observed_concentrations(self) -> np.ndarray:
        """Compact alias for assay-facing observed concentrations."""

        return self.observed_concentration_particles_per_ml


# Explicit long-form alias for callers that name the extracellular role.
SizeResolvedExtracellularKineticsResult = SizeResolvedKineticsResult


def lognormal_bin_probabilities(
    bin_edges_nm: Sequence[float] | np.ndarray,
    median_diameter_nm: float,
    geometric_sd: float,
) -> np.ndarray:
    """Integrate and normalize a lognormal kernel over diameter bins.

    The returned probabilities sum to one over the represented latent domain.
    This conditional normalization keeps pathway release mass-conserving when
    a reduced finite size range is modeled; callers should choose that range
    wide enough that omitted physical tail mass is negligible.
    """

    edges = _validate_bin_edges(bin_edges_nm, name="bin_edges_nm")
    if not np.isfinite(median_diameter_nm) or median_diameter_nm <= 0.0:
        raise ValueError("median_diameter_nm must be finite and positive")
    if not np.isfinite(geometric_sd) or geometric_sd <= 1.0:
        raise ValueError("geometric_sd must be finite and greater than one")

    log_sigma = float(np.log(geometric_sd))
    z = (np.log(edges) - np.log(float(median_diameter_nm))) / log_sigma
    probabilities = np.diff(ndtr(z))
    probabilities = np.clip(probabilities, 0.0, None)
    represented_mass = float(np.sum(probabilities))
    if not np.isfinite(represented_mass) or represented_mass <= 0.0:
        raise ValueError("lognormal kernel has no resolvable mass inside bin edges")
    return probabilities / represented_mass


def build_size_observation_matrix(
    latent_bin_edges_nm: Sequence[float] | np.ndarray,
    observed_bin_edges_nm: Sequence[float] | np.ndarray | None = None,
    instrument_log_diameter_sd: float = 0.0,
) -> np.ndarray:
    """Build an observed-bin by latent-bin size measurement operator.

    With zero instrument width, each entry is the linear-diameter overlap of
    an observed bin with a latent bin divided by the latent-bin width.  With a
    positive width, a particle at each latent bin's geometric center is
    measured with lognormal diameter error.  Columns may sum below one when
    the observed diameter range does not cover the latent range or measurement
    error scatters particles outside it.
    """

    latent_edges = _validate_bin_edges(latent_bin_edges_nm, name="latent_bin_edges_nm")
    observed_edges = (
        latent_edges.copy()
        if observed_bin_edges_nm is None
        else _validate_bin_edges(observed_bin_edges_nm, name="observed_bin_edges_nm")
    )
    if not np.isfinite(instrument_log_diameter_sd) or instrument_log_diameter_sd < 0.0:
        raise ValueError("instrument_log_diameter_sd must be finite and nonnegative")

    latent_lower = latent_edges[:-1]
    latent_upper = latent_edges[1:]
    observed_lower = observed_edges[:-1]
    observed_upper = observed_edges[1:]

    if instrument_log_diameter_sd == 0.0:
        overlap = np.maximum(
            0.0,
            np.minimum(observed_upper[:, None], latent_upper[None, :])
            - np.maximum(observed_lower[:, None], latent_lower[None, :]),
        )
        return overlap / (latent_upper - latent_lower)[None, :]

    latent_centers = np.sqrt(latent_lower * latent_upper)
    sigma = float(instrument_log_diameter_sd)
    upper_z = (
        np.log(observed_upper[:, None]) - np.log(latent_centers[None, :])
    ) / sigma
    lower_z = (
        np.log(observed_lower[:, None]) - np.log(latent_centers[None, :])
    ) / sigma
    matrix = np.clip(ndtr(upper_z) - ndtr(lower_z), 0.0, 1.0)
    if np.any(~np.isfinite(matrix)):
        raise ValueError("instrument observation matrix is non-finite")
    return matrix


def simulate_size_resolved_extracellular_kinetics(
    time_s: Sequence[float] | np.ndarray,
    cumulative_release_per_cell: Mapping[str, Sequence[float] | np.ndarray],
    latent_bin_edges_nm: Sequence[float] | np.ndarray,
    pathway_params: Mapping[str, PathwaySizeKernelParams | Mapping[str, float]],
    kinetics_params: SizeResolvedKineticsParams | Mapping[str, float] | None = None,
    *,
    initial_cell_density_per_ml: float,
    initial_concentration_particles_per_ml: Sequence[float] | np.ndarray | None = None,
    viable_fraction: Sequence[float] | np.ndarray | None = None,
    state_signals: Mapping[str, Sequence[float] | np.ndarray] | None = None,
    observed_bin_edges_nm: Sequence[float] | np.ndarray | None = None,
    observed_background_concentration_particles_per_ml: float
    | Sequence[float]
    | np.ndarray = 0.0,
) -> SizeResolvedKineticsResult:
    """Simulate a reduced, size-resolved extracellular concentration model.

    Cumulative pathway release is interpreted in model units per representative
    cell.  Within each interval, its increment is distributed by a lognormal
    kernel evaluated at the midpoint state.  The resulting constant interval
    source and first-order bin loss are integrated analytically.

    ``initial_concentration_particles_per_ml`` is tracked separately as an
    ambient contribution, while pathway-resolved histories retain attribution
    of newly released particles.  State signals default to zero independently
    for each pathway.  Background may be scalar, one value per observed bin, or
    an ``(n_time, n_observed_bins)`` array.
    """

    t = _validate_time(time_s)
    latent_edges = _validate_bin_edges(latent_bin_edges_nm, name="latent_bin_edges_nm")
    observed_edges = (
        latent_edges.copy()
        if observed_bin_edges_nm is None
        else _validate_bin_edges(observed_bin_edges_nm, name="observed_bin_edges_nm")
    )
    params = _coerce_kinetics_params(kinetics_params)
    kernels = _coerce_pathway_params(pathway_params)
    cumulative = _coerce_cumulative_release(cumulative_release_per_cell, len(t))
    viability = _coerce_viability(viable_fraction, len(t))
    signals = _coerce_state_signals(state_signals, len(t))

    if (
        not np.isfinite(initial_cell_density_per_ml)
        or initial_cell_density_per_ml <= 0.0
    ):
        raise ValueError("initial_cell_density_per_ml must be finite and positive")

    n_time = len(t)
    n_latent = len(latent_edges) - 1
    n_observed = len(observed_edges) - 1
    initial = _coerce_initial_concentration(
        initial_concentration_particles_per_ml, n_latent
    )
    background = _coerce_background(
        observed_background_concentration_particles_per_ml,
        n_time,
        n_observed,
    )

    latent_centers = np.sqrt(latent_edges[:-1] * latent_edges[1:])
    if params.effective_loss_rate_s == 0.0:
        loss_rates = np.zeros(n_latent, dtype=float)
    else:
        loss_rates = (
            params.effective_loss_rate_s
            * (latent_centers / params.loss_reference_diameter_nm)
            ** params.loss_size_exponent
        )
    if np.any(~np.isfinite(loss_rates)) or np.any(loss_rates < 0.0):
        raise ValueError("size-dependent loss rates must be finite and nonnegative")

    ambient_history = np.zeros((n_time, n_latent), dtype=float)
    ambient_history[0] = initial
    pathway_history = {
        pathway: np.zeros((n_time, n_latent), dtype=float) for pathway in SIZE_PATHWAYS
    }
    kernel_history = {
        pathway: np.zeros((n_time - 1, n_latent), dtype=float)
        for pathway in SIZE_PATHWAYS
    }

    for interval in range(1, n_time):
        dt = float(t[interval] - t[interval - 1])
        decay = np.exp(-loss_rates * dt)
        input_factor = np.ones_like(loss_rates)
        positive_loss = loss_rates > 0.0
        input_factor[positive_loss] = -np.expm1(-loss_rates[positive_loss] * dt) / (
            loss_rates[positive_loss] * dt
        )

        ambient_history[interval] = ambient_history[interval - 1] * decay
        midpoint_viability = 0.5 * (viability[interval - 1] + viability[interval])

        for pathway in SIZE_PATHWAYS:
            pathway_parameter = kernels[pathway]
            midpoint_state = 0.5 * (
                signals[pathway][interval - 1] + signals[pathway][interval]
            )
            shifted_median = pathway_parameter.median_diameter_nm * np.exp(
                pathway_parameter.state_shift_coefficient * midpoint_state
            )
            if not np.isfinite(shifted_median) or shifted_median <= 0.0:
                raise ValueError(
                    f"state-shifted median for {pathway} must be finite and positive"
                )
            probability = lognormal_bin_probabilities(
                latent_edges,
                float(shifted_median),
                pathway_parameter.geometric_sd,
            )
            kernel_history[pathway][interval - 1] = probability

            release_increment = max(
                float(
                    cumulative[pathway][interval] - cumulative[pathway][interval - 1]
                ),
                0.0,
            )
            source_increment = (
                release_increment
                * float(initial_cell_density_per_ml)
                * pathway_parameter.source_scale_particles_per_model_unit
                * midpoint_viability
                * probability
            )
            pathway_history[pathway][interval] = (
                pathway_history[pathway][interval - 1] * decay
                + source_increment * input_factor
            )

    total_latent = ambient_history.copy()
    for pathway in SIZE_PATHWAYS:
        total_latent += pathway_history[pathway]

    observation_matrix = build_size_observation_matrix(
        latent_edges,
        observed_edges,
        params.instrument_log_diameter_sd,
    )
    observed = (
        params.assay_recovery_fraction * (total_latent @ observation_matrix.T)
        + background
    )
    if (
        np.any(~np.isfinite(total_latent))
        or np.any(total_latent < 0.0)
        or np.any(~np.isfinite(observed))
        or np.any(observed < 0.0)
    ):
        raise ValueError("size-resolved simulation produced invalid concentrations")

    return SizeResolvedKineticsResult(
        time_s=t.copy(),
        latent_bin_edges_nm=latent_edges.copy(),
        observed_bin_edges_nm=observed_edges.copy(),
        ambient_latent_concentration_particles_per_ml=ambient_history,
        pathway_latent_concentration_particles_per_ml=pathway_history,
        total_latent_concentration_particles_per_ml=total_latent,
        observed_concentration_particles_per_ml=observed,
        kernel_probabilities=kernel_history,
        loss_rates_s=loss_rates,
        observation_matrix=observation_matrix,
    )


def _validate_time(values: Sequence[float] | np.ndarray) -> np.ndarray:
    time = np.asarray(values, dtype=float)
    if time.ndim != 1 or len(time) < 2 or np.any(~np.isfinite(time)):
        raise ValueError(
            "time_s must be a finite one-dimensional array with at least two points"
        )
    if np.any(np.diff(time) <= 0.0):
        raise ValueError("time_s must be strictly increasing")
    return time


def _validate_bin_edges(
    values: Sequence[float] | np.ndarray, *, name: str
) -> np.ndarray:
    edges = np.asarray(values, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or np.any(~np.isfinite(edges)):
        raise ValueError(
            f"{name} must be a finite one-dimensional array with at least two edges"
        )
    if np.any(edges <= 0.0):
        raise ValueError(f"{name} must be strictly positive for log-diameter models")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing with nonzero widths")
    return edges


def _coerce_pathway_params(
    values: Mapping[str, PathwaySizeKernelParams | Mapping[str, float]],
) -> dict[str, PathwaySizeKernelParams]:
    if not isinstance(values, Mapping):
        raise TypeError("pathway_params must be a mapping")
    missing = set(SIZE_PATHWAYS) - set(values)
    unknown = set(values) - set(SIZE_PATHWAYS)
    if missing or unknown:
        raise ValueError(
            f"pathway_params keys must be exactly {SIZE_PATHWAYS}; "
            f"missing={sorted(missing)}, unknown={sorted(unknown)}"
        )
    coerced: dict[str, PathwaySizeKernelParams] = {}
    for pathway in SIZE_PATHWAYS:
        value = values[pathway]
        if isinstance(value, PathwaySizeKernelParams):
            coerced[pathway] = value
        elif isinstance(value, Mapping):
            coerced[pathway] = PathwaySizeKernelParams(**dict(value))
        else:
            raise TypeError(
                f"pathway_params[{pathway!r}] must be a parameter object or mapping"
            )
    return coerced


def _coerce_kinetics_params(
    values: SizeResolvedKineticsParams | Mapping[str, float] | None,
) -> SizeResolvedKineticsParams:
    if values is None:
        return SizeResolvedKineticsParams()
    if isinstance(values, SizeResolvedKineticsParams):
        return values
    if isinstance(values, Mapping):
        return SizeResolvedKineticsParams(**dict(values))
    raise TypeError("kinetics_params must be a parameter object or mapping")


def _coerce_cumulative_release(
    values: Mapping[str, Sequence[float] | np.ndarray], expected_length: int
) -> dict[str, np.ndarray]:
    if not isinstance(values, Mapping):
        raise TypeError("cumulative_release_per_cell must be a mapping")
    missing = set(SIZE_PATHWAYS) - set(values)
    unknown = set(values) - set(SIZE_PATHWAYS)
    if missing or unknown:
        raise ValueError(
            f"cumulative release keys must be exactly {SIZE_PATHWAYS}; "
            f"missing={sorted(missing)}, unknown={sorted(unknown)}"
        )
    result: dict[str, np.ndarray] = {}
    for pathway in SIZE_PATHWAYS:
        cumulative = np.asarray(values[pathway], dtype=float)
        if (
            cumulative.shape != (expected_length,)
            or np.any(~np.isfinite(cumulative))
            or np.any(cumulative < 0.0)
        ):
            raise ValueError(
                f"cumulative release for {pathway} has invalid shape or values"
            )
        differences = np.diff(cumulative)
        if np.any(differences < -1.0e-10):
            raise ValueError(f"cumulative release for {pathway} must be non-decreasing")
        result[pathway] = cumulative.copy()
    return result


def _coerce_viability(
    values: Sequence[float] | np.ndarray | None, expected_length: int
) -> np.ndarray:
    if values is None:
        return np.ones(expected_length, dtype=float)
    viability = np.asarray(values, dtype=float)
    if viability.shape != (expected_length,) or np.any(~np.isfinite(viability)):
        raise ValueError("viable_fraction has invalid shape or values")
    if np.any((viability < 0.0) | (viability > 1.0)):
        raise ValueError("viable_fraction must remain between zero and one")
    return viability


def _coerce_state_signals(
    values: Mapping[str, Sequence[float] | np.ndarray] | None,
    expected_length: int,
) -> dict[str, np.ndarray]:
    if values is None:
        return {
            pathway: np.zeros(expected_length, dtype=float) for pathway in SIZE_PATHWAYS
        }
    if not isinstance(values, Mapping):
        raise TypeError("state_signals must be a mapping")
    unknown = set(values) - set(SIZE_PATHWAYS)
    if unknown:
        raise ValueError(f"unknown state-signal pathways: {sorted(unknown)}")
    signals: dict[str, np.ndarray] = {}
    for pathway in SIZE_PATHWAYS:
        signal = np.asarray(
            values.get(pathway, np.zeros(expected_length, dtype=float)),
            dtype=float,
        )
        if signal.shape != (expected_length,) or np.any(~np.isfinite(signal)):
            raise ValueError(f"state signal for {pathway} has invalid shape or values")
        signals[pathway] = signal
    return signals


def _coerce_initial_concentration(
    values: Sequence[float] | np.ndarray | None, n_latent_bins: int
) -> np.ndarray:
    if values is None:
        return np.zeros(n_latent_bins, dtype=float)
    initial = np.asarray(values, dtype=float)
    if (
        initial.shape != (n_latent_bins,)
        or np.any(~np.isfinite(initial))
        or np.any(initial < 0.0)
    ):
        raise ValueError(
            "initial_concentration_particles_per_ml must contain one finite, "
            "nonnegative value per latent bin"
        )
    return initial.copy()


def _coerce_background(
    values: float | Sequence[float] | np.ndarray,
    n_time: int,
    n_observed_bins: int,
) -> np.ndarray:
    background = np.asarray(values, dtype=float)
    if background.ndim == 0:
        background = np.full((n_time, n_observed_bins), float(background), dtype=float)
    elif background.shape == (n_observed_bins,):
        background = np.broadcast_to(
            background[None, :], (n_time, n_observed_bins)
        ).copy()
    elif background.shape == (n_time, n_observed_bins):
        background = background.copy()
    else:
        raise ValueError(
            "observed background must be scalar, one value per observed bin, "
            "or an (n_time, n_observed_bins) array"
        )
    if np.any(~np.isfinite(background)) or np.any(background < 0.0):
        raise ValueError("observed background must be finite and nonnegative")
    return background
