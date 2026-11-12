"""Sensitivity, local identifiability, and experiment-priority utilities."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import qmc

from electro_exocytosis.io.readers import load_default_parameters
from electro_exocytosis.models.ev_release import (
    coerce_ev_release_params,
    get_ev_initial_conditions,
    get_ev_state_names,
)
from electro_exocytosis.models.ion_transport import (
    coerce_ion_transport_params,
    get_ion_initial_conditions,
    get_ion_state_names,
)
from electro_exocytosis.parameter_registry import build_parameter_snapshot

ArrayModel = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True, slots=True)
class TargetDefinition:
    """One scalar quantity of interest tied to a condition and time window."""

    identifier: str
    output: str
    condition: str
    time_window: str
    unit: str
    description: str


_STATE_METADATA: dict[str, tuple[str, str, str, float, float, str]] = {
    "Ca_i": (
        "ion_transport",
        "calcium influx and stores",
        "Cytosolic free calcium",
        0.0,
        20.0,
        "uM",
    ),
    "Ca_ER": (
        "ion_transport",
        "calcium influx and stores",
        "Endoplasmic-reticulum calcium",
        0.0,
        2000.0,
        "uM",
    ),
    "Ca_mito": (
        "ion_transport",
        "mitochondrial calcium",
        "Mitochondrial calcium",
        0.0,
        20.0,
        "uM",
    ),
    "mitochondrial_potential": (
        "ion_transport",
        "mitochondrial calcium",
        "Normalized mitochondrial potential",
        0.0,
        1.0,
        "dimensionless",
    ),
    "ROS": (
        "ion_transport",
        "metabolic stress",
        "Normalized reactive-oxygen state",
        0.0,
        10.0,
        "model units",
    ),
    "ATP": (
        "ion_transport",
        "metabolic stress",
        "Normalized ATP state",
        0.0,
        1.0,
        "dimensionless",
    ),
    "Na_i": (
        "ion_transport",
        "monovalent ions and osmosis",
        "Intracellular sodium",
        0.0,
        200.0,
        "mM",
    ),
    "K_i": (
        "ion_transport",
        "monovalent ions and osmosis",
        "Intracellular potassium",
        0.0,
        200.0,
        "mM",
    ),
    "Cl_i": (
        "ion_transport",
        "monovalent ions and osmosis",
        "Intracellular chloride",
        0.0,
        200.0,
        "mM",
    ),
    "osmotic_stress": (
        "ion_transport",
        "monovalent ions and osmosis",
        "Normalized osmotic stress",
        0.0,
        10.0,
        "model units",
    ),
    "MVB_pool": (
        "ev_release",
        "endosomal sEV pathway",
        "Multivesicular-body pool",
        0.0,
        np.inf,
        "model units",
    ),
    "ILV_load": (
        "ev_release",
        "endosomal sEV pathway",
        "Intraluminal-vesicle load",
        0.0,
        np.inf,
        "model units",
    ),
    "docked_MVB_pool": (
        "ev_release",
        "endosomal sEV pathway",
        "Docked MVB pool",
        0.0,
        np.inf,
        "model units",
    ),
    "budding_pool": (
        "ev_release",
        "plasma-membrane mlEV pathway",
        "Plasma-membrane budding pool",
        0.0,
        np.inf,
        "model units",
    ),
    "apoptotic_commitment": (
        "ev_release",
        "apoptotic-body pathway",
        "Apoptotic commitment",
        0.0,
        np.inf,
        "model units",
    ),
    "cumulative_sEV": (
        "simulation",
        "cumulative release",
        "Cumulative small-EV release",
        0.0,
        np.inf,
        "model units cell^-1",
    ),
    "cumulative_mlEV": (
        "simulation",
        "cumulative release",
        "Cumulative medium/large-EV release",
        0.0,
        np.inf,
        "model units cell^-1",
    ),
    "cumulative_AB": (
        "simulation",
        "cumulative release",
        "Cumulative apoptotic-body release",
        0.0,
        np.inf,
        "model units cell^-1",
    ),
    "damage": (
        "injury_quality",
        "cell injury and fate",
        "Integrated cellular damage",
        0.0,
        np.inf,
        "model units",
    ),
    "extracellular_sEV_stock": (
        "extracellular_kinetics",
        "extracellular stocks",
        "Small-EV extracellular particle stock",
        0.0,
        np.inf,
        "particles",
    ),
    "extracellular_mlEV_stock": (
        "extracellular_kinetics",
        "extracellular stocks",
        "Medium/large-EV extracellular particle stock",
        0.0,
        np.inf,
        "particles",
    ),
    "extracellular_AB_stock": (
        "extracellular_kinetics",
        "extracellular stocks",
        "Apoptotic-body extracellular particle stock",
        0.0,
        np.inf,
        "particles",
    ),
    "extracellular_medium_volume": (
        "extracellular_kinetics",
        "sampling and assay",
        "Culture-medium volume changed by sampling events",
        0.0,
        np.inf,
        "mL",
    ),
    "ambient_latent_size_stock[*]": (
        "ev_size_observation",
        "size-dependent loss",
        "Ambient extracellular concentration in each latent size bin",
        0.0,
        np.inf,
        "particles mL^-1",
    ),
    "sEV_latent_size_stock[*]": (
        "ev_size_observation",
        "pathway-to-size kernel",
        "sEV-attributed concentration in each latent size bin",
        0.0,
        np.inf,
        "particles mL^-1",
    ),
    "mlEV_latent_size_stock[*]": (
        "ev_size_observation",
        "pathway-to-size kernel",
        "m/lEV-attributed concentration in each latent size bin",
        0.0,
        np.inf,
        "particles mL^-1",
    ),
    "AB_latent_size_stock[*]": (
        "ev_size_observation",
        "pathway-to-size kernel",
        "AB-attributed concentration in each latent size bin",
        0.0,
        np.inf,
        "particles mL^-1",
    ),
}


def build_state_inventory() -> pd.DataFrame:
    """Return dynamic states separately from model parameters."""

    defaults = load_default_parameters()
    ion = coerce_ion_transport_params(defaults)
    ev = coerce_ev_release_params(defaults)
    initial = dict(
        zip(get_ion_state_names(), get_ion_initial_conditions(ion), strict=True)
    )
    initial.update(
        dict(zip(get_ev_state_names(), get_ev_initial_conditions(ev), strict=True))
    )
    initial.update(
        cumulative_sEV=0.0,
        cumulative_mlEV=0.0,
        cumulative_AB=0.0,
        damage=0.0,
        extracellular_sEV_stock=0.0,
        extracellular_mlEV_stock=0.0,
        extracellular_AB_stock=0.0,
        extracellular_medium_volume=1.0,
        **{
            "ambient_latent_size_stock[*]": 0.0,
            "sEV_latent_size_stock[*]": 0.0,
            "mlEV_latent_size_stock[*]": 0.0,
            "AB_latent_size_stock[*]": 0.0,
        },
    )
    rows = []
    for name, (
        module,
        submodule,
        description,
        lower,
        upper,
        unit,
    ) in _STATE_METADATA.items():
        rows.append(
            {
                "state_id": f"state.{module}.{name}",
                "module": module,
                "submodule": submodule,
                "name": name,
                "symbol": name,
                "description": description,
                "unit": unit,
                "initial_condition": initial[name],
                "biological_lower_bound": lower,
                "biological_upper_bound": upper if np.isfinite(upper) else pd.NA,
                "optimizer_initial": pd.NA,
                "optimizer_lower_bound": pd.NA,
                "optimizer_upper_bound": pd.NA,
                "fit_status": "state_not_parameter",
                "fitted_value": pd.NA,
                "evidence_source": "model equation and canonical default",
                "evidence_confidence": "model-defined; biological bound provisional",
                "observability_notes": "Dynamic state; requires a state-specific assay or proxy and is not ranked as a parameter.",
            }
        )
    return pd.DataFrame(rows)


def build_parameter_inventory(
    *,
    fit_parameters: pd.DataFrame | None = None,
    selected_fit_variant: str | None = None,
) -> pd.DataFrame:
    """Enrich the canonical snapshot without replacing its numerical truth."""

    snapshot = build_parameter_snapshot(
        fit_parameters=(
            fit_parameters.drop(columns=["module"], errors="ignore")
            if fit_parameters is not None
            else None
        ),
        fit_module="ev_size_observation" if fit_parameters is not None else None,
        selected_fit_variant=selected_fit_variant,
    ).copy()
    biological_lower: list[Any] = []
    biological_upper: list[Any] = []
    screening_lower: list[Any] = []
    screening_upper: list[Any] = []
    range_sources: list[str] = []
    for row in snapshot.itertuples(index=False):
        value = row.fit_final if not pd.isna(row.fit_final) else row.effective_value
        if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(
            float(value)
        ):
            lo, hi, source = default_screening_bounds(float(value), str(row.parameter))
            screening_lower.append(lo)
            screening_upper.append(hi)
            # Unknown literature domains remain blank; do not relabel robustness ranges.
            biological_lower.append(pd.NA)
            biological_upper.append(pd.NA)
            range_sources.append(source)
        else:
            screening_lower.append(pd.NA)
            screening_upper.append(pd.NA)
            biological_lower.append(pd.NA)
            biological_upper.append(pd.NA)
            range_sources.append("not_numeric_or_not_curated")
    snapshot.insert(
        0, "parameter_id", "parameter." + snapshot["qualified_name"].astype(str)
    )
    snapshot["name_symbol"] = snapshot["parameter"]
    snapshot["description"] = snapshot["submodule"].astype(str) + " model parameter"
    snapshot["biological_lower_bound"] = biological_lower
    snapshot["biological_upper_bound"] = biological_upper
    snapshot["screening_lower_bound"] = screening_lower
    snapshot["screening_upper_bound"] = screening_upper
    snapshot["screening_range_source"] = range_sources
    snapshot["evidence_source"] = snapshot["source"]
    snapshot["evidence_confidence"] = np.where(
        snapshot["runtime_status"].eq("declared_not_consumed"),
        "inactive in current equations",
        "preliminary",
    )
    snapshot["observability_notes"] = np.where(
        snapshot["fit_status"].eq("fitted"),
        "Engaged by the current size-resolved fit; identifiability still requires diagnostics.",
        "Not directly estimated by the current longitudinal size-resolved dataset.",
    )
    return snapshot


def default_screening_bounds(value: float, name: str) -> tuple[float, float, str]:
    """Return a conservative *provisional robustness* range around a default.

    These ranges are not literature-derived biological validity intervals.
    Positive quantities span 0.8 to 1.2 times the default; fractions and
    efficiencies are additionally clipped to [0, 1].  Zero-valued optional
    rates use [0, 1e-4] s^-1 and other zero-valued coefficients use [-0.1,
    0.1].  Every generated registry row labels this provenance explicitly.
    """

    lower_name = name.lower()
    bounded_fraction = any(
        token in lower_name
        for token in ("fraction", "efficiency", "purity", "viability", "probability")
    )
    if value > 0.0:
        lower, upper = 0.8 * value, 1.2 * value
        if bounded_fraction:
            lower, upper = max(0.0, lower), min(1.0, upper)
    elif value < 0.0:
        lower, upper = 1.2 * value, 0.8 * value
    elif lower_name.endswith("_rate_s"):
        lower, upper = 0.0, 1.0e-4
    elif any(
        token in lower_name
        for token in (
            "fraction",
            "efficiency",
            "background",
            "density",
            "concentration",
        )
    ):
        lower, upper = 0.0, 0.1
    else:
        lower, upper = -0.1, 0.1
    if not lower < upper:
        lower = max(0.0, value * 0.5)
        upper = value + max(abs(value) * 0.1, 1.0e-9)
    return float(lower), float(upper), "provisional_default_centered_robustness_range"


def morris_trajectories(
    n_parameters: int,
    trajectories: int,
    *,
    levels: int = 6,
    seed: int = 0,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Create randomized one-at-a-time Morris paths on the unit hypercube.

    Returns ``(points, changed_parameter, signed_step)`` per trajectory.
    The construction uses a grid step of ``levels / (2*(levels-1))`` and a
    random orientation/order.  It is suitable for elementary-effect screening
    and intentionally does not claim an optimized trajectory design.
    """

    if n_parameters < 1 or trajectories < 1 or levels < 4:
        raise ValueError("n_parameters, trajectories, and levels must be valid")
    rng = np.random.default_rng(seed)
    delta = levels / (2.0 * (levels - 1.0))
    grid = np.arange(levels, dtype=float) / (levels - 1.0)
    paths: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for _ in range(trajectories):
        orientation = rng.choice(np.array([-1.0, 1.0]), size=n_parameters)
        start = np.empty(n_parameters, dtype=float)
        for j, sign in enumerate(orientation):
            candidates = grid[grid <= 1.0 - delta] if sign > 0 else grid[grid >= delta]
            start[j] = float(rng.choice(candidates))
        order = rng.permutation(n_parameters)
        points = [start.copy()]
        signed_steps: list[float] = []
        current = start.copy()
        for index in order:
            step = orientation[index] * delta
            current = current.copy()
            current[index] += step
            points.append(current)
            signed_steps.append(float(step))
        paths.append((np.asarray(points), order, np.asarray(signed_steps)))
    return paths


def morris_screen(
    model: ArrayModel,
    lower: Sequence[float],
    upper: Sequence[float],
    *,
    trajectories: int = 8,
    levels: int = 6,
    seed: int = 0,
    parameter_ids: Sequence[str] | None = None,
    target_ids: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate Morris elementary effects and return summary and raw tables."""

    lo = np.asarray(lower, dtype=float)
    hi = np.asarray(upper, dtype=float)
    if lo.ndim != 1 or hi.shape != lo.shape or np.any(~np.isfinite(lo + hi)):
        raise ValueError("bounds must be finite one-dimensional arrays")
    if np.any(hi <= lo):
        raise ValueError("every upper bound must exceed its lower bound")
    pids = list(parameter_ids or [f"p{j}" for j in range(len(lo))])
    raw: list[dict[str, Any]] = []
    inferred_targets: list[str] | None = list(target_ids) if target_ids else None
    for trajectory_index, (unit_points, order, steps) in enumerate(
        morris_trajectories(len(lo), trajectories, levels=levels, seed=seed)
    ):
        physical = lo + unit_points * (hi - lo)
        outputs = np.asarray([model(point) for point in physical], dtype=float)
        if outputs.ndim == 1:
            outputs = outputs[:, None]
        if np.any(~np.isfinite(outputs)):
            raise FloatingPointError("model returned non-finite Morris outputs")
        if inferred_targets is None:
            inferred_targets = [f"y{k}" for k in range(outputs.shape[1])]
        if len(inferred_targets) != outputs.shape[1]:
            raise ValueError("target_ids do not match model output width")
        for step_index, (parameter_index, signed_step) in enumerate(
            zip(order, steps, strict=True)
        ):
            effect = (outputs[step_index + 1] - outputs[step_index]) / signed_step
            for target_index, target in enumerate(inferred_targets):
                raw.append(
                    {
                        "parameter_id": pids[int(parameter_index)],
                        "target_id": target,
                        "trajectory": trajectory_index,
                        "elementary_effect": float(effect[target_index]),
                    }
                )
    raw_frame = pd.DataFrame(raw)
    summary = (
        raw_frame.groupby(["parameter_id", "target_id"], sort=False)[
            "elementary_effect"
        ]
        .agg(
            morris_mu="mean",
            morris_mu_star=lambda x: float(np.mean(np.abs(x))),
            morris_sigma=lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else np.nan,
            n_effects="count",
        )
        .reset_index()
    )
    return summary, raw_frame


def sobol_saltelli(
    model: ArrayModel,
    lower: Sequence[float],
    upper: Sequence[float],
    *,
    base_samples: int = 256,
    seed: int = 0,
    parameter_ids: Sequence[str] | None = None,
    target_ids: Sequence[str] | None = None,
    bootstrap_replicates: int = 200,
) -> pd.DataFrame:
    """Estimate Saltelli/Jansen first- and total-order Sobol indices.

    The returned estimates are Monte Carlo estimates (not analytic/exact
    indices). ``base_samples`` model rows are used for each A, B, and A_Bi
    block, for ``base_samples * (p + 2)`` total evaluations.
    """

    lo, hi = np.asarray(lower, float), np.asarray(upper, float)
    p = len(lo)
    if base_samples < 16 or hi.shape != lo.shape or np.any(hi <= lo):
        raise ValueError("invalid Sobol configuration")
    sampler = qmc.Sobol(d=2 * p, scramble=True, seed=seed)
    exponent = int(np.ceil(np.log2(base_samples)))
    design = sampler.random_base2(exponent)[:base_samples]
    a = lo + design[:, :p] * (hi - lo)
    b = lo + design[:, p:] * (hi - lo)
    fa = np.asarray([model(x) for x in a], float)
    fb = np.asarray([model(x) for x in b], float)
    if fa.ndim == 1:
        fa, fb = fa[:, None], fb[:, None]
    pids = list(parameter_ids or [f"p{j}" for j in range(p)])
    tids = list(target_ids or [f"y{k}" for k in range(fa.shape[1])])
    variance = np.var(np.vstack([fa, fb]), axis=0, ddof=1)
    rows: list[dict[str, Any]] = []
    bootstrap_rng = np.random.default_rng(seed + 7919)
    for j in range(p):
        hybrid = a.copy()
        hybrid[:, j] = b[:, j]
        fab = np.asarray([model(x) for x in hybrid], float)
        if fab.ndim == 1:
            fab = fab[:, None]
        for k, target in enumerate(tids):
            if variance[k] <= np.finfo(float).eps:
                s1 = st = np.nan
            else:
                s1 = np.mean(fb[:, k] * (fab[:, k] - fa[:, k])) / variance[k]
                st = 0.5 * np.mean((fa[:, k] - fab[:, k]) ** 2) / variance[k]
            s1_boot: list[float] = []
            st_boot: list[float] = []
            if variance[k] > np.finfo(float).eps and bootstrap_replicates > 0:
                for _ in range(bootstrap_replicates):
                    sample = bootstrap_rng.integers(0, base_samples, size=base_samples)
                    sampled_variance = np.var(
                        np.r_[fa[sample, k], fb[sample, k]], ddof=1
                    )
                    if sampled_variance <= np.finfo(float).eps:
                        continue
                    s1_boot.append(
                        float(
                            np.mean(fb[sample, k] * (fab[sample, k] - fa[sample, k]))
                            / sampled_variance
                        )
                    )
                    st_boot.append(
                        float(
                            0.5
                            * np.mean((fa[sample, k] - fab[sample, k]) ** 2)
                            / sampled_variance
                        )
                    )
            rows.append(
                {
                    "parameter_id": pids[j],
                    "target_id": target,
                    "sobol_S1": float(s1),
                    "sobol_ST": float(st),
                    "sobol_S1_ci_low": float(np.quantile(s1_boot, 0.025))
                    if s1_boot
                    else np.nan,
                    "sobol_S1_ci_high": float(np.quantile(s1_boot, 0.975))
                    if s1_boot
                    else np.nan,
                    "sobol_ST_ci_low": float(np.quantile(st_boot, 0.025))
                    if st_boot
                    else np.nan,
                    "sobol_ST_ci_high": float(np.quantile(st_boot, 0.975))
                    if st_boot
                    else np.nan,
                    "base_samples": base_samples,
                    "bootstrap_replicates": bootstrap_replicates,
                    "estimator": "Saltelli S1 / Jansen ST Monte Carlo",
                }
            )
    return pd.DataFrame(rows)


def local_fisher_diagnostics(
    residual: Callable[[np.ndarray], np.ndarray],
    estimate: Sequence[float],
    lower: Sequence[float],
    upper: Sequence[float],
    *,
    parameter_ids: Sequence[str] | None = None,
    relative_step: float = 1.0e-4,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    """Finite-difference local Fisher diagnostics for practical identifiability."""

    x, lo, hi = map(lambda a: np.asarray(a, float), (estimate, lower, upper))
    if not (x.shape == lo.shape == hi.shape) or np.any(hi <= lo):
        raise ValueError("estimate and bounds must have matching valid shapes")
    baseline = np.asarray(residual(x), float)
    jacobian = np.empty((len(baseline), len(x)), dtype=float)
    for j in range(len(x)):
        step = max(relative_step * max(abs(x[j]), hi[j] - lo[j], 1.0), 1.0e-8)
        plus, minus = x.copy(), x.copy()
        plus[j], minus[j] = min(hi[j], x[j] + step), max(lo[j], x[j] - step)
        denominator = plus[j] - minus[j]
        jacobian[:, j] = (residual(plus) - residual(minus)) / denominator
    scale = np.maximum(np.abs(x), 0.1 * (hi - lo))
    scaled = jacobian * scale[None, :]
    fisher = scaled.T @ scaled
    covariance = np.linalg.pinv(fisher, rcond=1.0e-10)
    sd = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = covariance / np.outer(sd, sd)
    pids = list(parameter_ids or [f"p{j}" for j in range(len(x))])
    information = np.diag(fisher)
    information_floor = max(float(np.max(information)) * 1.0e-10, 1.0e-14)
    condition = float(np.linalg.cond(fisher)) if np.any(fisher) else np.inf
    rows = []
    for j, pid in enumerate(pids):
        relative_sd = sd[j] if information[j] >= information_floor else np.inf
        rows.append(
            {
                "parameter_id": pid,
                "scaled_fisher_information": float(information[j]),
                "local_relative_sd": float(relative_sd),
                "identifiability_status": (
                    "practically_nonidentifiable"
                    if not np.isfinite(relative_sd)
                    else "locally_identifiable"
                    if relative_sd < 0.5
                    else "weak"
                    if relative_sd < 2.0
                    else "practically_nonidentifiable"
                ),
            }
        )
    correlation = pd.DataFrame(corr, index=pids, columns=pids)
    singular = np.linalg.svd(scaled, compute_uv=False)
    summary = {
        "fisher_condition_number": condition,
        "effective_rank": float(np.linalg.matrix_rank(scaled, tol=1.0e-8)),
        "largest_scaled_jacobian_singular_value": float(singular[0])
        if len(singular)
        else 0.0,
        "smallest_scaled_jacobian_singular_value": float(singular[-1])
        if len(singular)
        else 0.0,
    }
    return pd.DataFrame(rows), correlation, summary


def priority_scores(
    components: pd.DataFrame,
    *,
    weights: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    """Compute a configurable geometric decision metric from named components."""

    columns = (
        "sensitivity",
        "uncertainty",
        "identifiability_gap",
        "measurability",
        "feasibility",
    )
    missing = set(columns) - set(components.columns)
    if missing:
        raise ValueError(f"missing priority components: {sorted(missing)}")
    chosen = {column: 1.0 for column in columns}
    chosen.update(dict(weights or {}))
    result = components.copy()
    score = np.ones(len(result), dtype=float)
    total_weight = 0.0
    for column in columns:
        values = np.clip(result[column].to_numpy(float), 0.0, 1.0)
        weight = float(chosen[column])
        score *= np.maximum(values, 1.0e-12) ** weight
        total_weight += weight
    result["decision_score"] = score ** (1.0 / max(total_weight, 1.0e-12))
    result["decision_metric"] = "weighted geometric mean of normalized components"
    return result.sort_values("decision_score", ascending=False, ignore_index=True)
