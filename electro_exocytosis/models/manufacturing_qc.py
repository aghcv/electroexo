from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(slots=True)
class ManufacturingParams:
    cell_count: float = 1.0e6
    harvest_time_h: float = 24.0
    isolation_efficiency: float = 0.3
    isolation_method_factor: float = 1.0
    purity_factor: float = 0.7
    protein_contamination_factor: float = 0.15
    batch_consistency: float = 0.9
    batch_variability_fraction: float = 0.12
    scalability_factor: float = 1.0
    potency_weight: float = 0.45
    yield_weight: float = 0.30
    purity_weight: float = 0.15
    viability_weight: float = 0.10


# TODO-literature-review: replace fixed isolation and purity transforms with workflow-specific QC models.
def compute_manufacturing_outputs(
    cumulative_sEV: float,
    cumulative_mlEV: float,
    cumulative_AB: float,
    viability_fraction: float,
    params: ManufacturingParams,
    *,
    potency_score: float = 1.0,
    purity_score: float = 1.0,
) -> dict[str, float]:
    """Compute Layer 8 yield, recovery, purity, batch, and optimization outputs."""
    recovery = _clip01(params.isolation_efficiency * params.isolation_method_factor)
    scale = max(float(params.scalability_factor), 0.0)
    viability = _clip01(float(viability_fraction))
    isolated_sev = float(cumulative_sEV) * recovery * scale
    isolated_mlev = float(cumulative_mlEV) * recovery * scale
    isolated_ab = float(cumulative_AB) * recovery * max(0.5, viability)
    total_particles = isolated_sev + isolated_mlev + isolated_ab
    cell_normalized_yield = total_particles / max(float(params.cell_count), 1.0)
    harvest_rate = total_particles / max(float(params.harvest_time_h), 1e-12)
    process_purity = _clip01(
        float(purity_score)
        * params.purity_factor
        * viability
        / (1.0 + params.protein_contamination_factor * max(isolated_ab, 0.0))
    )
    batch_adjusted_yield = total_particles * _clip01(params.batch_consistency) * (1.0 - _clip01(params.batch_variability_fraction))
    yield_score = total_particles / (1.0 + total_particles)
    objective_score = (
        params.potency_weight * _clip01(float(potency_score))
        + params.yield_weight * yield_score
        + params.purity_weight * process_purity
        + params.viability_weight * viability
    )
    return {
        "isolated_yield_sEV": float(isolated_sev),
        "isolated_yield_mlEV": float(isolated_mlev),
        "isolated_yield_AB": float(isolated_ab),
        "total_measured_particles": float(total_particles),
        "cell_normalized_yield": float(cell_normalized_yield),
        "harvest_rate": float(harvest_rate),
        "process_recovery": float(recovery),
        "purity_score": float(process_purity),
        "batch_consistency": float(params.batch_consistency),
        "batch_adjusted_yield": float(batch_adjusted_yield),
        "optimization_objective": float(objective_score),
    }


def _clip01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))



def manufacturing_defaults() -> dict:
    """Return manufacturing defaults as a dict."""
    return asdict(ManufacturingParams())
