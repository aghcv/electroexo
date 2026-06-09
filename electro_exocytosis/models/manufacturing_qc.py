from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(slots=True)
class ManufacturingParams:
    isolation_efficiency: float = 0.3
    purity_factor: float = 0.7
    batch_consistency: float = 0.9
    scalability_factor: float = 1.0


# TODO-literature-review: replace fixed isolation and purity transforms with workflow-specific QC models.
def compute_manufacturing_outputs(
    cumulative_sEV: float,
    cumulative_mlEV: float,
    cumulative_AB: float,
    viability_fraction: float,
    params: ManufacturingParams,
) -> dict[str, float]:
    """Compute placeholder Layer 8 manufacturing and QC outputs."""
    isolated_sev = cumulative_sEV * params.isolation_efficiency * params.scalability_factor
    isolated_mlev = cumulative_mlEV * params.isolation_efficiency * params.scalability_factor
    isolated_ab = cumulative_AB * params.isolation_efficiency * max(0.5, viability_fraction)
    total_particles = isolated_sev + isolated_mlev + isolated_ab
    purity_score = params.purity_factor * viability_fraction
    return {
        "isolated_yield_sEV": float(isolated_sev),
        "isolated_yield_mlEV": float(isolated_mlev),
        "isolated_yield_AB": float(isolated_ab),
        "total_measured_particles": float(total_particles),
        "purity_score": float(purity_score),
        "batch_consistency": float(params.batch_consistency),
    }



def manufacturing_defaults() -> dict:
    """Return manufacturing defaults as a dict."""
    return asdict(ManufacturingParams())
