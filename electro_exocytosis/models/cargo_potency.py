from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass(slots=True)
class CargoPotencyParams:
    protein_enrichment_baseline: float = 1.0
    RNA_enrichment_baseline: float = 1.0
    lipid_enrichment_baseline: float = 1.0
    Ca_protein_coupling: float = 0.5
    ROS_RNA_coupling: float = 0.3
    stress_lipid_coupling: float = 0.2
    subtype_sEV_weight: float = 1.0
    subtype_mlEV_weight: float = 0.65
    subtype_AB_weight: float = 0.20
    ESCRT_protein_weight: float = 0.30
    rbp_RNA_sorting_weight: float = 0.35
    ceramide_lipid_weight: float = 0.30
    antigen_stress_weight: float = 0.40
    antigen_sorting_weight: float = 0.25
    direct_loading_efficiency: float = 0.35
    direct_loading_leak_fraction: float = 0.10
    recipient_dose_half_max: float = 1.0
    potency_saturation: float = 2.5
    potency_weights: dict[str, float] = field(
        default_factory=lambda: {
            "protein": 0.30,
            "RNA": 0.25,
            "lipid": 0.15,
            "antigen": 0.20,
            "direct_load": 0.10,
        }
    )


def compute_cargo_state(
    Ca_i: float,
    ROS: float,
    ATP: float,
    params: CargoPotencyParams,
    *,
    cumulative_sEV: float = 0.0,
    cumulative_mlEV: float = 0.0,
    cumulative_AB: float = 0.0,
    escrt_signal: float = 0.0,
    ceramide_signal: float = 0.0,
    secretory_bias: float = 0.0,
    direct_loading_drive: float = 0.0,
) -> dict[str, float]:
    """Compute reduced Layer 6 cargo, composition, direct-loading, and potency metrics.

    The form follows the report's computational table: subtype-weighted cargo,
    stress-modulated sorting, antigen enrichment, lipid state, potency mapping,
    and direct EV loading are kept as distinct outputs rather than collapsed into
    one opaque score.
    """
    ca_gate = _clip01(max(float(Ca_i) - 0.1, 0.0))
    ros_gate = _clip01(max(float(ROS) - 0.1, 0.0))
    atp_stress = _clip01(max(1.0 - float(ATP), 0.0))
    escrt = _clip01(float(escrt_signal))
    ceramide = _clip01(float(ceramide_signal))
    secretory = _clip01(float(secretory_bias))

    subtype_total = max(float(cumulative_sEV) + float(cumulative_mlEV) + float(cumulative_AB), 0.0)
    if subtype_total > 0.0:
        subtype_weighted_output = (
            params.subtype_sEV_weight * float(cumulative_sEV)
            + params.subtype_mlEV_weight * float(cumulative_mlEV)
            + params.subtype_AB_weight * float(cumulative_AB)
        ) / subtype_total
    else:
        subtype_weighted_output = 1.0

    protein = params.protein_enrichment_baseline * (
        1.0
        + params.Ca_protein_coupling * ca_gate
        + params.ESCRT_protein_weight * escrt
    )
    rna = params.RNA_enrichment_baseline * (
        1.0
        + params.ROS_RNA_coupling * ros_gate
        + params.rbp_RNA_sorting_weight * secretory
    )
    lipid = params.lipid_enrichment_baseline * (
        1.0
        + params.stress_lipid_coupling * atp_stress
        + params.ceramide_lipid_weight * ceramide
    )
    antigen = 1.0 + params.antigen_stress_weight * ros_gate + params.antigen_sorting_weight * escrt
    direct_loaded = max(
        0.0,
        params.direct_loading_efficiency * max(float(direct_loading_drive), 0.0)
        - params.direct_loading_leak_fraction,
    )
    cargo_vector_score = (
        params.potency_weights.get("protein", 0.0) * protein
        + params.potency_weights.get("RNA", 0.0) * rna
        + params.potency_weights.get("lipid", 0.0) * lipid
        + params.potency_weights.get("antigen", 0.0) * antigen
        + params.potency_weights.get("direct_load", 0.0) * direct_loaded
    )
    dose_gate = subtype_total / (params.recipient_dose_half_max + subtype_total) if subtype_total > 0 else 1.0
    potency = params.potency_saturation * cargo_vector_score * subtype_weighted_output * dose_gate / (
        params.potency_saturation + cargo_vector_score
    )
    return {
        "protein_enrichment": float(protein),
        "RNA_enrichment": float(rna),
        "lipid_enrichment": float(lipid),
        "antigen_enrichment": float(antigen),
        "direct_loaded_cargo": float(direct_loaded),
        "subtype_weighted_cargo": float(subtype_weighted_output),
        "cargo_vector_score": float(cargo_vector_score),
        "potency_score": float(potency),
    }


def _clip01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))



def cargo_potency_defaults() -> dict:
    """Return cargo/potency defaults as a dict."""
    return asdict(CargoPotencyParams())
