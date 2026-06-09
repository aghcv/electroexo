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
    potency_weights: dict[str, float] = field(default_factory=lambda: {"protein": 0.5, "RNA": 0.3, "lipid": 0.2})


# TODO-literature-review: replace proxy enrichment and potency equations with cargo-specific evidence.
def compute_cargo_state(Ca_i: float, ROS: float, ATP: float, params: CargoPotencyParams) -> dict[str, float]:
    """Compute placeholder Layer 6 cargo enrichment and potency metrics."""
    protein = params.protein_enrichment_baseline * (1.0 + params.Ca_protein_coupling * max(Ca_i - 0.1, 0.0))
    rna = params.RNA_enrichment_baseline * (1.0 + params.ROS_RNA_coupling * max(ROS - 0.1, 0.0))
    lipid = params.lipid_enrichment_baseline * (1.0 + params.stress_lipid_coupling * max(1.0 - ATP, 0.0))
    potency = (
        params.potency_weights["protein"] * protein
        + params.potency_weights["RNA"] * rna
        + params.potency_weights["lipid"] * lipid
    )
    return {
        "protein_enrichment": float(protein),
        "RNA_enrichment": float(rna),
        "lipid_enrichment": float(lipid),
        "potency_score": float(potency),
    }



def cargo_potency_defaults() -> dict:
    """Return cargo/potency defaults as a dict."""
    return asdict(CargoPotencyParams())
