from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(slots=True)
class RemodelingParams:
    K_PS_uM: float = 1.0
    n_PS: float = 2.0
    PS_max: float = 1.0
    K_calpain_uM: float = 2.0
    n_calpain: float = 2.0
    K_annex_uM: float = 0.5
    n_annex: float = 1.5
    tau_repair_s: float = 300.0


# TODO-literature-review: replace generic Hill functions with experimentally grounded Ca2+-remodeling relations.
def compute_remodeling_state(Ca_i: float, params: RemodelingParams) -> dict[str, float]:
    """Compute placeholder Layer 4 remodeling and repair observables from cytosolic calcium."""
    ps_exposure = (Ca_i ** params.n_PS) / (params.K_PS_uM ** params.n_PS + Ca_i ** params.n_PS) * params.PS_max
    calpain = (Ca_i ** params.n_calpain) / (params.K_calpain_uM ** params.n_calpain + Ca_i ** params.n_calpain)
    annexin = (Ca_i ** params.n_annex) / (params.K_annex_uM ** params.n_annex + Ca_i ** params.n_annex)
    actin_disruption = min(1.0, 0.7 * calpain + 0.3 * ps_exposure)
    repair_state = max(0.0, min(1.0, annexin * (1.0 - 0.5 * calpain)))
    return {
        "PS_exposure": float(ps_exposure),
        "calpain_activity": float(calpain),
        "annexin_activity": float(annexin),
        "actin_disruption": float(actin_disruption),
        "repair_state": float(repair_state),
    }



def remodeling_defaults() -> dict:
    """Return remodeling defaults as a dict."""
    return asdict(RemodelingParams())
