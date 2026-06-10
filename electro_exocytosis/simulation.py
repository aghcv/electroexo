from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import pandas as pd

from electro_exocytosis.config import SimulationScenario
from electro_exocytosis.io.readers import load_default_parameters, merge_parameters
from electro_exocytosis.models.cargo_potency import CargoPotencyParams, compute_cargo_state
from electro_exocytosis.models.cell_state import apply_cell_state_modifiers
from electro_exocytosis.models.dosimetry import compute_dosimetry
from electro_exocytosis.models.electrodynamics import compute_electrodynamics_state
from electro_exocytosis.models.injury_quality import InjuryParams, compute_quality_gate
from electro_exocytosis.models.manufacturing_qc import ManufacturingParams, compute_manufacturing_outputs
from electro_exocytosis.models.pulse import compute_pulse_descriptors
from electro_exocytosis.models.remodeling_repair import RemodelingParams, compute_remodeling_state
from electro_exocytosis.numerics.multiscale import MultiscaleScheduler
from electro_exocytosis.numerics.solvers import solve_ode_system


@dataclass(slots=True)
class SimulationResult:
    scenario_name: str
    mode: str
    t_array: np.ndarray
    state_timeseries: pd.DataFrame
    ev_timeseries: pd.DataFrame
    summary: dict[str, Any]
    parameters_used: dict[str, Any]
    warnings: list[str]


class Simulation:
    def __init__(self, scenario: SimulationScenario, params_override: dict | None = None):
        self.scenario = scenario
        self.params_nested = merge_parameters(load_default_parameters(), params_override or {})
        self.params_flat = self._flatten_parameters(self.params_nested)
        self.params_flat = apply_cell_state_modifiers(self.params_flat, self.scenario.cell_state)
        self.warnings: list[str] = [
            "All scientific modules use placeholder equations/parameters in v0.1.0.",
        ]

    @staticmethod
    def _flatten_parameters(params_nested: dict[str, Any]) -> dict[str, Any]:
        flat: dict[str, Any] = {}
        for value in params_nested.values():
            if isinstance(value, dict):
                for inner_key, inner_value in value.items():
                    flat[inner_key] = inner_value
        return flat

    def run(self) -> SimulationResult:
        descriptors = compute_pulse_descriptors(self.scenario.pulse, self.scenario.exposure)
        dosimetry = compute_dosimetry(descriptors, self.scenario.exposure)
        electro_state = compute_electrodynamics_state(
            descriptors,
            dosimetry,
            self.scenario.cell_state,
            self.params_flat,
        )

        scheduler = MultiscaleScheduler(descriptors, self.scenario.simulation.t_start_s, self.scenario.simulation.t_end_s)
        t_eval = scheduler.get_time_array(self.scenario.simulation.output_dt_s)
        t_span = (float(t_eval[0]), float(t_eval[-1]))
        y0 = [
            self.params_flat["Ca_baseline_uM"],
            self.params_flat["Ca_ER_uM"],
            self.params_flat["ROS_baseline"],
            self.params_flat["ATP_baseline"],
            0.0,
            0.0,
            0.0,
            0.0,
        ]

        def full_rhs(t: float, y: np.ndarray) -> list[float]:
            Ca_i, Ca_ER, ROS, ATP, sEV_cum, mlEV_cum, AB_cum, damage = [float(v) for v in y]
            Ca_i = max(Ca_i, 0.0)
            Ca_ER = max(Ca_ER, 0.0)
            ROS = max(ROS, 0.0)
            ATP = max(ATP, 0.0)
            damage = max(damage, 0.0)
            p = self.params_flat

            # Ca2+ flux (Layer 3 / Table A5-A6)
            decay_factor = np.exp(-t / max(p["tau_Ca_release_s"], 1e-9))
            J_Ca_in = (
                p["J_Ca_pore_factor"]
                * electro_state.membrane_permeability
                * max(Ca_ER - Ca_i, 0.0)
                * decay_factor
            )
            J_Ca_out = max(Ca_i - p["Ca_baseline_uM"], 0.0) / max(p["tau_Ca_uptake_s"], 1e-9)
            dCa_i = J_Ca_in - J_Ca_out
            dCa_ER = -J_Ca_in

            Ca_excess = max(Ca_i / max(p["Ca_baseline_uM"], 1e-9) - 1.0, 0.0)
            dROS = (
                p["ROS_production_factor"] * Ca_excess
                - (ROS - p["ROS_baseline"]) / max(p["tau_ROS_s"], 1e-9)
            )
            dATP = (
                -p["ATP_depletion_factor"] * max(ROS - p["ROS_baseline"], 0.0)
                / max(p["tau_ATP_s"], 1e-9)
                - (ATP - p["ATP_baseline"]) / max(p["tau_ATP_s"], 1e-9)
            )

            K_PS = p["K_PS_uM"]
            n_PS = p["n_PS"]
            PS = (Ca_i ** n_PS) / (K_PS ** n_PS + Ca_i ** n_PS) * p["PS_max"] if Ca_i > 0 else 0.0

            K_sEV = p["K_sEV_Ca_uM"]
            n_sEV = p["n_sEV_Ca"]
            K_mlEV = p["K_mlEV_PS"]
            n_mlEV = p["n_mlEV_PS"]
            K_AB = p["K_AB_damage"]
            n_AB = p["n_AB_damage"]

            sEV_base = p["baseline_sEV_rate"] * self.scenario.cell_state.baseline_EV_release_modifier
            Ca_hill = (Ca_i ** n_sEV) / (K_sEV ** n_sEV + Ca_i ** n_sEV)
            R_sEV = sEV_base * (1.0 + p["sEV_Ca_scale"] * Ca_hill)

            mlEV_base = p["baseline_mlEV_rate"] * self.scenario.cell_state.baseline_EV_release_modifier
            PS_hill = (PS ** n_mlEV) / (K_mlEV ** n_mlEV + PS ** n_mlEV)
            R_mlEV = mlEV_base * (1.0 + p["mlEV_PS_scale"] * PS_hill)

            damage_input = max(ROS - p["ROS_baseline"], 0.0) + max(Ca_i - p["Ca_max_uM"], 0.0)
            stress_mod = self.scenario.cell_state.stress_sensitivity_modifier
            if damage_input > 0:
                AB_denom = K_AB ** n_AB + damage_input ** n_AB
                AB_hill = (damage_input ** n_AB) / AB_denom
            else:
                AB_hill = 0.0
            R_AB = p["baseline_AB_rate"] * (1.0 + p["AB_damage_scale"] * stress_mod * AB_hill)

            dDamage = p["damage_rate"] * damage_input - p["repair_rate"] * damage / (1.0 + damage)
            return [float(dCa_i), float(dCa_ER), float(dROS), float(dATP), float(R_sEV), float(R_mlEV), float(R_AB), float(dDamage)]

        result = solve_ode_system(full_rhs, y0, t_span=t_span, t_eval=t_eval, method="RK45")
        y = result.y
        ca_i = np.clip(y[0], 0.0, None)
        ca_er = np.clip(y[1], 0.0, None)
        ros = np.clip(y[2], 0.0, None)
        atp = np.clip(y[3], 0.0, None)
        sev_cum = np.clip(y[4], 0.0, None)
        mlev_cum = np.clip(y[5], 0.0, None)
        ab_cum = np.clip(y[6], 0.0, None)
        damage = np.clip(y[7], 0.0, None)

        remodeling_params = RemodelingParams(
            K_PS_uM=self.params_flat["K_PS_uM"],
            n_PS=self.params_flat["n_PS"],
            PS_max=self.params_flat["PS_max"],
            K_calpain_uM=self.params_flat["K_calpain_uM"],
            n_calpain=self.params_flat["n_calpain"],
            K_annex_uM=self.params_flat["K_annex_uM"],
            n_annex=self.params_flat["n_annex"],
            tau_repair_s=self.params_flat["tau_repair_s"],
        )
        remodeling = [compute_remodeling_state(float(value), remodeling_params) for value in ca_i]
        ps = np.array([item["PS_exposure"] for item in remodeling], dtype=float)

        sev_rate = np.gradient(sev_cum, t_eval, edge_order=1)
        mlev_rate = np.gradient(mlev_cum, t_eval, edge_order=1)
        ab_rate = np.gradient(ab_cum, t_eval, edge_order=1)

        cargo_params = CargoPotencyParams(
            protein_enrichment_baseline=self.params_flat["protein_enrichment_baseline"],
            RNA_enrichment_baseline=self.params_flat["RNA_enrichment_baseline"],
            lipid_enrichment_baseline=self.params_flat["lipid_enrichment_baseline"],
            Ca_protein_coupling=self.params_flat["Ca_protein_coupling"],
            ROS_RNA_coupling=self.params_flat["ROS_RNA_coupling"],
            stress_lipid_coupling=self.params_flat["stress_lipid_coupling"],
            potency_weights=self.params_flat["potency_weights"],
        )
        cargo_state = compute_cargo_state(float(ca_i[-1]), float(ros[-1]), float(atp[-1]), cargo_params)

        injury_params = InjuryParams(
            K_apoptosis_damage=self.params_flat["K_apoptosis_damage"],
            K_necrosis_damage=self.params_flat["K_necrosis_damage"],
            n_apoptosis=self.params_flat["n_apoptosis"],
            n_necrosis=self.params_flat["n_necrosis"],
            debris_fraction_scale=self.params_flat["debris_fraction_scale"],
            contamination_threshold=self.params_flat["contamination_threshold"],
            damage_rate=self.params_flat["damage_rate"],
            repair_rate=self.params_flat["repair_rate"],
        )
        # Use peak damage for quality gate: captures worst-case acute cellular stress
        # rather than the post-recovery end-state, which is more biologically relevant
        # for assessing the EV-producing cell population during the harvest window.
        peak_damage = float(np.max(damage))
        debris = min(1.0, injury_params.debris_fraction_scale * peak_damage)
        quality = compute_quality_gate(peak_damage, debris, injury_params)

        manufacturing_params = ManufacturingParams(
            isolation_efficiency=self.params_flat["isolation_efficiency"],
            purity_factor=self.params_flat["purity_factor"],
            batch_consistency=self.params_flat["batch_consistency"],
            scalability_factor=self.params_flat["scalability_factor"],
        )
        manufacturing = compute_manufacturing_outputs(float(sev_cum[-1]), float(mlev_cum[-1]), float(ab_cum[-1]), float(quality["viability_fraction"]), manufacturing_params)

        state_timeseries = pd.DataFrame(
            {
                "t": t_eval,
                "Ca_i": ca_i,
                "Ca_ER": ca_er,
                "ROS": ros,
                "ATP": atp,
                "damage": damage,
            }
        )
        ev_timeseries = pd.DataFrame(
            {
                "t": t_eval,
                "sEV_rate": sev_rate,
                "mlEV_rate": mlev_rate,
                "AB_rate": ab_rate,
                "sEV_cumulative": sev_cum,
                "mlEV_cumulative": mlev_cum,
                "AB_cumulative": ab_cum,
            }
        )

        summary = {
            "scenario_name": self.scenario.scenario.name,
            "mode": self.scenario.scenario.mode,
            "dose_index": float(descriptors.dose_index),
            "peak_ca_i": float(np.max(ca_i)),
            "cumulative_small_EV": float(sev_cum[-1]),
            "cumulative_medium_large_EV": float(mlev_cum[-1]),
            "cumulative_apoptotic_body": float(ab_cum[-1]),
            "total_measured_particles": float(manufacturing["total_measured_particles"]),
            "purity_score": float(min(1.0, manufacturing["purity_score"] * float(quality["purity_score"]))),
            "viability_fraction": float(quality["viability_fraction"]),
            "potency_score": float(cargo_state["potency_score"]),
            "placeholder_fraction": 1.0,
            "warnings": self.warnings,
        }

        parameters_used = {
            "scenario": self.scenario.model_dump(mode="python"),
            "parameters_nested": self.params_nested,
            "parameters_flat": self.params_flat,
            "pulse_descriptors": asdict(descriptors),
            "dosimetry": asdict(dosimetry),
            "electrodynamics": asdict(electro_state),
            "terminal_remodeling": remodeling[-1],
            "terminal_cargo": cargo_state,
            "terminal_quality": quality,
            "manufacturing": manufacturing,
            "mean_ps_exposure": float(np.mean(ps)),
        }

        return SimulationResult(
            scenario_name=self.scenario.scenario.name,
            mode=self.scenario.scenario.mode,
            t_array=t_eval,
            state_timeseries=state_timeseries,
            ev_timeseries=ev_timeseries,
            summary=summary,
            parameters_used=parameters_used,
            warnings=self.warnings,
        )
