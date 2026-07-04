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
from electro_exocytosis.models.ev_release import (
    compute_ev_release_derivatives,
    compute_ev_release_fluxes,
    coerce_ev_release_params,
    get_ev_initial_conditions,
    get_ev_state_names,
)
from electro_exocytosis.models.injury_quality import InjuryParams, compute_quality_gate
from electro_exocytosis.models.ion_transport import (
    build_ion_transport_rhs,
    coerce_ion_transport_params,
    compute_ion_transport_fluxes,
    get_ion_initial_conditions,
    get_ion_state_names,
    ion_state_to_dict,
)
from electro_exocytosis.models.manufacturing_qc import ManufacturingParams, compute_manufacturing_outputs
from electro_exocytosis.models.pulse import compute_pulse_descriptors
from electro_exocytosis.models.remodeling_repair import coerce_remodeling_params, compute_remodeling_state
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
            "Scientific modules use reduced, uncalibrated equations and literature-derived default parameters.",
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
        ion_params = coerce_ion_transport_params(self.params_flat)
        ev_params = coerce_ev_release_params(self.params_flat)
        remodeling_params = coerce_remodeling_params(self.params_flat)
        ion_state_names = get_ion_state_names()
        ev_state_names = get_ev_state_names()
        effective_ion_perturbation_s = min(
            descriptors.train_duration_s,
            ion_params.tau_pore_reseal_s * np.sqrt(max(descriptors.pulse_number, 1)),
        )
        ion_rhs = build_ion_transport_rhs(
            ion_params,
            electro_state,
            t_pulse_end=effective_ion_perturbation_s,
        )

        scheduler = MultiscaleScheduler(descriptors, self.scenario.simulation.t_start_s, self.scenario.simulation.t_end_s)
        t_eval = scheduler.get_time_array(self.scenario.simulation.output_dt_s)
        t_span = (float(t_eval[0]), float(t_eval[-1]))
        y0 = [
            *get_ion_initial_conditions(ion_params),
            *get_ev_initial_conditions(ev_params),
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        ion_state_count = len(ion_state_names)
        ev_state_count = len(ev_state_names)

        def full_rhs(t: float, y: np.ndarray) -> list[float]:
            ion_y = y[:ion_state_count]
            ev_y = y[ion_state_count : ion_state_count + ev_state_count]
            ion_state = ion_state_to_dict(ion_y, ion_params)
            ion_derivatives = ion_rhs(t, ion_y)
            sEV_cum, mlEV_cum, AB_cum, damage = [
                float(v) for v in y[ion_state_count + ev_state_count :]
            ]
            Ca_i = ion_state["Ca_i"]
            ROS = ion_state["ROS"]
            ATP = ion_state["ATP"]
            osmotic_stress = ion_state["osmotic_stress"]
            mitochondrial_potential = ion_state["mitochondrial_potential"]
            damage = max(damage, 0.0)
            p = self.params_flat
            fluxes = compute_ion_transport_fluxes(
                ion_params,
                electro_state,
                t,
                ion_y,
                effective_ion_perturbation_s,
            )
            remodeling_state = compute_remodeling_state(
                Ca_i,
                remodeling_params,
                osmotic_stress=float(osmotic_stress),
                mitochondrial_potential=float(mitochondrial_potential),
                pore_activation=float(fluxes["pore_activation"]),
            )

            damage_input = (
                max(ROS - ion_params.ROS_baseline, 0.0)
                + max(Ca_i - ion_params.Ca_max_uM, 0.0)
                + ion_params.osmotic_damage_scale * osmotic_stress
                + ion_params.mitochondrial_damage_scale * max(1.0 - mitochondrial_potential, 0.0)
                + ion_params.ATP_damage_scale * max(ion_params.ATP_damage_threshold - ATP, 0.0)
            )
            ev_fluxes = compute_ev_release_fluxes(
                ev_params,
                ev_y,
                Ca_i=Ca_i,
                Ca_submembrane=float(remodeling_state["Ca_submembrane"]),
                ROS=ROS,
                ATP=ATP,
                damage_state=damage,
                delta_V_MVB=float(electro_state.delta_V_MVB),
                pore_activation=float(fluxes["pore_activation"]),
                PS_exposure=float(remodeling_state["PS_exposure"]),
                calpain_activity=float(remodeling_state["calpain_activity"]),
                annexin_activity=float(remodeling_state["annexin_activity"]),
                actomyosin_tension=float(remodeling_state["actomyosin_tension"]),
                actin_disruption=float(remodeling_state["actin_disruption"]),
                repair_state=float(remodeling_state["repair_state"]),
                repair_shedding_rate=float(remodeling_state["repair_shedding_rate"]),
            )
            ev_derivatives = compute_ev_release_derivatives(ev_params, ev_y, ev_fluxes)

            dDamage = p["damage_rate"] * damage_input - p["repair_rate"] * damage / (1.0 + damage)
            return [
                *ion_derivatives,
                *ev_derivatives,
                float(ev_fluxes["sEV_rate"]),
                float(ev_fluxes["mlEV_rate"]),
                float(ev_fluxes["AB_rate"]),
                float(dDamage),
            ]

        result = solve_ode_system(full_rhs, y0, t_span=t_span, t_eval=t_eval, method="RK45")
        y = result.y
        ion_series = {
            name: np.clip(y[index], 0.0, None)
            for index, name in enumerate(ion_state_names)
        }
        ev_series = {
            name: np.clip(y[ion_state_count + index], 0.0, None)
            for index, name in enumerate(ev_state_names)
        }
        ion_series["Ca_ER"] = np.clip(ion_series["Ca_ER"], 0.0, ion_params.Ca_ER_uM)
        ion_series["mitochondrial_potential"] = np.clip(ion_series["mitochondrial_potential"], 0.0, 1.0)
        ion_series["ATP"] = np.clip(ion_series["ATP"], 0.0, ion_params.ATP_baseline)
        ca_i = ion_series["Ca_i"]
        ca_er = ion_series["Ca_ER"]
        ca_mito = ion_series["Ca_mito"]
        mito_potential = ion_series["mitochondrial_potential"]
        ros = ion_series["ROS"]
        atp = ion_series["ATP"]
        osmotic_stress = ion_series["osmotic_stress"]
        sev_cum = np.clip(y[ion_state_count + ev_state_count], 0.0, None)
        mlev_cum = np.clip(y[ion_state_count + ev_state_count + 1], 0.0, None)
        ab_cum = np.clip(y[ion_state_count + ev_state_count + 2], 0.0, None)
        damage = np.clip(y[ion_state_count + ev_state_count + 3], 0.0, None)

        state_timeseries = pd.DataFrame({"t": t_eval, **ion_series, "damage": damage})
        flux_records = [
            compute_ion_transport_fluxes(
                ion_params,
                electro_state,
                float(t),
                [ion_series[name][index] for name in ion_state_names],
                effective_ion_perturbation_s,
            )
            for index, t in enumerate(t_eval)
        ]
        flux_timeseries = pd.DataFrame(flux_records)
        remodeling = [
            compute_remodeling_state(
                float(ca_i[index]),
                remodeling_params,
                osmotic_stress=float(osmotic_stress[index]),
                mitochondrial_potential=float(mito_potential[index]),
                pore_activation=float(flux_timeseries["pore_activation"].iloc[index]),
            )
            for index in range(len(t_eval))
        ]
        remodeling_timeseries = pd.DataFrame(remodeling)
        ev_release = [
            compute_ev_release_fluxes(
                ev_params,
                [ev_series[name][index] for name in ev_state_names],
                Ca_i=float(ca_i[index]),
                Ca_submembrane=float(remodeling_timeseries["Ca_submembrane"].iloc[index]),
                ROS=float(ros[index]),
                ATP=float(atp[index]),
                damage_state=float(damage[index]),
                delta_V_MVB=float(electro_state.delta_V_MVB),
                pore_activation=float(flux_timeseries["pore_activation"].iloc[index]),
                PS_exposure=float(remodeling_timeseries["PS_exposure"].iloc[index]),
                calpain_activity=float(remodeling_timeseries["calpain_activity"].iloc[index]),
                annexin_activity=float(remodeling_timeseries["annexin_activity"].iloc[index]),
                actomyosin_tension=float(remodeling_timeseries["actomyosin_tension"].iloc[index]),
                actin_disruption=float(remodeling_timeseries["actin_disruption"].iloc[index]),
                repair_state=float(remodeling_timeseries["repair_state"].iloc[index]),
                repair_shedding_rate=float(remodeling_timeseries["repair_shedding_rate"].iloc[index]),
            )
            for index in range(len(t_eval))
        ]
        ev_release_timeseries = pd.DataFrame(ev_release)
        ps = remodeling_timeseries["PS_exposure"].to_numpy(dtype=float)
        ev_state_frame = pd.DataFrame(ev_series)
        state_timeseries = pd.concat(
            [
                state_timeseries,
                flux_timeseries,
                remodeling_timeseries,
                ev_state_frame,
                ev_release_timeseries,
            ],
            axis=1,
        )

        cargo_params = CargoPotencyParams(
            protein_enrichment_baseline=self.params_flat["protein_enrichment_baseline"],
            RNA_enrichment_baseline=self.params_flat["RNA_enrichment_baseline"],
            lipid_enrichment_baseline=self.params_flat["lipid_enrichment_baseline"],
            Ca_protein_coupling=self.params_flat["Ca_protein_coupling"],
            ROS_RNA_coupling=self.params_flat["ROS_RNA_coupling"],
            stress_lipid_coupling=self.params_flat["stress_lipid_coupling"],
            subtype_sEV_weight=self.params_flat["subtype_sEV_weight"],
            subtype_mlEV_weight=self.params_flat["subtype_mlEV_weight"],
            subtype_AB_weight=self.params_flat["subtype_AB_weight"],
            ESCRT_protein_weight=self.params_flat["ESCRT_protein_weight"],
            rbp_RNA_sorting_weight=self.params_flat["rbp_RNA_sorting_weight"],
            ceramide_lipid_weight=self.params_flat["ceramide_lipid_weight"],
            antigen_stress_weight=self.params_flat["antigen_stress_weight"],
            antigen_sorting_weight=self.params_flat["antigen_sorting_weight"],
            direct_loading_efficiency=self.params_flat["direct_loading_efficiency"],
            direct_loading_leak_fraction=self.params_flat["direct_loading_leak_fraction"],
            recipient_dose_half_max=self.params_flat["recipient_dose_half_max"],
            potency_saturation=self.params_flat["potency_saturation"],
            potency_weights=self.params_flat["potency_weights"],
        )
        cargo_state = compute_cargo_state(
            float(ca_i[-1]),
            float(ros[-1]),
            float(atp[-1]),
            cargo_params,
            cumulative_sEV=float(sev_cum[-1]),
            cumulative_mlEV=float(mlev_cum[-1]),
            cumulative_AB=float(ab_cum[-1]),
            escrt_signal=float(ev_release_timeseries["escrt_dependent_signal"].iloc[-1]),
            ceramide_signal=float(ev_release_timeseries["ceramide_signal"].iloc[-1]),
            secretory_bias=float(ev_release_timeseries["secretory_bias"].iloc[-1]),
            direct_loading_drive=float(self.scenario.scenario.mode == "direct_EV_engineering"),
        )

        injury_params = InjuryParams(
            K_apoptosis_damage=self.params_flat["K_apoptosis_damage"],
            K_necrosis_damage=self.params_flat["K_necrosis_damage"],
            K_stress_damage=self.params_flat["K_stress_damage"],
            n_apoptosis=self.params_flat["n_apoptosis"],
            n_necrosis=self.params_flat["n_necrosis"],
            debris_fraction_scale=self.params_flat["debris_fraction_scale"],
            aggregate_fraction_scale=self.params_flat["aggregate_fraction_scale"],
            apoptotic_body_contamination_weight=self.params_flat["apoptotic_body_contamination_weight"],
            necrotic_debris_weight=self.params_flat["necrotic_debris_weight"],
            marker_panel_weight=self.params_flat["marker_panel_weight"],
            contamination_threshold=self.params_flat["contamination_threshold"],
            viability_threshold=self.params_flat["viability_threshold"],
            damage_rate=self.params_flat["damage_rate"],
            repair_rate=self.params_flat["repair_rate"],
        )
        # Use peak damage for quality gate: captures worst-case acute cellular stress
        # rather than the post-recovery end-state, which is more biologically relevant
        # for assessing the EV-producing cell population during the harvest window.
        peak_damage = float(np.max(damage))
        debris = min(1.0, injury_params.debris_fraction_scale * peak_damage)
        quality = compute_quality_gate(
            peak_damage,
            debris,
            injury_params,
            cumulative_sEV=float(sev_cum[-1]),
            cumulative_mlEV=float(mlev_cum[-1]),
            cumulative_AB=float(ab_cum[-1]),
        )

        manufacturing_params = ManufacturingParams(
            cell_count=self.params_flat["cell_count"],
            harvest_time_h=self.params_flat["harvest_time_h"],
            isolation_efficiency=self.params_flat["isolation_efficiency"],
            isolation_method_factor=self.params_flat["isolation_method_factor"],
            purity_factor=self.params_flat["purity_factor"],
            protein_contamination_factor=self.params_flat["protein_contamination_factor"],
            batch_consistency=self.params_flat["batch_consistency"],
            batch_variability_fraction=self.params_flat["batch_variability_fraction"],
            scalability_factor=self.params_flat["scalability_factor"],
            potency_weight=self.params_flat["potency_weight"],
            yield_weight=self.params_flat["yield_weight"],
            purity_weight=self.params_flat["purity_weight"],
            viability_weight=self.params_flat["viability_weight"],
        )
        manufacturing = compute_manufacturing_outputs(
            float(sev_cum[-1]),
            float(mlev_cum[-1]),
            float(ab_cum[-1]),
            float(quality["viability_fraction"]),
            manufacturing_params,
            potency_score=float(cargo_state["potency_score"]),
            purity_score=float(quality["purity_score"]),
        )
        ev_timeseries = pd.DataFrame(
            {
                "t": t_eval,
                "sEV_rate": ev_release_timeseries["sEV_rate"].to_numpy(dtype=float),
                "mlEV_rate": ev_release_timeseries["mlEV_rate"].to_numpy(dtype=float),
                "AB_rate": ev_release_timeseries["AB_rate"].to_numpy(dtype=float),
                "sEV_cumulative": sev_cum,
                "mlEV_cumulative": mlev_cum,
                "AB_cumulative": ab_cum,
                **ev_series,
                "rab_conversion_signal": ev_release_timeseries["rab_conversion_signal"].to_numpy(dtype=float),
                "escrt_dependent_signal": ev_release_timeseries["escrt_dependent_signal"].to_numpy(dtype=float),
                "ceramide_signal": ev_release_timeseries["ceramide_signal"].to_numpy(dtype=float),
                "acidification_signal": ev_release_timeseries["acidification_signal"].to_numpy(dtype=float),
                "lysosomal_routing": ev_release_timeseries["lysosomal_routing"].to_numpy(dtype=float),
                "secretory_bias": ev_release_timeseries["secretory_bias"].to_numpy(dtype=float),
                "rab_docking_signal": ev_release_timeseries["rab_docking_signal"].to_numpy(dtype=float),
                "fusion_signal": ev_release_timeseries["fusion_signal"].to_numpy(dtype=float),
                "budding_signal": ev_release_timeseries["budding_signal"].to_numpy(dtype=float),
                "scission_signal": ev_release_timeseries["scission_signal"].to_numpy(dtype=float),
                "apoptotic_blebbing_signal": ev_release_timeseries["apoptotic_blebbing_signal"].to_numpy(dtype=float),
            }
        )

        summary = {
            "scenario_name": self.scenario.scenario.name,
            "mode": self.scenario.scenario.mode,
            "dose_index": float(descriptors.dose_index),
            "peak_ca_i": float(np.max(ca_i)),
            "peak_ca_mito": float(np.max(ca_mito)),
            "peak_ros": float(np.max(ros)),
            "min_atp": float(np.min(atp)),
            "min_mitochondrial_potential": float(np.min(mito_potential)),
            "peak_osmotic_stress": float(np.max(osmotic_stress)),
            "peak_ps_exposure": float(np.max(ps)),
            "peak_calpain_activity": float(remodeling_timeseries["calpain_activity"].max()),
            "peak_annexin_activity": float(remodeling_timeseries["annexin_activity"].max()),
            "peak_actin_disruption": float(remodeling_timeseries["actin_disruption"].max()),
            "peak_repair_state": float(remodeling_timeseries["repair_state"].max()),
            "peak_secretory_bias": float(ev_release_timeseries["secretory_bias"].max()),
            "peak_lysosomal_routing": float(ev_release_timeseries["lysosomal_routing"].max()),
            "peak_docked_MVB_pool": float(ev_state_frame["docked_MVB_pool"].max()),
            "peak_budding_pool": float(ev_state_frame["budding_pool"].max()),
            "peak_apoptotic_commitment": float(ev_state_frame["apoptotic_commitment"].max()),
            "cumulative_repair_shedding": float(np.trapezoid(remodeling_timeseries["repair_shedding_rate"], t_eval)),
            "cumulative_small_EV": float(sev_cum[-1]),
            "cumulative_medium_large_EV": float(mlev_cum[-1]),
            "cumulative_apoptotic_body": float(ab_cum[-1]),
            "total_measured_particles": float(manufacturing["total_measured_particles"]),
            "purity_score": float(min(1.0, manufacturing["purity_score"] * float(quality["purity_score"]))),
            "viability_fraction": float(quality["viability_fraction"]),
            "stressed_viable_fraction": float(quality["stressed_viable_fraction"]),
            "apoptosis_fraction": float(quality["apoptosis_fraction"]),
            "necrosis_fraction": float(quality["necrosis_fraction"]),
            "bona_fide_EV_fraction": float(quality["bona_fide_EV_fraction"]),
            "protein_enrichment": float(cargo_state["protein_enrichment"]),
            "RNA_enrichment": float(cargo_state["RNA_enrichment"]),
            "lipid_enrichment": float(cargo_state["lipid_enrichment"]),
            "antigen_enrichment": float(cargo_state["antigen_enrichment"]),
            "direct_loaded_cargo": float(cargo_state["direct_loaded_cargo"]),
            "potency_score": float(cargo_state["potency_score"]),
            "cell_normalized_yield": float(manufacturing["cell_normalized_yield"]),
            "batch_adjusted_yield": float(manufacturing["batch_adjusted_yield"]),
            "optimization_objective": float(manufacturing["optimization_objective"]),
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
            "ion_transport": asdict(ion_params),
            "remodeling_repair": asdict(remodeling_params),
            "ev_release": asdict(ev_params),
            "effective_ion_perturbation_s": float(effective_ion_perturbation_s),
            "terminal_remodeling": remodeling[-1],
            "terminal_ev_release": ev_release[-1],
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
