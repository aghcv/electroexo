from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from electro_exocytosis.config import SimulationScenario
from electro_exocytosis.io.readers import load_default_parameters
from electro_exocytosis.simulation import Simulation


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = REPO_ROOT / "results" / "ffrci_parameter_screen"


@dataclass(frozen=True)
class CandidateParameter:
    section: str
    name: str
    evidence_role: str


CANDIDATES = (
    CandidateParameter("electrodynamics", "cell_radius_m", "measure_or_literature_prior"),
    CandidateParameter("electrodynamics", "membrane_charging_tau_s", "literature_prior"),
    CandidateParameter("electrodynamics", "permeability_threshold_V", "literature_prior"),
    CandidateParameter("electrodynamics", "permeability_slope", "literature_prior"),
    CandidateParameter("electrodynamics", "delta_V_ER_fraction", "literature_prior"),
    CandidateParameter("electrodynamics", "delta_V_MVB_fraction", "literature_prior"),
    CandidateParameter("ion_transport", "tau_pore_reseal_s", "literature_prior"),
    CandidateParameter("ion_transport", "J_Ca_pore_factor", "literature_prior"),
    CandidateParameter("ion_transport", "J_ER_release_factor", "literature_prior"),
    CandidateParameter("ion_transport", "ER_activation_threshold_V", "literature_prior"),
    CandidateParameter("ion_transport", "PMCA_Vmax_uM_s", "literature_prior"),
    CandidateParameter("ion_transport", "tau_ROS_s", "literature_prior"),
    CandidateParameter("ion_transport", "ROS_pulse_factor_s", "literature_prior"),
    CandidateParameter("ion_transport", "ATP_depletion_factor", "literature_prior"),
    CandidateParameter("remodeling_repair", "K_scramblase_uM", "literature_prior"),
    CandidateParameter("remodeling_repair", "K_calpain_uM", "literature_prior"),
    CandidateParameter("remodeling_repair", "tau_repair_s", "literature_prior"),
    CandidateParameter("remodeling_repair", "microdomain_pore_gain", "literature_prior"),
    CandidateParameter("remodeling_repair", "resealing_annexin_weight", "literature_prior"),
    CandidateParameter("remodeling_repair", "shedding_rate_scale", "fit_after_upstream_calibration"),
    CandidateParameter("ev_release", "baseline_sEV_rate", "fit_with_ev_data"),
    CandidateParameter("ev_release", "baseline_mlEV_rate", "fit_with_ev_data"),
    CandidateParameter("ev_release", "baseline_AB_rate", "fit_only_with_subtype_or_injury_data"),
    CandidateParameter("ev_release", "k_MVB_maturation_s", "fit_with_ev_time_course"),
    CandidateParameter("ev_release", "k_MVB_docking_s", "fit_with_ev_time_course"),
    CandidateParameter("ev_release", "k_ILV_release_s", "fit_with_ev_time_course"),
    CandidateParameter("ev_release", "k_budding_s", "fit_with_subtype_data"),
    CandidateParameter("ev_release", "k_apoptotic_commitment_s", "fit_only_with_injury_data"),
    CandidateParameter("ev_release", "K_fusion_Ca_uM", "literature_prior"),
    CandidateParameter("ev_release", "K_PS_budding", "literature_prior"),
    CandidateParameter("ev_release", "repair_shedding_mlEV_weight", "fit_with_subtype_data"),
    CandidateParameter("injury_quality", "K_apoptosis_damage", "fit_only_with_injury_data"),
    CandidateParameter("injury_quality", "K_necrosis_damage", "fit_only_with_injury_data"),
    CandidateParameter("injury_quality", "damage_rate", "fit_only_with_viability_data"),
    CandidateParameter("injury_quality", "repair_rate", "fit_only_with_viability_data"),
)

OBSERVABLES = (
    "total_ev_output_proxy",
    "small_ev_fraction",
    "medium_large_ev_fraction",
    "apoptotic_body_fraction",
    "viability_fraction",
    "peak_ca_i_uM",
    "peak_ros",
    "atp_depletion_fraction",
    "peak_ps_exposure",
    "peak_repair_state",
)


def scenario_payload(amplitude_label_kv: int, pulse_count: int, harvest_time_h: float) -> dict[str, object]:
    return {
        "scenario": {
            "name": f"provisional_{pulse_count}p{amplitude_label_kv}kV_{harvest_time_h:g}h",
            "mode": "cell_based_electro_exocytosis",
        },
        "pulse": {
            # Provisional interpretation only. The experimental label may be generator voltage,
            # not electric field in kV/cm.
            "amplitude_kV_cm": float(amplitude_label_kv),
            "pulse_width_ns": 60.0,
            "pulse_number": pulse_count,
            "repetition_rate_Hz": 1.0,
            "waveform": "square",
        },
        "exposure": {
            "geometry": "cuvette",
            "medium_conductivity_S_m": 1.6,
            "temperature_C": 37.0,
            "cell_density_per_ml": 1_000_000.0,
            "dosimetry_model": "legacy",
        },
        "cell_state": {
            "cell_type": "CD4",
            "membrane_modifier": 1.0,
            "calcium_handling_modifier": 1.0,
            "baseline_EV_release_modifier": 1.0,
            "stress_sensitivity_modifier": 1.0,
        },
        "simulation": {
            "t_start_s": 0.0,
            "t_end_s": harvest_time_h * 3600.0,
            "output_dt_s": 60.0,
            "numerical_method": "solve_ivp",
        },
    }


def model_outputs(result: object) -> dict[str, float]:
    summary = result.summary
    small = float(summary["cumulative_small_EV"])
    medium_large = float(summary["cumulative_medium_large_EV"])
    apoptotic = float(summary["cumulative_apoptotic_body"])
    total = small + medium_large + apoptotic
    denominator = max(total, 1e-12)
    return {
        "total_ev_output_proxy": total,
        "small_ev_fraction": small / denominator,
        "medium_large_ev_fraction": medium_large / denominator,
        "apoptotic_body_fraction": apoptotic / denominator,
        "viability_fraction": float(summary["viability_fraction"]),
        "peak_ca_i_uM": float(summary["peak_ca_i"]),
        "peak_ros": float(summary["peak_ros"]),
        "atp_depletion_fraction": max(1.0 - float(summary["min_atp"]), 1e-12),
        "peak_ps_exposure": float(summary["peak_ps_exposure"]),
        "peak_repair_state": float(summary["peak_repair_state"]),
    }


def run_scenarios(params_override: dict[str, dict[str, float]] | None = None) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for amplitude_label_kv, pulse_count in ((20, 1), (40, 3), (40, 5)):
        for harvest_time_h in (0.5, 1.0, 3.0):
            scenario = SimulationScenario.model_validate(
                scenario_payload(amplitude_label_kv, pulse_count, harvest_time_h)
            )
            outputs = model_outputs(Simulation(scenario, params_override=params_override).run())
            rows.append(
                {
                    "condition": f"{pulse_count}p{amplitude_label_kv}kV",
                    "pulse_count": pulse_count,
                    "amplitude_label_kV": amplitude_label_kv,
                    "harvest_time_h": harvest_time_h,
                    **outputs,
                }
            )
    return pd.DataFrame(rows)


def screen_parameters(
    log_step: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    defaults = load_default_parameters()
    baseline = run_scenarios()
    records: list[dict[str, object]] = []
    sensitivity_columns: dict[str, np.ndarray] = {}
    row_labels = [
        f"{row.condition}|{row.harvest_time_h:g}h|{observable}"
        for row in baseline.itertuples()
        for observable in OBSERVABLES
    ]

    for candidate in CANDIDATES:
        default_value = float(defaults[candidate.section][candidate.name])
        low_value = default_value * np.exp(-log_step)
        high_value = default_value * np.exp(log_step)
        low = run_scenarios({candidate.section: {candidate.name: low_value}})
        high = run_scenarios({candidate.section: {candidate.name: high_value}})
        sensitivity_vector: list[float] = []
        for row_index in range(len(baseline)):
            for observable in OBSERVABLES:
                low_output = max(abs(float(low.iloc[row_index][observable])), 1e-12)
                high_output = max(abs(float(high.iloc[row_index][observable])), 1e-12)
                elasticity = (np.log(high_output) - np.log(low_output)) / (2.0 * log_step)
                sensitivity_vector.append(float(elasticity))
                records.append(
                    {
                        "section": candidate.section,
                        "parameter": candidate.name,
                        "evidence_role": candidate.evidence_role,
                        "default_value": default_value,
                        "condition": baseline.iloc[row_index]["condition"],
                        "harvest_time_h": baseline.iloc[row_index]["harvest_time_h"],
                        "observable": observable,
                        "local_log_elasticity": float(elasticity),
                    }
                )
        sensitivity_columns[f"{candidate.section}.{candidate.name}"] = np.asarray(sensitivity_vector)

    long_frame = pd.DataFrame.from_records(records)
    ranking_rows = []
    for candidate in CANDIDATES:
        subset = long_frame[
            (long_frame["section"] == candidate.section) & (long_frame["parameter"] == candidate.name)
        ]
        ev_subset = subset[
            subset["observable"].isin(
                [
                    "total_ev_output_proxy",
                    "small_ev_fraction",
                    "medium_large_ev_fraction",
                    "apoptotic_body_fraction",
                ]
            )
        ]
        state_subset = subset[~subset.index.isin(ev_subset.index)]
        values = subset["local_log_elasticity"].to_numpy(dtype=float)
        ranking_rows.append(
            {
                "section": candidate.section,
                "parameter": candidate.name,
                "evidence_role": candidate.evidence_role,
                "default_value": float(subset["default_value"].iloc[0]),
                "overall_rms_log_elasticity": float(np.sqrt(np.mean(np.square(values)))),
                "maximum_absolute_log_elasticity": float(np.max(np.abs(values))),
                "ev_observable_rms_log_elasticity": float(
                    np.sqrt(np.mean(np.square(ev_subset["local_log_elasticity"].to_numpy(dtype=float))))
                ),
                "state_observable_rms_log_elasticity": float(
                    np.sqrt(np.mean(np.square(state_subset["local_log_elasticity"].to_numpy(dtype=float))))
                ),
            }
        )
    ranking = pd.DataFrame(ranking_rows).sort_values(
        ["overall_rms_log_elasticity", "maximum_absolute_log_elasticity"], ascending=False
    )
    matrix = pd.DataFrame(sensitivity_columns, index=row_labels)

    # The Gram matrix identifies stiff and sloppy combinations of log parameters.
    gram = matrix.to_numpy(dtype=float).T @ matrix.to_numpy(dtype=float)
    eigenvalues = np.clip(np.linalg.eigvalsh(gram)[::-1], 0.0, None)
    spectrum = pd.DataFrame(
        {
            "eigenvalue_rank": np.arange(1, len(eigenvalues) + 1),
            "sensitivity_gram_eigenvalue": eigenvalues,
        }
    )
    return baseline, long_frame, ranking, matrix, spectrum


def write_readme(path: Path, ranking: pd.DataFrame, spectrum: pd.DataFrame, log_step: float) -> None:
    top = ranking.head(10)
    active_threshold = max(float(spectrum["sensitivity_gram_eigenvalue"].max()) * 1e-6, 1e-12)
    active_direction_count = int((spectrum["sensitivity_gram_eigenvalue"] > active_threshold).sum())
    lines = [
        "# Provisional FFRCI parameter importance screen",
        "",
        "This is a local sensitivity screen of the current model, not a calibration result. "
        "It interprets the 20 kV and 40 kV labels as electric field in kV/cm and assumes 60 ns square pulses, 1 Hz, cuvette geometry, 1.6 S/m conductivity, and 37 C. "
        "Those assumptions must be replaced after the experimental protocol is confirmed.",
        "",
        f"Each parameter was varied symmetrically by exp(+/-{log_step:g}) around its default. "
        "The reported values are local log elasticities, d log(output) / d log(parameter), across the three labeled exposure regimes and three harvest times.",
        "",
        "## Highest local importance scores",
        "",
        "| Parameter | Evidence role | Overall RMS elasticity | Maximum absolute elasticity |",
        "|---|---|---:|---:|",
    ]
    for row in top.itertuples():
        lines.append(
            f"| `{row.section}.{row.parameter}` | {row.evidence_role} | "
            f"{row.overall_rms_log_elasticity:.3f} | {row.maximum_absolute_log_elasticity:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Identifiability warning",
            "",
            f"The local sensitivity Gram matrix has {active_direction_count} numerically active directions among {len(CANDIDATES)} screened parameters at a relative eigenvalue threshold of 1e-6. "
            "This does not prove practical identifiability because the RNA-seq measurements are indirect state proxies and the Exoid size bins do not identify EV biogenesis classes. "
            "Use this screen to reduce the parameter set before fitting, then apply profile likelihood or posterior diagnostics to the reduced model.",
            "",
            "The current model predicts nearly identical cumulative EV output for the three exposure regimes and monotonically increasing output over time. "
            "The measured condition-level concentration peaks at 1 h for 3p40kV and declines by 3 h for 5p40kV. "
            "This mismatch indicates a missing observation layer such as vesicle clearance, uptake, degradation, isolation recovery, or time-varying viable-cell number. "
            "A parameter fit should not proceed until the measured quantity and collection protocol are clarified and this observation layer is specified.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Screen local parameter importance for the FFRCI experiment design.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--log-step", type=float, default=np.log(1.25))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    baseline, long_frame, ranking, matrix, spectrum = screen_parameters(args.log_step)
    baseline.to_csv(args.out / "provisional_model_predictions.csv", index=False)
    long_frame.to_csv(args.out / "parameter_local_sensitivity_long.csv", index=False)
    ranking.to_csv(args.out / "parameter_importance_ranking.csv", index=False)
    matrix.to_csv(args.out / "parameter_sensitivity_matrix.csv", index_label="scenario_observable")
    spectrum.to_csv(args.out / "sensitivity_gram_eigenvalue_spectrum.csv", index=False)
    write_readme(args.out / "README.md", ranking, spectrum, args.log_step)
    print(f"Wrote provisional parameter screen to {args.out}")


if __name__ == "__main__":
    main()
