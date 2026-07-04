from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
from matplotlib.figure import Figure


@dataclass(frozen=True)
class AbbreviationEntry:
    short: str
    long: str
    plot_label: str | None = None
    csv_name: str | None = None


class AbbreviationRegistry:
    """Central abbreviation and standard-language export policy."""

    def __init__(
        self,
        entries: dict[str, AbbreviationEntry],
        *,
        exact_column_aliases: dict[str, str] | None = None,
    ) -> None:
        self._entries = dict(entries)
        self._exact_column_aliases = dict(exact_column_aliases or {})
        self._token_aliases = (
            ("ca_submembrane", "submembrane_calcium"),
            ("ca_mito", "mitochondrial_calcium"),
            ("ca_er", "endoplasmic_reticulum_calcium"),
            ("ca_i", "cytosolic_calcium"),
            ("mlev", "medium_large_extracellular_vesicle"),
            ("sev", "small_extracellular_vesicle"),
            ("mvb", "multivesicular_body"),
            ("ilv", "intraluminal_vesicle"),
            ("ros", "reactive_oxygen_species"),
            ("atp", "adenosine_triphosphate"),
            ("ps", "phosphatidylserine"),
            ("ab", "apoptotic_body"),
        )

    def entry(self, key: str) -> AbbreviationEntry | None:
        return self._entries.get(key)

    def plot_label(self, key: str) -> str:
        entry = self.entry(key)
        if entry and entry.plot_label:
            return entry.plot_label
        if entry:
            return entry.short
        return self._humanize_label(self.csv_name(key))

    def csv_name(self, key: str) -> str:
        entry = self.entry(key)
        if entry and entry.csv_name:
            return entry.csv_name
        if key in self._exact_column_aliases:
            return self._exact_column_aliases[key]

        normalized = key.lower()
        normalized = normalized.replace("2+", "2plus")
        normalized = normalized.replace("/", "_")
        normalized = normalized.replace("-", "_")
        normalized = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", normalized).lower()
        for token, replacement in self._token_aliases:
            normalized = re.sub(
                rf"(?<![a-z0-9]){re.escape(token)}(?![a-z0-9])",
                replacement,
                normalized,
            )
        normalized = re.sub(r"[^a-z0-9]+", "_", normalized)
        normalized = re.sub(r"_+", "_", normalized).strip("_")
        return normalized or key

    def rename_columns(self, frame: pd.DataFrame) -> pd.DataFrame:
        return frame.rename(columns={column: self.csv_name(column) for column in frame.columns})

    def note(self, keys: Iterable[str], *, latex: bool = False) -> str:
        items: list[str] = []
        seen: set[str] = set()
        for key in keys:
            entry = self.entry(key)
            if entry is None or entry.short in seen:
                continue
            seen.add(entry.short)
            if latex:
                items.append(f"{entry.short}, {entry.long}")
            else:
                items.append(f"{entry.short}, {entry.long}")
        if not items:
            return ""
        return "Abbreviations: " + "; ".join(items) + "."

    def add_figure_note(self, fig: Figure, keys: Iterable[str]) -> None:
        note = self.note(keys)
        if not note:
            return
        fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, 0.18))
        fig.text(
            0.01,
            0.015,
            note,
            ha="left",
            va="bottom",
            fontsize=7,
            color="#222222",
            wrap=True,
        )

    def write_bundle(self, outdir: Path, *, keys: Iterable[str] | None = None) -> None:
        outdir.mkdir(parents=True, exist_ok=True)
        entries = self._selected_entries(keys)
        payload = [asdict(entry) for entry in entries]
        (outdir / "abbreviations.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        (outdir / "abbreviations.md").write_text(self._markdown(entries), encoding="utf-8")
        (outdir / "abbreviations.tex").write_text(self._latex(entries), encoding="utf-8")

    def _selected_entries(self, keys: Iterable[str] | None) -> list[AbbreviationEntry]:
        if keys is None:
            return sorted(self._entries.values(), key=lambda entry: entry.short.lower())
        selected: list[AbbreviationEntry] = []
        seen: set[str] = set()
        for key in keys:
            entry = self.entry(key)
            if entry is None or entry.short in seen:
                continue
            seen.add(entry.short)
            selected.append(entry)
        return selected

    @staticmethod
    def _humanize_label(key: str) -> str:
        key = key.replace("_", " ").strip()
        return re.sub(r"\s+", " ", key).title()

    @staticmethod
    def _markdown(entries: list[AbbreviationEntry]) -> str:
        lines = [
            "# Abbreviations",
            "",
            "| Short form | Expanded form |",
            "| --- | --- |",
        ]
        lines.extend(f"| {entry.short} | {entry.long} |" for entry in entries)
        return "\n".join(lines) + "\n"

    @staticmethod
    def _latex(entries: list[AbbreviationEntry]) -> str:
        lines = [
            "% Generated by electro_exocytosis.abbreviations",
            "\\begin{tabular}{ll}",
            "\\toprule",
            "\\textbf{Short form} & \\textbf{Expanded form} \\\\",
            "\\midrule",
        ]
        lines.extend(f"{entry.short} & {entry.long} \\\\" for entry in entries)
        lines.extend(["\\bottomrule", "\\end{tabular}", ""])
        return "\n".join(lines)


STANDARD_ABBREVIATIONS = AbbreviationRegistry(
    entries={
        "AB": AbbreviationEntry("AB", "apoptotic body"),
        "AB_rate": AbbreviationEntry(
            "AB",
            "apoptotic body",
            plot_label="Apoptotic-body release rate",
            csv_name="apoptotic_body_rate",
        ),
        "ATP": AbbreviationEntry(
            "ATP",
            "adenosine triphosphate",
            plot_label="ATP",
            csv_name="adenosine_triphosphate_relative_state",
        ),
        "Ca_ER": AbbreviationEntry(
            "ER",
            "endoplasmic reticulum",
            plot_label="Endoplasmic-reticulum calcium",
            csv_name="endoplasmic_reticulum_calcium_uM",
        ),
        "Ca_i": AbbreviationEntry(
            "Ca2+",
            "cytosolic calcium",
            plot_label="Cytosolic calcium",
            csv_name="cytosolic_calcium_uM",
        ),
        "Ca_mito": AbbreviationEntry(
            "Mito",
            "mitochondrial compartment",
            plot_label="Mitochondrial calcium",
            csv_name="mitochondrial_calcium_uM",
        ),
        "Ca_submembrane": AbbreviationEntry(
            "Ca2+",
            "submembrane calcium",
            plot_label="Submembrane calcium",
            csv_name="submembrane_calcium_uM",
        ),
        "ER": AbbreviationEntry("ER", "endoplasmic reticulum"),
        "ESCRT": AbbreviationEntry("ESCRT", "endosomal sorting complex required for transport"),
        "EV": AbbreviationEntry("EV", "extracellular vesicle"),
        "ILV": AbbreviationEntry("ILV", "intraluminal vesicle"),
        "ILV_load": AbbreviationEntry(
            "ILV",
            "intraluminal vesicle",
            plot_label="Intraluminal-vesicle load",
            csv_name="intraluminal_vesicle_load",
        ),
        "MVB": AbbreviationEntry("MVB", "multivesicular body"),
        "MVB_pool": AbbreviationEntry(
            "MVB",
            "multivesicular body",
            plot_label="Multivesicular-body pool",
            csv_name="multivesicular_body_pool",
        ),
        "PM": AbbreviationEntry("PM", "plasma membrane"),
        "PS": AbbreviationEntry("PS", "phosphatidylserine"),
        "PS_exposure": AbbreviationEntry(
            "PS",
            "phosphatidylserine",
            plot_label="Phosphatidylserine exposure",
            csv_name="phosphatidylserine_exposure_fraction",
        ),
        "ROS": AbbreviationEntry(
            "ROS",
            "reactive oxygen species",
            plot_label="Reactive oxygen species",
            csv_name="reactive_oxygen_species_relative_state",
        ),
        "m/lEV": AbbreviationEntry("m/lEV", "medium/large extracellular vesicle"),
        "mlEV_rate": AbbreviationEntry(
            "m/lEV",
            "medium/large extracellular vesicle",
            plot_label="Medium/large-EV release rate",
            csv_name="medium_large_extracellular_vesicle_rate",
        ),
        "nsPEF": AbbreviationEntry("nsPEF", "nanosecond pulsed electric field"),
        "sEV": AbbreviationEntry("sEV", "small extracellular vesicle"),
        "sEV_rate": AbbreviationEntry(
            "sEV",
            "small extracellular vesicle",
            plot_label="Small-EV release rate",
            csv_name="small_extracellular_vesicle_rate",
        ),
    },
    exact_column_aliases={
        "t": "time_s",
        "AB_cumulative": "apoptotic_body_cumulative_output",
        "AB_rate": "apoptotic_body_rate",
        "ATP": "adenosine_triphosphate_relative_state",
        "Ca_ER": "endoplasmic_reticulum_calcium_uM",
        "Ca_i": "cytosolic_calcium_uM",
        "Ca_mito": "mitochondrial_calcium_uM",
        "Ca_submembrane": "submembrane_calcium_uM",
        "Cl_i": "intracellular_chloride_mM",
        "ILV_load": "intraluminal_vesicle_load",
        "J_Ca_pore": "electropore_calcium_flux_uM_s",
        "J_ER_release": "endoplasmic_reticulum_calcium_release_flux_uM_s",
        "K_i": "intracellular_potassium_mM",
        "MVB_pool": "multivesicular_body_pool",
        "Na_i": "intracellular_sodium_mM",
        "PS_exposure": "phosphatidylserine_exposure_fraction",
        "ROS": "reactive_oxygen_species_relative_state",
        "abbreviations": "abbreviations",
        "acidification_signal": "acidification_signal",
        "actin_disruption": "actin_disruption_state",
        "actomyosin_tension": "actomyosin_tension_state",
        "annexin_activity": "annexin_activity_state",
        "apoptotic_blebbing_signal": "apoptotic_blebbing_signal",
        "apoptotic_commitment": "apoptotic_commitment_state",
        "budding_pool": "budding_precursor_pool",
        "budding_signal": "budding_signal",
        "calpain_activity": "calpain_activity_state",
        "ceramide_signal": "ceramide_signal",
        "damage": "damage_index",
        "docked_MVB_pool": "docked_multivesicular_body_pool",
        "escrt_dependent_signal": "escrt_dependent_signal",
        "fusion_signal": "fusion_signal",
        "lysosomal_repair_activity": "lysosomal_repair_activity_state",
        "lysosomal_routing": "lysosomal_routing_bias",
        "membrane_permeability": "membrane_permeability_proxy",
        "mitochondrial_potential": "mitochondrial_membrane_potential",
        "mlEV_cumulative": "medium_large_extracellular_vesicle_cumulative_output",
        "mlEV_rate": "medium_large_extracellular_vesicle_rate",
        "osmotic_stress": "osmotic_stress_index",
        "pore_activation": "pore_activation_state",
        "rab_conversion_signal": "rab_conversion_signal",
        "rab_docking_signal": "rab_docking_signal",
        "repair_shedding_rate": "repair_associated_shedding_rate",
        "repair_state": "resealing_state",
        "sEV_cumulative": "small_extracellular_vesicle_cumulative_output",
        "sEV_rate": "small_extracellular_vesicle_rate",
        "scission_signal": "scission_signal",
        "scramblase_activity": "scramblase_activity_state",
        "secretory_bias": "secretory_routing_bias",
        "viability_fraction": "viability_fraction",
        "delta_Vm_V": "plasma_membrane_voltage_V",
        "delta_V_ER_V": "endoplasmic_reticulum_voltage_V",
        "delta_V_mito_V": "mitochondrial_voltage_V",
        "delta_V_MVB_V": "multivesicular_body_voltage_V",
        "E_peak_kV_cm": "peak_electric_field_kV_cm",
        "adiabatic_temp_rise_C": "adiabatic_temperature_rise_C",
        "absorbed_energy_density_mJ_mm3": "absorbed_energy_density_mJ_mm3",
        "amplitude_kV_cm": "amplitude_kV_cm",
        "cell_radius_um": "cell_radius_um",
        "cumulative_repair_shedding": "cumulative_repair_associated_shedding",
        "peak_J_Ca_pore_uM_s": "peak_electropore_calcium_flux_uM_s",
        "peak_J_ER_release_uM_s": "peak_endoplasmic_reticulum_calcium_release_flux_uM_s",
        "dose_index": "dose_index",
        "end_temp_rise_C": "end_temperature_rise_C",
        "end_temperature_C": "end_temperature_C",
        "field_uniformity_factor": "field_uniformity_factor",
        "geometry": "geometry",
        "geometry_heat_density_mJ_mm3": "geometry_adjusted_heat_density_mJ_mm3",
        "mean_E2_factor": "mean_electric_field_squared_factor",
        "medium_conductivity_S_m": "medium_conductivity_S_m",
        "membrane_charging_factor": "membrane_charging_factor",
        "membrane_charging_tau_ns": "membrane_charging_time_constant_ns",
        "min_mitochondrial_potential": "minimum_mitochondrial_membrane_potential",
        "model_note": "model_note",
        "peak_actin_disruption": "peak_actin_disruption",
        "peak_annexin_activity": "peak_annexin_activity",
        "peak_apoptotic_commitment": "peak_apoptotic_commitment",
        "peak_ca_i": "peak_cytosolic_calcium_uM",
        "peak_ca_i_uM": "peak_cytosolic_calcium_uM",
        "peak_ca_mito_uM": "peak_mitochondrial_calcium_uM",
        "peak_ca_submembrane_uM": "peak_submembrane_calcium_uM",
        "peak_calpain_activity": "peak_calpain_activity",
        "peak_docked_MVB_pool": "peak_docked_multivesicular_body_pool",
        "peak_lysosomal_routing": "peak_lysosomal_routing_bias",
        "peak_osmotic_stress": "peak_osmotic_stress_index",
        "peak_pore_activation": "peak_pore_activation_state",
        "peak_ps_exposure": "peak_phosphatidylserine_exposure",
        "peak_repair_state": "peak_resealing_state",
        "peak_ros": "peak_reactive_oxygen_species",
        "peak_secretory_bias": "peak_secretory_routing_bias",
        "pore_density_m2": "pore_density_per_m2",
        "pore_density_um2": "pore_density_per_um2",
        "pulse_number": "pulse_number",
        "pulse_width_ns": "pulse_width_ns",
        "repetition_rate_Hz": "repetition_rate_Hz",
        "scenario": "scenario",
        "scenario_label": "scenario_label",
        "schwan_limit_V": "schwan_limit_voltage_V",
        "thermal_retention_factor": "thermal_retention_factor",
        "temp_rise_C": "temperature_rise_C",
        "temperature_C": "temperature_C",
        "time_s": "time_s",
        "train_duration_s": "pulse_train_duration_s",
        "waveform": "waveform",
        "waveform_energy_factor": "waveform_energy_factor",
        "min_atp": "minimum_adenosine_triphosphate_state",
        "cumulative_small_EV": "cumulative_small_extracellular_vesicle_output",
        "cumulative_medium_large_EV": "cumulative_medium_large_extracellular_vesicle_output",
        "cumulative_apoptotic_body": "cumulative_apoptotic_body_output",
    },
)
