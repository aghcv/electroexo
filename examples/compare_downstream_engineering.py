from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from electro_exocytosis.abbreviations import STANDARD_ABBREVIATIONS
from electro_exocytosis.config import (
    CellStateConfig,
    ExposureConfig,
    PulseConfig,
    ScenarioConfig,
    SimulationConfig,
    SimulationScenario,
)
from electro_exocytosis.simulation import Simulation
from electro_exocytosis.visualization.style import (
    MANUSCRIPT_DOUBLE_COLUMN_WIDTH_IN,
    bar_colors,
    bar_hatch,
    manuscript_style_context,
    save_manuscript_figure,
)


@dataclass(frozen=True)
class DownstreamScenarioSpec:
    name: str
    label: str
    pulse_overrides: dict = field(default_factory=dict)
    cell_state_overrides: dict = field(default_factory=dict)
    parameter_overrides: dict = field(default_factory=dict)


SCENARIOS = (
    DownstreamScenarioSpec(
        name="productive_cargo_window",
        label="Productive cargo",
        pulse_overrides={"amplitude_kV_cm": 15.0, "pulse_number": 20},
        parameter_overrides={
            "ev_release": {
                "ceramide_baseline": 0.50,
                "acidification_ceramide_relief": 0.75,
                "baseline_sEV_rate": 1.25,
            },
            "cargo_potency": {
                "ESCRT_protein_weight": 0.45,
                "rbp_RNA_sorting_weight": 0.45,
                "ceramide_lipid_weight": 0.45,
            },
        },
    ),
    DownstreamScenarioSpec(
        name="direct_loading_mode",
        label="Direct loading",
        pulse_overrides={"amplitude_kV_cm": 12.0, "pulse_number": 12},
        parameter_overrides={
            "cargo_potency": {
                "direct_loading_efficiency": 0.70,
                "direct_loading_leak_fraction": 0.05,
                "potency_weights": {
                    "protein": 0.20,
                    "RNA": 0.20,
                    "lipid": 0.15,
                    "antigen": 0.15,
                    "direct_load": 0.30,
                },
            },
            "manufacturing_qc": {"isolation_efficiency": 0.45},
        },
    ),
    DownstreamScenarioSpec(
        name="injury_quality_penalty",
        label="Injury penalty",
        pulse_overrides={"amplitude_kV_cm": 24.0, "pulse_number": 40},
        cell_state_overrides={"stress_sensitivity_modifier": 1.6},
        parameter_overrides={
            "ev_release": {"baseline_AB_rate": 0.035, "apoptosis_damage_weight": 0.70},
            "injury_quality": {"contamination_threshold": 0.25},
        },
    ),
)


def build_response_table() -> pd.DataFrame:
    rows: list[dict[str, float | str | bool]] = []
    for spec in SCENARIOS:
        result = Simulation(
            _build_scenario(spec), params_override=spec.parameter_overrides
        ).run()
        terminal_quality = result.parameters_used["terminal_quality"]
        manufacturing = result.parameters_used["manufacturing"]
        rows.append(
            {
                "scenario": spec.name,
                "label": spec.label,
                "protein_enrichment": result.summary["protein_enrichment"],
                "RNA_enrichment": result.summary["RNA_enrichment"],
                "lipid_enrichment": result.summary["lipid_enrichment"],
                "antigen_enrichment": result.summary["antigen_enrichment"],
                "direct_loaded_cargo": result.summary["direct_loaded_cargo"],
                "potency_score": result.summary["potency_score"],
                "viability_fraction": result.summary["viability_fraction"],
                "apoptosis_fraction": result.summary["apoptosis_fraction"],
                "necrosis_fraction": result.summary["necrosis_fraction"],
                "bona_fide_EV_fraction": result.summary["bona_fide_EV_fraction"],
                "quality_pass": bool(terminal_quality["quality_pass"]),
                "cell_normalized_yield": result.summary["cell_normalized_yield"],
                "process_recovery": manufacturing["process_recovery"],
                "purity_score": result.summary["purity_score"],
                "batch_adjusted_yield": result.summary["batch_adjusted_yield"],
                "optimization_objective": result.summary["optimization_objective"],
            }
        )
    return pd.DataFrame(rows)


def write_outputs(summary: pd.DataFrame, outdir: Path, make_plots: bool = True) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    STANDARD_ABBREVIATIONS.rename_columns(summary).to_csv(
        outdir / "downstream_engineering_summary.csv",
        index=False,
    )
    STANDARD_ABBREVIATIONS.write_bundle(
        outdir, keys=("EV", "RNA", "sEV", "m/lEV", "AB")
    )
    if make_plots:
        _plot_downstream_panel(summary, outdir)


def _build_scenario(spec: DownstreamScenarioSpec) -> SimulationScenario:
    pulse_kwargs = {
        "amplitude_kV_cm": 15.0,
        "pulse_width_ns": 200.0,
        "pulse_number": 20,
        "repetition_rate_Hz": 5.0,
    }
    pulse_kwargs.update(spec.pulse_overrides)
    cell_state_kwargs = {"cell_type": "generic"}
    cell_state_kwargs.update(spec.cell_state_overrides)
    mode = (
        "direct_EV_engineering"
        if spec.name == "direct_loading_mode"
        else "cell_based_electro_exocytosis"
    )
    return SimulationScenario(
        scenario=ScenarioConfig(name=spec.name, mode=mode),
        pulse=PulseConfig(**pulse_kwargs),
        exposure=ExposureConfig(dosimetry_model="joule_lumped_thermal"),
        cell_state=CellStateConfig(**cell_state_kwargs),
        simulation=SimulationConfig(t_start_s=0.0, t_end_s=1800.0, output_dt_s=5.0),
    )


def _plot_downstream_panel(summary: pd.DataFrame, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_labels = ("Productive", "Direct loading", "Injury shifted")
    colors = bar_colors(len(summary))
    hatches = [bar_hatch(index) for index in range(len(summary))]

    panel_specs = (
        (
            "Layer 6: cargo state",
            (
                ("protein_enrichment", "Protein"),
                ("RNA_enrichment", "RNA"),
                ("lipid_enrichment", "Lipid"),
                ("antigen_enrichment", "Antigen"),
                ("direct_loaded_cargo", "Direct loading"),
            ),
        ),
        (
            "Layer 7: quality gate",
            (
                ("viability_fraction", "Viability"),
                ("bona_fide_EV_fraction", "Bona fide EV"),
                ("purity_score", "Purity"),
            ),
        ),
        (
            "Layer 8: process objective",
            (
                ("optimization_objective", "Optimization"),
                ("normalized_batch_yield", "Batch yield"),
            ),
        ),
    )
    plot_frame = summary.copy()
    batch_yield = plot_frame["batch_adjusted_yield"].to_numpy(dtype=float)
    plot_frame["normalized_batch_yield"] = batch_yield / max(
        float(batch_yield.max()), 1e-12
    )

    with manuscript_style_context():
        fig, axes = plt.subplots(
            3,
            1,
            figsize=(MANUSCRIPT_DOUBLE_COLUMN_WIDTH_IN, 3.65),
            gridspec_kw={"height_ratios": (1.45, 1.0, 0.78)},
        )
        bar_height = 0.22
        for axis, (title, metric_specs) in zip(axes, panel_specs, strict=True):
            y_positions = np.arange(len(metric_specs), dtype=float)
            for scenario_index, (_, row) in enumerate(plot_frame.iterrows()):
                offsets = y_positions + (scenario_index - 1) * bar_height
                values = [float(row[metric]) for metric, _ in metric_specs]
                axis.barh(
                    offsets,
                    values,
                    height=bar_height,
                    color=colors[scenario_index],
                    hatch=hatches[scenario_index],
                    edgecolor="#222222",
                    linewidth=0.55,
                    label=plot_labels[scenario_index],
                )
            axis.set_title(title, loc="left", pad=3)
            axis.set_yticks(
                y_positions,
                [label for _, label in metric_specs],
            )
            axis.invert_yaxis()
            panel_max = max(
                float(plot_frame[metric].max()) for metric, _ in metric_specs
            )
            axis.set_xlim(0.0, panel_max * 1.08)
            axis.grid(axis="x", color="#D9D9D9", linewidth=0.6, alpha=0.65)
            axis.tick_params(labelsize=8.5)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=3,
            fontsize=8,
            frameon=True,
            facecolor="white",
            edgecolor="#555555",
            framealpha=0.96,
        )
        fig.supxlabel("Relative score", fontsize=9.5, y=0.015)
        fig.subplots_adjust(left=0.17, right=0.985, bottom=0.10, top=0.83, hspace=0.72)
        save_manuscript_figure(fig, outdir / "downstream_engineering_panel.png")
        plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=Path("results/downstream_engineering_comparison")
    )
    parser.add_argument(
        "--no-plots", action="store_true", help="Skip PNG plot generation."
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary = build_response_table()
    write_outputs(summary, args.out, make_plots=not args.no_plots)
    print(f"Wrote downstream engineering comparison outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
