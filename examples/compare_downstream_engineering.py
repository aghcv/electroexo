from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

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
    MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE,
    bar_colors,
    bar_hatch,
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
        result = Simulation(_build_scenario(spec), params_override=spec.parameter_overrides).run()
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
    STANDARD_ABBREVIATIONS.write_bundle(outdir, keys=("EV", "RNA", "sEV", "m/lEV", "AB"))
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
    mode = "direct_EV_engineering" if spec.name == "direct_loading_mode" else "cell_based_electro_exocytosis"
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

    fig, axes = plt.subplots(1, 3, figsize=MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE)
    labels = summary["label"].tolist()
    x = range(len(summary))
    colors = bar_colors(5)

    cargo_metrics = ["protein_enrichment", "RNA_enrichment", "lipid_enrichment", "antigen_enrichment", "direct_loaded_cargo"]
    width = 0.16
    for index, metric in enumerate(cargo_metrics):
        values = summary[metric].to_numpy(dtype=float)
        positions = [pos + (index - 2) * width for pos in x]
        axes[0].bar(positions, values, width=width, color=colors[index], hatch=bar_hatch(index), label=STANDARD_ABBREVIATIONS.plot_label(metric))
    axes[0].set_title("Layer 6 cargo state")
    axes[0].set_ylabel("Relative score")
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(labels, rotation=20, ha="right")
    axes[0].legend(fontsize=7)

    quality_metrics = ["viability_fraction", "bona_fide_EV_fraction", "purity_score"]
    width = 0.22
    for index, metric in enumerate(quality_metrics):
        values = summary[metric].to_numpy(dtype=float)
        positions = [pos + (index - 1) * width for pos in x]
        axes[1].bar(positions, values, width=width, color=colors[index], hatch=bar_hatch(index), label=STANDARD_ABBREVIATIONS.plot_label(metric))
    axes[1].set_title("Layer 7 quality gate")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(labels, rotation=20, ha="right")
    axes[1].legend(fontsize=7)

    objective_values = summary["optimization_objective"].to_numpy(dtype=float)
    yield_values = summary["batch_adjusted_yield"].to_numpy(dtype=float)
    yield_values = yield_values / max(float(yield_values.max()), 1e-12)
    axes[2].bar([pos - 0.12 for pos in x], objective_values, width=0.24, color=colors[0], hatch=bar_hatch(0), label="Optimization objective")
    axes[2].bar([pos + 0.12 for pos in x], yield_values, width=0.24, color=colors[1], hatch=bar_hatch(1), label="Batch-adjusted yield")
    axes[2].set_title("Layer 8 process objective")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_xticks(list(x))
    axes[2].set_xticklabels(labels, rotation=20, ha="right")
    axes[2].legend(fontsize=7)

    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "downstream_engineering_panel.png", abbreviation_keys=("EV", "RNA"))
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/downstream_engineering_comparison"))
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary = build_response_table()
    write_outputs(summary, args.out, make_plots=not args.no_plots)
    print(f"Wrote downstream engineering comparison outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
