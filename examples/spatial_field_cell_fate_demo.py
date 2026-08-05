from __future__ import annotations

import argparse
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib")
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import yaml

from electro_exocytosis.config import SimulationScenario
from electro_exocytosis.simulation import Simulation, SimulationResult
from electro_exocytosis.spatial import (
    FieldGeometry2D,
    FieldSolution2D,
    sample_field_kV_cm,
    solve_quasistatic_field,
)


DEFAULT_OUTDIR = Path(
    "results/bioelectrics_seed_2026_aim2/spatial_cell_fate_demo"
)


@dataclass(frozen=True, slots=True)
class SelectedCell:
    label: str
    x_mm: float
    y_mm: float
    color: str
    linestyle: str


SELECTED_CELLS = (
    SelectedCell("A", -1.25, 0.35, "#D55E00", "-"),
    SelectedCell("B", 0.30, 0.80, "#0072B2", "--"),
    SelectedCell("C", 3.70, 2.70, "#009E73", "-."),
    SelectedCell("D", 3.80, -2.80, "#CC79A7", ":"),
)


def build_scenario(cell: SelectedCell, local_field_kV_cm: float) -> dict[str, Any]:
    return {
        "scenario": {
            "name": f"spatial_cell_{cell.label}_{local_field_kV_cm:.3f}_kV_cm",
            "mode": "cell_based_electro_exocytosis",
        },
        "pulse": {
            "amplitude_kV_cm": max(float(local_field_kV_cm), 1e-6),
            "pulse_width_ns": 60.0,
            "pulse_number": 10,
            "repetition_rate_Hz": 1.0,
            "waveform": "square",
        },
        "exposure": {
            "geometry": "dish",
            "medium_conductivity_S_m": 0.20,
            "temperature_C": 37.0,
            "cell_density_per_ml": 1_000_000,
            "dosimetry_model": "joule_lumped_thermal",
        },
        "cell_state": {
            "cell_type": "N1S1-HCC provisional",
            "membrane_modifier": 1.0,
            "calcium_handling_modifier": 1.0,
            "baseline_EV_release_modifier": 1.0,
            "stress_sensitivity_modifier": 1.0,
        },
        "simulation": {
            "t_start_s": 0.0,
            "t_end_s": 1_800.0,
            "output_dt_s": 5.0,
            "numerical_method": "solve_ivp",
        },
    }


def run_location_scenarios(
    solution: FieldSolution2D,
    outdir: Path,
) -> tuple[pd.DataFrame, dict[str, SimulationResult]]:
    points = [(cell.x_mm, cell.y_mm) for cell in SELECTED_CELLS]
    local_fields = sample_field_kV_cm(solution, points)
    scenario_dir = outdir / "scenarios"
    scenario_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    results: dict[str, SimulationResult] = {}
    for cell, local_field in zip(SELECTED_CELLS, local_fields, strict=True):
        x_index = int(np.argmin(np.abs(solution.x_mm - cell.x_mm)))
        y_index = int(np.argmin(np.abs(solution.y_mm - cell.y_mm)))
        if solution.electrode_mask[y_index, x_index] or solution.vessel_mask[y_index, x_index]:
            raise ValueError(f"Selected Cell {cell.label} is not located in tissue.")
        if not np.isfinite(local_field) or local_field <= 0:
            raise ValueError(f"Selected Cell {cell.label} has an invalid local field.")
        payload = build_scenario(cell, float(local_field))
        scenario_path = scenario_dir / f"cell_{cell.label}.yaml"
        with scenario_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False)

        result = Simulation(SimulationScenario.model_validate(payload)).run()
        results[cell.label] = result
        rows.append(
            {
                "cell": cell.label,
                "x_mm": cell.x_mm,
                "y_mm": cell.y_mm,
                "local_field_kV_cm": float(local_field),
                "peak_cytosolic_calcium_uM": float(result.summary["peak_ca_i"]),
                "peak_reactive_oxygen_species": float(result.summary["peak_ros"]),
                "minimum_ATP_state": float(result.summary["min_atp"]),
                "peak_repair_state": float(result.summary["peak_repair_state"]),
                "peak_damage_state": float(result.state_timeseries["damage"].max()),
                "terminal_damage_state": float(result.state_timeseries["damage"].iloc[-1]),
                "viability_fraction": float(result.summary["viability_fraction"]),
                "cumulative_small_EV": float(result.summary["cumulative_small_EV"]),
                "color": cell.color,
            }
        )
    return pd.DataFrame(rows), results


def write_numerical_outputs(
    solution: FieldSolution2D,
    summary: pd.DataFrame,
    results: dict[str, SimulationResult],
    outdir: Path,
) -> None:
    summary.to_csv(outdir / "selected_cell_summary.csv", index=False)

    state_frames: list[pd.DataFrame] = []
    for cell in SELECTED_CELLS:
        frame = results[cell.label].state_timeseries.copy()
        frame.insert(0, "cell", cell.label)
        state_frames.append(frame)
    pd.concat(state_frames, ignore_index=True).to_csv(
        outdir / "selected_cell_state_timeseries.csv", index=False
    )

    np.savez_compressed(
        outdir / "quasistatic_field_solution.npz",
        x_mm=solution.x_mm,
        y_mm=solution.y_mm,
        conductivity_S_m=solution.conductivity_S_m,
        potential_V=solution.potential_V,
        field_x_kV_cm=solution.field_x_kV_cm,
        field_y_kV_cm=solution.field_y_kV_cm,
        field_magnitude_kV_cm=solution.field_magnitude_kV_cm,
        electrode_mask=solution.electrode_mask,
        vessel_mask=solution.vessel_mask,
    )

    manifest = {
        "analysis": "2D quasi-static field to location-specific electroexo scenarios",
        "field_equation": "div(sigma grad(phi)) = 0",
        "pulse_protocol": {
            "pulse_width_ns": 60.0,
            "pulse_number": 10,
            "repetition_rate_Hz": 1.0,
            "waveform": "square",
        },
        "applied_voltage_V": solution.applied_voltage_V,
        "reference_field_kV_cm": solution.geometry.reference_field_kV_cm,
        "tissue_conductivity_S_m": solution.geometry.tissue_conductivity_S_m,
        "vessel_conductivity_S_m": solution.geometry.vessel_conductivity_S_m,
        "selected_cells": summary.to_dict(orient="records"),
        "model_status": (
            "Computed proof of concept. Field solve is quasi-static and the biological "
            "model remains exploratory and uncalibrated."
        ),
    }
    (outdir / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


def plot_linked_figure(
    solution: FieldSolution2D,
    summary: pd.DataFrame,
    results: dict[str, SimulationResult],
    outdir: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelcolor": "#202124",
            "axes.titlecolor": "#202124",
            "xtick.color": "#202124",
            "ytick.color": "#202124",
        }
    )

    fig = plt.figure(figsize=(12.6, 6.8), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.38, 1.0),
        height_ratios=(0.82, 1.18),
        left=0.065,
        right=0.965,
        bottom=0.12,
        top=0.955,
        wspace=0.29,
        hspace=0.34,
    )
    field_ax = fig.add_subplot(grid[:, 0])
    pulse_ax = fig.add_subplot(grid[0, 1])
    fate_ax = fig.add_subplot(grid[1, 1])

    _plot_field_map(field_ax, solution, summary)
    _plot_local_pulses(pulse_ax, summary)
    _plot_damage_trajectories(fate_ax, summary, results)

    fig.text(
        0.5,
        0.035,
        "Computed proof of concept: the 2D solve is quasi-static and the current cell model is uncalibrated. "
        "Colors link each mapped cell to its local pulse and predicted injury trajectory.",
        ha="center",
        va="center",
        fontsize=7.3,
        color="#687078",
    )
    fig.savefig(
        outdir / "spatial_field_cell_fate_linked.png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        outdir / "spatial_field_cell_fate_linked.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def _plot_field_map(
    ax: plt.Axes,
    solution: FieldSolution2D,
    summary: pd.DataFrame,
) -> None:
    xx_mm, yy_mm = np.meshgrid(solution.x_mm, solution.y_mm)
    field = np.ma.masked_less_equal(solution.field_magnitude_kV_cm, 0.12)
    filled_levels = np.geomspace(0.12, 100.0, 28)
    contour = ax.contourf(
        xx_mm,
        yy_mm,
        field,
        levels=filled_levels,
        norm=LogNorm(vmin=filled_levels[0], vmax=filled_levels[-1]),
        cmap="viridis",
        extend="both",
    )
    line_levels = (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0)
    lines = ax.contour(
        xx_mm,
        yy_mm,
        solution.field_magnitude_kV_cm,
        levels=line_levels,
        colors="white",
        linewidths=0.72,
        alpha=0.78,
    )
    ax.clabel(lines, fmt=lambda value: f"{value:g}", fontsize=6.4, inline=True)

    ax.contourf(
        xx_mm,
        yy_mm,
        solution.vessel_mask.astype(float),
        levels=(0.5, 1.5),
        colors=("#DDE4E9",),
        alpha=0.58,
        zorder=4,
    )
    ax.contour(
        xx_mm,
        yy_mm,
        solution.vessel_mask.astype(float),
        levels=(0.5,),
        colors=("#4F5A62",),
        linewidths=1.0,
        zorder=5,
    )
    electrode_x, electrode_y = solution.geometry.electrode_center_mm
    electrode = plt.Circle(
        (electrode_x, electrode_y),
        solution.geometry.electrode_radius_mm,
        facecolor="#111111",
        edgecolor="white",
        linewidth=0.8,
        zorder=7,
    )
    ax.add_patch(electrode)

    legend_handles: list[Line2D | Patch] = []
    for cell in SELECTED_CELLS:
        record = summary.loc[summary["cell"] == cell.label].iloc[0]
        marker_size = 65.0 + 460.0 * float(record["peak_damage_state"])
        ax.scatter(
            [cell.x_mm],
            [cell.y_mm],
            s=marker_size,
            c=[cell.color],
            edgecolors="white",
            linewidths=1.45,
            zorder=8,
            clip_on=False,
        )
        ax.text(
            cell.x_mm,
            cell.y_mm,
            cell.label,
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="bold",
            color="white",
            zorder=9,
            clip_on=False,
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=cell.color,
                marker="o",
                linestyle=cell.linestyle,
                linewidth=1.6,
                markersize=6,
                label=(
                    f"Cell {cell.label}: {float(record['local_field_kV_cm']):.2f} kV/cm; "
                    f"viability {float(record['viability_fraction']):.2f}"
                ),
            )
        )
    legend_handles.append(
        Patch(
            facecolor="#DDE4E9",
            edgecolor="#4F5A62",
            label="Conductive vessel (5x tissue conductivity)",
        )
    )

    ax.annotate(
        f"needle electrode\n{solution.applied_voltage_V / 1000.0:.2f} kV",
        xy=(electrode_x, electrode_y),
        xytext=(-1.75, -0.68),
        ha="center",
        va="center",
        fontsize=6.7,
        color="#202124",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        arrowprops={"arrowstyle": "-", "color": "#202124", "linewidth": 0.8},
        zorder=9,
    )
    ax.text(
        1.35,
        1.28,
        "conductive vessel",
        ha="center",
        va="center",
        fontsize=7.0,
        color="#202124",
        zorder=9,
    )
    ax.legend(
        handles=legend_handles,
        loc="lower left",
        fontsize=6.6,
        frameon=True,
        framealpha=0.94,
        borderpad=0.7,
        handlelength=2.1,
    )
    ax.set_aspect("equal")
    ax.set_xlim(*solution.geometry.x_limits_mm)
    ax.set_ylim(*solution.geometry.y_limits_mm)
    ax.set_xlabel("Distance (mm)")
    ax.set_ylabel("Distance (mm)")
    ax.set_title(
        "A   Computed 2D field and sampled cell locations",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.tick_params(labelsize=7.4)
    cbar = plt.colorbar(contour, ax=ax, fraction=0.046, pad=0.035)
    cbar.set_label("Local field magnitude (kV/cm; log scale)", fontsize=7.6)
    cbar.set_ticks((0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0))
    cbar.set_ticklabels(("0.2", "0.5", "1", "2", "5", "10", "20", "40", "80"))
    cbar.ax.tick_params(labelsize=6.8)
    ax.text(
        0.5,
        -0.105,
        "Cell marker area is proportional to peak predicted damage; color links location to panels B-C.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=6.8,
        color="#687078",
    )


def _plot_local_pulses(ax: plt.Axes, summary: pd.DataFrame) -> None:
    time_ns = np.linspace(-15.0, 105.0, 900)
    on_mask = (time_ns >= 0.0) & (time_ns <= 60.0)
    for cell in SELECTED_CELLS:
        record = summary.loc[summary["cell"] == cell.label].iloc[0]
        amplitude = float(record["local_field_kV_cm"])
        waveform = np.zeros_like(time_ns)
        waveform[on_mask] = amplitude
        ax.step(
            time_ns,
            waveform,
            where="post",
            color=cell.color,
            linestyle=cell.linestyle,
            linewidth=1.65,
            label=f"{cell.label}: {amplitude:.2f} kV/cm",
        )
    ax.axvspan(0.0, 60.0, color="#E5E8EA", alpha=0.36, zorder=0)
    ax.set_yscale("symlog", linthresh=0.20, linscale=1.0)
    ax.set_xlim(-10.0, 100.0)
    ax.set_ylim(0.0, 70.0)
    ax.set_xlabel("Time within representative pulse (ns)")
    ax.set_ylabel("Local field (kV/cm)")
    ax.set_title(
        "B   Local nsPEF signal sampled at each cell",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
    )
    ax.grid(True, color="#E5E8EA", linewidth=0.7)
    ax.tick_params(labelsize=7.2)
    ax.legend(fontsize=6.4, loc="upper right", frameon=True, ncol=2)


def _plot_damage_trajectories(
    ax: plt.Axes,
    summary: pd.DataFrame,
    results: dict[str, SimulationResult],
) -> None:
    for cell in SELECTED_CELLS:
        result = results[cell.label]
        record = summary.loc[summary["cell"] == cell.label].iloc[0]
        time_min = result.state_timeseries["t"].to_numpy(dtype=float) / 60.0
        damage = result.state_timeseries["damage"].to_numpy(dtype=float)
        ax.plot(
            time_min,
            damage,
            color=cell.color,
            linestyle=cell.linestyle,
            linewidth=2.0,
            label=(
                f"{cell.label}: peak {float(record['peak_damage_state']):.2f}; "
                f"viability {float(record['viability_fraction']):.2f}"
            ),
        )
    pulse_train_end_min = 10.0 / 60.0
    ax.axvline(
        pulse_train_end_min,
        color="#687078",
        linewidth=0.9,
        linestyle=(0, (2, 3)),
    )
    ax.text(
        pulse_train_end_min + 0.3,
        0.315,
        "10-pulse train ends",
        fontsize=6.7,
        color="#687078",
        va="top",
    )
    ax.set_xlim(0.0, 30.0)
    ax.set_ylim(0.0, 0.33)
    ax.set_xlabel("Time after pulse-train onset (min)")
    ax.set_ylabel("Predicted damage state (0-1)")
    ax.set_title(
        "C   Location-specific cell-injury trajectories",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
    )
    ax.grid(True, color="#E5E8EA", linewidth=0.7)
    ax.tick_params(labelsize=7.2)
    ax.legend(fontsize=6.5, loc="upper right", frameon=True)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    geometry = FieldGeometry2D(nx=201, ny=151)
    solution = solve_quasistatic_field(geometry)
    summary, results = run_location_scenarios(solution, args.out)
    write_numerical_outputs(solution, summary, results, args.out)
    plot_linked_figure(solution, summary, results, args.out)
    print(f"Wrote spatial field/cell-fate demo to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
