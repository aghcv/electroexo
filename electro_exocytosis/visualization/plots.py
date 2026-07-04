from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from electro_exocytosis.abbreviations import STANDARD_ABBREVIATIONS
from electro_exocytosis.visualization.style import (
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    line_styles,
    save_manuscript_figure,
)

if TYPE_CHECKING:
    from electro_exocytosis.simulation import SimulationResult


def plot_calcium_timeseries(result: SimulationResult, outdir: Path) -> None:
    """Plot cytosolic, ER, and mitochondrial calcium trajectories."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    ax_er = ax.twinx()
    styles = line_styles(3)
    ax.plot(
        result.state_timeseries["t"],
        result.state_timeseries["Ca_i"],
        label=STANDARD_ABBREVIATIONS.plot_label("Ca_i"),
        **styles[0],
    )
    ax.plot(
        result.state_timeseries["t"],
        result.state_timeseries["Ca_mito"],
        label=STANDARD_ABBREVIATIONS.plot_label("Ca_mito"),
        **styles[2],
    )
    ax_er.plot(
        result.state_timeseries["t"],
        result.state_timeseries["Ca_ER"],
        label=STANDARD_ABBREVIATIONS.plot_label("Ca_ER"),
        **styles[1],
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Cytosolic and mitochondrial calcium (uM)")
    ax_er.set_ylabel("Endoplasmic-reticulum calcium (uM)")
    lines = ax.get_lines() + ax_er.get_lines()
    ax.legend(lines, [line.get_label() for line in lines])
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "calcium_timeseries.png", abbreviation_keys=("ER",))
    plt.close(fig)


def plot_layer3_state_panel(result: SimulationResult, outdir: Path) -> None:
    """Plot ion transport, organelle stress, ROS, and ATP diagnostics."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    t = result.state_timeseries["t"]

    calcium_styles = line_styles(3)
    ca_ax = axes[0, 0]
    er_ax = ca_ax.twinx()
    ca_ax.plot(t, result.state_timeseries["Ca_i"], label=STANDARD_ABBREVIATIONS.plot_label("Ca_i"), **calcium_styles[0])
    ca_ax.plot(
        t,
        result.state_timeseries["Ca_mito"],
        label=STANDARD_ABBREVIATIONS.plot_label("Ca_mito"),
        **calcium_styles[2],
    )
    er_ax.plot(t, result.state_timeseries["Ca_ER"], label=STANDARD_ABBREVIATIONS.plot_label("Ca_ER"), **calcium_styles[1])
    ca_ax.set_ylabel("Cytosolic and mitochondrial calcium (uM)")
    er_ax.set_ylabel("Endoplasmic-reticulum calcium (uM)")
    calcium_lines = ca_ax.get_lines() + er_ax.get_lines()
    ca_ax.legend(calcium_lines, [line.get_label() for line in calcium_lines])

    ion_styles = line_styles(4)
    ion_ax = axes[0, 1]
    osmotic_ax = ion_ax.twinx()
    ion_ax.plot(t, result.state_timeseries["Na_i"], label="Intracellular sodium", **ion_styles[0])
    ion_ax.plot(t, result.state_timeseries["K_i"], label="Intracellular potassium", **ion_styles[1])
    ion_ax.plot(t, result.state_timeseries["Cl_i"], label="Intracellular chloride", **ion_styles[2])
    osmotic_ax.plot(t, result.state_timeseries["osmotic_stress"], label="osmotic stress", **ion_styles[3])
    ion_ax.set_ylabel("Ion concentration (mM)")
    osmotic_ax.set_ylabel("Osmotic stress")
    ion_lines = ion_ax.get_lines() + osmotic_ax.get_lines()
    ion_ax.legend(ion_lines, [line.get_label() for line in ion_lines])

    mito_styles = line_styles(2)
    axes[1, 0].plot(t, result.state_timeseries["mitochondrial_potential"], label="Mitochondrial membrane potential", **mito_styles[0])
    axes[1, 0].plot(t, result.state_timeseries["ATP"], label=STANDARD_ABBREVIATIONS.plot_label("ATP"), **mito_styles[1])
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("Normalized state")
    axes[1, 0].legend()

    stress_styles = line_styles(2)
    axes[1, 1].plot(t, result.state_timeseries["ROS"], label=STANDARD_ABBREVIATIONS.plot_label("ROS"), **stress_styles[0])
    axes[1, 1].plot(t, result.state_timeseries["damage"], label="damage", **stress_styles[1])
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("Normalized state")
    axes[1, 1].legend()

    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "layer3_state_panel.png", abbreviation_keys=("ER", "ATP", "ROS"))
    plt.close(fig)


def plot_remodeling_repair_panel(result: SimulationResult, outdir: Path) -> None:
    """Plot Ca2+-dependent remodeling and repair diagnostics."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    t = result.state_timeseries["t"]

    calcium_styles = line_styles(2)
    axes[0, 0].plot(t, result.state_timeseries["Ca_i"], label=STANDARD_ABBREVIATIONS.plot_label("Ca_i"), **calcium_styles[0])
    axes[0, 0].plot(
        t,
        result.state_timeseries["Ca_submembrane"],
        label=STANDARD_ABBREVIATIONS.plot_label("Ca_submembrane"),
        **calcium_styles[1],
    )
    axes[0, 0].set_ylabel("Calcium (uM)")
    axes[0, 0].legend()

    lipid_styles = line_styles(3)
    axes[0, 1].plot(
        t,
        result.state_timeseries["PS_exposure"],
        label=STANDARD_ABBREVIATIONS.plot_label("PS_exposure"),
        **lipid_styles[0],
    )
    axes[0, 1].plot(t, result.state_timeseries["scramblase_activity"], label="scramblase", **lipid_styles[1])
    axes[0, 1].plot(t, result.state_timeseries["flippase_activity"], label="flippase", **lipid_styles[2])
    axes[0, 1].set_ylabel("Normalized state")
    axes[0, 1].legend()

    effector_styles = line_styles(3)
    axes[1, 0].plot(t, result.state_timeseries["calpain_activity"], label="calpain", **effector_styles[0])
    axes[1, 0].plot(t, result.state_timeseries["annexin_activity"], label="annexin", **effector_styles[1])
    axes[1, 0].plot(t, result.state_timeseries["lysosomal_repair_activity"], label="lysosomal repair", **effector_styles[2])
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("Effector activity")
    axes[1, 0].legend()

    repair_styles = line_styles(4)
    axes[1, 1].plot(t, result.state_timeseries["actin_disruption"], label="actin disruption", **repair_styles[0])
    axes[1, 1].plot(t, result.state_timeseries["actomyosin_tension"], label="actomyosin", **repair_styles[1])
    axes[1, 1].plot(t, result.state_timeseries["repair_state"], label="resealing", **repair_styles[2])
    axes[1, 1].plot(t, result.state_timeseries["repair_shedding_rate"], label="shedding rate", **repair_styles[3])
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("Normalized state / rate")
    axes[1, 1].legend()

    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "remodeling_repair_panel.png", abbreviation_keys=("PS",))
    plt.close(fig)



def plot_ev_biogenesis_panel(result: SimulationResult, outdir: Path) -> None:
    """Plot Layer 5 EV biogenesis and subtype-release diagnostics."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    t = result.ev_timeseries["t"]

    pool_styles = line_styles(3)
    axes[0, 0].plot(t, result.ev_timeseries["MVB_pool"], label=STANDARD_ABBREVIATIONS.plot_label("MVB_pool"), **pool_styles[0])
    axes[0, 0].plot(t, result.ev_timeseries["ILV_load"], label=STANDARD_ABBREVIATIONS.plot_label("ILV_load"), **pool_styles[1])
    axes[0, 0].plot(t, result.ev_timeseries["docked_MVB_pool"], label="Docked multivesicular-body pool", **pool_styles[2])
    axes[0, 0].set_ylabel("Relative pool size")
    axes[0, 0].legend()

    fate_styles = line_styles(4)
    axes[0, 1].plot(t, result.ev_timeseries["escrt_dependent_signal"], label="ESCRT-dependent signal", **fate_styles[0])
    axes[0, 1].plot(t, result.ev_timeseries["ceramide_signal"], label="ceramide", **fate_styles[1])
    axes[0, 1].plot(t, result.ev_timeseries["secretory_bias"], label="secretory bias", **fate_styles[2])
    axes[0, 1].plot(t, result.ev_timeseries["lysosomal_routing"], label="lysosomal routing", **fate_styles[3])
    axes[0, 1].set_ylabel("Normalized signal")
    axes[0, 1].legend()

    budding_styles = line_styles(4)
    axes[1, 0].plot(t, result.ev_timeseries["budding_pool"], label="budding pool", **budding_styles[0])
    axes[1, 0].plot(t, result.ev_timeseries["budding_signal"], label="budding", **budding_styles[1])
    axes[1, 0].plot(t, result.ev_timeseries["scission_signal"], label="scission", **budding_styles[2])
    axes[1, 0].plot(t, result.ev_timeseries["apoptotic_commitment"], label="apoptotic commitment", **budding_styles[3])
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylabel("Relative state")
    axes[1, 0].legend()

    rate_styles = line_styles(3)
    axes[1, 1].plot(t, result.ev_timeseries["sEV_rate"], label=STANDARD_ABBREVIATIONS.plot_label("sEV_rate"), **rate_styles[0])
    axes[1, 1].plot(t, result.ev_timeseries["mlEV_rate"], label=STANDARD_ABBREVIATIONS.plot_label("mlEV_rate"), **rate_styles[1])
    axes[1, 1].plot(t, result.ev_timeseries["AB_rate"], label=STANDARD_ABBREVIATIONS.plot_label("AB_rate"), **rate_styles[2])
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("Release rate")
    axes[1, 1].legend()

    fig.tight_layout()
    save_manuscript_figure(
        fig,
        outdir / "ev_biogenesis_panel.png",
        abbreviation_keys=("MVB", "ILV", "ESCRT", "sEV", "m/lEV", "AB"),
    )
    plt.close(fig)


def plot_ev_release_rates(result: SimulationResult, outdir: Path) -> None:
    """Plot EV release rates."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(3)
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["sEV_rate"], label=STANDARD_ABBREVIATIONS.plot_label("sEV"), **styles[0])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["mlEV_rate"], label=STANDARD_ABBREVIATIONS.plot_label("m/lEV"), **styles[1])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["AB_rate"], label=STANDARD_ABBREVIATIONS.plot_label("AB"), **styles[2])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rate")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "ev_release_rates.png", abbreviation_keys=("EV", "sEV", "m/lEV", "AB"))
    plt.close(fig)



def plot_cumulative_ev_yield(result: SimulationResult, outdir: Path) -> None:
    """Plot cumulative EV yields."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(3)
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["sEV_cumulative"], label=STANDARD_ABBREVIATIONS.plot_label("sEV"), **styles[0])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["mlEV_cumulative"], label=STANDARD_ABBREVIATIONS.plot_label("m/lEV"), **styles[1])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["AB_cumulative"], label=STANDARD_ABBREVIATIONS.plot_label("AB"), **styles[2])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Cumulative output")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "cumulative_ev_yield.png", abbreviation_keys=("EV", "sEV", "m/lEV", "AB"))
    plt.close(fig)



def plot_quality_viability(result: SimulationResult, outdir: Path) -> None:
    """Plot damage and viability proxy."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(2)
    ax.plot(result.state_timeseries["t"], result.state_timeseries["damage"], label="Damage", **styles[0])
    ax.axhline(result.summary["viability_fraction"], label="viability fraction", **line_styles(2, include_markers=False)[1])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalized value")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "quality_viability.png")
    plt.close(fig)



def generate_all_plots(result: SimulationResult, outdir: Path) -> None:
    """Generate all standard plots."""
    plot_calcium_timeseries(result, outdir)
    plot_layer3_state_panel(result, outdir)
    plot_remodeling_repair_panel(result, outdir)
    plot_ev_biogenesis_panel(result, outdir)
    plot_ev_release_rates(result, outdir)
    plot_cumulative_ev_yield(result, outdir)
    plot_quality_viability(result, outdir)
