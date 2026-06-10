from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from electro_exocytosis.visualization.style import (
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    line_styles,
    save_manuscript_figure,
)

if TYPE_CHECKING:
    from electro_exocytosis.simulation import SimulationResult


def plot_calcium_timeseries(result: SimulationResult, outdir: Path) -> None:
    """Plot calcium and ER calcium trajectories."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(2)
    ax.plot(result.state_timeseries["t"], result.state_timeseries["Ca_i"], label="Ca_i", **styles[0])
    ax.plot(result.state_timeseries["t"], result.state_timeseries["Ca_ER"], label="Ca_ER", **styles[1])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Calcium (uM)")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "calcium_timeseries.png")
    plt.close(fig)



def plot_ev_release_rates(result: SimulationResult, outdir: Path) -> None:
    """Plot EV release rates."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(3)
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["sEV_rate"], label="sEV", **styles[0])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["mlEV_rate"], label="mlEV", **styles[1])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["AB_rate"], label="AB", **styles[2])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rate")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "ev_release_rates.png")
    plt.close(fig)



def plot_cumulative_ev_yield(result: SimulationResult, outdir: Path) -> None:
    """Plot cumulative EV yields."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(3)
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["sEV_cumulative"], label="sEV", **styles[0])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["mlEV_cumulative"], label="mlEV", **styles[1])
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["AB_cumulative"], label="AB", **styles[2])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Cumulative output")
    ax.legend()
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "cumulative_ev_yield.png")
    plt.close(fig)



def plot_quality_viability(result: SimulationResult, outdir: Path) -> None:
    """Plot damage and viability proxy."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)
    styles = line_styles(2)
    ax.plot(result.state_timeseries["t"], result.state_timeseries["damage"], label="damage", **styles[0])
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
    plot_ev_release_rates(result, outdir)
    plot_cumulative_ev_yield(result, outdir)
    plot_quality_viability(result, outdir)
