from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from electro_exocytosis.simulation import SimulationResult


def plot_calcium_timeseries(result: SimulationResult, outdir: Path) -> None:
    """Plot calcium and ER calcium trajectories."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    ax.plot(result.state_timeseries["t"], result.state_timeseries["Ca_i"], label="Ca_i")
    ax.plot(result.state_timeseries["t"], result.state_timeseries["Ca_ER"], label="Ca_ER")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Calcium (uM)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "calcium_timeseries.png", dpi=150)
    plt.close(fig)



def plot_ev_release_rates(result: SimulationResult, outdir: Path) -> None:
    """Plot EV release rates."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["sEV_rate"], label="sEV")
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["mlEV_rate"], label="mlEV")
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["AB_rate"], label="AB")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rate")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "ev_release_rates.png", dpi=150)
    plt.close(fig)



def plot_cumulative_ev_yield(result: SimulationResult, outdir: Path) -> None:
    """Plot cumulative EV yields."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["sEV_cumulative"], label="sEV")
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["mlEV_cumulative"], label="mlEV")
    ax.plot(result.ev_timeseries["t"], result.ev_timeseries["AB_cumulative"], label="AB")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Cumulative output")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "cumulative_ev_yield.png", dpi=150)
    plt.close(fig)



def plot_quality_viability(result: SimulationResult, outdir: Path) -> None:
    """Plot damage and viability proxy."""
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    ax.plot(result.state_timeseries["t"], result.state_timeseries["damage"], label="damage")
    ax.axhline(result.summary["viability_fraction"], color="tab:green", linestyle="--", label="viability fraction")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalized value")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "quality_viability.png", dpi=150)
    plt.close(fig)



def generate_all_plots(result: SimulationResult, outdir: Path) -> None:
    """Generate all standard plots."""
    plot_calcium_timeseries(result, outdir)
    plot_ev_release_rates(result, outdir)
    plot_cumulative_ev_yield(result, outdir)
    plot_quality_viability(result, outdir)
