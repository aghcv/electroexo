from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from electro_exocytosis.config import (
    CellStateConfig,
    ExposureConfig,
    PulseConfig,
    ScenarioConfig,
    SimulationConfig,
    SimulationScenario,
)
from electro_exocytosis.simulation import Simulation
from electro_exocytosis.visualization.style import MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE, save_manuscript_figure


DEFAULT_AMPLITUDES = (6.0, 10.0, 14.0, 18.0, 22.0, 26.0)
DEFAULT_PULSE_NUMBERS = (5, 10, 20, 30, 40, 60)


def build_regime_table(
    amplitudes_kV_cm: Sequence[float] = DEFAULT_AMPLITUDES,
    pulse_numbers: Sequence[int] = DEFAULT_PULSE_NUMBERS,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for amplitude_kV_cm in amplitudes_kV_cm:
        for pulse_number in pulse_numbers:
            scenario = SimulationScenario(
                scenario=ScenarioConfig(name=f"a{amplitude_kV_cm:g}_n{pulse_number}"),
                pulse=PulseConfig(
                    amplitude_kV_cm=amplitude_kV_cm,
                    pulse_width_ns=200.0,
                    pulse_number=pulse_number,
                    repetition_rate_Hz=5.0,
                ),
                exposure=ExposureConfig(dosimetry_model="joule_lumped_thermal"),
                cell_state=CellStateConfig(),
                simulation=SimulationConfig(t_start_s=0.0, t_end_s=1800.0, output_dt_s=10.0),
            )
            result = Simulation(scenario).run()
            sev = float(result.summary["cumulative_small_EV"])
            mlev = float(result.summary["cumulative_medium_large_EV"])
            ab = float(result.summary["cumulative_apoptotic_body"])
            dominant = max(
                [("sEV", sev), ("m/lEV", mlev), ("AB", ab)],
                key=lambda item: item[1],
            )[0]
            rows.append(
                {
                    "amplitude_kV_cm": amplitude_kV_cm,
                    "pulse_number": pulse_number,
                    "dose_index": result.summary["dose_index"],
                    "cumulative_small_EV": sev,
                    "cumulative_medium_large_EV": mlev,
                    "cumulative_apoptotic_body": ab,
                    "viability_fraction": float(result.summary["viability_fraction"]),
                    "dominant_subtype": dominant,
                }
            )
    return pd.DataFrame(rows)


def write_outputs(summary: pd.DataFrame, outdir: Path, make_plots: bool = True) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "ev_release_regime_map.csv", index=False)
    if make_plots:
        plot_regime_map(summary, outdir / "ev_release_regime_map.png")


def plot_regime_map(summary: pd.DataFrame, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    amplitudes = sorted(summary["amplitude_kV_cm"].unique())
    pulse_numbers = sorted(summary["pulse_number"].unique())
    fig, axes = plt.subplots(2, 2, figsize=MANUSCRIPT_PANEL_LANDSCAPE_FIGSIZE)

    metrics = [
        ("cumulative_small_EV", "sEV yield"),
        ("cumulative_medium_large_EV", "m/lEV yield"),
        ("cumulative_apoptotic_body", "AB yield"),
        ("viability_fraction", "Viability"),
    ]
    for ax, (metric, title) in zip(axes.flat, metrics, strict=True):
        grid = _pivot_grid(summary, metric, amplitudes, pulse_numbers)
        im = ax.imshow(grid, origin="lower", aspect="auto", cmap="viridis")
        ax.set_title(title, fontsize=10)
        ax.set_xticks(range(len(pulse_numbers)))
        ax.set_xticklabels([str(v) for v in pulse_numbers], fontsize=7)
        ax.set_yticks(range(len(amplitudes)))
        ax.set_yticklabels([f"{v:g}" for v in amplitudes], fontsize=7)
        ax.set_xlabel("Pulse number", fontsize=8)
        ax.set_ylabel("Amplitude (kV/cm)", fontsize=8)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()
    save_manuscript_figure(fig, output_path)
    plt.close(fig)


def _pivot_grid(summary: pd.DataFrame, metric: str, amplitudes: Sequence[float], pulse_numbers: Sequence[int]) -> np.ndarray:
    table = summary.pivot(index="amplitude_kV_cm", columns="pulse_number", values=metric)
    table = table.reindex(index=amplitudes, columns=pulse_numbers)
    return table.to_numpy(dtype=float)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/ev_release_regime_map"))
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG plot generation.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary = build_regime_table()
    write_outputs(summary, args.out, make_plots=not args.no_plots)
    print(f"Wrote EV release regime map outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
