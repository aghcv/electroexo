from __future__ import annotations

import argparse
import os
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "electroexo_matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from electro_exocytosis.config import CellStateConfig, ExposureConfig, PulseConfig
from electro_exocytosis.models.dosimetry import compute_dosimetry
from electro_exocytosis.models.electrodynamics import (
    ElectrodynamicsParams,
    compute_electrodynamics_state,
)
from electro_exocytosis.models.pulse import compute_pulse_descriptors
from electro_exocytosis.visualization.style import (
    MANUSCRIPT_LANDSCAPE_FIGSIZE,
    line_styles,
    save_manuscript_figure,
)

PULSE_WIDTHS_NS = (5.0, 10.0, 30.0, 100.0, 300.0, 1000.0)
FIELD_STRENGTHS_KV_CM = (0.5, 1.0, 3.0)
CELL_RADII_UM = (5.0, 10.0, 20.0)
MEMBRANE_CHARGING_TAU_S = 100e-9
PERMEABILIZATION_REFERENCE_V = 0.25


def build_response_table() -> pd.DataFrame:
    """Compute a compact grid of reduced membrane-electrodynamics responses."""
    rows: list[dict[str, float | str]] = []
    for amplitude_kV_cm in FIELD_STRENGTHS_KV_CM:
        for pulse_width_ns in PULSE_WIDTHS_NS:
            for radius_um in CELL_RADII_UM:
                rows.append(_response_row(amplitude_kV_cm, pulse_width_ns, radius_um))
    return pd.DataFrame(rows)


def _response_row(amplitude_kV_cm: float, pulse_width_ns: float, radius_um: float) -> dict[str, float | str]:
    pulse = PulseConfig(
        amplitude_kV_cm=amplitude_kV_cm,
        pulse_width_ns=pulse_width_ns,
        pulse_number=1,
        repetition_rate_Hz=1,
    )
    exposure = ExposureConfig()
    params = ElectrodynamicsParams(
        cell_radius_m=radius_um * 1e-6,
        membrane_charging_tau_s=MEMBRANE_CHARGING_TAU_S,
    )
    descriptors = compute_pulse_descriptors(pulse, exposure)
    dosimetry = compute_dosimetry(descriptors, exposure)
    state = compute_electrodynamics_state(descriptors, dosimetry, CellStateConfig(), params)
    return {
        "amplitude_kV_cm": amplitude_kV_cm,
        "pulse_width_ns": pulse_width_ns,
        "cell_radius_um": radius_um,
        "membrane_charging_tau_ns": MEMBRANE_CHARGING_TAU_S * 1e9,
        "schwan_limit_V": state.schwan_limit_V,
        "membrane_charging_factor": state.membrane_charging_factor,
        "delta_Vm_V": state.delta_Vm,
        "delta_V_ER_V": state.delta_V_ER,
        "delta_V_mito_V": state.delta_V_mito,
        "delta_V_MVB_V": state.delta_V_MVB,
        "membrane_permeability": state.membrane_permeability,
        "pore_density_m2": state.pore_density,
        "pore_density_um2": state.pore_density / 1.0e12,
    }


def plot_pulse_width_response(summary: pd.DataFrame, outdir: Path) -> None:
    """Plot charging and membrane voltage response versus pulse duration."""
    radius_um = 10.0
    subset = summary[summary["cell_radius_um"] == radius_um]
    styles = line_styles(len(FIELD_STRENGTHS_KV_CM))
    fig, axes = plt.subplots(1, 2, figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)

    charging = subset[subset["amplitude_kV_cm"] == FIELD_STRENGTHS_KV_CM[0]]
    axes[0].semilogx(
        charging["pulse_width_ns"],
        charging["membrane_charging_factor"],
        **styles[0],
    )
    axes[0].set_xlabel("Pulse width (ns)")
    axes[0].set_ylabel("Membrane charging factor")
    axes[0].set_ylim(0.0, 1.05)

    for index, amplitude_kV_cm in enumerate(FIELD_STRENGTHS_KV_CM):
        field_subset = subset[subset["amplitude_kV_cm"] == amplitude_kV_cm]
        axes[1].semilogx(
            field_subset["pulse_width_ns"],
            field_subset["delta_Vm_V"],
            label=f"{amplitude_kV_cm:g} kV/cm",
            **styles[index],
        )
    axes[1].axhline(PERMEABILIZATION_REFERENCE_V, color="0.35", linestyle=":", linewidth=1)
    axes[1].set_xlabel("Pulse width (ns)")
    axes[1].set_ylabel("Induced membrane voltage (V)")
    axes[1].legend(title="Peak field", fontsize=8)
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "membrane_voltage_by_pulse_width.png")
    plt.close(fig)


def plot_radius_sensitivity(summary: pd.DataFrame, outdir: Path) -> None:
    """Plot cell-radius sensitivity at a fixed threshold-scale field."""
    amplitude_kV_cm = 1.0
    pulse_widths_ns = (10.0, 100.0, 300.0)
    subset = summary[summary["amplitude_kV_cm"] == amplitude_kV_cm]
    styles = line_styles(len(pulse_widths_ns))
    fig, ax = plt.subplots(figsize=MANUSCRIPT_LANDSCAPE_FIGSIZE)

    for index, pulse_width_ns in enumerate(pulse_widths_ns):
        width_subset = subset[subset["pulse_width_ns"] == pulse_width_ns]
        ax.plot(
            width_subset["cell_radius_um"],
            width_subset["delta_Vm_V"],
            label=f"{pulse_width_ns:g} ns",
            **styles[index],
        )
    ax.axhline(PERMEABILIZATION_REFERENCE_V, color="0.35", linestyle=":", linewidth=1)
    ax.set_xlabel("Cell radius (um)")
    ax.set_ylabel("Induced membrane voltage (V)")
    ax.legend(title="Pulse width", fontsize=8)
    fig.tight_layout()
    save_manuscript_figure(fig, outdir / "cell_radius_sensitivity.png")
    plt.close(fig)


def write_outputs(summary: pd.DataFrame, outdir: Path, make_plots: bool) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "electrodynamics_response_summary.csv", index=False)
    if make_plots:
        plot_pulse_width_response(summary, outdir)
        plot_radius_sensitivity(summary, outdir)


def print_console_summary(summary: pd.DataFrame) -> None:
    columns = [
        "amplitude_kV_cm",
        "pulse_width_ns",
        "cell_radius_um",
        "membrane_charging_factor",
        "delta_Vm_V",
        "membrane_permeability",
    ]
    compact = summary[
        (summary["cell_radius_um"] == 10.0)
        & (summary["amplitude_kV_cm"].isin((0.5, 1.0)))
        & (summary["pulse_width_ns"].isin((10.0, 100.0, 300.0)))
    ]
    print(compact[columns].to_string(index=False, float_format=lambda value: f"{value:0.3f}"))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare reduced membrane/organelle electrodynamics responses across nsPEF scenarios."
    )
    parser.add_argument(
        "--out",
        default="results/electrodynamics_model_comparison",
        help="Directory for CSV and PNG outputs.",
    )
    parser.add_argument("--no-plots", action="store_true", help="Write CSV outputs only.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    summary = build_response_table()
    outdir = Path(args.out)
    write_outputs(summary, outdir, make_plots=not args.no_plots)
    print_console_summary(summary)
    print(f"\nWrote electrodynamics comparison outputs to {outdir}")


if __name__ == "__main__":
    main()
