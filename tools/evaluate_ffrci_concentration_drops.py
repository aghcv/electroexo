from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tools.analyze_ffrci_data import (  # noqa: E402
    EXOID_FILE,
    parse_exoid_file,
    parse_sample_label,
)


DEFAULT_DATA_DIR = REPO_ROOT / "data" / "experimental" / "ffrci_data_sharing"
DEFAULT_OUT_DIR = REPO_ROOT / "results" / "ffrci_concentration_drop_evaluation"

METRICS = {
    "total": "All measured sizes",
    "below_100_nm": "<100 nm",
    "100_to_below_200_nm": "100–<200 nm",
    "at_or_above_200_nm": "≥200 nm",
}
CONDITIONS = ["1p20kV", "3p40kV", "5p40kV"]
TIMES_H = [0.5, 1.0, 3.0]
INTERVALS = [(0.5, 1.0), (1.0, 3.0), (0.5, 3.0)]
SURFACE_BIN_EDGES_NM = np.arange(80.0, 381.0, 20.0)
SURFACE_BIN_CENTERS_NM = (SURFACE_BIN_EDGES_NM[:-1] + SURFACE_BIN_EDGES_NM[1:]) / 2.0


def _metric_values(group: pd.DataFrame) -> dict[str, float]:
    diameter = group["particle_diameter_nm"].to_numpy(dtype=float)
    concentration = group["concentration_particles_per_ml"].to_numpy(dtype=float)
    return {
        "total": float(concentration.sum()),
        "below_100_nm": float(concentration[diameter < 100].sum()),
        "100_to_below_200_nm": float(
            concentration[(diameter >= 100) & (diameter < 200)].sum()
        ),
        "at_or_above_200_nm": float(concentration[diameter >= 200].sum()),
    }


def build_sample_metrics(raw: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_number, group in raw.groupby("dataset_number", sort=True):
        label = str(group["label"].iloc[0])
        metrics = _metric_values(group)
        total = metrics["total"]
        rows.append(
            {
                "dataset_number": int(dataset_number),
                "label": label,
                **parse_sample_label(label),
                **{f"{key}_particles_per_ml": value for key, value in metrics.items()},
                "fraction_below_100_nm": metrics["below_100_nm"] / total,
                "fraction_100_to_below_200_nm": metrics["100_to_below_200_nm"] / total,
                "fraction_at_or_above_200_nm": metrics["at_or_above_200_nm"] / total,
            }
        )
    return pd.DataFrame(rows)


def to_long(sample_metrics: pd.DataFrame) -> pd.DataFrame:
    id_columns = [
        "dataset_number",
        "label",
        "sample_type",
        "condition",
        "control_label",
        "pulse_count",
        "amplitude_label_kV",
        "harvest_time_h",
        "replicate",
    ]
    value_columns = [f"{metric}_particles_per_ml" for metric in METRICS]
    long = sample_metrics.melt(
        id_vars=id_columns,
        value_vars=value_columns,
        var_name="metric",
        value_name="concentration_particles_per_ml",
    )
    long["metric"] = long["metric"].str.removesuffix("_particles_per_ml")
    long["size_range"] = long["metric"].map(METRICS)
    return long


def summarize_timecourses(long: pd.DataFrame) -> pd.DataFrame:
    treated = long[long["sample_type"] == "treatment"]
    summary = (
        treated.groupby(["condition", "harvest_time_h", "metric", "size_range"], sort=True)[
            "concentration_particles_per_ml"
        ]
        .agg(n_nominal_p_labels="count", mean="mean", sample_sd="std", minimum="min", maximum="max")
        .reset_index()
    )
    summary["cv"] = summary["sample_sd"] / summary["mean"]
    summary["relative_range"] = (summary["maximum"] - summary["minimum"]) / summary["mean"]
    return summary.rename(
        columns={
            "mean": "mean_particles_per_ml",
            "sample_sd": "sample_sd_particles_per_ml",
            "minimum": "minimum_particles_per_ml",
            "maximum": "maximum_particles_per_ml",
        }
    )


def summarize_controls(long: pd.DataFrame) -> pd.DataFrame:
    controls = long[long["sample_type"] == "control"]
    summary = (
        controls.groupby(["condition", "metric", "size_range"], sort=True)[
            "concentration_particles_per_ml"
        ]
        .agg(n_nominal_p_labels="count", mean="mean", sample_sd="std", minimum="min", maximum="max")
        .reset_index()
    )
    summary["cv"] = summary["sample_sd"] / summary["mean"]
    summary["relative_range"] = (summary["maximum"] - summary["minimum"]) / summary["mean"]
    summary.insert(1, "harvest_time_available", False)
    return summary.rename(
        columns={
            "mean": "mean_particles_per_ml",
            "sample_sd": "sample_sd_particles_per_ml",
            "minimum": "minimum_particles_per_ml",
            "maximum": "maximum_particles_per_ml",
        }
    )


def _direction(value: float, tolerance: float = 1e-12) -> str:
    if value < -tolerance:
        return "decrease"
    if value > tolerance:
        return "increase"
    return "unchanged"


def interval_changes(long: pd.DataFrame) -> pd.DataFrame:
    treated = long[long["sample_type"] == "treatment"]
    rows: list[dict[str, object]] = []
    for condition in CONDITIONS:
        condition_data = treated[treated["condition"] == condition]
        for metric, size_range in METRICS.items():
            metric_data = condition_data[condition_data["metric"] == metric]
            pivot = metric_data.pivot(index="replicate", columns="harvest_time_h", values="concentration_particles_per_ml")
            means = metric_data.groupby("harvest_time_h")["concentration_particles_per_ml"].mean()
            for start_h, end_h in INTERVALS:
                paired = pivot[[start_h, end_h]].dropna()
                paired_percent = 100.0 * (paired[end_h] - paired[start_h]) / paired[start_h]
                start = float(means.loc[start_h])
                end = float(means.loc[end_h])
                change = end - start
                percent = 100.0 * change / start
                decreases = int((paired_percent < 0).sum())
                increases = int((paired_percent > 0).sum())
                unchanged = int((paired_percent == 0).sum())
                if len(paired) and decreases == len(paired):
                    consistency = "all paired p labels decrease"
                elif len(paired) and increases == len(paired):
                    consistency = "all paired p labels increase"
                else:
                    consistency = "mixed directions across paired p labels"
                rows.append(
                    {
                        "condition": condition,
                        "metric": metric,
                        "size_range": size_range,
                        "start_time_h": start_h,
                        "end_time_h": end_h,
                        "start_mean_particles_per_ml": start,
                        "end_mean_particles_per_ml": end,
                        "absolute_change_particles_per_ml": change,
                        "percent_change_of_group_mean": percent,
                        "direction_of_group_mean": _direction(percent),
                        "n_paired_nominal_p_labels": int(len(paired)),
                        "n_paired_decreasing": decreases,
                        "n_paired_increasing": increases,
                        "n_paired_unchanged": unchanged,
                        "paired_direction_consistency": consistency,
                        "median_paired_percent_change": float(paired_percent.median()),
                        "minimum_paired_percent_change": float(paired_percent.min()),
                        "maximum_paired_percent_change": float(paired_percent.max()),
                    }
                )
    return pd.DataFrame(rows)


def _time_label(time_h: float) -> str:
    return "30 min" if time_h == 0.5 else f"{time_h:g} h"


def plot_timecourses(summary: pd.DataFrame, output: Path) -> None:
    colors = {"1p20kV": "#2878B5", "3p40kV": "#E07A1F", "5p40kV": "#B23A48"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
    for ax, (metric, title) in zip(axes.flat, METRICS.items(), strict=True):
        for condition in CONDITIONS:
            subset = summary[(summary["condition"] == condition) & (summary["metric"] == metric)].sort_values(
                "harvest_time_h"
            )
            x = subset["harvest_time_h"].to_numpy(dtype=float)
            y = subset["mean_particles_per_ml"].to_numpy(dtype=float)
            sd = subset["sample_sd_particles_per_ml"].to_numpy(dtype=float)
            ax.errorbar(
                x,
                y / 1e9,
                yerr=sd / 1e9,
                marker="o",
                linewidth=2,
                capsize=3,
                color=colors[condition],
                label=condition,
            )
        ax.set_title(title)
        ax.set_xticks(TIMES_H, [_time_label(value) for value in TIMES_H])
        ax.set_ylabel(r"Concentration ($10^9$ particles/mL)")
        ax.grid(alpha=0.25)
    axes[0, 0].legend(frameon=False, ncol=3, fontsize=9)
    fig.suptitle(
        "Stimulated CD4+ extracellular-particle trajectories\n"
        "Error bars are spread across nominal p labels, not confirmed biological-replicate uncertainty",
        fontsize=14,
    )
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_interval_heatmap(changes: pd.DataFrame, output: Path) -> None:
    interval_order = ["30 min→1 h", "1 h→3 h", "30 min→3 h"]
    interval_map = {
        (0.5, 1.0): interval_order[0],
        (1.0, 3.0): interval_order[1],
        (0.5, 3.0): interval_order[2],
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    for ax, (metric, title) in zip(axes.flat, METRICS.items(), strict=True):
        subset = changes[changes["metric"] == metric].copy()
        subset["interval"] = [interval_map[(a, b)] for a, b in zip(subset["start_time_h"], subset["end_time_h"])]
        table = subset.pivot(index="condition", columns="interval", values="percent_change_of_group_mean").loc[
            CONDITIONS, interval_order
        ]
        image = ax.imshow(table.to_numpy(), cmap="RdBu_r", vmin=-80, vmax=80, aspect="auto")
        for row in range(table.shape[0]):
            for column in range(table.shape[1]):
                value = table.iloc[row, column]
                ax.text(column, row, f"{value:+.1f}%", ha="center", va="center", fontsize=10)
        ax.set_title(title)
        ax.set_xticks(range(len(interval_order)), interval_order, rotation=15, ha="right")
        ax.set_yticks(range(len(CONDITIONS)), CONDITIONS)
    colorbar = fig.colorbar(image, ax=axes, shrink=0.78, label="Change in group mean (%)")
    colorbar.ax.axhline(0, color="black", linewidth=0.5)
    fig.suptitle("Direction and magnitude of stimulated-sample concentration changes", fontsize=14)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_controls(summary: pd.DataFrame, output: Path) -> None:
    order = ["sham2", "sham_media"]
    display = {"sham2": "Sham2", "sham_media": "Sham media"}
    colors = ["#5B8E7D", "#8A7090"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), constrained_layout=True)

    total = summary[summary["metric"] == "total"].set_index("condition").loc[order]
    axes[0].bar(
        [display[value] for value in order],
        total["mean_particles_per_ml"] / 1e9,
        yerr=total["sample_sd_particles_per_ml"] / 1e9,
        capsize=4,
        color=colors,
    )
    axes[0].set_title("All measured sizes")
    axes[0].set_ylabel(r"Concentration ($10^9$ particles/mL)")
    axes[0].grid(axis="y", alpha=0.25)

    bottom = np.zeros(len(order))
    band_colors = {"below_100_nm": "#4C78A8", "100_to_below_200_nm": "#F2CF5B", "at_or_above_200_nm": "#E45756"}
    for metric in ["below_100_nm", "100_to_below_200_nm", "at_or_above_200_nm"]:
        subset = summary[summary["metric"] == metric].set_index("condition").loc[order]
        values = subset["mean_particles_per_ml"].to_numpy(dtype=float) / 1e9
        axes[1].bar(
            [display[value] for value in order],
            values,
            bottom=bottom,
            color=band_colors[metric],
            label=METRICS[metric],
        )
        bottom += values
    axes[1].set_title("Broad size-band composition")
    axes[1].set_ylabel(r"Concentration ($10^9$ particles/mL)")
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].grid(axis="y", alpha=0.25)
    fig.suptitle("Control samples are not a time course: the file contains no sham harvest times", fontsize=14)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def build_sham_normalized_surface_values(
    raw: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata_rows: list[dict[str, object]] = []
    for dataset_number, group in raw.groupby("dataset_number", sort=True):
        metadata_rows.append(
            {
                "dataset_number": int(dataset_number),
                **parse_sample_label(str(group["label"].iloc[0])),
            }
        )
    metadata = pd.DataFrame(metadata_rows)

    binned = raw.copy()
    binned["size_bin_center_nm"] = pd.cut(
        binned["particle_diameter_nm"],
        bins=SURFACE_BIN_EDGES_NM,
        labels=SURFACE_BIN_CENTERS_NM,
        right=False,
    ).astype(float)
    binned = (
        binned.dropna(subset=["size_bin_center_nm"])
        .groupby(["dataset_number", "size_bin_center_nm"], observed=True)[
            "concentration_particles_per_ml"
        ]
        .sum()
        .reset_index()
        .merge(metadata, on="dataset_number", validate="many_to_one")
    )
    binned["size_bin_lower_nm"] = binned["size_bin_center_nm"] - 10.0
    binned["size_bin_upper_nm"] = binned["size_bin_center_nm"] + 10.0

    treated = binned[binned["sample_type"] == "treatment"].copy()
    sham2 = binned[binned["condition"] == "sham2"][
        [
            "replicate",
            "size_bin_center_nm",
            "concentration_particles_per_ml",
        ]
    ].rename(columns={"concentration_particles_per_ml": "matched_sham2_particles_per_ml"})
    matched = treated.merge(
        sham2,
        on=["replicate", "size_bin_center_nm"],
        how="left",
        validate="many_to_one",
    )
    matched = matched.rename(
        columns={"concentration_particles_per_ml": "treatment_particles_per_ml"}
    )
    if matched["matched_sham2_particles_per_ml"].isna().any():
        raise ValueError("A treated size bin is missing its p-matched Sham2 reference")
    if (matched["matched_sham2_particles_per_ml"] <= 0).any():
        raise ValueError("A p-matched Sham2 surface denominator is not positive")
    matched["fold_of_matched_sham2"] = (
        matched["treatment_particles_per_ml"]
        / matched["matched_sham2_particles_per_ml"]
    )
    matched = matched[
        [
            "condition",
            "harvest_time_h",
            "replicate",
            "size_bin_lower_nm",
            "size_bin_upper_nm",
            "size_bin_center_nm",
            "treatment_particles_per_ml",
            "matched_sham2_particles_per_ml",
            "fold_of_matched_sham2",
        ]
    ].sort_values(["replicate", "condition", "harvest_time_h", "size_bin_center_nm"])

    treatment_average = (
        matched.groupby(
            [
                "condition",
                "harvest_time_h",
                "size_bin_lower_nm",
                "size_bin_upper_nm",
                "size_bin_center_nm",
            ],
            sort=True,
        )["treatment_particles_per_ml"]
        .agg(treatment_mean_particles_per_ml="mean", n_treatment_p_labels="count")
        .reset_index()
    )
    sham_average = (
        sham2.groupby("size_bin_center_nm", sort=True)["matched_sham2_particles_per_ml"]
        .agg(sham2_mean_particles_per_ml="mean", n_sham2_p_labels="count")
        .reset_index()
    )
    average = treatment_average.merge(
        sham_average,
        on="size_bin_center_nm",
        how="left",
        validate="many_to_one",
    )
    average["fold_of_mean_sham2"] = (
        average["treatment_mean_particles_per_ml"]
        / average["sham2_mean_particles_per_ml"]
    )
    return matched, average


def _complete_row_blocks(z_values: np.ndarray) -> list[np.ndarray]:
    complete = np.isfinite(z_values).all(axis=1)
    blocks: list[np.ndarray] = []
    start: int | None = None
    for index, valid in enumerate(complete):
        if valid and start is None:
            start = index
        if start is not None and (not valid or index == len(complete) - 1):
            stop = index if not valid else index + 1
            if stop - start >= 2:
                blocks.append(np.arange(start, stop))
            start = None
    return blocks


def plot_sham_normalized_surfaces(
    values: pd.DataFrame,
    *,
    value_column: str,
    title: str,
    output: Path,
) -> None:
    finite_values = values[value_column].replace([np.inf, -np.inf], np.nan).dropna()
    z_max = max(2.0, float(np.ceil(finite_values.max())))
    norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=z_max)
    figure = plt.figure(figsize=(16, 5.8))
    axes = [figure.add_subplot(1, 3, index + 1, projection="3d") for index in range(3)]

    for axis, condition in zip(axes, CONDITIONS, strict=True):
        subset = values[values["condition"] == condition]
        pivot = subset.pivot(
            index="harvest_time_h",
            columns="size_bin_center_nm",
            values=value_column,
        ).reindex(index=TIMES_H, columns=SURFACE_BIN_CENTERS_NM)
        z_values = pivot.to_numpy(dtype=float)
        x_values, y_values = np.meshgrid(SURFACE_BIN_CENTERS_NM, TIMES_H)
        for row_indices in _complete_row_blocks(z_values):
            axis.plot_surface(
                x_values[row_indices],
                y_values[row_indices],
                z_values[row_indices],
                cmap="coolwarm",
                norm=norm,
                linewidth=0.25,
                edgecolor=(0.15, 0.15, 0.15, 0.35),
                antialiased=True,
                alpha=0.92,
            )
        for row_index, time_h in enumerate(TIMES_H):
            finite = np.isfinite(z_values[row_index])
            if finite.any():
                axis.plot(
                    x_values[row_index, finite],
                    y_values[row_index, finite],
                    z_values[row_index, finite],
                    color="#202020",
                    linewidth=0.9,
                    alpha=0.8,
                )

        missing_times = [
            _time_label(TIMES_H[index])
            for index in range(len(TIMES_H))
            if not np.isfinite(z_values[index]).any()
        ]
        panel_title = condition
        if missing_times:
            panel_title += f"\n{', '.join(missing_times)} missing"
        axis.set_title(panel_title, pad=5, fontsize=12)
        axis.set_xlabel("Size-bin center (nm)", labelpad=4, fontsize=9)
        axis.set_ylabel("Time (h)", labelpad=4, fontsize=9)
        axis.set_zlabel("")
        axis.set_xlim(float(SURFACE_BIN_CENTERS_NM.min()), float(SURFACE_BIN_CENTERS_NM.max()))
        axis.set_ylim(0.5, 3.0)
        axis.set_zlim(0.0, z_max)
        axis.set_xticks([90, 150, 210, 270, 330, 370])
        axis.set_yticks(TIMES_H, ["0.5", "1", "3"])
        axis.tick_params(axis="both", which="major", labelsize=8, pad=1)
        axis.view_init(elev=28, azim=-124)
        axis.set_box_aspect((1.75, 1.05, 0.95), zoom=0.88)

    scalar = matplotlib.cm.ScalarMappable(norm=norm, cmap="coolwarm")
    scalar.set_array([])
    color_axis = figure.add_axes([0.925, 0.20, 0.015, 0.52])
    figure.colorbar(
        scalar,
        cax=color_axis,
        label="Treated / matched Sham2",
    )
    figure.suptitle(
        f"{title}\n20-nm bands over 80–380 nm; a ratio of 1 equals the Sham2 reference",
        fontsize=14,
        y=0.98,
    )
    figure.subplots_adjust(left=0.015, right=0.90, bottom=0.05, top=0.80, wspace=0.02)
    figure.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(figure)


def build_data_quality(
    sample_metrics: pd.DataFrame,
    timecourse_summary: pd.DataFrame,
) -> dict[str, object]:
    total_summary = timecourse_summary[timecourse_summary["metric"] == "total"]
    identical = total_summary["relative_range"] < 1e-8
    expected = pd.MultiIndex.from_product([CONDITIONS, TIMES_H, [1, 2, 3]])
    treated = sample_metrics[sample_metrics["sample_type"] == "treatment"]
    observed = pd.MultiIndex.from_frame(treated[["condition", "harvest_time_h", "replicate"]])
    missing = expected.difference(observed)
    return {
        "n_datasets": int(len(sample_metrics)),
        "n_treatment_datasets": int((sample_metrics["sample_type"] == "treatment").sum()),
        "n_control_datasets": int((sample_metrics["sample_type"] == "control").sum()),
        "control_harvest_times_present": False,
        "missing_expected_treatment_labels": [
            {"condition": str(condition), "harvest_time_h": float(time_h), "replicate": int(replicate)}
            for condition, time_h, replicate in missing
        ],
        "n_treatment_condition_time_cells": int(len(total_summary)),
        "n_cells_with_numerically_identical_total_across_p_labels": int(identical.sum()),
        "identical_total_relative_range_tolerance": 1e-8,
        "interpretive_warning": (
            "The p-labeled records have different size distributions but effectively identical summed total "
            "concentration within every condition/time cell. They must not be treated as independent biological "
            "replicates for inference about total concentration without provenance metadata."
        ),
    }


def write_report(
    path: Path,
    changes: pd.DataFrame,
    controls: pd.DataFrame,
    data_quality: dict[str, object],
) -> None:
    total = changes[changes["metric"] == "total"].copy()

    def pct(condition: str, start: float, end: float) -> float:
        return float(
            total[
                (total["condition"] == condition)
                & (total["start_time_h"] == start)
                & (total["end_time_h"] == end)
            ]["percent_change_of_group_mean"].iloc[0]
        )

    sham_total = controls[controls["metric"] == "total"].set_index("condition")
    sham_ratio = float(
        sham_total.loc["sham2", "mean_particles_per_ml"]
        / sham_total.loc["sham_media", "mean_particles_per_ml"]
    )
    report = f"""# FFRCI Exoid concentration-drop evaluation

## Bottom line

The exported stimulated trajectories do contain decreases, but not one common monotonic pattern:

- **1p20kV:** {_time_label(0.5)}→{_time_label(1.0)} {pct('1p20kV', 0.5, 1.0):+.1f}%, {_time_label(1.0)}→{_time_label(3.0)} {pct('1p20kV', 1.0, 3.0):+.1f}%, and net {_time_label(0.5)}→{_time_label(3.0)} {pct('1p20kV', 0.5, 3.0):+.1f}%.
- **3p40kV:** {_time_label(0.5)}→{_time_label(1.0)} {pct('3p40kV', 0.5, 1.0):+.1f}%, {_time_label(1.0)}→{_time_label(3.0)} {pct('3p40kV', 1.0, 3.0):+.1f}%, and net {_time_label(0.5)}→{_time_label(3.0)} {pct('3p40kV', 0.5, 3.0):+.1f}%.
- **5p40kV:** {_time_label(0.5)}→{_time_label(1.0)} {pct('5p40kV', 0.5, 1.0):+.1f}%, {_time_label(1.0)}→{_time_label(3.0)} {pct('5p40kV', 1.0, 3.0):+.1f}%, and net {_time_label(0.5)}→{_time_label(3.0)} {pct('5p40kV', 0.5, 3.0):+.1f}%.

Thus, **5p40kV is the only condition with a strong net fall from 30 min to 3 h**. The 1p20kV series has an early trough followed by rebound; the 3p40kV series has a 1 h peak followed by partial decline.

These are measured extracellular-particle concentrations, not direct measurements of secretion rate. A decrease can reflect reduced production, extracellular removal, sample handling, or measurement processing; it is not evidence of a negative release rate by itself. Exoid particle counts also do not establish that every detected particle is an EV.

No sham time trend can be evaluated from this file. `Sham2` and `sham media` each have p1–p3 records, but neither has a harvest-time label. They are distinct control labels, not ordered time points. Their total concentrations differ by {sham_ratio:.2f}-fold, which reinforces that they should not be interpreted as successive sham measurements without protocol metadata.

## Size-resolved evidence

- For **5p40kV**, the 1→3 h decline occurs in all three broad bands (`<100`, `100–<200`, and `≥200 nm`) and in all three paired p labels. This makes a simple shift across the 200 nm cutoff an insufficient explanation.
- For **3p40kV**, the 1→3 h total decline is driven mainly by the `100–<200 nm` band. The `<100 nm` direction is mixed across the two p labels available at 1 h.
- For **1p20kV**, all p labels show the early fall and later rebound in both `<100` and `100–<200 nm`; the `≥200 nm` direction is less uniform.

See `stimulated_interval_changes.csv` for all effect sizes and direction counts.

## Size by time surfaces normalized to matched sham

The surface figures use common 20-nm display bands from 80 to 380 nm. Each p-specific surface divides a treated concentration by the `Sham2` concentration with the same p label: treated p1 by Sham2 p1, treated p2 by Sham2 p2, and treated p3 by Sham2 p3. The average surface is the mean treated concentration divided by the mean Sham2 concentration in each size band; it is not the mean of the three ratios.

These ratios are descriptive because the file gives no harvest time for `Sham2`. The same static p-matched Sham2 size distribution is used at all treated time points. The `3p40kV 1 h p3` record is absent, so that row is deliberately not imputed in the p3 surface.

## What can and cannot be claimed

This is **descriptive evidence of trajectory shape**, not a significance test. Every stimulated condition/time cell has an effectively identical all-size total across its p labels ({data_quality['n_cells_with_numerically_identical_total_across_p_labels']}/{data_quality['n_treatment_condition_time_cells']} cells), despite different size distributions. That pattern suggests the totals may have been copied, normalized, or used to scale multiple distributions. Consequently, the p labels do not provide independent variance for testing total-concentration change. The size bins are subdivisions of the same observation and are also not replicates.

The conclusion also depends on the measurements being comparable endpoint concentrations. It changes if the time points represent interval collections, cumulative collections, serial sampling with volume removal, different wells, or medium replacement.

## Hypotheses to investigate

1. **Export or normalization behavior.** Confirm whether p1–p3 are biological cultures, technical repeat scans, merged distributions, or software-generated curves, and obtain the raw per-sample concentration totals. The within-cell total duplication is the first issue to resolve.
2. **Sampling and collection protocol.** Medium replacement, serial withdrawal, dilution, different collection intervals, or unequal recovered volume can create apparent falls. Reconstruct volume added/removed and whether each time point came from the same or separate culture.
3. **Loss of viable producer cells after the stronger exposure.** The broad, consistent late decline for 5p40kV is compatible with fewer viable/secreting cells or delayed injury. Time-matched viable-cell counts, viability/apoptosis, LDH, and total protein would test this.
4. **Extracellular clearance or loss.** Cellular re-uptake, particle degradation, adsorption to plastic, sedimentation, aggregation, or loss during storage/handling could reduce measured supernatant concentration. Cell-free spike/recovery controls and matched supernatant/cell-associated measurements would help separate these mechanisms.
5. **Transient release and adaptation.** The 3p40kV peak and the 1p20kV rebound could reflect pulse-triggered release followed by clearance or replenishment, but protocol and raw-replicate validation are required before treating these shapes as biological kinetics.
6. **Size redistribution or detection efficiency.** This may contribute to condition-specific band changes, especially 1p20kV and 3p40kV, but it does not by itself explain the 5p40kV decline across all broad bands. Instrument detection limits, concentration range, dilution, and QC flags should be reviewed.

## Scale and cell count

The analysis stays in the reported particles/mL units. If all samples truly use the same 5 million starting cells and the same effective medium/recovery volume, conversion to particles per initial cell multiplies every point by the same constant and does not change any percent change or direction reported here. If volume, viable-cell count, or recovery differs by time/condition, those values are required for a valid per-cell comparison.

## Files

- `sample_size_band_concentrations.csv`: one row per exported dataset.
- `stimulated_timecourse_summary.csv`: mean and spread across nominal p labels.
- `stimulated_interval_changes.csv`: effect sizes and paired direction counts.
- `sham_control_summary.csv`: standalone control summaries; explicitly not a time course.
- `data_quality_summary.json`: missing-label and total-duplication audit.
- `particle_size_time_surface_values_by_p.csv`: 20-nm concentrations and p-matched Sham2 ratios.
- `particle_size_time_surface_values_average.csv`: ratio of treatment and Sham2 means.
- `surface_relative_to_sham2_p1.png`, `p2.png`, and `p3.png`: separate p-matched surfaces.
- `surface_relative_to_sham2_average.png`: average surface across available p labels.

## Reproduce

From the repository root, run:

```bash
MPLCONFIGDIR=/tmp/electroexo-mpl python tools/evaluate_ffrci_concentration_drops.py
```
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate concentration decreases in the FFRCI CD4+ Exoid export.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw = parse_exoid_file(args.data_dir / EXOID_FILE)
    sample_metrics = build_sample_metrics(raw)
    long = to_long(sample_metrics)
    timecourse_summary = summarize_timecourses(long)
    control_summary = summarize_controls(long)
    changes = interval_changes(long)
    matched_surface, average_surface = build_sham_normalized_surface_values(raw)
    data_quality = build_data_quality(sample_metrics, timecourse_summary)

    sample_metrics.to_csv(args.out_dir / "sample_size_band_concentrations.csv", index=False)
    timecourse_summary.to_csv(args.out_dir / "stimulated_timecourse_summary.csv", index=False)
    changes.to_csv(args.out_dir / "stimulated_interval_changes.csv", index=False)
    control_summary.to_csv(args.out_dir / "sham_control_summary.csv", index=False)
    matched_surface.to_csv(
        args.out_dir / "particle_size_time_surface_values_by_p.csv", index=False
    )
    average_surface.to_csv(
        args.out_dir / "particle_size_time_surface_values_average.csv", index=False
    )
    (args.out_dir / "data_quality_summary.json").write_text(
        json.dumps(data_quality, indent=2) + "\n", encoding="utf-8"
    )
    write_report(args.out_dir / "README.md", changes, control_summary, data_quality)
    plot_timecourses(timecourse_summary, args.out_dir / "stimulated_concentration_timecourses.png")
    plot_interval_heatmap(changes, args.out_dir / "stimulated_interval_change_heatmap.png")
    plot_controls(control_summary, args.out_dir / "sham_control_context.png")
    for replicate in (1, 2, 3):
        plot_sham_normalized_surfaces(
            matched_surface[matched_surface["replicate"] == replicate],
            value_column="fold_of_matched_sham2",
            title=f"Stimulated CD4+ particle concentration relative to Sham2 p{replicate}",
            output=args.out_dir / f"surface_relative_to_sham2_p{replicate}.png",
        )
    plot_sham_normalized_surfaces(
        average_surface,
        value_column="fold_of_mean_sham2",
        title="Mean stimulated CD4+ particle concentration relative to mean Sham2",
        output=args.out_dir / "surface_relative_to_sham2_average.png",
    )


if __name__ == "__main__":
    main()
