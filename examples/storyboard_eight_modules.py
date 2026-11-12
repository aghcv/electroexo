"""Five controlled, eight-module storyboards for a presentation.

These figures deliberately use a presentation-size canvas (17.78 by 6.4 in),
not the repository's 7.1-in double-column manuscript layout.  Each panel shows
one quantity.  Axes are fixed across all studies so unchanged upstream panels
remain unchanged.  Module 6--8 quantities are endpoint proxies, not ODE states.

Run from the repository root::

    python -m examples.storyboard_eight_modules --out results/eight_module_storyboards

The manifest contains decimated plotting data for editable slide charts.  The
CSV files retain the complete 0.25-s output grid used to resolve the acute phase.
"""

from __future__ import annotations

import argparse
import json
import hashlib
import math
import platform
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

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
from electro_exocytosis.io.readers import merge_parameters
from electro_exocytosis.simulation import Simulation, SimulationResult
from electro_exocytosis.visualization.style import (
    COLORBLIND_PALETTE,
    manuscript_style_context,
    save_manuscript_figure,
)
from examples.storyboard_multilayer_overview import SCENARIOS as REGIME_SCENARIOS


PRESENTATION_FIGSIZE = (17.78, 6.4)
END_TIME_S = 1800.0
OUTPUT_DT_S = 0.25
SERIES_COLOR = COLORBLIND_PALETTE[0]

COMMON_PULSE = {
    "amplitude_kV_cm": 15.0,
    "pulse_width_ns": 200.0,
    "pulse_number": 20,
    "repetition_rate_Hz": 5.0,
    "waveform": "square",
}

COMMON_LIMITATIONS = (
    "These are illustrative simulations with provisional, uncalibrated parameters, "
    "not validated treatment recommendations or predictions for a named cell type.",
    "Deposited energy is an exposure descriptor. Heating is not coupled back into "
    "the injury equations, so the model cannot establish heat-caused injury.",
    "The reduced model uses a single-pulse charging descriptor and an effective "
    "perturbation duration; it does not resolve every nanosecond pulse or bipolar cancellation.",
    "The remodeling/repair score is algebraic and does not feed back to close pores "
    "or repair the damage state in this implementation.",
    "Release and recovered yield are relative model units, not calibrated particle counts. "
    "The protein cargo score is a terminal enrichment proxy, not measured cargo or function.",
    "The viability score is derived from peak modeled damage over the simulation; "
    "it is not a measured survival percentage or the recovered end-state.",
    "Recovered small-EV yield is the Layer 8 transform of cumulative release and "
    "assumed isolation efficiency. It is not extracellular medium stock or a sampling model.",
)

PANEL_DESCRIPTION = (
    "The eight columns show, respectively: deposited energy (J/mL); single-pulse "
    "plasma-membrane voltage (V); cytosolic calcium (µM, 0–60 s); repair activity "
    "(0–1 score, 0–60 s); cumulative small-EV release (relative model units, 0–30 min); "
    "protein cargo enrichment (terminal model score at 30 min); viability score "
    "(%, derived from peak modeled damage); and recovered small-EV yield "
    "(relative model units at 30 min). Energy, voltage, cargo, viability and recovery "
    "are scalar bars; calcium, repair and release are time courses. Axes match across "
    "all five figures. Input labels identify dependencies of the selected plotted readout, "
    "rather than every output of the module. The displayed recovered small-EV yield depends "
    "on cumulative release and assumed recovery; other process objectives also use cargo "
    "and quality. All runs use a 30-min horizon and 0.25-s output interval. "
)


@dataclass(frozen=True)
class StoryRow:
    id: str
    label: str
    detail: str
    pulse_overrides: dict[str, Any] = field(default_factory=dict)
    exposure_overrides: dict[str, Any] = field(default_factory=dict)
    cell_state_overrides: dict[str, Any] = field(default_factory=dict)
    parameter_overrides: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class StoryStudy:
    id: str
    title: str
    question: str
    takeaway: str
    caption: str
    script: str
    rows: tuple[StoryRow, ...]


@dataclass
class StoryboardDataset:
    manifest: dict[str, Any]
    summary: pd.DataFrame
    outputs: dict[str, SimulationResult]


MODULES = (
    {
        "number": 1,
        "title": "Pulse delivery, exposure geometry and dosimetry",
        "short_title": "Exposure and\ndosimetry",
        "metric_label": "Deposited energy",
        "unit": "J/mL",
        "kind": "bar",
        "x_label": "Exposure",
        "x_min": 0.0,
        "x_max": 1.0,
        "y_min": 0.0,
        "inputs": ["settings"],
    },
    {
        "number": 2,
        "title": "Plasma membrane and organelle electrodynamics",
        "short_title": "Membrane\nelectrodynamics",
        "metric_label": "Membrane voltage",
        "unit": "V",
        "kind": "bar",
        "x_label": "Per pulse",
        "x_min": 0.0,
        "x_max": 1.0,
        "y_min": 0.0,
        "inputs": [1],
    },
    {
        "number": 3,
        "title": "Ion transport, Ca²⁺ mobilization, ROS and bioenergetics",
        "short_title": "Calcium and\ncell energy",
        "metric_label": "Calcium",
        "unit": "µM",
        "kind": "line",
        "x_label": "Time (s)",
        "x_min": 0.0,
        "x_max": 60.0,
        "y_min": 0.0,
        "inputs": [2],
    },
    {
        "number": 4,
        "title": "Ca²⁺-dependent remodeling and repair",
        "short_title": "Remodeling\nand repair",
        "metric_label": "Repair activity",
        "unit": "0–1 score",
        "kind": "line",
        "x_label": "Time (s)",
        "x_min": 0.0,
        "x_max": 60.0,
        "y_min": 0.0,
        "y_max": 1.0,
        "inputs": [3],
    },
    {
        "number": 5,
        "title": "EV biogenesis and subtype release",
        "short_title": "EV formation\nand release",
        "metric_label": "Small-EV release",
        "unit": "relative units",
        "kind": "line",
        "x_label": "Time (min)",
        "x_min": 0.0,
        "x_max": 30.0,
        "y_min": 0.0,
        "inputs": [2, 3, 4, 7],
    },
    {
        "number": 6,
        "title": "Cargo sorting, EV composition and potency",
        "short_title": "Cargo and\npotency",
        "metric_label": "Protein cargo",
        "unit": "model score",
        "kind": "bar",
        "x_label": "At 30 min",
        "x_min": 0.0,
        "x_max": 1.0,
        "y_min": 0.0,
        "inputs": [3, 5],
    },
    {
        "number": 7,
        "title": "Injury, debris and quality gate",
        "short_title": "Injury and\nquality",
        "metric_label": "Viability score",
        "unit": "%",
        "kind": "bar",
        "x_label": "Peak injury",
        "x_min": 0.0,
        "x_max": 1.0,
        "y_min": 0.0,
        "y_max": 100.0,
        "inputs": [3],
    },
    {
        "number": 8,
        "title": "Manufacturing, isolation and QC interface",
        "short_title": "Isolation and\nquality control",
        "metric_label": "Recovered small EVs",
        "unit": "relative units",
        "kind": "bar",
        "x_label": "At 30 min",
        "x_min": 0.0,
        "x_max": 1.0,
        "y_min": 0.0,
        "inputs": [5, "recovery"],
    },
)


def _regime_rows() -> tuple[StoryRow, ...]:
    labels = ("Lower exposure", "Secretory example", "Injury example")
    narratives = (
        "Lower-exposure example with a shorter calcium response.",
        "Coordinated calcium, repair and secretory-pathway assumptions favor small-EV release.",
        "Higher modeled stress and apoptotic-body release.",
    )
    return tuple(
        StoryRow(
            id=spec.name,
            label=label,
            detail=(
                f"{spec.pulse_overrides['amplitude_kV_cm']:g} kV/cm; "
                f"{spec.pulse_overrides['pulse_width_ns']:g} ns; "
                f"{spec.pulse_overrides['pulse_number']} pulses; "
                f"{spec.pulse_overrides['repetition_rate_Hz']:g} Hz. "
                f"{narrative} This reproduces the original illustrative scenario, "
                "including its cell-state and pathway overrides; these rows are not a "
                "controlled dose-response comparison."
            ),
            pulse_overrides=spec.pulse_overrides,
            exposure_overrides=spec.exposure_overrides,
            cell_state_overrides=spec.cell_state_overrides,
            parameter_overrides=spec.parameter_overrides,
        )
        for spec, label, narrative in zip(REGIME_SCENARIOS, labels, narratives, strict=True)
    )


STUDIES = (
    StoryStudy(
        id="01_regime_examples",
        title="Follow the response from exposure to recovered vesicles",
        question="How can the same framework describe different response patterns?",
        takeaway="Different cell responses can lead to different vesicle outcomes.",
        caption=(
            "Original illustrative lower-exposure, secretory and injury-dominant "
            "scenarios, extended from five to eight modules. Several inputs differ "
            "between rows; this figure demonstrates the framework, not a controlled "
            "dose-response experiment."
        ),
        script=(
            "Read each row from left to right. We start with the electrical exposure, "
            "then follow membrane charging, the calcium response and repair activity. "
            "The last four columns show small-vesicle release, a simple protein cargo "
            "score, a cell viability score and what an assumed isolation step would "
            "recover. These are the same three illustrative examples as our earlier "
            "figure, with the remaining modules now visible. The secretory example "
            "produces more small vesicles in this model. The injury example shows why "
            "release must be considered alongside cell condition. Several settings "
            "change between these examples, so we use the next slides to vary one "
            "aspect at a time. The numerical scores are provisional model outputs."
        ),
        rows=_regime_rows(),
    ),
    StoryStudy(
        id="02_pulse_number",
        title="More pulses can change the response after membrane charging",
        question="What changes when we increase pulse number alone?",
        takeaway="A larger pulse train can prolong calcium signaling and increase injury.",
        caption=(
            "Controlled comparison of 5, 20 and 80 square pulses, each at 15 kV/cm, "
            "200 ns and 5 Hz. Only pulse number changes. Single-pulse membrane "
            "voltage is therefore identical. The effective perturbation duration "
            "is a reduced pulse-train approximation, not explicit pulse-by-pulse kinetics."
        ),
        script=(
            "Here we change only the number of pulses. Each pulse has the same field, "
            "width and spacing. Deposited energy rises, while the voltage reached "
            "during one pulse stays the same. The model then shows a larger and longer "
            "calcium response as the train gets larger. That does not translate into "
            "a simple increase in small-vesicle release: the viability score also "
            "falls. For an experiment, this comparison suggests measuring early "
            "calcium recovery and cell injury alongside vesicle output when choosing "
            "pulse number. Heating is shown only indirectly through the energy "
            "descriptor; the current equations do not establish that heating caused "
            "the injury. These patterns are hypotheses to test, not a calibrated dose recommendation."
        ),
        rows=tuple(
            StoryRow(
                id=f"pulses_{number:02d}",
                label=f"{number} pulses",
                detail=(
                    f"Only pulse number is changed to {number}; 15 kV/cm, 200 ns, "
                    "5 Hz square pulses, default cell and pathway settings."
                ),
                pulse_overrides={"pulse_number": number},
            )
            for number in (5, 20, 80)
        ),
    ),
    StoryStudy(
        id="03_calcium_clearance",
        title="Recovery speed can change the outcome of the same exposure",
        question="What if cells clear calcium at different rates?",
        takeaway="A prolonged calcium response can increase injury even with the same pulses.",
        caption=(
            "Controlled calcium-clearance comparison. The SERCA, PMCA and NCX maximum "
            "transport rates are scaled together by 0.25, 1 and 2, with all calcium influx "
            "settings and the 15 kV/cm, 200 ns, 20-pulse, 5-Hz exposure held fixed. "
            "These are hypothetical calcium-clearance states, not named or measured cell types."
        ),
        script=(
            "The first two columns are identical because the exposure is identical. "
            "This time we change only the capacities of three calcium-removal pathways. "
            "With lower clearance capacity, calcium remains elevated longer and the model "
            "assigns more injury. The repair activity also lasts longer, which should "
            "be read as a sustained demand for repair rather than evidence that the "
            "cell is healthier. A larger release signal can therefore accompany a "
            "less favorable cell state. This is useful for planning experiments: "
            "the decay of the calcium response may tell us more than its peak alone. "
            "The protein cargo score changes very little here because it uses the "
            "terminal cell state. The model has not identified these rows with actual cell types."
        ),
        rows=tuple(
            StoryRow(
                id=f"pump_{str(scale).replace('.', 'p')}",
                label=label,
                detail=(
                    f"Calcium-removal transport capacity {scale:g}× reference: "
                    f"SERCA={0.15 * scale:g}, PMCA={0.20 * scale:g} and "
                    f"NCX={0.06 * scale:g} µM/s. Calcium influx, injury thresholds, "
                    "all other parameters and exposure are unchanged."
                ),
                parameter_overrides={
                    "ion_transport": {
                        "SERCA_Vmax_uM_s": 0.15 * scale,
                        "PMCA_Vmax_uM_s": 0.20 * scale,
                        "NCX_Vmax_uM_s": 0.06 * scale,
                    }
                },
            )
            for scale, label in (
                (0.25, "Low clearance"),
                (1.0, "Reference"),
                (2.0, "High clearance"),
            )
        ),
    ),
    StoryStudy(
        id="04_ev_machinery",
        title="Vesicle machinery offers a different way to change small-EV output",
        question="Can vesicle pathways change output without increasing cell injury?",
        takeaway="Vesicle pathway settings change small-EV release and modeled protein cargo.",
        caption=(
            "Controlled EV-machinery comparison. The basal Rab-conversion, ESCRT and "
            "Rab-docking coefficients are scaled together by 0.25, 1 and 2. Basal "
            "small-EV release rate, all other parameters and exposure are held fixed. "
            "This is a coordinated pathway perturbation, not a validated drug or inhibitor model."
        ),
        script=(
            "In this comparison, exposure, calcium clearance and repair settings all "
            "stay the same. We change three parts of the vesicle machinery together: "
            "maturation, vesicle formation and docking. The calcium and repair traces "
            "therefore overlap. Farther to the right, the small-vesicle output and "
            "protein cargo score rise, while the viability score is almost unchanged. "
            "The point is that increasing electrical exposure is not the only possible "
            "way to alter output. The model can help formulate pathway experiments "
            "that ask whether vesicle production can be adjusted separately from "
            "cell stress. The cargo column is a simple model score; it does not "
            "establish a measured protein change or a functional benefit."
        ),
        rows=tuple(
            StoryRow(
                id=f"machinery_{str(scale).replace('.', 'p')}",
                label=label,
                detail=(
                    f"Rab-conversion, ESCRT and Rab-docking baselines {scale:g}× "
                    f"reference: {0.35 * scale:g}, {0.40 * scale:g} and "
                    f"{0.30 * scale:g}, respectively (dimensionless model coefficients). "
                    "The basal small-EV release rate is unchanged."
                ),
                parameter_overrides={
                    "ev_release": {
                        "rab_conversion_baseline": 0.35 * scale,
                        "escrt_baseline": 0.40 * scale,
                        "rab_docking_baseline": 0.30 * scale,
                    }
                },
            )
            for scale, label in (
                (0.25, "Lower pathway capacity"),
                (1.0, "Reference"),
                (2.0, "Higher pathway capacity"),
            )
        ),
    ),
    StoryStudy(
        id="05_isolation_recovery",
        title="Recovery can change the final yield without changing the cells",
        question="How much of the final difference could arise during isolation?",
        takeaway="The same biological release can give different recovered yields.",
        caption=(
            "Isolation-only negative control. Assumed isolation efficiency is 10%, "
            "30% or 60%, with all upstream inputs unchanged. Modules 1–7 therefore "
            "remain identical. The final column applies that efficiency to cumulative "
            "small-EV release; it is not extracellular medium stock, assay concentration "
            "or a comparison of named isolation methods."
        ),
        script=(
            "Every biological setting is identical in these three rows. The only "
            "change is the fraction assumed to be recovered during isolation. The "
            "first seven columns should therefore be identical, and they are. Only "
            "the final yield changes: a process that recovers sixty percent gives "
            "six times the small-vesicle yield of one that recovers ten percent. "
            "For an experimentalist, this is a reminder to measure process recovery "
            "before assigning an output difference to the cell. It also shows how "
            "the framework separates a biological hypothesis from a process assumption. "
            "These efficiencies are examples, not claims about particular isolation "
            "methods. This final column uses cumulative release rather than the "
            "separate medium-stock and sampling model."
        ),
        rows=tuple(
            StoryRow(
                id=f"recovery_{int(efficiency * 100):02d}",
                label=f"{int(efficiency * 100)}% recovery",
                detail=(
                    f"Only isolation_efficiency changes to {efficiency:g}; every "
                    "biological setting and the isolation method factor of 1 are "
                    "unchanged. The modeled harvest horizon is 0.5 h."
                ),
                parameter_overrides={
                    "manufacturing_qc": {"isolation_efficiency": efficiency}
                },
            )
            for efficiency in (0.10, 0.30, 0.60)
        ),
    ),
)


def build_scenario(study: StoryStudy, row: StoryRow) -> SimulationScenario:
    return SimulationScenario(
        scenario=ScenarioConfig(name=f"{study.id}__{row.id}"),
        pulse=PulseConfig(**(COMMON_PULSE | row.pulse_overrides)),
        exposure=ExposureConfig(
            **({"dosimetry_model": "joule_lumped_thermal"} | row.exposure_overrides)
        ),
        cell_state=CellStateConfig(**row.cell_state_overrides),
        simulation=SimulationConfig(t_end_s=END_TIME_S, output_dt_s=OUTPUT_DT_S),
    )


def _parameters(row: StoryRow) -> dict[str, Any]:
    # The simulation horizon and process harvest time must agree in every row.
    return merge_parameters(
        {"manufacturing_qc": {"harvest_time_h": END_TIME_S / 3600.0}},
        row.parameter_overrides,
    )


def _acute_plot_times() -> np.ndarray:
    return np.unique(
        np.concatenate(
            (
                np.arange(0.0, 10.0 + 0.125, 0.25),
                np.arange(11.0, 31.0, 1.0),
                np.arange(35.0, 61.0, 5.0),
            )
        )
    )


def _release_plot_times() -> np.ndarray:
    return np.unique(
        np.concatenate(
            (
                np.array([0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]),
                np.arange(60.0, END_TIME_S + 1.0, 60.0),
            )
        )
    )


def _bar_panel(value: float) -> dict[str, Any]:
    return {"x": [0.5], "y": [float(value)], "value": float(value)}


def _line_panel(
    raw_t: np.ndarray,
    raw_y: np.ndarray,
    plot_t_s: np.ndarray,
    *,
    time_divisor: float = 1.0,
) -> dict[str, Any]:
    return {
        "x": (plot_t_s / time_divisor).tolist(),
        "y": np.interp(plot_t_s, raw_t, raw_y).tolist(),
    }


def _panels(result: SimulationResult) -> list[dict[str, Any]]:
    state = result.state_timeseries
    release = result.ev_timeseries
    acute_t = _acute_plot_times()
    return [
        _bar_panel(result.parameters_used["dosimetry"]["joule_heat_density_J_m3"] / 1e6),
        _bar_panel(result.parameters_used["electrodynamics"]["delta_Vm"]),
        _line_panel(result.t_array, state["Ca_i"].to_numpy(), acute_t),
        _line_panel(result.t_array, state["repair_state"].to_numpy(), acute_t),
        _line_panel(
            result.t_array,
            release["sEV_cumulative"].to_numpy(),
            _release_plot_times(),
            time_divisor=60.0,
        ),
        _bar_panel(result.summary["protein_enrichment"]),
        _bar_panel(result.summary["viability_fraction"] * 100.0),
        _bar_panel(result.parameters_used["manufacturing"]["isolated_yield_sEV"]),
    ]


def _global_axis_max(value: float) -> float:
    """Round an axis ceiling to two significant digits with a little headroom."""
    if value <= 0.0:
        return 1.0
    step = 10.0 ** (math.floor(math.log10(value)) - 1)
    return float(math.ceil(1.12 * value / step) * step)


def _provenance() -> dict[str, Any]:
    source_root = Path(__file__).resolve().parents[1]
    default_path = source_root / "electro_exocytosis" / "parameters" / "default_parameters.yaml"
    model_paths = sorted((source_root / "electro_exocytosis" / "models").glob("*.py"))
    model_paths += [source_root / "electro_exocytosis" / "simulation.py"]
    return {
        "calibration_status": "Unfitted illustrative comparisons; literature-derived provisional defaults.",
        "t_end_s": END_TIME_S,
        "dt_s": OUTPUT_DT_S,
        "default_parameters_sha256": hashlib.sha256(default_path.read_bytes()).hexdigest(),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_source_sha256": {
            str(path.relative_to(source_root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in model_paths
        },
    }


def build_storyboards(studies: Sequence[StoryStudy] = STUDIES) -> StoryboardDataset:
    study_manifests: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    outputs: dict[str, SimulationResult] = {}
    for study in studies:
        row_manifests: list[dict[str, Any]] = []
        for row in study.rows:
            scenario = build_scenario(study, row)
            parameters = _parameters(row)
            result = Simulation(scenario, params_override=parameters).run()
            outputs[scenario.scenario.name] = result
            panels = _panels(result)
            row_manifests.append(
                {
                    "id": row.id,
                    "label": row.label,
                    "detail": row.detail,
                    "scenario": scenario.model_dump(mode="json"),
                    "parameters": parameters,
                    "panels": panels,
                }
            )
            summary_rows.append(
                {
                    "study": study.id,
                    "row": row.id,
                    "label": row.label,
                    **{key: value for key, value in result.summary.items() if key != "warnings"},
                    "deposited_energy_J_ml": panels[0]["value"],
                    "delta_Vm_V": panels[1]["value"],
                    "protein_cargo_score": panels[5]["value"],
                    "viability_score_percent": panels[6]["value"],
                    "recovered_small_EV_model_units": panels[7]["value"],
                    "calcium_above_1uM_duration_s": float(
                        np.trapezoid(
                            (result.state_timeseries["Ca_i"].to_numpy() > 1.0).astype(float),
                            result.t_array,
                        )
                    ),
                }
            )
        study_manifests.append(
            {
                "id": study.id,
                "title": study.title,
                "question": study.question,
                "takeaway": study.takeaway,
                "caption": (
                    f"{study.caption} {PANEL_DESCRIPTION}"
                    + "Module titles: "
                    + "; ".join(f"{module['number']}. {module['title']}" for module in MODULES)
                    + ". " + " ".join(COMMON_LIMITATIONS)
                ),
                "script": study.script,
                "rows": row_manifests,
            }
        )

    if not outputs:
        raise ValueError("At least one study is required.")
    modules = [dict(module) for module in MODULES]
    for index, module in enumerate(modules):
        if "y_max" not in module:
            maximum = max(
                max(row["panels"][index]["y"])
                for study in study_manifests
                for row in study["rows"]
            )
            module["y_max"] = _global_axis_max(float(maximum))

    manifest = {
        "schema_version": 1,
        "provenance": _provenance(),
        "unfitted_illustrative": True,
        "generator": "examples/storyboard_eight_modules.py",
        "python_version": platform.python_version(),
        "figure_purpose": "Presentation-size figures; not a 7.1-in manuscript layout.",
        "figure_width_in": PRESENTATION_FIGSIZE[0],
        "figure_height_in": PRESENTATION_FIGSIZE[1],
        "simulation_horizon_s": END_TIME_S,
        "raw_output_interval_s": OUTPUT_DT_S,
        "curve_decimation": {
            "calcium_and_repair": "0.25 s through 10 s; 1 s through 30 s; 5 s through 60 s",
            "small_EV_release": "1-min grid, supplemented by 0–30-s early points",
        },
        "limitations": list(COMMON_LIMITATIONS),
        "modules": modules,
        "studies": study_manifests,
    }
    return StoryboardDataset(manifest, pd.DataFrame(summary_rows), outputs)


def _selected_timeseries(result: SimulationResult) -> pd.DataFrame:
    state = result.state_timeseries
    release = result.ev_timeseries
    return pd.DataFrame(
        {
            "time_s": result.t_array,
            "cytosolic_calcium_uM": state["Ca_i"].to_numpy(),
            "repair_activity_score": state["repair_state"].to_numpy(),
            "small_EV_cumulative_model_units": release["sEV_cumulative"].to_numpy(),
            "medium_large_EV_cumulative_model_units": release["mlEV_cumulative"].to_numpy(),
            "apoptotic_body_cumulative_model_units": release["AB_cumulative"].to_numpy(),
            "ROS_score": state["ROS"].to_numpy(),
            "ATP_score": state["ATP"].to_numpy(),
            "damage_score": state["damage"].to_numpy(),
        }
    )


def _format_bar_value(value: float) -> str:
    if value == 0.0:
        return "0"
    if value >= 10.0:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.3g}"


def plot_study(
    study: dict[str, Any],
    modules: Sequence[dict[str, Any]],
    output_stem: Path,
    *,
    dpi: int = 600,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, ScalarFormatter

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    with manuscript_style_context(), plt.rc_context(
        {"font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"]}
    ):
        fig, axes = plt.subplots(3, 8, figsize=PRESENTATION_FIGSIZE, squeeze=False)
        fig.subplots_adjust(
            left=0.105,
            right=0.992,
            bottom=0.115,
            top=0.755,
            wspace=0.62,
            hspace=0.40,
        )
        for column, module in enumerate(modules):
            title = textwrap.fill(module["title"], width=23, break_long_words=False)
            axes[0, column].set_title(
                f"{module['number']}. {title}", fontsize=10.5, pad=29, linespacing=1.10
            )
            inputs = ", ".join(str(value) for value in module["inputs"])
            axes[0, column].text(
                0.5, 1.055, f"Inputs: {inputs}", transform=axes[0, column].transAxes,
                ha="center", va="bottom", fontsize=8.5, color="#555555",
            )
        for row_index, row in enumerate(study["rows"]):
            for column, (module, panel) in enumerate(zip(modules, row["panels"], strict=True)):
                axis = axes[row_index, column]
                axis.set_xlim(module["x_min"], module["x_max"])
                axis.set_ylim(module["y_min"], module["y_max"])
                axis.spines[["top", "right"]].set_visible(False)
                axis.tick_params(labelsize=9.0, width=0.8, length=3.0, pad=2.5)
                axis.yaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
                formatter = ScalarFormatter(useOffset=False)
                formatter.set_scientific(False)
                axis.yaxis.set_major_formatter(formatter)
                axis.grid(axis="y", color="#D9DDE3", linewidth=0.65, alpha=0.7)
                axis.set_axisbelow(True)
                if module["kind"] == "line":
                    axis.plot(panel["x"], panel["y"], color=SERIES_COLOR, linewidth=1.8)
                    axis.set_xticks([0, 30, 60] if column in (2, 3) else [0, 15, 30])
                    if row_index == 2:
                        axis.set_xlabel(module["x_label"], fontsize=10.0, labelpad=4)
                    else:
                        axis.tick_params(labelbottom=False)
                else:
                    value = panel["value"]
                    axis.bar([0.5], [value], width=0.43, color=SERIES_COLOR, zorder=3)
                    axis.set_xticks([0.5])
                    axis.tick_params(axis="x", length=0)
                    axis.set_xticklabels([module["x_label"]] if row_index == 2 else [""])
                    axis.text(
                        0.5,
                        min(value + module["y_max"] * 0.035, module["y_max"] * 0.96),
                        _format_bar_value(value),
                        fontsize=9.5,
                        ha="center",
                        va="bottom",
                    )
                if row_index == 1:
                    axis.set_ylabel(
                        f"{module['metric_label']}\n({module['unit']})",
                        fontsize=10.0,
                        labelpad=5,
                    )

            location = axes[row_index, 0].get_position()
            fig.text(
                0.011,
                (location.y0 + location.y1) / 2,
                textwrap.fill(row["label"], width=10, break_long_words=False),
                ha="left",
                va="center",
                fontsize=11.0,
                fontweight="bold",
            )

        fig.suptitle(study["title"], fontsize=15.0, y=0.993, fontweight="normal")
        fig.legend(
            [Line2D([0], [0], color=SERIES_COLOR, linewidth=2.0)],
            ["Modeled output"],
            loc="lower center",
            bbox_to_anchor=(0.55, 0.006),
            frameon=False,
            fontsize=9.0,
        )
        for suffix in ("png", "svg"):
            save_manuscript_figure(
                fig,
                output_stem.with_suffix(f".{suffix}"),
                dpi=dpi,
                bbox_inches=None,
            )
        plt.close(fig)


def write_outputs(
    dataset: StoryboardDataset,
    outdir: Path,
    *,
    make_plots: bool = True,
    dpi: int = 600,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for subdir in ("timeseries", "snapshots", "scenarios"):
        (outdir / subdir).mkdir(exist_ok=True)
    dataset.summary.to_csv(outdir / "storyboard_summary.csv", index=False)
    (outdir / "manifest.json").write_text(
        json.dumps(dataset.manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    for study in dataset.manifest["studies"]:
        for row in study["rows"]:
            name = row["scenario"]["scenario"]["name"]
            result = dataset.outputs[name]
            _selected_timeseries(result).to_csv(outdir / "timeseries" / f"{name}.csv", index=False)
            (outdir / "scenarios" / f"{name}.json").write_text(
                json.dumps(
                    {"scenario": row["scenario"], "parameter_overrides": row["parameters"]},
                    indent=2,
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )
            (outdir / "snapshots" / f"{name}.json").write_text(
                json.dumps(
                    {
                        "summary": result.summary,
                        "parameters_used": result.parameters_used,
                        "warnings": result.warnings,
                        "panels": row["panels"],
                    },
                    indent=2,
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )
        if make_plots:
            plot_study(
                study,
                dataset.manifest["modules"],
                outdir / "figures" / study["id"],
                dpi=dpi,
            )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/eight_module_storyboards"))
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.dpi <= 0:
        raise ValueError("dpi must be positive")
    dataset = build_storyboards()
    write_outputs(dataset, args.out, make_plots=not args.no_plots, dpi=args.dpi)
    print(f"Wrote {len(dataset.outputs)} scenarios to {args.out.resolve()}")
    print(f"Manifest: {(args.out / 'manifest.json').resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
