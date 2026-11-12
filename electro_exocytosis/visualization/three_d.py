from __future__ import annotations

import json
import shutil
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

if False:  # pragma: no cover - typing only without importing optional pyvista
    from electro_exocytosis.simulation import SimulationResult


DEFAULT_HUMAN_CELL_ASSET_DIR = Path(__file__).resolve().parent / "assets" / "human_cell_vtp"

HUMAN_CELL_COMPONENT_LABELS = {
    "dna": "DNA / chromatin",
    "golgin_laite": "Golgi apparatus",
    "lysosome": "Lysosome / peroxisome",
    "mitocondria": "Mitochondria",
    "nucleous": "Nucleus / nucleolus",
    "plasma": "Plasma membrane proteins",
    "ribosomes": "Ribosomes",
    "rough_endoplasmic_reticulum": "Rough ER",
    "secretory_vesicle": "Secretory vesicles",
    "smoothe_endoplasmic_reticulum": "Smooth ER",
    "surface": "Cell surface / plasma membrane",
    "sytoskeleton_filament": "Cytoskeleton filaments",
    "transport_vesicle": "Transport vesicles",
}

COMPONENT_RENDER_STYLES = {
    "surface": {"color": "#43d34f", "opacity": 0.30, "role": "membrane_signal"},
    "plasma": {"color": "#45e36c", "opacity": 0.42, "role": "membrane_signal"},
    "mitocondria": {"color": "#9e9cff", "opacity": 0.88, "role": "mitochondrial_stress"},
    "rough_endoplasmic_reticulum": {"color": "#20c4c8", "opacity": 0.58, "role": "er_release"},
    "smoothe_endoplasmic_reticulum": {"color": "#26b6c7", "opacity": 0.58, "role": "er_release"},
    "sytoskeleton_filament": {"color": "#ff6464", "opacity": 0.82, "role": "repair"},
    "secretory_vesicle": {"color": "#f4c34f", "opacity": 0.92, "role": "ev_release"},
    "transport_vesicle": {"color": "#f4d06f", "opacity": 0.85, "role": "ev_release"},
    "lysosome": {"color": "#65dd5f", "opacity": 0.86, "role": "repair"},
    "golgin_laite": {"color": "#c28b52", "opacity": 0.84, "role": "context"},
    "nucleous": {"color": "#4649d9", "opacity": 0.60, "role": "context"},
    "dna": {"color": "#f7f2b8", "opacity": 0.95, "role": "context"},
    "ribosomes": {"color": "#3a3a3a", "opacity": 0.55, "role": "context"},
}


@dataclass(frozen=True)
class HumanCellComponent:
    """One component extracted from a VTM wrapper into VTP files."""

    name: str
    label: str
    role: str
    files: tuple[Path, ...]


@dataclass(frozen=True)
class HumanCellAssetManifest:
    """Manifest for local human-cell VTP assets."""

    asset_dir: Path
    source_dir: Path | None
    components: tuple[HumanCellComponent, ...]

    @property
    def file_count(self) -> int:
        return sum(len(component.files) for component in self.components)


@dataclass(frozen=True)
class VisualizationEventFrame:
    """One scheduled 3D frame for a mechanism-focused storyboard."""

    event_id: str
    title: str
    time_s: float
    focus_components: tuple[str, ...]
    signal_key: str
    scalar_label: str
    cmap: str
    description: str


def prepare_human_cell_vtp_assets(
    source_dir: str | Path,
    dest_dir: str | Path = DEFAULT_HUMAN_CELL_ASSET_DIR,
    components: Iterable[str] | None = None,
) -> HumanCellAssetManifest:
    """Extract VTP blocks referenced by VTM component files into one asset folder.

    The supplied cell model stores each organelle as a small VTM wrapper that points
    to one or more VTP polydata files. This helper keeps the VTP files separate,
    because separate meshes are easier to color, fade, and animate from model
    variables than one merged surface.
    """

    source = Path(source_dir)
    destination = Path(dest_dir)
    component_names = tuple(components) if components is not None else tuple(HUMAN_CELL_COMPONENT_LABELS)
    extracted_components: list[HumanCellComponent] = []

    for component_name in component_names:
        vtm_path = source / f"{component_name}.vtm"
        if not vtm_path.exists():
            raise FileNotFoundError(f"Missing VTM component file: {vtm_path}")
        vtp_refs = _vtp_refs_from_vtm(vtm_path)
        if not vtp_refs:
            raise ValueError(f"No VTP DataSet entries found in {vtm_path}")

        copied_files: list[Path] = []
        component_dir = destination / component_name
        component_dir.mkdir(parents=True, exist_ok=True)
        for ref in vtp_refs:
            source_vtp = (source / ref).resolve()
            if not source_vtp.exists():
                raise FileNotFoundError(f"VTM reference does not exist: {source_vtp}")
            dest_vtp = component_dir / source_vtp.name
            if source_vtp != dest_vtp.resolve():
                shutil.copy2(source_vtp, dest_vtp)
            copied_files.append(dest_vtp.relative_to(destination))

        style = COMPONENT_RENDER_STYLES.get(component_name, {})
        extracted_components.append(
            HumanCellComponent(
                name=component_name,
                label=HUMAN_CELL_COMPONENT_LABELS.get(component_name, component_name.replace("_", " ").title()),
                role=str(style.get("role", "context")),
                files=tuple(copied_files),
            )
        )

    manifest = HumanCellAssetManifest(
        asset_dir=destination,
        source_dir=source,
        components=tuple(extracted_components),
    )
    _write_manifest(manifest)
    return manifest


def load_human_cell_manifest(
    asset_dir: str | Path = DEFAULT_HUMAN_CELL_ASSET_DIR,
) -> HumanCellAssetManifest:
    """Load a prepared human-cell VTP asset manifest."""

    destination = Path(asset_dir)
    manifest_path = destination / "manifest.json"
    with manifest_path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    components = []
    for item in raw["components"]:
        components.append(
            HumanCellComponent(
                name=item["name"],
                label=item["label"],
                role=item["role"],
                files=tuple(Path(file_name) for file_name in item["files"]),
            )
        )
    source_dir = raw.get("source_dir")
    return HumanCellAssetManifest(
        asset_dir=destination,
        source_dir=Path(source_dir) if source_dir else None,
        components=tuple(components),
    )


def scene_values_from_result(result: "SimulationResult", time_s: float | None = None) -> dict[str, float | str]:
    """Reduce a simulation result to normalized scalar controls for a 3D scene."""

    state = result.state_timeseries
    ev = result.ev_timeseries
    if state.empty:
        raise ValueError("Simulation result has no state time-series rows.")

    time_values = state["t"].to_numpy(dtype=float)
    if time_s is None:
        idx = _representative_frame_index(result)
    else:
        idx = int(np.argmin(np.abs(time_values - float(time_s))))

    ev_idx = int(np.argmin(np.abs(ev["t"].to_numpy(dtype=float) - float(time_values[idx])))) if not ev.empty else idx
    total_ev = _column_at(ev, "sEV_cumulative", ev_idx) + _column_at(ev, "mlEV_cumulative", ev_idx) + _column_at(ev, "AB_cumulative", ev_idx)
    peak_total_ev = max(
        _column_peak(ev, "sEV_cumulative") + _column_peak(ev, "mlEV_cumulative") + _column_peak(ev, "AB_cumulative"),
        1e-12,
    )

    mitochondrial_potential = _column_at(state, "mitochondrial_potential", idx, default=1.0)
    atp = _column_at(state, "ATP", idx, default=1.0)
    atp_peak = max(_column_peak(state, "ATP"), 1e-12)

    return {
        "scenario": result.scenario_name,
        "mode": result.mode,
        "time_s": float(time_values[idx]),
        "membrane_activation": _normalized_column_at(state, "pore_activation", idx),
        "membrane_permeability": _normalized_column_at(state, "membrane_permeability", idx),
        "calcium_signal": _normalized_column_at(state, "Ca_i", idx),
        "submembrane_calcium_signal": _normalized_column_at(state, "Ca_submembrane", idx),
        "osmotic_signal": _normalized_column_at(state, "osmotic_stress", idx),
        "er_release": float(np.clip(1.0 - _column_at(state, "Ca_ER", idx) / max(_column_peak(state, "Ca_ER"), 1e-12), 0.0, 1.0)),
        "ps_exposure": float(np.clip(_column_at(state, "PS_exposure", idx), 0.0, 1.0)),
        "repair_state": float(np.clip(_column_at(state, "repair_state", idx), 0.0, 1.0)),
        "actin_disruption": float(np.clip(_column_at(state, "actin_disruption", idx), 0.0, 1.0)),
        "ros_signal": _normalized_column_at(state, "ROS", idx),
        "mitochondrial_stress": float(np.clip(1.0 - mitochondrial_potential, 0.0, 1.0)),
        "atp_depletion": float(np.clip(1.0 - atp / atp_peak, 0.0, 1.0)),
        "ev_release_signal": float(np.clip(total_ev / peak_total_ev, 0.0, 1.0)),
        "viability_fraction": float(np.clip(result.summary.get("viability_fraction", 1.0), 0.0, 1.0)),
    }


def build_default_event_schedule(result: "SimulationResult") -> tuple[VisualizationEventFrame, ...]:
    """Build a compact mechanism-first sequence for presentation frames."""

    state = result.state_timeseries
    t0 = float(state["t"].iloc[0])
    final_t = float(state["t"].iloc[-1])

    membrane_t = _time_of_peak_metric(result, ("pore_activation", "membrane_permeability"), fallback=t0)
    calcium_t = _time_of_peak_metric(result, ("Ca_submembrane", "Ca_i", "Ca_mito"), fallback=membrane_t)
    mito_t = _time_of_peak_metric(result, ("ROS", "Ca_mito", "osmotic_stress", "mitochondrial_stress_proxy"), fallback=calcium_t)
    repair_t = _time_of_peak_metric(
        result,
        ("PS_exposure", "repair_state", "actin_disruption", "calpain_activity", "annexin_activity"),
        fallback=calcium_t,
    )

    return (
        VisualizationEventFrame(
            event_id="01_nsPEF_exposure",
            title="nsPEF exposure",
            time_s=t0,
            focus_components=("surface", "plasma"),
            signal_key="field_exposure",
            scalar_label="field coupling",
            cmap="cool",
            description="Electric-field axis and surface coupling before downstream state changes are emphasized.",
        ),
        VisualizationEventFrame(
            event_id="02_membrane_electrodynamics",
            title="Membrane electrodynamics",
            time_s=membrane_t,
            focus_components=("surface", "plasma"),
            signal_key="membrane_electrodynamics",
            scalar_label="pore/permeability signal",
            cmap="turbo",
            description="The plasma membrane and surface proteins carry the electroporation-like response.",
        ),
        VisualizationEventFrame(
            event_id="03_Ca_influx_ER_release",
            title="Ca2+ influx and ER release",
            time_s=calcium_t,
            focus_components=("surface", "plasma", "rough_endoplasmic_reticulum", "smoothe_endoplasmic_reticulum"),
            signal_key="calcium_transport",
            scalar_label="Ca2+ / ER signal",
            cmap="winter_r",
            description="Electropore influx, submembrane calcium, and ER calcium mobilization are emphasized.",
        ),
        VisualizationEventFrame(
            event_id="04_ionic_mito_ROS_ATP",
            title="Ionic, mitochondrial, ROS, and ATP stress",
            time_s=mito_t,
            focus_components=("mitocondria", "surface", "plasma"),
            signal_key="mitochondrial_bioenergetics",
            scalar_label="stress signal",
            cmap="magma",
            description="Ionic/osmotic perturbation is linked to mitochondrial stress, ROS, and ATP depletion.",
        ),
        VisualizationEventFrame(
            event_id="05_remodeling_repair",
            title="Ca2+-dependent remodeling and repair",
            time_s=repair_t,
            focus_components=("surface", "plasma", "sytoskeleton_filament", "lysosome"),
            signal_key="remodeling_repair",
            scalar_label="repair/remodeling signal",
            cmap="plasma",
            description="PS exposure, cytoskeletal remodeling, lysosomal repair, and resealing are emphasized.",
        ),
        VisualizationEventFrame(
            event_id="06_EV_release",
            title="EV budding and release",
            time_s=final_t,
            focus_components=("surface", "plasma", "secretory_vesicle", "transport_vesicle"),
            signal_key="ev_release",
            scalar_label="EV release signal",
            cmap="viridis",
            description="Vesicle compartments and membrane shedding/release are emphasized at the harvest-scale endpoint.",
        ),
    )


def write_event_schedule_metadata(
    result: "SimulationResult",
    path: str | Path,
) -> tuple[VisualizationEventFrame, ...]:
    """Write the scheduled 3D storyboard frame metadata."""

    schedule = build_default_event_schedule(result)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump([_event_to_dict(event) for event in schedule], handle, indent=2)
        handle.write("\n")
    return schedule


def write_scene_metadata(
    result: "SimulationResult",
    path: str | Path,
    time_s: float | None = None,
) -> dict[str, float | str]:
    """Write normalized 3D scene controls for review or downstream rendering."""

    values = scene_values_from_result(result, time_s=time_s)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(values, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return values


def render_human_cell_snapshot(
    result: "SimulationResult",
    asset_dir: str | Path = DEFAULT_HUMAN_CELL_ASSET_DIR,
    out_png: str | Path | None = None,
    time_s: float | None = None,
    show: bool = False,
    show_colorbar: bool = False,
) -> dict[str, float | str]:
    """Render a first PyVista snapshot of the simulation projected on the cell model.

    PyVista is an optional visualization dependency. Install the `viz` extra before
    calling this function in a fresh environment:

        python -m pip install -e ".[viz]"
    """

    try:
        import pyvista as pv
    except ModuleNotFoundError as exc:  # pragma: no cover - depends on optional dependency
        raise ModuleNotFoundError(
            "PyVista is required for 3D rendering. Install it with "
            '`python -m pip install -e ".[viz]"` or run the example with '
            "`--prepare-only` to only organize the VTP assets."
        ) from exc

    manifest = load_human_cell_manifest(asset_dir)
    values = scene_values_from_result(result, time_s=time_s)
    out_path = Path(out_png) if out_png is not None else None

    plotter = pv.Plotter(
        off_screen=out_path is not None and not show,
        window_size=(2400, 1600),
    )
    plotter.set_background("#111217")

    bounds = _add_human_cell_components(plotter, pv, manifest, values, show_colorbar=show_colorbar)
    _add_field_arrow(plotter, pv, bounds)
    _add_ev_particles(plotter, pv, bounds, values)
    _add_scene_text(plotter, values, title="Composite simulation state")

    _set_camera(plotter, bounds, zoom=1.55)

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plotter.show(screenshot=str(out_path), auto_close=True)
    elif show:
        plotter.show()
    else:
        plotter.close()
    return values


def render_human_cell_storyboard_frames(
    result: "SimulationResult",
    asset_dir: str | Path = DEFAULT_HUMAN_CELL_ASSET_DIR,
    outdir: str | Path = "results/human_cell_3d/storyboard_frames",
    show_colorbar: bool = True,
) -> tuple[Path, ...]:
    """Render one presentation frame per scheduled electro-exocytosis event."""

    try:
        import pyvista as pv
    except ModuleNotFoundError as exc:  # pragma: no cover - depends on optional dependency
        raise ModuleNotFoundError(
            "PyVista is required for 3D rendering. Install it with "
            '`python -m pip install -e ".[viz]"` or run the example with '
            "`--prepare-only` to only organize the VTP assets."
        ) from exc

    manifest = load_human_cell_manifest(asset_dir)
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rendered_paths: list[Path] = []

    for frame_number, event in enumerate(build_default_event_schedule(result), start=1):
        values = scene_values_from_result(result, time_s=event.time_s)
        plotter = pv.Plotter(off_screen=True, window_size=(2400, 1600))
        plotter.set_background("#111217")

        bounds = _add_human_cell_components(
            plotter,
            pv,
            manifest,
            values,
            event=event,
            show_colorbar=show_colorbar,
        )
        if event.signal_key in {"field_exposure", "membrane_electrodynamics"}:
            _add_field_arrow(plotter, pv, bounds)
        if event.signal_key in {"calcium_transport", "mitochondrial_bioenergetics"}:
            _add_event_glow(plotter, pv, bounds, values, event)
        if event.signal_key == "ev_release":
            _add_ev_particles(plotter, pv, bounds, values)

        _add_scene_text(plotter, values, title=event.title, caption=event.scalar_label)
        _set_camera(plotter, bounds, zoom=1.65)

        out_path = output_dir / f"{frame_number:02d}_{event.event_id}.png"
        plotter.show(screenshot=str(out_path), auto_close=True)
        rendered_paths.append(out_path)

    return tuple(rendered_paths)


def _vtp_refs_from_vtm(vtm_path: Path) -> tuple[Path, ...]:
    tree = ET.parse(vtm_path)
    refs: list[Path] = []
    for element in tree.iter():
        if element.tag.split("}")[-1] != "DataSet":
            continue
        file_name = element.attrib.get("file")
        if file_name and file_name.endswith(".vtp"):
            refs.append(Path(file_name))
    return tuple(refs)


def _write_manifest(manifest: HumanCellAssetManifest) -> None:
    manifest.asset_dir.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "description": "VTP blocks extracted from the human cell VTM component wrappers for 3D electro-exocytosis visualization.",
        "source_dir": str(manifest.source_dir) if manifest.source_dir is not None else None,
        "components": [
            {
                "name": component.name,
                "label": component.label,
                "role": component.role,
                "files": [str(file_name) for file_name in component.files],
                "style": COMPONENT_RENDER_STYLES.get(component.name, {}),
            }
            for component in manifest.components
        ],
    }
    with (manifest.asset_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _representative_frame_index(result: "SimulationResult") -> int:
    state = result.state_timeseries
    scores = np.zeros(len(state), dtype=float)
    for column in (
        "pore_activation",
        "Ca_submembrane",
        "Ca_i",
        "PS_exposure",
        "repair_state",
        "ROS",
        "repair_shedding_rate",
    ):
        if column in state:
            scores += _normalized_series(state[column].to_numpy(dtype=float))
    if result.ev_timeseries is not None and not result.ev_timeseries.empty:
        ev_total = (
            result.ev_timeseries.get("sEV_cumulative", 0.0)
            + result.ev_timeseries.get("mlEV_cumulative", 0.0)
            + result.ev_timeseries.get("AB_cumulative", 0.0)
        )
        scores += np.interp(
            state["t"].to_numpy(dtype=float),
            result.ev_timeseries["t"].to_numpy(dtype=float),
            _normalized_series(np.asarray(ev_total, dtype=float)),
        )
    return int(np.argmax(scores))


def _normalized_series(values: np.ndarray) -> np.ndarray:
    values = np.nan_to_num(np.asarray(values, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    shifted = values - min(float(values.min()), 0.0)
    peak = max(float(shifted.max()), 1e-12)
    return np.clip(shifted / peak, 0.0, 1.0)


def _column_at(frame: Any, column: str, index: int, default: float = 0.0) -> float:
    if frame is None or column not in frame or len(frame) == 0:
        return float(default)
    bounded_index = min(max(index, 0), len(frame) - 1)
    return float(frame[column].iloc[bounded_index])


def _column_peak(frame: Any, column: str, default: float = 0.0) -> float:
    if frame is None or column not in frame or len(frame) == 0:
        return float(default)
    return float(np.nanmax(np.abs(frame[column].to_numpy(dtype=float))))


def _normalized_column_at(frame: Any, column: str, index: int) -> float:
    if frame is None or column not in frame or len(frame) == 0:
        return 0.0
    series = _normalized_series(frame[column].to_numpy(dtype=float))
    bounded_index = min(max(index, 0), len(series) - 1)
    return float(series[bounded_index])


def _time_of_peak_metric(result: "SimulationResult", columns: tuple[str, ...], fallback: float) -> float:
    state = result.state_timeseries
    if state.empty:
        return fallback
    score = np.zeros(len(state), dtype=float)
    for column in columns:
        if column == "mitochondrial_stress_proxy" and "mitochondrial_potential" in state:
            score += _normalized_series(1.0 - state["mitochondrial_potential"].to_numpy(dtype=float))
        elif column in state:
            score += _normalized_series(state[column].to_numpy(dtype=float))
    if float(score.max()) <= 0.0:
        return fallback
    return float(state["t"].iloc[int(np.argmax(score))])


def _event_to_dict(event: VisualizationEventFrame) -> dict[str, Any]:
    return {
        "event_id": event.event_id,
        "title": event.title,
        "time_s": event.time_s,
        "focus_components": list(event.focus_components),
        "signal_key": event.signal_key,
        "scalar_label": event.scalar_label,
        "cmap": event.cmap,
        "description": event.description,
    }


def _scalar_bar_args(title: str) -> dict[str, Any]:
    return {
        "title": title,
        "title_font_size": 8,
        "label_font_size": 7,
        "n_labels": 3,
        "position_x": 0.86,
        "position_y": 0.11,
        "width": 0.08,
        "height": 0.24,
        "fmt": "%.1f",
    }


def _set_camera(plotter: Any, bounds: tuple[float, float, float, float, float, float], zoom: float = 1.55) -> None:
    center = ((bounds[0] + bounds[1]) / 2.0, (bounds[2] + bounds[3]) / 2.0, (bounds[4] + bounds[5]) / 2.0)
    extent = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4], 1.0)
    plotter.camera_position = [
        (center[0] + 0.25 * extent, center[1] - 2.9 * extent, center[2] + 1.35 * extent),
        center,
        (0.0, 0.0, 1.0),
    ]
    plotter.camera.zoom(zoom)


def _add_human_cell_components(
    plotter: Any,
    pv: Any,
    manifest: HumanCellAssetManifest,
    values: dict[str, float | str],
    event: VisualizationEventFrame | None = None,
    show_colorbar: bool = False,
) -> tuple[float, float, float, float, float, float]:
    bounds: list[float] | None = None
    scalar_bar_added = False
    for component in manifest.components:
        style = COMPONENT_RENDER_STYLES.get(component.name, {"color": "#b8c0cc", "opacity": 0.7, "role": "context"})
        is_focused = event is None or component.name in event.focus_components
        opacity = float(style.get("opacity", 0.7))
        if event is not None and not is_focused:
            opacity = min(opacity, 0.13)
        for rel_path in component.files:
            mesh = pv.read(manifest.asset_dir / rel_path)
            mesh = mesh.copy(deep=True)
            _add_component_signal(mesh, component.name, values, event=event)
            role = str(style.get("role", "context"))
            show_this_colorbar = bool(show_colorbar and is_focused and not scalar_bar_added)
            if event is not None and is_focused:
                plotter.add_mesh(
                    mesh,
                    scalars="event_signal",
                    cmap=event.cmap,
                    clim=(0.0, 1.0),
                    opacity=opacity,
                    smooth_shading=True,
                    show_scalar_bar=show_this_colorbar,
                    scalar_bar_args=_scalar_bar_args(event.scalar_label) if show_this_colorbar else None,
                )
                scalar_bar_added = scalar_bar_added or show_this_colorbar
            elif role == "membrane_signal":
                plotter.add_mesh(
                    mesh,
                    scalars="nsPEF_membrane_signal",
                    cmap="turbo",
                    clim=(0.0, 1.0),
                    opacity=opacity,
                    smooth_shading=True,
                    show_scalar_bar=show_this_colorbar,
                    scalar_bar_args=_scalar_bar_args("membrane signal") if show_this_colorbar else None,
                )
                scalar_bar_added = scalar_bar_added or show_this_colorbar
            elif role == "mitochondrial_stress":
                plotter.add_mesh(
                    mesh,
                    scalars="mitochondrial_stress",
                    cmap="magma",
                    clim=(0.0, 1.0),
                    opacity=opacity,
                    smooth_shading=True,
                    show_scalar_bar=show_this_colorbar,
                    scalar_bar_args=_scalar_bar_args("mitochondrial stress") if show_this_colorbar else None,
                )
                scalar_bar_added = scalar_bar_added or show_this_colorbar
            elif role == "er_release":
                plotter.add_mesh(
                    mesh,
                    scalars="er_release",
                    cmap="winter_r",
                    clim=(0.0, 1.0),
                    opacity=opacity,
                    smooth_shading=True,
                    show_scalar_bar=show_this_colorbar,
                    scalar_bar_args=_scalar_bar_args("ER release") if show_this_colorbar else None,
                )
                scalar_bar_added = scalar_bar_added or show_this_colorbar
            elif role in {"repair", "ev_release"}:
                scalar_name = "repair_signal" if role == "repair" else "ev_release_signal"
                plotter.add_mesh(
                    mesh,
                    scalars=scalar_name,
                    cmap="plasma",
                    clim=(0.0, 1.0),
                    opacity=opacity,
                    smooth_shading=True,
                    show_scalar_bar=show_this_colorbar,
                    scalar_bar_args=_scalar_bar_args(scalar_name.replace("_", " ")) if show_this_colorbar else None,
                )
                scalar_bar_added = scalar_bar_added or show_this_colorbar
            else:
                plotter.add_mesh(
                    mesh,
                    color=str(style.get("color", "#b8c0cc")) if is_focused else "#6f7480",
                    opacity=opacity,
                    smooth_shading=True,
                )
            bounds = _merge_bounds(bounds, mesh.bounds)
    if bounds is None:
        raise ValueError("No VTP meshes were loaded from the human-cell manifest.")
    return tuple(bounds)  # type: ignore[return-value]


def _add_component_signal(
    mesh: Any,
    component_name: str,
    values: dict[str, float | str],
    event: VisualizationEventFrame | None = None,
) -> None:
    n_points = int(mesh.n_points)
    if n_points <= 0:
        return
    polar_weight = np.ones(n_points, dtype=float)
    if component_name in {"surface", "plasma"}:
        points = np.asarray(mesh.points, dtype=float)
        center_x = float(points[:, 0].mean())
        span_x = max(float(np.max(np.abs(points[:, 0] - center_x))), 1e-12)
        polar_weight = np.abs(points[:, 0] - center_x) / span_x
        polar_weight = polar_weight**1.8
        membrane_activation = float(values["membrane_activation"])
        ps_exposure = float(values["ps_exposure"])
        repair_state = float(values["repair_state"])
        mesh.point_data["nsPEF_membrane_signal"] = np.clip(
            0.72 * membrane_activation * polar_weight + 0.18 * ps_exposure + 0.10 * repair_state,
            0.0,
            1.0,
        )
    mesh.point_data["mitochondrial_stress"] = np.full(n_points, float(values["mitochondrial_stress"]))
    mesh.point_data["er_release"] = np.full(n_points, float(values["er_release"]))
    mesh.point_data["repair_signal"] = np.full(
        n_points,
        max(float(values["repair_state"]), float(values["actin_disruption"])),
    )
    mesh.point_data["ev_release_signal"] = np.full(n_points, float(values["ev_release_signal"]))
    if event is not None:
        mesh.point_data["event_signal"] = _event_signal_for_component(component_name, polar_weight, values, event)


def _event_signal_for_component(
    component_name: str,
    polar_weight: np.ndarray,
    values: dict[str, float | str],
    event: VisualizationEventFrame,
) -> np.ndarray:
    if event.signal_key == "field_exposure":
        if component_name in {"surface", "plasma"}:
            return np.clip(0.18 + 0.82 * polar_weight, 0.0, 1.0)
        return np.full_like(polar_weight, 0.10, dtype=float)
    if event.signal_key == "membrane_electrodynamics":
        signal = max(float(values["membrane_activation"]), float(values["membrane_permeability"]))
        if component_name in {"surface", "plasma"}:
            return np.clip(0.10 + 0.90 * signal * polar_weight, 0.0, 1.0)
        return np.full_like(polar_weight, 0.12 * signal, dtype=float)
    if event.signal_key == "calcium_transport":
        if component_name in {"surface", "plasma"}:
            return np.clip(float(values["submembrane_calcium_signal"]) * (0.35 + 0.65 * polar_weight), 0.0, 1.0)
        if component_name in {"rough_endoplasmic_reticulum", "smoothe_endoplasmic_reticulum"}:
            return np.full_like(polar_weight, max(float(values["calcium_signal"]), float(values["er_release"])), dtype=float)
        return np.full_like(polar_weight, 0.10 * float(values["calcium_signal"]), dtype=float)
    if event.signal_key == "mitochondrial_bioenergetics":
        stress = max(
            float(values["mitochondrial_stress"]),
            float(values["ros_signal"]),
            float(values["atp_depletion"]),
            float(values["osmotic_signal"]),
        )
        if component_name == "mitocondria":
            return np.full_like(polar_weight, stress, dtype=float)
        if component_name in {"surface", "plasma"}:
            return np.clip(float(values["osmotic_signal"]) * (0.25 + 0.75 * polar_weight), 0.0, 1.0)
        return np.full_like(polar_weight, 0.10 * stress, dtype=float)
    if event.signal_key == "remodeling_repair":
        repair = max(float(values["repair_state"]), float(values["actin_disruption"]), float(values["ps_exposure"]))
        if component_name in {"surface", "plasma"}:
            return np.clip(max(float(values["ps_exposure"]), float(values["repair_state"])) * (0.25 + 0.75 * polar_weight), 0.0, 1.0)
        if component_name in {"sytoskeleton_filament", "lysosome"}:
            return np.full_like(polar_weight, repair, dtype=float)
        return np.full_like(polar_weight, 0.10 * repair, dtype=float)
    if event.signal_key == "ev_release":
        ev_signal = float(values["ev_release_signal"])
        if component_name in {"surface", "plasma"}:
            return np.clip(ev_signal * (0.20 + 0.80 * polar_weight), 0.0, 1.0)
        if component_name in {"secretory_vesicle", "transport_vesicle"}:
            return np.full_like(polar_weight, ev_signal, dtype=float)
        return np.full_like(polar_weight, 0.08 * ev_signal, dtype=float)
    return np.zeros_like(polar_weight, dtype=float)


def _merge_bounds(bounds: list[float] | None, next_bounds: tuple[float, float, float, float, float, float]) -> list[float]:
    if bounds is None:
        return list(next_bounds)
    bounds[0] = min(bounds[0], next_bounds[0])
    bounds[1] = max(bounds[1], next_bounds[1])
    bounds[2] = min(bounds[2], next_bounds[2])
    bounds[3] = max(bounds[3], next_bounds[3])
    bounds[4] = min(bounds[4], next_bounds[4])
    bounds[5] = max(bounds[5], next_bounds[5])
    return bounds


def _add_field_arrow(plotter: Any, pv: Any, bounds: tuple[float, float, float, float, float, float]) -> None:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    dx = max(xmax - xmin, 1.0)
    y_mid = (ymin + ymax) / 2.0
    z_mid = zmax + 0.16 * max(zmax - zmin, 1.0)
    arrow = pv.Arrow(start=(xmin - 0.28 * dx, y_mid, z_mid), direction=(1.0, 0.0, 0.0), scale=0.22 * dx)
    plotter.add_mesh(arrow, color="#7bdcff", opacity=0.82)
    plotter.add_text("nsPEF", position=(0.49, 0.90), viewport=True, font_size=10, color="#b9efff")


def _add_event_glow(
    plotter: Any,
    pv: Any,
    bounds: tuple[float, float, float, float, float, float],
    values: dict[str, float | str],
    event: VisualizationEventFrame,
) -> None:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    center = ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0)
    radius = 0.28 * max(xmax - xmin, ymax - ymin, zmax - zmin, 1.0)
    if event.signal_key == "calcium_transport":
        signal = max(float(values["calcium_signal"]), float(values["submembrane_calcium_signal"]))
        color = "#5bd7ff"
    elif event.signal_key == "mitochondrial_bioenergetics":
        signal = max(float(values["ros_signal"]), float(values["mitochondrial_stress"]), float(values["osmotic_signal"]))
        color = "#ff9b4a"
    else:
        return
    if signal <= 0.02:
        return
    sphere = pv.Sphere(radius=radius, center=center, theta_resolution=64, phi_resolution=32)
    plotter.add_mesh(sphere, color=color, opacity=min(0.20, 0.05 + 0.14 * signal), smooth_shading=True)


def _add_ev_particles(plotter: Any, pv: Any, bounds: tuple[float, float, float, float, float, float], values: dict[str, float | str]) -> None:
    ev_signal = float(values["ev_release_signal"])
    if ev_signal <= 0.01:
        return
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    center = np.array([(xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0], dtype=float)
    radii = np.array([max(xmax - xmin, 1.0), max(ymax - ymin, 1.0), max(zmax - zmin, 1.0)], dtype=float) / 2.0
    rng = np.random.default_rng(17)
    particle_count = 10 + int(45 * ev_signal)
    particle_radius = 0.018 * max(radii)
    for index in range(particle_count):
        phi = rng.uniform(0.0, 2.0 * np.pi)
        costheta = rng.uniform(-0.62, 0.88)
        sintheta = np.sqrt(max(1.0 - costheta * costheta, 0.0))
        direction = np.array([sintheta * np.cos(phi), sintheta * np.sin(phi), costheta])
        distance = 1.06 + 0.22 * rng.random() + 0.20 * ev_signal
        position = center + distance * radii * direction
        sphere = pv.Sphere(radius=particle_radius * (0.72 + 0.8 * rng.random()), center=position)
        color = "#55d8ff" if index % 3 else "#ffd45d"
        plotter.add_mesh(sphere, color=color, opacity=0.82, smooth_shading=True)


def _add_scene_text(
    plotter: Any,
    values: dict[str, float | str],
    title: str,
    caption: str | None = None,
) -> None:
    label = (
        f"{title}\n"
        f"t = {float(values['time_s']):.1f} s | "
        f"Ca {float(values['calcium_signal']):.2f} | "
        f"PS {float(values['ps_exposure']):.2f} | "
        f"repair {float(values['repair_state']):.2f} | "
        f"EV {float(values['ev_release_signal']):.2f} | "
        f"viability {float(values['viability_fraction']):.2f}"
    )
    if caption:
        label = f"{label}\n{caption}"
    plotter.add_text(label, position=(0.015, 0.945), viewport=True, font_size=9, color="#f7f7f2")
