from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from electro_exocytosis.io.readers import load_scenario
from electro_exocytosis.simulation import Simulation
from electro_exocytosis.visualization.three_d import (
    DEFAULT_HUMAN_CELL_ASSET_DIR,
    prepare_human_cell_vtp_assets,
    render_human_cell_snapshot,
    render_human_cell_storyboard_frames,
    write_event_schedule_metadata,
    write_scene_metadata,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scenario",
        type=Path,
        default=Path("examples/scenario_layer3_ion_bioenergetics.yaml"),
        help="Scenario YAML used to drive the 3D scene controls.",
    )
    parser.add_argument(
        "--source-vtm-dir",
        type=Path,
        default=Path("docs"),
        help="Directory containing the source VTM component wrappers.",
    )
    parser.add_argument(
        "--assets",
        type=Path,
        default=DEFAULT_HUMAN_CELL_ASSET_DIR,
        help="Destination folder for extracted VTP assets.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/human_cell_3d"),
        help="Output directory for scene metadata and optional render.",
    )
    parser.add_argument(
        "--time-s",
        type=float,
        default=None,
        help="Simulation time for composite rendering. Storyboard mode chooses event-specific times.",
    )
    parser.add_argument(
        "--method",
        choices=("composite", "storyboard"),
        default="composite",
        help="Composite overlays several model states; storyboard renders one scheduled event per frame.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only extract VTP assets and write scene metadata; skip PyVista rendering.",
    )
    parser.add_argument(
        "--no-color-bars",
        action="store_true",
        help="Hide scalar bars in rendered 3D frames.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive PyVista window for composite mode instead of only writing a screenshot.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = prepare_human_cell_vtp_assets(args.source_vtm_dir, args.assets)
    scenario = load_scenario(args.scenario)
    result = Simulation(scenario).run()
    scene_values = write_scene_metadata(result, args.out / "scene_values.json", time_s=args.time_s)
    schedule = write_event_schedule_metadata(result, args.out / "storyboard_schedule.json")

    print(f"Prepared {manifest.file_count} VTP blocks in {manifest.asset_dir}")
    print(f"Wrote scene controls to {args.out / 'scene_values.json'}")
    print(f"Wrote storyboard schedule with {len(schedule)} events to {args.out / 'storyboard_schedule.json'}")

    if args.prepare_only:
        return 0

    show_colorbar = not args.no_color_bars
    try:
        if args.method == "composite":
            screenshot = args.out / "electro_exocytosis_human_cell_3d.png"
            render_human_cell_snapshot(
                result,
                asset_dir=args.assets,
                out_png=screenshot if not args.show else None,
                time_s=float(scene_values["time_s"]),
                show=args.show,
                show_colorbar=show_colorbar,
            )
        else:
            rendered_paths = render_human_cell_storyboard_frames(
                result,
                asset_dir=args.assets,
                outdir=args.out / "storyboard_frames",
                show_colorbar=show_colorbar,
            )
    except ModuleNotFoundError as exc:
        print(str(exc))
        print("VTP assets and scene metadata were still generated successfully.")
        return 0

    if args.method == "composite" and not args.show:
        print(f"Wrote PyVista snapshot to {screenshot}")
    elif args.method == "storyboard":
        print(f"Wrote {len(rendered_paths)} storyboard frames to {args.out / 'storyboard_frames'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
