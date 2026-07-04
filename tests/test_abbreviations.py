from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from electro_exocytosis.abbreviations import STANDARD_ABBREVIATIONS


def test_registry_renames_columns_to_standard_language() -> None:
    frame = pd.DataFrame(
        {
            "t": [0.0, 1.0],
            "Ca_i": [0.1, 0.2],
            "ROS": [0.0, 0.4],
            "sEV_rate": [1.0, 2.0],
        }
    )

    renamed = STANDARD_ABBREVIATIONS.rename_columns(frame)

    assert list(renamed.columns) == [
        "time_s",
        "cytosolic_calcium_uM",
        "reactive_oxygen_species_relative_state",
        "small_extracellular_vesicle_rate",
    ]


def test_registry_writes_abbreviation_bundle_and_figure_note(tmp_path: Path) -> None:
    STANDARD_ABBREVIATIONS.write_bundle(tmp_path, keys=("EV", "ROS", "ATP"))

    assert (tmp_path / "abbreviations.json").exists()
    assert (tmp_path / "abbreviations.md").read_text(encoding="utf-8").startswith("# Abbreviations")
    assert "EV" in (tmp_path / "abbreviations.tex").read_text(encoding="utf-8")

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    STANDARD_ABBREVIATIONS.add_figure_note(fig, ("EV", "ROS"))

    assert any("Abbreviations:" in text.get_text() for text in fig.texts)
    plt.close(fig)
