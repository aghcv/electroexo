from __future__ import annotations

from pathlib import Path

import pandas as pd


class EvidenceLoader:
    """
    Loads and indexes the nsPEF EV literature evidence workbook.
    Currently reads sheets into DataFrames and maps rows to model layers.
    Future versions will use this to override placeholder parameters.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._loaded = False
        self._module_map = pd.DataFrame()
        self._literature_tracker = pd.DataFrame()

    def load(self) -> bool:
        """Load workbook sheets, returning False if the file is unavailable."""
        if not self.path.exists():
            self._loaded = False
            self._module_map = pd.DataFrame()
            self._literature_tracker = pd.DataFrame()
            return False
        self._module_map = pd.read_excel(self.path, sheet_name="Module Map")
        self._literature_tracker = pd.read_excel(self.path, sheet_name="Literature Tracker")
        self._loaded = True
        return True

    def get_module_map(self) -> pd.DataFrame:
        """Return the module map sheet."""
        return self._module_map.copy()

    def get_literature_tracker(self) -> pd.DataFrame:
        """Return the literature tracker sheet."""
        return self._literature_tracker.copy()

    def get_placeholder_report(self) -> dict[str, object]:
        """Report module placeholder coverage."""
        layers = []
        if self._loaded and not self._module_map.empty and "Layer" in self._module_map.columns:
            layers = [str(value) for value in self._module_map["Layer"].dropna().tolist()]
        return {
            "loaded": self._loaded,
            "placeholder_fraction": 1.0,
            "placeholder_modules": [
                "pulse",
                "dosimetry",
                "electrodynamics",
                "ion_transport",
                "remodeling_repair",
                "ev_release",
                "cargo_potency",
                "injury_quality",
                "manufacturing_qc",
                "cell_state",
            ],
            "layers": layers,
        }

    def summarize(self) -> str:
        """Return a short text summary of workbook status."""
        if not self._loaded:
            return f"Evidence workbook not loaded: {self.path}"
        return (
            f"Evidence workbook loaded with {len(self._literature_tracker)} literature rows "
            f"and {len(self._module_map)} module-map rows."
        )
