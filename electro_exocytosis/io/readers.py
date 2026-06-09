from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from electro_exocytosis.config import SimulationScenario



def load_scenario(path: str | Path) -> SimulationScenario:
    """Load and validate a YAML scenario file."""
    with Path(path).open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    return SimulationScenario.model_validate(raw)



def load_default_parameters() -> dict[str, Any]:
    """Load default_parameters.yaml from package data."""
    path = Path(__file__).resolve().parents[1] / "parameters" / "default_parameters.yaml"
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}



def merge_parameters(defaults: dict[str, Any], overrides: dict[str, Any] | None) -> dict[str, Any]:
    """Deep-merge parameter dicts."""
    merged = deepcopy(defaults)
    if not overrides:
        return merged
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_parameters(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged
