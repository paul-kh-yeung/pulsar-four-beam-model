"""JSON configuration helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict


def read_json(path: str | Path) -> Dict[str, Any]:
    """Read a JSON file into a dictionary."""
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(data: Dict[str, Any], path: str | Path) -> None:
    """Write a dictionary to a JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
