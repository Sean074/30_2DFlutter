"""Input parsing: JSON / CSV → dicts."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def load_json(path: str | Path) -> dict:
    """
    Load a JSON input file.

    Parameters
    ----------
    path : str or Path
        Path to the JSON file.

    Returns
    -------
    dict
        Parsed contents of the file.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file is not valid JSON.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc


def load_csv(path: str | Path) -> list[dict]:
    """
    Load a CSV input file with a header row.

    Parameters
    ----------
    path : str or Path
        Path to the CSV file.

    Returns
    -------
    list[dict]
        One dict per data row, keyed by header names.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))
