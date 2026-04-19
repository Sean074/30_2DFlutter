"""Results serialisation: complex arrays → JSON / CSV."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


def save_flutter_results(
    velocity_range_m_s: NDArray[np.float64],
    flutter_p_modes: dict[str, NDArray[np.complex128]],
    input_path: str,
    output_dir: str | Path = "data/outputs",
) -> Path:
    """
    Serialise flutter results to a timestamped JSON file.

    Parameters
    ----------
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep in m/s.
    flutter_p_modes : dict[str, NDArray[np.complex128]]
        Mapping of mode label → complex p array.
    input_path : str
        Path of the input file used (recorded in metadata).
    output_dir : str or Path
        Directory to write results into.

    Returns
    -------
    Path
        Path to the written JSON file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    out_path = output_dir / f"pk_flutter_{timestamp}.json"

    payload: dict = {
        "run_metadata": {
            "method": "p-k flutter",
            "input_file": str(input_path),
            "timestamp_utc": timestamp,
        },
        "velocity_range_m_s": velocity_range_m_s.tolist(),
        "modes": {
            label: {
                "real_damping_1_s": np.real(p_arr).tolist(),
                "imag_freq_rad_s": np.imag(p_arr).tolist(),
            }
            for label, p_arr in flutter_p_modes.items()
        },
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    return out_path
