# Engineering Code Standard

## 1. Scope

This standard applies to all Python-based engineering analysis tools in this project. It defines conventions for structure, naming, data I/O, user interaction, and testing. It is intended to produce code that is readable, reproducible, and maintainable by engineers who did not write it.

---

## 2. Project Structure

```
project/
├── main.py                  # Entry point — TUI and orchestration only
├── defaults/
│   ├── constants.json       # Universal physical and mathematical constants
│   └── default_flutter.json # Default parameter values for p-k flutter analysis
├── data/
│   ├── inputs/              # Input files: JSON or CSV
│   └── outputs/             # Generated results: JSON, CSV, PNG
├── analysis/
│   ├── __init__.py
│   └── <method_name>.py     # One module per analysis method
├── io/
│   ├── __init__.py
│   ├── loader.py            # Input parsing (JSON / CSV → dicts / dataclasses)
│   └── exporter.py          # Results serialisation
├── ui/
│   ├── __init__.py
│   └── tui.py               # All TUI logic (menus, prompts, validation)
├── plotting/
│   ├── __init__.py
│   └── plots.py             # All matplotlib / visualisation logic
├── tests/
│   └── test_<module>.py
├── docs/
│   └── engineering_code_standard.md
├── requirements.txt
└── CLAUDE.md
```

### Rules

- `main.py` **must not** contain analysis logic, I/O parsing, or plotting code.
- Each analysis method lives in its own module under `analysis/`.
- Tests mirror the source structure: `tests/test_flutter.py` tests `analysis/flutter.py`.

---

## 3. Entry Point (`main.py`)

`main.py` is the sole runnable entry point.

```python
"""Entry point for the engineering analysis tool."""
from ui.tui import run_tui

if __name__ == "__main__":
    run_tui()
```

`run_tui()` handles all user interaction and dispatches to analysis modules.

---

## 4. Text User Interface (TUI)

See `docs/user_interface.md` for the full TUI style guide. Key rules:

- Text-only, 1990s-style ANSI terminal output. No `rich`, `questionary`, `curses`, or animation.
- Four mandatory stages in order: INPUT → CHECK INPUT → OUTPUT → CHECK OUTPUT.
- Never `print()` inside analysis or I/O modules — terminal output only in `ui/tui.py`.
- Provide a `--batch` CLI flag (via `argparse`) so scripts can bypass the interactive TUI.

---

## 5. Defaults and Constants

### Rule: no hardcoded values in program code

**Parameter defaults and physical/mathematical constants must never be hardcoded in Python source files.** They must be loaded from the JSON files in `defaults/` at runtime via `file_io/loader.py`.

This applies to:
- Default input values shown as prompts in the TUI (e.g. default semi-chord, default air density)
- Mathematical constants (π)
- Physical constants (standard gravity, sea-level ISA density, speed of sound)
- Solver parameters: convergence tolerances, iteration limits, and numerical seed values (e.g. `det_tol_nd`, `p_init_damp_1_nd`)
- Any value that would otherwise appear as a bare numeric literal in analysis or UI code

### `defaults/constants.json`

Stores universal physical and mathematical constants. Keys follow the `<name>_<unit>` convention.

```json
{
  "schema_version": 1,
  "pi_nd": 3.141592653589793,
  "g_m_s2": 9.80665,
  "rho_sl_kg_m3": 1.225,
  "speed_of_sound_sl_m_s": 340.29,
  "gamma_air_nd": 1.4
}
```

### `defaults/default_flutter.json`

Stores default parameter values for the p-k flutter analysis. This includes **all** defaults: values pre-filled in TUI prompts, solver tolerances, and numerical seed values for the iterative method.

```json
{
  "schema_version": 1,
  "b_m": 0.5,
  "a_h_nd": 0.25,
  "x_alpha_nd": -0.5,
  "mass_kg": 10.4,
  "inertia_alpha_kg_m2": 1.12,
  "omega_alpha_rad_s": 26.0,
  "omega_heave_rad_s": 20.0,
  "rho_kg_m3": 1.21,
  "velocity_range_m_s": [0.01, 50.0, 4990],
  "det_tol_nd": 0.001,
  "p_init_damp_1_nd": -0.01,
  "p_init_damp_2_nd": 0.0
}
```

`det_tol_nd` is the secant-method convergence tolerance on `|det(F)|`. `p_init_damp_1_nd` and `p_init_damp_2_nd` are the damping scale factors for the two initial secant seeds: `real(p) = p_init_damp_*_nd * k`.

### Loading pattern

Load defaults once at TUI startup via `file_io/loader.py`. Pass the loaded dict into the prompting functions — do not re-open the files inside analysis modules.

```python
# In ui/tui.py — load once at startup
from file_io.loader import load_json

_DEFAULTS = load_json("defaults/default_flutter.json")
_CONSTANTS = load_json("defaults/constants.json")

# Use loaded values as prompt defaults
b_m = _prompt_float("SEMI-CHORD b", "m", default=_DEFAULTS["b_m"], gt=0.0)
```

```python
# In analysis/flutter.py — receive constants as arguments, never import directly
def compute_flutter(
    ...,
    pi_nd: float,  # passed in from _CONSTANTS["pi_nd"]
) -> NDArray[np.complex128]:
    ...
```

### What is forbidden

```python
# FORBIDDEN — hardcoded numeric literal
PI = 3.14159

# FORBIDDEN — numpy constant used directly in analysis module
import numpy as np
pi = np.pi

# FORBIDDEN — hardcoded default in a prompt
b_m = _prompt_float("SEMI-CHORD b", "m", default=0.5)

# CORRECT — value loaded from defaults/
b_m = _prompt_float("SEMI-CHORD b", "m", default=_DEFAULTS["b_m"])
```

> **Exception:** local arithmetic intermediates inside a ≤10-line function (e.g. `2.0` in `2.0 * pi_nd * freq`) are permitted. What is forbidden is defining named constants or default-value literals in source code.

---

## 6. Variable Naming Convention

### Format

```
<name>_<unit>
```

Unit suffixes are appended with an underscore. Use SI base and coherent derived units only.

| Quantity             | Unit              | Suffix         | Example                          |
|----------------------|-------------------|----------------|----------------------------------|
| Length               | metre             | `_m`           | `chord_m`, `span_m`              |
| Mass                 | kilogram          | `_kg`          | `mass_kg`, `m_total_kg`          |
| Time                 | second            | `_s`           | `period_s`, `dt_s`               |
| Velocity             | metre/second      | `_m_s`         | `velocity_m_s`, `v_flutter_m_s`  |
| Frequency (rad)      | radian/second     | `_rad_s`       | `omega_h_rad_s`, `omega_a_rad_s` |
| Frequency (Hz)       | hertz             | `_hz`          | `freq_hz`                        |
| Density              | kilogram/metre³   | `_kg_m3`       | `rho_kg_m3`                      |
| Pressure             | pascal            | `_pa`          | `q_dyn_pa`                       |
| Stiffness            | newton/metre      | `_n_m`         | `k_h_n_m`                        |
| Rotational stiffness | newton·metre/rad  | `_nm_rad`      | `k_a_nm_rad`                     |
| Damping              | dimensionless     | `_nd`          | `g_nd`, `damping_nd`             |
| Reduced frequency    | dimensionless     | `_nd`          | `k_nd`                           |
| Angle                | radian            | `_rad`         | `alpha_rad`, `theta_rad`         |

### Additional rules

- **No abbreviations without a suffix.** `k` alone is forbidden; `k_nd` or `k_h_n_m` is required.
- **NumPy arrays:** append `_arr` before the unit suffix only if the shape is non-obvious from context, e.g. `p_arr_rad_s`.
- **Dimensionless ratios** use `_nd` ("non-dimensional").
- **Loop indices** are exempt: `i`, `j`, `n` are acceptable.
- **Local intermediate variables** inside a ≤10-line function are exempt from the unit suffix if the unit is apparent from the surrounding named variables.

---

## 7. Input Data

### Formats

Input files are **JSON** (preferred for structured/nested data) or **CSV** (preferred for tabular parameter sweeps).

### JSON schema example

```json
{
  "schema_version": 1,
  "chord_m": 1.0,
  "span_m": 6.0,
  "mass_kg": 10.5,
  "omega_h_rad_s": 10.0,
  "omega_a_rad_s": 25.0,
  "rho_kg_m3": 1.225,
  "velocity_range_m_s": [10.0, 200.0, 100]
}
```

- Keys **must** follow the same `<name>_<unit>` convention as Python variables.
- All values are SI. No unit conversion inside analysis modules.
- Include a `"schema_version"` field so loaders can detect incompatible files.

### CSV convention

- Row 1: header row with `<name>_<unit>` column names.
- No blank rows or comment rows.
- Numeric values only; strings only for categorical columns (e.g. `"material"`).

### Loading pattern

```python
# file_io/loader.py
import json
import csv
from pathlib import Path

def load_json(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)

def load_csv(path: str | Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))
```

---

## 8. Type Annotations

All public functions must have full type annotations.

```python
import numpy as np
from numpy.typing import NDArray

def compute_flutter(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    velocity_range_m_s: NDArray[np.float64],
    rho_kg_m3: float,
    chord_m: float,
) -> NDArray[np.complex128]:
    ...
```

Use `from __future__ import annotations` at the top of each file for forward-reference support.

---

## 9. Docstrings

Use **NumPy-style** docstrings for all public functions and classes.

```python
def compute_flutter(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    velocity_range_m_s: NDArray[np.float64],
    rho_kg_m3: float,
    chord_m: float,
) -> NDArray[np.complex128]:
    """
    Compute p-k flutter eigenvalues over a velocity sweep.

    Parameters
    ----------
    mass_matrix_kg : NDArray[np.float64], shape (n, n)
        Structural mass matrix in kg.
    stiffness_matrix_n_m : NDArray[np.float64], shape (n, n)
        Structural stiffness matrix in N/m.
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Sweep velocities in m/s.
    rho_kg_m3 : float
        Air density in kg/m³.
    chord_m : float
        Airfoil chord length in metres.

    Returns
    -------
    NDArray[np.complex128], shape (n_modes, nv)
        Complex p values: real = damping (1/s), imag = frequency (rad/s).
    """
```

---

## 10. Error Handling

- Raise descriptive exceptions at **system boundaries** (file I/O, user input).
- Do not use bare `except:` or `except Exception:` — catch specific types.
- Analysis modules raise `ValueError` for invalid physical inputs; the TUI catches and displays them cleanly.

```python
# Good
if chord_m <= 0:
    raise ValueError(f"chord_m must be positive; got {chord_m}")

# Bad
try:
    result = compute(...)
except Exception:
    pass
```

---

## 11. Code Style

- Follow **PEP 8**.
- Maximum line length: **100 characters**.
- Use `black` for auto-formatting and `ruff` for linting.
- No `import *`.
- Imports ordered: stdlib → third-party → local, separated by blank lines.

```
pip install black ruff
black .
ruff check .
```

---

## 12. Testing

- Use `pytest`.
- One test file per module: `tests/test_flutter.py` for `analysis/flutter.py`.
- Tests use only hardcoded SI values — no file I/O in unit tests.
- Numerical results checked with `np.testing.assert_allclose(actual, expected, rtol=1e-4)`.

```python
def test_flutter_speed_matches_theodorsen():
    p_arr = compute_flutter(
        mass_matrix_kg=...,
        stiffness_matrix_n_m=...,
        velocity_range_m_s=np.linspace(10, 200, 100),
        rho_kg_m3=1.225,
        chord_m=1.0,
    )
    np.testing.assert_allclose(np.real(p_arr[0, 50]), 0.0, atol=1e-3)
```

---

## 13. Dependencies (`requirements.txt`)

Pin exact versions for reproducibility. Minimum required packages:

```
numpy==2.2.4
scipy==1.15.2
matplotlib==3.10.1
pytest==8.3.5
black==25.1.0
ruff==0.11.5
```

Remove unused packages (e.g. `flask` if not serving a web interface).

---

## 14. Output / Results

- Save results to `data/outputs/` as JSON or CSV, never overwriting inputs.
- Filename convention: `<method>_<timestamp_iso8601>.json`, e.g. `pk_flutter_20260419T142301.json`.
- Plots saved as PNG at 150 dpi minimum alongside the result file.
- Log run parameters (input file path, method, timestamp) in the output JSON under a `"run_metadata"` key.

---

## 15. Quick-Reference Checklist

Before committing, verify:

- [ ] All public function parameters and return values carry `_<unit>` suffixes.
- [ ] All values are SI with no in-module unit conversion.
- [ ] `main.py` contains only entry-point code.
- [ ] Input loaded via `file_io/loader.py`, not inline `open()` calls in analysis modules.
- [ ] Default parameter values (including solver tolerances and seed values) loaded from `defaults/default_flutter.json`.
- [ ] Physical and mathematical constants loaded from `defaults/constants.json`.
- [ ] No hardcoded numeric constants or default values in Python source files.
- [ ] `black` and `ruff` pass with no errors.
- [ ] New functions have NumPy-style docstrings.
- [ ] At least one `pytest` test covers each new public function.
- [ ] Unused packages removed from `requirements.txt`.
