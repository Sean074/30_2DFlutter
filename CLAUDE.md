# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A 2D aeroelastic flutter analysis tool implementing the **p-k method** for determining flutter onset speed. Given structural and aerodynamic properties of a 2D airfoil section, it sweeps a velocity range and computes the complex eigenvalue `p` (where `real(p)` = damping, `imag(p)` = frequency) at each speed.

## Running

```bash
pip install -r requirements.txt
python main.py
```

`main.py` starts the terminal UI (`ui/tui.py`); results go to `data/outputs/` as JSON/CSV/PNG. To build a single-file executable (`flutter2d`): `bash build.sh` on macOS, `build.bat` on Windows (both run PyInstaller on `build.spec`).

## Architecture

- **`analysis/flutter.py`** — the p-k solver: `build_mass_matrix`, `build_stiffness_matrix`, `build_aero_matrix` (Theodorsen via the Jones two-pole approximation), `compute_structural_frequencies` and `compute_flutter`, which sweeps velocity and, for every mode, solves the eigenproblem of `F = (V/b)²p²M + K − ρπb V² A(k)` at fixed `k`, keeps the root nearest the mode's estimate and iterates `k = imag(p)`. Modes are seeded from the coupled structural frequencies and extrapolated speed to speed, so identity is kept; it returns shape `(n_modes, nv)` and raises if two modes collapse onto one root. `find_flutter_speed` interpolates the first g = 0 crossing of one mode (speed and frequency).
- **`file_io/`** — `loader.py` (JSON/CSV input), `exporter.py` (results out).
- **`plotting/plots.py`** — V-g / V-f figures: top subplot `imag(p)·V/b` (frequency, rad/s) vs V; bottom `real(p)/imag(p)` (damping `g`) vs V, with a `g = 0.03` limit line. Flutter onset is where `g` crosses zero.
- **`ui/tui.py`** — the four-stage terminal workflow; the only module that prints.
- **`tests/`** — mirrors the source tree (`tests/test_flutter.py` ↔ `analysis/flutter.py`).

The pre-refactor scripts (`flut.py`, `h1_aeroelastic_calc.py`) were removed on 21 Oct 2026; their history is in git.

## Code Standards

All code generation MUST conform to:

- **`docs/engineering_code_standard.md`** — project structure, variable naming (`<name>_<unit>` SI suffixes), type annotations, docstrings, error handling, testing, and output conventions.
- **`docs/user_interface.md`** — TUI style: 1990s plain ANSI terminal, no `rich`/curses/animation, strict 4-stage workflow (INPUT → CHECK INPUT → OUTPUT → CHECK OUTPUT).

**MUST be followed for all code in this project.**

## Project Structure

```
project/
├── main.py                  # Entry point — TUI and orchestration only
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
└── requirements.txt
```

- `main.py` **must not** contain analysis logic, I/O parsing, or plotting code.
- Each analysis method lives in its own module under `analysis/`.
- Tests mirror the source structure: `tests/test_flutter.py` tests `analysis/flutter.py`.

## Entry Point (`main.py`)

`main.py` is the sole runnable entry point and must only call `run_tui()` from `ui/tui.py`.

## TUI Style

**This is a 1990s-style command-line engineering program. Text-only, functional, direct. No decoration for its own sake.**

### Do NOT use
- `rich`, `textual`, `curses`, or any library requiring a full terminal
- Spinners, progress bars with ASCII animation
- Box-drawing characters (`┌`, `│`, `─`, etc.)
- Emoji or unicode beyond standard ASCII
- 256-colour codes, colour gradients, background colours, blinking
- Interactive menus or arrow-key navigation

### Terminal output rules
- Use only standard 16-colour ANSI codes: cyan (bold) for headers, red for errors, yellow for warnings, green for OK, cyan for calculated output.
- All caps for section titles and field labels.
- Pad labels to 28 chars so values align. Width: target 72 chars, never exceed 80.
- Use `=` for major dividers, `-` for minor sub-dividers.
- Never `print()` inside analysis or I/O modules — terminal output only in `ui/`.

### Workflow stages (in order, never combined or skipped)
1. **INPUT** — geometry, mass, structural, aero parameters (one group at a time)
2. **CHECK INPUT** — echo all values; user confirms with Y/Enter or re-enters with N
3. **OUTPUT** — print progress and results as solver runs; do not buffer
4. **CHECK OUTPUT** — summary block with flutter speed, warnings (`[W]`), errors (`[E]`)

### Input field format
```
  SEMI-CHORD b [m]             : _
```
- Label left-aligned, padded to 28 chars; unit in brackets; colon and cursor on same line.
- Default values shown in brackets: `[0.5]`. Press Enter to accept.
- `q`/`Q` exits cleanly at any prompt.

### Error/warning format
```
  [E] INVALID INPUT: b must be > 0. Re-enter.
  [W] x_a outside typical range [-0.5, 0.5]. Confirm value is correct.
```
- `[E]` = fatal for that field; re-prompt immediately. No stack traces shown to user.
- `[W]` = advisory; shown at check-input and check-output stages.

## Variable Naming Convention

Format: `<name>_<unit>` — unit suffixes appended with underscore, SI units only.

| Quantity             | Suffix         | Example                          |
|----------------------|----------------|----------------------------------|
| Length [m]           | `_m`           | `chord_m`, `span_m`              |
| Mass [kg]            | `_kg`          | `mass_kg`                        |
| Time [s]             | `_s`           | `period_s`, `dt_s`               |
| Velocity [m/s]       | `_m_s`         | `velocity_m_s`, `v_flutter_m_s`  |
| Frequency [rad/s]    | `_rad_s`       | `omega_h_rad_s`, `omega_a_rad_s` |
| Frequency [Hz]       | `_hz`          | `freq_hz`                        |
| Density [kg/m³]      | `_kg_m3`       | `rho_kg_m3`                      |
| Pressure [Pa]        | `_pa`          | `q_dyn_pa`                       |
| Stiffness [N/m]      | `_n_m`         | `k_h_n_m`                        |
| Rot. stiffness [Nm/rad] | `_nm_rad`   | `k_a_nm_rad`                     |
| Dimensionless        | `_nd`          | `g_nd`, `k_nd`, `damping_nd`     |
| Angle [rad]          | `_rad`         | `alpha_rad`, `theta_rad`         |

- `k` alone is forbidden; `k_nd` or `k_h_n_m` is required.
- NumPy arrays: append `_arr` before unit suffix if shape is non-obvious, e.g. `p_arr_rad_s`.
- Loop indices (`i`, `j`, `n`) and local intermediates in ≤10-line functions are exempt.

## Input Data

- JSON (preferred for structured data) or CSV (preferred for tabular sweeps).
- JSON keys must follow `<name>_<unit>` convention. Include `"schema_version"` field.
- All values are SI. No unit conversion inside analysis modules.
- Load via `io/loader.py` — no inline `open()` calls in analysis modules.

## Type Annotations

All public functions must have full type annotations. Use `from __future__ import annotations` at the top of each file.

## Docstrings

Use **NumPy-style** docstrings for all public functions and classes.

## Error Handling

- Raise descriptive exceptions at system boundaries (file I/O, user input).
- No bare `except:` or `except Exception:` — catch specific types.
- Analysis modules raise `ValueError` for invalid physical inputs; TUI catches and displays them.

## Code Style

- Follow PEP 8. Max line length: 100 characters.
- Use `black` for formatting, `ruff` for linting.
- No `import *`. Imports ordered: stdlib → third-party → local.

## Testing

- Use `pytest`. One test file per module.
- Tests use only hardcoded SI values — no file I/O in unit tests.
- Numerical results checked with `np.testing.assert_allclose(actual, expected, rtol=1e-4)`.

## Output / Results

- Save to `data/outputs/` as JSON or CSV. Never overwrite inputs.
- Filename: `<method>_<timestamp_iso8601>.json`, e.g. `pk_flutter_20260419T142301.json`.
- Plots saved as PNG at 150 dpi minimum alongside result file.
- Log run parameters under a `"run_metadata"` key in the output JSON.

## Pre-commit Checklist

- [ ] All public function parameters and return values carry `_<unit>` suffixes.
- [ ] All values are SI with no in-module unit conversion.
- [ ] `main.py` contains only entry-point code.
- [ ] Input loaded via `io/loader.py`, not inline `open()` calls.
- [ ] `black` and `ruff` pass with no errors.
- [ ] New functions have NumPy-style docstrings.
- [ ] At least one `pytest` test covers each new public function.
- [ ] Unused packages removed from `requirements.txt`.
