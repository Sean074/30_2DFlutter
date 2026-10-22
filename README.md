# 2D Flutter Analysis (p-k method)

Flutter onset speed of a two-degree-of-freedom (heave + pitch) airfoil section in
incompressible flow. The tool sweeps a velocity range, solves the p-k flutter equation at
each speed for every structural mode, and reports where the damping of a mode crosses zero.

It is a 1990s-style terminal program: plain ANSI text, no curses, four fixed stages
(INPUT, CHECK INPUT, OUTPUT, CHECK OUTPUT).

## Install

Python 3.12 is used for development.

```bash
pip install -r requirements.txt
```

## Run

Interactive (prompts for every input, with defaults; `q` quits at any prompt):

```bash
python main.py
```

Batch, from a JSON input file:

```bash
python main.py --batch data/inputs/default.json
python main.py --batch data/inputs/default.json --no-save   # print only, write nothing
```

Both modes print a V-g / V-f table, a summary with the flutter speed and frequency of each
mode that flutters, and show the V-g / V-f plot.

## Input file

JSON, SI units, keys named `<name>_<unit>`. Example (`data/inputs/default.json`):

```json
{
  "schema_version": 1,
  "b_m": 0.5,
  "a_h_nd": 0.25,
  "x_alpha_nd": -0.5,
  "mass_kg": 10.4,
  "inertia_alpha_kg_m2": 1.12,
  "omega_alpha_rad_s": 21.0,
  "omega_heave_rad_s": 20.0,
  "rho_kg_m3": 1.21,
  "velocity_range_m_s": [0.01, 50.0, 4990]
}
```

| Key | Meaning |
|---|---|
| `b_m` | Semi-chord b [m] |
| `a_h_nd` | Elastic axis position aft of mid-chord, in semi-chords (`a`) |
| `x_alpha_nd` | CG position aft of the elastic axis, in semi-chords (`x_a`) |
| `mass_kg` | Section mass [kg] |
| `inertia_alpha_kg_m2` | Pitch inertia about the elastic axis [kg m2] |
| `omega_heave_rad_s`, `omega_alpha_rad_s` | Uncoupled heave and pitch frequencies [rad/s] |
| `rho_kg_m3` | Air density [kg/m3] |
| `velocity_range_m_s` | `[v_min, v_max, n_points]` of the sweep |

The heave and pitch stiffnesses are derived as `omega_heave^2 * mass` and
`omega_alpha^2 * inertia`.

This example has the CG ahead of the elastic axis (`x_alpha_nd < 0`), which is stable: it
correctly reports no flutter up to 50 m/s. The interactive defaults (`a = -0.2`,
`x_a = 0.1`, `omega_h = 10`, `omega_a = 25` rad/s, 10 kg, 0.5 kg m2) flutter at 19.39 m/s.

## Output

With saving enabled, each run writes to `data/outputs/` (git-ignored):

- `pk_flutter_<UTC timestamp>.json`: `run_metadata`, the velocity sweep, and per mode the
  real and imaginary parts of p at every speed.
- `pk_flutter_<UTC timestamp>.png`: V-f (top) and V-g (bottom) plot at 150 dpi, with a
  `g = 0.03` limit line.

Modes are labelled by their coupled in-vacuo frequency, e.g. `MODE 1 (15.43 rad/s)`.

## Method

For each mode and speed V the p-k equation

```
[ (V/b)^2 p^2 M + K - rho pi b V^2 A(k) ] x = 0,      k = imag(p)
```

is solved by fixed-point iteration on k: with k fixed it is a 2x2 eigenproblem in p; the
root nearest the mode's current estimate is kept and k is updated until it converges.
`A(k)` is Theodorsen's oscillatory aerodynamics with the Jones two-pole approximation of
`C(k)`. Modes are seeded from the coupled structural frequencies and extrapolated from speed
to speed so each mode keeps its identity.

Damping is `g = real(p) / imag(p)` and frequency is `imag(p) V / b`. Flutter onset is the
first speed where g goes from negative to non-negative; the speed and frequency are
interpolated linearly in g between the two bracketing sweep points, so the result does not
snap to the grid.

The solver raises an error rather than return a wrong answer if the iteration does not
converge, if two modes collapse onto one root (refine the sweep), or if a mode's frequency
reaches zero (static divergence; lower the maximum speed). A warning is printed when the
reduced frequency drops below 0.01, where the quasi-steady limit is approached.

### Verification

The textbook section of Hodges & Pierce, *Introduction to Structural Dynamics and
Aeroelasticity*, ch. 5 (`a = -1/5`, `x_a = 1/10`, `mu = 20`, `r^2 = 6/25`,
`omega_h / omega_a = 2/5`) is a regression test: `V_F / (b omega_a) = 2.1536`,
`omega_F / omega_a = 0.650`. An independent k-method solve gives the same values; with the
exact Theodorsen function it gives 2.184, so the Jones approximation costs about 1.4 %.

## Tests

```bash
python -m pytest
```

## Build a single-file executable

```bash
bash build.sh     # macOS  -> dist/flutter2d
build.bat         # Windows -> dist\flutter2d.exe
```

Both run PyInstaller on `build.spec`. Build on each platform for its own binary.

## Layout

```
main.py              entry point (calls ui.tui.run_tui)
analysis/flutter.py  matrix builders, p-k solver, flutter-speed interpolation
file_io/             loader.py (JSON / CSV in), exporter.py (results out)
plotting/plots.py    V-g / V-f figure
ui/tui.py            four-stage terminal workflow and batch mode
defaults/            default values and constants (JSON)
data/inputs/         example input files
tests/               pytest, mirrors the source tree
docs/                engineering code standard and TUI style guide
```

Code conventions are in `docs/engineering_code_standard.md` and `docs/user_interface.md`.
