# 2D Flutter Analysis (p-k and k methods)

Flutter onset speed of a two-degree-of-freedom (heave + pitch) airfoil section in
incompressible flow, optionally with a trailing-edge flap as a third degree of freedom. The tool sweeps a velocity range, solves the flutter equation at each
speed for every structural mode with the p-k method (default) or the k-method, and reports
where the damping of a mode crosses zero.

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
  "method": "pk",
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
| `method` | Solver: `"pk"` (p-k, default if omitted) or `"k"` (k-method) |
| `b_m` | Semi-chord b [m] |
| `a_h_nd` | Elastic axis position aft of mid-chord, in semi-chords (`a`) |
| `x_alpha_nd` | CG position aft of the elastic axis, in semi-chords (`x_a`) |
| `mass_kg` | Section mass [kg] |
| `inertia_alpha_kg_m2` | Pitch inertia about the elastic axis [kg m2] |
| `omega_heave_rad_s`, `omega_alpha_rad_s` | Uncoupled heave and pitch frequencies [rad/s] |
| `rho_kg_m3` | Air density [kg/m3] |
| `velocity_range_m_s` | `[v_min, v_max, n_points]` of the sweep |
| `c_hinge_nd` | Flap hinge position aft of mid-chord, in semi-chords (`c`, between -1 and 1) |
| `x_beta_nd` | Flap CG aft of the hinge, in semi-chords, normalised by the total mass (`x_b`) |
| `inertia_beta_kg_m2` | Flap inertia about the hinge [kg m2] |
| `omega_beta_rad_s` | Uncoupled flap (hinge-spring) frequency [rad/s] |

The heave and pitch stiffnesses are derived as `omega_heave^2 * mass` and
`omega_alpha^2 * inertia`.

The four flap keys are optional: give all of them for the three-DOF section (heave, pitch,
flap) or none for two DOF; a partial set is an error. The flap stiffness is
`omega_beta^2 * inertia_beta`, the flap static moment is `mass * b * x_beta`, and
`inertia_alpha_kg_m2` is for the whole section, flap included.
`data/inputs/flap_example.json` is the shipped-defaults section with a flap; it flutters at
20.50 m/s. With three modes the results table uses narrower gaps so it stays within 80
columns.

This example has the CG ahead of the elastic axis (`x_alpha_nd < 0`), which is stable: it
correctly reports no flutter up to 50 m/s.

## Defaults

The interactive prompts are pre-filled from `defaults/default_flutter.json` (same keys as an
input file, plus the solver settings `k_tol_nd` and `max_iter_nd` (p-k) and `n_k_nd`
(number of reduced frequencies in the k-method sweep), and `flap` (whether the INCLUDE FLAP
prompt defaults to yes)); `pi_nd` comes from
`defaults/constants.json`. Edit those files to change the defaults; both are bundled into the
executable. The shipped defaults (`a = -0.2`, `x_a = 0.1`, `omega_h = 10`, `omega_a = 25`
rad/s, 10 kg, 0.5 kg m2, 1-100 m/s, p-k method) flutter at 19.56 m/s.

## Output

With saving enabled, each run writes to `data/outputs/` (git-ignored):

- `<method>_flutter_<UTC timestamp>.json` (`pk_` or `k_`): `run_metadata` (including the
  method), the velocity sweep, and per mode the real and imaginary parts of p at every speed
  (`null` where a k-method branch does not reach that speed).
- `<method>_flutter_<UTC timestamp>.png`: V-f (top) and V-g (bottom) plot at 150 dpi, with
  a `g = 0.03` limit line.

Modes are labelled by their coupled in-vacuo frequency, e.g. `MODE 1 (15.43 rad/s)`.

## Method

For each mode and speed V the p-k equation

```
[ (V/b)^2 p^2 M + K - rho pi b V^2 A(k) ] x = 0,      k = imag(p)
```

is solved by fixed-point iteration on k: with k fixed it is a 2x2 eigenproblem in p; the
root nearest the mode's current estimate is kept and k is updated until it converges.
`A(k)` is Theodorsen's oscillatory aerodynamics with R. T. Jones' two-pole approximation of
`C(k)` (NACA Report 681: poles at `k = 0.0455` and `0.300`). Modes are seeded from the
coupled structural frequencies and extrapolated from speed to speed so each mode keeps its
identity.

Damping is the equivalent structural damping `g = 2 real(p) / imag(p)` (twice the decay-rate
ratio `gamma`; `g = 0.03` is `gamma = 1.5 %`) and frequency is `imag(p) V / b`. Flutter onset is the
first speed where g goes from negative to non-negative; the speed and frequency are
interpolated linearly in g between the two bracketing sweep points, so the result does not
snap to the grid.

### k-method

With `"method": "k"` the classical k-method is used instead (`analysis/k_method.py`). It
adds artificial structural damping `K(1 + ig)`, sweeps the reduced frequency k and solves one
eigenproblem per k for `Z = (1 + ig) / omega^2`, giving `omega`, `V = b omega / k` and the
*required* damping g of each branch. The branches are interpolated onto the velocity sweep, so
the table, plot and output look the same as for p-k. The flutter point (g = 0) is the same as
the p-k one; elsewhere the k-method g is not physical damping and the V-g curves differ, so
read damping margins from the p-k method. Where a branch's speed folds back (near static
divergence) the higher speeds are left empty and a warning is printed. Details: section 7 of
the theory document.

### Errors and warnings

The p-k solver raises an error rather than return a wrong answer if two modes collapse onto
one root (refine the sweep). A mode the p-k iteration can no longer follow (its frequency
falls to zero, near static divergence or because it is so heavily damped it stops
oscillating, or the iteration does not converge) ends at that speed: the table shows `-`
beyond it, the other modes carry on, and a warning gives the speed. A warning is also
printed when the reduced frequency drops below 0.01, where the quasi-steady limit is
approached.

Full derivation, the old and new p-k solvers, the k-method, the flap degree of freedom, and
verification:
[docs/theory_pk_flutter.md](docs/theory_pk_flutter.md).

### Verification

The textbook section of Hodges & Pierce, *Introduction to Structural Dynamics and
Aeroelasticity*, ch. 5 (`a = -1/5`, `x_a = 1/10`, `mu = 20`, `r^2 = 6/25`,
`omega_h / omega_a = 2/5`) is a regression test: `V_F / (b omega_a) = 2.1704`,
`omega_F / omega_a = 0.644`. An independent k-method solve gives the same values; with the
exact Theodorsen function it gives 2.184, so the Jones approximation costs about 0.6 %.

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
analysis/k_method.py k-method solver
file_io/             loader.py (JSON / CSV in), exporter.py (results out)
plotting/plots.py    V-g / V-f figure
ui/tui.py            four-stage terminal workflow and batch mode
defaults/            default values and constants (JSON)
data/inputs/         example input files (default.json, flap_example.json)
tests/               pytest, mirrors the source tree
docs/                theory (theory_pk_flutter.md + figures/), code standard, TUI guide
```

Code conventions are in `docs/engineering_code_standard.md` and `docs/user_interface.md`.
