# p-k Flutter Analysis of a Typical Section: Theory and Solver

This document explains the physics and numerics behind `analysis/flutter.py` (p-k) and
`analysis/k_method.py` (k-method, selectable since issue #10), and records the change of p-k
solver made for issue #2 (commit `d147b89`): what the old secant solver did, why it failed,
what replaced it, and how the result was verified.

Figures are regenerated with `python docs/figures/make_figures.py`, which also prints the
verification table in section 9.

---

## 1. Purpose and scope

The tool finds the **flutter speed** of a two-degree-of-freedom airfoil section (the
"typical section"): the lowest airspeed at which one of its aeroelastic modes stops being
damped and starts to grow. It sweeps a velocity range and, at each speed, solves the p-k
flutter equation for every structural mode.

Assumptions:

- 2D section, rigid in its own plane, on a heave spring and a torsion spring.
- Small motions about a zero-lift equilibrium (linear structure and aerodynamics).
- Incompressible, inviscid, attached flow over a thin airfoil (Theodorsen's theory).
- No structural damping in the model; damping is the *output* (section 4).

Outside this envelope (high subsonic Mach number, stall, large amplitude, 3D effects) the
results should not be used without corrections.

---

## 2. Structural model

![Typical section geometry](figures/section_geometry.png)

All lengths are measured along the chord in units of the semi-chord `b`, positive aft,
from mid-chord:

| Symbol | Meaning | Code |
|---|---|---|
| `b` | semi-chord (chord `c = 2b`) | `b_m` |
| `a` | elastic axis (EA) position aft of mid-chord, in semi-chords | `a_h_nd` |
| `x_α` | CG position aft of the EA, in semi-chords | `x_alpha_nd` |
| `m` | section mass | `mass_kg` |
| `I_α` | pitch inertia about the EA | `inertia_alpha_kg_m2` |
| `ω_h`, `ω_α` | uncoupled heave and pitch frequencies | `omega_heave_rad_s`, `omega_alpha_rad_s` |

The degrees of freedom are the heave `h` of the EA (**positive down**) and the pitch `α`
(**positive nose up**). Lift `L` is positive up and the moment `M_EA` about the EA is
positive nose up, following Theodorsen.

The equations of motion are

$$
m\ddot h + S_\alpha\ddot\alpha + k_h h = -L, \qquad
S_\alpha\ddot h + I_\alpha\ddot\alpha + k_\alpha\alpha = M_{EA},
\qquad S_\alpha = m\,b\,x_\alpha ,
$$

so the mass and stiffness matrices (`build_mass_matrix`, `build_stiffness_matrix`) are

$$
\mathbf M = \begin{bmatrix} m & m b x_\alpha \\ m b x_\alpha & I_\alpha \end{bmatrix},
\qquad
\mathbf K = \begin{bmatrix} k_h & 0 \\ 0 & k_\alpha \end{bmatrix},
\qquad k_h = m\omega_h^2,\; k_\alpha = I_\alpha\omega_\alpha^2 .
$$

The usual dimensionless groups are the mass ratio `μ = m/(πρb²)`, the radius of gyration
`r_α² = I_α/(m b²)` and the frequency ratio `σ = ω_h/ω_α`. `M` is positive definite only if
`r_α² > x_α²`; `compute_structural_frequencies` raises `ValueError` otherwise.

### Uncoupled vs coupled frequencies

`ω_h` and `ω_α` describe each spring with the other degree of freedom locked. When the CG is
off the EA (`x_α ≠ 0`) the off-diagonal mass term couples the two, and the actual in-vacuo
natural frequencies are the roots of `det(K − ω²M) = 0`, returned by
`compute_structural_frequencies`.

For `data/inputs/default.json` (`m = 10.4 kg`, `I_α = 1.12 kg m²`, `b = 0.5 m`,
`x_α = −0.5`) the coupling is strong: `m b² x_α² / I_α = 0.58`. The uncoupled frequencies of
20 and 21 rad/s become **15.43 and 42.02 rad/s** once coupled. That gap is the root cause of
issue #2 (section 5).

---

## 3. Unsteady aerodynamics

### Theodorsen's lift and moment

For harmonic motion of reduced frequency `k = bω/V`, Theodorsen's lift and moment about the
EA are

$$
L = \pi\rho b^2\left(\ddot h + V\dot\alpha - b a\ddot\alpha\right)
  + 2\pi\rho V b\, C(k)\left[\dot h + V\alpha + b\left(\tfrac12 - a\right)\dot\alpha\right]
$$

$$
M_{EA} = \pi\rho b^2\left[b a\ddot h - V b\left(\tfrac12 - a\right)\dot\alpha
  - b^2\left(\tfrac18 + a^2\right)\ddot\alpha\right]
  + 2\pi\rho V b^2\left(a + \tfrac12\right) C(k)\left[\dot h + V\alpha
  + b\left(\tfrac12 - a\right)\dot\alpha\right]
$$

The first bracket of each is non-circulatory (apparent mass and damping). The second is
circulatory: the lift of the quasi-steady downwash at the three-quarter chord,
`ḣ + Vα + b(½ − a)α̇`, acting at the quarter chord and reduced and lagged by the wake through
the Theodorsen function `C(k)`.

### The aerodynamic matrix

With `h = h̄ e^{iωt}`, `α = ᾱ e^{iωt}`, the loads are linear in `x̄ = [h̄, ᾱ]ᵀ`.
`build_aero_matrix` returns the 2×2 complex matrix `A(k)` defined by

$$
\rho\pi b V^2\,\mathbf A(k)\,\bar{\mathbf x} = \begin{bmatrix} -\bar L \\ \bar M_{EA} \end{bmatrix},
$$

so that, for example, `A₁₁ = k²/b − 2iC(k)k/b` is the heave-on-heave term of `−L` divided by
`ρπbV²`. The minus sign on `L` comes from `h` being positive down while `L` is positive up.

### The Theodorsen function and the Jones approximation

$$
C(k) = F(k) + iG(k) = \frac{H_1^{(2)}(k)}{H_1^{(2)}(k) + i H_0^{(2)}(k)},
$$

with `C → 1` as `k → 0` (quasi-steady) and `C → ½` as `k → ∞`. The code (`theo_c`) uses
R. T. Jones' two-pole rational approximation (NACA Report 681, from his fit to the Wagner
function):

$$
C(k) \approx 1 - \frac{0.165}{1 - 0.0455\,i/k} - \frac{0.335}{1 - 0.300\,i/k}.
$$

![Theodorsen function, exact vs Jones](figures/theodorsen_jones_vs_exact.png)

The largest error against the exact function is `|ΔC| = 0.0145` for `k ≤ 2`. Until
issue #9 the code used `0.041` and `0.32` in place of `0.0455` and `0.300`; that fit is worse
(largest `|ΔC| = 0.021`, near `k = 0.18`) and made the flutter speed noticeably low:

| Case | Old `0.041 / 0.32` | R. T. Jones (code) | Exact Theodorsen |
|---|---|---|---|
| Hodges & Pierce, `V_F/(bω_α)` | 2.1536 | 2.1704 | 2.1839 |
| Shipped defaults, `V_F` [m/s] | 19.393 | 19.558 | 19.575 |

So most of the gap between a two-pole fit and exact theory comes from the choice of
coefficients, not from the two-pole form itself.

---

## 4. The p-k flutter equation

### Derivation

Substituting the aerodynamic loads into the equations of motion and writing the motion as
`x̄ e^{st}` with complex `s = σ + iω` gives

$$
\left[s^2\mathbf M + \mathbf K - \rho\pi b V^2\,\mathbf A(k)\right]\bar{\mathbf x} = 0 .
$$

With the non-dimensional eigenvalue `p = s b / V` this is the form solved in the code:

$$
\mathbf F(p)\,\bar{\mathbf x} =
\left[\left(\tfrac{V}{b}\right)^2 p^2\,\mathbf M + \mathbf K - \rho\pi b V^2\,\mathbf A(k)\right]
\bar{\mathbf x} = 0,
\qquad k = \operatorname{Im}(p).
$$

The **p-k approximation** is the condition `k = Im(p)`: the aerodynamics are those of
harmonic motion at the mode's own frequency, even though the mode is actually growing or
decaying. This is exact only when `Re(p) = 0`, which is exactly the flutter point.

### Meaning of p, and the damping convention

| Quantity | Expression | Code |
|---|---|---|
| reduced frequency | `k = Im(p) = bω/V` | `np.imag(p)` |
| frequency [rad/s] | `ω = Im(p) V / b` | plotted on the V-f axis |
| decay-rate ratio | `γ = Re(p)/Im(p) = σ/ω` | |
| **damping** | **`g = 2γ = 2 Re(p)/Im(p)`** | `damping_g_nd` |

`g` is the flutter-analysis convention: the equivalent **structural damping** a mode has (if
negative) or would need (if positive) to be neutrally stable. For a lightly damped mode
`g ≈ −2ζ`, where `ζ` is the viscous damping ratio. A structure is commonly assumed to have
about `g = 0.03` (`γ = 1.5 %`) of structural damping, which is why the V-g plot carries a
`g = 0.03` line; the speed at which a mode crosses it is the flutter speed with that damping
allowed for.

### Flutter onset

Flutter onset is the first speed at which a mode's `g` goes from negative to non-negative.
`find_flutter_speed` finds the bracketing pair of sweep speeds and interpolates both speed
and frequency linearly in `g`, so the answer does not snap to the velocity grid.

### Relation to the k-method

The older **k-method** adds artificial structural damping, `K(1 + ig)`, and solves the
eigenproblem at fixed `k` for `Z = (1 + ig)/ω²`. Its `g` is the damping *required* for
harmonic motion. It has no physical meaning away from `g = 0` and its modes are hard to
track, but **at `g = 0` it solves exactly the same equation as the p-k method**. That makes
it a clean independent check of the flutter point (section 9), and it is available as an
alternative solver (section 7).

---

## 5. Old solver: secant iteration on det F (removed in `d147b89`)

### Algorithm

The old `compute_flutter` solved one mode per call, from one starting frequency `ω₀`:

1. At the first speed, `k₀ = ω₀ b / V₁`, with two starting points
   `p₁ = k₀(d₁ + i)` and `p₂ = k₀(d₂ + i)` (`d₁ = 0`, `d₂ = 0.01`).
2. Secant step on the scalar determinant, each `F` evaluated with `k = Im(p)`:
   `p₃ = (p₂ D₁ − p₁ D₂)/(D₁ − D₂)`, where `D_i = det F(p_i)`.
3. Stop when `|det F(p₃)| ≤ 0.001` (`det_tol_nd`).
4. Seed the next speed by scaling both points by `V/(V + ΔV)`, which keeps the physical
   eigenvalue `s = pV/b` fixed.

The TUI called it twice, with `ω₀ = ω_h` and `ω₀ = ω_α` (the *uncoupled* frequencies), and
labelled the results HEAVE and PITCH.

This is close to Hassig's original p-k formulation, which also iterated on the
determinant. The method is not wrong in principle; the failures were in seeding and
bookkeeping.

### Why it failed

- **Root selection.** A secant search on `det F(p)` converges to whichever root's basin the
  starting points are in. Nothing tied a call to a particular mode.
- **Seeding from uncoupled frequencies.** For `default.json` both starting frequencies
  (20 and 21 rad/s) are far from the upper coupled mode (42 rad/s) and both lie in the basin
  of the lower one (15.4 rad/s). Both calls returned the **same** mode, so the second mode
  was never solved and could never show flutter.
- **Tolerance with units.** `det F` has units of `k_h k_α` (N²/rad). The absolute tolerance
  `0.001` was therefore a relative tolerance of about `5 × 10⁻¹⁰` for `default.json`
  (`det K ≈ 2 × 10⁶`) but about `7 × 10⁻⁶` for the Hodges & Pierce section
  (`det K ≈ 150`). Its strictness depended on the units and size of the model.
- **No iteration cap.** If the secant step did not converge, the loop never ended.
- **No mode check.** Nothing detected two calls returning the same root, or a call jumping
  to another mode between speeds. The step `ΔV` was also taken from the first two points,
  which assumes a uniform grid.

![Old vs new solver on default.json](figures/old_vs_new_default_json.png)

In the left column of the figure the two old calls lie exactly on top of each other. In the
right column the new solver finds both modes (neither flutters for this input, which is
correct: the CG is ahead of the EA).

---

## 6. New solver: eigenvalue tracking

### The k iteration

For fixed `k` the p-k equation is a quadratic eigenproblem in `p` with no first-order term:

$$
\left(\tfrac{V}{b}\right)^2 p^2\,\mathbf M\,\bar{\mathbf x} = -\mathbf Q(k)\,\bar{\mathbf x},
\qquad \mathbf Q(k) = \mathbf K - \rho\pi b V^2 \mathbf A(k)
\;\;\Rightarrow\;\;
p^2 = -\left(\tfrac{b}{V}\right)^2 \operatorname{eig}\!\left(\mathbf M^{-1}\mathbf Q(k)\right).
$$

Each of the two eigenvalues `p²` gives a pair `±p`; the root with `Im(p) ≥ 0` is kept (the
other is its mirror). For each mode, `_solve_mode_at_speed`:

1. takes the current estimate `p`, sets `k = Im(p)` and computes the two candidate roots;
2. keeps the candidate **nearest the current estimate**;
3. updates `k` and repeats until `|Δk| ≤ k_tol · |k|` (`k_tol_nd`, default `10⁻¹⁰`).

### Seeding and mode identity

- At the first speed, each mode starts from its **coupled** structural frequency,
  `s = iω_i` (`compute_structural_frequencies`), so every mode has its own starting point.
- At the second speed the start is the first converged `s = pV/b`; from the third speed on it
  is extrapolated linearly in `V` from the previous two. `s` is the physical eigenvalue and varies smoothly with speed,
  whereas `p` grows like `1/V` at low speed.
- All modes are solved together at each speed. After each speed, if any two modes are
  closer than `10⁻⁶ max|s|`, the solver raises `ValueError` rather than return duplicates.

`compute_flutter` returns an array of shape `(n_modes, nv)`, modes ordered by ascending
structural frequency.

### Failure behaviour

The solver raises `ValueError` instead of returning a wrong answer when:

- `M` or `K` is not positive definite;
- the velocity sweep is empty, non-positive or not strictly increasing;
- the k iteration does not converge in `max_iter_nd` steps (default 200);
- `k` falls to zero, which is static divergence (the remedy is a lower maximum speed);
- two modes collapse onto one root (the remedy is more velocity points).

---

## 7. k-method solver (`analysis/k_method.py`)

Selected with `"method": "k"` in the input file or `k` at the METHOD prompt; `"pk"` (the
default) selects section 6.

### Equation

With harmonic motion, `p = ik`, and artificial structural damping `K(1 + ig)`, the flutter
equation of section 4 becomes, after using `ρπbV² = ρπb³ω²/k²`,

$$
\left[ K^{-1}\left( M + \frac{\rho\pi b^3}{k^2} A(k) \right) \right] x = Z x,
\qquad Z = \frac{1 + ig}{\omega^2}.
$$

For each `k` this is a plain eigenproblem with no iteration. Every eigenvalue gives

$$
\omega = \frac{1}{\sqrt{\operatorname{Re} Z}}, \qquad V = \frac{b\omega}{k}, \qquad
g = \frac{\operatorname{Im} Z}{\operatorname{Re} Z}.
$$

### Sweep and branches

`compute_k_branches` sweeps `k` downwards (speed upwards) and keeps each branch by matching the
eigenvalues at one `k` to the nearest ones at the previous `k`. Branches start in ascending
frequency, so mode numbers match the p-k solver. `compute_flutter_k` chooses the `k` range from
the velocity range and the coupled frequencies (from `2bω_max/V_min` down to
`0.25bω_min/V_max`, `n_k_nd` = 4000 points, geometric). It then interpolates each branch's `ω`
and `g` linearly in `V` onto the requested velocity grid and returns them as an equivalent
`p = k(g/2 + i)`. That has the shape and meaning of the p-k output, so the table, plot,
exporter and `find_flutter_speed` are shared.

`V(k)` is not always monotonic. Approaching static divergence, the lower branch's speed peaks
and then falls again as `k → 0`. A branch is used from its lowest speed up to that peak; at
higher speeds the result is NaN (`-` in the table, `null` in the JSON, a gap in the plot), and
CHECK OUTPUT warns where the branch ends.

### What the k-method g means

The k-method `g` is the structural damping that would be *needed* for neutral harmonic
motion at that speed. It is not the damping the section actually has. Only at `g = 0` is the
motion truly harmonic, and there the equation is the p-k equation, so the two methods give
the same flutter speed and frequency. Away from `g = 0` the V-g curves differ, and in the
k-method the branches can behave differently from the p-k roots (here the lower branch folds
back instead of heading towards divergence):

![p-k vs k-method, shipped defaults](figures/pk_vs_k_shipped_defaults.png)

Use the k-method to cross-check a flutter point. Read damping margins (the `g = 0.03` line)
from the p-k method.

---

## 8. Comparison

| | Old solver | New solver |
|---|---|---|
| Equation | `det F(p) = 0`, `k = Im(p)` | `F(p) x̄ = 0`, `k = Im(p)` |
| Iteration | secant on the scalar determinant | fixed point on `k`, 2×2 eigenproblem per step |
| Root selection | whichever root the secant reaches | nearest to the tracked estimate |
| Start, first speed | uncoupled `ω_h` or `ω_α` | coupled structural frequencies |
| Start, next speed | previous `p` × `V/(V+ΔV)` | linear extrapolation of `s = pV/b` |
| Modes per call | one | all, together |
| Convergence test | `\|det F\| ≤ 0.001` (absolute, has units) | `\|Δk\| ≤ 10⁻¹⁰ \|k\|` (relative) |
| Iteration cap | none | 200 |
| Duplicate modes | not detected | `ValueError` |
| Output | shape `(nv,)` per call | shape `(n_modes, nv)` |

Both solve the same equation, so wherever the old solver happened to land on the right
root the two agree. The difference is robustness and keeping track of which mode is which.

---

## 9. Verification

An independent **k-method** solver, written directly from Theodorsen's `L` and `M_EA`
(section 3) and sharing no code with `analysis/flutter.py`, is in
`docs/figures/make_figures.py`. At `g = 0` it must agree with the p-k solver exactly
(section 4). It is kept separate from the program's own k-method (section 7), which reuses
`build_aero_matrix` and so would not catch an error there; `tests/test_k_method.py` checks
that one against the same flutter points.

| Case | p-k `V_F` | p-k `ω_F` | k-method `V_F` | k-method `ω_F` | Exact `C(k)` `V_F` |
|---|---|---|---|---|---|
| Hodges & Pierce (`V_F/bω_α`, `ω_F/ω_α`) | 2.17036 | 0.64433 | 2.17036 | 0.64433 | 2.18389 |
| Shipped defaults (m/s, rad/s) | 19.55820 | 16.01076 | 19.55813 | 16.01074 | 19.57463 |
| `data/inputs/default.json` | none | | none | | none |

The Hodges & Pierce section (*Introduction to Structural Dynamics and Aeroelasticity*,
ch. 5) has `a = −1/5`, `x_α = 1/10`, `μ = 20`, `r_α² = 6/25`, `σ = 2/5`. The shipped
defaults are those of `defaults/default_flutter.json`. p-k and k-method agree to 5–6
significant figures. With the exact Theodorsen function the flutter speed is 0.6 % (Hodges &
Pierce) and 0.1 % (shipped defaults) higher, which is the cost of the Jones approximation
(section 3).

![V-g / V-f, Hodges & Pierce](figures/vg_hodges_pierce.png)

![V-g / V-f, shipped defaults](figures/vg_shipped_defaults.png)

Both cases show classic bending-torsion (coalescence) flutter: as speed rises the two
frequencies approach, and the upper mode's damping turns positive. The lower mode becomes
heavily damped and heads towards static divergence.

![V-g / V-f, default.json](figures/vg_default_json.png)

With the CG ahead of the EA (`x_α < 0`) the section does not flutter in the range.

### Tests

| Claim | Test (`tests/test_flutter.py`) |
|---|---|
| Hodges & Pierce flutter speed and frequency | `test_compute_flutter_hodges_pierce_textbook_case` |
| Shipped-defaults flutter speed | `test_compute_flutter_detects_flutter_speed` |
| Modes stay distinct (issue #2) | `test_compute_flutter_keeps_modes_distinct` |
| Coupled frequencies | `test_compute_structural_frequencies_uncoupled` |
| Interpolated onset | `test_find_flutter_speed_interpolates_between_grid_points` |
| `g = 2γ` | `test_damping_g_nd_is_twice_decay_rate_ratio` |
| `C(k) → ½` | `test_theo_c_approaches_half_at_high_k` |
| k-method flutter point = p-k (both cases) | `tests/test_k_method.py` |
| k-method branch ends at its fold | `test_k_method_branch_is_nan_beyond_its_fold` |

---

## 10. Limitations and possible extensions

- **Jones approximation.** The two-pole `C(k)` makes the flutter speed 0.1–0.6 % low in
  the verified cases. The exact Theodorsen function (Hankel functions) would remove this
  but needs Bessel functions the code does not otherwise depend on.
- **Incompressible, thin-airfoil aerodynamics.** No Mach effects; no thickness, camber or
  viscous effects; no stall.
- **No structural damping input.** `g` is an output; a damping margin is read off the plot
  against the `g = 0.03` line rather than included in the model.
- **p-k damping away from `g = 0` is approximate**, because the aerodynamics assume
  harmonic motion.
- **Mode swaps are not detected.** The collapse check catches two modes on one root, but if
  a coarse grid near frequency coalescence lets two modes exchange identities, it is not
  reported. More velocity points avoid it.
- **k-method branches can fold back** (`V(k)` not monotonic); speeds beyond the fold are
  not reported for that branch (section 7).
- **Divergence is not computed.** The solver stops with an error when `k → 0`; it does not
  report the divergence speed.
- **Extensions:** a control-surface (third) degree of freedom, a structural damping input,
  compressibility corrections, a divergence-speed calculation.

---

## 11. Nomenclature and references

| Symbol | Meaning | Code |
|---|---|---|
| `a` | EA position aft of mid-chord / `b` | `a_h_nd` |
| `A(k)` | aerodynamic matrix | `build_aero_matrix` |
| `b` | semi-chord | `b_m` |
| `C(k) = F + iG` | Theodorsen function | `theo_c` |
| `g` | damping, `2 Re(p)/Im(p)` | `damping_g_nd` |
| `γ` | decay-rate ratio, `Re(p)/Im(p)` | |
| `h` | heave of the EA, positive down | |
| `I_α` | pitch inertia about the EA | `inertia_alpha_kg_m2` |
| `k` | reduced frequency `bω/V = Im(p)` | |
| `k_h`, `k_α` | heave and pitch stiffness | `k_heave_n_m`, `k_alpha_nm_rad` |
| `L`, `M_EA` | lift (up), moment about the EA (nose up) | |
| `m` | section mass | `mass_kg` |
| `M`, `K` | mass and stiffness matrices | `build_mass_matrix`, `build_stiffness_matrix` |
| `p` | non-dimensional eigenvalue `s b / V` | rows of `compute_flutter` output |
| `r_α²` | `I_α/(m b²)` | |
| `s = σ + iω` | physical eigenvalue | |
| `V`, `V_F` | airspeed, flutter speed | `velocity_range_m_s`, `find_flutter_speed` |
| `x_α` | CG position aft of the EA / `b` | `x_alpha_nd` |
| `α` | pitch, positive nose up | |
| `μ` | mass ratio `m/(πρb²)` | |
| `ρ` | air density | `rho_kg_m3` |
| `σ` | frequency ratio `ω_h/ω_α` (section 2); real part of `s` (section 4) | |
| `ω_h`, `ω_α` | uncoupled heave and pitch frequencies | `omega_heave_rad_s`, `omega_alpha_rad_s` |

1. T. Theodorsen, *General Theory of Aerodynamic Instability and the Mechanism of
   Flutter*, NACA Report 496, 1935.
2. R. T. Jones, *The Unsteady Lift of a Wing of Finite Aspect Ratio*, NACA Report 681, 1940.
3. H. J. Hassig, "An Approximate True Damping Solution of the Flutter Equation by
   Determinant Iteration", *Journal of Aircraft* 8(11), 885–889, 1971.
4. D. H. Hodges and G. A. Pierce, *Introduction to Structural Dynamics and Aeroelasticity*,
   2nd ed., Cambridge University Press, 2011, ch. 5.
5. R. L. Bisplinghoff, H. Ashley and R. L. Halfman, *Aeroelasticity*, Addison-Wesley, 1955
   (Dover reprint, 1996).
