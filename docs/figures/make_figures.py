"""
Regenerate the figures and verification numbers for docs/theory_pk_flutter.md.

Run from the repository root:

    python docs/figures/make_figures.py

Writes PNGs (150 dpi) next to this script and prints the verification table. Needs only the
project's runtime dependencies (numpy, matplotlib).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402
from numpy.typing import NDArray  # noqa: E402

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from analysis.flutter import (  # noqa: E402
    build_aero_matrix,
    build_mass_matrix,
    build_stiffness_matrix,
    compute_flutter,
    compute_structural_frequencies,
    damping_g_nd,
    find_flutter_speed,
    theo_c,
)
from plotting.plots import plot_vg_vf  # noqa: E402

_OUT_DIR = Path(__file__).resolve().parent
_DPI = 150

# Test cases: (mass_kg, inertia_alpha_kg_m2, b_m, x_alpha_nd, a_h_nd, omega_h, omega_a, rho)
_MU_HP_ND = 20.0
_MASS_HP_KG = _MU_HP_ND * np.pi
CASES = {
    "hodges_pierce": (_MASS_HP_KG, 0.24 * _MASS_HP_KG, 1.0, 0.1, -0.2, 0.4, 1.0, 1.0),
    "shipped_defaults": (10.0, 0.5, 0.5, 0.1, -0.2, 10.0, 25.0, 1.225),
    "default_json": (10.4, 1.12, 0.5, -0.5, 0.25, 20.0, 21.0, 1.21),
}


# ---------------------------------------------------------------------------
# Exact Theodorsen function from Bessel functions (numpy-only quadrature)
# ---------------------------------------------------------------------------


def _bessel_j(n: int, x_nd: float) -> float:
    """J_n(x) = (1/pi) int_0^pi cos(n tau - x sin tau) d tau."""
    tau = np.linspace(0.0, np.pi, 4001)
    return float(np.trapezoid(np.cos(n * tau - x_nd * np.sin(tau)), tau) / np.pi)


def _bessel_y(n: int, x_nd: float) -> float:
    """Y_n(x) from its integral representation (valid for x > 0)."""
    tau = np.linspace(0.0, np.pi, 4001)
    part_1 = np.trapezoid(np.sin(x_nd * np.sin(tau) - n * tau), tau)
    t_max = np.arcsinh(60.0 / x_nd) + 1.0
    t = np.linspace(0.0, t_max, 200001)
    kernel = (np.exp(n * t) + (-1.0) ** n * np.exp(-n * t)) * np.exp(-x_nd * np.sinh(t))
    return float((part_1 - np.trapezoid(kernel, t)) / np.pi)


def theodorsen_exact(k_nd: float) -> complex:
    """C(k) = H1(2)(k) / (H1(2)(k) + i H0(2)(k)), H(2) = J - iY."""
    h0 = _bessel_j(0, k_nd) - 1j * _bessel_y(0, k_nd)
    h1 = _bessel_j(1, k_nd) - 1j * _bessel_y(1, k_nd)
    return h1 / (h1 + 1j * h0)


# ---------------------------------------------------------------------------
# Independent k-method, written from Theodorsen's L and M_ea (shares no code with the solver)
# ---------------------------------------------------------------------------


def _aero_force_matrix(omega, v, rho, b, a, c_k):
    """Harmonic aerodynamic force matrix: rows (+L in heave eq., -M_ea in pitch eq.)."""
    iw = 1j * omega
    lift_h = np.pi * rho * b**2 * (-(omega**2)) + 2 * np.pi * rho * v * b * c_k * iw
    lift_a = np.pi * rho * b**2 * (v * iw + b * a * omega**2) + 2 * np.pi * rho * v * b * c_k * (
        v + b * (0.5 - a) * iw
    )
    mom_h = np.pi * rho * b**2 * (-b * a * omega**2) + 2 * np.pi * rho * v * b**2 * (
        a + 0.5
    ) * c_k * iw
    mom_a = np.pi * rho * b**2 * (
        -v * b * (0.5 - a) * iw + b**2 * (1 / 8 + a**2) * omega**2
    ) + 2 * np.pi * rho * v * b**2 * (a + 0.5) * c_k * (v + b * (0.5 - a) * iw)
    return np.array([[lift_h, lift_a], [-mom_h, -mom_a]])


def k_method_flutter(case, c_func, k_arr_nd=None):
    """Return (V_F, omega_F) from the k-method g = 0 crossing, or None."""
    m, ia, b, xa, a, wh, wa, rho = case
    mass = np.array([[m, m * b * xa], [m * b * xa, ia]])
    stiff = np.diag([wh**2 * m, wa**2 * ia])
    if k_arr_nd is None:
        k_arr_nd = np.geomspace(3.0, 0.02, 3000)
    rows = []
    for k_nd in k_arr_nd:
        omega = 1.0
        v = b * omega / k_nd
        z_arr = np.linalg.eigvals(
            np.linalg.solve(stiff, mass - _aero_force_matrix(omega, v, rho, b, a, c_func(k_nd)))
        )
        row = []
        for z in z_arr:
            if z.real <= 0:
                row.append((np.nan, np.nan, np.nan))
                continue
            om = 1.0 / np.sqrt(z.real)
            row.append((b * om / k_nd, om, z.imag / z.real))
        rows.append(sorted(row, key=lambda r: -r[1] if np.isfinite(r[1]) else 0.0))
    out = np.array(rows)
    best = None
    for r in range(out.shape[1]):
        v_arr, om_arr, g_arr = out[:, r, 0], out[:, r, 1], out[:, r, 2]
        for j in range(len(k_arr_nd) - 1):
            if np.isfinite(g_arr[j]) and np.isfinite(g_arr[j + 1]) and g_arr[j] < 0 <= g_arr[j + 1]:
                t = -g_arr[j] / (g_arr[j + 1] - g_arr[j])
                cand = (
                    v_arr[j] + t * (v_arr[j + 1] - v_arr[j]),
                    om_arr[j] + t * (om_arr[j + 1] - om_arr[j]),
                )
                if best is None or cand[0] < best[0]:
                    best = cand
    return best


# ---------------------------------------------------------------------------
# Old secant solver, frozen copy of analysis/flutter.py at d147b89^ (for the comparison figure)
# ---------------------------------------------------------------------------


def old_secant_solver(mass, stiff, v_arr, rho, b, a, start_omega, det_tol=0.001, d1=0.0, d2=0.01):
    """The pre-#2 solver: secant on det F(p), seeded from one (uncoupled) frequency."""
    dv = v_arr[1] - v_arr[0]
    k_2 = start_omega * b / v_arr[0]
    k_1 = k_2
    p_1 = complex(d1 * k_2, k_2)
    p_2 = complex(d2 * k_2, k_2)
    out = np.empty(len(v_arr), dtype=np.complex128)

    def det_f(p, k, v):
        f = (v / b) ** 2 * p**2 * mass + stiff - rho * np.pi * b * v**2 * build_aero_matrix(
            k, theo_c(k), b, a
        )
        return np.linalg.det(f)

    for i, v in enumerate(v_arr):
        while True:
            det_1, det_2 = det_f(p_1, k_1, v), det_f(p_2, k_2, v)
            p_3 = (p_2 * det_1 - p_1 * det_2) / (det_1 - det_2)
            if abs(det_f(p_3, np.imag(p_3), v)) <= det_tol:
                break
            p_1, p_2 = p_2, p_3
            k_1, k_2 = np.imag(p_1), np.imag(p_2)
        out[i] = p_3
        p_1 = p_1 * v / (v + dv)
        p_2 = p_2 * v / (v + dv)
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def _matrices(case):
    m, ia, b, xa, a, wh, wa, rho = case
    return (
        build_mass_matrix(m, ia, b, xa),
        build_stiffness_matrix(wh**2 * m, wa**2 * ia),
        b,
        a,
        rho,
    )


def _naca0012(x_c: NDArray[np.float64]) -> NDArray[np.float64]:
    return 0.6 * (
        0.2969 * np.sqrt(x_c) - 0.1260 * x_c - 0.3516 * x_c**2 + 0.2843 * x_c**3 - 0.1015 * x_c**4
    )


def _dim(ax, x0, x1, y, text, ha="center", text_x=None, colour="k"):
    """Dimension arrow between x0 and x1 at height y, with its label just above."""
    ax.add_patch(
        FancyArrowPatch(
            (x0, y), (x1, y), arrowstyle="<|-|>", mutation_scale=9, lw=0.9, color=colour
        )
    )
    tx = (x0 + x1) / 2 if text_x is None else text_x
    ax.text(tx, y + 0.015, text, ha=ha, va="bottom", fontsize=10, color=colour)


def fig_section_geometry(path: Path) -> None:
    """Typical section: chord 2b, mid-chord, elastic axis at a b, CG at x_a b aft of the EA."""
    a_nd, x_alpha_nd = -0.2, 0.1  # shipped defaults
    ea, cg, ac = a_nd, a_nd + x_alpha_nd, -0.5
    blue, red, green, purple = "tab:blue", "tab:red", "tab:green", "tab:purple"
    x_c = np.linspace(0.0, 1.0, 400)
    x = -1.0 + 2.0 * x_c
    t = _naca0012(x_c) * 2.0

    fig, ax = plt.subplots(figsize=(11, 6.6))
    ax.fill_between(x, -t, t, color="0.9", zorder=1)
    ax.plot(x, t, "k", x, -t, "k", lw=1.2, zorder=2)
    ax.plot([-1.12, 1.12], [0, 0], color="0.35", lw=0.8, ls="-.", zorder=2)

    # Points on the chord line, labelled on one row under the section
    ax.plot(ac, 0, "^", ms=8, color=green, zorder=5)
    ax.plot(0, 0, "o", ms=5, color="k", zorder=5)
    ax.plot(cg, 0, "o", ms=10, color=red, zorder=6)
    ax.plot(ea, 0, "o", ms=10, mfc="white", mec=blue, mew=2.2, zorder=7)
    row_y = -0.19
    ax.text(ac + 0.04, row_y, "1/4 chord\n(aero. centre)", ha="right", va="top", fontsize=9,
            color=green)
    ax.text(ea - 0.02, row_y, "elastic\naxis (EA)", ha="right", va="top", fontsize=9, color=blue)
    ax.text(cg + 0.012, row_y, "CG", ha="left", va="top", fontsize=9, color=red)
    ax.text(0.02, row_y, "mid-chord", ha="left", va="top", fontsize=9)

    # Extension lines and dimensions (one row each)
    for xp, y_end in ((-1.0, -0.78), (1.0, -0.64), (0.0, -0.64), (ea, -0.52), (cg, -0.40)):
        ax.plot([xp, xp], [-0.02, y_end], color="0.45", lw=0.7, ls=":", zorder=0)
    _dim(ax, ea, cg, -0.40, r"$x_\alpha b$  (CG aft of EA)", ha="left", text_x=cg + 0.02,
         colour=red)
    _dim(ax, ea, 0.0, -0.52, r"$a\,b$  ($a<0$: EA ahead of mid-chord)", ha="left",
         text_x=0.02, colour=blue)
    _dim(ax, 0.0, 1.0, -0.64, r"$b$ (semi-chord)")
    _dim(ax, -1.0, 1.0, -0.78, r"chord $c = 2b$", text_x=-0.55)

    # Heave spring (zig-zag) from ground to the EA
    zz_y = np.linspace(0.70, 0.14, 13)
    zz_x = ea + np.r_[0.0, np.tile([0.035, -0.035], 5), 0.0, 0.0][: len(zz_y)]
    ax.plot(zz_x, zz_y, color=blue, lw=1.4, zorder=3)
    ax.plot([ea, ea], [0.14, 0.0], color=blue, lw=1.4, zorder=3)
    ax.plot([ea - 0.08, ea + 0.08], [0.70, 0.70], color="k", lw=2)
    for xh in np.linspace(ea - 0.08, ea + 0.06, 6):
        ax.plot([xh, xh + 0.03], [0.70, 0.75], color="k", lw=0.8)
    ax.text(ea + 0.055, 0.44, r"$k_h$ (heave spring)", fontsize=11, color=blue, va="center")

    # Torsion spring (spiral) around the EA
    th = np.linspace(0.0, 3.0 * np.pi, 200)
    r = 0.018 + 0.011 * th / np.pi
    ax.plot(ea + r * np.cos(th), r * np.sin(th), color=blue, lw=1.0, zorder=6)
    ax.text(ea + 0.075, 0.035, r"$k_\alpha$", fontsize=12, color=blue, zorder=8)

    # Heave DOF, beside the spring
    ax.add_patch(FancyArrowPatch((ea - 0.11, 0.52), (ea - 0.11, 0.26), arrowstyle="-|>",
                                 mutation_scale=14, lw=1.7, color=purple))
    ax.text(ea - 0.13, 0.55, r"$h$ (+ down)", ha="right", va="bottom", fontsize=11,
            color=purple)

    # Pitch DOF and moment about the EA: same sense (nose up = clockwise in this view)
    ax.add_patch(FancyArrowPatch((ea + 0.12, 0.36), (ea + 0.36, 0.14), arrowstyle="-|>",
                                 connectionstyle="arc3,rad=-0.35", mutation_scale=14, lw=1.7,
                                 color=purple))
    ax.text(ea + 0.33, 0.30, r"$\alpha$, $M_{EA}$ (+ nose up, about EA)", fontsize=11,
            color=purple)

    # Lift at the quarter chord, positive up
    ax.add_patch(FancyArrowPatch((ac, 0.11), (ac, 0.50), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.7, color=green))
    ax.text(ac - 0.03, 0.44, r"$L$ (+ up)", ha="right", fontsize=11, color=green)

    # Free stream
    ax.add_patch(FancyArrowPatch((-1.55, 0.0), (-1.18, 0.0), arrowstyle="-|>",
                                 mutation_scale=16, lw=1.8, color="k"))
    ax.text(-1.40, 0.04, r"$V$", fontsize=13, ha="center")

    ax.text(
        -1.55,
        -0.92,
        r"Shown to scale: $a=-0.2$, $x_\alpha=0.1$ (shipped defaults).  Mass $m$, pitch "
        r"inertia $I_\alpha$ about the EA,  $r_\alpha^2 = I_\alpha/(m b^2)$,  "
        r"$\mu = m/(\pi\rho b^2)$.",
        fontsize=9,
    )
    ax.set_xlim(-1.6, 1.2)
    ax.set_ylim(-0.97, 0.80)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Typical section: geometry, degrees of freedom and sign conventions")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)


def fig_theodorsen(path: Path) -> None:
    k_arr_nd = np.geomspace(0.01, 10.0, 120)
    exact = np.array([theodorsen_exact(k) for k in k_arr_nd])
    jones = np.array([theo_c(k) for k in k_arr_nd])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogx(k_arr_nd, exact.real, "k", label="F(k) exact")
    ax.semilogx(k_arr_nd, jones.real, "k--", label="F(k) Jones")
    ax.semilogx(k_arr_nd, exact.imag, "tab:red", label="G(k) exact")
    ax.semilogx(k_arr_nd, jones.imag, color="tab:red", ls="--", label="G(k) Jones")
    ax.set_xlabel("Reduced frequency k = b omega / V")
    ax.set_ylabel("C(k) = F + iG")
    ax.set_title("Theodorsen function: exact vs Jones two-pole approximation")
    ax.grid(True, which="both", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)


def fig_old_vs_new(path: Path) -> None:
    case = CASES["default_json"]
    mass, stiff, b, a, rho = _matrices(case)
    v_arr = np.linspace(0.01, 50.0, 500)
    old = {
        f"old, seeded at {w:.0f} rad/s": old_secant_solver(mass, stiff, v_arr, rho, b, a, w)
        for w in (case[5], case[6])
    }
    new_p = compute_flutter(mass, stiff, v_arr, rho, b, a)
    omega_arr = compute_structural_frequencies(mass, stiff)
    new = {f"new, MODE {i + 1} ({omega_arr[i]:.1f} rad/s)": new_p[i] for i in range(2)}

    fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for col, (title, modes) in enumerate((("Old secant solver", old), ("New solver", new))):
        for (label, p), style in zip(modes.items(), ("-", "--")):
            ax[0, col].plot(v_arr, np.imag(p) * v_arr / b, style, lw=2, label=label)
            ax[1, col].plot(v_arr, damping_g_nd(p), style, lw=2, label=label)
        ax[0, col].set_title(title)
        ax[0, col].set_ylim(0, 48)
        ax[1, col].set_ylim(-1.0, 0.5)
        ax[1, col].axhline(0.0, color="k", lw=1.5)
        ax[1, col].axhline(0.03, color="k", lw=1.5, ls=":")
        ax[1, col].set_xlabel("Velocity [m/s]")
        for row in (0, 1):
            ax[row, col].grid(True)
            ax[row, col].legend(fontsize=9)
    ax[0, 0].set_ylabel("Frequency [rad/s]")
    ax[1, 0].set_ylabel("Damping g")
    fig.suptitle("data/inputs/default.json: both old calls return the same mode")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)


def fig_vg(name: str, v_max_m_s: float, path: Path) -> NDArray[np.complex128]:
    mass, stiff, b, a, rho = _matrices(CASES[name])
    v_arr = np.linspace(v_max_m_s / 600, v_max_m_s, 600)
    p = compute_flutter(mass, stiff, v_arr, rho, b, a)
    omega_arr = compute_structural_frequencies(mass, stiff)
    modes = {f"MODE {i + 1} ({omega_arr[i]:.2f} rad/s)": p[i] for i in range(len(omega_arr))}
    plot_vg_vf(v_arr, modes, b, save_path=path)
    plt.close("all")
    return p, v_arr, b


def main() -> None:
    fig_section_geometry(_OUT_DIR / "section_geometry.png")
    fig_theodorsen(_OUT_DIR / "theodorsen_jones_vs_exact.png")
    fig_old_vs_new(_OUT_DIR / "old_vs_new_default_json.png")

    print(f"{'case':18s} {'p-k V_F':>10s} {'p-k w_F':>10s} {'k-meth V_F':>11s} "
          f"{'k-meth w_F':>11s} {'exact-C V_F':>12s}")
    for name, v_max in (("hodges_pierce", 3.0), ("shipped_defaults", 40.0), ("default_json", 50.0)):
        p, v_arr, b = fig_vg(name, v_max, _OUT_DIR / f"vg_{name}.png")
        pk = [find_flutter_speed(row, v_arr, b) for row in p]
        pk = min((x for x in pk if x is not None), default=None)
        km = k_method_flutter(CASES[name], theo_c)
        ex = k_method_flutter(CASES[name], theodorsen_exact, np.geomspace(3.0, 0.05, 400))
        fmt = lambda x, i: f"{x[i]:.5f}" if x is not None else "none"  # noqa: E731
        print(f"{name:18s} {fmt(pk, 0):>10s} {fmt(pk, 1):>10s} {fmt(km, 0):>11s} "
              f"{fmt(km, 1):>11s} {fmt(ex, 0):>12s}")


if __name__ == "__main__":
    main()
