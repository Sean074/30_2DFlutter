"""TUI and orchestration for the 2D aeroelastic analysis tool."""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from analysis.flutter import (
    build_mass_matrix,
    build_stiffness_matrix,
    compute_flutter,
    compute_structural_frequencies,
    damping_g_nd,
    find_flutter_speed,
)
from analysis.k_method import compute_flutter_k
from file_io.loader import load_json
from file_io.exporter import save_flutter_results
from plotting.plots import plot_vg_vf

# ANSI colour codes — 16-colour ANSI only
_R = "\033[0m"
_CYAN_B = "\033[1;36m"
_WHITE_B = "\033[1;37m"
_RED = "\033[0;31m"
_YELLOW = "\033[0;33m"
_GREEN = "\033[0;32m"
_CYAN = "\033[0;36m"

_WIDTH = 72

# defaults/ sits next to the package in a source checkout and in the unpacked PyInstaller bundle.
_DEFAULTS_DIR = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent)) / "defaults"
_DEFAULT_KEYS = (
    "b_m",
    "a_h_nd",
    "x_alpha_nd",
    "mass_kg",
    "inertia_alpha_kg_m2",
    "omega_alpha_rad_s",
    "omega_heave_rad_s",
    "rho_kg_m3",
    "velocity_range_m_s",
    "k_tol_nd",
    "max_iter_nd",
    "method",
    "n_k_nd",
    "flap",
    "c_hinge_nd",
    "x_beta_nd",
    "inertia_beta_kg_m2",
    "omega_beta_rad_s",
)
# Flap (third DOF) inputs: an input file gives all of them or none.
_FLAP_KEYS = ("c_hinge_nd", "x_beta_nd", "inertia_beta_kg_m2", "omega_beta_rad_s")

# Selectable solvers: input-file / prompt code -> display name.
_METHODS = {"pk": "p-k METHOD", "k": "k METHOD"}


def _load_defaults(defaults_dir: Path = _DEFAULTS_DIR) -> tuple[dict, dict]:
    """Load default_flutter.json and constants.json. Raises ValueError if a key is missing."""
    defaults = load_json(defaults_dir / "default_flutter.json")
    constants = load_json(defaults_dir / "constants.json")
    missing = [key for key in _DEFAULT_KEYS if key not in defaults]
    if "pi_nd" not in constants:
        missing.append("pi_nd")
    if missing:
        raise ValueError(f"missing keys in {defaults_dir}: {', '.join(missing)}")
    _check_method(defaults["method"])
    return defaults, constants


def _check_method(method: str) -> str:
    """Return method if it is a known solver code; ValueError otherwise."""
    if method not in _METHODS:
        raise ValueError(f"unknown method '{method}'; use one of {', '.join(_METHODS)}")
    return method


# ---------------------------------------------------------------------------
# Low-level terminal helpers
# ---------------------------------------------------------------------------

def _major() -> None:
    print(_CYAN_B + "=" * _WIDTH + _R)


def _minor() -> None:
    print("  " + "-" * 20)


def _header(left: str, right: str = "") -> None:
    _major()
    if right:
        content_w = _WIDTH - 2
        line = f"  {left:<{content_w - len(right)}}{right}"
    else:
        line = f"  {left}"
    print(_CYAN_B + line + _R)
    _major()


def _prompt(label_padded: str, default: str | None = None) -> str | None:
    """Read one line from stdin. Returns None on q/Q or EOF; '' on empty with no default."""
    hint = f"[{default}] " if default is not None else ""
    try:
        raw = input(f"{_WHITE_B}{label_padded}: {hint}{_R}").strip()
    except (EOFError, KeyboardInterrupt):
        print()
        return None
    if raw.lower() == "q":
        return None
    if raw == "" and default is not None:
        return str(default)
    return raw


def _prompt_float(
    label: str,
    unit: str,
    *,
    default: float | None = None,
    gt: float | None = None,
    ge: float | None = None,
    lt: float | None = None,
    warn_range: tuple[float, float] | None = None,
) -> float | None:
    """Prompt for a float with validation. Returns None on quit."""
    label_str = f"{label} [{unit}]" if unit else label
    padded = f"  {label_str:<28}"
    default_str = str(default) if default is not None else None

    while True:
        raw = _prompt(padded, default_str)
        if raw is None:
            return None
        if raw == "":
            print(f"{_RED}  [E] Value required. Re-enter.{_R}")
            continue
        try:
            val = float(raw)
        except ValueError:
            print(f"{_RED}  [E] INVALID INPUT: must be a number. Re-enter.{_R}")
            continue
        if gt is not None and val <= gt:
            print(f"{_RED}  [E] INVALID INPUT: must be > {gt}. Re-enter.{_R}")
            continue
        if ge is not None and val < ge:
            print(f"{_RED}  [E] INVALID INPUT: must be >= {ge}. Re-enter.{_R}")
            continue
        if lt is not None and val >= lt:
            print(f"{_RED}  [E] INVALID INPUT: must be < {lt}. Re-enter.{_R}")
            continue
        if warn_range is not None:
            lo, hi = warn_range
            if not lo <= val <= hi:
                print(f"{_YELLOW}  [W] Value outside typical range [{lo}, {hi}]. Confirm value is correct.{_R}")
        return val


def _prompt_choice(label: str, choices: tuple[str, ...], default: str) -> str | None:
    """Prompt for one of choices (case-insensitive). Returns None on quit."""
    padded = f"  {label + ' [' + '/'.join(choices) + ']':<28}"
    while True:
        raw = _prompt(padded, default)
        if raw is None:
            return None
        if raw.lower() in choices:
            return raw.lower()
        print(f"{_RED}  [E] INVALID INPUT: enter one of {', '.join(choices)}. Re-enter.{_R}")


def _confirm(prompt: str, default: bool = True) -> bool | None:
    """Y/n prompt. Returns None on quit."""
    hint = "[Y/n]" if default else "[y/N]"
    padded = f"  {prompt:<28}"
    while True:
        try:
            raw = input(f"{_WHITE_B}{padded}: {hint} {_R}").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print()
            return None
        if raw == "q":
            return None
        if raw in ("y", "yes"):
            return True
        if raw in ("n", "no"):
            return False
        if raw == "":
            return default
        print(f"{_RED}  [E] Enter Y or N.{_R}")


# ---------------------------------------------------------------------------
# Stage 1 — INPUT
# ---------------------------------------------------------------------------

def _stage1_input(defaults: dict) -> dict | None:
    """Collect all input parameters interactively, pre-filled from defaults. None on quit."""
    print()
    _header("INPUT")
    print()

    v_min_default_m_s, v_max_default_m_s, n_pts_default = defaults["velocity_range_m_s"]

    print(f"{_WHITE_B}  SOLVER{_R}")
    _minor()
    method = _prompt_choice("METHOD", tuple(_METHODS), defaults["method"])
    if method is None:
        return None

    print()
    print(f"{_WHITE_B}  GEOMETRY{_R}")
    _minor()
    b_m = _prompt_float("SEMI-CHORD b", "m", default=defaults["b_m"], gt=0.0)
    if b_m is None:
        return None
    a_h_nd = _prompt_float(
        "ELASTIC AXIS a", "-1 to 1", default=defaults["a_h_nd"], warn_range=(-1.0, 1.0)
    )
    if a_h_nd is None:
        return None
    x_alpha_nd = _prompt_float(
        "CG OFFSET x_a", "-1 to 1", default=defaults["x_alpha_nd"], warn_range=(-0.5, 0.5)
    )
    if x_alpha_nd is None:
        return None

    print()
    print(f"{_WHITE_B}  MASS PROPERTIES{_R}")
    _minor()
    mass_kg = _prompt_float("MASS", "kg", default=defaults["mass_kg"], gt=0.0)
    if mass_kg is None:
        return None
    inertia_alpha_kg_m2 = _prompt_float(
        "MOMENT OF INERTIA Ia", "kg.m2", default=defaults["inertia_alpha_kg_m2"], gt=0.0
    )
    if inertia_alpha_kg_m2 is None:
        return None

    print()
    print(f"{_WHITE_B}  STRUCTURAL{_R}")
    _minor()
    omega_heave_rad_s = _prompt_float(
        "HEAVE FREQ omega_h", "rad/s", default=defaults["omega_heave_rad_s"], gt=0.0
    )
    if omega_heave_rad_s is None:
        return None
    omega_alpha_rad_s = _prompt_float(
        "PITCH FREQ omega_a", "rad/s", default=defaults["omega_alpha_rad_s"], gt=0.0
    )
    if omega_alpha_rad_s is None:
        return None

    print()
    print(f"{_WHITE_B}  FLAP{_R}")
    _minor()
    flap = _stage1_flap(defaults)
    if flap is None:
        return None

    print()
    print(f"{_WHITE_B}  AERODYNAMIC{_R}")
    _minor()
    rho_kg_m3 = _prompt_float("AIR DENSITY rho", "kg/m3", default=defaults["rho_kg_m3"], gt=0.0)
    if rho_kg_m3 is None:
        return None
    v_min_m_s = _prompt_float("VELOCITY MIN", "m/s", default=v_min_default_m_s, gt=0.0)
    if v_min_m_s is None:
        return None

    while True:
        v_max_m_s = _prompt_float("VELOCITY MAX", "m/s", default=v_max_default_m_s, gt=0.0)
        if v_max_m_s is None:
            return None
        if v_max_m_s > v_min_m_s:
            break
        print(f"{_RED}  [E] INVALID INPUT: must be > v_min ({v_min_m_s}). Re-enter.{_R}")

    n_pts_raw = _prompt_float("N VELOCITY POINTS", "-", default=float(n_pts_default), ge=2.0)
    if n_pts_raw is None:
        return None

    return {
        "method": method,
        "b_m": b_m,
        "a_h_nd": a_h_nd,
        "x_alpha_nd": x_alpha_nd,
        "mass_kg": mass_kg,
        "inertia_alpha_kg_m2": inertia_alpha_kg_m2,
        "omega_heave_rad_s": omega_heave_rad_s,
        "omega_alpha_rad_s": omega_alpha_rad_s,
        "rho_kg_m3": rho_kg_m3,
        "v_min_m_s": v_min_m_s,
        "v_max_m_s": v_max_m_s,
        "n_pts": int(n_pts_raw),
        **flap,
    }


def _stage1_flap(defaults: dict) -> dict | None:
    """Ask whether to model the flap and, if so, its inputs. Returns None on quit."""
    include = _confirm("INCLUDE FLAP (3-DOF)?", default=bool(defaults["flap"]))
    if include is None:
        return None
    if not include:
        return {"flap": False}
    c_hinge_nd = _prompt_float(
        "HINGE POSITION c", "-1 to 1", default=defaults["c_hinge_nd"], gt=-1.0, lt=1.0
    )
    if c_hinge_nd is None:
        return None
    x_beta_nd = _prompt_float(
        "FLAP CG OFFSET x_b", "-", default=defaults["x_beta_nd"], warn_range=(-0.1, 0.1)
    )
    if x_beta_nd is None:
        return None
    inertia_beta_kg_m2 = _prompt_float(
        "FLAP INERTIA Ib", "kg.m2", default=defaults["inertia_beta_kg_m2"], gt=0.0
    )
    if inertia_beta_kg_m2 is None:
        return None
    omega_beta_rad_s = _prompt_float(
        "FLAP FREQ omega_b", "rad/s", default=defaults["omega_beta_rad_s"], gt=0.0
    )
    if omega_beta_rad_s is None:
        return None
    return {
        "flap": True,
        "c_hinge_nd": c_hinge_nd,
        "x_beta_nd": x_beta_nd,
        "inertia_beta_kg_m2": inertia_beta_kg_m2,
        "omega_beta_rad_s": omega_beta_rad_s,
    }


# ---------------------------------------------------------------------------
# Stage 2 — CHECK INPUT
# ---------------------------------------------------------------------------

def _stage2_check(data: dict) -> bool | None:
    """Echo all inputs and ask user to confirm. Returns None on quit, False to re-enter."""
    print()
    _header("INPUT SUMMARY", "CONFIRM BEFORE RUNNING")
    print()

    def _row(label: str, unit: str, val: float) -> None:
        label_str = f"{label} [{unit}]" if unit else label
        print(f"  {label_str:<28}: {val:.4f}")

    print(f"{_WHITE_B}  SOLVER{_R}")
    _minor()
    print(f"  {'METHOD':<28}: {_METHODS[data['method']]}")

    print()
    print(f"{_WHITE_B}  GEOMETRY{_R}")
    _minor()
    _row("SEMI-CHORD b", "m", data["b_m"])
    _row("ELASTIC AXIS a", "-1 to 1", data["a_h_nd"])
    _row("CG OFFSET x_a", "-1 to 1", data["x_alpha_nd"])

    print()
    print(f"{_WHITE_B}  MASS PROPERTIES{_R}")
    _minor()
    _row("MASS", "kg", data["mass_kg"])
    _row("MOMENT OF INERTIA Ia", "kg.m2", data["inertia_alpha_kg_m2"])

    print()
    print(f"{_WHITE_B}  STRUCTURAL{_R}")
    _minor()
    _row("HEAVE FREQ omega_h", "rad/s", data["omega_heave_rad_s"])
    _row("PITCH FREQ omega_a", "rad/s", data["omega_alpha_rad_s"])

    print()
    print(f"{_WHITE_B}  FLAP{_R}")
    _minor()
    print(f"  {'INCLUDE FLAP (3-DOF)':<28}: {'YES' if data['flap'] else 'NO'}")
    if data["flap"]:
        _row("HINGE POSITION c", "-1 to 1", data["c_hinge_nd"])
        _row("FLAP CG OFFSET x_b", "-", data["x_beta_nd"])
        _row("FLAP INERTIA Ib", "kg.m2", data["inertia_beta_kg_m2"])
        _row("FLAP FREQ omega_b", "rad/s", data["omega_beta_rad_s"])

    print()
    print(f"{_WHITE_B}  AERODYNAMIC{_R}")
    _minor()
    _row("AIR DENSITY rho", "kg/m3", data["rho_kg_m3"])
    _row("VELOCITY MIN", "m/s", data["v_min_m_s"])
    _row("VELOCITY MAX", "m/s", data["v_max_m_s"])
    print(f"  {'N VELOCITY POINTS':<28}: {data['n_pts']}")

    print()
    proceed = _confirm("Proceed?", default=True)
    if proceed is None:
        return None
    if not proceed:
        return False

    save = _confirm("Save results to data/outputs?", default=True)
    if save is None:
        return None

    data["_save"] = save
    return True


# ---------------------------------------------------------------------------
# Stage 3 — OUTPUT
# ---------------------------------------------------------------------------

def _print_results_table(
    flutter_p_modes: dict[str, np.ndarray],
    velocity_range_m_s: np.ndarray,
    b_m: float,
) -> None:
    """Print velocity-damping / velocity-frequency table. Highlights flutter rows in cyan."""
    col_w = 10
    mode_labels = list(flutter_p_modes.keys())
    # Two modes fit with 4-space gaps; three need 1-space gaps to stay within 80 columns.
    gap = " " * (4 if len(mode_labels) <= 2 else 1)

    header = f"  {'V [m/s]':>{col_w}}"
    for i, lbl in enumerate(mode_labels, start=1):
        header += f"{gap}{'g (MODE ' + str(i) + ')':>{col_w}}{gap}{'f [rad/s]':>{col_w}}"
    print(_WHITE_B + header + _R)

    sep = "  " + "-" * col_w + (gap + "-" * col_w + gap + "-" * col_w) * len(mode_labels)
    print(sep)

    for i, v_m_s in enumerate(velocity_range_m_s):
        flutter_row = False
        row = f"  {v_m_s:>{col_w}.3f}"
        for lbl in mode_labels:
            p = flutter_p_modes[lbl][i]
            if np.isnan(p):
                # The mode's branch does not reach this speed.
                row += f"{gap}{'-':>{col_w}}{gap}{'-':>{col_w}}"
                continue
            im = np.imag(p)
            g_nd = float(damping_g_nd(p)) if im != 0.0 else 0.0
            f_rad_s = im * v_m_s / b_m
            if g_nd >= 0:
                flutter_row = True
            row += f"{gap}{g_nd:>{col_w}.4f}{gap}{f_rad_s:>{col_w}.3f}"
        print((_CYAN + row + _R) if flutter_row else row)


def _solve_modes(
    data: dict,
    velocity_range_m_s: np.ndarray,
    solver: dict,
) -> dict[str, np.ndarray]:
    """Build matrices, run the chosen solver, print progress. Returns flutter_p_modes."""
    b_m = data["b_m"]
    k_alpha_nm_rad = data["omega_alpha_rad_s"] ** 2 * data["inertia_alpha_kg_m2"]
    k_heave_n_m = data["omega_heave_rad_s"] ** 2 * data["mass_kg"]

    if data["flap"]:
        c_hinge_nd = data["c_hinge_nd"]
        mass_matrix_kg = build_mass_matrix(
            data["mass_kg"],
            data["inertia_alpha_kg_m2"],
            b_m,
            data["x_alpha_nd"],
            data["x_beta_nd"],
            data["inertia_beta_kg_m2"],
            c_hinge_nd,
            data["a_h_nd"],
        )
        k_beta_nm_rad = data["omega_beta_rad_s"] ** 2 * data["inertia_beta_kg_m2"]
    else:
        c_hinge_nd = None
        mass_matrix_kg = build_mass_matrix(
            data["mass_kg"], data["inertia_alpha_kg_m2"], b_m, data["x_alpha_nd"]
        )
        k_beta_nm_rad = None
    stiffness_matrix_n_m = build_stiffness_matrix(k_heave_n_m, k_alpha_nm_rad, k_beta_nm_rad)

    omega_arr_rad_s = compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)
    print(f"  Solving {len(omega_arr_rad_s)} modes ...", end="", flush=True)
    common = {
        "mass_matrix_kg": mass_matrix_kg,
        "stiffness_matrix_n_m": stiffness_matrix_n_m,
        "velocity_range_m_s": velocity_range_m_s,
        "rho_kg_m3": data["rho_kg_m3"],
        "b_m": b_m,
        "a_h_nd": data["a_h_nd"],
        "pi_nd": solver["pi_nd"],
        "c_hinge_nd": c_hinge_nd,
    }
    if data["method"] == "k":
        p_arr_modes = compute_flutter_k(**common, n_k_nd=solver["n_k_nd"])
    else:
        p_arr_modes = compute_flutter(
            **common, k_tol_nd=solver["k_tol_nd"], max_iter_nd=solver["max_iter_nd"]
        )
    print(" done")

    return {
        f"MODE {idx} ({omega_rad_s:.2f} rad/s)": p_arr_modes[idx - 1]
        for idx, omega_rad_s in enumerate(omega_arr_rad_s, start=1)
    }


def _stage3_output(data: dict, solver: dict) -> tuple[dict[str, np.ndarray], np.ndarray] | None:
    """Run analysis and print results table. Returns (flutter_p_modes, velocity_range) or None."""
    v_min_m_s: float = data["v_min_m_s"]
    v_max_m_s: float = data["v_max_m_s"]
    n_pts: int = data["n_pts"]
    velocity_range_m_s = np.linspace(v_min_m_s, v_max_m_s, n_pts)

    print()
    _header("RUNNING ANALYSIS", _METHODS[data["method"]])
    print()
    print(f"  Velocity sweep : {v_min_m_s:.1f} to {v_max_m_s:.1f} m/s  ({n_pts} steps)")
    print()

    t_start = time.perf_counter()
    try:
        flutter_p_modes = _solve_modes(data, velocity_range_m_s, solver)
    except (ValueError, KeyError) as exc:
        print(f"\n  {_RED}[E] {exc}{_R}")
        return None
    elapsed_s = time.perf_counter() - t_start

    print()
    print(f"  Elapsed : {elapsed_s:.2f} s")

    print()
    _header("RESULTS")
    print()
    _print_results_table(flutter_p_modes, velocity_range_m_s, data["b_m"])

    return flutter_p_modes, velocity_range_m_s


# ---------------------------------------------------------------------------
# Stage 4 — CHECK OUTPUT
# ---------------------------------------------------------------------------

def _stage4_summary(
    flutter_p_modes: dict[str, np.ndarray],
    velocity_range_m_s: np.ndarray,
    b_m: float,
    save_path: Path | None,
    method: str = "pk",
) -> None:
    """Print analysis summary, flutter speed estimates, and warnings."""
    print()
    _header("ANALYSIS SUMMARY")
    print()

    mode_labels = list(flutter_p_modes.keys())
    for i, lbl in enumerate(mode_labels, start=1):
        flutter_point = find_flutter_speed(flutter_p_modes[lbl], velocity_range_m_s, b_m)

        if flutter_point is not None:
            v_flutter_m_s, omega_flutter_rad_s = flutter_point
            print(f"  {('FLUTTER DETECTED (MODE ' + str(i) + ')'):<28}: {_GREEN}YES{_R}")
            print(f"  {('FLUTTER SPEED (MODE ' + str(i) + ')'):<28}: "
                  f"{_CYAN}{v_flutter_m_s:.2f} m/s{_R}  (g=0 crossing)")
            print(f"  {('FLUTTER FREQ (MODE ' + str(i) + ')'):<28}: "
                  f"{_CYAN}{omega_flutter_rad_s:.2f} rad/s{_R}")
        else:
            print(f"  {('FLUTTER DETECTED (MODE ' + str(i) + ')'):<28}: NO")

    print()

    # imag(p) is the reduced frequency k; report the lowest speed where any mode drops below 0.01.
    k_arr_nd = np.imag(np.array(list(flutter_p_modes.values())))
    low_k_arr = np.flatnonzero(np.any(k_arr_nd < 0.01, axis=0))
    if len(low_k_arr) > 0:
        print(f"  {_YELLOW}[W] Reduced frequency k < 0.01 at V >= "
              f"{velocity_range_m_s[low_k_arr[0]]:.1f} m/s -- quasi-steady limit approached.{_R}")

    if method == "k":
        print(f"  {_YELLOW}[W] k-method g is the damping required for neutral motion; only "
              f"g = 0 is physical.{_R}")
    # k-method: V(k) folds back. p-k: frequency -> 0 (overdamped or divergence) or no convergence.
    end_reason = "V(k) folds back" if method == "k" else "p-k cannot follow it"
    for i, lbl in enumerate(mode_labels, start=1):
        known_arr = np.flatnonzero(~np.isnan(flutter_p_modes[lbl]))
        if len(known_arr) == 0:
            print(f"  {_YELLOW}[W] MODE {i} not found in the velocity range.{_R}")
        elif known_arr[-1] < len(velocity_range_m_s) - 1:
            print(f"  {_YELLOW}[W] MODE {i} ends at V = "
                  f"{velocity_range_m_s[known_arr[-1]]:.1f} m/s ({end_reason}).{_R}")

    if save_path is not None:
        print(f"  {'PLOT SAVED':<28}: {_GREEN}{save_path}{_R}")

    print()
    _major()


# ---------------------------------------------------------------------------
# Main runners
# ---------------------------------------------------------------------------

def _run_pk_flutter_batch(input_path: Path, save_outputs: bool, solver: dict) -> None:
    """Load inputs from JSON and run p-k flutter (batch / headless mode)."""
    raw = load_json(input_path)

    v_min, v_max, n_pts = raw["velocity_range_m_s"]
    data = {
        "method": _check_method(raw.get("method", "pk")),
        "b_m": float(raw["b_m"]),
        "a_h_nd": float(raw["a_h_nd"]),
        "x_alpha_nd": float(raw["x_alpha_nd"]),
        "mass_kg": float(raw["mass_kg"]),
        "inertia_alpha_kg_m2": float(raw["inertia_alpha_kg_m2"]),
        "omega_alpha_rad_s": float(raw["omega_alpha_rad_s"]),
        "omega_heave_rad_s": float(raw["omega_heave_rad_s"]),
        "rho_kg_m3": float(raw["rho_kg_m3"]),
        "v_min_m_s": float(v_min),
        "v_max_m_s": float(v_max),
        "n_pts": int(n_pts),
        **_flap_from_input(raw),
    }

    velocity_range_m_s = np.linspace(data["v_min_m_s"], data["v_max_m_s"], data["n_pts"])

    print()
    _header("RUNNING ANALYSIS", f"BATCH MODE  {_METHODS[data['method']]}")
    print()
    print(f"  Input          : {input_path}")
    print(f"  Velocity sweep : {data['v_min_m_s']:.1f} to {data['v_max_m_s']:.1f} m/s"
          f"  ({data['n_pts']} steps)")
    print()

    t_start = time.perf_counter()
    flutter_p_modes = _solve_modes(data, velocity_range_m_s, solver)
    elapsed_s = time.perf_counter() - t_start

    print()
    print(f"  Elapsed : {elapsed_s:.2f} s")

    print()
    _header("RESULTS")
    print()
    _print_results_table(flutter_p_modes, velocity_range_m_s, data["b_m"])

    save_path: Path | None = None
    if save_outputs:
        out_path = save_flutter_results(
            velocity_range_m_s=velocity_range_m_s,
            flutter_p_modes=flutter_p_modes,
            input_path=str(input_path),
            method=data["method"],
        )
        print(f"\n  {_GREEN}Results saved to {out_path}{_R}")
        save_path = out_path.with_suffix(".png")

    plot_shown = plot_vg_vf(
        velocity_range_m_s=velocity_range_m_s,
        flutter_p_modes=flutter_p_modes,
        b_m=data["b_m"],
        save_path=save_path,
        title=f"Flutter V-g / V-f ({_METHODS[data['method']]})",
    )
    if not plot_shown:
        hint = "see the saved PNG" if save_path is not None else "the plot was not saved"
        print(f"\n  {_YELLOW}[W] No plot window available (no GUI backend); {hint}.{_R}")

    _stage4_summary(
        flutter_p_modes, velocity_range_m_s, data["b_m"], save_path, data["method"]
    )


def _flap_from_input(raw: dict) -> dict:
    """Flap entries of an input file: all of _FLAP_KEYS (3-DOF) or none (2-DOF)."""
    present = [key for key in _FLAP_KEYS if key in raw]
    if not present:
        return {"flap": False}
    missing = [key for key in _FLAP_KEYS if key not in raw]
    if missing:
        raise ValueError(f"flap input incomplete; missing {', '.join(missing)}")
    return {"flap": True} | {key: float(raw[key]) for key in _FLAP_KEYS}


def _run_pk_flutter_interactive(defaults: dict, solver: dict) -> None:
    """Run the full 4-stage interactive workflow."""
    print()
    _header("2D AEROELASTIC FLUTTER ANALYSIS", "p-k / k METHOD")
    print()
    print("  Type q at any prompt to quit.")
    print()

    # Stage 1: INPUT
    data = _stage1_input(defaults)
    if data is None:
        print(f"\n  {_GREEN}Goodbye.{_R}\n")
        return

    # Stage 2: CHECK INPUT  (loop until confirmed or quit)
    while True:
        result = _stage2_check(data)
        if result is None:
            print(f"\n  {_GREEN}Goodbye.{_R}\n")
            return
        if result is True:
            break
        data = _stage1_input(defaults)
        if data is None:
            print(f"\n  {_GREEN}Goodbye.{_R}\n")
            return

    save_outputs: bool = data.get("_save", True)

    # Stage 3: OUTPUT
    outcome = _stage3_output(data, solver)
    if outcome is None:
        return
    flutter_p_modes, velocity_range_m_s = outcome

    save_path: Path | None = None
    if save_outputs:
        out_path = save_flutter_results(
            velocity_range_m_s=velocity_range_m_s,
            flutter_p_modes=flutter_p_modes,
            input_path="interactive",
            method=data["method"],
        )
        print(f"\n  {_GREEN}Results saved to {out_path}{_R}")
        save_path = out_path.with_suffix(".png")

    plot_shown = plot_vg_vf(
        velocity_range_m_s=velocity_range_m_s,
        flutter_p_modes=flutter_p_modes,
        b_m=data["b_m"],
        save_path=save_path,
        title=f"Flutter V-g / V-f ({_METHODS[data['method']]})",
    )
    if not plot_shown:
        hint = "see the saved PNG" if save_path is not None else "the plot was not saved"
        print(f"\n  {_YELLOW}[W] No plot window available (no GUI backend); {hint}.{_R}")

    # Stage 4: CHECK OUTPUT
    _stage4_summary(
        flutter_p_modes, velocity_range_m_s, data["b_m"], save_path, data["method"]
    )


def run_tui(args: list[str] | None = None) -> None:
    """
    Entry point for the interactive TUI and batch mode.

    Parameters
    ----------
    args : list[str] or None
        CLI arguments (defaults to sys.argv if None).
    """
    parser = argparse.ArgumentParser(description="2D Aeroelastic Flutter Analysis")
    parser.add_argument("--batch", metavar="INPUT_FILE", help="Run headless with this input file")
    parser.add_argument("--no-save", action="store_true", help="Skip saving outputs to disk")
    parsed = parser.parse_args(args)

    save_outputs = not parsed.no_save

    try:
        defaults, constants = _load_defaults()
    except (FileNotFoundError, ValueError) as exc:
        print(f"  {_RED}[E] Cannot load defaults: {exc}{_R}")
        return
    solver = {
        "k_tol_nd": float(defaults["k_tol_nd"]),
        "max_iter_nd": int(defaults["max_iter_nd"]),
        "pi_nd": float(constants["pi_nd"]),
        "n_k_nd": int(defaults["n_k_nd"]),
    }

    if parsed.batch:
        input_path = Path(parsed.batch)
        if not input_path.exists():
            print(f"  {_RED}[E] File not found: {input_path}{_R}")
            return
        try:
            _run_pk_flutter_batch(input_path, save_outputs, solver)
        except (ValueError, KeyError) as exc:
            print(f"  {_RED}[E] {exc}{_R}")
        return

    _run_pk_flutter_interactive(defaults, solver)
