"""TUI and orchestration for the 2D aeroelastic analysis tool."""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from analysis.flutter import (
    build_mass_matrix,
    build_stiffness_matrix,
    compute_flutter,
)
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
        if warn_range is not None:
            lo, hi = warn_range
            if not lo <= val <= hi:
                print(f"{_YELLOW}  [W] Value outside typical range [{lo}, {hi}]. Confirm value is correct.{_R}")
        return val


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

def _stage1_input() -> dict | None:
    """Collect all input parameters interactively. Returns None on quit."""
    print()
    _header("INPUT", "p-k METHOD")
    print()

    print(f"{_WHITE_B}  GEOMETRY{_R}")
    _minor()
    b_m = _prompt_float("SEMI-CHORD b", "m", default=0.5, gt=0.0)
    if b_m is None:
        return None
    a_h_nd = _prompt_float("ELASTIC AXIS a", "-1 to 1", default=-0.2, warn_range=(-1.0, 1.0))
    if a_h_nd is None:
        return None
    x_alpha_nd = _prompt_float("CG OFFSET x_a", "-1 to 1", default=0.1, warn_range=(-0.5, 0.5))
    if x_alpha_nd is None:
        return None

    print()
    print(f"{_WHITE_B}  MASS PROPERTIES{_R}")
    _minor()
    mass_kg = _prompt_float("MASS", "kg", default=10.0, gt=0.0)
    if mass_kg is None:
        return None
    inertia_alpha_kg_m2 = _prompt_float("MOMENT OF INERTIA Ia", "kg.m2", default=0.5, gt=0.0)
    if inertia_alpha_kg_m2 is None:
        return None

    print()
    print(f"{_WHITE_B}  STRUCTURAL{_R}")
    _minor()
    omega_heave_rad_s = _prompt_float("HEAVE FREQ omega_h", "rad/s", default=10.0, gt=0.0)
    if omega_heave_rad_s is None:
        return None
    omega_alpha_rad_s = _prompt_float("PITCH FREQ omega_a", "rad/s", default=25.0, gt=0.0)
    if omega_alpha_rad_s is None:
        return None

    print()
    print(f"{_WHITE_B}  AERODYNAMIC{_R}")
    _minor()
    rho_kg_m3 = _prompt_float("AIR DENSITY rho", "kg/m3", default=1.225, gt=0.0)
    if rho_kg_m3 is None:
        return None
    v_min_m_s = _prompt_float("VELOCITY MIN", "m/s", default=1.0, gt=0.0)
    if v_min_m_s is None:
        return None

    while True:
        v_max_m_s = _prompt_float("VELOCITY MAX", "m/s", default=100.0, gt=0.0)
        if v_max_m_s is None:
            return None
        if v_max_m_s > v_min_m_s:
            break
        print(f"{_RED}  [E] INVALID INPUT: must be > v_min ({v_min_m_s}). Re-enter.{_R}")

    n_pts_raw = _prompt_float("N VELOCITY POINTS", "-", default=50.0, ge=2.0)
    if n_pts_raw is None:
        return None

    return {
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

    header = f"  {'V [m/s]':>{col_w}}"
    for i, lbl in enumerate(mode_labels, start=1):
        header += f"    {'g (MODE ' + str(i) + ')':>{col_w}}    {'f [rad/s]':>{col_w}}"
    print(_WHITE_B + header + _R)

    sep = "  " + "-" * col_w + ("    " + "-" * col_w + "    " + "-" * col_w) * len(mode_labels)
    print(sep)

    for i, v_m_s in enumerate(velocity_range_m_s):
        flutter_row = False
        row = f"  {v_m_s:>{col_w}.3f}"
        for lbl in mode_labels:
            p = flutter_p_modes[lbl][i]
            im = np.imag(p)
            g_nd = np.real(p) / im if im != 0.0 else 0.0
            f_rad_s = im * v_m_s / b_m
            if g_nd >= 0:
                flutter_row = True
            row += f"    {g_nd:>{col_w}.4f}    {f_rad_s:>{col_w}.3f}"
        print((_CYAN + row + _R) if flutter_row else row)


def _solve_modes(
    data: dict,
    velocity_range_m_s: np.ndarray,
) -> dict[str, np.ndarray]:
    """Build matrices, solve p-k for each mode, print progress. Returns flutter_p_modes."""
    b_m = data["b_m"]
    k_alpha_nm_rad = data["omega_alpha_rad_s"] ** 2 * data["inertia_alpha_kg_m2"]
    k_heave_n_m = data["omega_heave_rad_s"] ** 2 * data["mass_kg"]

    mass_matrix_kg = build_mass_matrix(
        data["mass_kg"], data["inertia_alpha_kg_m2"], b_m, data["x_alpha_nd"]
    )
    stiffness_matrix_n_m = build_stiffness_matrix(k_heave_n_m, k_alpha_nm_rad)

    modes = {
        f"HEAVE ({data['omega_heave_rad_s']:.2f} rad/s)": data["omega_heave_rad_s"],
        f"PITCH ({data['omega_alpha_rad_s']:.2f} rad/s)": data["omega_alpha_rad_s"],
    }

    flutter_p_modes: dict[str, np.ndarray] = {}
    n_modes = len(modes)
    for idx, (label, omega_rad_s) in enumerate(modes.items(), start=1):
        print(f"  Mode {idx} of {n_modes} : {label} ...", end="", flush=True)
        flutter_p_modes[label] = compute_flutter(
            mass_matrix_kg=mass_matrix_kg,
            stiffness_matrix_n_m=stiffness_matrix_n_m,
            velocity_range_m_s=velocity_range_m_s,
            rho_kg_m3=data["rho_kg_m3"],
            b_m=b_m,
            a_h_nd=data["a_h_nd"],
            start_omega_rad_s=omega_rad_s,
        )
        print(" done")

    return flutter_p_modes


def _stage3_output(data: dict) -> tuple[dict[str, np.ndarray], np.ndarray] | None:
    """Run analysis and print results table. Returns (flutter_p_modes, velocity_range) or None."""
    v_min_m_s: float = data["v_min_m_s"]
    v_max_m_s: float = data["v_max_m_s"]
    n_pts: int = data["n_pts"]
    velocity_range_m_s = np.linspace(v_min_m_s, v_max_m_s, n_pts)

    print()
    _header("RUNNING ANALYSIS")
    print()
    print(f"  Velocity sweep : {v_min_m_s:.1f} to {v_max_m_s:.1f} m/s  ({n_pts} steps)")
    print()

    t_start = time.perf_counter()
    try:
        flutter_p_modes = _solve_modes(data, velocity_range_m_s)
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
) -> None:
    """Print analysis summary, flutter speed estimates, and warnings."""
    print()
    _header("ANALYSIS SUMMARY")
    print()

    mode_labels = list(flutter_p_modes.keys())
    for i, lbl in enumerate(mode_labels, start=1):
        p_arr = flutter_p_modes[lbl]
        im_arr = np.imag(p_arr)
        g_arr = np.where(im_arr != 0.0, np.real(p_arr) / im_arr, 0.0)

        flutter_speed_m_s: float | None = None
        for j in range(len(g_arr) - 1):
            if g_arr[j] < 0 and g_arr[j + 1] >= 0:
                v1, v2 = velocity_range_m_s[j], velocity_range_m_s[j + 1]
                g1, g2 = g_arr[j], g_arr[j + 1]
                flutter_speed_m_s = v1 + (0.0 - g1) * (v2 - v1) / (g2 - g1)
                break

        short = lbl.split(" ")[0]
        if flutter_speed_m_s is not None:
            print(f"  {('FLUTTER DETECTED (MODE ' + str(i) + ')'):<28}: {_GREEN}YES{_R}")
            speed_str = f"{flutter_speed_m_s:.2f} m/s"
            print(f"  {('FLUTTER SPEED (MODE ' + str(i) + ')'):<28}: {_CYAN}{speed_str:<16}{_R}  ({short}, g=0 crossing)")
        else:
            print(f"  {('FLUTTER DETECTED (MODE ' + str(i) + ')'):<28}: NO")

    print()

    v_max_m_s = float(velocity_range_m_s[-1])
    for lbl, p_arr in flutter_p_modes.items():
        k_min = float(np.min(np.imag(p_arr))) * b_m / v_max_m_s
        if k_min < 0.01:
            print(f"  {_YELLOW}[W] Reduced frequency k < 0.01 at V > {0.8 * v_max_m_s:.0f} m/s "
                  f"-- quasi-steady limit approached.{_R}")
            break

    if save_path is not None:
        print(f"  {'PLOT SAVED':<28}: {_GREEN}{save_path}{_R}")

    print()
    _major()


# ---------------------------------------------------------------------------
# Main runners
# ---------------------------------------------------------------------------

def _run_pk_flutter_batch(input_path: Path, save_outputs: bool) -> None:
    """Load inputs from JSON and run p-k flutter (batch / headless mode)."""
    raw = load_json(input_path)

    v_min, v_max, n_pts = raw["velocity_range_m_s"]
    data = {
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
    }

    velocity_range_m_s = np.linspace(data["v_min_m_s"], data["v_max_m_s"], data["n_pts"])

    print()
    _header("RUNNING ANALYSIS", "BATCH MODE")
    print()
    print(f"  Input          : {input_path}")
    print(f"  Velocity sweep : {data['v_min_m_s']:.1f} to {data['v_max_m_s']:.1f} m/s"
          f"  ({data['n_pts']} steps)")
    print()

    t_start = time.perf_counter()
    flutter_p_modes = _solve_modes(data, velocity_range_m_s)
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
        )
        print(f"\n  {_GREEN}Results saved to {out_path}{_R}")
        save_path = out_path.with_suffix(".png")

    plot_vg_vf(
        velocity_range_m_s=velocity_range_m_s,
        flutter_p_modes=flutter_p_modes,
        b_m=data["b_m"],
        speed_max_m_s=data["v_max_m_s"],
        save_path=save_path,
    )

    _stage4_summary(flutter_p_modes, velocity_range_m_s, data["b_m"], save_path)


def _run_pk_flutter_interactive() -> None:
    """Run the full 4-stage interactive workflow."""
    print()
    _header("2D AEROELASTIC FLUTTER ANALYSIS", "p-k METHOD")
    print()
    print("  Type q at any prompt to quit.")
    print()

    # Stage 1: INPUT
    data = _stage1_input()
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
        data = _stage1_input()
        if data is None:
            print(f"\n  {_GREEN}Goodbye.{_R}\n")
            return

    save_outputs: bool = data.get("_save", True)

    # Stage 3: OUTPUT
    outcome = _stage3_output(data)
    if outcome is None:
        return
    flutter_p_modes, velocity_range_m_s = outcome

    save_path: Path | None = None
    if save_outputs:
        out_path = save_flutter_results(
            velocity_range_m_s=velocity_range_m_s,
            flutter_p_modes=flutter_p_modes,
            input_path="interactive",
        )
        print(f"\n  {_GREEN}Results saved to {out_path}{_R}")
        save_path = out_path.with_suffix(".png")

    plot_vg_vf(
        velocity_range_m_s=velocity_range_m_s,
        flutter_p_modes=flutter_p_modes,
        b_m=data["b_m"],
        speed_max_m_s=data["v_max_m_s"],
        save_path=save_path,
    )

    # Stage 4: CHECK OUTPUT
    _stage4_summary(flutter_p_modes, velocity_range_m_s, data["b_m"], save_path)


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

    if parsed.batch:
        input_path = Path(parsed.batch)
        if not input_path.exists():
            print(f"  {_RED}[E] File not found: {input_path}{_R}")
            return
        try:
            _run_pk_flutter_batch(input_path, save_outputs)
        except (ValueError, KeyError) as exc:
            print(f"  {_RED}[E] {exc}{_R}")
        return

    _run_pk_flutter_interactive()
