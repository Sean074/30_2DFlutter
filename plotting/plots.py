"""Matplotlib visualisation for flutter V-g / V-f diagrams."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends import BackendFilter, backend_registry
from numpy.typing import NDArray

from analysis.flutter import damping_g_nd, find_flutter_speed

# A frozen (PyInstaller) build cannot rely on matplotlib's backend auto-detection, so pick the
# first window backend that imports, falling back to Agg (file output only) if none does.
_FROZEN_BACKENDS = ("TkAgg", "MacOSX") if sys.platform == "darwin" else ("TkAgg",)


def _select_frozen_backend(candidates: tuple[str, ...]) -> str:
    """Switch to the first backend in candidates that imports; Agg if none. Returns its name."""
    for name in candidates:
        try:
            plt.switch_backend(name)
            return name
        except ImportError:
            continue
    plt.switch_backend("Agg")
    return "Agg"


if getattr(sys, "frozen", False):
    _select_frozen_backend(_FROZEN_BACKENDS)


# Damping window [-L, +L * _G_POS_TO_NEG_RATIO_ND], with L clamped to this range so the g = 0
# crossing stays readable even when one mode is very heavily damped.
_G_POS_TO_NEG_RATIO_ND = 0.5
_G_NEG_LIM_RANGE_ND = (0.2, 1.0)
_FREQ_MARGIN_ND = 1.1


def _axis_limits(
    freq_arr_rad_s: NDArray[np.float64],
    g_arr_nd: NDArray[np.float64],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return ((freq_min, freq_max) in rad/s, (g_min, g_max)) for the V-f and V-g axes."""
    freq_max_rad_s = float(np.nanmax(freq_arr_rad_s)) * _FREQ_MARGIN_ND
    if not freq_max_rad_s > 0.0:
        freq_max_rad_s = 1.0

    g_neg_needed_nd = max(
        -float(np.nanmin(g_arr_nd)),
        float(np.nanmax(g_arr_nd)) / _G_POS_TO_NEG_RATIO_ND,
    )
    g_neg_lim_nd = float(np.clip(g_neg_needed_nd, *_G_NEG_LIM_RANGE_ND))
    return (0.0, freq_max_rad_s), (-g_neg_lim_nd, g_neg_lim_nd * _G_POS_TO_NEG_RATIO_ND)


def plot_vg_vf(
    velocity_range_m_s: NDArray[np.float64],
    flutter_p_modes: dict[str, NDArray[np.complex128]],
    b_m: float,
    save_path: str | Path | None = None,
    freq_lim_rad_s: tuple[float, float] | None = None,
    g_lim_nd: tuple[float, float] | None = None,
    title: str = "Flutter V-g / V-f",
) -> bool:
    """
    Produce V-f (top) and V-g (bottom) flutter diagrams.

    Axis limits follow the data: frequency from 0 to 1.1 x the highest modal frequency; damping
    g = 2 real(p) / imag(p) in a window [-L, +L/2] around g = 0, with L between 0.2 and 1.0.
    The flutter point of each mode that flutters is marked on both plots.

    Parameters
    ----------
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep in m/s.
    flutter_p_modes : dict[str, NDArray[np.complex128]]
        Mapping of mode label → complex p array (real=damping, imag=freq).
    b_m : float
        Semi-chord in metres (used to convert p to frequency in rad/s).
    save_path : str, Path, or None
        If provided, save the figure to this path at 150 dpi.
    freq_lim_rad_s : tuple[float, float] or None
        Frequency axis limits in rad/s, overriding the automatic ones.
    g_lim_nd : tuple[float, float] or None
        Damping axis limits, overriding the automatic ones.
    title : str
        Figure title.

    Returns
    -------
    bool
        True if the figure was shown in a window; False if the active matplotlib backend is
        non-interactive (e.g. Agg when no display or GUI toolkit is available).
    """
    fig, ax = plt.subplots(2, 1, figsize=(9, 7))

    ax[0].set_title(title)
    ax[0].set_ylabel("Frequency [rad/s]")
    ax[1].set_ylabel("Structural damping g")
    ax[1].set_xlabel("Velocity [m/s]")

    p_arr_modes = np.array(list(flutter_p_modes.values()))
    freq_arr_rad_s = np.imag(p_arr_modes) * velocity_range_m_s / b_m
    g_arr_nd = damping_g_nd(p_arr_modes)
    auto_freq_lim_rad_s, auto_g_lim_nd = _axis_limits(freq_arr_rad_s, g_arr_nd)

    ax[0].set_ylim(freq_lim_rad_s if freq_lim_rad_s is not None else auto_freq_lim_rad_s)
    ax[1].set_ylim(g_lim_nd if g_lim_nd is not None else auto_g_lim_nd)
    for axis in ax:
        axis.set_xlim(velocity_range_m_s[0], velocity_range_m_s[-1])
        axis.grid(True)

    for mode, label in enumerate(flutter_p_modes):
        lines = ax[0].plot(velocity_range_m_s, freq_arr_rad_s[mode], label=label)
        ax[1].plot(velocity_range_m_s, g_arr_nd[mode], label=label)

        flutter_point = find_flutter_speed(p_arr_modes[mode], velocity_range_m_s, b_m)
        if flutter_point is not None:
            v_flutter_m_s, omega_flutter_rad_s = flutter_point
            colour = lines[0].get_color()
            ax[0].plot(v_flutter_m_s, omega_flutter_rad_s, "o", color=colour)
            ax[1].plot(v_flutter_m_s, 0.0, "o", color=colour)
            for axis in ax:
                axis.axvline(v_flutter_m_s, color=colour, linestyle="--", linewidth=1)
            ax[1].plot([], [], "--", color=colour, label=f"V_F = {v_flutter_m_s:.2f} m/s")

    ax[1].axhline(0.03, linestyle="dotted", color="black", linewidth=2)
    ax[1].axhline(0.0, color="black", linewidth=2)

    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150)

    non_interactive = backend_registry.list_builtin(BackendFilter.NON_INTERACTIVE)
    if plt.get_backend().lower() in (name.lower() for name in non_interactive):
        plt.close(fig)
        return False
    plt.show()
    return True
