"""Matplotlib visualisation for flutter V-g / V-f diagrams."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends import BackendFilter, backend_registry
from numpy.typing import NDArray

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


def plot_vg_vf(
    velocity_range_m_s: NDArray[np.float64],
    flutter_p_modes: dict[str, NDArray[np.complex128]],
    b_m: float,
    speed_max_m_s: float,
    save_path: str | Path | None = None,
) -> bool:
    """
    Produce V-f (top) and V-g (bottom) flutter diagrams.

    Parameters
    ----------
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep in m/s.
    flutter_p_modes : dict[str, NDArray[np.complex128]]
        Mapping of mode label → complex p array (real=damping, imag=freq).
    b_m : float
        Semi-chord in metres (used to convert p to frequency in rad/s).
    speed_max_m_s : float
        Upper x-axis limit in m/s.
    save_path : str, Path, or None
        If provided, save the figure to this path at 150 dpi.

    Returns
    -------
    bool
        True if the figure was shown in a window; False if the active matplotlib backend is
        non-interactive (e.g. Agg when no display or GUI toolkit is available).
    """
    fig, ax = plt.subplots(2, 1, figsize=(9, 7))

    ax[0].set_title("Flutter V-g / V-f")
    ax[0].set_ylabel("Frequency [rad/s]")
    ax[1].set_ylabel("Structural damping g")
    ax[1].set_xlabel("Velocity [m/s]")

    ax[0].set_ylim(0, 50)
    ax[1].set_ylim(-0.2, 0.3)
    ax[0].set_xlim(0, speed_max_m_s)
    ax[1].set_xlim(0, speed_max_m_s)

    ax[0].grid(True)
    ax[1].grid(True)

    for label, p_arr in flutter_p_modes.items():
        freq_rad_s = np.imag(p_arr) * velocity_range_m_s / b_m
        damping_nd = np.real(p_arr) / np.imag(p_arr)
        ax[0].plot(velocity_range_m_s, freq_rad_s, label=label)
        ax[1].plot(velocity_range_m_s, damping_nd, label=label)

    zeros = np.zeros_like(velocity_range_m_s)
    ax[1].plot(velocity_range_m_s, zeros + 0.03, linestyle="dotted", color="black", linewidth=2)
    ax[1].plot(velocity_range_m_s, zeros, color="black", linewidth=2)

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
