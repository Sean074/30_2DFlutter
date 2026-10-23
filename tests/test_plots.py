"""Tests for plotting/plots.py."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from plotting import plots


@pytest.fixture(autouse=True)
def _agg_backend():
    """Run every test on the non-interactive Agg backend and restore it afterwards."""
    previous = plt.get_backend()
    plt.switch_backend("Agg")
    yield
    plt.close("all")
    plt.switch_backend(previous)


def test_plot_vg_vf_reports_not_shown_on_agg():
    velocity_range_m_s = np.array([5.0, 10.0, 20.0])
    p_arr = np.array([1.0, 0.8, 0.6]) * (-0.05 + 1j)
    shown = plots.plot_vg_vf(
        velocity_range_m_s, {"MODE 1 (10.00 rad/s)": p_arr}, b_m=0.5
    )
    assert shown is False
    assert plt.get_fignums() == []


def test_select_frozen_backend_falls_through_to_first_that_imports(monkeypatch):
    attempted = []

    def fake_switch(name: str) -> None:
        attempted.append(name)
        if name == "TkAgg":
            raise ImportError("No module named '_tkinter'")

    monkeypatch.setattr(plots.plt, "switch_backend", fake_switch)
    assert plots._select_frozen_backend(("TkAgg", "MacOSX")) == "MacOSX"
    assert attempted == ["TkAgg", "MacOSX"]


def test_select_frozen_backend_falls_back_to_agg(monkeypatch):
    def fake_switch(name: str) -> None:
        if name != "Agg":
            raise ImportError(name)

    monkeypatch.setattr(plots.plt, "switch_backend", fake_switch)
    assert plots._select_frozen_backend(("TkAgg", "MacOSX")) == "Agg"


def test_freq_axis_limits_follow_data():
    """Frequency axis runs from 0 to 1.1 x the highest frequency, not a fixed 50 rad/s."""
    freq_arr_rad_s = np.array([[0.2, 0.4], [0.9, 1.0]])
    np.testing.assert_allclose(
        plots._freq_axis_limits(freq_arr_rad_s), (0.0, 1.1), rtol=1e-4
    )


def test_plot_vg_vf_damping_axis_is_fixed(monkeypatch):
    """A heavily damped mode (g = -9.4) does not stretch the damping axis beyond [-0.4, +0.2]."""
    captured = {}
    monkeypatch.setattr(plots.plt, "close", lambda fig: captured.setdefault("fig", fig))
    velocity_range_m_s = np.array([5.0, 10.0, 20.0])
    p_arr = np.array([1.0, 0.8, 0.6]) * (np.array([-4.7, -0.1, 0.2]) + 1j)

    plots.plot_vg_vf(velocity_range_m_s, {"MODE 1": p_arr}, b_m=0.5)
    np.testing.assert_allclose(
        captured["fig"].axes[1].get_ylim(), (-0.4, 0.2), rtol=1e-4
    )
    plt.close(captured["fig"])


def test_plot_vg_vf_marks_flutter_speed(monkeypatch):
    """gamma goes -0.1 -> +0.3 between 10 and 20 m/s: flutter at 12.5 m/s gets a legend entry."""
    captured = {}
    monkeypatch.setattr(plots.plt, "close", lambda fig: captured.setdefault("fig", fig))
    velocity_range_m_s = np.array([5.0, 10.0, 20.0, 30.0])
    p_arr = np.array([2.0, 1.0, 0.6, 0.4]) * (np.array([-0.2, -0.1, 0.3, 0.5]) + 1j)

    plots.plot_vg_vf(velocity_range_m_s, {"MODE 1": p_arr}, b_m=0.5)

    ax_damping = captured["fig"].axes[1]
    labels = [text.get_text() for text in ax_damping.get_legend().get_texts()]
    assert "V_F = 12.50 m/s" in labels
    np.testing.assert_allclose(ax_damping.get_xlim(), (5.0, 30.0), rtol=1e-4)
