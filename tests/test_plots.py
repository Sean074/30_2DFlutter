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
        velocity_range_m_s, {"MODE 1 (10.00 rad/s)": p_arr}, b_m=0.5, speed_max_m_s=20.0
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
