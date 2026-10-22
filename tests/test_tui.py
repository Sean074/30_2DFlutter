"""Tests for ui/tui.py."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ui import tui

# Deliberately unlike the shipped defaults, to prove the TUI does not use hardcoded values.
_DEFAULTS = {
    "b_m": 0.7,
    "a_h_nd": -0.3,
    "x_alpha_nd": 0.2,
    "mass_kg": 12.0,
    "inertia_alpha_kg_m2": 0.8,
    "omega_alpha_rad_s": 30.0,
    "omega_heave_rad_s": 12.0,
    "rho_kg_m3": 1.1,
    "velocity_range_m_s": [2.0, 80.0, 40],
    "k_tol_nd": 1e-9,
    "max_iter_nd": 150,
}
_CONSTANTS = {"pi_nd": 3.141592653589793}


def _fake_load_json(defaults: dict, constants: dict):
    def load(path: Path) -> dict:
        return dict(defaults if path.name == "default_flutter.json" else constants)

    return load


def test_load_defaults_returns_both_files(monkeypatch):
    monkeypatch.setattr(tui, "load_json", _fake_load_json(_DEFAULTS, _CONSTANTS))
    defaults, constants = tui._load_defaults(Path("defaults"))
    assert defaults == _DEFAULTS
    assert constants == _CONSTANTS


def test_load_defaults_reports_missing_keys(monkeypatch):
    incomplete = {k: v for k, v in _DEFAULTS.items() if k != "k_tol_nd"}
    monkeypatch.setattr(tui, "load_json", _fake_load_json(incomplete, {}))
    with pytest.raises(ValueError, match="k_tol_nd, pi_nd"):
        tui._load_defaults(Path("defaults"))


def test_stage1_input_prefills_from_defaults(monkeypatch, capsys):
    """Pressing Enter at every prompt returns exactly the loaded defaults."""
    monkeypatch.setattr("builtins.input", lambda prompt="": "")
    data = tui._stage1_input(_DEFAULTS)
    v_min_m_s, v_max_m_s, n_pts = _DEFAULTS["velocity_range_m_s"]
    for key in (
        "b_m",
        "a_h_nd",
        "x_alpha_nd",
        "mass_kg",
        "inertia_alpha_kg_m2",
        "omega_heave_rad_s",
        "omega_alpha_rad_s",
        "rho_kg_m3",
    ):
        assert data[key] == _DEFAULTS[key]
    assert (data["v_min_m_s"], data["v_max_m_s"], data["n_pts"]) == (v_min_m_s, v_max_m_s, n_pts)


def test_solve_modes_passes_solver_settings(monkeypatch, capsys):
    captured = {}

    def fake_compute_flutter(**kwargs):
        captured.update(kwargs)
        return np.ones((2, len(kwargs["velocity_range_m_s"])), dtype=np.complex128)

    monkeypatch.setattr(tui, "compute_flutter", fake_compute_flutter)
    data = {key: _DEFAULTS[key] for key in _DEFAULTS if key in tui._DEFAULT_KEYS[:8]}
    solver = {"k_tol_nd": 1e-9, "max_iter_nd": 150, "pi_nd": 3.14159}
    tui._solve_modes(data, np.linspace(2.0, 80.0, 40), solver)
    assert {key: captured[key] for key in solver} == solver
