"""Tests for ui/tui.py."""
from __future__ import annotations

import re
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
    "method": "k",
    "n_k_nd": 1234,
    "flap": False,
    "c_hinge_nd": 0.55,
    "x_beta_nd": 0.02,
    "inertia_beta_kg_m2": 0.03,
    "omega_beta_rad_s": 40.0,
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
        "method",
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
    data["method"] = "pk"
    data["flap"] = False
    solver = {"k_tol_nd": 1e-9, "max_iter_nd": 150, "pi_nd": 3.14159}
    tui._solve_modes(data, np.linspace(2.0, 80.0, 40), solver | {"n_k_nd": 1234})
    assert {key: captured[key] for key in solver} == solver


def test_solve_modes_dispatches_k_method(monkeypatch, capsys):
    captured = {}

    def fake_compute_flutter_k(**kwargs):
        captured.update(kwargs)
        return np.ones((2, len(kwargs["velocity_range_m_s"])), dtype=np.complex128)

    def fail_compute_flutter(**kwargs):
        raise AssertionError("p-k solver called for method 'k'")

    monkeypatch.setattr(tui, "compute_flutter_k", fake_compute_flutter_k)
    monkeypatch.setattr(tui, "compute_flutter", fail_compute_flutter)
    data = {key: _DEFAULTS[key] for key in _DEFAULTS if key in tui._DEFAULT_KEYS[:8]}
    data["method"] = "k"
    data["flap"] = False
    solver = {"k_tol_nd": 1e-9, "max_iter_nd": 150, "pi_nd": 3.14159, "n_k_nd": 1234}
    tui._solve_modes(data, np.linspace(2.0, 80.0, 40), solver)
    assert (captured["n_k_nd"], captured["pi_nd"]) == (1234, 3.14159)


def test_stage1_input_rejects_unknown_method(monkeypatch, capsys):
    """An unknown method is re-prompted; 'PK' is accepted case-insensitively."""
    answers = iter(["xyz", "PK"])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers, ""))
    data = tui._stage1_input(_DEFAULTS)
    assert data["method"] == "pk"
    assert "[E] INVALID INPUT" in capsys.readouterr().out


def test_load_defaults_rejects_unknown_method(monkeypatch):
    monkeypatch.setattr(tui, "load_json", _fake_load_json(_DEFAULTS | {"method": "x"}, _CONSTANTS))
    with pytest.raises(ValueError, match="unknown method 'x'"):
        tui._load_defaults(Path("defaults"))


def test_print_results_table_shows_dash_where_branch_missing(capsys):
    p_arr = np.array([0.5 * (-0.02 + 1j), np.nan])
    tui._print_results_table({"MODE 1": p_arr}, np.array([5.0, 10.0]), 0.5)
    last_row = capsys.readouterr().out.strip().splitlines()[-1]
    assert last_row.split() == ["10.000", "-", "-"]


def test_stage1_input_flap_prompts_when_included(monkeypatch, capsys):
    """Answering Y at INCLUDE FLAP prompts for the four flap inputs, pre-filled from defaults."""
    monkeypatch.setattr(
        "builtins.input", lambda prompt="": "y" if "INCLUDE FLAP" in prompt else ""
    )
    data = tui._stage1_input(_DEFAULTS)
    assert data["flap"] is True
    for key in tui._FLAP_KEYS:
        assert data[key] == _DEFAULTS[key]


def test_stage1_input_rejects_hinge_at_trailing_edge(monkeypatch, capsys):
    answers = iter(["y", "1.0"])

    def fake_input(prompt=""):
        if "INCLUDE FLAP" in prompt or "HINGE" in prompt:
            return next(answers, "")
        return ""

    monkeypatch.setattr("builtins.input", fake_input)
    data = tui._stage1_input(_DEFAULTS)
    assert data["c_hinge_nd"] == _DEFAULTS["c_hinge_nd"]
    assert "must be < 1.0" in capsys.readouterr().out


def test_flap_from_input_all_or_none():
    assert tui._flap_from_input({"b_m": 0.5}) == {"flap": False}
    full = {key: _DEFAULTS[key] for key in tui._FLAP_KEYS}
    assert tui._flap_from_input(full) == {"flap": True} | full
    with pytest.raises(ValueError, match="omega_beta_rad_s"):
        tui._flap_from_input({key: 1.0 for key in tui._FLAP_KEYS[:3]})


def test_solve_modes_builds_three_dof_system_with_flap(monkeypatch, capsys):
    captured = {}

    def fake_compute_flutter(**kwargs):
        captured.update(kwargs)
        return np.ones((3, len(kwargs["velocity_range_m_s"])), dtype=np.complex128)

    monkeypatch.setattr(tui, "compute_flutter", fake_compute_flutter)
    data = {key: _DEFAULTS[key] for key in _DEFAULTS if key in tui._DEFAULT_KEYS[:8]}
    data |= {key: _DEFAULTS[key] for key in tui._FLAP_KEYS}
    data |= {"method": "pk", "flap": True}
    solver = {"k_tol_nd": 1e-9, "max_iter_nd": 150, "pi_nd": 3.14159, "n_k_nd": 1234}
    modes = tui._solve_modes(data, np.linspace(2.0, 80.0, 40), solver)
    assert len(modes) == 3
    assert captured["mass_matrix_kg"].shape == (3, 3)
    assert captured["c_hinge_nd"] == _DEFAULTS["c_hinge_nd"]
    np.testing.assert_allclose(
        captured["stiffness_matrix_n_m"][2, 2], 40.0**2 * 0.03, rtol=1e-4
    )


def test_print_results_table_three_modes_fits_80_columns(capsys):
    p_arr = np.array([0.5 * (-0.02 + 1j), 0.4 * (-0.01 + 1j)])
    tui._print_results_table({f"MODE {i}": p_arr for i in (1, 2, 3)}, np.array([5.0, 10.0]), 0.5)
    lines = re.sub(r"\x1b\[[0-9;]*m", "", capsys.readouterr().out).splitlines()
    assert max(len(line) for line in lines) <= 80


def test_stage4_summary_warns_where_pk_mode_ends(capsys):
    p_arr = np.array([0.5 * (-0.02 + 1j), complex(np.nan, np.nan)])
    tui._stage4_summary({"MODE 1": p_arr}, np.array([5.0, 10.0]), 0.5, None, "pk")
    assert "MODE 1 ends at V = 5.0 m/s (p-k cannot follow it)" in capsys.readouterr().out
