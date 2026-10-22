"""Tests for analysis/flutter.py."""
from __future__ import annotations

import numpy as np
import pytest

from analysis.flutter import (
    build_mass_matrix,
    build_stiffness_matrix,
    compute_flutter,
    compute_structural_frequencies,
    theo_c,
)


def test_theo_c_approaches_half_at_high_k():
    """C(k) → 0.5 as k → ∞ (Jones two-pole limit matches exact Theodorsen)."""
    result = theo_c(1000.0)
    np.testing.assert_allclose(abs(result), 0.5, atol=0.01)


def test_build_mass_matrix_shape_and_symmetry():
    m = build_mass_matrix(mass_kg=10.4, inertia_alpha_kg_m2=1.12, b_m=0.5, x_alpha_nd=-0.5)
    assert m.shape == (2, 2)
    np.testing.assert_allclose(m[0, 1], m[1, 0])


def test_build_stiffness_matrix_diagonal():
    k = build_stiffness_matrix(k_heave_n_m=8000.0, k_alpha_nm_rad=757.12)
    assert k[0, 1] == 0.0
    assert k[1, 0] == 0.0
    np.testing.assert_allclose(k[0, 0], 8000.0)
    np.testing.assert_allclose(k[1, 1], 757.12)


def test_compute_flutter_returns_correct_shape():
    mass_kg = 10.4
    inertia_alpha_kg_m2 = 1.12
    b_m = 0.5
    x_alpha_nd = -0.5
    omega_alpha_rad_s = 26.0
    omega_heave_rad_s = 20.0

    mass_matrix_kg = build_mass_matrix(mass_kg, inertia_alpha_kg_m2, b_m, x_alpha_nd)
    stiffness_matrix_n_m = build_stiffness_matrix(
        k_heave_n_m=omega_heave_rad_s**2 * mass_kg,
        k_alpha_nm_rad=omega_alpha_rad_s**2 * inertia_alpha_kg_m2,
    )
    velocity_range_m_s = np.linspace(0.01, 50.0, 50)

    p_arr = compute_flutter(
        mass_matrix_kg=mass_matrix_kg,
        stiffness_matrix_n_m=stiffness_matrix_n_m,
        velocity_range_m_s=velocity_range_m_s,
        rho_kg_m3=1.21,
        b_m=b_m,
        a_h_nd=0.25,
    )
    assert p_arr.shape == (2, 50)
    assert p_arr.dtype == np.complex128


def test_compute_flutter_raises_on_empty_velocity():
    mass_matrix_kg = build_mass_matrix(10.0, 1.0, 0.5, -0.5)
    stiffness_matrix_n_m = build_stiffness_matrix(4000.0, 676.0)
    with pytest.raises(ValueError, match="empty"):
        compute_flutter(
            mass_matrix_kg=mass_matrix_kg,
            stiffness_matrix_n_m=stiffness_matrix_n_m,
            velocity_range_m_s=np.array([]),
            rho_kg_m3=1.225,
            b_m=0.5,
            a_h_nd=0.0,
        )


def test_compute_structural_frequencies_uncoupled():
    """With the CG on the elastic axis the modes are the uncoupled heave and pitch frequencies."""
    mass_matrix_kg = build_mass_matrix(10.0, 0.5, 0.5, 0.0)
    stiffness_matrix_n_m = build_stiffness_matrix(25.0**2 * 10.0, 10.0**2 * 0.5)
    omega_arr_rad_s = compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)
    np.testing.assert_allclose(omega_arr_rad_s, [10.0, 25.0], rtol=1e-4)


def test_compute_structural_frequencies_raises_on_singular_mass():
    """x_a^2 m b^2 > I_a makes the mass matrix indefinite."""
    mass_matrix_kg = build_mass_matrix(10.0, 0.1, 0.5, 0.5)
    stiffness_matrix_n_m = build_stiffness_matrix(4000.0, 676.0)
    with pytest.raises(ValueError, match="mass matrix"):
        compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)


def test_compute_flutter_raises_on_non_increasing_velocity():
    mass_matrix_kg = build_mass_matrix(10.0, 1.0, 0.5, -0.5)
    stiffness_matrix_n_m = build_stiffness_matrix(4000.0, 676.0)
    with pytest.raises(ValueError, match="increasing"):
        compute_flutter(
            mass_matrix_kg=mass_matrix_kg,
            stiffness_matrix_n_m=stiffness_matrix_n_m,
            velocity_range_m_s=np.array([10.0, 5.0]),
            rho_kg_m3=1.225,
            b_m=0.5,
            a_h_nd=0.0,
        )


def test_compute_flutter_keeps_modes_distinct():
    """Regression for issue #2: input whose uncoupled frequencies are 20 and 21 rad/s."""
    mass_kg = 10.4
    inertia_alpha_kg_m2 = 1.12
    b_m = 0.5
    mass_matrix_kg = build_mass_matrix(mass_kg, inertia_alpha_kg_m2, b_m, -0.5)
    stiffness_matrix_n_m = build_stiffness_matrix(20.0**2 * mass_kg, 21.0**2 * inertia_alpha_kg_m2)
    velocity_range_m_s = np.linspace(0.5, 50.0, 100)

    p_arr = compute_flutter(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.21, b_m, a_h_nd=0.25
    )
    freq_arr_rad_s = np.imag(p_arr) * velocity_range_m_s / b_m
    assert np.all(np.abs(freq_arr_rad_s[1] - freq_arr_rad_s[0]) > 20.0)
    omega_arr_rad_s = compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)
    # Air added mass lowers the low-speed frequencies a few percent below the in-vacuo values.
    np.testing.assert_allclose(freq_arr_rad_s[:, 0], omega_arr_rad_s, rtol=0.1)


def test_compute_flutter_detects_flutter_speed():
    """TUI default section: flutter at 19.3934 m/s, from an independent k-method (g = 0) solve."""
    mass_kg = 10.0
    inertia_alpha_kg_m2 = 0.5
    b_m = 0.5
    mass_matrix_kg = build_mass_matrix(mass_kg, inertia_alpha_kg_m2, b_m, 0.1)
    stiffness_matrix_n_m = build_stiffness_matrix(10.0**2 * mass_kg, 25.0**2 * inertia_alpha_kg_m2)
    velocity_range_m_s = np.linspace(0.1, 40.0, 400)

    p_arr = compute_flutter(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.225, b_m, a_h_nd=-0.2
    )
    g_arr_nd = np.real(p_arr) / np.imag(p_arr)
    assert np.all(g_arr_nd[0, 1:] < 0.0)
    j = np.flatnonzero((g_arr_nd[1, :-1] < 0.0) & (g_arr_nd[1, 1:] >= 0.0))[0]
    v1, v2 = velocity_range_m_s[j], velocity_range_m_s[j + 1]
    g1, g2 = g_arr_nd[1, j], g_arr_nd[1, j + 1]
    v_flutter_m_s = v1 - g1 * (v2 - v1) / (g2 - g1)
    np.testing.assert_allclose(v_flutter_m_s, 19.3934, rtol=1e-4)
