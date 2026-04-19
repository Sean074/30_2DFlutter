"""Tests for analysis/flutter.py."""
from __future__ import annotations

import numpy as np
import pytest

from analysis.flutter import (
    build_mass_matrix,
    build_stiffness_matrix,
    compute_flutter,
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
        start_omega_rad_s=omega_heave_rad_s,
    )
    assert p_arr.shape == (50,)
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
            start_omega_rad_s=26.0,
        )
