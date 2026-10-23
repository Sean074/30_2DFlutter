"""Tests for analysis/k_method.py."""

from __future__ import annotations

import numpy as np
import pytest

from analysis.flutter import (
    build_mass_matrix,
    build_stiffness_matrix,
    compute_flutter,
    compute_structural_frequencies,
    find_flutter_speed,
)
from analysis.k_method import compute_flutter_k, compute_k_branches


def _hodges_pierce() -> tuple[np.ndarray, np.ndarray]:
    """Hodges & Pierce ch. 5 section (see tests/test_flutter.py), b = 1 m, w_a = 1 rad/s."""
    mass_kg = 20.0 * np.pi
    inertia_alpha_kg_m2 = 6.0 / 25.0 * mass_kg
    return (
        build_mass_matrix(mass_kg, inertia_alpha_kg_m2, 1.0, 0.1),
        build_stiffness_matrix(0.4**2 * mass_kg, inertia_alpha_kg_m2),
    )


def _shipped_defaults() -> tuple[np.ndarray, np.ndarray]:
    """defaults/default_flutter.json section: b = 0.5 m, a = -0.2, x_a = 0.1."""
    return (
        build_mass_matrix(10.0, 0.5, 0.5, 0.1),
        build_stiffness_matrix(10.0**2 * 10.0, 25.0**2 * 0.5),
    )


def test_k_method_hodges_pierce_matches_pk_flutter_point():
    """Same flutter point as the p-k solver: V_F / (b w_a) = 2.17036, w_F / w_a = 0.6443."""
    mass_matrix_kg, stiffness_matrix_n_m = _hodges_pierce()
    velocity_range_m_s = np.linspace(0.01, 3.0, 600)

    p_arr = compute_flutter_k(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.0, 1.0, a_h_nd=-0.2
    )

    assert p_arr.shape == (2, 600)
    assert find_flutter_speed(p_arr[0], velocity_range_m_s, 1.0) is None
    v_flutter_m_s, omega_flutter_rad_s = find_flutter_speed(p_arr[1], velocity_range_m_s, 1.0)
    np.testing.assert_allclose(v_flutter_m_s, 2.17036, rtol=1e-4)
    np.testing.assert_allclose(omega_flutter_rad_s, 0.6443, rtol=1e-3)


def test_k_method_shipped_defaults_matches_pk_flutter_speed():
    mass_matrix_kg, stiffness_matrix_n_m = _shipped_defaults()
    velocity_range_m_s = np.linspace(0.1, 40.0, 400)
    args = (mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.225, 0.5, -0.2)

    v_k_m_s, _ = find_flutter_speed(compute_flutter_k(*args)[1], velocity_range_m_s, 0.5)
    v_pk_m_s, _ = find_flutter_speed(compute_flutter(*args)[1], velocity_range_m_s, 0.5)

    np.testing.assert_allclose(v_k_m_s, 19.5582, rtol=1e-4)
    np.testing.assert_allclose(v_k_m_s, v_pk_m_s, rtol=1e-4)


def test_k_method_no_flutter_with_cg_ahead_of_elastic_axis():
    """data/inputs/default.json section (x_a < 0): no g = 0 crossing up to 50 m/s."""
    mass_matrix_kg = build_mass_matrix(10.4, 1.12, 0.5, -0.5)
    stiffness_matrix_n_m = build_stiffness_matrix(20.0**2 * 10.4, 21.0**2 * 1.12)
    velocity_range_m_s = np.linspace(0.5, 50.0, 200)

    p_arr = compute_flutter_k(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.21, 0.5, a_h_nd=0.25
    )

    for p_mode_arr in p_arr:
        assert find_flutter_speed(p_mode_arr, velocity_range_m_s, 0.5) is None


def test_k_method_branch_is_nan_beyond_its_fold():
    """Shipped defaults: the lower branch's V(k) peaks near 25 m/s, so 40 m/s is not reached."""
    mass_matrix_kg, stiffness_matrix_n_m = _shipped_defaults()
    velocity_range_m_s = np.array([5.0, 10.0, 40.0])

    p_arr = compute_flutter_k(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.225, 0.5, a_h_nd=-0.2
    )

    assert np.all(np.isfinite(p_arr[:, :2]))
    assert np.isnan(p_arr[0, 2].real) and np.isnan(p_arr[0, 2].imag)
    assert np.isfinite(p_arr[1, 2])


def test_k_branches_start_in_structural_frequency_order_without_air():
    """With rho = 0, Z = 1 / omega^2 exactly: coupled frequencies, ascending, and g = 0."""
    mass_matrix_kg, stiffness_matrix_n_m = _shipped_defaults()
    k_arr_nd = np.array([2.0, 1.0])

    v_arr_m_s, omega_arr_rad_s, g_arr_nd = compute_k_branches(
        mass_matrix_kg, stiffness_matrix_n_m, k_arr_nd, 0.0, 0.5, -0.2
    )

    expected_rad_s = compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)
    np.testing.assert_allclose(omega_arr_rad_s[:, 0], expected_rad_s, rtol=1e-4)
    np.testing.assert_allclose(v_arr_m_s[:, 1], 0.5 * expected_rad_s / 1.0, rtol=1e-4)
    np.testing.assert_allclose(g_arr_nd, 0.0, atol=1e-12)


@pytest.mark.parametrize("k_arr_nd", [np.array([]), np.array([1.0, 2.0]), np.array([1.0, 0.0])])
def test_k_branches_reject_bad_k_sweep(k_arr_nd):
    mass_matrix_kg, stiffness_matrix_n_m = _shipped_defaults()
    with pytest.raises(ValueError, match="k_arr_nd"):
        compute_k_branches(mass_matrix_kg, stiffness_matrix_n_m, k_arr_nd, 1.225, 0.5, -0.2)


def test_compute_flutter_k_rejects_too_few_k_points():
    mass_matrix_kg, stiffness_matrix_n_m = _shipped_defaults()
    with pytest.raises(ValueError, match="n_k_nd"):
        compute_flutter_k(
            mass_matrix_kg, stiffness_matrix_n_m, np.array([1.0, 2.0]), 1.225, 0.5, -0.2, n_k_nd=1
        )


def test_k_method_three_dof_flap_example_matches_pk():
    """data/inputs/flap_example.json: same flutter point as p-k, 20.4965 m/s at 16.859 rad/s."""
    mass_matrix_kg = build_mass_matrix(
        10.0,
        0.5,
        0.5,
        0.1,
        x_beta_nd=0.0125,
        inertia_beta_kg_m2=0.0156,
        c_hinge_nd=0.6,
        a_h_nd=-0.2,
    )
    stiffness_matrix_n_m = build_stiffness_matrix(1000.0, 312.5, 30.0**2 * 0.0156)
    velocity_range_m_s = np.linspace(0.5, 60.0, 600)

    p_arr = compute_flutter_k(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.225, 0.5, -0.2, c_hinge_nd=0.6
    )

    assert p_arr.shape == (3, 600)
    v_flutter_m_s, omega_flutter_rad_s = find_flutter_speed(p_arr[1], velocity_range_m_s, 0.5)
    np.testing.assert_allclose(v_flutter_m_s, 20.4965, rtol=1e-4)
    np.testing.assert_allclose(omega_flutter_rad_s, 16.859, rtol=1e-4)


def test_k_branches_reject_flap_matrices_without_hinge():
    mass_matrix_kg, stiffness_matrix_n_m = np.eye(3), np.eye(3)
    with pytest.raises(ValueError, match="hinge"):
        compute_k_branches(
            mass_matrix_kg, stiffness_matrix_n_m, np.array([1.0, 0.5]), 1.2, 0.5, 0.0
        )
