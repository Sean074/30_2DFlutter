"""Tests for flutter2d/analysis/flutter.py."""

from __future__ import annotations

import numpy as np
import pytest

from flutter2d.analysis.flutter import (
    build_aero_matrix,
    build_mass_matrix,
    build_stiffness_matrix,
    check_flap_dofs,
    compute_flutter,
    compute_structural_frequencies,
    damping_g_nd,
    find_flutter_speed,
    theo_c,
    validate_velocity_range,
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
    """TUI default section: flutter at 19.5582 m/s, from an independent k-method (g = 0) solve."""
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
    np.testing.assert_allclose(v_flutter_m_s, 19.5582, rtol=1e-4)


def test_find_flutter_speed_interpolates_between_grid_points():
    """gamma = -0.1 at 10 m/s and +0.3 at 20 m/s crosses zero a quarter of the way: 12.5 m/s."""
    b_m = 0.5
    velocity_range_m_s = np.array([5.0, 10.0, 20.0, 30.0])
    k_arr_nd = np.array([2.0, 1.0, 0.6, 0.4])
    gamma_arr_nd = np.array([-0.2, -0.1, 0.3, 0.5])
    p_arr_nd = k_arr_nd * (gamma_arr_nd + 1j)

    v_flutter_m_s, omega_flutter_rad_s = find_flutter_speed(p_arr_nd, velocity_range_m_s, b_m)

    np.testing.assert_allclose(v_flutter_m_s, 12.5, rtol=1e-4)
    # frequency (k V / b) is 20 rad/s at 10 m/s and 24 rad/s at 20 m/s
    np.testing.assert_allclose(omega_flutter_rad_s, 21.0, rtol=1e-4)


def test_find_flutter_speed_returns_none_without_crossing():
    velocity_range_m_s = np.array([5.0, 10.0, 20.0])
    p_arr_nd = np.array([1.0, 0.8, 0.6]) * (-0.05 + 1j)
    assert find_flutter_speed(p_arr_nd, velocity_range_m_s, 0.5) is None


def test_compute_flutter_hodges_pierce_textbook_case():
    """
    Textbook two-DOF section (Hodges & Pierce, Introduction to Structural Dynamics and
    Aeroelasticity, ch. 5): a = -1/5, e = -1/10 (x_a = 1/10), mu = 20, r^2 = 6/25,
    sigma = w_h / w_a = 2/5. Non-dimensionalised with b = 1 m, w_a = 1 rad/s, rho = 1 kg/m3.

    With the R. T. Jones approximation of C(k) the flutter point is V_F / (b w_a) = 2.17036,
    w_F / w_a = 0.6443, confirmed by an independent k-method (g = 0) solve. With exact
    Theodorsen the same solve gives 2.184, consistent with the book's ~2.2; the 0.6 %
    difference is the Jones approximation.
    """
    b_m = 1.0
    rho_kg_m3 = 1.0
    omega_alpha_rad_s = 1.0
    mass_kg = 20.0 * np.pi * rho_kg_m3 * b_m**2
    inertia_alpha_kg_m2 = 6.0 / 25.0 * mass_kg * b_m**2
    mass_matrix_kg = build_mass_matrix(mass_kg, inertia_alpha_kg_m2, b_m, 0.1)
    stiffness_matrix_n_m = build_stiffness_matrix(
        (0.4 * omega_alpha_rad_s) ** 2 * mass_kg,
        omega_alpha_rad_s**2 * inertia_alpha_kg_m2,
    )
    velocity_range_m_s = np.linspace(0.01, 3.0, 600)

    p_arr = compute_flutter(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, rho_kg_m3, b_m, a_h_nd=-0.2
    )

    assert find_flutter_speed(p_arr[0], velocity_range_m_s, b_m) is None
    v_flutter_m_s, omega_flutter_rad_s = find_flutter_speed(p_arr[1], velocity_range_m_s, b_m)
    np.testing.assert_allclose(v_flutter_m_s / (b_m * omega_alpha_rad_s), 2.17036, rtol=1e-4)
    np.testing.assert_allclose(omega_flutter_rad_s / omega_alpha_rad_s, 0.6443, rtol=1e-3)


@pytest.mark.parametrize(
    "velocity_range_m_s", [np.array([]), np.array([0.0, 1.0]), np.array([2.0, 1.0])]
)
def test_validate_velocity_range_rejects_bad_sweeps(velocity_range_m_s):
    with pytest.raises(ValueError, match="velocity_range_m_s"):
        validate_velocity_range(velocity_range_m_s)


def test_damping_g_nd_is_twice_decay_rate_ratio():
    """p = k (gamma + i) with gamma = -1.5 % gives g = -0.03."""
    p_arr_nd = np.array([0.4 * (-0.015 + 1j), 2.0 * (0.01 + 1j)])
    np.testing.assert_allclose(damping_g_nd(p_arr_nd), [-0.03, 0.02], rtol=1e-4)


# Shipped defaults plus the flap of data/inputs/flap_example.json.
_FLAP_SECTION = {"c_hinge_nd": 0.6, "x_beta_nd": 0.0125, "inertia_beta_kg_m2": 0.0156}


def _flap_matrices(omega_beta_rad_s: float) -> tuple[np.ndarray, np.ndarray]:
    mass_matrix_kg = build_mass_matrix(10.0, 0.5, 0.5, 0.1, a_h_nd=-0.2, **_FLAP_SECTION)
    stiffness_matrix_n_m = build_stiffness_matrix(
        10.0**2 * 10.0, 25.0**2 * 0.5, omega_beta_rad_s**2 * _FLAP_SECTION["inertia_beta_kg_m2"]
    )
    return mass_matrix_kg, stiffness_matrix_n_m


@pytest.mark.parametrize("k_nd", [0.05, 0.4, 1.5])
def test_aero_matrix_flap_hinged_at_leading_edge_equals_pitch_about_it(k_nd):
    """With c = a = -1 the flap rotation is the whole-chord pitch: rows and columns coincide."""
    aero_matrix = build_aero_matrix(k_nd, theo_c(k_nd), 0.7, -1.0, -1.0 + 1e-12)
    np.testing.assert_allclose(aero_matrix[:, 2], aero_matrix[:, 1], rtol=1e-4, atol=1e-9)
    np.testing.assert_allclose(aero_matrix[2, :], aero_matrix[1, :], rtol=1e-4, atol=1e-9)


def test_aero_matrix_flap_steady_lift_matches_thin_airfoil_theory():
    """Steady lift per flap angle: dCl/dbeta = 2 (pi - theta_h + sin theta_h), cos theta_h = -c."""
    c_hinge_nd = 0.45
    theta_h_rad = np.arccos(-c_hinge_nd)
    aero_matrix = build_aero_matrix(1e-9, 1.0 + 0j, 0.8, -0.3, c_hinge_nd)
    # Row 0 is the downward force over rho pi b V^2: L = rho V^2 b Cl = -rho pi b V^2 a13 beta.
    np.testing.assert_allclose(
        -np.pi * np.real(aero_matrix[0, 2]),
        2.0 * (np.pi - theta_h_rad + np.sin(theta_h_rad)),
        rtol=1e-4,
    )


def test_aero_matrix_flap_added_mass_is_symmetric_positive_definite():
    """The k^2 (apparent mass) part of A with C = 0 is symmetric and positive definite."""

    def aero(k_nd: float) -> np.ndarray:
        return build_aero_matrix(k_nd, 0.0, 1.0, -0.3, 0.5)

    # A = X k^2 + Y i k + Z with C = 0, so X = (A(2) - 2 A(1) + A(0)) / 2.
    added_mass = np.real(aero(2.0) - 2.0 * aero(1.0) + aero(0.0)) / 2.0
    np.testing.assert_allclose(added_mass, added_mass.T, atol=1e-12)
    assert np.all(np.linalg.eigvalsh(added_mass) > 0.0)


def test_aero_matrix_rejects_hinge_outside_chord():
    with pytest.raises(ValueError, match="hinge"):
        build_aero_matrix(0.5, theo_c(0.5), 0.5, -0.2, 1.0)


def test_build_mass_matrix_with_flap():
    mass_matrix_kg = build_mass_matrix(10.0, 0.5, 0.5, 0.1, a_h_nd=-0.2, **_FLAP_SECTION)
    s_beta_kg_m = 10.0 * 0.5 * 0.0125
    pitch_flap_kg_m2 = 0.0156 + 0.5 * (0.6 + 0.2) * s_beta_kg_m
    np.testing.assert_allclose(
        mass_matrix_kg,
        [
            [10.0, 0.5, s_beta_kg_m],
            [0.5, 0.5, pitch_flap_kg_m2],
            [s_beta_kg_m, pitch_flap_kg_m2, 0.0156],
        ],
        rtol=1e-4,
    )


def test_build_mass_matrix_rejects_partial_flap():
    with pytest.raises(ValueError, match="flap"):
        build_mass_matrix(10.0, 0.5, 0.5, 0.1, x_beta_nd=0.01)


def test_build_stiffness_matrix_with_flap_is_diagonal():
    np.testing.assert_allclose(
        build_stiffness_matrix(1000.0, 300.0, 14.0), np.diag([1000.0, 300.0, 14.0]), rtol=1e-4
    )


def test_check_flap_dofs_matches_hinge_to_matrix_size():
    with pytest.raises(ValueError, match="hinge"):
        check_flap_dofs(np.eye(3), None)
    with pytest.raises(ValueError, match="no flap"):
        check_flap_dofs(np.eye(2), 0.6)


def test_compute_flutter_stiff_flap_recovers_two_dof_speed():
    """A flap at 1000 rad/s barely moves: the 2-DOF flutter speed 19.5582 m/s is recovered."""
    mass_matrix_kg, stiffness_matrix_n_m = _flap_matrices(1000.0)
    velocity_range_m_s = np.linspace(0.5, 40.0, 400)
    p_arr = compute_flutter(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.225, 0.5, -0.2, c_hinge_nd=0.6
    )
    v_flutter_m_s, _ = find_flutter_speed(p_arr[1], velocity_range_m_s, 0.5)
    np.testing.assert_allclose(v_flutter_m_s, 19.5582, rtol=1e-4)


def test_compute_flutter_three_dof_flap_example():
    """
    data/inputs/flap_example.json: flutter of mode 2 at 20.4965 m/s, 16.859 rad/s (the k-method
    gives the same point, see tests/test_k_method.py). The heavily damped heave mode stops
    oscillating near 40 m/s; its branch ends there (NaN) instead of stopping the sweep.
    """
    mass_matrix_kg, stiffness_matrix_n_m = _flap_matrices(30.0)
    velocity_range_m_s = np.linspace(0.5, 60.0, 600)
    p_arr = compute_flutter(
        mass_matrix_kg, stiffness_matrix_n_m, velocity_range_m_s, 1.225, 0.5, -0.2, c_hinge_nd=0.6
    )
    assert p_arr.shape == (3, 600)
    v_flutter_m_s, omega_flutter_rad_s = find_flutter_speed(p_arr[1], velocity_range_m_s, 0.5)
    np.testing.assert_allclose(v_flutter_m_s, 20.4965, rtol=1e-4)
    np.testing.assert_allclose(omega_flutter_rad_s, 16.859, rtol=1e-4)
    ended_arr = np.isnan(p_arr[0])
    assert ended_arr[-1] and not ended_arr[0]
    assert np.all(np.isfinite(p_arr[1:]))
