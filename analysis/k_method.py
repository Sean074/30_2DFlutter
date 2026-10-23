"""k-method (V-g) flutter analysis: required structural damping over a reduced-frequency sweep."""

from __future__ import annotations

from itertools import permutations

import numpy as np
from numpy.typing import NDArray

from analysis.flutter import (
    build_aero_matrix,
    check_flap_dofs,
    compute_structural_frequencies,
    theo_c,
    validate_velocity_range,
)

# The k sweep runs from K_MAX_MARGIN * b * omega_max / v_min down to
# K_MIN_MARGIN * b * omega_min / v_max, wide enough that every branch spans the velocity range
# (air added mass lowers frequencies, flutter usually raises the lower one).
_K_MAX_MARGIN_ND = 2.0
_K_MIN_MARGIN_ND = 0.25


def compute_k_branches(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    k_arr_nd: NDArray[np.float64],
    rho_kg_m3: float,
    b_m: float,
    a_h_nd: float,
    pi_nd: float = 3.141592653589793,
    c_hinge_nd: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """
    Solve the k-method eigenproblem at each reduced frequency and track the branches.

    With harmonic motion and structural damping K(1 + ig) the flutter equation becomes
    Z x = K^-1 (M + rho pi b^3 A(k) / k^2) x, Z = (1 + ig) / omega^2, so each eigenvalue gives
    omega = 1 / sqrt(real(Z)), V = b omega / k and g = imag(Z) / real(Z). Branches are matched
    from one k to the next by the nearest set of eigenvalues, starting in ascending frequency.

    Parameters
    ----------
    mass_matrix_kg : NDArray[np.float64], shape (n, n), n = 2 or 3
        Generalised structural mass matrix.
    stiffness_matrix_n_m : NDArray[np.float64], shape (n, n), n = 2 or 3
        Structural stiffness matrix in N/m (or N·m/rad on pitch diagonal).
    k_arr_nd : NDArray[np.float64], shape (nk,)
        Reduced frequencies, positive and strictly decreasing (i.e. increasing speed).
    rho_kg_m3 : float
        Air density in kg/m³.
    b_m : float
        Semi-chord length in metres.
    a_h_nd : float
        Elastic axis location aft of semi-chord mid-point (dimensionless).
    pi_nd : float
        Mathematical constant π, loaded from defaults/constants.json.
    c_hinge_nd : float or None
        Flap hinge location aft of mid-chord in semi-chords for 3x3 (heave, pitch, flap)
        matrices; None (default) for the two-DOF section.

    Returns
    -------
    tuple of three NDArray[np.float64], each shape (n_modes, nk)
        (velocity in m/s, frequency in rad/s, required structural damping g) per branch and k.
        Where real(Z) <= 0 the branch has no real frequency and all three are NaN.

    Raises
    ------
    ValueError
        If k_arr_nd is empty, not positive or not strictly decreasing, or c_hinge_nd does not
        match the matrix size.
    """
    check_flap_dofs(mass_matrix_kg, c_hinge_nd)
    if len(k_arr_nd) == 0:
        raise ValueError("k_arr_nd must not be empty")
    if np.any(k_arr_nd <= 0.0) or np.any(np.diff(k_arr_nd) >= 0.0):
        raise ValueError("k_arr_nd must be positive and strictly decreasing")

    n_modes = mass_matrix_kg.shape[0]
    z_arr = np.empty((n_modes, len(k_arr_nd)), dtype=np.complex128)
    for i, k_nd in enumerate(k_arr_nd):
        aero_matrix = build_aero_matrix(k_nd, theo_c(k_nd), b_m, a_h_nd, c_hinge_nd)
        z_col = np.linalg.eigvals(
            np.linalg.solve(
                stiffness_matrix_n_m,
                mass_matrix_kg + rho_kg_m3 * pi_nd * b_m**3 / k_nd**2 * aero_matrix,
            )
        )
        if i == 0:
            # Largest real(Z) is the lowest frequency: order as compute_structural_frequencies.
            z_arr[:, 0] = z_col[np.argsort(-np.real(z_col))]
            continue
        z_prev = z_arr[:, i - 1]
        order = min(
            permutations(range(n_modes)),
            key=lambda o: np.sum(np.abs(z_col[list(o)] - z_prev)),
        )
        z_arr[:, i] = z_col[list(order)]

    valid_arr = np.real(z_arr) > 0.0
    z_real_arr = np.where(valid_arr, np.real(z_arr), np.nan)
    omega_arr_rad_s = 1.0 / np.sqrt(z_real_arr)
    velocity_arr_m_s = b_m * omega_arr_rad_s / k_arr_nd[None, :]
    g_arr_nd = np.imag(z_arr) / z_real_arr
    return velocity_arr_m_s, omega_arr_rad_s, g_arr_nd


def compute_flutter_k(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    velocity_range_m_s: NDArray[np.float64],
    rho_kg_m3: float,
    b_m: float,
    a_h_nd: float,
    n_k_nd: int = 4000,
    pi_nd: float = 3.141592653589793,
    c_hinge_nd: float | None = None,
) -> NDArray[np.complex128]:
    """
    Compute k-method flutter results for every mode, on a velocity sweep.

    The reduced frequency is swept downwards on a geometric grid (`compute_k_branches`) and each
    branch's frequency and required damping g are interpolated linearly in V onto
    velocity_range_m_s. A branch is used from its lowest speed up to its highest; beyond that
    (where V(k) folds back, typically towards static divergence) and below its lowest speed the
    result is NaN.

    The output has the shape and encoding of `compute_flutter`: p = k (g / 2 + i), with
    k = b omega / V, so `damping_g_nd`, `find_flutter_speed` and the plots apply unchanged. g is
    the damping the structure would need for neutral oscillation; it equals the p-k damping only
    at g = 0, which is why the two methods give the same flutter point but different V-g curves.

    Parameters
    ----------
    mass_matrix_kg : NDArray[np.float64], shape (n, n), n = 2 or 3
        Generalised structural mass matrix.
    stiffness_matrix_n_m : NDArray[np.float64], shape (n, n), n = 2 or 3
        Structural stiffness matrix in N/m (or N·m/rad on pitch diagonal).
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep points in m/s. Must be positive and strictly increasing.
    rho_kg_m3 : float
        Air density in kg/m³.
    b_m : float
        Semi-chord length in metres.
    a_h_nd : float
        Elastic axis location aft of semi-chord mid-point (dimensionless).
    n_k_nd : int
        Number of reduced frequencies in the sweep.
    pi_nd : float
        Mathematical constant π, loaded from defaults/constants.json.
    c_hinge_nd : float or None
        Flap hinge location aft of mid-chord in semi-chords for 3x3 (heave, pitch, flap)
        matrices; None (default) for the two-DOF section.

    Returns
    -------
    NDArray[np.complex128], shape (n_modes, nv)
        Equivalent p per mode and speed, modes ordered by ascending structural frequency; NaN
        where a branch does not reach that speed.

    Raises
    ------
    ValueError
        If the velocity sweep is invalid, n_k_nd < 2, the matrices are not positive definite,
        or c_hinge_nd does not match the matrix size.
    """
    validate_velocity_range(velocity_range_m_s)
    if n_k_nd < 2:
        raise ValueError(f"n_k_nd must be at least 2; got {n_k_nd}")

    omega_arr_rad_s = compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)
    k_max_nd = _K_MAX_MARGIN_ND * b_m * omega_arr_rad_s[-1] / velocity_range_m_s[0]
    k_min_nd = _K_MIN_MARGIN_ND * b_m * omega_arr_rad_s[0] / velocity_range_m_s[-1]
    k_arr_nd = np.geomspace(k_max_nd, k_min_nd, n_k_nd)
    v_branch_m_s, omega_branch_rad_s, g_branch_nd = compute_k_branches(
        mass_matrix_kg,
        stiffness_matrix_n_m,
        k_arr_nd,
        rho_kg_m3,
        b_m,
        a_h_nd,
        pi_nd,
        c_hinge_nd,
    )

    # Both parts NaN: np.nan alone would give nan + 0j, i.e. a frequency of zero.
    p_arr_nd = np.full(
        (len(omega_arr_rad_s), len(velocity_range_m_s)),
        complex(np.nan, np.nan),
        np.complex128,
    )
    for mode in range(len(omega_arr_rad_s)):
        v_arr_m_s = v_branch_m_s[mode]
        # Keep the leading stretch where V rises with falling k (and is finite).
        rising_arr = np.isfinite(v_arr_m_s[1:]) & (np.diff(v_arr_m_s) > 0.0)
        n_keep = 1 + (np.argmin(rising_arr) if not np.all(rising_arr) else len(rising_arr))
        if n_keep < 2 or not np.isfinite(v_arr_m_s[0]):
            continue
        v_keep_m_s = v_arr_m_s[:n_keep]
        inside_arr = (velocity_range_m_s >= v_keep_m_s[0]) & (velocity_range_m_s <= v_keep_m_s[-1])
        v_out_m_s = velocity_range_m_s[inside_arr]
        omega_out_rad_s = np.interp(v_out_m_s, v_keep_m_s, omega_branch_rad_s[mode, :n_keep])
        g_out_nd = np.interp(v_out_m_s, v_keep_m_s, g_branch_nd[mode, :n_keep])
        k_out_nd = b_m * omega_out_rad_s / v_out_m_s
        p_arr_nd[mode, inside_arr] = k_out_nd * (0.5 * g_out_nd + 1j)
    return p_arr_nd
