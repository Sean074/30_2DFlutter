"""p-k flutter analysis: matrix builders and secant-method solver."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def theo_c(k_nd: float) -> complex:
    """
    Theodorsen aerodynamic transfer function (Jones two-pole approximation).

    Parameters
    ----------
    k_nd : float
        Reduced frequency (dimensionless), k = b * omega / V.

    Returns
    -------
    complex
        Theodorsen function C(k).
    """
    return 1 - 0.165 / (1 - 0.041 / k_nd * 1j) - 0.335 / (1 - 0.32 / k_nd * 1j)


def build_aero_matrix(
    k_nd: float,
    c_nd: complex,
    b_m: float,
    a_h_nd: float,
) -> NDArray[np.complex128]:
    """
    Build the 2D oscillatory aerodynamic influence matrix.

    Parameters
    ----------
    k_nd : float
        Reduced frequency (dimensionless).
    c_nd : complex
        Theodorsen function value at k_nd.
    b_m : float
        Semi-chord length in metres.
    a_h_nd : float
        Elastic axis location aft of semi-chord mid-point (dimensionless, +ve aft).

    Returns
    -------
    NDArray[np.complex128], shape (2, 2)
        Aerodynamic influence matrix.
    """
    a11 = (k_nd**2) / b_m - 2 * c_nd * 1j * k_nd / b_m
    a12 = -(k_nd**2) * a_h_nd - 1j * k_nd - 2 * c_nd * ((0.5 - a_h_nd) * 1j * k_nd + 1)
    a21 = 2 * c_nd * (0.5 + a_h_nd) * 1j * k_nd - (k_nd**2) * a_h_nd
    a22 = (
        2 * c_nd * (0.5 + a_h_nd) * ((0.5 - a_h_nd) * 1j * k_nd + 1) * b_m
        + (k_nd**2) * (a_h_nd**2) * b_m
        - (0.5 - a_h_nd) * 1j * k_nd * b_m
        + (k_nd**2) * b_m / 8
    )
    return np.array([[a11, a12], [a21, a22]], dtype=np.complex128)


def build_mass_matrix(
    mass_kg: float,
    inertia_alpha_kg_m2: float,
    b_m: float,
    x_alpha_nd: float,
) -> NDArray[np.float64]:
    """
    Build the 2D generalised structural mass matrix [heave, pitch].

    Parameters
    ----------
    mass_kg : float
        Section mass in kg.
    inertia_alpha_kg_m2 : float
        Pitch moment of inertia about elastic axis in kg·m².
    b_m : float
        Semi-chord length in metres.
    x_alpha_nd : float
        CG offset aft of elastic axis (dimensionless, +ve aft).

    Returns
    -------
    NDArray[np.float64], shape (2, 2)
        Generalised mass matrix.
    """
    m12 = mass_kg * b_m * x_alpha_nd
    return np.array(
        [[mass_kg, m12], [m12, inertia_alpha_kg_m2]],
        dtype=np.float64,
    )


def build_stiffness_matrix(
    k_heave_n_m: float,
    k_alpha_nm_rad: float,
) -> NDArray[np.float64]:
    """
    Build the 2D uncoupled structural stiffness matrix [heave, pitch].

    Parameters
    ----------
    k_heave_n_m : float
        Heave (bending) stiffness in N/m.
    k_alpha_nm_rad : float
        Pitch (torsional) stiffness in N·m/rad.

    Returns
    -------
    NDArray[np.float64], shape (2, 2)
        Stiffness matrix.
    """
    return np.array(
        [[k_heave_n_m, 0.0], [0.0, k_alpha_nm_rad]],
        dtype=np.float64,
    )


def compute_flutter(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    velocity_range_m_s: NDArray[np.float64],
    rho_kg_m3: float,
    b_m: float,
    a_h_nd: float,
    start_omega_rad_s: float,
    det_tol_nd: float = 0.001,
    p_init_damp_1_nd: float = 0.0,
    p_init_damp_2_nd: float = 0.01,
    pi_nd: float = 3.141592653589793,
) -> NDArray[np.complex128]:
    """
    Compute p-k flutter eigenvalues over a velocity sweep using the secant method.

    Parameters
    ----------
    mass_matrix_kg : NDArray[np.float64], shape (2, 2)
        Generalised structural mass matrix.
    stiffness_matrix_n_m : NDArray[np.float64], shape (2, 2)
        Structural stiffness matrix in N/m (or N·m/rad on pitch diagonal).
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep points in m/s. Must be monotonically increasing.
    rho_kg_m3 : float
        Air density in kg/m³.
    b_m : float
        Semi-chord length in metres.
    a_h_nd : float
        Elastic axis location aft of semi-chord mid-point (dimensionless).
    start_omega_rad_s : float
        Initial frequency estimate for this mode in rad/s.
    det_tol_nd : float
        Convergence tolerance on |det(F)|; iteration stops when below this value.
    p_init_damp_1_nd : float
        Damping scale for first secant seed: real(p_1) = p_init_damp_1_nd * k.
    p_init_damp_2_nd : float
        Damping scale for second secant seed: real(p_2) = p_init_damp_2_nd * k.
    pi_nd : float
        Mathematical constant π, loaded from defaults/constants.json.

    Returns
    -------
    NDArray[np.complex128], shape (nv,)
        Complex p values per speed: real(p) = damping (1/s), imag(p) = frequency (rad/s).

    Raises
    ------
    ValueError
        If velocity_range_m_s is empty or non-positive.
    """
    if len(velocity_range_m_s) == 0:
        raise ValueError("velocity_range_m_s must not be empty")
    if velocity_range_m_s[0] <= 0:
        raise ValueError(f"velocity_range_m_s must be positive; got {velocity_range_m_s[0]}")

    speed_inc_m_s = velocity_range_m_s[1] - velocity_range_m_s[0]

    k_2 = start_omega_rad_s * b_m / velocity_range_m_s[0]
    k_1 = k_2
    p_1 = complex(p_init_damp_1_nd * k_2, k_2)
    p_2 = complex(p_init_damp_2_nd * k_2, k_2)

    flutter_p = np.empty(len(velocity_range_m_s), dtype=np.complex128)

    for i, v_m_s in enumerate(velocity_range_m_s):
        converged = False

        while not converged:
            c_1 = theo_c(k_1)
            a1 = build_aero_matrix(k_1, c_1, b_m, a_h_nd)
            f_1 = (v_m_s / b_m) ** 2 * p_1**2 * mass_matrix_kg + stiffness_matrix_n_m - rho_kg_m3 * pi_nd * b_m * v_m_s**2 * a1
            det_f1 = np.linalg.det(f_1)

            c_2 = theo_c(k_2)
            a2 = build_aero_matrix(k_2, c_2, b_m, a_h_nd)
            f_2 = (v_m_s / b_m) ** 2 * p_2**2 * mass_matrix_kg + stiffness_matrix_n_m - rho_kg_m3 * pi_nd * b_m * v_m_s**2 * a2
            det_f2 = np.linalg.det(f_2)

            p_3 = (p_2 * det_f1 - p_1 * det_f2) / (det_f1 - det_f2)

            k_3 = np.imag(p_3)
            c_3 = theo_c(k_3)
            a3 = build_aero_matrix(k_3, c_3, b_m, a_h_nd)
            f_3 = (v_m_s / b_m) ** 2 * p_3**2 * mass_matrix_kg + stiffness_matrix_n_m - rho_kg_m3 * pi_nd * b_m * v_m_s**2 * a3
            det_f3 = np.linalg.det(f_3)

            if abs(det_f3) <= det_tol_nd:
                converged = True
            else:
                p_1, p_2 = p_2, p_3
                k_1, k_2 = np.imag(p_1), np.imag(p_2)

        flutter_p[i] = p_3
        p_1 = p_1 * v_m_s / (v_m_s + speed_inc_m_s)
        p_2 = p_2 * v_m_s / (v_m_s + speed_inc_m_s)

    return flutter_p
