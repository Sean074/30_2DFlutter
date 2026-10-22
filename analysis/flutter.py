"""p-k flutter analysis: matrix builders and eigenvalue-tracking p-k solver."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# R. T. Jones two-pole fit of C(k) (NACA Report 681): C = 1 - sum(A_i k / (k - i B_i)).
_JONES_A_ND = (0.165, 0.335)
_JONES_B_ND = (0.0455, 0.300)


def theo_c(k_nd: float) -> complex:
    """
    Theodorsen aerodynamic transfer function (R. T. Jones two-pole approximation).

    C(k) ~ 1 - 0.165 / (1 - 0.0455 i / k) - 0.335 / (1 - 0.300 i / k); largest error against
    the exact function is |dC| = 0.0145 for k <= 2.

    Parameters
    ----------
    k_nd : float
        Reduced frequency (dimensionless), k = b * omega / V.

    Returns
    -------
    complex
        Theodorsen function C(k).
    """
    return 1 - sum(a / (1 - b / k_nd * 1j) for a, b in zip(_JONES_A_ND, _JONES_B_ND))


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


def compute_structural_frequencies(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Compute the coupled in-vacuo natural frequencies of the section.

    Parameters
    ----------
    mass_matrix_kg : NDArray[np.float64], shape (2, 2)
        Generalised structural mass matrix.
    stiffness_matrix_n_m : NDArray[np.float64], shape (2, 2)
        Structural stiffness matrix.

    Returns
    -------
    NDArray[np.float64], shape (2,)
        Natural frequencies in rad/s, ascending.

    Raises
    ------
    ValueError
        If the mass or stiffness matrix is not positive definite.
    """
    for name, mat in (("mass", mass_matrix_kg), ("stiffness", stiffness_matrix_n_m)):
        if np.any(np.linalg.eigvalsh(mat) <= 0.0):
            raise ValueError(
                f"{name} matrix is not positive definite; check mass, inertia and CG offset"
            )
    omega_sq_arr = np.linalg.eigvals(np.linalg.solve(mass_matrix_kg, stiffness_matrix_n_m))
    return np.sort(np.sqrt(np.real(omega_sq_arr)))


def _solve_mode_at_speed(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    v_m_s: float,
    rho_kg_m3: float,
    b_m: float,
    a_h_nd: float,
    p_guess_nd: complex,
    k_tol_nd: float,
    max_iter_nd: int,
    pi_nd: float,
) -> complex:
    """Iterate k = imag(p) to convergence, following the eigenvalue nearest p_guess_nd."""
    p_nd = p_guess_nd
    for _ in range(max_iter_nd):
        k_nd = np.imag(p_nd)
        if k_nd <= 0.0:
            raise ValueError(
                f"reduced frequency fell to zero at V = {v_m_s:.3f} m/s "
                "(static divergence); reduce the maximum velocity"
            )
        aero_matrix = build_aero_matrix(k_nd, theo_c(k_nd), b_m, a_h_nd)
        q_matrix = stiffness_matrix_n_m - rho_kg_m3 * pi_nd * b_m * v_m_s**2 * aero_matrix
        # (V/b)^2 p^2 M + Q = 0  ->  p^2 = -(b/V)^2 eig(M^-1 Q); keep the oscillatory root.
        p_sq_arr = -((b_m / v_m_s) ** 2) * np.linalg.eigvals(
            np.linalg.solve(mass_matrix_kg, q_matrix)
        )
        root_arr = np.sqrt(p_sq_arr.astype(np.complex128))
        root_arr = np.where(np.imag(root_arr) < 0.0, -root_arr, root_arr)
        p_new_nd = complex(root_arr[np.argmin(np.abs(root_arr - p_nd))])
        if abs(np.imag(p_new_nd) - k_nd) <= k_tol_nd * abs(np.imag(p_new_nd)):
            return p_new_nd
        p_nd = p_new_nd
    raise ValueError(f"p-k iteration did not converge at V = {v_m_s:.3f} m/s")


def compute_flutter(
    mass_matrix_kg: NDArray[np.float64],
    stiffness_matrix_n_m: NDArray[np.float64],
    velocity_range_m_s: NDArray[np.float64],
    rho_kg_m3: float,
    b_m: float,
    a_h_nd: float,
    k_tol_nd: float = 1e-10,
    max_iter_nd: int = 200,
    pi_nd: float = 3.141592653589793,
) -> NDArray[np.complex128]:
    """
    Compute p-k flutter eigenvalues for every structural mode over a velocity sweep.

    At each speed and for each mode, the eigenvalues of F(p) = (V/b)^2 p^2 M + K - rho pi b V^2 A(k)
    are found for fixed k, the root nearest the mode's current estimate is kept, and k = imag(p)
    is updated until it converges. Modes are seeded from the coupled structural frequencies and
    each speed is seeded by linear extrapolation from the previous two, which keeps mode identity.

    Parameters
    ----------
    mass_matrix_kg : NDArray[np.float64], shape (2, 2)
        Generalised structural mass matrix.
    stiffness_matrix_n_m : NDArray[np.float64], shape (2, 2)
        Structural stiffness matrix in N/m (or N·m/rad on pitch diagonal).
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep points in m/s. Must be positive and strictly increasing.
    rho_kg_m3 : float
        Air density in kg/m³.
    b_m : float
        Semi-chord length in metres.
    a_h_nd : float
        Elastic axis location aft of semi-chord mid-point (dimensionless).
    k_tol_nd : float
        Relative convergence tolerance on the reduced frequency k = imag(p).
    max_iter_nd : int
        Maximum p-k iterations per mode per speed.
    pi_nd : float
        Mathematical constant π, loaded from defaults/constants.json.

    Returns
    -------
    NDArray[np.complex128], shape (n_modes, nv)
        Complex p per mode and speed, modes ordered by ascending structural frequency
        (see `compute_structural_frequencies`): real(p) = decay rate, imag(p) = reduced
        frequency k, so the physical frequency is imag(p) * V / b in rad/s.

    Raises
    ------
    ValueError
        If the velocity sweep is empty, non-positive or not increasing, the matrices are not
        positive definite, the iteration fails to converge, or two modes collapse onto one root.
    """
    if len(velocity_range_m_s) == 0:
        raise ValueError("velocity_range_m_s must not be empty")
    if velocity_range_m_s[0] <= 0:
        raise ValueError(f"velocity_range_m_s must be positive; got {velocity_range_m_s[0]}")
    if np.any(np.diff(velocity_range_m_s) <= 0):
        raise ValueError("velocity_range_m_s must be strictly increasing")

    omega_arr_rad_s = compute_structural_frequencies(mass_matrix_kg, stiffness_matrix_n_m)
    n_modes = len(omega_arr_rad_s)
    # Track s = p V / b (rad/s, complex): it varies smoothly with V, unlike p itself.
    s_arr_rad_s = np.empty((n_modes, len(velocity_range_m_s)), dtype=np.complex128)

    for i, v_m_s in enumerate(velocity_range_m_s):
        for mode in range(n_modes):
            if i == 0:
                s_guess_rad_s = 1j * omega_arr_rad_s[mode]
            elif i == 1:
                s_guess_rad_s = s_arr_rad_s[mode, 0]
            else:
                dv_ratio = (v_m_s - velocity_range_m_s[i - 1]) / (
                    velocity_range_m_s[i - 1] - velocity_range_m_s[i - 2]
                )
                s_guess_rad_s = s_arr_rad_s[mode, i - 1] + dv_ratio * (
                    s_arr_rad_s[mode, i - 1] - s_arr_rad_s[mode, i - 2]
                )
            p_nd = _solve_mode_at_speed(
                mass_matrix_kg,
                stiffness_matrix_n_m,
                v_m_s,
                rho_kg_m3,
                b_m,
                a_h_nd,
                s_guess_rad_s * b_m / v_m_s,
                k_tol_nd,
                max_iter_nd,
                pi_nd,
            )
            s_arr_rad_s[mode, i] = p_nd * v_m_s / b_m

        s_col = s_arr_rad_s[:, i]
        gap_arr_rad_s = np.abs(s_col[:, None] - s_col[None, :])
        np.fill_diagonal(gap_arr_rad_s, np.inf)
        if np.min(gap_arr_rad_s) <= 1e-6 * np.max(np.abs(s_col)):
            raise ValueError(
                f"modes collapsed onto one root at V = {v_m_s:.3f} m/s; "
                "increase the number of velocity points"
            )

    return s_arr_rad_s * b_m / velocity_range_m_s[None, :]


def damping_g_nd(p_arr_nd: NDArray[np.complex128]) -> NDArray[np.float64]:
    """
    Equivalent structural damping g = 2 * gamma, with gamma = real(p) / imag(p).

    gamma is the p-k decay-rate ratio (sigma / omega); g = 2 * gamma is the flutter convention,
    so g = 0.03 corresponds to gamma = 1.5 %. Both cross zero at the same speed.

    Parameters
    ----------
    p_arr_nd : NDArray[np.complex128]
        Complex p values, any shape.

    Returns
    -------
    NDArray[np.float64]
        Damping g, same shape as p_arr_nd.
    """
    return 2.0 * np.real(p_arr_nd) / np.imag(p_arr_nd)


def find_flutter_speed(
    p_arr_nd: NDArray[np.complex128],
    velocity_range_m_s: NDArray[np.float64],
    b_m: float,
) -> tuple[float, float] | None:
    """
    Find the first flutter onset of one mode by interpolating the g = 0 crossing.

    Damping g = 2 real(p) / imag(p) (`damping_g_nd`). Onset is the first interval where g goes from
    negative to non-negative; speed and frequency are interpolated linearly in g within it,
    so the result does not snap to the velocity grid.

    Parameters
    ----------
    p_arr_nd : NDArray[np.complex128], shape (nv,)
        Complex p for one mode, as one row of `compute_flutter` output.
    velocity_range_m_s : NDArray[np.float64], shape (nv,)
        Velocity sweep points in m/s.
    b_m : float
        Semi-chord length in metres.

    Returns
    -------
    tuple[float, float] or None
        (flutter speed in m/s, flutter frequency in rad/s), or None if g never crosses zero.
    """
    g_arr_nd = damping_g_nd(p_arr_nd)
    freq_arr_rad_s = np.imag(p_arr_nd) * velocity_range_m_s / b_m
    crossing_arr = np.flatnonzero((g_arr_nd[:-1] < 0.0) & (g_arr_nd[1:] >= 0.0))
    if len(crossing_arr) == 0:
        return None
    j = crossing_arr[0]
    frac_nd = -g_arr_nd[j] / (g_arr_nd[j + 1] - g_arr_nd[j])
    v_flutter_m_s = velocity_range_m_s[j] + frac_nd * (
        velocity_range_m_s[j + 1] - velocity_range_m_s[j]
    )
    omega_flutter_rad_s = freq_arr_rad_s[j] + frac_nd * (freq_arr_rad_s[j + 1] - freq_arr_rad_s[j])
    return float(v_flutter_m_s), float(omega_flutter_rad_s)
