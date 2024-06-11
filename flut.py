# Flutter Analysis Functions

# theo_c: Theoderson Function approximation

# twod_aero: 
# twod_mass
# twod_stiffness

# flut


import numpy as np

def theo_c (k):
# Theoderson Function approximation
    c_of_k = 1 - 0.165/(1-0.041/k*1j) - 0.335/(1-0.32/k*1j)

    return c_of_k


def twod_aero (k, theo_c, geometry):

    b = geometry['b']
    a_h = geometry['a_h']

# 2D Aerodynamics assumed the oscillatory aerodynamics so position and velocity
# are related.
    a11 = (k**2)/b - 2*theo_c*1j*k/b
    a12 = -(k**2)*a_h - 1j*k - 2*theo_c*((0.5-a_h)*1j*k+1)
    a21 = 2*theo_c*(0.5+a_h)*1j*k - (k**2)*a_h
    a22 = 2*theo_c*(0.5+a_h)*((0.5-a_h)*1j*k+1)*b + (k**2)*(a_h**2)*b - (0.5-a_h)*1j*k*b + (k**2)*b/8

    a_matrix = np.matrix([[a11, a12], [a21, a22]])

    return a_matrix


def twod_mass (mass, inertia_alpha, geometry):

    b = geometry['b']
    x_alpha = geometry['x_alpha']

# 2D Mass matrix
    m11 = mass
    m12 = mass*b*x_alpha
    m21 = mass*b*x_alpha
    m22 = inertia_alpha

    m_matrix = np.matrix([[m11, m12], [m21, m22]])

    return m_matrix


def twod_stiffness (k_alpha, k_heave):
# 2D Stiffness matrix
    k11 = k_heave
    k12 = 0
    k21 = 0
    k22 = k_alpha

    k_matrix = np.matrix([[k11, k12], [k21, k22]])

    return k_matrix


def flut(start_mode_omega,velocity_vector, density, geometry, stiffness, mass):
    # Flutter solver p-k method with determinate iteration.  Note aero embedded, currentlimitation.
    # TODO make the A matrix a function of of k and interpolate inside of flut.

    import flut
    import numpy as np

    # Stuff you need from problem
 
    b = geometry['b']

    mass_matrix = mass
    stiff_matrix = stiffness

    # Function Defined Constants
    DETERMINATE_ITERATIONS = 10
    K_TOLLERANCE = 0.0001

    PI = 3.141
    SPEED_INC = velocity_vector[1]-velocity_vector[0]

    # Initalize output
    flutter_freq = np.array([])
    flutter_damp = np.array([])
    flutter_p = np.array([])

    # FOR a given start point

    # Initial estimate of the differential operator p at low speed
    k_n = start_mode_omega*b/velocity_vector[0]
    k_np1 = start_mode_omega*b/velocity_vector[0]

    p_n = 0.01*k_n + k_n*1j
    p_np1 = 0 + k_np1*1j

    for speed in velocity_vector:
        # For each speed determine the flutter solution
            
        tollerance=False
        iterations = 1

        while (tollerance == False) & (iterations <= DETERMINATE_ITERATIONS):
            # Iterate the det(F) until a converged solution
            
            c_kn = flut.theo_c(k_n)
            A_matrix_n = flut.twod_aero(k=k_n, theo_c=c_kn, geometry=geometry)       
            F_V_p_n = (speed/b)**2*p_n**2*mass_matrix + stiff_matrix - density*PI*b*speed**2*A_matrix_n
            F_n = np.linalg.det(F_V_p_n)

            c_knp1 = flut.theo_c(k_np1)
            A_matrix_np1 = flut.twod_aero(k=k_np1, theo_c=c_knp1, geometry=geometry)
            F_V_p_np1 = (speed/b)**2*p_np1**2*mass_matrix + stiff_matrix - density*PI*b*speed**2*A_matrix_np1
            F_np1 = np.linalg.det(F_V_p_np1)

            p_np2 = (p_np1*F_n - p_n*F_np1)/(F_n - F_np1)

            if (np.imag(p_np2) - np.imag(p_np1))**2 < K_TOLLERANCE**2:
                tollerance = True
                    
            iterations = iterations + 1

            # New values of k
            k_n = np.imag(p_np1)
            k_np1 = np.imag(p_np2)

            # Freq and Damping
            omega = speed * np.imag(p_np2) / b
            damping_g = np.real(p_np2) / np.imag(p_np2)

            # New values of p
            p_n = p_np1
            p_np1 = p_np2

        if iterations == DETERMINATE_ITERATIONS:
            print(f"WARNING: Max iterations reached at speed {speed} delta k = {k_np1-k_n}")

        flutter_freq = np.append(flutter_freq,omega)
        flutter_damp = np.append(flutter_damp,damping_g)
        flutter_p = np.append(flutter_p,p_np1)

        p_n = p_n*speed/(speed+SPEED_INC)
        p_np1 = p_np1*speed/(speed+SPEED_INC)

    return flutter_p # flutter_freq, flutter_damp
