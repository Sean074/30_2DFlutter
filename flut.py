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


def flut(start_mode_omega, velocity_vector, density, geometry, stiffness, mass):
    # Flutter solver p-k method with determinate iteration.  Note aero embedded, current limitation.
    # TODO make the A matrix a function of k and interpolate inside of flut.

    import flut
    import numpy as np

    # Stuff you need from problem
 
    b = geometry['b']

    mass_matrix = mass
    stiff_matrix = stiffness

    # Function Defined Constants
    DETERMINATE_ITERATIONS = 100
    detF_TOLL = 0.001

    PI = 3.1415926
    SPEED_INC = velocity_vector[1]-velocity_vector[0]

    # Initalize output
    flutter_p = np.array([])

    # FOR a given start point

    # Initial estimate of the differential operator p at low speed
    k_2 = start_mode_omega*b/velocity_vector[0]
    k_1 = k_2

    p_1 = -0.01*k_2 + k_2*1j
    p_2 = 0.0 + k_2*1j

    for speed in velocity_vector:
        # For each speed determine the flutter solution
            
        tollerance = False
        iterations = 0

        while (tollerance == False):
            iterations = iterations + 1

            # Iterate the det(F) until a converged solution
            c_1 = flut.theo_c(k_1)
            A_matrix_1 = flut.twod_aero(k=k_1, theo_c=c_1, geometry=geometry)       
            F_1 = (speed/b)**2*p_1**2*mass_matrix + stiff_matrix - density*PI*b*speed**2*A_matrix_1
            detF_1 = np.linalg.det(F_1)

            #evals, evecs = np.linalg.eig(F_V_p_n)

            c_2 = flut.theo_c(k_2)
            A_matrix_2 = flut.twod_aero(k=k_2, theo_c=c_2, geometry=geometry)
            F_2 = (speed/b)**2*p_2**2*mass_matrix + stiff_matrix - density*PI*b*speed**2*A_matrix_2
            detF_2 = np.linalg.det(F_2)
            
            p_3 = (p_2*detF_1 - p_1*detF_2)/(detF_1 - detF_2)

            # Doing this in this step and not using in in the next is a computational waste
            # TODO refactor to use this in the next iteration.
            k_3 = np.imag(p_3)
            c_3 = flut.theo_c(k_3)
            A_matrix_3 = flut.twod_aero(k=k_3, theo_c=c_3, geometry=geometry)       
            F_3 = (speed/b)**2*p_3**2*mass_matrix + stiff_matrix - density*PI*b*speed**2*A_matrix_3
            detF_3 = np.linalg.det(F_3)
            detF3_length = np.sqrt((np.real(detF_3)**2 + np.imag(detF_3)**2))

            if detF3_length <= detF_TOLL:
                #print(f"{iterations}: P_3 {p_3}, LF_3 = {detF3_length}")
                tollerance = True
            else:            
                # New values of p
                p_1 = p_2
                p_2 = p_3

                k_1 = np.imag(p_1)
                k_2 = np.imag(p_2)


        flutter_p = np.append(flutter_p, p_3)

        #print(f"Speed {speed}, omega {omega}, k_n {k_n}")

        p_1 = p_1*speed/(speed+SPEED_INC)
        p_2 = p_2*speed/(speed+SPEED_INC)

    return flutter_p # flutter_freq, flutter_damp

