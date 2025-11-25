# H.1 Dynamics Aeroelastic Calculation

# NOT WORKING!

# Imports
import numpy as np
from numpy import linalg as npla

# CONSTANTS
PI = 3.141

# System Definition Geometry and Structure
s = 7.5            # semi span [m]
c = 2.0            # chord [m]
m = 100.0          # unit mass/area of wing
kappa_freq = 5     # heave frequecy [Hz]
thete_freq = 10    # torsion freq [Hz]
xcm = 0.5 * c      # posiion of cg from leading edge
xf = 0.48 * c      # position of the elastic axis from leading edge

e = xf/c-0.25      # eccentricity between the flexual axis and the aero center at (0.25c)
mass_le = (m*c**2 - 2*m*c*xcm) / (2*xcm) # leading edge mass term

damping = True

if damping:  # structural proportional damping C = alpha*M+beta*K
    z1 = 0.0   # critical damping 1st freq
    z2 = 0.0   # critical damping 2nd freq
    w1 = 2*2*PI    # first freq
    w2 = 14*2*PI   # 2nd freq
    alpha = 2*w1*w2*(-z2*w1 + z1*w2)/(w1**2*w2**2)
    beta = 2*(z2*w2-z1*w1) / (w2**2-w1**2)

# Aerodynamic Definition
vel_start = 1.0      # initial starting velocity
vel_max   = 18.0    # velocity range
vel_inc   = 1.0    # velicity increment

rho = 1.225        # air density

a = 2*PI      # lift curve slope
m_thete_dot = -1.2 # Unsteady aero damping

# Set up system matrices

# Mass Matrix A
a11 = (m*s**3*c)/3 + mass_le*s**3/3  # I_kappa
a22 = m*s*(c**3/3 - c*c*xf + xf**2*c) + mass_le*xf**2*s # I theta
a12 = m*s**2/2*(c**2/2 - c*xf) - mass_le*xf*s**2/2
a21 = a12

A = np.array([[a11, a12],
              [a21, a22]])

# Structral stiffness matrix E
k11 = (kappa_freq*PI*2)**2*a11  # k_kappa heave stiffness
k22 = (thete_freq*PI*2)**2*a11  # k_theta torsion stiffness

E = np.array([[k11, 0],
              [0, k22]])

# Aero and strucutral damping matrix C
def aero_damp_matrix(damping, velocity, density):

    if damping:
        c11 = c*s**3*a/6.0
        c12 = 0.0
        c21 = -c**2*s**2*e*a/4
        c22 = -c**3*s*m_thete_dot/8
        C = density*velocity*np.array([[c11, c12],
                                       [c21, c22]])
    else:
        C = np.array([[0, 0],
                      [0, 0]])

    return C

# Set up the eigen value problem 
# Loop for each velocity
for velocity in np.arange(vel_start, vel_max, vel_inc):
          
    C = aero_damp_matrix(damping = damping, velocity = velocity, density = rho)
   
    K = rho*velocity**2*np.array([[0, c*s**2*a/4], [0, -c**2*s*e*a/2]]) + E
   
    matrix11 = np.array([[0, 0], [0, 0]])
    matrix12 = np.array([[1, 0], [0, 1]])
    A_inv = npla.inv(A)
    matrix21 = npla.solve(A_inv,K)
    matrix22 = npla.solve(A_inv,C)

    matrix1 = np.concatenate((matrix11, matrix12),axis=1)
    matrix2 = np.concatenate((matrix21, matrix22),axis=1)
    matrix = np.concatenate((matrix1, matrix2),axis=0)

    landa = npla.eigvals(matrix)

    # Natural Freq and Damping Ratio
    for mode in landa:
        eigen_val_real = np.real(mode)
        eigen_val_imag = np.imag(mode)
        freq = np.sqrt(eigen_val_real**2+eigen_val_imag**2)/(2*PI)  # Hz
        damp = -100*eigen_val_real/freq

        print(velocity, freq, damp)