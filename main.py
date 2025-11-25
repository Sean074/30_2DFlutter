# 2D Flutter Analysis
# Created: Sean O'Meara
# Date: Nov 24, 2025
# Version: Pre Basic and Stable

# Import
import flut
import numpy as np
import matplotlib.pyplot as plt

## USER INPUT
# TODO make this a septate function/file
# TODO make user commd line input and output to file

# Geometry
geometry = {'b': 0.5,       # Semi chord length [m]
            'a_h': 0.25,   # Location aft of the semi-chord of the elastic axis [% semi chord, +ve aft]
            'x_alpha': 0.2, # Location of the cg aft of the elastic axis [% semi chord, +ve aft]
            }

# Mass Properties
mass = 10.4     # [kg]
i_alpha = 1.12
i_beta = 0.0607

# Stiffness
omega_alpha = 30 # [rad/sec]
omega_heave = 20 # [rad/sec]

modes = {'pitch': omega_alpha,
         'heave': omega_heave}

k_alpha = (omega_alpha**2)*i_alpha  # [N/rad]
k_heave = (omega_heave**2)*mass     # [N/m]

# Aerodynamics
density = 1.21      # Air density [kg/m^3]

## SOLUTION PARAMETER
PI = 3.1415926
SPEED_MIN = .01   #[m/s]
SPEED_MAX = 50.0   #[m/s]
SPEED_INC = .01   #[m/s]

# BUILD VELOCITY ARRAY
velocity = np.array([])
speed = SPEED_MIN

while speed <= (SPEED_MAX):
    velocity = np.append(velocity, speed)
    speed = speed + SPEED_INC

## BUILD MATRIX
# Build Mass Matrix
mass_matrix = flut.twod_mass(inertia_alpha=i_alpha, mass=mass, geometry=geometry)

# Build Stiffness Matrix
stiff_matrix = flut.twod_stiffness(k_alpha=k_alpha,k_heave=k_heave)


# The better function version.
flutter_results = {'TAS': velocity}

for mode in modes.values():
    print(f'Starting mode {mode}')

    p_result = flut.flut(start_mode_omega= mode,
                             velocity_vector=velocity,
                             density=density,
                             geometry=geometry,
                             stiffness=stiff_matrix,
                             mass=mass_matrix,
                             )
    
    flutter_results[str(mode)+'_p'] = p_result

# Make some plots
# TODO make this a septate function/file
# TODO add some data about the case
# TODO add some engineering controls
fig1, ax = plt.subplots(2,1)

# set the basic properties
ax[0].set_title('Flutter V-g')
ax[0].set_ylabel('Frequency [rad/sec]')
ax[0].set_xlabel('Velocity [TAS]')
ax[1].set_ylabel('Damping g')

# set the limits
ax[0].set_ylim(0,50)
ax[1].set_ylim(-0.2,0.3)
ax[0].set_xlim(0,SPEED_MAX)
ax[1].set_xlim(0,SPEED_MAX)


# set the grid on
ax[0].grid('on')
ax[1].grid('on')


# Plot all modes in the modes array
for mode in modes.values():
    ax[0].plot(flutter_results['TAS'], np.imag(flutter_results[str(mode)+'_p']) * flutter_results['TAS']/ geometry['b'])
    ax[1].plot(flutter_results['TAS'], np.real(flutter_results[str(mode)+'_p']) / np.imag(flutter_results[str(mode)+'_p']))

# Add the g=0.03 line
ax[1].plot(flutter_results['TAS'], flutter_results['TAS']*0+0.03,linestyle='dotted', color='black', linewidth='2')
ax[1].plot(flutter_results['TAS'], flutter_results['TAS']*0+0.0, color='black', linewidth='2')

plt.show()

# TODO Find a target value
# search P-target to find index where damping went positive (changed sign?).
# 0.0 then occurs between index and index-1.
# linear interpolate to find value of speed.

