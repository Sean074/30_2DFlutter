# 2D Flutter Analysis
# Created: Sean O'Meara
# Date: June 11, 2024
# Version: Pre Basic and Stable

# Import
import flut
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## USER INPUT

# Geometry
geometry = {'b': 0.6,
            'a_h': -0.28,
            'x_alpha': 0.2}

# Mass Properties
mass = 10.4     # [kg]
i_alpha = 1.12
i_beta = 0.0607

# Stiffness
omega_alpha = 160 # [rad/sec]
omega_heave = 83.3 # [rad/sec]

modes = {'pitch': omega_alpha,
         'heave': omega_heave}

k_alpha = (omega_alpha**2)*i_alpha  # [N/rad]
k_heave = (omega_heave**2)*mass     # [N/m]

# Aerodynamics
density = 1.21      # Air density [kg/m^3]

## SOLUTION PARAMETERS
PI = 3.141
SPEED_MIN = 0.001   #[m/s]
SPEED_MAX = 200.0   #[m/s]
SPEED_INC = 1.0   #[m/s]

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


# The function version of the flutter solver
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
fig1, ax = plt.subplots(2,1)

# set the basic properties
ax[0].set_title('Flutter V-g')
ax[0].set_ylabel('Frequency [rad/sec]')
ax[0].set_xlabel('Velocity [TAS]')
ax[1].set_ylabel('Damping g')


# set the limits
ax[0].set_ylim(0,200)
ax[1].set_ylim(-0.1,0.1)
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

# Check if any of the modes cross the taget g level 0 and 0.03.

flutter_targets = []
growth_target = 0.0
output_targets = pd.DataFrame({})
speed_target = {}

for mode in modes.values():

    growth_rate = np.real(flutter_results[str(mode)+'_p']) / np.imag(flutter_results[str(mode)+'_p'])

    index = 0
    growth_0 = 0
    speed_0 = 0

    for speed_1 in flutter_results['TAS']:
        growth_1 = growth_rate[index]

        #print(f"{growth_0} {growth_1}")

        if growth_0 < growth_target and growth_1 > growth_target:
            print('found a crossing')
            print(f"Between {speed_0} and {speed_1}")

            # Linear interp to find crossing
            speed_target = speed_0 + (growth_target-growth_0)*(speed_1-speed_0)/(growth_1-growth_0)
            found_targets = pd.DataFrame({'Mode': mode,
                                         'Speed': speed_target,
                                         'Growth': growth_target,
                                        }, index=[0])
            
            output_targets = pd.concat([output_targets, found_targets], ignore_index=True)

            print(speed_target)

        growth_0 = growth_1
        speed_0 = speed_1
        
        index = index+1

print("Targets found:")
print(found_targets)

        # TODO print output to file
        # TODO make the plot a function and add to flut.py
        # TODO make the target and output functions and create flut_plot.py and flut_output.py      