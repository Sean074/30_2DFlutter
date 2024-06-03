
# Import
import flut
import numpy as np
import matplotlib.pyplot as plt

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

modes = {'pitch': 160,
         'heave': 83.3}

k_alpha = (omega_alpha**2)*i_alpha  # [N/rad]
k_heave = (omega_heave**2)*mass     # [N/m]

# Aerodynamics
density = 1.21      # Air density [kg/m^3]

## SOLUTION PARAMETERS
PI = 3.141
SPEED_MIN = 1.0   #[m/s]
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
    
    #flutter_results[str(mode)+'_freq'] = freq
    #flutter_results[str(mode)+'_damp'] = damp
    flutter_results[str(mode)+'_p'] = p_result


print(flutter_results.keys())

print(np.imag(flutter_results['160_p']))

fig1, ax = plt.subplots()
# ax.plot(flutter_results['TAS'], flutter_results['160_freq'])
ax.plot(flutter_results['TAS'], np.imag(flutter_results['160_p']) * flutter_results['TAS']/ geometry['b'])
plt.show()

