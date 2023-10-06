import numpy as np
import matplotlib.pyplot as plt
pi = np.pi

NPRO_wavelength = 1064e-9 #meters

def gaussian_beam(z, z0, w0):
    # Function shift
    z = z - z0

    # This is the Rayleigh length
    zR = pi * w0**2 / NPRO_wavelength

    # Waist as a function of propagation w(z)
    wz = w0 * np.sqrt(1 + (z/zR)**2)

    # Get the beam parameter at a given z position
    #qz = z - 1j*zR
    return wz


z_arr = np.arange(start=-100, stop=100, step=5)
print("this is beam waist in meters:")
print("when z=", z_arr, ", z0=0.0m, w0=2mm")
print(gaussian_beam(z_arr, z0=0.0, w0=2e-3))

### CHALLENGE ### 

# Make a plot of the beam waist as a function of z_arr
# HINT: The beam waist can be computed by solving for it from the
# beam parameter.

###Answer###

plt.plot(z_arr, gaussian_beam(z_arr, z0=0.0, w0=2e-3), 'ro')
plt.suptitle('Beam waist plot')
plt.xlabel('z position (m)')
plt.ylabel('Beam Radius (m)')
plt.show()