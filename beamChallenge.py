import numpy as np
import matplotlib.pyplot as plt
pi = np.pi

from physunits import *

NPRO_wavelength = 1064*nm

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
print(gaussian_beam(z_arr, z0=0.0, w0=2*mm))

### CHALLENGE ### 

# Make a plot of the beam waist as a function of z_arr
# HINT: The beam waist can be computed by solving for it from the
# beam parameter.

###Answer###

def shade_gaussian(axis, z, waist_profile, color):
    axis.plot(z, waist_profile, c=color)
    axis.plot(z, -waist_profile, c=color)
    axis.fill_between(z, waist_profile, -waist_profile, color=color, alpha=0.5)

plt.figure()
ax = plt.subplot(111)
shade_gaussian(ax, z_arr, gaussian_beam(z_arr, z0=0.0, w0=2e-3)/mm, color='green')
ax.set_title('Beam waist plot')
ax.set_xlabel('z position (m)')
ax.set_ylabel('Beam Radius (mm)')

#plt.show()

### Propagate gaussian beams using ray transfer matrices

# Beam parameter is
# q = z0 + 1j * zR

# ABCD matrix transforms q as:
# [[A, B], 
#  [C, D]]

# This is how a beam parameter transforms under an ABCD matrix
# qout = (qin * A + B) / (qin * C + D)

# For a series of ABCD matrices, we can combine the net effect as
# Mpath = Mzlaser0 @ Mlens1 @ Mm1 @ Mlens2 ... @ METM

def free_space(qin, d):

    A = 1
    B = d 
    C = 0
    D = 1

    qout = (qin * A + B) / (qin * C + D)
    return qout

# HERE ADD thin_lens, refraction
# TEST the following path example:
# 
# w(-z) = 2*cm------(f=100mm)------30*cm------(f=200mm)--------w(z=?)
# 
# What happened?

# Test
w0 = 2*mm
k0 = 2 * np.pi * 1 / NPRO_wavelength
zR0 = 0.5 * k0 * w0**2
qinit = 0*m + 1j * zR0

# As a function of z_arr, propagate and plot! 
qout = []
for z_pos in z_arr:
    qout.append(free_space(qinit, z_pos))
qout = np.array(qout)

plt.figure()
ax = plt.subplot(111)
shade_gaussian(ax, z_arr, w0*np.sqrt(1 + (qout.real/qout.imag)**2), color='goldenrod')
plt.show()

# Thought provoking ideas: Challenge

# q = (z-z0) + 1j * (0.5 * k * w0**2)

# q.real = z - z0
# q.imag = (0.5 * k * w0**2)

# wz = np.sqrt(q.imag / ( np.pi / NPRO_wavelength)) * np.sqrt(1 + (q.real/q.imag)**2)
