import numpy as np
import matplotlib.pyplot as plt
import functools
from physunits import *
pi = np.pi

NPRO_wavelength_RED= 1064.50*nm
NPRO_wavelength_GREEN= 532.25*nm


def beam_param_at_z(z, z0, w0):
    z = z - z0
    # This is the Rayleigh length
    zR = pi * w0**2 / NPRO_wavelength_RED
    # Waist as a function of propagation w(z)
    # Get the beam parameter at a given z position
    qz = z - 1j*zR
    return qz


def thin_lens_abcd_matrix(f):
    return np.matrix([[1, 0], [-1/f, 1]])


def free_space_abcd_matrix(d):
    return np.matrix([[1, d], [0, 1]])


def refraction_curved_abcd_matrix(R, n1, n2):
    return np.matrix([[1, 0], [(n1-n2)/(R*n2), n1/n2]])


def beam_q_param(qin, abcd_matrix):
    A = abcd_matrix.item(0,0)
    B = abcd_matrix.item(0, 1)
    C = abcd_matrix.item(1, 0)
    D = abcd_matrix.item(1, 1)
    return (A + (B / qin)) / (C + (D / qin))

def refraction_flat_mirror():
    return np.matrix([[1, 0], [0, 1]])



def combine_beam_path_components(abcd_arr):
    return functools.reduce(np.dot, reversed(abcd_arr))


def propagate(z_arr, abcd_matrix, z0, w0):
    waist = []
    for z_pos in z_arr:
        qz = beam_q_param(beam_param_at_z(z_pos, z0, w0), abcd_matrix)
        wz = w0 * np.sqrt(1 + (qz.real / qz.imag) ** 2)
        waist.append(wz)
    return np.array(waist)


# plot of the beam waist as a function of z_arr
def shade_gaussian(axis, z, waist_profile, color):
    axis.plot(z, waist_profile, c=color)
    axis.plot(z, -waist_profile, c=color)
    axis.fill_between(z, waist_profile, -waist_profile, color=color, alpha=0.5)

# TEST the following path example:
#
# w(-z) = 2*cm------(f=100mm)------30*cm------(f=200mm)--------w(z=?)
#
# What happened?


z_arr = np.arange(start=-100, stop=100, step=5)
ax = plt.subplot(111)
beam_path = [free_space_abcd_matrix(2*cm),
             thin_lens_abcd_matrix(100*mm),
             free_space_abcd_matrix(30*cm),
             thin_lens_abcd_matrix(200*mm)]
path_abcd = combine_beam_path_components(beam_path)
shade_gaussian(ax, z_arr, propagate(z_arr, path_abcd, z0=0.0, w0=2*mm) / mm, color='blue')

ax.set_title('Beam waist plot after optical path')
ax.set_xlabel('z position (m)')
ax.set_ylabel('Beam Radius (mm)')
plt.show()

# LIGO Laser Beam Path #





# Thought-provoking ideas: Challenge

# q = (z-z0) + 1j * (0.5 * k * w0**2)

# q.real = z - z0
# q.imag = (0.5 * k * w0**2)

# wz = np.sqrt(q.imag / ( np.pi / NPRO_wavelength)) * np.sqrt(1 + (q.real/q.imag)**2)

