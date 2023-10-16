import numpy as np
import matplotlib.pyplot as plt
import functools
import math
from physunits import *
pi = np.pi

# IMPORTANT - everything is in cm

NPRO_wavelength = 1064.50*nm / cm
NPRO_wavelength_green = NPRO_wavelength/2

W0 = 217*um / cm
Z0 = 5

def q_at_z(z, z0, w0):
    z = z - z0
    # This is the Rayleigh length
    zR = pi * w0**2 / NPRO_wavelength
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


def combine_abcd_beam_path_components(abcd_arr):
    return functools.reduce(np.dot, reversed(abcd_arr))


def initial_waist_at_z(z, z0, w0):
    # Function shift
    z = z - z0
    # This is the Rayleigh length
    zR = pi * w0**2 / NPRO_wavelength
    # Waist as a function of propagation w(z)
    wz = w0 * np.sqrt(1 + (z/zR)**2)
    return np.array(wz)


def propagate(z_values, abcd_matrix, z0, w0):
    waist = []
    for z_pos in z_values:
        qz = beam_q_param(q_at_z(z_pos, z0, w0), abcd_matrix)
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


#z_arr = np.arange(start=-100, stop=100, step=5)
ax = plt.subplot(111)

z_arrs = []
beam_path_partial_propagation_abcd = []

free_space_segments = [7, 30, 20]
beam_path = [free_space_abcd_matrix(free_space_segments[0]),
             thin_lens_abcd_matrix(10),
             free_space_abcd_matrix(free_space_segments[1]),
             thin_lens_abcd_matrix(20)]

# taking pairs of free space and device
for partial_path_end_index in range(2, len(beam_path)+1, 2):
    partial_bean_path = beam_path[0:partial_path_end_index]
    combined_partial_abcd = combine_abcd_beam_path_components(partial_bean_path)
    beam_path_partial_propagation_abcd.append(combined_partial_abcd)

accumulated_zs = []

# creating an array of the free spaces
for indx, free_space_length in enumerate(free_space_segments):
    z_arrs.append(np.arange(start=0, stop=free_space_length, step=0.1))
    accumulated_zs.append(np.arange(start=sum(free_space_segments[0:indx]),
                                    stop=sum(free_space_segments[0:indx+1]),
                                    step=0.1))

# calculating the propagation throughout the entire beam:
# calculate the initial propagation in free space
# breaking the beam to segments of free space + optical device, and calculate the propagation in the next free space

beam_waist = initial_waist_at_z(z_arrs[0], z0=Z0, w0=W0)

for i in range(len(beam_path_partial_propagation_abcd)):
    curr_path = beam_path_partial_propagation_abcd[i]  # path until (incl.) an optical device
    # question - does Z always start from 0, or from the z we got up till this point?
    curr_z_arr = z_arrs[i+1]  # z into the beam propagates (after the device)
    curr_accumulated_z = accumulated_zs[i+1]
    curr_propagation = propagate(curr_accumulated_z, curr_path, Z0, W0)  # propagation
    beam_waist = np.concatenate((beam_waist, curr_propagation))  # accumulating the beam waist

all_z = np.arange(start=0, stop=sum(free_space_segments), step=0.1)
path_abcd = combine_abcd_beam_path_components(beam_path)

shade_gaussian(ax, all_z, beam_waist, color='red')

ax.set_title('Beam waist plot after optical path')
ax.set_xlabel('z position (m)')
ax.set_ylabel('Beam Radius (mm)')
plt.show()


# Thought-provoking ideas: Challenge

# q = (z-z0) + 1j * (0.5 * k * w0**2)

# q.real = z - z0
# q.imag = (0.5 * k * w0**2)

# wz = np.sqrt(q.imag / ( np.pi / NPRO_wavelength)) * np.sqrt(1 + (q.real/q.imag)**2)
