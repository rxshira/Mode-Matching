import numpy as np
import matplotlib.pyplot as plt
from physunits import *
pi = np.pi

# IMPORTANT - everything is in cm

NPRO_wavelength = 1064.50*nm
NPRO_wavelength_green = NPRO_wavelength/2

W0 = 217*um
Z0 = 5*cm


def thin_lens_abcd_matrix(f):
    return np.matrix([[1.0, 0.0], [1.0/-f, 1.0]])


def refraction_curved_abcd_matrix(R, n1, n2):
    return np.matrix([[1, 0], [(n1-n2)/(R*n2), n1/n2]])


def beam_q_param(qin, abcd_matrix):
    A = abcd_matrix.item(0, 0)
    B = abcd_matrix.item(0, 1)
    C = abcd_matrix.item(1, 0)
    D = abcd_matrix.item(1, 1)
    return (A + (B / qin)) / (C + (D / qin))


def q_at_min_waist(z0, w0):
    # This is the Rayleigh length
    zR = pi * w0**2 / NPRO_wavelength
    # Get the beam parameter at a given z position
    q0 = z0 + 1j*zR  # !! SHOULD BE -, not +
    return q0


def waist_at_z(z, z0, w0):
    q0 = q_at_min_waist(z0, w0)
    q = q0-z
    wz = w0 * np.sqrt(1 + (q.real/q.imag)**2)
    return np.array(wz), q0


# plot of the beam waist as a function of z_arr
def shade_gaussian(axis, z, waist_profile, color):
    axis.plot(z, waist_profile, c=color)
    axis.plot(z, -waist_profile, c=color)
    axis.fill_between(z, waist_profile/um, -waist_profile/um, color=color, alpha=0.2)

# TEST the following path example:
#
# w(-z) = 2*cm------(f=100mm)------30*cm------(f=200mm)--------w(z=?)
#
# What happened?

ax = plt.subplot(111)

free_space_segments = [2*cm, 30*cm, 30*cm]
devices = [thin_lens_abcd_matrix(100*mm), thin_lens_abcd_matrix(200*mm)]
beam_path = [val for pair in zip(free_space_segments, devices) for val in pair]

propagation_segments_z = []
for indx, free_space_length in enumerate(free_space_segments):
    segment_range = np.arange(start=sum(free_space_segments[0:indx]),
                              stop=sum(free_space_segments[0:indx+1]),
                              step=1*mm)
    propagation_segments_z.append(segment_range)

beam_waist, q0_of_propagated_beam = waist_at_z(propagation_segments_z[0], z0=Z0, w0=W0)
all_z = propagation_segments_z[0]
for i in range(len(devices)):
    curr_device_abcd = devices[i]
    segment_z = propagation_segments_z[i+1]  # the segment after the device
    device_location = segment_z[0]
    q_at_device_location = q0_of_propagated_beam - device_location  # beam params at the device location
    q_new_beam_segment = beam_q_param(q_at_device_location, curr_device_abcd)  # beam params for the next segment
    z0 = q_new_beam_segment.real + device_location  # z0 value as it should be based on the new beam params
    w0 = np.sqrt(NPRO_wavelength * q_new_beam_segment.imag / np.pi)
    curr_propagation, q0_of_propagated_beam = waist_at_z(segment_z,
                                                         z0=z0,
                                                         w0=w0)
    beam_waist = np.concatenate((beam_waist, curr_propagation))  # accumulating the beam waist
    all_z = np.concatenate((all_z, segment_z))

shade_gaussian(ax, all_z, beam_waist, color='red')
ax.set_title('Beam waist plot after optical path')
ax.set_xlabel('z position (m)')
ax.set_ylabel('Beam Radius (um)')
devices_locs = [seg[-1] for seg in propagation_segments_z]
devices_locs.pop()
pos = [int(loc/mm) for loc in devices_locs]
heights = np.array([p/um for p in beam_waist[pos]])
plt.vlines(x=devices_locs, ymax=heights, ymin=-heights,  colors='k')
plt.show()

# q.real = z - z0
# q.imag = (0.5 * k * w0**2)

# wz = np.sqrt(q.imag / ( np.pi / NPRO_wavelength)) * np.sqrt(1 + (q.real/q.imag)**2)
