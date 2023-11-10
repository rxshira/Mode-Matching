import math

import numpy as np
import matplotlib.pyplot as plt
from physunits import *
pi = np.pi


class BeamPathElement:
    def __init__(self, x, y, label, w0=0, z0=0, wavelength=None):
        self.w0 = w0
        self.z0 = z0
        self.x = x
        self.y = y
        self.label = label
        self.wavelength = wavelength
        self.q0 = None if wavelength is None else self.q_at_min_waist()

    def set_element_properties(self, w0, z0):
        self.w0 = w0
        self.z0 = z0
        self.q0 = self.q_at_min_waist()

    def distance_to(self, other_device):
        dx = self.x - other_device.x
        dy = self.y - other_device.y
        return math.sqrt(dx*dx + dy*dy)

    def q_at_min_waist(self):
        # This is the Rayleigh length
        zR = pi * self.w0 ** 2 / self.wavelength
        # Get the beam parameter at a given z position
        q0 = self.z0 + 1j * zR  # !! SHOULD BE -, not +
        return q0

    def waist_at_z(self, z):
        q = self.q0 - z
        wz = self.w0 * np.sqrt(1 + (q.real / q.imag) ** 2)
        return np.array(wz), self.outout_wavelength()

    def outout_wavelength(self):
        return self.wavelength

    def __str__(self):
        return self.label

class OpticalDevice(BeamPathElement):
    def __init__(self, abcd_matrix, x, y, label):
        super().__init__(x, y, label)
        self.abcd_matrix = abcd_matrix

    def beam_q_param(self, qin):
        A = self.abcd_matrix.item(0, 0)
        B = self.abcd_matrix.item(0, 1)
        C = self.abcd_matrix.item(1, 0)
        D = self.abcd_matrix.item(1, 1)
        return (A + (B / qin)) / (C + (D / qin))

    def prepare_device_params(self, q, location_on_path, wavelength):
        self.wavelength = wavelength
        q_at_device_location = q-location_on_path
        q_new_beam_segment = self.beam_q_param(q_at_device_location)  # beam params for the next segment
        z0 = q_new_beam_segment.real + location_on_path  # z0 value as it should be based on the new beam params
        w0 = np.sqrt(self.wavelength * q_new_beam_segment.imag / np.pi)
        self.set_element_properties(w0, z0)

# IMPORTANT - everything is in cm

NPRO_wavelength = 1064.50*nm
NPRO_wavelength_green = NPRO_wavelength/2

def thin_lens_abcd_matrix(f):
    return np.matrix([[1.0, 0.0], [1.0/-f, 1.0]])

def refraction_curved_abcd_matrix(R, n1, n2):
    return np.matrix([[1, 0], [(n1-n2)/(R*n2), n1/n2]])

def ref_prog_ref_abcd_matrix(n1, n2, d):
    in_ref = np.matrix([[1,0], [0, n1/n2]])
    prop = np.matrix([[1, d], [0, 1]])
    out_ref = np.matrix([[1, 0], [0, n2/n1]])
    return np.matmul(np.matmul(out_ref, prop), in_ref)

# plot of the beam waist as a function of z_arr
def shade_gaussian(axis, z, waist_profile, color):
    axis.plot(z, waist_profile, c=color)
    axis.plot(z, -waist_profile, c=color)
    axis.fill_between(z, waist_profile/um, -waist_profile/um, color=color, alpha=0.2)

devices = [BeamPathElement(0, 0, "Start", w0=217*um, z0=5*cm, wavelength=NPRO_wavelength),
           OpticalDevice(thin_lens_abcd_matrix(-10*cm), 7*cm, 0, "L1"),
           OpticalDevice(thin_lens_abcd_matrix(-20*cm), 37*cm, 5*cm, "L2"),
           OpticalDevice(thin_lens_abcd_matrix(-40*cm), 67*cm, 0, "L3"),
           BeamPathElement(117 * cm, 0, "End")]

propagation_segments_z = []

distance_so_far = 0
for index, device in enumerate(devices[:-1]):
    next_device = devices[index+1]
    distance_to_next = device.distance_to(next_device)
    next_device_z = distance_so_far+distance_to_next
    segment_range = np.arange(start=distance_so_far, stop=next_device_z, step=1*mm)
    distance_so_far = next_device_z
    propagation_segments_z.append(segment_range)

beam_waist = np.array([])
all_z = np.array([])
for i, segment in enumerate(propagation_segments_z):
    current_device = devices[i]
    curr_propagation, output_wavelength = current_device.waist_at_z(segment)
    beam_waist = np.concatenate((beam_waist, curr_propagation))
    all_z = np.concatenate((all_z, segment))
    if i+1 != len(propagation_segments_z):
        devices[i+1].prepare_device_params(current_device.q0, segment[-1], output_wavelength)

ax = plt.subplot(111)
shade_gaussian(ax, all_z, beam_waist,color='red')
ax.set_title('Beam waist plot after optical path')
ax.set_xlabel('z position (m)')
ax.set_ylabel('Beam Radius (um)')
devices_locs = [seg[-1] for seg in propagation_segments_z]
devices_locs.pop()
device_labels = [d.label for d in devices[1:-1]]
pos = [int(loc/mm) for loc in devices_locs]
heights = np.array([p/um for p in beam_waist[pos]])
plt.vlines(x=devices_locs, ymax=heights, ymin=-heights,  colors='k')
for i, label in enumerate(device_labels):
    plt.text(s=" "+label, x=devices_locs[i], y=10, c="blue")
plt.show()
