import math
from physunits import *
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
pi = np.pi


class LaserBeamSetup:
    def __init__(self, devices):
        self.devices = devices
        self.path_segments_z = self.calculate_segments()
        self.propagation_segments = []

    def calculate_segments(self):
        segments = []
        distance_so_far = 0
        for index, device in enumerate(self.devices[:-1]):
            next_device = self.devices[index + 1]
            distance_to_next = device.distance_to(next_device)
            next_device_z = distance_so_far + distance_to_next
            segment_range = np.arange(start=distance_so_far, stop=next_device_z, step=1 * mm)
            distance_so_far = next_device_z
            segments.append(segment_range)
        return segments

    def propagate(self):
        self.propagation_segments = []
        for i, segment in enumerate(self.path_segments_z):
            current_device = self.devices[i]
            curr_propagation, output_wavelength = current_device.waist_at_z(segment)
            self.propagation_segments.append(PropagationSegment(segment, curr_propagation, current_device.color))
            if i + 1 != len(self.path_segments_z):
                self.devices[i + 1].prepare_device_params(current_device.q0, segment[-1], output_wavelength)

    def plot_beam_path(self):
        self.propagate()
        axis = plt.subplot(111)
        for seg in self.propagation_segments:
            axis.plot(seg.z_values, seg.waist_values, c='black')
            axis.plot(seg.z_values, -seg.waist_values, c='black')
            axis.fill_between(seg.z_values, seg.waist_values / um, -seg.waist_values / um, color=seg.color, alpha=0.2)
        axis.set_title('Beam waist plot after optical path')
        axis.set_xlabel('z position (m)')
        axis.set_ylabel('Beam Radius (um)')
        devices_location = [seg[-1] for seg in self.path_segments_z[0:-1]]
        device_labels = [d.label for d in self.devices[1:-1]]
        heights = np.array([p.waist_values[-1] / um for p in self.propagation_segments[0:-1]])
        plt.vlines(x=devices_location, ymax=heights, ymin=-heights, colors='k')
        for device_loc, lbl in zip(devices_location, device_labels):
            plt.text(s=" " + lbl, x=device_loc, y=10, c="blue")
        plt.show()

    def waist_at_end(self):
        self.propagate()
        return self.propagation_segments[-1].waist_values[-1]


class BeamPathElement:
    def __init__(self, x, y, label, w0=0, z0=0, wavelength=None, color=""):
        self.w0 = w0
        self.z0 = z0
        self.x = x
        self.y = y
        self.label = label
        self.wavelength = wavelength
        self.color = color
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
        z_r = pi * self.w0 ** 2 / self.wavelength
        # Get the beam parameter at a given z position
        q0 = self.z0 + 1j * z_r  # !! SHOULD BE -, not +
        return q0

    def waist_at_z(self, z):
        q = self.q0 - z
        wz = self.w0 * np.sqrt(1 + (q.real / q.imag) ** 2)
        return np.array(wz), self.output_wavelength()

    def output_wavelength(self):
        return self.wavelength

    def __str__(self):
        return self.label


class OpticalDevice(BeamPathElement):
    def __init__(self, abcd_matrix, x, y, label, color):
        super().__init__(x, y, label, color=color)
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


class SHG(OpticalDevice):
    def __init__(self, n1, n2, n3, d1, d2, x, y, label, color):
        super().__init__(ABCD_Matrix_of.shg(n1, n2, n3, d1, d2), x, y, label, color=color)

    def prepare_device_params(self, q, location_on_path, wavelength):
        super().prepare_device_params(q, location_on_path, wavelength/2.0)


class Mirror(OpticalDevice):
    def __init__(self, x, y, label, color):
        super().__init__(ABCD_Matrix_of.mirror(), x, y, label, color)


class ThinLens(OpticalDevice):
    def __init__(self, f, x, y, label, color):
        super().__init__(ABCD_Matrix_of.thin_lens(f), x, y, label, color)


class FaradayIsolator(OpticalDevice):
    def __init__(self, n1, n2, d, x, y, label, color):
        super().__init__(ABCD_Matrix_of.faraday_isolator(n1, n2, d), x, y, label, color)


class FabryPerotCavity(OpticalDevice):
    def __init__(self, rc, d, x, y,  label, color):
        super().__init__(ABCD_Matrix_of.fabry_perot_cavity(rc, d), x, y, label, color)
        self.d = d
        self.rc = rc

    def stable_waist(self, wave_length):
        return ((wave_length/pi)**2 * ((self.rc-self.d)**2) * (self.rc+self.rc - self.d)) / (2*(self.rc-self.d))**2


class PropagationSegment:
    def __init__(self, z_values, waist_values, color):
        self.z_values = z_values
        self.waist_values = waist_values
        self.color = color


class ABCD_Matrix_of:
    @staticmethod
    def thin_lens(f):
        return np.matrix(np.array([[1.0, 0.0], [1.0 / f, 1.0]]))

    @staticmethod
    def refraction_curved(r, n1, n2):
        return np.matrix(np.array([[1, 0], [(n1 - n2) / (r * n2), n1 / n2]]))

    @staticmethod
    def mirror():
        return np.matrix(np.array([[1, 0], [0, 1]]))

    @staticmethod
    def shg(n1, n2, n3, d1, d2):  # single harmonic generator
        input_refraction = np.matrix(np.array([[1, 0], [0, n1 / n2]]))
        propagation1 = np.matrix(np.array([[1, d1], [0, 1]]))
        inside_refraction = np.matrix(np.array([[1, 0], [0, n2 / n3]]))
        propagation2 = np.matrix(np.array([[1, d2], [0, 1]]))
        output_refraction = np.matrix(np.array([[1, 0], [0, n3 / n1]]))
        reverse_path = [output_refraction, propagation2, inside_refraction, propagation1, input_refraction]
        return reduce(np.matmul, reverse_path)

    @staticmethod
    def faraday_isolator(n1, n2, d):
        entering_refraction = np.matrix(np.array([[1, 0], [0, n1 / n2]]))
        propagation = np.matrix(np.array([[1, d], [0, 1]]))
        exiting_refraction = np.matrix(np.array([[1, 0], [0, n2 / n1]]))
        reverse_path = [exiting_refraction, propagation, entering_refraction]
        return reduce(np.matmul, reverse_path)

    @staticmethod
    def fabry_perot_cavity(Rc, d):
        curved_mirror_refraction = np.matrix(np.array([[1, 0], [-2/Rc, 1]]))
        propagation = np.matrix(np.array([[1, d], [0, 1]]))
        reverse_path = [curved_mirror_refraction, propagation, curved_mirror_refraction,
                        propagation, curved_mirror_refraction]
        return reduce(np.matmul, reverse_path)
