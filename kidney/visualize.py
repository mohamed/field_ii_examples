#!/usr/bin/env python3

#%matplotlib inline

import os
import sys
import time
import numpy as np
import matplotlib
import matplotlib.pyplot
import matplotlib.animation
import scipy
import scipy.io
import scipy.signal
import scipy.interpolate

#matplotlib.use("Agg")

####### CONSTANTS #######
speed_of_sound = 1540.0         # Speed of sound (m/s)
sampling_frequency = 10.0e6     # Sampling frequency (Hz)
Ft = 7.0e6                      # Transmit frequency of the PMUT (Hz)
d = speed_of_sound/(2.0*Ft)     # Distance between every two sensors in the array (m)
depth = 0.1289                  # Scan depth (m)
#########################


class RxParams(object):
    """Container class to store receiver parameters"""
    def __init__(self, c, Ft, Fs, N, d, angles, depth, nr_samples):
        self.c = c
        self.Ft = Ft
        self.Tt = 1.0 / self.Ft
        self.Fs = Fs
        self.Ts = 1.0 / self.Fs
        self.N = N
        self.d = d
        self.angles = angles
        self.depth = depth
        self.nr_samples = nr_samples
        self.rs = np.linspace(0.0, depth, nr_samples)
        self.tgc = self.compute_tgc()
        self.elements = self.get_elements_positions()

    def get_elements_positions(self):
        # Elements positions in the x-axis
        dx = [-(self.d/2.0 + i*self.d) for i in range(self.N//2-1, -1, -1)]
        dx.extend([(self.d/2.0 + i*self.d) for i in range(0, self.N//2, 1)])
        return np.array(dx)

    def compute_tgc(self):
        # Time-gain control
        alpha = 4.0e-6
        tgc = np.exp(2.0 * alpha * self.rs * self.Ft)
        return tgc


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def generate_bmode(steering_angles, depths, scan_lines, resolution=[512, 512]):
    x_resolution = resolution[0]
    y_resolution = resolution[-1]

    dep_min = depths[0]
    dep_max = depths[-1]
    lat_min = depths[-1] * np.sin(np.deg2rad(steering_angles[0]))
    lat_max = depths[-1] * np.sin(np.deg2rad(steering_angles[-1]))

    step_x = 1.0 / (x_resolution - 1.)
    step_y = 1.0 / (y_resolution - 1.)
    pos_vec_y_new = np.arange(0.0,1.0 + step_x, step_x) * (lat_max - lat_min) + lat_min
    pos_vec_x_new = np.arange(0.0,1.0 + step_y, step_y) * (dep_max - dep_min) + dep_min

    pos_mat_y_new, pos_mat_x_new = np.meshgrid(pos_vec_y_new, pos_vec_x_new)
    r_cart, th_cart = cart2pol(pos_mat_x_new, pos_mat_y_new)

    F = scipy.interpolate.RegularGridInterpolator((depths, np.deg2rad(steering_angles)), scan_lines, bounds_error=False, fill_value=np.nan)
    pts = np.column_stack((np.ravel(r_cart), np.ravel(th_cart)))
    b_mode = F(pts)
    b_mode = b_mode.reshape(x_resolution, y_resolution)

    mask = np.ones(b_mode.shape)
    mask[np.isnan(b_mode)] = 0
    b_mode[np.isnan(b_mode)] = np.amin(b_mode)

    return (b_mode, pos_vec_y_new, pos_vec_x_new)


def get_angles(start, delta, nr_angles):
    angles = np.array([start + delta * i for i in range(nr_angles)])
    return angles


def beamform_delayline(params, rf_data):
    nr_angles = len(params.angles)
    Beamformed = np.zeros((params.nr_samples, nr_angles))

    for a in range(nr_angles):
        tstart = rf_data[a]['tstart'][0]
        data = rf_data[a]['rf_data']
        delay = np.round(tstart * params.Fs).astype(int)[0]
        dim = delay if delay > 0 else 0
        delays = np.zeros((dim, 1))
        bf = np.concatenate((delays, data))
        Beamformed[0:bf.shape[0],a] = bf[:,0]
    return Beamformed


def post_processing(Beamformed, nr_angles, tgc):
    # DC component removal
    mean_value = np.mean(Beamformed, axis=0)
    bf_data = np.zeros(Beamformed.shape)
    for i in range(nr_angles):
        bf_data[:,i] = Beamformed[:,i] - mean_value[i]

    # TGC
    #bf_data = np.multiply(bf_data.T, tgc).T

    # Envelope detection
    bf_data = np.abs(scipy.signal.hilbert(bf_data, axis=0))

    # Log compression
    norm_value = np.amax(bf_data)
    bf_data /= norm_value
    bf_data = 20.0 * np.log10(bf_data)

    return bf_data


def read_rf_data(dirname, N):
    # Load data from Matlab dataset
    rf_data = {}
    for i in range(N):
        rf_data[i] = scipy.io.loadmat(os.path.join(dirname, 'rf_ln' + str(i + 1)))
    return rf_data


def prepare_frame(idx):
    N = 128
    rf_data = read_rf_data('rf_data', N)
    dtheta = 90.0 / N
    nr_angles = N
    nr_samples = 0
    for v in rf_data.values():
        if len(v['rf_data']) > nr_samples:
            nr_samples = len(v['rf_data'])

    if idx == 0:
        print("Number of angles = %d" % nr_angles)
        print("Number of elements = %d" % N)
        print("Speed of sound = %f" % speed_of_sound)
        print("Sampling Frequency = %f" % sampling_frequency)
        print("Number of samples = %d" % nr_samples)

    params = RxParams(speed_of_sound, Ft, sampling_frequency, N, d, get_angles(-N/2*dtheta, dtheta, nr_angles), depth, nr_samples)
    return (params, rf_data)


def main():
    nr_frames = 1
    print("Number of frames = %d" % nr_frames)
    images = list()
    # Setup visualization stuff here
    fig = matplotlib.pyplot.figure()
    fig.suptitle('B Mode', fontsize=14)
    matplotlib.pyplot.xlabel('Azimuth (mm)', fontsize=12)
    matplotlib.pyplot.ylabel('Depth (mm)', fontsize=12)

    for idx in range(nr_frames):
        tic = time.time()
        print("Processing frame %d" % idx)
        params, rf_data = prepare_frame(idx)
        bf_data = beamform_delayline(params, rf_data)
        bf_data = post_processing(bf_data, len(params.angles), params.tgc)
        b_mode, az_lims, dep_lims = generate_bmode(params.angles, params.rs, bf_data)
        az_lims *= 1000
        dep_lims *= 1000
        print("Processing framed finished. It took %s seconds" % (time.time() - tic))
        im = matplotlib.pyplot.imshow(b_mode, cmap='gray', animated=True, interpolation=None, aspect='equal', vmin=-60, vmax=0, extent=[az_lims[0], az_lims[-1], dep_lims[-1], dep_lims[0]])

    matplotlib.pyplot.show()


if "__main__" == __name__:
    tic = time.time()
    main()
    print("%s took %f seconds" % (sys.argv[0], (time.time() - tic)))
