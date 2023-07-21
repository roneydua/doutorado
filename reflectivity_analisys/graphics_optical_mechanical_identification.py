#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   graphics_optical_mechanical_identification.py
@Time    :   2023/07/17 18:34:34
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

from modeling.mathModelAccel import AccelModel
import locale
from debugpy import trace_this_thread

import h5py
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.linalg import eig

from common_functions.common_functions import *

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# plt.style.use("default")
plt.style.use("./common_functions/roney3.mplstyle")
FIG_L = 6.29
FIG_A = (90.0) / 25.4
am = AccelModel()
max_traction_10g = am.seismic_mass * 10.0 * 9.89 / 4.0
# load data


def plot_colleted_data(fbg_number: str):
    # laod hdf5 file with source and fbg
    f = h5py.File('./../data/phd_data.hdf5',
                  'r')['optical_mechanical_identification']
    fbg_data = []
    source_data = []
    for key in f.keys():
        if key.startswith("fbg_2"):
            fbg_data.append(key)
        else:
            source_data.append(key)
    fig, ax = plt.subplots(1, 1, num=1, figsize=(FIG_L, FIG_A))
    for fbg in fbg_data:
        ax.plot(
            f[fbg]['wavelength'][:],
            f[fbg]['power_dbm'][:],
            lw=0.5,
            label=locale.format_string('T=%.2f', f[fbg].attrs['traction_N']) + r"$\si{\newton}$, " + locale.format_string('T=%.2f', f[fbg].attrs['micrometer_position_um']) + r"$\si{\um}$")
    ax.set_xlim(1530, 1560)
    ax.set_ylim(bottom=-42, top=-32)
    ax.legend(ncols=2)
    ax.set_xlabel(r"$\lambda,\unit{\nm}$")
    ax.set_ylabel(r"$\unit{\dbm\per\nm}$")
    plt.savefig(
        "../dissertacao/images/plot_colleted_data_fbg_2.pdf",
        format="pdf")
    plt.close(fig=1)
    f.close()


def plot_reflectivity_of_colleted_data_fbg_3():
    f = h5py.File('./../data/phd_data.hdf5',
                  'r')['optical_mechanical_identification']
    fbg_data = []
    source_data = []
    for key in f.keys():
        if key.startswith("fbg_2"):
            fbg_data.append(key)
        elif key.startswith("fonte_2"):
            source_data.append(key)
    fig.clear()
    fig, ax = plt.subplots(1, 1, num=1, figsize=(FIG_L, FIG_A))
    for fbg in fbg_data:
        diff = (f[source_data[0]]['power_dbm'][2600:3000] -
                f[fbg]['power_dbm'][2600:3000]).mean()
        # fix bias between source and fbg collected data
        source_fix = f[source_data[0]]['power_dbm'] - diff
        # compute reflectivity
        r = 1.0 - 10.0**(0.1 * (f[fbg]['power_dbm'][:] - source_fix))
        ax.plot(f[source_data[0]]['wavelength'][:], r, lw=0.5,
                label=locale.format_string('T=%.2f', f[fbg].attrs['traction_N']) + r"$\si{\newton}$")
        ax.set_xlim(1520, 1560)
        ax.set_ylim(bottom=-0.1, top=1)
        # ax.legend(ncols=2)
        ax.set_xlabel(r"$\lambda,\unit{\nm}$")
        ax.set_ylabel(r"$\unit{\dbm\per\nm}$")
        plt.savefig(
            "../dissertacao/images/plot_reflectivity_of_colleted_data_fbg_3.pdf",
            format="pdf")
        plt.close(fig=1)


def generate_toy_date():
    # estimation of E and l0
    f = h5py.File(
        './../data/phd_data.hdf5',
        'r')['optical_mechanical_identification/test_fibers/fbg_2/test_002']

    fbg_data = []
    for key in f.keys():
        if key.startswith("fbg_2"):
            fbg_data.append(key)

    _N = len(fbg_data)
    mat_a = np.zeros((_N, 2))
    fo_area = (125e-6)**2 * np.pi * 0.25
    tractions = np.zeros(_N)
    deformations = np.zeros(_N)
    _E = 60e9
    _l0 = 21e-3
    for i in range(len(fbg_data)):
        deformations[i] = f[fbg_data[i]].attrs['micrometer_position_um']
        tractions[i] = _E * fo_area / _l0 * deformations[i] * 1e-7
        mat_a[i, 0] = tractions[i]
        mat_a[i, 1] = -fo_area * deformations[i] * 1e-6 * 1e9
    mat_a_T_mat_a = mat_a.T @ mat_a
    eigen_values, eigen_vectors = np.linalg.eig(mat_a_T_mat_a)
    eigen_vector_min = eigen_vectors[:, np.argmin(eigen_values)]
    eigen_vector_min / min(eigen_values)


def calc_optical_mechanical():
    # estimation of E and l0
    f = h5py.File(
        './../data/phd_data.hdf5',
        'r')['optical_mechanical_identification/test_fibers/fbg_2/test_002']

    fbg_data = []
    for key in f.keys():
        if key.startswith("fbg_2"):
            fbg_data.append(key)
    _N = len(fbg_data)
    mat_a = np.zeros((_N, 2))

    tractions = np.zeros(_N)
    deformations = np.zeros(_N)
    for i in range(len(fbg_data)):
        mat_a[i, 0] = -f[fbg_data[i]].attrs['traction_N']
        mat_a[i, 1] = -fo_area * f[fbg_data[i]
                                   ].attrs['micrometer_position_um'] * 1e-6 * 1e9
        tractions[i] = f[fbg_data[i]].attrs['traction_N']
        deformations[i] = f[fbg_data[i]].attrs['micrometer_position_um']
    mat_a_T_mat_a = mat_a.T @ mat_a
    eigen_values, eigen_vectors = np.linalg.eig(mat_a_T_mat_a)
    eigen_vector_min = eigen_vectors[:, np.argmin(eigen_values)]
    mat_a @ eigen_vector_min * min(eigen_values)
    a = 70.0 / eigen_vector_min[1]
    eigen_vector_min[0] * a


def traction_vs_deformation():
    data = h5py.File('./../data/phd_data.hdf5', 'r')['optical_mechanical_identification/test_fibers/fbg_6/']
    tests = ['test_003']
    fo_area = (125.4e-6)**2 * np.pi * 0.25
    l0 = 39.7e-3
    for t in tests:
        traction = np.array([])
        delta_l = np.array([])
        mat_a = np.array([])
        vec_y = np.array([])
        _N = 0
        for key in data[t].keys():
            if key.startswith('fbg_6'):
                traction = np.append(
                    traction, -data[t][key].attrs['traction_N'])
                delta_l = np.append(
                    delta_l, data[t][key].attrs['micrometer_position_um'])
                _N += 1
            mat_a = np.append(mat_a, fo_area*(delta_l[-1])*1e-6/l0)
            vec_y = np.append(vec_y,traction[-1])
        E = 1./(mat_a.T @ mat_a) * mat_a.T @ vec_y
        print(1e-9*E)

        plt.plot(1e-6*delta_l/l0, traction)
        plt.plot(1e-6 * delta_l / l0, fo_area * E * delta_l * 1e-6 / l0, '*')
        
        
