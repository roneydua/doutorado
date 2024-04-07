#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   accel_1.py
@Time    :   2024/01/23 20:10:54
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""

import numpy as np
import sympy as sp
import locale
import matplotlib.pyplot as plt
import h5py
import allantools

# from IPython.core.interactiveshell import InteractiveShell
# from ipywidgets import interactive, fixed
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4


def plot_accel_1_pm_components():
    f = h5py.File("phd_data.hdf5", "r")
    ff = f["accel_1/test_of_pm_components"]
    fig, ax = plt.subplots(1, 1, num=1, figsize=(FIG_L, FIG_A * 0.75))
    ax.plot(ff["sld/time"][:], ff["sld/v_trans"][:], label="SLD")
    ax.plot(
        ff["sld_isolator/time"][:], ff["sld_isolator/v_trans"][:], label="SLD+Isolator"
    )
    ax.plot(
        ff["sld_isolator_coupler/time"][:],
        ff["sld_isolator_coupler/v_trans"][:],
        label="SLD+Isolator+Acoplador PM 99/1",
    )
    ax.legend()
    ax.set_xlabel(r"Tempo [\unit{\s}]")
    ax.set_ylabel(r"Tensão $v^{\text{tr}}[\unit{\volt}]$ ")
    plt.savefig("../tese/images/accel_1_pm_components.pdf", format="pdf")
    plt.close(fig=1)
    f.close()


def calibration_date_accel_1():
    def calc_mean_std(data: np.ndarray):
        mean = np.mean(data)
        std = np.std(data)
        str_data = (
            locale.format_string("%.3f", mean)
            + r"$\pm$"
            + locale.format_string("%.3f", std)
        )
        return str_data

    f = h5py.File("phd_data.hdf5", "r")
    ff = f.require_group("accel_1/allan_variance")
    fig, ax = plt.subplots(2, 2, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    plt.show()

    ax[0, 0].plot(
        ff["y_up/time"][:],
        ff["y_up/y"][:],
        label="y cima" + "(" + calc_mean_std(ff["y_up/y"][:]) + ")",
    )

    ax[0, 1].plot(
        ff["y_up/time"][:],
        ff["y_up/y"][:] / ff["y_up/tap"][:],
        label="y cima normalizado"
        + "("
        + calc_mean_std(ff["y_up/y"][:] / ff["y_up/tap"][:])
        + ")",
    )

    ax[1, 0].plot(
        ff["y_down/time"][:],
        ff["y_down/y"][:],
        label="y baixo" + "(" + calc_mean_std(ff["y_down/y"][:]) + ")",
    )
    ax[1, 1].plot(
        ff["y_down/time"][:],
        ff["y_down/y"][:] / ff["y_down/tap"][:],
        label="y baixo normalizado"
        + "("
        + calc_mean_std(ff["y_down/y"][:] / ff["y_down/tap"][:])
        + ")",
    )
    ax[0, 0].legend()
    ax[0, 1].legend()
    ax[1, 0].legend()
    ax[1, 1].legend()
    ax[0, 0].set_title("Medidas em Volt")
    ax[0, 1].set_title("Medidas normalizadas")
    fig.supxlabel(r"Tempo[\unit{\second}]")
    plt.savefig("../tese/images/medidas_calibracao.pdf", format="pdf")
    plt.close(fig=1)
    f.close()


def plot_calibrated_accel_1_data():
    a = np.array([[0.572, 1], [0.715, 1]])
    b = np.array([-9.8, 9.8])
    coef = np.linalg.inv(a) @ b

    def calibrate_data(x):
        return coef[0] * x + coef[1]

    f = h5py.File("phd_data.hdf5", "r")
    ff = f.require_group("accel_1/allan_variance")
    fig, ax = plt.subplots(1, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    plt.show()

    ax.plot(
        ff["y_up/time"][:],
        calibrate_data(ff["y_up/y"][:] / ff["y_up/tap"][:]),
        label="y cima",
    )
    ax.plot(
        ff["y_down/time"][:],
        calibrate_data(ff["y_down/y"][:] / ff["y_down/tap"][:]),
        label="y baixo",
    )
    ax.legend()
    ax.set_ylabel(r"Aceleração[\unit{\meter\per\second\squared}]")
    fig.supxlabel(r"Tempo[\unit{\second}]")
    plt.savefig("../tese/images/medidas_calibrada.pdf", format="pdf")
    plt.close(fig=1)
    f.close()


def plot_allan_data_calibrated():
    a = np.array([[0.572, 1], [0.715, 1]])
    b = np.array([-9.8, 9.8])
    coef = np.linalg.inv(a) @ b

    def calibrate_data(x):
        return coef[0] * x + coef[1]

    f = h5py.File("phd_data.hdf5", "r")
    ff = f.require_group("accel_1/allan_variance")
    fig, ax = plt.subplots(3, 1, num=1, sharex=True, figsize=(FIG_L, 1.25 * FIG_A))
    plt.show()
    ax[0].plot(ff["y_up_long/time"][:] / 3600, ff["y_up_long/y"][:])
    ax[0].set_ylabel(r"Tensão y[\unit{\volt}]")
    ax[1].plot(ff["y_up_long/time"][:] / 3600, ff["y_up_long/tap"][:])
    ax[2].plot(
        ff["y_up_long/time"][:] / 3600,
        ff["y_up_long/y_calibrated"][:],
    )
    ax[1].set_ylabel(r"Tensão Tap[\unit{\volt}]")
    ax[2].set_ylabel(r"Aceleração [\unit{\meter\per\second\squared}]")
    plt.savefig("../tese/images/allan_data_calibrated.pdf", format="pdf")
    plt.close(fig=1)
    f.close()


def plot_allan_variace_accel1():
    f = h5py.File("phd_data.hdf5", "a")
    ff = f.require_group("accel_1/allan_variance")
    fig, ax = plt.subplots(1, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    # taus_y, adevs_y, errors, ns = allantools.adev(
    #     ff["y_up_long/y_calibrated"][:], 10, data_type="freq", taus="all"
    # )
    # ff["y_up_long/adevs_y_long"] = adevs_y
    # ff["y_up_long/taus_y_long"] = taus_y
    ax.loglog(ff["y_up_long/taus_y_long"][:], ff["y_up_long/adevs_y_long"][:],label="Completo")
    # taus_y, adevs_y, errors, ns = allantools.adev(
    #     ff["y_up_long/y_calibrated"][2_000_000:2_500_000], 10, data_type="freq", taus="all"
    # )
    # ff["y_up_long/adevs_y"] = adevs_y
    # ff["y_up_long/taus_y"] = taus_y
    ax.loglog(ff["y_up_long/taus_y"][:], ff["y_up_long/adevs_y"][:],label=r"20\unit{\hour}")
    ax.set_ylabel(r"Variância de Allan [\unit{\meter\per\second\squared}]")
    ax.set_xlabel(r"Tempo de correlação [\unit{\second}]")
    ax.legend()
    ax.set_ylim(top=0.1, bottom=0.01)
    ax.set_xlim(left=0.5, right=1000)
    ax.grid(which="both")
    plt.savefig("../tese/images/allan_variace_accel1.pdf", format="pdf")
    plt.close(fig=1)
    f.close()
