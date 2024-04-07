#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   graphics_accel_v1.py
@Time    :   2023/12/21 09:56:48
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""

import numpy as np
import sympy as sp
import pandas as pd
import locale
import matplotlib.pyplot as plt
import h5py
import lvm_read
import allantools
import scipy.signal

FIG_L = 6.29
FIG_A = (90.0) / 25.4

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
# plt.style.use("fast")

plt.style.use("common_functions/roney3.mplstyle")
# my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]e.


# End of header
def plot_reflectivity_data_accel_1():
    f = h5py.File("acquisitions/file_temp.hdf5", "r")
    ff = f.require_group("reflection_with_seismic_mass_unlocked")
    fig, ax = plt.subplots(1, 1, num=1, figsize=(FIG_L, FIG_A))
    for i in ff.keys():
        if i == "y_minus" or i == "y_plus":
            ax.plot(ff[i]["wavelength_m"][:], ff[i]["optical_power_dbm"][:], label=i)
    ax.legend()
    plt.show()


def lvm_to_hdf5_accel_1():
    d_down = pd.read_csv(
        "non_calibrated_y_down.lvm", names=["y", "tap"], sep=r"\t", engine="python"
    )
    d_up = pd.read_csv(
        "non_calibrated_y_up.lvm", names=["y", "tap"], sep=r"\t", engine="python"
    )
    fig, ax = plt.subplots(2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A), dpi=72)
    dt = 10e-3
    _N = len(d_up["y"]) + len(d_down["y"])
    t_down = np.arange(0, len(d_down["y"]) * dt, dt)
    t_up = np.arange(0, len(d_up["y"]) * dt, dt)

    ax[0].plot(t_down, d_down["y"] / d_down["tap"])
    ax[0].plot(t_up, d_up["y"] / d_up["tap"])
    mat_A = np.ones((_N, 2))
    vec_y = np.zeros(_N)
    mat_A[:, 0] = pd.concat((d_up["y"] / d_up["tap"], d_down["y"] / d_down["tap"]))
    vec_y = pd.concat((9.80 + 0.0 * d_up["y"], -9.80 + 0.0 * d_down["y"]))
    x = np.linalg.inv(mat_A.T @ mat_A) @ mat_A.T @ vec_y
    ax[1].plot(t_down, x[0] * d_down["y"] / d_down["tap"] + x[1])
    ax[1].plot(t_up, x[0] * d_up["y"] / d_up["tap"] + x[1])
    ax[1].set_xlabel(r"Tempo [$\si{\second}$]")
    # e = (x[0] * d_down["y"] + x[1]).mean()
    #     (x[0] * d_up["y"] + x[1]).mean()
    d = pd.read_csv(
        "non_calibrated_y_down_allan.lvm",
        names=["y", "tap"],
        sep=r"\t",
        engine="python",
    )
    d_calibrated = x[0] * (d["y"] / d["tap"]) + x[1]
    plt.plot(d_calibrated)
    plt.plot(d_calibrated * 0 + d_calibrated.mean())
    # np.savetxt("data_accel_1.txt",d_calibrated)
    taus, adevs, errors, ns = allantools.adev(d_calibrated, 100, data_type="freq")
    plt.loglog(taus, adevs)
    # project filter
    # b, a = scipy.signal.iirfilter(N=4,Wn=10,fs=100,btype='low')
    # d_calibrated_filtered = scipy.signal.lfilter(b,a,d_calibrated)
    # plt.plot(d_calibrated_filtered)
    # taus, adevs, errors, ns = allantools.adev(d_calibrated_filtered,100,data_type="freq")
    # plt.loglog(taus,adevs)


def plot_allan_2():
    d = pd.read_csv(
        "20231228/accel_y_allan.txt",
        names=["y", "tap"],
        sep=r"\t",
        engine="python",
    )
    fig, ax = plt.subplots(2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    ax[0].plot(d["y"])
    ax[1].plot(d["tap"])
    fig2, ax2 = plt.subplots(2, 1, num=2, sharex=True, figsize=(FIG_L, FIG_A))
    taus_tap, adevs_tap, errors, ns = allantools.adev(
        d["tap"][:], 100, data_type="freq"
    )
    ax2[0].semilogx(taus_tap, adevs_tap)
    taus_y, adevs_y, errors, ns = allantools.adev(d["y"][:], 100, data_type="freq")
    ax2[1].semilogx(taus_y, adevs_y)


def plot_allan_20240117():
    y_down = pd.read_csv(
        "20240117/y_baixo",
        names=["y", "tap"],
        sep=r"\t",
        engine="python",
    )
    y_up = pd.read_csv(
        "20240117/y_cima",
        names=["y", "tap"],
        sep=r"\t",
        engine="python",
    )
    y_allan = pd.read_csv(
        "20240117/allan_variance_longo",
        names=["time", "y", "tap"],
        sep=r",",
        engine="python",
    )
    a = 19.6 / (y_up["y"].mean() - 0.01 - y_down["y"].mean())
    b = -0.5 * a * (y_up["y"].mean() - 0.01 + y_down["y"].mean())
    y_allan_calibrated = a * y_allan["y"] + b
    # time_in_hours = pd.to_datetime(y_allan["time"], unit="s")
    fig, ax = plt.subplots(2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    ax[0].plot(y_up["y"])
    ax[0].plot(y_down["y"])
    # ax.plot(y_allan["time"] / 3600.0 / 24.0, y_allan["y"])
    ax[0].plot(y_allan["time"] / 3600.0 / 24.0, y_allan_calibrated)
    ax[1].plot(y_allan["y"])
    ax[1].plot(y_allan["y"] / y_allan["tap"])
    ax[1].plot(y_up["y"])
    ax[1].plot(y_down["y"])

    fig2, ax2 = plt.subplots(1, 1, num=2, sharex=True, figsize=(FIG_L, FIG_A))
    taus_tap, adevs_tap, errors, ns = allantools.adev(
        y_allan_calibrated[2_000_000:2_500_000].to_numpy(), 10, data_type="freq",taus='all'
    )
    ax2.loglog(taus_tap, adevs_tap)
    # taus_y, adevs_y, errors, ns = allantools.adev(d["y"][:], 100, data_type="freq")
    # ax2[1].semilogx(taus_y, adevs_y)
    # plt.plot(y_allan["y"])
    # plt.savefig('y_allan.png')


plot_allan_20240117()
