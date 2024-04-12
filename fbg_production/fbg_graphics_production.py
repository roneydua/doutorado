#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   fbg_graphics_production.py
@Time    :   2023/10/02 10:08:38
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""

import locale

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import regex
from natsort import natsorted
from scipy.signal import butter, lfilter

from common_functions.generic_functions import find_index_of_x_span

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = FIG_L / 1.6
# header


def plot_production_one_column(date: str) -> None:
    f = h5py.File("production_file.hdf5", "r")
    ff = f["fbg_production/" + date]
    # ff = f[date]
    # find index of fbgs
    fbg_keys = list(ff.keys())
    number_of_graphics = len(fbg_keys)
    print(number_of_graphics)
    if number_of_graphics > 5:
        row_graphic = 5
        col_graphics = number_of_graphics // 5
    else:
        row_graphic = number_of_graphics
        col_graphics = 1
    fig, ax = plt.subplots(
        row_graphic, col_graphics, num=1, sharex=True, figsize=(FIG_L, 11)
    )
    i = 0
    for ax_count in ax.flat:
        index_min, index_max = find_index_of_x_span(
            1540e-9, 1560e-9, ff[fbg_keys[i] + "/wavelength_m"][:]
        )
        ax_count.set_title(fbg_keys[i])
        ax_count.plot(
            ff[fbg_keys[i] + "/wavelength_m"][:] * 1e9,
            ff[fbg_keys[i] + "/reflectivity"][:, -1],
        )
        ax_count.set_xlim(
            ff[fbg_keys[i] + "/wavelength_m"][index_min] * 1e9,
            ff[fbg_keys[i] + "/wavelength_m"][index_max] * 1e9,
        )
        ax_count.set_ylim(0, 1)
        i += 1

    fig.supylabel("Refletividade")

    _date = pd.to_datetime(date).strftime("%d de %B de %Y")
    fig.suptitle("FBG produzidas em " + _date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$")
    f.close()
    plt.savefig("fbg_production/" + date + ".pdf", format="pdf")
    plt.close(fig=1)


# plot_production_one_column('20231010')


def loop_graphics():
    f = h5py.File("production_files.hdf5", "r")
    ff = f["fbg_production/"]
    fbg_keys = list(ff.keys())
    for i in fbg_keys:
        plot_production_one_column(date=i)
    f.close()


def plot_one_graphic(date: str, fbg_keys=[]):
    f = h5py.File("production_files.hdf5", "r")
    # ff = f[date]
    ff = f["fbg_production/" + date]
    # find index of fbgs
    append_name = ""
    if len(fbg_keys) == 0:
        fbg_keys = natsorted(list(ff.keys()))
    else:
        for i in fbg_keys:
            append_name += "_" + i + "_"
    number_of_graphics = len(fbg_keys)
    fig, ax = plt.subplots(1, 1, num=1, sharex=True, figsize=(FIG_L, FIG_L / 2))
    i = 0
    for i in fbg_keys:
        if i not in [
            "fbg5",
            "fbg15",
            "fbg1",
            "fbg13",
            "fbg14",
            "fbg8",  #
            "fbg6",
            "fbg9",
            "fbg2",
            "fbg16",
            "fbg3",
            "fbg17",
        ]:
            index_min, index_max = find_index_of_x_span(
                1540e-9, 1560e-9, ff[i + "/wavelength_m"][:]
            )
            # ax_count.set_title(i)
            ax.plot(
                ff[i + "/wavelength_m"][:] * 1e9,
                ff[i + "/reflectivity"][:, -1],
                label=i,
            )
            ax.set_xlim(
                ff[i + "/wavelength_m"][index_min] * 1e9,
                ff[i + "/wavelength_m"][index_max] * 1e9,
            )
            ax.set_ylim(0, 1)
    fig.supylabel("Refletividade")
    _date = pd.to_datetime(date).strftime("%d de %B de %Y")
    # fig.suptitle("FBG produzidas em " + _date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$")
    ax.legend(ncols=2)
    # ax.legend([r"${FBG}_1$",r"${FBG}_2$"],ncols=2)
    f.close()
    # plt.show()
    plt.savefig(
        "fbg_production/" + date + append_name + ".png", format="png", transparent=False
    )
    plt.close(fig=1)


# plot_one_graphic("20231130",fbg_keys=["fbg4","fbg9"])
plot_one_graphic(date="20240207")


def plot_one_graphic_2(date="20240328", fbg_keys=[]):
    f = h5py.File("fbg_production/dado_davi_20240304.hdf5", "r")
    # ff = f[date]
    ff = f["20240328"]
    # find index of fbgs
    append_name = ""
    if len(fbg_keys) == 0:
        fbg_keys = natsorted(list(ff.keys()))
    else:
        for i in fbg_keys:
            append_name += "_" + i + "_"
    number_of_graphics = len(fbg_keys)
    fig, ax = plt.subplots(1, 1, num=1, sharex=True, figsize=(FIG_L, FIG_L / 2))
    i = 0
    for i in fbg_keys:
        index_min, index_max = find_index_of_x_span(
            1540e-9, 1560e-9, ff[i + "/wavelength_m"][:]
        )
        # ax_count.set_title(i)
        ax.plot(
            ff[i + "/wavelength_m"][:] * 1e9,
            ff[i + "/reflectivity"][:, -1],
            label=i,
        )
        ax.set_xlim(
            ff[i + "/wavelength_m"][index_min] * 1e9,
            ff[i + "/wavelength_m"][index_max] * 1e9,
        )
        ax.set_ylim(0, 1)
    fig.supylabel("Refletividade")
    _date = pd.to_datetime(date).strftime("%d de %B de %Y")
    # fig.suptitle("FBG produzidas em " + _date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$")
    ax.legend(ncols=2)
    # ax.legend([r"${FBG}_1$",r"${FBG}_2$"],ncols=2)
    f.close()
    plt.show()
    plt.savefig(
        "fbg_production/" + date + append_name + ".png", format="png", transparent=False
    )
    plt.close(fig=1)


def remove_wrong_data(file_name: str, last_index: int, data_name: str):
    f = h5py.File(file_name, mode="a")
    new_op = f[data_name + "/optical_power_dbm"][:]
    new_r = f[data_name + "/reflectivity"][:]
    del f[data_name + "/optical_power_dbm"]
    del f[data_name + "/reflectivity"]
    f[data_name + "/optical_power_dbm"] = new_op[:, :last_index]
    f[data_name + "/reflectivity"] = new_r[:, :last_index]
    f.close()


# remove_wrong_data("production_files.hdf5",25,"fbg_production/20240207/fbg3")


def plot_graphics_with_pairs_acc_3():
    """plot_graphics_with_pairs_acc_3 plot graphics of fbgs used on acc 3."""

    def plot_fbg(ax: plt.axes, fbg_number: str):
        ax.plot(
            ff[fbg_number + "/wavelength_m"][:] * 1e6,
            ff[fbg_number + "/reflectivity"][:, -1],
            label=fbg_number,
        )
        ax.legend()

    f = h5py.File(name="./production_files.hdf5", mode="r")
    ff = f["fbg_production/20240207"]
    fig, ax = plt.subplots(3, 1, num=2, sharex=True, figsize=(FIG_L, FIG_L))
    # fbg5 fbg15
    fig.supxlabel(r"$\lambda [\unit{\nm}]$")
    fig.supylabel(r"Refletividade $[\unit{\percent}]$")
    plot_fbg(ax[0], fbg_number="fbg5")
    plot_fbg(ax[0], fbg_number="fbg15")
    plot_fbg(ax[1], fbg_number="fbg8")
    plot_fbg(ax[1], fbg_number="fbg14")
    plot_fbg(ax[2], fbg_number="fbg1")
    plot_fbg(ax[2], fbg_number="fbg13")
    plt.show()
    plt.savefig("../tese/images/fbg_acc_3.pdf", format="pdf")
    plt.close(fig=2)


def plot_graphics_with_pairs_acc_4():
    """plot_graphics_with_pairs_acc_3 plot graphics of fbgs used on acc 3."""

    def plot_fbg(ax: plt.axes, fbg_number: str, ff):
        _fbg_name = regex.findall(r"[a-zA-Z]+", fbg_number)[0]
        _fbg_number = regex.findall(r"(\d+)", fbg_number)[0]

        ax.plot(
            ff[fbg_number + "/wavelength_m"][:] * 1e6,
            ff[fbg_number + "/reflectivity"][:, -1],
            label=_fbg_name.upper() + " " + _fbg_number,
        )
        ax.legend()

    f = h5py.File("./production_files.hdf5", "r")
    ff = f["fbg_production/20240207"]
    fff = f["fbg_production/20231130"]

    plt.close(fig=2)
    fig, ax = plt.subplots(3, 1, num=2, sharex=True, figsize=(FIG_L, FIG_L * 0.75))
    # fbg5 fbg15
    fig.supxlabel(r"$\lambda [\unit{\nm}]$")
    fig.supylabel(r"Refletividade $[\unit{\percent}]$")
    plot_fbg(ax[0], fbg_number="fbg3", ff=ff)
    plot_fbg(ax[0], fbg_number="fbg11", ff=ff)
    plot_fbg(ax[1], fbg_number="fbg7", ff=fff)
    plot_fbg(ax[1], fbg_number="fbg9", ff=ff)
    plot_fbg(ax[2], fbg_number="fbg2", ff=ff)
    plot_fbg(ax[2], fbg_number="fbg17", ff=ff)
    ax[0].set_ylabel("x")
    ax[1].set_ylabel("y")
    ax[2].set_ylabel("z")
    # plt.show()
    plt.savefig("../tese/images/fbg_acc_4.pdf", format="pdf")
    plt.close()


def plot_graphics_with_pairs_acc_5():
    """plot_graphics_with_pairs_acc_3 plot graphics of fbgs used on acc 3."""

    def plot_fbg(ax: plt.axes, fbg_number: str, ff):
        _fbg_name = regex.findall(r"[a-zA-Z]+", fbg_number)[0]
        _fbg_number = regex.findall(r"(\d+)", fbg_number)[0]

        ax.plot(
            ff[fbg_number + "/wavelength_m"][:] * 1e9,
            ff[fbg_number + "/reflectivity"][:, -1],
            label=_fbg_name.upper() + " " + _fbg_number,
        )
        ax.legend()

    f = h5py.File("./production_files.hdf5", "r")
    # ff = f["fbg_production/20240207"]
    # fff = f["fbg_production/20231130"]
    ff = f["fbg_production/20240328"]
    plt.close(fig=2)
    # fig.clear()
    fig, ax = plt.subplots(3, 1, num=5, sharex=True, figsize=(FIG_L, FIG_A))
    # fbg5 fbg15

    fig.supxlabel(r"$\lambda [\unit{\nm}]$")
    fig.supylabel(r"Refletividade $[\unit{\percent}]$")
    plot_fbg(ax[0], fbg_number="fbg3", ff=ff)
    plot_fbg(ax[0], fbg_number="fbg17", ff=ff)
    plot_fbg(ax[1], fbg_number="fbg2", ff=ff)
    plot_fbg(ax[1], fbg_number="fbg16", ff=ff)
    plot_fbg(ax[2], fbg_number="fbg4", ff=ff)
    plot_fbg(ax[2], fbg_number="fbg14", ff=ff)
    ax[0].set_ylabel("x")
    ax[1].set_ylabel("y")
    ax[2].set_ylabel("z")
    ax[2].set_xlim(left=1545, right=1560)
    plt.show()
    plt.savefig("../tese/images/fbg_acc_5.pdf", format="pdf")
    plt.close()


def plot_graphics_with_pairs_acc_6():
    """plot_graphics_with_pairs_acc_3 plot graphics of fbgs used on acc 3."""

    def plot_fbg(ax: plt.axes, fbg_number: str, ff):
        _fbg_name = regex.findall(r"[a-zA-Z]+", fbg_number)[0]
        _fbg_number = regex.findall(r"(\d+)", fbg_number)[0]
        # argmax = ff[fbg_number]["reflectivity"][:, -1].argmax()
        # w_peak = ff[fbg_number]["wavelength_m"][argmax]
        # ax.vlines(w_peak*1e9,1,0)
        ax.plot(
            ff[fbg_number + "/wavelength_m"][:] * 1e9,
            ff[fbg_number + "/reflectivity"][:, -1],
            label=_fbg_name.upper() + " " + _fbg_number,
        )
        ax.legend()

    f = h5py.File("./production_files.hdf5", "r")
    ff = f["fbg_production/20240328"]
    fff = f["fbg_production/20231130"]
    ffff = f["fbg_production/20240207"]
    plt.close(fig=5)
    fig1.clear()
    fig1.set_dpi(144)
    fig1, ax1 = plt.subplots(1, 1, num=7, sharex=True, figsize=(FIG_L, FIG_A))
    # fbg5 fbg15
    for i in [
        5,
        # 6,
        7,
        8,
        9,
        # 10,
        # 11,
        # 12,
        # 13,
        # 15
        18,
    ]:
        plot_fbg(ax1, fbg_number="fbg" + str(i), ff=ff)
    plot_fbg(ax1, fbg_number="fbg14", ff=fff)
    for i in [4, 7, 10, 12]:
        plot_fbg(ax1, fbg_number="fbg" + str(i), ff=ffff)

    fig.clear()
    # sn.set_palette("muted")
    # sn.set_palette("husl", 4)
    fig, ax = plt.subplots(3, 1, num=6, sharex=True, figsize=(FIG_L, FIG_A))
    fig.supxlabel(r"$\lambda [\unit{\nm}]$")
    fig.supylabel(r"Refletividade $[\unit{\percent}]$")
    plot_fbg(ax=ax[0], fbg_number="fbg15", ff=ff)
    plot_fbg(ax=ax[0], fbg_number="fbg12", ff=ff)
    plot_fbg(ax=ax[1], fbg_number="fbg10", ff=ff)
    plot_fbg(ax=ax[1], fbg_number="fbg13", ff=ff)
    plot_fbg(ax=ax[2], fbg_number="fbg6", ff=ff)
    plot_fbg(ax=ax[2], fbg_number="fbg11", ff=ff)
    ax[0].set_ylabel("x")
    ax[1].set_ylabel("y")
    ax[2].set_ylabel("z")
    ax[2].set_xlim(left=1545, right=1560)
    # plt.show()
    plt.savefig("../tese/images/fbg_acc_6.pdf", format="pdf")
    plt.close()
