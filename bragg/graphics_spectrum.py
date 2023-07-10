#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   graphics_spectrum.py
@Time    :   2023/03/17 17:43:02
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''
# type: ignore
from matplotlib import ticker
import numpy as np
import locale
import matplotlib.pyplot as plt
import os
from acquisitionAnritsu import dataAcquisition
from matplotlib import units
from scipy.signal import find_peaks_cwt
import pandas as pd
from pathlib import Path

from tol_colors import tol_cset
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("./../../../../programasComuns/roney3.mplstyle")

cores = plt.rcParams["axes.prop_cycle"].by_key()["color"]

figL = 6.29
figA = (90.0) / 25.4
# End of header

path_locked = "./data/03202023/seismicMass_locked/"
# files2Plot_locked = os.listdir(path_locked)
path_unlocked = "./data/03312023/seismicMass_unlocked/z_up/"
path_unlocked_z_down = "./data/03312023/seismicMass_unlocked/z_down/"
# files2Plot_unlocked = os.listdir(path_unlocked)
graphic = ['A', 'B', 'C', 'D', 'E', 'F', 'tap']
graphic_axis = [ r'$y_-$', r'$z_{+}$', r'$z_-$', r'$y_{+}$', r'$x_{+}$', r'$x_-$', r'\text{\emph{Tap}}']
acquisition_label = [
    r'\emph{Locked}', r'\emph{Unlocked} $\hat{\boldsymbol{z}}\uparrow$',
    r'\emph{Unlocked} $\hat{\boldsymbol{z}}\downarrow$'
]
da_unlocked = []
da_unlocked_z_down = []
da_locked = []
da = [da_locked, da_unlocked, da_unlocked_z_down]

def read_data_from_folder(name, path_to_read):
    # print(index_graphic,graf)
    _da = 0
    for i in os.listdir(path_to_read):
        _da = dataAcquisition(test_name=path_to_read+i)
        if i.split('_')[0] == name:
            break
    return _da


for i in range(len(graphic)):
    da_locked.append(read_data_from_folder(graphic[i],path_locked))
    da_unlocked.append(read_data_from_folder(graphic[i],path_unlocked))
    da_unlocked_z_down.append(read_data_from_folder(graphic[i],path_unlocked_z_down))


def dbm2W(_a):
    return 10.0**(_a * 0.1)



def plot1():
    plt.close(1)
    fig, ax = plt.subplots(7, 1, num=1, sharex='col', figsize=(figL, 1.75*figA),dpi=576)
    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(2))
    def plotGraphics(_ax, _x, _y, _label,line='-'):
        _ax.plot(_x,_y,line,label=_label)
        # _ax.set_ylabel(_name)


    for row in range(7):
        for graf_index in range(3):
            plotGraphics(ax[row], da[graf_index][row].wave_length,
                         da[graf_index][row].power, acquisition_label[graf_index])
            # plotGraphics(ax[row], da[graf_index][row].wave_length,
            # da[graf_index][row].power,
            # graphic_axis[row])


            ax[-1].set_xlabel("Comprimento de onda, " +r"$\si{\nano\meter}$")
            ax[-1].set_xlim(1535, 1555)  # type: ignore

    plt.legend(ncols=3, loc='center',bbox_to_anchor=(0.5, -.75))
    fig.supylabel("Potência, " + r"$\si{\decibelm}$")

    # ax[row,0].set_ylabel(r'Potência,$\si{\decibel m}$')
    plt.savefig("./../../images/acquisition_dbm.pdf", format="pdf",dpi=576)
    plt.close(fig=1)

plot1()


def plot_tap_correction_test(normalized=False):
    plt.close(2)
    fig, ax = plt.subplots(6,
                           1,
                           num=2,
                           sharex='col',
                           figsize=(figL, 1.75 * figA),
                           dpi=144)
    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(2))
    def plotGraphics(_ax, _x, _y, _label, _name, line='-'):
        _ax.plot(_x, _y, line, label=_label)
        _ax.set_ylabel(_name)

    x_tap_normalized = []
    y_tap_normalized = []

    for i in range(3):
        x_tap_normalized.append(da[i][6].wave_length)
        y_tap_normalized.append(dbm2W(da[i][6].power))
        #normalized tap
        # y_tap_normalized[i] -= y_tap_normalized[i].max()
    for row in range(6):
        for graf_index in range(3):
            # find index of wave length with 1535 nm.
            index = np.where(da[graf_index][row].wave_length>1535.0)[0][0]
            # print(index,da[graf_index][row].wave_length[index])
            x = da[graf_index][row].wave_length[index:]
            y = dbm2W(da[graf_index][row].power[index:])
            # normalize all values by tap
            y_normalized = y / y_tap_normalized[graf_index][index:]

            if normalized == True:
                y_normalized /= y_normalized.max()
            plotGraphics(ax[row], x, y_normalized,
                         acquisition_label[graf_index], graphic_axis[row])
            ax[-1].set_xlabel(r"$\lambda,\si{\nano\meter}$")
            ax[-1].set_xlim(1535, 1555)  # type: ignore
        # ax[row].vlines(1546.1,-1,2)
        # ax[row].set_ylim(0.45,1.25)

    plt.legend(ncols=3, loc='center', bbox_to_anchor=(0.5, -.75))
    # fig.supylabel(r'$\frac{p_{\text{measured}}}{p_{\text{tap}}}$')
    # ax[row,0].set_ylabel(r'Potência,$\si{\decibel m}$')
    if normalized:
        plt.savefig("./../../images/acquisition_tap_normalized.eps", format="eps")
    else:

        plt.savefig("./../../images/acquisition_tap.eps", format="eps")
    plt.close(fig=2)


# plot_tap_correction_test(True)
# plot_tap_correction_test()



def plot_tap_correction(normalized=False):
    plt.close(2)
    fig, ax = plt.subplots(6,
                           1,
                           num=2,
                           sharex='col',
                           figsize=(figL, 1.75 * figA),
                           dpi=144)
    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(2))

    def plotGraphics(_ax, _x, _y, _label, _name, line='-'):
        _ax.plot(_x, _y, line, label=_label)
        _ax.set_ylabel(_name)

    x_tap = []
    y_tap = []

    for i in range(3):
        x_tap.append(da[i][6].wave_length)
        y_tap.append(da[i][6].power_watt)
        #normalized tap
        # y_tap[i] -= y_tap[i].max()
    for row in range(6):
        for graf_index in range(3):
            # find index of wave length with 1535 nm.
            index = np.where(da[graf_index][row].wave_length > 1535)[0][0]
            # print(index,da[graf_index][row].wave_length[index])
            x = da[graf_index][row].wave_length[index:]
            y = da[graf_index][row].power_watt[index:]
            y_normalized = y / y_tap[graf_index][index:]

            if normalized == True:
                if graf_index > 0:
                    # now we use the point on index to make same power to compare
                    y_normalized /= y_normalized[0]
                    y_normalized *= (da[0][row].power_watt[index] / y_tap[0][0])
            plotGraphics(ax[row], x, y_normalized,
                         acquisition_label[graf_index], graphic_axis[row])
            ax[-1].set_xlabel(r"$\lambda,\si{\nano\meter}$")
            ax[-1].set_xlim(1535, 1555)  # type: ignore
        # ax[row].vlines(1546.1,-1,2)
        # ax[row].set_ylim(0.45,1.25)

    plt.legend(ncols=3, loc='center', bbox_to_anchor=(0.5, -.75))
    # fig.supylabel(r'$\frac{p_{\text{measured}}}{p_{\text{tap}}}$')
    # ax[row,0].set_ylabel(r'Potência,$\si{\decibel m}$')
    if normalized:
        plt.savefig("./../../images/acquisition_tap_normalized.eps",
                    format="eps")
    else:

        plt.savefig("./../../images/acquisition_tap.eps", format="eps")
    plt.close(fig=2)


# plot_tap_correction()
# plot_tap_correction(True)


## Plot Graphics with pandas

r_e = pd.read_csv(Path('./../../../../experimentos/20042023/e_reflection_042023__142133.csv'))
f_4_percent = pd.read_csv(Path('./../../../../experimentos/20042023/e_reflection_4_percent042023__142133.csv'))

fig, ax = plt.subplots(2,1 , num=1, sharex=True)
ax[0].plot(r_e['wave_length'],f_4_percent['power'],lw=0.75)
ax[0].plot(r_e['wave_length'],r_e['power'],lw=0.75)
ax[0].plot(r_e['wave_length'], r_e['power'] - f_4_percent['power'], lw=0.75)


ax[1].plot(r_e['wave_length'],4.0/(10**((f_4_percent['power'] - r_e['power'])*0.1-1.0)))
