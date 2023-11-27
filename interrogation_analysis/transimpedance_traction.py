#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   plot_linearity_fbg_transimpedance.py
@Time    :   2023/09/26 11:22:35
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import h5py
import numpy as np
import sympy as sp
import locale
import pandas as pd
import matplotlib.pyplot as plt
import lvm_read
from modeling.math_model_accel import AccelModel
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4

# get mass of seismic mass
seismic_mass = AccelModel().seismic_mass


def traction_and_tension():
    _f = h5py.File('phd_data.hdf5', 'r')
    f = _f['transimpedance_analysis/saturation_test']

    fig, ax = plt.subplots(2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    ax[0].plot(f['time'][:], f['dyn'][:], '.')
    ax[0].set_ylabel(r'Tração$[\si{\newton}]$')
    ax[1].plot(f['time'][:], f['fbg'][:], '.')
    ax[1].set_ylabel(r'Tensão$[\si{\volt}]$')

    fig.supxlabel('Tempo$[\\si{\\second}]$')
    plt.savefig("../dissertacao/images/traction_and_tension.pdf", format="pdf")
    plt.close(fig=1)
    # ax.plot(4.0*(-f['dyn']-0.5)/(seismic_mass*9.8), f['fbg'], ':')
    # ax.set_ylabel(r'$v^\text{tr}[\si{\volt}]$')
    # ax.set_xlabel('Aceleração[g]')
    # plt.show()
