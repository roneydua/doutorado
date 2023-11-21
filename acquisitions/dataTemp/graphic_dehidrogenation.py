#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   graphic_dehidrogenation.py
@Time    :   2023/08/28 17:29:37
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import pandas as pd
import locale
import matplotlib.pyplot as plt

from common_functions.generic_functions import *

# from IPython.core.interactiveshell import InteractiveShell
# from ipywidgets import interactive, fixed
# InteractiveShell.ast_node_interactivity = "all"

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4

# 04/09/2023


def plot230912():
    s_1_before = pd.read_csv('aquisitions/dataTemp/source.txt',
                             delimiter=r'\s+', names=['wavelength', 'power_dbm'])
    fbg_1_before_dehydrogenated = pd.read_csv('aquisitions/dataTemp/fbg_1.txt',
                                              delimiter=r'\s+', names=['wavelength', 'power_dbm'])

    s_1_after_dehydrogenated = pd.read_csv('aquisitions/dataTemp/source_dehydrogenated.txt',
                                           delimiter=r'\s+', names=['wavelength', 'power_dbm'])
    fbg_1_after_dehydrogenated = pd.read_csv(
        'aquisitions/dataTemp/fbg_1_dehydrogenated.txt', delimiter=r'\s+', names=['wavelength', 'power_dbm'])

    # fbg_2 = pd.read_csv('aquisitions/dataTemp/fbg_2.txt',
    # delimiter=r'\s+')
    data_on_production = pd.read_csv(
        '../data/danteAlex/20230912/FBG#1.txt', sep=r'\t', engine='python')
    # fig.clear()
    fig, ax = plt.subplots(1, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))

    ax.set_xlabel(r'$\lambda,\unit{\m}$')
    ax.set_ylim(-0.1, 1)
    # before dehydrogenation
    r = calc_reflectivity_by_transmission(
        s_1_before.power_dbm, fbg_1_before_dehydrogenated.power_dbm, fbg_1_before_dehydrogenated.wavelength, normalize_source=True, min_wavelength=1500e-9, max_wavelength=1520e-9)
    ax.plot(s_1_before.wavelength, r, label='Hidrogenada pós produção')
    # after dehydrogenation
    r = calc_reflectivity_by_transmission(
        s_1_after_dehydrogenated.power_dbm, fbg_1_after_dehydrogenated.power_dbm, fbg_1_after_dehydrogenated.wavelength, normalize_source=True, min_wavelength=1565e-9, max_wavelength=1580e-9)
    ax.plot(s_1_after_dehydrogenated.wavelength, r, label='Desidrogenada')

    r_on_production = calc_reflectivity_by_transmission(
        data_on_production.iloc[:, 1], data_on_production.iloc[:, -1], data_on_production.iloc[:, 0])

    ax.plot(data_on_production.iloc[:, 0],
            r_on_production, label='Durante a produção')
    ax.legend()
    ax.vlines(x=1549.8*1e-9,ymin=-0.1,ymax=1)





s_1_before = pd.read_csv('aquisitions/dataTemp/source.txt',
                         delimiter=r'\s+', names=['wavelength', 'power_dbm'])
fbg_1_before_dehydrogenated = pd.read_csv('aquisitions/dataTemp/fbg_1.txt',
                                          delimiter=r'\s+', names=['wavelength', 'power_dbm'])

s_1_after_dehydrogenated = pd.read_csv('aquisitions/dataTemp/source_dehydrogenated.txt',
                                       delimiter=r'\s+', names=['wavelength', 'power_dbm'])
fbg_1_after_dehydrogenated = pd.read_csv(
    'aquisitions/dataTemp/fbg_1_dehydrogenated.txt', delimiter=r'\s+', names=['wavelength', 'power_dbm'])

fig, ax = plt.subplots(1, 2, num=1 ,sharex=True, figsize=(FIG_L, FIG_A))
fig.clear()
ax[0].plot(s_1_after_dehydrogenated.wavelength,s_1_after_dehydrogenated.power_dbm)
ax[0].plot(fbg_1_after_dehydrogenated.wavelength,fbg_1_after_dehydrogenated.power_dbm)

ax[1].plot(s_1_before.wavelength,s_1_before.power_dbm)

ax[1].plot(fbg_1_before_dehydrogenated.wavelength,fbg_1_before_dehydrogenated.power_dbm)


gabriel()