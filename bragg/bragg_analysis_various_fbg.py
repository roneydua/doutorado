#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   bragg_analysis_various_fbg.py
@Time    :   2023/08/15 18:49:47
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import sympy as sp
import locale
import matplotlib.pyplot as plt
from bragg.bragg import Bragg
from common_functions.generic_functions import calc_laser
# from IPython.core.interactiveshell import InteractiveShell
# from ipywidgets import interactive, fixed
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4
# End of header


def test_with_multiples_fbg():
    fbg1 = Bragg(2e-3,1550.0,5e-4,delta_span_wavelength=1000,diff_of_peak=10)
    fbg2 = Bragg(fbg1.fbg_size,1550.0,fbg1.delta_n,wavelength_span=fbg1.wavelength_span)
    fbg3 = Bragg(fbg1.fbg_size, 1550.0, fbg1.delta_n, wavelength_span=fbg1.wavelength_span,
                 number_of_grating_period_forced=2*fbg1.number_of_grating_period)
    s = calc_laser(fbg1.wavelength_span,fbg1.wavelength_peak,0.001)
    fig.clear()
    fig, ax = plt.subplots(1, 1, num=1 , figsize=(FIG_L, FIG_A))
    ax.plot(fbg1.wavelength_span_nm,fbg1.r0,label=r'$\text{FBG}_1$')
    ax.plot(fbg1.wavelength_span_nm, fbg1.r0*fbg2.r0,
            label=r'$\text{FBG}_1*\text{FBG}_2$')
    ax.plot(fbg1.wavelength_span_nm,fbg3.r0,label=r'$\text{FBG}_3$')
    # ax.plot(fbg1.wavelength_span_nm,r,label=r'$Dupla \text{FBG}$')
    ax.legend()