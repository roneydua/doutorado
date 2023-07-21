#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   Untitled-1
@Time    :   2023/05/29 10:11:01
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import sympy as sp
from pathlib import Path
import pandas as pd
import locale
import matplotlib.pyplot as plt
from IPython.core.interactiveshell import InteractiveShell
from ipywidgets import interactive, fixed

# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# plt.style.use("default")
plt.style.use("./common_functions/roney3.mplstyle")
FIG_L = 6.29
FIG_A = (90.0) / 25.4
# End of header

lf = pd.read_csv(Path(
    "../../../../experimentos/24042023/reflectivity_approximations_and_peaks.csv"
),
                 index_col=0)

p11 = 0.113
p12 = 0.252
v = 0.16
n_eff = 1.482
pe = n_eff**2 * (p12 - v * (p11 + p12)) / 2.

# lf.r_a[0]*lf.ref

k_tr, e = sp.symbols(r'\transimpedanceGain, \deformation')
leg = [r'y_{-}', r'z_{+}', r'z_{-}', r'y_{+}', r'x_{+}', r'x_{-}']

s_i = lambda _i: 2.790 * _i - 0.393
s = s_i(.470)
k_pd = 0.85
l_laser = 1549.3

max_deformation = 3.4265e-5

v_tr = lambda ind: k_tr * k_pd * s * (
    lf.r_a[ind] * lf.reflectivity_peak_wavelength[ind] * (1.0 - pe) * e +
    (lf.r_a[ind] * l_laser + lf.r_l[ind])) / 24.0

f = []
for i in range(6):
    f.append(v_tr(i).expand())
    print(f[-1])



k_tr_values = []
str_to_latex = ""
for i in range(6):
    f[i] = f[i].subs(e, 0)

    k_tr_values.append(sp.solve(f[i] - 1500))
    str_to_latex += "\\transimpedanceGain_{" + leg[i] + r"} & =\num{" + str(k_tr_values[-1][0]) + "} \\unit{\\volt\\per\\ampere}, \\\\"


print(str_to_latex)


def calc_transimpedance_tension_vs_strain_coefficients():
    v_tr_with_gain = lambda ind: k_tr_values[ind][0] * k_pd * s * (
        lf.r_a[ind] * lf.reflectivity_peak_wavelength[ind] * (1.0 - pe) * e +
        (lf.r_a[ind] * l_laser + lf.r_l[ind])) / 24.0

    f_with_gain = []
    str_to_latex = ""

    transimpedance_tension_vs_strain_coefficients = pd.DataFrame()
    transimpedance_tension_vs_strain_coefficients['coef_ang'] = np.zeros(6,dtype=np.float32)
    transimpedance_tension_vs_strain_coefficients['coef_lin'] = np.zeros(6,dtype=np.float32)


    for i in range(6):
        f_with_gain.append(v_tr_with_gain(i).expand())
        # print(f_with_gain[-1])
        _t = sp.collect(f_with_gain[-1],e,evaluate=False)
        transimpedance_tension_vs_strain_coefficients['coef_ang'][i] = _t[e]
        transimpedance_tension_vs_strain_coefficients['coef_lin'][i] = _t[1]
        str_to_latex += "\\transimpedanceTension_{" + leg[i] + r"} & =\num{" + str(
            _t[e]) + "} \\deformation+\\num{"+str(_t[1])+"}, \\\\"
    print(str_to_latex)
    transimpedance_tension_vs_strain_coefficients.to_csv(
        'data/transimpedance_tension_vs_strain_coefficients.csv')


# Constant Gain

def calc_transimpedance_tension_vs_strain_coefficients_constant_gain():
    v_tr_with_gain = lambda ind: 230300 * k_pd * s * (
        lf.r_a[ind] * lf.reflectivity_peak_wavelength[ind] * (1.0 - pe) * e +
        (lf.r_a[ind] * l_laser + lf.r_l[ind])) / 24.0
    f_with_gain = []

    transimpedance_tension_vs_strain_coefficients_constant_gain = pd.DataFrame()
    transimpedance_tension_vs_strain_coefficients_constant_gain['coef_ang'] = np.zeros(6,dtype=np.float32)
    transimpedance_tension_vs_strain_coefficients_constant_gain['coef_lin'] = np.zeros(6,dtype=np.float32)


    str_to_latex = ""
    for i in range(6):
        f_with_gain.append(v_tr_with_gain(i).expand())
        # print(f_with_gain[-1])
        _t = sp.collect(f_with_gain[-1],e,evaluate=False)
        transimpedance_tension_vs_strain_coefficients_constant_gain['coef_ang'][i] = _t[e]
        transimpedance_tension_vs_strain_coefficients_constant_gain['coef_lin'][i] = _t[1]
        str_to_latex += "\\transimpedanceTension_{" + leg[i] + r"} & =\num{" + str(
            _t[e]*max_deformation) + "} \\cdot \\text{g}+\\num{"+str(_t[1])+"},\\unit{\\milli\\volt}, \\\\"
    print(str_to_latex)
    transimpedance_tension_vs_strain_coefficients_constant_gain.to_csv(
        'data/transimpedance_tension_vs_strain_coefficients_constant_gain.csv')


def calc_transimpedance_tension_vs_strain_coefficients_constant_gain_transmission():
    s_t = s_i(0.1957)

    print((4500./(230300*k_pd*(-lf.r_a[0] * l_laser + 1.0 - lf.r_l[0])/6.)+.393)/2.79)

    v_tr_with_gain = lambda ind: 230300 * k_pd * s_t * (
        -lf.r_a[ind] * lf.reflectivity_peak_wavelength[ind] * (1.0 - pe) * e +
        (-lf.r_a[ind] * l_laser + 1.0-lf.r_l[ind])) / 6.0
    f_with_gain = []

    transimpedance_tension_vs_strain_coefficients_constant_gain = pd.DataFrame(
    )
    transimpedance_tension_vs_strain_coefficients_constant_gain[
        'coef_ang'] = np.zeros(6, dtype=np.float32)
    transimpedance_tension_vs_strain_coefficients_constant_gain[
        'coef_lin'] = np.zeros(6, dtype=np.float32)

    str_to_latex = ""
    for i in range(6):
        f_with_gain.append(v_tr_with_gain(i).expand())
        # print(f_with_gain[-1])
        _t = sp.collect(f_with_gain[-1], e, evaluate=False)
        transimpedance_tension_vs_strain_coefficients_constant_gain[
            'coef_ang'][i] = _t[e]
        transimpedance_tension_vs_strain_coefficients_constant_gain[
            'coef_lin'][i] = _t[1]
        str_to_latex += "\\transimpedanceTension_{" + leg[
            i] + r"} & =\num{" + str(
                _t[e] * max_deformation) + "} \\cdot \\text{g}+\\num{" + str(
                    _t[1]) + "},\\unit{\\milli\\volt},  \\\\"
    print(str_to_latex)
    transimpedance_tension_vs_strain_coefficients_constant_gain.to_csv(
        'data/transimpedance_tension_vs_strain_coefficients_constant_gain_transmission.csv')
