#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   push_pull_interrogation_analysis.py
@Time    :   2024/01/25 14:42:05
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""

import numpy as np
import sympy as sp
from sympy.integrals.risch import risch_integrate
import locale
import matplotlib.pyplot as plt
import h5py

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
from modeling.math_model_accel import AccelModelInertialFrame
from common_functions.generic_functions import calc_laser

FIG_L = 6.29
FIG_A = (90.0) / 25.4
from scipy.ndimage import shift


class fbg_simulation(object):
    def __init__(self):
        f = h5py.File("./production_files_.hdf5", "r")
        w_fbg_e_ = f["fbg_production/20231130/fbg4/wavelength_m"][:]
        "wavelength in meters"
        fbg_e_ = f["fbg_production/20231130/fbg4/reflectivity"][:, -1]
        w_fbg_d_ = f["fbg_production/20231130/fbg9/wavelength_m"][:]
        "wavelength in meters"
        fbg_d_ = f["fbg_production/20231130/fbg9/reflectivity"][:, -1]
        f.close()

        self.w_fbg_d_interp, self.w_fbg_d = self.extend_vector(w_fbg_d_, "wavelength")
        "wavelength in meters"
        self.fbg_d = self.extend_vector(fbg_d_, "reflectivity")
        self.fbg_d_interp = np.interp(self.w_fbg_d_interp, self.w_fbg_d, self.fbg_d)

        self.w_fbg_e_interp, self.w_fbg_e = self.extend_vector(w_fbg_e_, "wavelength")
        "wavelength in meters"
        self.fbg_e = self.extend_vector(fbg_e_, "reflectivity")
        self.fbg_e_interp = np.interp(self.w_fbg_e_interp, self.w_fbg_e, self.fbg_e)
        self.power_of_second_reflected_spectrum_mW = 0.0
        self.step_of_w_fbg = np.diff(self.w_fbg_d_interp).mean()
        "wavelength in meters"

        self.amplitude_mw_by_nm = 20 / 20  # 20mW distributed in  40nm
        self.translate_fbgs(0)

    def extend_vector(self, v: np.ndarray, type_of_data: str):
        if type_of_data == "wavelength":
            _w = np.hstack(
                (v - (v.size * np.diff(v).mean()), v, v + (v.size * np.diff(v).mean()))
            )
            return np.linspace(_w[0], _w[-1], 100_000), _w

        else:
            return np.hstack((np.zeros(v.size), v, np.zeros(v.size)))

    def translate_fbgs(self, delta_lambda_d: float, delta_lambda_e=None):
        index_to_shift_d = int(delta_lambda_d // self.step_of_w_fbg)
        if delta_lambda_e == None:
            index_to_shift_e = -index_to_shift_d
        else:
            index_to_shift_e = int(delta_lambda_e // self.step_of_w_fbg)

        self.fbg_d_shifted = shift(
            self.fbg_d_interp, index_to_shift_d, mode="grid-wrap"
        )
        self.fbg_e_shifted = shift(
            self.fbg_e_interp, index_to_shift_e, mode="grid-wrap"
        )
        self.make_dot_product()

    def make_dot_product(self):
        self.first_reflected_spectrum = self.amplitude_mw_by_nm * self.fbg_e_shifted

        self.second_reflected_spectrum = (
            self.first_reflected_spectrum * self.fbg_d_shifted
        )
        self.power_of_second_reflected_spectrum_mW = (
            np.trapz(y=self.second_reflected_spectrum) * self.step_of_w_fbg * 1e9
        )


def power_vs_delta_lambda_animation():
    fbgs = fbg_simulation()
    alpha = np.linspace(1, 0, 10)
    fig, ax = plt.subplots(
        2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A * 1.25), dpi=72
    )

    index = 0
    for i in np.arange(-1, 1, 0.2):
        #
        # ax.plot(fbgs.w_fbg_d * 1e9, fbgs.fbg_d_shifted)
        # ax.plot(fbgs.w_fbg_e * 1e9, fbgs.fbg_e_shifted)
        fig.clear()
        fig, ax = plt.subplots(
            2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A * 1.25), dpi=72
        )
        fbgs.translate_fbgs(i * 1e-9)
        ax[0].plot(
            fbgs.w_fbg_e_interp * 1e9,
            fbgs.fbg_e_shifted,
            lw=0.5,
            ls="-",
            color=my_colors[0],
            alpha=alpha[index],
        )
        ax[0].plot(
            fbgs.w_fbg_d_interp * 1e9,
            fbgs.fbg_d_shifted,
            lw=0.5,
            ls="-",
            color=my_colors[1],
            alpha=alpha[index],
        )
        ax[0].set_ylabel("Refletividade")
        ax[1].plot(
            fbgs.w_fbg_d_interp * 1e9,
            fbgs.second_reflected_spectrum,
            lw=0.5,
            ls="-",
            color=my_colors[2],
            label=locale.format_string(
                r"%.3f \unit{\nano\watt}",
                fbgs.power_of_second_reflected_spectrum_mW * 1e9,
            ),
            alpha=alpha[index],
        )
        # index += 1e
        # ax[1].legend()
        ax[1].set_ylabel(r"Potência óptica [\unit{\milli\watt\per\nano\meter}]")
        ax[1].set_xlabel(r"Comprimento de onda [\unit{\nano\meter}]")
        ax[0].set_xlim(left=1540, right=1560)
        plt.pause(0.01)

    # plt.savefig(
    #     "../tese/images/power_vs_delta_lambda_animation.pdf", format="pdf"
    # )
    # plt.close(fig=1)


def plot_power_vs_accel_accel1():
    fbgs = fbg_simulation()
    accel = AccelModelInertialFrame()

    def compute_delta_lambda(dd_r: float, delta_T=0.0):
        lambda_e = 1548e-9
        lambda_d = 1552e-9
        pe = 0.23
        alpha = 0.55e-6
        zeta = 8.60e-6
        temp_var = (
            1.0 - pe
        ) * accel.seismic_mass * accel.fiber_length * dd_r / accel.E / np.pi / (
            accel.fiber_diameter
        ) ** 2 + (
            alpha + zeta
        ) * delta_T
        delta_l_e = lambda_e * temp_var
        delta_l_d = lambda_d * temp_var
        return delta_l_e, delta_l_d

    dd_r_vec = np.arange(-200, 200, 10)
    pot_vec = 0.0 * dd_r_vec
    delta_l_d = 0.0 * dd_r_vec
    delta_l_e = 0.0 * dd_r_vec

    for i in range(len(dd_r_vec)):
        delta_l_e[i], delta_l_d[i] = compute_delta_lambda(dd_r_vec[i], 0)
        # print(delta_l_e[i], delta_l_d[i])
        fbgs.translate_fbgs(delta_lambda_d=delta_l_d[i], delta_lambda_e=-delta_l_e[i])
        pot_vec[i] = fbgs.power_of_second_reflected_spectrum_mW
    fig, ax = plt.subplots(2, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))
    ax[0].plot(dd_r_vec, pot_vec, "-*")
    ax[1].plot(dd_r_vec, -delta_l_e * 1e12, "-*")
    ax[1].plot(dd_r_vec, delta_l_d * 1e12, "-*")
    ax[1].legend([r"$\Delta\lambda_{g,1}$", r"$\Delta\lambda_{g,2}$"])
    ax[0].set_ylabel(r"Potência óptica [\unit{\nano\watt}]")
    ax[1].set_ylabel("Variação do $\\lambda_{g} [\\unit{\\pico\\meter}]$")
    ax[1].set_xlabel(
        r"Aceleração ${}^\mathcal{B}\ddot{\mathbf{r}}$ [\unit{\meter\per\second\squared}]"
    )
    plt.savefig("../tese/images/plot_power_vs_accel_accel1.pdf", format="pdf")
    plt.close(fig=1)
    print(
        "Accel is null in",
        np.where(dd_r_vec == 0)[0][0],
        "with power in [nW] equal ",
        pot_vec[np.where(dd_r_vec == 0)[0][0]],
    )
    # plt.pause(0.01)


def simulation_push_pull_symbolic():
    def gaussian_func(wavelength, wavelength_center, standart_deviation):
        calc_a = 1 / (standart_deviation * sp.sqrt(2 * sp.pi))
        return calc_a * sp.exp(
            -(((wavelength - wavelength_center) / standart_deviation) ** 2) / 2
        )

    l = sp.symbols(r"\lambda")
    l_center_1 = sp.symbols(r"\lambda_{g\,1}", real=True)
    l_center_2 = sp.symbols(r"\lambda_{g\,2}", real=True)
    standart_deviation = sp.symbols(r"\sigma", real=True)
    f_1 = gaussian_func(l, l_center_1, standart_deviation)
    f_2 = gaussian_func(l, l_center_2, standart_deviation)
    f_1 * f_2
    # integral_f1_f2 = sp.integrate(f_1 * f_2, l)
    integral_f1_f2 = risch_integrate(f_1 * f_2, l).doit()
    print(sp.latex(integral_f1_f2))
    print(sp.latex(f_2))
