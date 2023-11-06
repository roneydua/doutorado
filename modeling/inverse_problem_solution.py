#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   inverse_problem_solution.py
@Time    :   2023/10/19 15:15:52
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""
from modeling.mathModelAccel import InverseProblem, SimpleSolution
import locale

import h5py
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
plt.rcParams["figure.dpi"] = 100
FIG_L = 6.29
FIG_A = (90.0) / 25.4


def recover_acceleration():
    # accel = AccelModelInertialFrame()
    fibers_with_length_info = np.array([1, 5, 9, 11])
    deformation = np.zeros((fibers_with_length_info.size, 1))
    ip = InverseProblem(fibers_with_length_info)
    ss = SimpleSolution(np.array([1, 5, 9]))
    f = h5py.File("modeling_data.hdf5", "a")
    ff = f["pure_translational_movement"]
    # ff.create_dataset('recover_accel_simple', (3, ff['t'].size))
    # ff.create_dataset('recover_accel_ls_simple', (3, ff['t'].size))
    # Solve inverse problem
    temp_fiber_l = np.take(ff["fiber_len"][:], fibers_with_length_info - 1, axis=0)

    for idx_t, _t in enumerate(tqdm(ff["t"])):
        # for idx_t, _t in enumerate(tqdm(range(10))):
        # ip.compute_inverse_problem_solution(temp_fiber_l[:, idx_t])
        # ip.estimate_f_vector()
        # ff['recover_accel_ls_simple'][:, idx_t] = ip.estimate_ddrm_B()
        ff["recover_accel_simple"][:, idx_t] = ss.estimated_ddrm_B(
            ff["fiber_len"][:, idx_t]
        )

    f.close()


# recover_acceleration()


def print_recover_acceleration():
    f = h5py.File("modeling_data_4.hdf5", "r")
    ff = f["pure_translational_movement"]
    plt.close(1)
    fig, ax = plt.subplots(3, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))

    for i in range(3):
        ax[i].plot(ff["t"][:], ff["true_accel"][i, :], label="Referência")
        ax[i].plot(ff["t"][:-1], ff["recover_accel_simple"][i, :-1], label="Simples")
        ax[i].plot(
            ff["t"][:],
            ff["recover_accel_ls_simple"][i, :],
            label="Cruzado de translação",
        )

    ax[2].legend()
    fig.supylabel(r"Aceleração $[\unit{\meter\per\second\squared}]$")
    fig.supxlabel(r"Tempo $[\unit{\second}]$")
    # fig.legend()
    plt.show()
    f.close()


def print_defromations():
    f = h5py.File("modeling_data_4.hdf5", "r")
    ff = f["pure_translational_movement"]
    plt.close(2)
    fig, ax = plt.subplots(3, 4, num=2, sharex=True, figsize=(FIG_L, FIG_A))

    for i, _ax in enumerate(ax.flat):
        print(i)
        _ax.plot(ff["t"][:], ff["fiber_len"][i, :])

    # fig.supylabel(r'Aceleração $[\unit{\meter\per\second\squared}]$')
    fig.supxlabel(r"Tempo $[\unit{\second}]$")
    # fig.legend()
    plt.show()
    f.close()


def print_inertial_states():
    f = h5py.File("modeling_data_4.hdf5", "r")
    ff = f["pure_translational_movement"]
    plt.close(3)
    fig, ax = plt.subplots(4, 3, num=3, sharex=True, figsize=(FIG_L, FIG_A))
    indx_x = [0, 6, 13, 20]
    ax[0][0].set_ylabel(r"$\dot{\mathbf{r}},\si{\meter\per\second}$")
    ax[1][0].set_ylabel(r"$\mathbf{r},\si{\meter}$")
    ax[2][0].set_ylabel(r"$\mathbf{q}$")
    ax[3][0].set_ylabel(r"$\boldsymbol{\omega},\si{\radian\per\second}$")
    for lin in range(4):
        for col in range(3):
            ax[lin, col].plot(ff["t"][:], ff["x"][indx_x[lin] + col, :])
            ax[lin, col].plot(ff["t"][:], ff["x"][indx_x[lin] + col + 3, :])

    # fig.supylabel(r'Aceleração $[\unit{\meter\per\second\squared}]$')
    fig.supxlabel(r"Tempo $[\unit{\second}]$")
    # fig.legend()
    plt.show()
    f.close()
