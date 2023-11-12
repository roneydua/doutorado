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

from modeling.math_model_accel import InverseProblem, SimpleSolution

locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
plt.rcParams["figure.dpi"] = 100
FIG_L = 6.29
FIG_A = (90.0) / 25.4

test_name = "translational_movement"
# test_name = "complete_movement"
h5py_file_name = "modeling_data.hdf5"


graphics = [
            # 'differential_aligned',
            # 'differential_cross',
            'ls_translational',
            'ls_translational_angular_full',
            'ls_translational_angular_reduced',
            # 'one_fiber'
            ]


def recover_acceleration():
    # accel = AccelModelInertialFrame()
    fibers_with_length_info = np.array([1, 5, 9, 11])
    deformation = np.zeros((fibers_with_length_info.size, 1))
    ip_trans = InverseProblem(fibers_with_info=np.array([1, 5, 9, 11]))
    ip_trans_ang_full = InverseProblem(
        fibers_with_info=np.arange(1, 13), recover_angular_accel=True, estimation="full"
    )
    ip_trans_ang_reduced = InverseProblem(
        fibers_with_info=np.array([1, 4, 5, 8, 9, 11, 12]),
        recover_angular_accel=True,
        estimation="reduced",
    )
    ss = SimpleSolution(np.array([1, 5, 9]))
    f = h5py.File(h5py_file_name, "a")
    _ff = f[test_name]
    if "accel_recover" in _ff.keys():
        del _ff["accel_recover"]
    ff = _ff.require_group("accel_recover")
    # Simple methods
    ff.create_dataset("one_fiber", (3, _ff["t"].size))
    ff["one_fiber"].attrs["method_name"] = "Simples (1 FO)"
    ff.create_dataset("differential_aligned", (3, _ff["t"].size))
    ff["differential_aligned"].attrs["method_name"] = "Diferencial alinhado"
    ff.create_dataset("differential_cross", (3, _ff["t"].size))
    ff["differential_cross"].attrs["method_name"] = "Diferencial cruzado"
    # Least squared methods
    ff.create_dataset("ls_translational", (3, _ff["t"].size))
    ff["ls_translational"].attrs["method_name"] = "MMQ translacional"
    ff["ls_translational"].attrs["fibers_with_info"] = ip_trans.fibers_with_info
    ff.create_dataset("ls_translational_angular_reduced", (6, _ff["t"].size))
    ff["ls_translational_angular_reduced"].attrs["method_name"] = "MMQ reduzido"
    ff["ls_translational_angular_reduced"].attrs[
        "fibers_with_info"
    ] = ip_trans_ang_reduced.fibers_with_info
    ff.create_dataset("ls_translational_angular_full", (6, _ff["t"].size))
    ff["ls_translational_angular_full"].attrs["method_name"] = "MMQ completo"
    ff["ls_translational_angular_full"].attrs[
        "fibers_with_info"
    ] = ip_trans_ang_full.fibers_with_info
    for idx in tqdm(range(_ff["t"].size)):
        ff["one_fiber"][:, idx] = ss.estimated_ddrm_B(
            _ff["fiber_len"][:, idx], method="one_fiber"
        )
        ff["differential_aligned"][:, idx] = ss.estimated_ddrm_B(
            _ff["fiber_len"][:, idx], method="differential_aligned"
        )
        ff["differential_cross"][:, idx] = ss.estimated_ddrm_B(
            _ff["fiber_len"][:, idx], method="differential_cross"
        )
        # Least squared matrix
        ff["ls_translational"][:, idx] = ip_trans.compute_inverse_problem_solution(
            np.take(_ff["fiber_len"][:, idx], ip_trans.fibers_with_info_index)
        )
        (
            ff["ls_translational_angular_reduced"][:3, idx],
            ff["ls_translational_angular_reduced"][3:, idx],
        ) = ip_trans_ang_reduced.compute_inverse_problem_solution(
            np.take(
                _ff["fiber_len"][:, idx], ip_trans_ang_reduced.fibers_with_info_index
            )
        )
        (
            ff["ls_translational_angular_full"][:3, idx],
            ff["ls_translational_angular_full"][3:, idx],
        ) = ip_trans_ang_full.compute_inverse_problem_solution(
            np.take(_ff["fiber_len"][:, idx], ip_trans_ang_full.fibers_with_info_index)
        )

    f.close()


def print_recover_acceleration():
    f = h5py.File(h5py_file_name, "r")
    ff = f[test_name]
    fff = ff["accel_recover"]
    plt.close(1)
    fig, ax = plt.subplots(3, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))

    for i in range(3):
        ax[i].plot(ff["t"][:] * 1e3, ff["true_accel_i"][i, :], "--", label="Referência")
        lw = 4
        for j in fff.keys():
            if j in graphics:
                ax[i].plot(
                    ff["t"][:] * 1e3, fff[j][i, :], label=fff[j].attrs["method_name"],
                lw=lw)
                lw-=.5
    ax[2].legend()
    fig.supylabel(r"Aceleração $[\unit{\meter\per\second\squared}]$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    # fig.legend()
    plt.show()
    f.close()


def print_defromations():
    f = h5py.File(h5py_file_name, "r")
    ff = f[test_name]
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
    f = h5py.File(h5py_file_name, "r")
    ff = f[test_name]
    plt.close(3)
    fig, ax = plt.subplots(4, 3, num=3, sharex=True, figsize=(FIG_L, FIG_A))
    ind_x = [0, 6, 13, 20]
    ax[0][0].set_ylabel(r"$\dot{\mathbf{r}},\si{\meter\per\second}$")
    ax[1][0].set_ylabel(r"$\mathbf{r},\si{\meter}$")
    ax[2][0].set_ylabel(r"$\mathbf{q}$")
    ax[3][0].set_ylabel(r"$\boldsymbol{\omega},\si{\radian\per\second}$")
    for lin in range(4):
        for col in range(3):
            ax[lin, col].plot(ff["t"][:], ff["x"][ind_x[lin] + col, :])
            ax[lin, col].plot(ff["t"][:], ff["x"][ind_x[lin] + col + 3, :])

    # fig.supylabel(r'Aceleração $[\unit{\meter\per\second\squared}]$')
    fig.supxlabel(r"Tempo $[\unit{\second}]$")
    # fig.legend()
    plt.show()
    f.close()


if __name__ == "__main__":
    recover_acceleration()
