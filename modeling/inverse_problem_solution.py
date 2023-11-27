#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   inverse_problem_solution.py
@Time    :   2023/10/19 15:15:52
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""
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

# TEST_NAME = "translational_movement"
TEST_NAME = "complete_movement"
# TEST_NAME = "angular_movement"
H5PY_FILE_NAME = "modeling_data.hdf5"


# ONLY_SHOW = True
ONLY_SHOW = False

# if ONLY_SHOW:
# plt.rcParams["figure.dpi"] = 100

graphics = [
    # 'differential_aligned',
    "differential_cross",
    # "ls_translational",
    "ls_translational_angular_full",
    # "ls_translational_angular_reduced",
    # "one_fiber",
]


def recover_acceleration():
    # accel = AccelModelInertialFrame()
    fibers_with_length_info = np.array([1, 5, 9, 11])
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
    f = h5py.File(H5PY_FILE_NAME, "a")
    _ff = f[TEST_NAME]
    if "accel_recover" in _ff.keys():
        del _ff["accel_recover"]
    ff = _ff.require_group("accel_recover")
    if "ls_solution" in _ff.keys():
        del _ff["ls_solution"]
    fff = _ff.require_group("ls_solution")
    # _delta_time_size_new = int(1e-2 / (_ff["t"][1] - _ff["t"][0]))
    _delta_n = int(1e-3 / (_ff["t"][1] - _ff["t"][0]))
    print(_delta_n)
    _time_size = int(_ff["t"].size * (_ff["t"][1] - _ff["t"][0]) / 1e-3)
    if _time_size < 1:
        print("_time_size muito pequeno")
        _delta_n = 1
        _time_size = _ff["t"].size
    ff.create_dataset("t", _time_size, dtype=np.float64)
    # Simple methods
    ff.create_dataset("one_fiber", (3, _time_size), dtype=np.float64)
    ff["one_fiber"].attrs["method_name"] = "Simples (1 FO)"
    ff.create_dataset("differential_aligned", (3, _time_size), dtype=np.float64)
    ff["differential_aligned"].attrs["method_name"] = "Diferencial alinhado"
    ff.create_dataset("differential_cross", (3, _time_size), dtype=np.float64)
    ff["differential_cross"].attrs["method_name"] = "Diferencial cruzado"
    # Least squared methods
    ff.create_dataset("ls_translational", (3, _time_size), dtype=np.float64)
    ff["ls_translational"].attrs["method_name"] = "MMQ translacional"
    ff["ls_translational"].attrs["fibers_with_info"] = ip_trans.fibers_with_info
    fff.create_dataset("ls_translational", (4, _time_size), dtype=np.float64)
    fff["ls_translational"].attrs["method_name"] = "MMQ translacional"

    ff.create_dataset(
        "ls_translational_angular_reduced", (6, _time_size), dtype=np.float64
    )
    ff["ls_translational_angular_reduced"].attrs["method_name"] = "MMQ reduzido"
    ff["ls_translational_angular_reduced"].attrs[
        "fibers_with_info"
    ] = ip_trans_ang_reduced.fibers_with_info
    fff.create_dataset(
        "ls_translational_angular_reduced", (7, _time_size), dtype=np.float64
    )
    fff["ls_translational_angular_reduced"].attrs["method_name"] = "MMQ reduzido"

    ff.create_dataset(
        "ls_translational_angular_full", (6, _time_size), dtype=np.float64
    )
    ff["ls_translational_angular_full"].attrs["method_name"] = "MMQ completo"
    ff["ls_translational_angular_full"].attrs[
        "fibers_with_info"
    ] = ip_trans_ang_full.fibers_with_info
    fff.create_dataset(
        "ls_translational_angular_full", (10, _time_size), dtype=np.float64
    )
    fff["ls_translational_angular_full"].attrs["method_name"] = "MMQ completo"
    idx_new = 0
    for idx in tqdm(range(0, int(_ff["t"].size) - 1, _delta_n)):
        ff["t"][idx_new] = _ff["t"][idx]

        ff["one_fiber"][:, idx_new] = ss.estimated_ddrm_B(
            _ff["fiber_len"][:, idx], method="one_fiber"
        )
        ff["differential_aligned"][:, idx_new] = ss.estimated_ddrm_B(
            _ff["fiber_len"][:, idx], method="differential_aligned"
        )
        ff["differential_cross"][:, idx_new] = ss.estimated_ddrm_B(
            _ff["fiber_len"][:, idx], method="differential_cross"
        )
        # Least squared matrix
        ff["ls_translational"][:, idx_new] = ip_trans.compute_inverse_problem_solution(
            np.take(_ff["fiber_len"][:, idx], ip_trans.fibers_with_info_index)
        )
        fff["ls_translational"][:, idx_new] = ip_trans.var_gamma

        (
            ff["ls_translational_angular_reduced"][:3, idx_new],
            ff["ls_translational_angular_reduced"][3:, idx_new],
        ) = ip_trans_ang_reduced.compute_inverse_problem_solution(
            np.take(
                _ff["fiber_len"][:, idx], ip_trans_ang_reduced.fibers_with_info_index
            )
        )
        fff["ls_translational_angular_reduced"][
            :, idx_new
        ] = ip_trans_ang_reduced.var_gamma
        (
            ff["ls_translational_angular_full"][:3, idx_new],
            ff["ls_translational_angular_full"][3:, idx_new],
        ) = ip_trans_ang_full.compute_inverse_problem_solution(
            np.take(_ff["fiber_len"][:, idx], ip_trans_ang_full.fibers_with_info_index)
        )
        fff["ls_translational_angular_full"][:, idx_new] = ip_trans_ang_full.var_gamma
        idx_new += 1
    f.close()


def plot_recover_acceleration():
    f = h5py.File(H5PY_FILE_NAME, "r")
    ff = f[TEST_NAME]
    fff = ff["accel_recover"]
    plt.close(1)
    fig, ax = plt.subplots(3, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))

    for i in range(3):
        ax[i].plot(ff["t"][:] * 1e3, ff["true_accel_b"][i, :], "--", label="Referência")
        # lw = 6
        for j in fff.keys():
            if j in graphics:
                ax[i].plot(
                    fff["t"][:] * 1e3, fff[j][i, :], label=fff[j].attrs["method_name"]
                )
    ax[2].legend()
    fig.supylabel(r"Aceleração $[\unit{\meter\per\second\squared}]$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    # fig.legend()
    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/recover_translational_acceleration_"
            + TEST_NAME
            + ".pdf",
            format="pdf",
        )
        plt.close(fig=3)
    else:
        plt.show()
    f.close()


def plot_recover_angular_acceleration():
    f = h5py.File(H5PY_FILE_NAME, "r")
    ff = f[TEST_NAME]
    fff = ff["accel_recover"]
    plt.close(1)
    y_axis_name = [r"$\dot{\omega}_{x}$", r"$\dot{\omega}_{y}$", r"$\dot{\omega}_{z}$"]
    fig, ax = plt.subplots(3, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))

    for i in range(3):
        ax[i].plot(
            ff["t"][:] * 1e3,
            ff["true_angular_acceleration_b"][i, :],
            "--",
            label="Referência",
        )
        ax[i].set_ylabel(y_axis_name[i])
        # lw = 4
        for j in ["ls_translational_angular_full", "ls_translational_angular_reduced"]:
            ax[i].plot(
                fff["t"][:] * 1e3,
                fff[j][i + 3, :],
                label=fff[j].attrs["method_name"],
                # lw=lw,
            )
            # lw /= 2
    ax[2].legend()
    fig.supylabel(r"Aceleração angular $[\unit{\radian\per\second\squared}]$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    # fig.legend()
    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/recover_angular_acceleration_" + TEST_NAME + ".pdf",
            format="pdf",
        )
        plt.close(fig=1)
    else:
        plt.show()
    f.close()


def plot_defromations():
    f = h5py.File(H5PY_FILE_NAME, "r")
    ff = f[TEST_NAME]
    plt.close(2)
    fig, ax = plt.subplots(3, 4, num=2, sharex=True, figsize=(FIG_L, FIG_A))

    for i, _ax in enumerate(ax.flat):
        _ax.plot(
            1e3 * ff["t"][:],
            1e6 * (ff["fiber_len"][i, :] / ff.attrs["fiber_length"] - 1.0),
        )
        _ax.set_title(r"$\text{FBG}_{" + str(i + 1) + "}$")
        # _ax.set_xlim(50,55)
        # _ax.set_ylim(-50, 50)

    fig.supylabel(r"Deformação $\frac{\ell-\ell_0}{\ell_0}\times{10}^{6}$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    # fig.legend()

    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/deformations_" + TEST_NAME + ".pdf", format="pdf"
        )
        plt.close(fig=3)
    else:
        plt.show()
    f.close()


def plot_inertial_states():
    f = h5py.File(H5PY_FILE_NAME, "r")
    ff = f[TEST_NAME]
    plt.close(3)
    fig, ax = plt.subplots(4, 3, num=3, sharex=True, figsize=(FIG_L, FIG_A * 1.25))
    ind_x = [0, 6, 12, 20]
    ax[0][0].set_ylabel(r"$\dot{\mathbf{r}}[\si{\meter\per\second}]$")
    ax[1][0].set_ylabel(r"$\mathbf{r}[\si{\meter}]$")
    ax[2][0].set_ylabel(r"$\mathbf{q}$")
    ax[3][0].set_ylabel(r"$\boldsymbol{\omega}[\si{\radian\per\second}]$")
    for lin in range(4):
        if lin == 2:
            # pass
            for col in range(3):
                ax[lin, col].plot(
                    1e3 * ff["t"][:],
                    np.float32(ff["x"][ind_x[lin] + col + 1, :]),
                    lw=3,
                    label=r"$q_{" + str(col + 1) + "}$ Base",
                )
                ax[lin, col].plot(
                    1e3 * ff["t"][:],
                    np.float32(ff["x"][ind_x[lin] + col + 5, :]),
                    label=r"$q_{" + str(col + 1) + "}$ Massa sísmica",
                )
            ax[lin, 0].plot(
                1e3 * ff["t"][:],
                np.float32(ff["x"][ind_x[lin], :]),
                label=r"$q_{0}$ Base",
            )
            ax[lin, 0].plot(
                1e3 * ff["t"][:],
                np.float32(ff["x"][ind_x[lin] + 4, :]),
                label=r"$q_{0}$ Massa sísmica",
            )
            ax[lin, 0].legend()
        else:
            for col in range(3):
                ax[lin, col].plot(
                    1e3 * ff["t"][:], np.float32(ff["x"][ind_x[lin] + col, :]), lw=3
                )
                ax[lin, col].plot(
                    1e3 * ff["t"][:], np.float32(ff["x"][ind_x[lin] + col + 3, :])
                )
    ax[0, 0].legend(
        [r"Base", "Massa sísmica"],
    )
    # fig.supylabel(r'Aceleração $[\unit{\meter\per\second\squared}]$')
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    # fig.legend()
    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/all_states_" + TEST_NAME + ".pdf", format="pdf"
        )
        plt.close(fig=3)
    else:
        plt.show()

    f.close()


def plot_length_of_fo():
    f = h5py.File(H5PY_FILE_NAME, "r")
    ff = f[TEST_NAME]
    plt.close(2)
    fig, ax = plt.subplots(3, 4, num=2, sharex=True, figsize=(FIG_L, FIG_A))

    for i, _ax in enumerate(ax.flat):
        _ax.plot(
            1e3 * ff["t"][:], 1e3* ff["fiber_len"][i, :],
        )
        _ax.set_title(r"$\text{FBG}_{" + str(i + 1) + "}$")
        # _ax.set_xlim(50,55)
        _ax.set_ylim(2, 4)

    fig.supylabel(r"$\ell$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    # fig.legend()

    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/length_of_fo_" + TEST_NAME + ".pdf", format="pdf"
        )
        plt.close(fig=2)
    else:
        plt.show()
    f.close()


def plot_ls_solution():
    f = h5py.File(H5PY_FILE_NAME, "r")
    ff = f[TEST_NAME]
    fff = ff["accel_recover"]
    plt.close(4)
    fig, ax = plt.subplots(3, 1, num=4, sharex=True, figsize=(FIG_L, FIG_A))

    for i in range(3):
        ax[i].plot(
            1e3 * ff["t"][:],
            1e6 * ff["true_relative_position"][i, :],
            label="Referência",
        )
        for j in ff["ls_solution"].keys():
            ax[i].plot(
                1e3 * fff["t"][:],
                1e6 * ff["ls_solution/" + j][i + 1, :],
                label=ff["ls_solution/" + j].attrs["method_name"],
            )

    fig.supylabel(r"Posição relativa ${\mathcal{B}}^{}\mathbf{r}_{m}$ $[\unit{\um}]$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    ax[0].legend()
    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/estimated_relative_position_" + TEST_NAME + ".pdf",
            format="pdf",
        )
        plt.close(fig=4)
    else:
        plt.show()

    plt.close(fig=5)
    fig, ax = plt.subplots(4, 1, num=5, sharex=True, figsize=(FIG_L, FIG_A))

    for i in range(4):
        ax[i].plot(
            1e3 * ff["t"][:],
            ff["true_relative_orientation"][i, :],
            label="Referência",
        )
        for j in ["ls_translational_angular_reduced", "ls_translational_angular_full"]:
            if i == 0:
                q0 = ff["ls_solution/" + j][-4, :]
                for k in range(fff["t"].size):
                    q0[k] = np.sqrt(
                        1.0
                        - ff["ls_solution/" + j][-3, k] ** 2
                        - ff["ls_solution/" + j][-2, k] ** 2
                        - ff["ls_solution/" + j][-1, k] ** 2
                    )
                ax[i].plot(
                    1e3 * fff["t"][:],
                    q0,
                    label=ff["ls_solution/" + j].attrs["method_name"],
                )
            else:
                ax[i].plot(
                    1e3 * fff["t"][:],
                    ff["ls_solution/" + j][-4 + i, :],
                    label=ff["ls_solution/" + j].attrs["method_name"],
                )

    fig.supylabel(r"Orientação relativa $\,_{\mathcal{M}}^{\mathcal{B}}\mathbf{q}$")
    fig.supxlabel(r"Tempo $[\unit{\ms}]$")
    ax[0].legend()
    if not ONLY_SHOW:
        plt.savefig(
            "../dissertacao/images/estimated_relative_orientation_"
            + TEST_NAME
            + ".pdf",
            format="pdf",
        )
        plt.close(fig=5)
    else:
        plt.show()
    f.close()


if __name__ == "__main__":
    recover_acceleration()
    plot_inertial_states()
    plot_recover_acceleration()
    plot_recover_angular_acceleration()
    plot_ls_solution()
    plot_defromations()
