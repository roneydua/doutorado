#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   Untitled-1
@Time    :   2023/07/08 17:15:59
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import locale

import matplotlib.pyplot as plt
import numpy as np

from common_functions.generic_functions import mW_dbm

# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
FIG_L = 6.29
FIG_A = (90.0) / 25.4


R0_REF = 0.04
TRUE_R = 0.85


def refletividade_case_1(error=0.0, diff_ref=3.0):
    return R0_REF / (10.0**(0.1 * (diff_ref + error)) - 1.0)


def refletividade_case_2(error=0.0, diff_ref=3.0):
    return R0_REF * (10.0**(0.1 * (diff_ref + error)))


def refletividade_case_3(error=0.0, diff_ref=3.0):
    return 1.0 - (10**(0.1 * (diff_ref + error)))


def plot_error_comparision():
    error_dbm = [-.25,-0.15, 0.15,0.25]
    reflectivity_span = np.arange(0.1, 1.0, 0.01)
    # fig.clear()
    def make_fig():
        fig, ax = plt.subplots(
            1, 1, num=1, sharex=True, sharey=True, figsize=(
                FIG_L, 0.75*FIG_A))
        ax.set_ylim(0, 1)
        ax.set_ylabel('Refletividade estimada')
        ax.set_xlabel('Refletividade real')
        return fig, ax

    def plot(ax, _y, **kwargs):
        ax.plot(
            reflectivity_span,
            _y, label=r'$\varrho$=' + locale.format_string('%.2f', kwargs['error']) +
            r"\unit{\dbm}")
        ax.legend(ncols=5)

    fig, ax = make_fig()
    for i in error_dbm:
        _y = refletividade_case_1(
            error=i, diff_ref=mW_dbm(reflectivity_span + R0_REF) - mW_dbm(reflectivity_span))
        plot(ax=ax, _y=_y, error=i)
    _y = refletividade_case_1(
        diff_ref=mW_dbm(
            reflectivity_span +
            R0_REF) -
        mW_dbm(reflectivity_span))
    ax.plot(reflectivity_span, _y, ls=':', label="Exato")
    ax.legend(ncols=3)
    plt.savefig("../tese/images/method_error_1.pdf", format="pdf")
    plt.close(fig=1)
    fig,ax = make_fig()
    for i in error_dbm:
        _y = refletividade_case_2(error=i, diff_ref=mW_dbm(reflectivity_span) -
                                  mW_dbm(R0_REF))
        plot(ax=ax, _y=_y, error=i)
    _y = refletividade_case_2(
        diff_ref=mW_dbm(reflectivity_span) -
        mW_dbm(R0_REF))
    ax.plot(reflectivity_span, _y, ls=':', label="Exato")
    ax.legend(ncols=3)
    plt.savefig("../tese/images/method_error_2.pdf", format="pdf")
    plt.close(fig=1)
    fig, ax = make_fig()
    for i in error_dbm:
        _y = refletividade_case_3(
            error=i, diff_ref=mW_dbm(
                1.0 - reflectivity_span) - mW_dbm(1.0))
        plot(ax=ax, _y=_y, error=i)
    _y = refletividade_case_3(
        diff_ref=mW_dbm(
            1.0 -
            reflectivity_span) -
        mW_dbm(1.0))
    ax.legend(ncols=3)

    ax.plot(reflectivity_span, _y, ls=':', label="Exato")
    plt.savefig("../tese/images/method_error_3.pdf", format="pdf")
    plt.close(fig=1)

