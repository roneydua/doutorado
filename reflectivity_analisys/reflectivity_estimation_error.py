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
import sympy as sp
from IPython.core.interactiveshell import InteractiveShell
from ipywidgets import fixed, interactive

from common_functions.common_functions import *

# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# plt.style.use("default")
# plt.style.use("~/Dropbox/pacotesPython/roney3.mplstyle")
figL = 6.29
figA = (90.0) / 25.4


r_ref = 0.04
r = 0.85


def refletividade_case_1(e=0.0, diff_ref=3.0):
    return r_ref / (10.0**(0.1 * (e + diff_ref)) - 1.0)


def refletividade_case_2(e=0.0, diff_ref=3.0):
    return r_ref * (10.0**(0.1 * (e + diff_ref)))


def refletividade_case_3(e=0.0, diff_ref=3.0):
    return 1.0 - 10**(0.1 * (e + diff_ref))


error_dbm = np.range()

refletividade_case_1(e=-1, diff_ref=mW_dbm(r + r_ref) - mW_dbm(r))
refletividade_case_2(e=-1, diff_ref=mW_dbm(r) - mW_dbm(r_ref))
refletividade_case_3(e=-1, diff_ref=mW_dbm(1.0 - r) - mW_dbm(1.0))
