#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   inveser_problem_solution.py
@Time    :   2023/10/19 15:15:52
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''
import numpy as np
import locale
import h5py
import matplotlib.pyplot as plt
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4

from modeling.mathModelAccel import InverseProblem
f = h5py.File('teste.hdf5','w')