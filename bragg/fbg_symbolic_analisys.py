#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   fbg_symbolic_analisys.py
@Time    :   2023/05/13 18:08:24
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import locale
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from IPython.core.interactiveshell import InteractiveShell
from ipywidgets import fixed, interactive




# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("../../../../programasComuns/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
figL = 6.29
figA = (90.0) / 25.4

n_source = sp.symbols('r_{s}')
"Noise of source "

