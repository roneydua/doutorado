#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   transimpedance_static_analysis.py
@Time    :   2023/10/04 13:33:54
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import sympy as sp
import pandas as pd
import locale
import matplotlib.pyplot as plt
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4


folder = "../data/tempData/20231004/"
data = pd.read_csv(folder+'20231004.lvm',names=['time','fbg','tap','traction'])


fig, ax = plt.subplots(3, 1, num=1 ,sharex=True, figsize=(FIG_L, FIG_A))    
ax[0].plot(data.time/60.,data.traction,label='Tração')
ax[1].plot(data.time/60.,data.fbg,label='FBG')
ax[2].plot(data.time/60.0,data.tap,label='Tap')
fig.supxlabel(r"Tempo, \si{\minute}")
ax[0].set_ylim(-.72,-.67)

