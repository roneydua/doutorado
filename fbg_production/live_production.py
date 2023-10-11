#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   live_production.py
@Time    :   2023/10/04 19:58:32
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import locale
import matplotlib.pyplot as plt
from time import sleep
from aquisitions.Q8347 import Q8347
from common_functions.generic_functions import reflectivity_transmition
from datetime import datetime
from pathlib import Path
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4


# t = np.arange(1,20,.1)
# y = np.sin(np.random.rand(1)[0]*t)
osa = Q8347(center=1550,span=20)
osa.read()


wavelength = osa.wavelength_nm
y = osa.optical_power_dbm
y_all = np.zeros((y.size,1))

y_all[:,0] = y
fig, ax = plt.subplots(2, 1, num=1, sharex=True,figsize=(FIG_L, FIG_A), dpi=144)
plt.show()
fig.supxlabel(r'$\lambda, [\si{\nm}]$')
ax[0].set_ylabel('Potência óptica [dBm]')
ax[1].set_ylabel('Refletividade')
iteration = 1
while True:
    try:
        ax[0].clear()
        ax[0].plot(wavelength*1e9, y)
        # ax[1].plot(wavelength*1e9, y_all[:, 0])
        r = reflectivity_transmition(d0=y_all[:, 0],di=y_all[:,iteration-1])
        ax[1].plot(wavelength*1e9, r)
        y_all = np.column_stack((y_all, y))
        plt.pause(.1)
        print('Aquisition...')
        osa.read()
        y = osa.optical_power_dbm
    except KeyboardInterrupt:
        print('Stop of aquisiton')
        break


data = np.column_stack((wavelength, y_all))
now = datetime.now()
# save_name = 'temp_data/'+now.strftime(r"%y%m%d/%H%M%S")+'.csv'
save_name = 'temp_data/'+now.strftime(r"%Y%m%d/%H%M%S")+'.csv'
print('save data on' +save_name )
np.savetxt(save_name, data)

