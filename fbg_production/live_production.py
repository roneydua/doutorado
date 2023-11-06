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
import h5py
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4

osa = Q8347(center=1550,span=20,high_resolution=True)
osa.read()


wavelength_m = osa.wavelength_m
y = osa.optical_power_dbm
y_all = np.zeros((y.size,1))
reflectivity = np.zeros((y.size,1))


y_all[:,0] = y
fig, ax = plt.subplots(2, 1, num=1, sharex=True,figsize=(FIG_L, 2*FIG_A), dpi=100)
plt.show()
fig.supxlabel(r'$\lambda, [\si{\nm}]$')
iteration=1
while True:
    try:
        print('Aquisition number '+str(iteration))
        r = reflectivity_transmition(d0=y_all[:, 0],di=y_all[:,iteration-1])
        ax[0].plot(wavelength_m*1e9, y)
        ax[1].clear()
        ax[1].set_ylabel('Refletividade')
        ax[1].plot(wavelength_m*1e9, r)
        ax[1].set_ylim(-.10,1)
        plt.pause(.01)
        osa.read()
        y = osa.optical_power_dbm
        y_all = np.column_stack((y_all, y))
        reflectivity = np.column_stack((reflectivity, r))
        iteration +=1
    except KeyboardInterrupt:
        print('Stop of aquisiton')
        osa.close()
        break

print("save data on hdf file")
f = h5py.File('fbg_production/production_file_temp.hdf5','a')
now = datetime.now()
ff = f.require_group(now.strftime(r"%Y%m%d"))
# check number of fbg on group
number_of_dataset = len(ff.keys())
fff = ff.require_group('fbg'+str(1+number_of_dataset))
fff.attrs['optical_fiber'] = 'sm1500(4.2/125)'
fff.create_dataset('wavelength_m', data=wavelength_m)
fff.create_dataset('optical_power_dbm',data=y_all)
fff.create_dataset('reflectivity',data=reflectivity)
f.close()
plt.savefig('fbg_production'+'/fig_production_figs/'+now.strftime(r"%Y%m%d") +
            '/fbg'+str(1+number_of_dataset)+".png", format="png", transparent=False)
plt.close(fig=1)
plt.close(fig=1)
print("ok")


