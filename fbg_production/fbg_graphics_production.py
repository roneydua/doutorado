#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   fbg_graphics_production.py
@Time    :   2023/10/02 10:08:38
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

from common_functions.generic_functions import find_index_of_x_span, reflectivity_transmition
import numpy as np
import locale
import h5py
import matplotlib.pyplot as plt
import pandas as pd
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29

def plot_production_one_collum(date:str):
    f = h5py.File('production_files.hdf5','r')
    ff = f['fbg_production/'+date]
    fbg_keys = list(ff.keys())
    number_of_graphics = len(fbg_keys)
    # find index of fbgs
    fig, ax = plt.subplots(number_of_graphics, 1, num=1 ,sharex=True, figsize=(FIG_L, 11))
    if number_of_graphics == 1:
        _ax =  ax
        ax = [_ax]
    for i in range(number_of_graphics):
        index_min, index_max = find_index_of_x_span(1540e-9,1560e-9,ff[fbg_keys[i]+'/wavelength_m'][:])
        ax[i].set_title(fbg_keys[i])
        ax[i].plot(ff[fbg_keys[i]+'/wavelength_m'][:] *
                   1e9, ff[fbg_keys[i]+'/reflectivity'][:,-1])
        ax[i].set_xlim(ff[fbg_keys[i]+'/wavelength_m'][index_min]*
                       1e9, ff[fbg_keys[i]+'/wavelength_m'][index_max] *
                       1e9)
        ax[i].set_ylim(0,1)
        ax[i].set_ylabel('Refletividade')

    _date = pd.to_datetime(date).strftime('%d de %B de %Y')
    fig.suptitle("FBG produzidas em "+_date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$") 
    f.close()
    plt.savefig("fbg_production/"+date+".pdf", format="pdf")
    plt.close(fig=1)




def loop_graphics():
    f = h5py.File('production_files.hdf5','r') 
    ff = f['fbg_production/']
    fbg_keys = list(ff.keys())
    for i in fbg_keys:
        plot_production_one_collum(date=i)
    f.close()
    



