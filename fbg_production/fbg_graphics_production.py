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
from scipy.signal import butter,lfilter
import locale
import h5py
import matplotlib.pyplot as plt
import pandas as pd
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29

def plot_production_one_column(date:str) -> None:
    f = h5py.File('production_files.hdf5','r')
    ff = f['fbg_production/'+date]
    # find index of fbgs
    fbg_keys = list(ff.keys())
    number_of_graphics = len(fbg_keys)
    print(number_of_graphics)
    if number_of_graphics >5:
        row_graphic = 5
        col_graphics = number_of_graphics // 5
    else:
        row_graphic = number_of_graphics
        col_graphics = 1
    fig, ax = plt.subplots(row_graphic, col_graphics, num=1 ,sharex=True, figsize=(FIG_L, 11))
    i = 0
    for ax_count in ax.flat:
        index_min, index_max = find_index_of_x_span(1540e-9,1560e-9,ff[fbg_keys[i]+'/wavelength_m'][:])
        ax_count.set_title(fbg_keys[i])
        ax_count.plot(ff[fbg_keys[i]+'/wavelength_m'][:] *
                   1e9, ff[fbg_keys[i]+'/reflectivity'][:,-1])
        ax_count.set_xlim(ff[fbg_keys[i]+'/wavelength_m'][index_min]*
                       1e9, ff[fbg_keys[i]+'/wavelength_m'][index_max] *
                       1e9)
        ax_count.set_ylim(0,1)
        i+=1
    
    fig.supylabel('Refletividade')

    _date = pd.to_datetime(date).strftime('%d de %B de %Y')
    fig.suptitle("FBG produzidas em "+_date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$") 
    f.close()
    plt.savefig("fbg_production/"+date+".pdf", format="pdf")
    plt.close(fig=1)


# plot_production_one_column('20231010')


def loop_graphics():
    f = h5py.File('production_files.hdf5','r') 
    ff = f['fbg_production/']
    fbg_keys = list(ff.keys())
    for i in fbg_keys:
        plot_production_one_column(date=i)
    f.close()
    




def plot_one_graphic(date:str):
    f = h5py.File('production_files.hdf5','r')
    ff = f['fbg_production/'+date]
    # find index of fbgs
    fbg_keys = list(ff.keys())
    number_of_graphics = len(fbg_keys)
    fig, ax = plt.subplots(1, 1, num=1 ,sharex=True, figsize=(FIG_L, 11))
    i = 0
    for i in fbg_keys:
        index_min, index_max = find_index_of_x_span(1540e-9,1560e-9,ff[i+'/wavelength_m'][:])
        # ax_count.set_title(i)
        ax.plot(ff[i+'/wavelength_m'][:] *
                      1e9, ff[i+'/reflectivity'][:, -1], label=i)
        ax.set_xlim(ff[i+'/wavelength_m'][index_min]*
                       1e9, ff[i+'/wavelength_m'][index_max] *
                       1e9)
        ax.set_ylim(0,1)    
    fig.supylabel('Refletividade')
    _date = pd.to_datetime(date).strftime('%d de %B de %Y')
    fig.suptitle("FBG produzidas em "+_date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$") 
    ax.legend(ncols=2)
    f.close()
    plt.savefig("fbg_production/"+date+".pdf", format="pdf")
    plt.close(fig=1)
    
    
def plot_one_graphic_filtered(date: str):
    f = h5py.File('production_files.hdf5', 'r')
    ff = f['fbg_production/'+date]
    # find index of fbgs
    fbg_keys = list(ff.keys())
    number_of_graphics = len(fbg_keys)
    fig, ax = plt.subplots(1, 1, num=1, sharex=True, figsize=(FIG_L, 11))
    i = 0
    for i in fbg_keys:
        index_min, index_max = find_index_of_x_span(
            1540e-9, 1560e-9, ff[i+'/wavelength_m'][:])
        # ax_count.set_title(i)
        ## Filtering
        # filtered_date = butter(5,2.0,f[i+'/reflectivity'][:, -1]
        
        ax.plot(ff[i+'/wavelength_m'][:] *
                1e9, ff[i+'/reflectivity'][:, -1], label=i)
        ax.set_xlim(ff[i+'/wavelength_m'][index_min] *
                    1e9, ff[i+'/wavelength_m'][index_max] *
                    1e9)
        ax.set_ylim(0, 1)
    fig.supylabel('Refletividade')
    _date = pd.to_datetime(date).strftime('%d de %B de %Y')
    fig.suptitle("FBG produzidas em "+_date)
    fig.supxlabel(r"$\lambda[\si{\nm}]$")
    ax.legend(ncols=2)
    f.close()
    plt.savefig("fbg_production/"+date+".pdf", format="pdf")
    plt.close(fig=1)
 
plot_one_graphic('20231031')
# plot_production_one_column('20231030')


# f = h5py.File('production_files.hdf5','a')
# ff = f['fbg_production/20231031']

# for i in list(ff.keys()):
#     # print(i,ff[i+'/optical_power'].size)
#     ff[i].create_dataset('optical_power_dbm', data=ff[i+'/optical_power'][:])
#     ff[i].create_dataset('wavelength_m', data=ff[i+'/wavelength'][:])
#     del ff[i+'/optical_power']
#     del ff[i+'/wavelength']
# f.close()





# f = h5py.File('production_files.hdf5','a')
# ff = f['fbg_production/20231031/fbg8']

# t = ff['optical_power_dbm'][:, :10]
# del ff['optical_power_dbm']
# ff['optical_power_dbm']=t
# t = ff['reflectivity'][:, :10]
# del ff['reflectivity']
# ff['reflectivity']=t

# f.close()


# f = h5py.File('production_files.hdf5','a')
# ff = f['fbg_production/20231030/fbg2']

# for i in range(4):
#     ff['reflectivity'][:,i] = reflectivity_transmition(ff['optical_power_dbm'][:,0],ff['optical_power_dbm'][:,i])
# f.close()