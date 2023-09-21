#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   graphics_fbg_production.py
@Time    :   2023/08/16 11:11:55
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import locale
import os

import h5py
import lvm_read
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd
import sympy as sp
from mpl_toolkits.mplot3d import axes3d

# from IPython.core.interactiveshell import InteractiveShell
# from ipywidgets import interactive, fixed
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
from common_functions.generic_functions import reflectivity_transmition

FIG_L = 6.29
FIG_A = (90.0) / 25.4

data_folder = './../data/danteAlex/20230814/'


def fix_time(t):
    t_fixed = 1.0*t
    count = 0
    for i in range(2, len(t)):
        if t[i]-t[i-1] < 0:
            count += 60
        t_fixed[i] = count+t[i]
    return t_fixed





def plot_all_fbgs():
    '''
    plot_all_fbgs Graphic Plot before and after after_dehydrogenation
    '''
    # fig.clear()
    # fig2.clear()
    fig, ax = plt.subplots(1, 1, num=2, sharex=True, figsize=(FIG_L, FIG_A))
    # fig, ax = plt.subplots(num=1,figsize=(FIG_L, FIG_A),subplot_kw={"projection":"3d"})
    f = h5py.File("./../data/phd_data.hdf5",'r')
    r_before = f['fbg_production/test1/fbg5/r/315'][...]
    wavelength_before = f['fbg_production/test1/fbg5/wavelength/315'][...]
    r_after = f['fbg_production/test1/after_dehydrogenation/fbg5/r'][...]
    wavelength_after = f['fbg_production/test1/after_dehydrogenation/fbg5/wavelength'][...]
    ax.plot(1e9*wavelength_before,r_before,label='Hidrogenada')
    ax.plot(1e9*wavelength_after, r_after, label='Desidrogenada')
    ax.set_ylabel("Refletividade")
    ax.set_xlabel(r"$\lambda[\unit{\nm}]$")
    ax.set_xlim(left=1550,right=1560)
    ax.set_ylim(0,0.8)
    ax.legend()
    # ax.legend(ncol=2)




# f = h5py.File('./../data/phd_data.hdf5', 'r')

# for i in [4,5,6]:
#     _files = os.popen('du -a '+data_folder+"FBG"+str(i) +
#                       "*/*.lvm").read().split('\n')[:-1]
#     files = natsort.natsorted(_files)
#     # meta = os.popen('du -a '+data_folder+"FBG"+str(i) +
#     # "*/meta*").read().split('\n')[:-1][0].split('\t')[-1]
#     time = pd.read_csv(os.popen('du -a '+data_folder+"FBG"+str(i) +
#                                 "*/*time.txt").read().split('\n')[:-1][0].split('\t')[-1], names=['time', 'none'], sep='\t')
#     _time = fix_time(time.time)
#     data = pd.read_csv(
#         files[1].split('\t')[-1], sep="\t", names=['none', 'wavelength', 'optical_power'])

#     # Reference for reflectivity calculation
#     index_start_bias = 120
#     """start index for source corrections"""
#     index_stop_bias = 180
#     """start index for source corrections"""
#     index_r_start = np.where(data.wavelength*1e9 > 1550)[0][0]
#     index_r_stop = np.where(data.wavelength*1e9 > 1570)[0][0]
#     d0 = data.optical_power[:]
#     ff = f.require_group('fbg_production/test1_/after_dwhydrogenation'+str(i))
#     l_bragg = []
#     r_peak = []
#     time = []
#     for fn in range(15, len(files)):
#         leg = files[fn].split('/')[-1]
#         data = pd.read_csv(
#             files[fn].split('\t')[-1], sep="\t", names=['none', 'wavelength', 'optical_power'])
#         bias = (d0[index_start_bias:index_stop_bias] -
#                 data.optical_power[index_start_bias:index_stop_bias]).mean()
#         # ax.plot(1e9*data.wavelength, data.optical_power, label=leg)
#         r = reflectivity_transmition(d0, data.optical_power[:]+bias)
#         ff['r/'+"{:03d}".format(fn)] = r
#         ff['optical_power/'+"{:03d}".format(fn)] = data.optical_power
#         ff['wavelength/'+"{:03d}".format(fn)] = data.wavelength
#         l_bragg.append(
#             1e9*data.wavelength[index_r_start+r[index_r_start:index_r_stop].argmax()])
#         r_peak.append(r[index_r_start:index_r_stop].max())
#         time.append(_time[fn])
#         ff['r/'+"{:03d}".format(fn)].attrs['time'] = _time[fn]
#         ff['r/'+"{:03d}".format(fn)].attrs['l_bragg'] = l_bragg[-1]
#         ff['r/'+"{:03d}".format(fn)].attrs['r_peak'] = r_peak[-1]
# f.close()







def analysis_laser_fbg():
    f = h5py.File('./../data/phd_data.hdf5', 'r')
    # fbg_power_1 = f['fbg_production/test2/fbg_1']
    fbg_power_2 = f['fbg_production/test2/fbg_2']
    fbg_power_3 = f['fbg_production/test2/fbg_3']
    fbg_power_4 = f['fbg_production/test2/fbg_4']
    fbg_power_5 = f['fbg_production/test2/fbg_5']

    ## Dehydrogenated
    fbg_power_2_dehydrogenated = f['fbg_production/test2/fbg_2/after_dehydrogenation']
    fbg_power_5_dehydrogenated = f['fbg_production/test2/fbg_5/after_dehydrogenation']

    laser_power = f['osa_comparision/advantest02/high_resolution/']
    
    plt.cla()
    # plt.plot(fbg_power_1['wavelength'][:],
            #  fbg_power_1['r'][:,-1],label='Fbg_1')
    plt.plot(fbg_power_2['wavelength'][:],
             fbg_power_2['r'][:,-1],label='Fbg_2')
    plt.plot(fbg_power_3['wavelength'][:],
             fbg_power_3['r'][:,-1],label='Fbg_3')
    plt.plot(fbg_power_4['wavelength'][:],
             fbg_power_4['r'][:,-1],label='Fbg_4')
    plt.plot(fbg_power_5['wavelength'][:],
             fbg_power_5['r'][:,-1],label='Fbg_5')



    plt.plot(fbg_power_2_dehydrogenated['wavelength'][:],
             fbg_power_2_dehydrogenated['r'][:], label='Fbg_2 Desidrogenada')
    plt.plot(fbg_power_5_dehydrogenated['wavelength'][:],
             fbg_power_5_dehydrogenated['r'][:], label='Fbg_5 Desidrogenada')
    
    # plt.plot(fbg_power_after_dehydrogenation['wavelength'][:],
    #          fbg_power_after_dehydrogenation['r'][:],label='Desidrogenada')
    # plt.plot(fbg_power['wavelength'][:], fbg_power['r'][:,8],label='Hidrogenada')
    plt.plot(laser_power['wavelength'][:],
             laser_power['power_dbm'][:]/laser_power['power_dbm'][:].max(),label='Laser')
    plt.ylim(0,1)
    plt.legend()
    plt.xlim(1545e-9,1555e-9)