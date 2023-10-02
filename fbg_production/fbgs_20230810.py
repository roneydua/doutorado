#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   fbgs_20230814.py
@Time    :   2023/10/01 22:39:22
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import sympy as sp
import h5py
import pandas as pd
import os
import locale
import matplotlib.pyplot as plt
import natsort
from common_functions.generic_functions import reflectivity_transmition
# from IPython.core.interactiveshell import InteractiveShell
# from ipywidgets import interactive, fixed
# InteractiveShell.ast_node_interactivity = "all"
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
# plt.style.use("common_function/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4


FOLDER = "../data/danteAlex/20230810/"

def put_data_production_on_hdf():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    for i in [1]:
        folder = FOLDER
        _files = os.popen('du -a '+folder+"/*.lvm").read().split('\n')[:-1]
        files = natsort.natsorted(_files)
        n = files.__len__()
        wavelength = pd.read_csv(files[0].split("\t")[1], sep="\t").iloc[:, 1]
        optical_power = np.zeros((wavelength.size, n))
        reflectivity = np.zeros((wavelength.size, n))
        for j in range(n):
            # read all data
            _data = pd.read_csv(files[j].split("\t")[1], sep="\t")
            optical_power[:, j] = _data.iloc[:, 2]
            reflectivity[:, j] = reflectivity_transmition(
                optical_power[:, 0], optical_power[:, j])
        f = h5py.File(production_files.hdf5","a")
        ff =f.require_group("fbg_production/20230810/fbg"+str(i))
        ff["reflectivity"] = reflectivity
        ff["wavelength_m"] = wavelength
        ff["optical_power_dbm"] = optical_power
        with open(folder+"/metadata.txt") as metadata:
            ff.attrs['metadata']=metadata.readlines()
        f.close()