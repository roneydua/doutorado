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


def put_data_production_on_hdf():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230911/"
    for i in [1,2,3]:
        file = FOLDER+"FBG#"+str(i)+".txt"
        data_pd =  pd.read_csv(file,sep="\t")
        f = h5py.File("production_files.hdf5","a")
        ff =f.require_group("fbg_production/20230911/fbg"+str(i))
        ff["wavelength_m"] = data_pd.iloc[:,0]
        ff["optical_power_dbm"] = data_pd.iloc[:,1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0],ff["optical_power_dbm"][:,j])
        ff["reflectivity"] = refletictivity
        with open(FOLDER+"metadata"+"_FBG#"+str(i)+".txt",errors='ignore') as metadata:
            ff.attrs['metadata']=metadata.readlines()
        f.close()