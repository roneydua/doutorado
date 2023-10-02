#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   fbgs_production_files.py
@Time    :   2023/10/02 09:36:43
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import sympy as sp
import locale
import matplotlib.pyplot as plt
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
plt.style.use("common_function/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4


def put_data_production_on_hdf_20230807():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230807/"
    for i in [1, 2, 3]:
        folder = FOLDER+"FBG"+str(i)+"_07Ago2023"
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
        f = h5py.File("production_files.hdf5", "a")
        ff = f.require_group("fbg_production/20230807/fbg"+str(i))
        ff["reflectivity"] = reflectivity
        ff["wavelength_m"] = wavelength
        ff["optical_power_dbm"] = optical_power
        with open(folder+"/metadata.txt") as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230817():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230817/"
    for i in [5]:
        file = FOLDER+"FBG"+str(i)+".txt"
        data_pd =  pd.read_csv(file,sep="\t")
        f = h5py.File("production_files.hdf5","a")
        ff =f.require_group("fbg_production/20230817/fbg"+str(i))
        ff["wavelength_m"] = data_pd.iloc[:,0]
        ff["optical_power_dbm"] = data_pd.iloc[:,1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0],ff["optical_power_dbm"][:,j])
        ff["reflectivity"] = refletictivity
        # with open(folder+"/metadata.txt") as metadata:
        #     ff.attrs['metadata']=metadata.readlines()
        f.close()
        
        
def put_data_production_on_hdf_20230810():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230810/"
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
        f = h5py.File("production_files.hdf5","a")
        ff = f.require_group("fbg_production/20230810/fbg"+str(i))
        ff["reflectivity"] = reflectivity
        ff["wavelength_m"] = wavelength
        ff["optical_power_dbm"] = optical_power
        with open(folder+"/metadata.txt") as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230814():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230814/"
    for i in [2, 3, 4, 5, 6]:
        folder = FOLDER+"FBG"+str(i)+"_14Ago2023_trans"
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
        f = h5py.File("production_files.hdf5", "a")
        ff = f.require_group("fbg_production/20230814/fbg"+str(i))
        ff["reflectivity"] = reflectivity
        ff["wavelength_m"] = wavelength
        ff["optical_power_dbm"] = optical_power
        with open(folder+"/metadata.txt") as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230921():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230921/"
    for i in [2, 3, 4, 5]:
        file = FOLDER+"FBG#"+str(i)+".txt"
        data_pd = pd.read_csv(file, sep="\t")
        f = h5py.File("production_files.hdf5", "a")
        ff = f.require_group("fbg_production/20230921/fbg"+str(i))
        ff["wavelength_m"] = data_pd.iloc[:, 0]
        ff["optical_power_dbm"] = data_pd.iloc[:, 1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0], ff["optical_power_dbm"][:, j])
        ff["reflectivity"] = refletictivity
        with open(FOLDER+"metadata"+"_FBG#"+str(i)+".txt", errors='ignore') as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230822():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230822/"
    for i in [1,2,3,4,5]:
        file = FOLDER+"FBG#"+str(i)+".txt"
        data_pd =  pd.read_csv(file,sep="\t")
        f = h5py.File("production_files.hdf5","a")
        ff =f.require_group("fbg_production/20230822/fbg"+str(i))
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


def put_data_production_on_hdf_20230904():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230904/"
    for i in ["2", "3_lente_cilindrica", "4_lente_cilindrica", "5_lente_cilindrica"]:
        file = FOLDER+"FBG#"+i+".txt"
        data_pd = pd.read_csv(file, sep="\t")
        f = h5py.File("temp.hdf5", "a")
        ff = f.require_group("fbg_production/20230904/fbg"+i)
        ff["wavelength_m"] = data_pd.iloc[:, 0]
        ff["optical_power_dbm"] = data_pd.iloc[:, 1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0], ff["optical_power_dbm"][:, j])
        ff["reflectivity"] = refletictivity
        with open(FOLDER+"metadata"+"_FBG#"+i+".txt", errors='ignore') as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230911():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230911/"
    for i in [1, 2, 3]:
        file = FOLDER+"FBG#"+str(i)+".txt"
        data_pd = pd.read_csv(file, sep="\t")
        f = h5py.File("production_files.hdf5", "a")
        ff = f.require_group("fbg_production/20230911/fbg"+str(i))
        ff["wavelength_m"] = data_pd.iloc[:, 0]
        ff["optical_power_dbm"] = data_pd.iloc[:, 1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0], ff["optical_power_dbm"][:, j])
        ff["reflectivity"] = refletictivity
        with open(FOLDER+"metadata"+"_FBG#"+str(i)+".txt", errors='ignore') as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230912():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230912/"
    for i in [1, 2]:
        file = FOLDER+"FBG#"+str(i)+".txt"
        data_pd = pd.read_csv(file, sep="\t")
        f = h5py.File("production_files.hdf5", "a")
        ff = f.require_group("fbg_production/20230912/fbg"+str(i))
        ff["wavelength_m"] = data_pd.iloc[:, 0]
        ff["optical_power_dbm"] = data_pd.iloc[:, 1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0], ff["optical_power_dbm"][:, j])
        ff["reflectivity"] = refletictivity
        with open(FOLDER+"metadata"+"_FBG#"+str(i)+".txt", errors='ignore') as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()


def put_data_production_on_hdf_20230913():
    '''
    put_data_production_on_hdf Funtion to append csv files collected on production to a hdf file. 
    '''
    FOLDER = "../data/danteAlex/20230913/"
    for i in [1, 2]:
        file = FOLDER+"FBG#"+str(i)+".txt"
        data_pd = pd.read_csv(file, sep="\t")
        f = h5py.File("production_files.hdf5", "a")
        ff = f.require_group("fbg_production/20230913/fbg"+str(i))
        ff["wavelength_m"] = data_pd.iloc[:, 0]
        ff["optical_power_dbm"] = data_pd.iloc[:, 1:]
        refletictivity = np.zeros(ff["optical_power_dbm"].shape)
        for j in range(ff["optical_power_dbm"].shape[1]):
            refletictivity[:, j] = reflectivity_transmition(
                ff["optical_power_dbm"][:, 0], ff["optical_power_dbm"][:, j])
        ff["reflectivity"] = refletictivity
        with open(FOLDER+"metadata"+"_FBG#"+str(i)+".txt", errors='ignore') as metadata:
            ff.attrs['metadata'] = metadata.readlines()
        f.close()
