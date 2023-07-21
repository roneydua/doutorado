#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   acquisitionAnritsu.py
@Time    :   2023/03/02 17:19:25
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import time
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymeasure.instruments.anritsu import AnritsuMS9710C
from pyvisa import ResourceManager

from aquisitions.auxiliary_classes import Dynamometer, EncoderMikeController
from common_functions.common_functions import *

mpl.rcParams['figure.dpi'] = 72
plt.style.use("common_functions/roney3.mplstyle")

ILX_CONTROL = False
MIKE_CONTROL = False
DYNAMOMETER_CONTROL = False
AXIS_MIKE_ENCODER = 3  # 1, 2 or 3 for x, y, or z respectively
# H - hour, M- minute, S - second
folder_save = "./../experimentos/20230720/test_3/"
# letter = 'f'
teste_save_name = 'fbg_6_'


class experiment_data_save():

    def __init__(self, size: int):
        self.df = pd.DataFrame()
        """Data frame with experiment data"""
        self.df['wavelength'] = np.zeros(size, dtype=np.float32)
        """Wavelength vector in nano memeter"""
        self.df['power_dbm'] = np.zeros(size, dtype=np.float16)
        """Power vector in dBm * actual_resolution. WARNING: ned correction with 'actual resolution'"""
        # Variable necessary for the correction of the power value in decibels
        self.df['actual_resolution'] = '0.05'
        self.df['ilx_current'] = 400.0
        self.df['resolution_nm'] = 0.05
        self.df['resolution_vbw'] = '100Hz'

        self.df['traction_N'] = np.zeros(size, dtype=np.float16)
        
        self.df['room_temperature_C'] = 22.5 * np.ones(size, dtype=np.float16)

    def save(self, name):
        self.df.to_csv(name, index=True)


def save_temporary_graphic(da, graphic_name: str):
    fig, ax2 = plt.subplots(1, 1, num='temp', dpi=36, facecolor="white")
    ax2.plot(da.df['wavelength'], da.df['power_dbm'])
    ax2.set_xlabel(r"$\lambda, \si{\nm}$")
    ax2.set_ylabel('dBm')
    plt.savefig(graphic_name, format="png", transparent=False)
    plt.close(fig='temp')


def test_setup():
    osa = AnritsuMS9710C("GPIB0::8::INSTR", timeout=15000)
    osa.clear()
    osa.write("DATE " + str(datetime.now().year)[2:] + ",0" +
              str(datetime.now().month) + "," + str(datetime.now().day))
    osa.write("time " + str(datetime.now().hour) + "," +
              str(datetime.now().minute))
    osa.write("SMT OFF")
    # osa.write("STA 1510")
    osa.write("STA 1450")
    # osa.write("STO 1600")
    osa.write("STO 1650")
    
    # osa.write("LOG 4") # Y scale in LOG
    # osa.write("RLV -30") # Ref level
    osa.resolution_vbw = '1kHz'
    # osa.write("STA 1525")
    osa.write("MPT 5001")
    # osa.write('SSI')
    osa.resolution = experiment_data_save(10).df['resolution_nm'][0]
    # Enabel true resulution
    osa.write('ARES ON')  # osa.write("SRT")
    # osa.clear()
    osa.write(r'ARED?')
    # ILX setup
    if ILX_CONTROL:
        rm = ResourceManager()
        ilx = rm.open_resource('GPIB0::1::INSTR')
        # put TEC on temperature mode
        ilx.write("TEC:MODE:T<NL>")  # type: ignore
        # enable tec
        ilx.write("TEC:OUT 1<NL>")  # type: ignore
        # enable laser output
        ilx.write("LAS:OUT 1")  # type: ignore
    # ilx.write("LAS:OUT 0")
    # set ILX current
    mike = EncoderMikeController(39.7e-3)
    # define velocity
    mike.sendCommand(AXIS_MIKE_ENCODER, "SV", 1)
    # instance of dynamometer
    dyn = Dynamometer()

    return osa, dyn, mike


def make_test(osa, dyn, mike, manual_micrometer: bool or float = False):
    # osa.single_sweep(n=10)
    # performe a singe sweep
    osa.write("ssi")
    # Ensure end of read
    while True:
        _t = osa.esr2
        time.sleep(1.)
        if _t != 0:
            break
    # update time to save data
    now = datetime.now()
    current_time = now.strftime("%m%d%y__%H%M%S")
    # da = dataAcquisition(n=osa.sampling_points, test_name="TEAP")
    test_name = Path(folder_save + teste_save_name + current_time + ".csv")
    # da = dataAcquisition(n=osa.sampling_points) # type: ignore
    sampling_points = int(osa.sampling_points)  # type: ignore
    da = experiment_data_save(size=sampling_points)
    # update Actual resolution, resolution and other variables
    da.df['actual_resolution'] = float(osa.ask("ARED?"))  # actual_resolution
    da.df['resolution_nm'] = float(osa.ask("RES?"))
    da.df['resolution_vbw'] = osa.resolution_vbw
    # osa.wait_for_sweep()
    da.df['wavelength'], da.df['power_dbm'] = osa.read_memory(slot='A')
    # Dynamomenter data
    da.df['traction_N'] = dyn.update_strain()
    # Encoder data
    da.df['initial_length_m'] = mike.initial_length_m
    if manual_micrometer == False:
        da.df['micrometer_position_um'] = mike.getPositionZ()
    else:
        da.df['micrometer_position_um'] = manual_micrometer
    # Save all data of test
    da.save(test_name)
    # Save figure
    save_temporary_graphic(
        da, folder_save + teste_save_name + current_time +
        ".png")
    plt.plot(da.df['wavelength'], da.df['power_dbm'], label=teste_save_name)
    plt.legend()
    # osa.write("SRT")


def main_with_oriel_encoder():
    osa, dyn, mike = test_setup()
    # change velocity temporarily
    
    positions = np.arange(20,121,20)
    # plt.plot(positions)
    # put low velocity
    mike.sendCommand(AXIS_MIKE_ENCODER, "SV", 10.0)
    # Put current on 100mA
    for i in positions:
        print(i)
        mike.goToZ(i)
        time.sleep(3)
        print("make_test")
        make_test(osa, dyn, mike,manual_micrometer=False)
        
    mike.goToZ(0)
    



def manual_with_vernier():
    '''
    manual_with_vernier Method used for test with the manual micrometer.
    '''    
    osa, dyn, mike = test_setup()
    # change velocity temporarily
    positions = np.arange(20,221,20)
    # put low velocity
    mike.sendCommand(AXIS_MIKE_ENCODER, "SV", 20.0)
    # Put current on 100mA
    for i in positions:
        _p = input("type positions...")
        print(i)
        # mike.goToZ(i)
        # time.sleep(3)
        print("make_test")
        make_test(osa, dyn, mike, manual_micrometer=i)
        
    mike.goToZ(0)