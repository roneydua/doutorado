#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   acquisitionAnritsu.py
@Time    :   2023/03/02 17:19:25
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import pyvisa
import time
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from pymeasure.instruments.anritsu import AnritsuMS9710C
from auxiliaryClasses import Dynamometer
import pandas as pd
from pathlib import Path
mpl.rcParams['figure.dpi'] = 72
plt.style.use("../../../../programasComuns/roney3.mplstyle")

# for save data
from datetime import datetime
ilx_control = True
mike_control =False
dynamometer_control = False
# H - hour, M- minute, S - second
folder_save = './data/12052023/'
# letter = 'f'
# teste_save_name = letter+"_reflection_"+"4_percent_"
# teste_save_name = letter+"_reflection_4_percent_"
teste_save_name = 'laser_bobina_21mm_M12_'

class dataAcquisition:

    def __init__(self, n=1, test_name='') -> None:
        self.test_name = test_name
        # With a test_name, the loading step is called
        if test_name != '':
            f = self.load_states(name=test_name)
            self.wave_length = f.wave_length
            self.power = f.power
            self.power_watt = f.power_watt
        else:
            self.wave_length = np.zeros(n)
            self.power = np.zeros(n)
            self.power_watt = np.zeros(n)
            self.enc_position = np.zeros(n)
            self.dyn = np.zeros(n)

    def dbm2W(self):
        self.power_watt = 10.0**(self.power * 0.1)

    def save_states(self, test_name):
        # save power in Watts
        self.test_name = test_name
        self.dbm2W()
        with open(self.test_name, 'wb') as handle:
            pickle.dump(self, handle)

    def load_states(self, name):
        with open(name, 'rb') as f:
            a = pickle.load(f)
        return a


def dbm2W(power):
    return 10.0**(power * 0.1)

class experiment_data_save():

    def __init__(self, size: int, actual_resolution=1.0):
        self.df = pd.DataFrame()
        """Data frame with experiment data"""

        """Wavelength vector in nano memeter"""
        self.df['wave_length'] = 1.0 + np.zeros(size, dtype=np.float32)
        """Power vector in dBm * actual_resolution. WARNING: ned correction with 'actual resolution'"""

        self.df['power_dbm'] = 2.0 + np.zeros(size, dtype=np.float16)
        "" "Variable necessary for the correction of the power value in decibels" ""
        self.df['actual_resolution'] = actual_resolution
        self.ilx_current = 0.0
        self.resolution_nm = 1.0
        self.resolution_vbw = '100Hz'

    def save(self, name):
        # fix dbm power with current actual resolution

        # self.df['power_dbm'] -= 10.0 * np.log10(self.df['actual_resolution'])
        self.df.to_csv(name, index=True)

def test_setup():
    osa = AnritsuMS9710C("GPIB0::8::INSTR", timeout=15000)
    osa.clear()
    osa.write("DATE " + str(datetime.now().year)[2:] + ",0" +
              str(datetime.now().month) + "," + str(datetime.now().day))
    osa.write("time " + str(datetime.now().hour) + "," +
              str(datetime.now().minute))
    return osa

def save_temporary_graphic(da, graphic_name:str):
    fig, ax2 = plt.subplots(1, 1, num='temp',dpi=36)
    ax2.plot(da.df['wave_length'], da.df['power'])
    ax2.set_xlabel(r"$\lambda, \si{\nm}$")
    ax2.set_ylabel('dBm')
    plt.savefig(graphic_name, format="pdf")
    plt.close(fig='temp')


def make_test(ilx_current:float):
        ## ILX setup
    ilx_current = '{:2.1f}'.format(ilx_control)
    if ilx_control:
        rm = pyvisa.ResourceManager()
        ilx = rm.open_resource('GPIB0::1::INSTR')
        # put TEC on temperature mode
        ilx.write("TEC:MODE:T<NL>")
        # enable tec
        ilx.write("TEC:OUT 1<NL>")
        # set ILX current

        ilx.write("LAS:LDI " + ilx_current)
        # enable laser output
        ilx.write("LAS:OUT 1")
        # ilx.write("LAS:OUT 0")
    if dynamometer_control:
        ## Dynamometer
        dyn = Dynamometer()
        dyn.set_zero()
    ## ORIEL Encoder Mike Controller.
    # rm = pyvisa.ResourceManager()
    if mike_control:
        enc_mike = rm.open_resource('GPIB0::3::INSTR', timeout=2 , write_termination = '\r', read_termination='\r')
        # enc_mike.close()
        print(rm.list_resources('GPIB*'))
        enc_mike.write('3PM')
        enc_mike.write('3SV1')
        enc_mike.write('3MR-10.0')
        enc_mike.write('3ST')
        enc_mike.write('3DH0')
        enc_mike.timeout = 5000
        # Set position mode
        # enc_mike.write('3TP')

    ## Anritsu OSA
    
    # osa.write("aut")
    # osa.center_at_peak()
    # osa.write("PKS PEAK")
    # turn on smooth filter
    osa.write("SMT OFF")
    # osa.write("STA 1548")
    osa.wavelength_span = 2.0
    osa.write("CNT 1549.0")
    osa.resolution_vbw = '1kHz'
    # osa.write("STA 1525")
    # osa.write("STO 1560")
    resolution = '.05'
    osa.write("RES "+resolution)
    osa.write("MPT 2001")
    # osa.write('SSI')
    # get current resolution
    osa.write('ARES ON')  # osa.write("SRT")
    # osa.clear()
    osa.write(r'ARED?')
    actual_resolution = float(osa.read())

    # teste_save_name = "TEAP_"
    osa.single_sweep(n=100)
    # update time to save data
    now = datetime.now()
    current_time = now.strftime("%m%d%y__%H%M%S")
    # da = dataAcquisition(n=osa.sampling_points, test_name="TEAP")
    test_name = Path(folder_save + teste_save_name + current_time+ ".csv")
    # da = dataAcquisition(n=osa.sampling_points) # type: ignore
    sampling_points = int(osa.sampling_points) # type: ignore
    da = experiment_data_save(size=sampling_points)
    ## update Actual resolution, resolution and other variables
    da.df['actual_resolution'] = actual_resolution
    da.df['ilx_current'] = ilx_current
    da.df['resolution_nm'] = float(resolution)
    da.df['resolution_vbw'] = osa.resolution_vbw
    # osa.wait_for_sweep()
    da.df['wave_length'], da.df['power'] = osa.read_memory(slot='A')

    da.save(test_name)
    save_temporary_graphic(
        da, folder_save + teste_save_name + current_time +
        ".pdf")
    # da.save_states(da)
    # To Load
    # da = dataAcquisition(test_name='data/A_032023__080225_500mA')
    plt.plot(da.df['wave_length'], da.df['power'], label=teste_save_name)
    plt.legend()
    osa.write("SRT")



if __name__ == "__main__":
    osa = setup_osa_anritsu()






# da = pd.read_csv(Path('/home/pegasus/Dropbox/doutorado/experimentos/25042023/laser_1549042523__131735.csv'))
# plt.plot(da['wave_length'],da['power'])
# plt.xlim([1548.5,1550])
# plt.ylim(bottom=-60,top=8)
# plt.savefig('laser_1594.jpg',format='jpg')


# test_name1 = Path("../../../../experimentos/20042023/a_reflection042423__130528.csv")
# test_name2 = Path("../../../../experimentos/20042023/a_reflection_4_percent_042423__130752.csv")

# da1 = pd.read_csv(test_name1)
# da2 = pd.read_csv(test_name2)

# plt.plot(da1['wave_length'].to_numpy(),da1['power'].to_numpy())
# plt.plot(da2['wave_length'].to_numpy(),da2['power'].to_numpy())

# plt.plot(da1['wave_length'],
#          (dbm2W(da2['power'])- dbm2W(da1['power']) ) / dbm2W(da2['power']))


