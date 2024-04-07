#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   AQ6370D.py
@Time    :   2024/02/06 20:49:41
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""
from inspect import Traceback
from pathlib import Path
import pyvisa
from numpy import asfarray, column_stack
import matplotlib.pyplot as plt
from time import sleep  # função para verificar o status byte do aparelho


class Trace(object):
    def __init__(self, trace_name: str):
        """__init__ class with trace's variables as wavelength and power
        Args:
            trace_name: a string of name of trace. TRA for trace A, TRB for trace B and etc. This OSA have 8 traces, so we have TRA, TRB, ... TRH.
        """
        super(Trace, self).__init__()
        self.trace_name = trace_name


# class Traces(Trace):
#     def __init__(self, trace_name: str):
#         """__init__ class with trace's variables as wavelength and power
#         Args:
#             trace_name: a string of name of trace. TRA for trace A, TRB for trace B and etc. This OSA have 8 traces, so we have TRA, TRB, ... TRH.
#         """
#         super(Traces, self).__init__(trace_name)
#         self.trace_name = trace_name
#         self = Trace(trace_name)





class AQ6370D(Trace):
    """
    AQ6370D Classe AQ6370D com GPIB

    _extended_summary_

    Args:
        object: _description_
    """

    def __init__(
        self,
        center=1550,
        gpib_address="GPIB0::1::INSTR",
        span=None
    ):
        """
        __init__ constructor of AQ6370D Advantest
        Args:
            gpib_address: Defaults to 'GPIB0::8::INSTR'.
        """
        # super(AQ6370D, self).__init__()
        rm = pyvisa.ResourceManager("@py")
        # rm = pyvisa.ResourceManager()
        # print(rm.list_resources())
        self.osa = rm.open_resource(gpib_address)  # osa livre
        # self.osa.chunk_size = 65535  # comunitaiton setup
        self.osa.timeout = 200_000  # comunitaiton setup
        # self.osa.read_termination = "\r\n"  # comunitaiton setup
        if center != None:
            self.set_center(center)
        if span != None:
            self.set_span(span)

        self.trace = {
            TRA: False,
            TRB: False,
            TRC: False,
            TRD: False,
            TRE: False,
            TRF: False,
            TRG: False,
            TRH: False,
        }
        # self.trace
        
        

    def checkSTB(self, t=1):  # FUNÇÃO BASEADA NO "wait of spectrometer" do Gabriel
        sleep(1)  
        i = True 
        while i:
            stb = self.osa.read_stb()
            if stb == 1:  # Quando o STB é igual a 1, o Advantest terminou
                i = False  # de fazer a ação que havia sido solicitada
            else:
                print("waiting")
                sleep(t)
        return

    def set_resolution(self, resolution='20pm'):
        '''set_resolution Set measurement resolution

        Args:
            resolution. Defaults to '20pm'.
        '''
        self.osa.write(":SENSE:BANDWIDTH:RESOLUTION "+resolution)
        self.resolution = self.osa.query(":sense:bandwidth?")
        print("current resolution ", self.resolution)

    def simple_sweep(self, trace="tra"):
        """
            simple_sweep() performs a data acquisition and updates the wavelength_NM variables with wavelength and optical_power_dbm with dbm power.
        Args:
                trace: a string with trance name trX with X a letter of trace. For example, trc for trace C.
                    Defaults to trace A.
        """
        
        self.osa.write(":INITIATE")  # make a sigle measurement
        # self.checkSTB(1)
        self.read(trace=trace)

    def read(self, trace='tra'):
        """
        Read () performs a data acquisition and updates the wavelength_NM variables with wavelength and optical_power_dbm with dbm power.
    Args:
            trace: a string with trance name trX with X a letter of trace. For example, trc for trace C.
                Defaults to trace A.
        """
        # self.checkSTB(1)
        self.x = self.osa.query(":TRACE:X? "+trace)  # get wavelength
        self.y = self.osa.query(":TRACE:Y? "+trace)  # get optical power
        # remove reader
        # x = x[10:]
        # y = y[10:]
        # convert to numpy arrays
        self.wavelength_m = asfarray(self.x.split(","))
        self.optical_power_dbm = asfarray(self.y.split(","))

    def set_span(self, span: float, unit="M"):
        """
        set_span Set span of reads

        Args:
            span: float of wavelength center in meters or hz
            unit: The unit of span. Can be [M, HZ]. Defaults to 'M'.
        """
        self.span = span
        self.osa.write(":SENSE:WAVELENGTH:SPAN " + str(span) + str(unit))

    def set_center(self, center: float, unit="M"):
        """
        set_center Set center of measurement


        Args:
            center: float of wavelength center in meters or hz
            unit: The unit of center. Can be [M, HZ]. Defaults to 'M'.
        """
        self.center = center
        self.osa.write(":SENSE:WAVELENGTH:CENTER " + str(center) + unit)

    def close(self):
        self.osa.close()


a = {"tra":Trace(trace_name="TRA")}
a["trb"] = Trace(trace_name="TRB")
'trc' in a.keys()

trace = {
    TRA: False,
    TRB: False,
    TRC: False,
    TRD: False,
    TRE: False,
    TRF: False,
    TRG: False,
    TRH: False,
}
