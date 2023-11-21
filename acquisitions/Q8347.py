#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   Q8347.py
@Time    :   2023/10/09 19:20:41
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''


from pathlib import Path
import pyvisa
from numpy import savetxt, asfarray, column_stack
import matplotlib.pyplot as plt
from time import sleep  # função para verificar o status byte do aparelho


class Q8347(object):
    '''
    Q8347 Classe Q8347 com GPIB

    _extended_summary_

    Args:
        object: _description_
    '''

    def __init__(self, center=1550, gpib_adress='GPIB0::8::INSTR', span=None, high_resolution=True):
        '''
        __init__ constructor of Q8347 Advantest

        _extended_summary_

        Args:
            gpib_adress: Defaults to 'GPIB0::8::INSTR'.
        '''
        # super(Q8347, self).__init__()
        rm = pyvisa.ResourceManager('@py')
        # rm.list_resources()
        self.osa = rm.open_resource('GPIB0::8::INSTR')  # osa livre
        self.osa.chunk_size = 65535  # configuraçoes da comunicaçao
        self.osa.timeout = 20000  # configuraçoes da comunicaçao
        self.osa.read_termination = '\n'  # configuraçoes da comunicaçao
        if center != None:
            self.set_center(center)
        if span != None:
            self.set_span(span)
        if high_resolution:
            self.set_resolution(high=True)
        else: 
            self.set_resolution(high=False)

    def checkSTB(self, t=1):  # FUNÇÃO BASEADA NO "wait of spectrometer" do Gabriel
        sleep(1)  # Como o Advantest nao tem a função do GPIB (SRQ),
        i = True  # é preciso usar uma funçõa de mais baixo nivel
        while i:
            stb = self.osa.read_stb()
            if stb == 1:  # Quando o STB é igual a 1, o Advantest terminou
                i = False  # de fazer a ação que havia sido solicitada
            else:
                sleep(t)
        return

    def read(self):
        '''
        Read () performs a data acquisition and updates the wavelength_NM variables with wavelength and optical_power_dbm with dbm power.

        '''
        self.osa.write('MEA1')  # make a sigle measurement
        self.checkSTB(0.5)  # cal the checkSTB() function and wait read ends
        x = self.osa.query('OSD1')  # get wavelength
        y = self.osa.query('OSD0')  # get optical power
        # remove reader
        x = x[5:]
        y = y[5:]
        # convert to numpy arrays
        self.wavelength_m = asfarray(x.split(','))
        self.optical_power_dbm = asfarray(y.split(','))
        
    def set_resolution(self, high=True):
        '''
        set_resolution Set resolution mode
        Args:
            high: mode of resolution. Defaults to True.
        '''
        if high:
            self.osa.write('RES 1')
        else:
            self.osa.write('RES 0')

    def set_span(self, span: float, unit='NM'):
        '''
        set_span Set span of reads

        Args:
            span: Span
            unit: string with values UM, NM
        '''
        self.span = span
        self.osa.write('SPA' + str(span) + unit)

    def set_center(self, center: float, unit='NM'):
        '''
        set_center Set center of measurement


        Args:
            center: flat of wavelength center
            unit: The unit of center. Defaults to 'NM'.
        '''
        self.center = center
        self.osa.write('CEN'+str(center)+unit)

    def close(self):
        self.osa.close()