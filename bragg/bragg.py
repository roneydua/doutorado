#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   bragg.py
@Time    :   2023/02/26 13:53:14
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''
import numpy as np
import sympy as sp
figL = 6.29
figA = (90.0) / 25.4


class Bragg(object):
    """docstring for Bragg."""
    # grating_period = 1e-6
    d_original = 0.0
    n = 0
    #! refraction index
    n = 1.45
    # refraction index perturbation
    delta_n = 0
    number_of_grating_period = 0
    S_even = np.zeros(shape=(2, 2), dtype=np.complex128)
    S_odd = np.zeros(shape=(2, 2), dtype=np.complex128)
    wavelength_peak = 1550e-9

    def __init__(self, fbg_size, wavelength_peak=1550.0, delta_n=1e-3,delta_span_wavelength=10, diff_of_peak=1,number_of_grating_period_forced=None, wavelength_span=None):
        self.delta_n = delta_n
        """Refraction index perturbation"""
        self.wavelength_peak = 1e-9 * wavelength_peak

        # To facilitate the interrogation of FBG can be passed as the wavelength vector argument.
        if type(wavelength_span) != np.ndarray:
            self.wavelength_span = 1e-9 * np.linspace(
                1e9 * self.wavelength_peak - diff_of_peak,
                1e9 * self.wavelength_peak + diff_of_peak, delta_span_wavelength)
        else:
            self.wavelength_span = wavelength_span
        
        self.wavelength_span_nm = self.wavelength_span * 1e9

        self.fbg_size = fbg_size
        self.n_eff = 2.0*self.n*(self.n+self.delta_n) / (2.0*self.n+self.delta_n)
        """Effective refractive index"""

        self.d_original = self.wavelength_peak / (4.0 * self.n_eff)
        # self.grating_period = self.fbg_size / self.d_original
        # self.d_original = 1e-9 * wavelength_peak * (self.n+self.delta_n*0.5) * 0.5
        if number_of_grating_period_forced is None:
            self.number_of_grating_period = int((int(self.fbg_size / self.d_original) + 1) / 2.0)
        else:
            self.number_of_grating_period  = number_of_grating_period_forced
        #  self.number_of_grating_period = 1500
        # reflection coefficient
        ## NOTE: For the even reflection coefficient the negative of self.r_odd is used.
        self.r_odd = self.delta_n / (2.0 * self.n + self.delta_n)
        self.r_even = -self.r_odd
        # transmition index
        self.t_odd = 2 * self.n / (2.0 * self.n + self.delta_n)
        # even %2 = 0
        self.t_even = 2 * (self.n + self.delta_n) / (2.0 * self.n +
                                                     self.delta_n)
        # compute reflection with no deformation
        self.r0 = self.calc_bragg()


    def calc_S(self, delta, r, t: float):
        return np.array([[np.exp(1j * delta), r * np.exp(1j * delta)],[r * np.exp(-1j * delta), np.exp(-1j * delta)]], dtype=np.complex128)/t

    def calc_bragg(self, deformation=0.0, wavelength_vector=None):
        d = (deformation+1.0)*self.d_original
        if wavelength_vector is None:
            wavelength_vector = self.wavelength_span
        reflectance = np.zeros(len(wavelength_vector))
        # even %2 = 0
        delta_even = (self.n + self.delta_n) * (2.0 * np.pi /wavelength_vector) * d
        delta_odd = self.n * (2.0 * np.pi / wavelength_vector) * d
        # Walk the self.wave length vector
        for i in range(len(wavelength_vector)):
            # even %2 = 0
            self.S_even = self.calc_S(delta_even[i],self.r_even,self.t_even)
            self.S_odd = self.calc_S(delta_odd[i],self.r_odd,self.t_odd)
            S_temp = self.S_even @ self.S_odd
            self.S = (np.array([[1.0, self.r_odd], [self.r_odd, 1.0]])@
            np.linalg.matrix_power(S_temp, self.number_of_grating_period - 1)) @ (self.S_even)
            reflectance[i] = self.calc_reflectance(self.S)
            # reflectance[i] = self.S[1, 0] / self.S[0, 0]
        return reflectance

    def calc_reflectance(self,_S):
        return np.linalg.norm(_S[1, 0] / _S[0, 0])**2

    def calc_transmisivity(self, _S):
        return np.linalg.norm(1.0 / _S[0, 0])**2

    # def reflection_of_transmition():

    def reflection_of_transmition(self, deformation: float, wavelength: float):
        d = (deformation+1.0)*self.d_original
        # even %2 = 0
        delta_even = (self.n + self.delta_n) * (2.0 * np.pi /wavelength) * d
        delta_odd = self.n * (2.0 * np.pi / wavelength) * d
        # even %2 = 0
        self.S_even = self.calc_S(delta_even,self.r_even,self.t_even)
        self.S_odd = self.calc_S(delta_odd,self.r_odd,self.t_odd)
        S_temp = self.S_even @ self.S_odd
        return (np.array([[1.0, self.r_odd], [self.r_odd, 1.0]])@ np.linalg.matrix_power(S_temp, self.number_of_grating_period - 1)) @ (self.S_even)





class OpticalCoupler(object):
    """docstring for opticalCoupler."""
    """ Transfer matrix of optical coupler"""
    tm = np.zeros((2,2),dtype=np.complex128)
    def __init__(self, e: float):
        self.tm[0,0] = np.sqrt(1.0-e)
        self.tm[0,1] = 1j*np.sqrt(e)
        self.tm[1, 0] = self.tm[0, 1]
        self.tm[1,1] = self.tm[0,0]
