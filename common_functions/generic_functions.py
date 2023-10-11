#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   generic_functions.py
@Time    :   2023/07/08 18:37:43
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import pandas as pd


# import sympy as sp


def dbm_mW(dbm:float or np.ndarray, mW=True):
    '''
    dbm2mW transform a dbm to mW or Watt
    Args:
        dbm: value to convert to mw. Can be float or numpy array.
        mw: _description_. Defaults to True.

    Returns:
        _description_
    '''    
    if mW==True:
        return 10.0**(0.1*dbm)
    else:
        return 1e-3 * (10.0**(0.1*dbm))
    

def mW_dbm(mW: float or np.ndarray):
    '''
    dbm2mW transform a dbm to mW or Watt
    Args:
        dbm: value to convert to mw. Can be float or numpy array.
        mw: _description_. Defaults to True.

    Returns:
        _description_
    '''
    return 10.0*np.log10(mW*1e3)


def calc_laser(w, center_w, std=0.01):
    '''
    calc_laser make a laser with Gaussian profile.
    Args:
        w: wavelength vector 
        center_w: wavelength peak
        std: Standart deviation of Gaussian profile. Defaults to 0.1*1e-9.
    Returns:
        A vector of laser
    '''
    return np.exp(-(w - center_w)**2 /
                  (2.0 * std**2)) / (std * np.sqrt(2.0 * np.pi))


def find_index_of_x_span(min:float, max:float, x:np.array):
    '''
    find_index_of_x_span gives the minimum and maximum indices of the wavelengths.
    Args:
        min: minimum index
        max: index of max
        x: n dimension vector
    Returns:
        index
    '''
    try:
        index_min = np.where(x>min)[0][0]
    except:
        index_min = 0
    try:
        index_max = np.where(x>max)[0][0]
    except:
        index_max = len(x)-1
    return index_min, index_max


def reflectivity_transmition(d0: np.array, di: np.array):
    '''
    reflectivity_transmition Calc reflectivity with transmition method
    Args:
        d0: Source in dbm
        di: transmited power in dbm

    Returns:
        reflectivity vector [np.array]
    '''
    return 1.0 - (10**(0.1 * (di-d0)))


def calc_reflectivity_by_transmission(source: np.array, power: np.array, wavelength: None or np.array, normalize_source =False, **kwargs):
    '''
    calc_reflectivity_by_transmission compute reflectivity by transmistion method.

    Args:
        wavelength: wavelength compriment in meters
        source: power of source in dbm
        power: power of transmited power
        normalize_source: If the source needs normalization with respecti a specific range.
            In this case, it is necessary to inform the minimum and maximum values of wavelengths. Defaults to False.

    Returns:
        _description_
    '''
    if normalize_source:
        min_w = kwargs.pop('min_wavelength')
        max_w = kwargs.pop('max_wavelength')
        #check with wavelength is in meters or in nm
        index_min, index_max = find_index_of_x_span(min_w,max_w,wavelength)
        bias_source = (source[index_min:index_max]-power[index_min:index_max]).mean()
        return reflectivity_transmition(source-bias_source, power)
    
    return reflectivity_transmition(source,power)


def calc_bias(self, min_w:float, max_w:float, wavelength:np.ndarray,source:np.ndarray,power:np.ndarray):
    index_min, index_max = find_index_of_x_span(min_w, max_w, wavelength)
    # bias_source = (source[index_min:index_max] -
                #    power[index_min:index_max]).mean()
    bias_source = np.cumsum(source[index_min:index_max] + power[index_min:index_max])*0.5
    return bias_source



def gabriel(x:float,y:float):
    return x*y