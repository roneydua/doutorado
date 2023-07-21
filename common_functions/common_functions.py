#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   common_functions.py
@Time    :   2023/07/08 18:37:43
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
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
