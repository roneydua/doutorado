#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 10:01:20 2021

@author: roney
"""

import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
import locale
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")


def printLatex(_b, printResult=True):
    def replaceAll(_a):
        _a = _a.replace('\\frac{d}{d t} \\phi{\\left(t \\right)}', '\\phiDerivadaEuler')
        _a = _a.replace('\\phi{\\left(t \\right)}', '\\phiEuler')
        _a = _a.replace('\\frac{d}{d t} \\theta{\\left(t \\right)}', '\\thetaDerivadaEuler')
        _a = _a.replace('\\theta{\\left(t \\right)}', '\\thetaEuler')
        _a = _a.replace('\\frac{d}{d t} \\psi{\\left(t \\right)}', '\\psiDerivadaEuler')
        _a = _a.replace('\\psi{\\left(t \\right)}', '\\psiEuler')
        _a = _a.replace('{\\left(\\psiEuler \\right)}', '\\psiEuler')
        _a = _a.replace('{\\left(\\thetaEuler \\right)}', '\\thetaEuler')
        _a = _a.replace('{\\left(\\phiEuler \\right)}', '\\phiEuler')
        _a = _a.replace('I_{x}', '\\Ix')
        _a = _a.replace('I_{y}', '\\Iy')
        _a = _a.replace('I_{z}', '\\Iz')
        _a = _a.replace('q_{1 }', '\\qi{1}')
        _a = _a.replace('q_{2 }', '\\qi{2}')
        _a = _a.replace('q_{3 }', '\\qi{3}')
        _a = _a.replace('q_{0 }', '\\qi{0}')
        _a = _a.replace('\\omega_{x}', '\\omegai{x}')
        _a = _a.replace('\\omega_{y}', '\\omegai{y}')
        _a = _a.replace('\\omega_{z}', '\\omegai{z}')
        _a = _a.replace('\\Delta t^{r}', '\\intervaloDiscretoRotacional')
        _a = _a.replace('\\alpha', '\\coeficineteSDCa')
        _a = _a.replace('\\beta', '\\coeficineteSDCb')
        return _a

    if printResult:
        print(replaceAll(sp.latex(_b)))
    else:
        return replaceAll(sp.latex(_b))