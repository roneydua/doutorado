#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   modeling_2d_case.py
@Time    :   2023/10/24 14:05:18
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''

import numpy as np
import locale
import matplotlib.pyplot as plt
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4
# %%

class Accel2D(object):
    """docstring for Accel2D."""
    def __init__(self, fibers_with_info):
        self.m_M = np.array([[-1,0], [0,-1], [1,0],[0,1]])
        self.b_B = np.array([[-2,0], [0,-2], [2,0],[0,2]])
        self.estimated_f_B = self.m_M-self.b_B
        self.fibers_with_info=fibers_with_info
        self.fibers_with_info_index = fibers_with_info-1
        self.var_xi = np.ones((self.fibers_with_info.size, 3))
        self.var_gamma = np.zeros(3)
        self.var_psi = np.zeros(self.fibers_with_info.size)
        _aux_vector = np.ones((self.fibers_with_info.size, 2))
        self.f = np.ones((self.fibers_with_info.size, 2))
        self.fiber_len = np.zeros(self.fibers_with_info.size)
        '''auxiliar vector to compute (m-b) with dimenstion fiber_with_sise by 3, used on var_xi and var_psi'''
        self.aux_var_psi_matrix = np.zeros(self.fibers_with_info.size)
        for i, j in enumerate(self.fibers_with_info_index):
            _aux_vector[i, :] = self.m_M[j, :] - self.b_B[j, :]
            self.aux_var_psi_matrix[i] = _aux_vector[i, :].dot(
                _aux_vector[i, :])
        self.var_xi[:, 1:] = 2.0*_aux_vector
        if self.var_xi.shape[0] == self.var_gamma.size:
            # in this case the least squared method use only the inverse matrix of var_gamma
            self.least_square_matrix = np.linalg.inv(self.var_xi)
        else:
            # It is necessary compute pseud inverse of matrix
            self.least_square_matrix = np.linalg.pinv(self.var_xi)
            
        self.diff_m_M_b_B = self.m_M-self.b_B
        
    def compute_f(self, _rb):
        for i, j in enumerate(self.fibers_with_info_index):
            self.f[i,:] = _rb+self.diff_m_M_b_B[j,:]
            self.fiber_len[i] = np.linalg.norm(self.f[i, :])
    def compute_inverse_problem_solution(self, fiber_len: np.ndarray):
        '''
        compute_inverse_problem_solution 
        Args:
            fiber_len: vector of fiber_len is ((f).dot(f))^2
        '''
        self.var_psi = np.square(
            fiber_len)-self.aux_var_psi_matrix
        self.var_gamma = self.least_square_matrix @ self.var_psi

    
accel= Accel2D(np.array([1,2,3]))
accel.compute_f(np.array([-5.0e-8,5.0e-8]))
accel.fiber_len[0]
accel.compute_inverse_problem_solution(accel.fiber_len)
accel.var_gamma[1:]
