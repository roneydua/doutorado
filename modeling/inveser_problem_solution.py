#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   inveser_problem_solution.py
@Time    :   2023/10/19 15:15:52
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
'''
from modeling.mathModelAccel import InverseProblem
import locale

import h5py
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
locale.setlocale(locale.LC_ALL, "pt_BR.UTF-8")
# plt.style.use("default")
plt.style.use("common_functions/roney3.mplstyle")
my_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FIG_L = 6.29
FIG_A = (90.0) / 25.4


def recover_acceleration():
    # accel = AccelModelInertialFrame()
    fibers_with_length_info = np.array([1, 5, 9, 11])
    deformation = np.zeros((fibers_with_length_info.size,1))
    ip = InverseProblem(fibers_with_length_info)
    f = h5py.File('modeling_data.hdf5', 'a')
    ff = f['test_of_group']
    ## Solver inverse problem
    temp_fiber_l = np.zeros((fibers_with_length_info.size, ff['t'].size))
    for idx, _f in enumerate(fibers_with_length_info):
        temp_fiber_l[idx,:] = ff['fiber_len'][_f-1,:]
    
    for idx_t,_t in enumerate(tqdm(ff['t'])):
        ip.compute_inverse_problem_solution(temp_fiber_l[:,idx_t])
        ip.estimate_f_vector()
        ff['recover_accel_simple'][:,idx_t] = ip.estimate_ddrm_B()
        
    f.close()


recover_acceleration()


