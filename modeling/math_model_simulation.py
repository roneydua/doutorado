import dill
import numpy as np
from tqdm import tqdm
import h5py
import pandas as pd
from modeling.mathModelAccel import AccelModelInertialFrame
from common_functions.RK import RungeKutta
from common_functions import funcoesQuaternion as fq

import matplotlib.pyplot as plt


class statesOfSimulation_object(object):

    def __init__(self, dt=1e-5, tf=1.0, hf=h5py.Group, accel=AccelModelInertialFrame()) -> None:
        '''
        __init__ states of simulation collect all states too simulate 
        Args:
            dt: Time step for integration. Defaults to 1e-5.NOTE: With dt<1e-5 the simulation show artificia damper.
            tf: Time of simulation. Defaults to 1.0.
        '''
        self.hf = hf
        self.hf.attrs['dt'] = dt
        ''' Step of time of simulation.'''
        self.hf['t'] = np.arange(0.0, tf, dt)
        self.hf['t'].attrs['about'] = 'Time of simulation in seconds'

        self.hf['x'] = np.zeros((26, self.hf['t'].size))
        self.hf['x'].attrs['about'] = 'States  d_rb[:3] d_rm[3:6] rb[6:9] rm[9:12] qb[12:16] qm[16:20] wb[20:23] wm[23:26]'
        # # if 0:
        # self.hf['f_norm'] = np.zeros(shape=(12, self.hf['t'].size - 1))
        # self.hf['f_norm'].attrs['about'] = 'Vector of length of fibers. This vector differs of initial length.'
        self.hf['u'] = np.zeros(shape=(6, self.hf['t'].size))
        self.hf['u'].attrs['about'] = 'translation and rotation generalized forces.'
        self.hf['true_accel'] = np.zeros(shape=(3, self.hf['t'].size))
        self.hf['true_accel'].attrs['about'] = 'Exact acceleration'
        """ Relative deformation with respect of the initial length (l-l0)/l. """

        self.hf.attrs['b_B'] = accel.b_B
        self.hf.attrs['m_M'] = accel.m_M
        self.hf.attrs['k'] = accel.k
        self.hf.attrs['seismic_mass'] = accel.seismic_mass
        self.hf.attrs['fiber_length'] = accel.fiber_length

        # """ Instance used on simulation"""

    # def saveData(self, filename):
    #     """Save data with Dill Package"""
    #     with open(filename, "wb") as fp:
    #         dill.dump(self, fp)
    #     print("data save as " + filename)


if __name__ == "__main__":
    ft = 5.0 * np.pi
    ft_y = 2.0 * np.pi
    at = 1.0
    at_y = 1.0
    at_z = 9.89
    tc = 3.0
    ft_w_x = 1.0
    ft_w_y = -1.0
    ft_w_z = -1.0

    def acceleration_value(time):
        y = np.array([0.0, 0.0, 0.0])
        y[0] = -at * (ft**2) * np.sin(ft * time)
        y[1] = -at_y * (ft_y**2) * np.sin(ft_y * time)
        y[2] = at_z * (1 - np.exp(-tc * time))
        return y

    def velocity_value(time):
        y = np.array([0.0, 0.0, 0.0])
        y[0] = at * ft * np.cos(ft * time)
        y[1] = at_y * ft_y * np.cos(ft_y * time)
        y[2] = at_z * (time + np.exp(-tc * time) / tc)
        return y

    def position_value(time):
        y = np.array([0.0, 0.0, 0.0])
        y[0] = at * np.sin(ft * time)
        y[1] = at_y * np.sin(ft_y * time)
        y[2] = at_z * (0.5 * time**2 - np.exp(-tc * time) / (tc**2))
        return y

    def angular_velocity_value(time):
        y = np.array([0.0, 0.0, 0.0])
        y[0] =  (np.exp(-time)-1)* np.sin(ft*time)
        y[1] =  (np.exp(-time)-1)* np.cos(ft_y * time)
        y[2] =   np.exp(-time)-1
        return y

    accel = AccelModelInertialFrame(damper_for_computation_simulations=1e-6)
    f = h5py.File('teste.hdf5', 'w')
    ff = f.create_group('test_of_group')
    s = statesOfSimulation_object(tf=1.0, dt=5e-6, hf=ff, accel=accel)
    # save instance accel used on model

    RK = RungeKutta(s.hf['x'].shape[0], accel.dd_x_forced_body_state)
    # initial conditions for all quaternions as [1, 0, 0, 0]
    s.hf['x'][12, :] = 1.0
    s.hf['x'][16, :] = 1.0
    # initial misalignment
    s.hf['x'][12:16, 0] = fq.eulerQuaternion(yaw=0, pitch=0, roll=0)  # qb
    s.hf['x'][16:20, 0] = fq.eulerQuaternion(yaw=0, pitch=0, roll=0)  # qm
    # Set body sensor and seismic mass initial conditions
    s.hf['x'][:3, 0] = velocity_value(s.hf['t'][0])  # body velocity
    s.hf['x'][3:6, 0] = velocity_value(s.hf['t'][0])  # seismic mass  velocity
    s.hf['x'][6:9, 0] = position_value(s.hf['t'][0])  # body position
    s.hf['x'][9:12, 0] = position_value(s.hf['t'][0])# seismic mass position

    s.hf['x'][20:23, 0] = angular_velocity_value(s.hf['t'][0])
    s.hf['x'][23:, 0] = angular_velocity_value(s.hf['t'][0])

    accel.update_states(rb=s.hf['x'][6:9, 0],
                        rm=s.hf['x'][9:12, 0],
                        qb=s.hf['x'][12:16, 0],
                        qm=s.hf['x'][16:20, 0])
    for i in tqdm(range(s.hf['t'].size - 1)):
        # integration of states
        s.hf['x'][:, i + 1] = RK.integrates_states(
            s.hf['x'][:, i], s.hf['u'][:, i], s.hf.attrs['dt'])
        # print(s.hf['x'][:, i + 1])
        # Put the body in specific trajectories
        s.hf['x'][:3, i + 1] = velocity_value(s.hf['t'][i + 1])
        s.hf['x'][6:9, i + 1] = position_value(s.hf['t'][i + 1])
        s.hf['x'][20:23, i + 1] = angular_velocity_value(s.hf['t'][i + 1])

        s.hf['x'][12:16, i + 1] = fq.mult_quat(p=s.hf['x'][12:16, i], q=fq.expMap(
            s.hf['x'][20:23, i + 1], dt=s.hf.attrs['dt']))

        s.hf['true_accel'][:, i + 1] = acceleration_value(s.hf['t'][i + 1])
        # print(s.hf['true_accel'][:, i + 1])
        # # update class structs to compute f vectors
        accel.update_states(rb=s.hf['x'][6:9, i + 1],
                            rm=s.hf['x'][9:12, i + 1],
                            qb=s.hf['x'][12:16, i + 1],
                            qm=s.hf['x'][16:20, i + 1])
        # s.f_norm[:, i] = np.linalg.norm(accel.f_B, axis=1) - accel.fiber_length
        # s.deformation[:, i] = s.f_norm[:, i] / accel.fiber_length
        # s.recover_accel_simple[
        #     0, i + 1] = 4.0 * accel.k * s.f_norm[0, i] / accel.seismic_mass
        # s.recover_accel_simple[
        #     1, i + 1] = 4.0 * accel.k * s.f_norm[4, i] / accel.seismic_mass
        # s.recover_accel_simple[
        #     2, i +
        #     1] = 4.0 * accel.k * s.f_norm[9, i] / accel.seismic_mass + accel.G
    # f.close()
    # s.saveData('data/vertical_acceleration.pickle')

    f.close()
