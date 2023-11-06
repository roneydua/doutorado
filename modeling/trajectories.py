#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   Trajectories.py
@Time    :   2023/11/05 23:55:44
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt

class Trajectory(object):
    """docstring for Trajectory."""

    def __init__(self, time_vector: np.ndarray, **kwargs):
        self.n = time_vector.size
        if kwargs["test"] == "linear":
            print("Linear pure trajectory")
            self.ft = 1.0 * 2.0 * np.pi
            self.ft_y = .50 * 2.0 * np.pi

            self.at = 150.0 / self.ft**2.0
            self.at_y = 150.0 / self.ft_y**2
            self.at_z = 10.0

            self.tc = 3000.0

            self.ft_w_x = 10.0 * 2.0 * np.pi
            self.ft_w_y = -1.0 * 2.0 * np.pi
            self.ft_w_z = -1.0 * 2.0 * np.pi

            self.aw_x = 2.0
            self.aw_y = 2.0
            self.aw_z = 3.0
            self.time_vector = time_vector
            self.acceleration_vector = np.zeros((3, self.n))
            self.velocity_vector = np.zeros((3, self.n))
            self.position_vector = np.zeros((3, self.n))
            self.angular_velocity_vector = np.zeros((3, self.n))
            # find time indexes
            try:
                self.idx_y = np.where(self.time_vector > 2.0 *(2.0*np.pi)/ self.ft)[0][0]
            except IndexError:
                self.idx_y = self.n
            try:
                self.idx_xy = np.where(self.time_vector > 3.0 *(2.0*np.pi)/ self.ft_y)[0][0]
            except IndexError:
                self.idx_xy = self.n
            try:
                self.idx_xyz = np.where(self.time_vector > 5.0 *(2.0*np.pi)/ self.ft_y)[0][0]
            except IndexError:
                self.idx_xyz = self.n
            for idx, t in enumerate(self.time_vector):
                self.acceleration_vector[:, idx] = self.acceleration_value(t)
                # self.angular_velocity_vector[:,idx] = self.angular_velocity_value(t)
            # set null values on firts block of y and z
            self.acceleration_vector[1:, : self.idx_y] *= 0.0
            # set null values on second block of
            self.acceleration_vector[0, self.idx_y : self.idx_xy] *= 0.0
            self.acceleration_vector[2, self.idx_y : self.idx_xy] *= 0.0
            # set null values on third block
            self.acceleration_vector[2, self.idx_xy : self.idx_xyz] *= 0.0
            # integration with trap
            for idx in range(self.n-1):
                self.velocity_vector[:,idx+1] = self.velocity_vector[:,idx]+(self.acceleration_vector[:,idx]+self.acceleration_vector[:,idx+1])*0.5*(self.time_vector[idx+1]-self.time_vector[idx])
                self.position_vector[:,idx+1] = self.position_vector[:,idx]+(self.velocity_vector[:,idx]+self.velocity_vector[:,idx+1])*0.5*(self.time_vector[idx+1]-self.time_vector[idx])

    def acceleration_value(self, time: float):
        y = np.array([0.0, 0.0, 0.0])
        y[0] = -self.at * (self.ft**2) * np.sin(self.ft * time)
        y[1] = -self.at_y * (self.ft_y**2) * np.sin(self.ft_y * time)
        y[2] = self.at_z * (1 - np.exp(-self.tc * time))
        return y

    def velocity_value(self, time: float):
        y = np.array([0.0, 0.0, 0.0])
        y[0] = self.at * self.ft * np.cos(self.ft * time)
        y[1] = self.at_y * self.ft_y * np.cos(self.ft_y * time)
        y[2] = self.at_z * (time + np.exp(-self.tc * time) / self.tc)
        return y

    def position_value(self, time: float):
        y = np.array([0.0, 0.0, 0.0])
        y[0] = self.at * np.sin(self.ft * time)
        y[1] = self.at_y * np.sin(self.ft_y * time)
        y[2] = self.at_z * (0.5 * time**2 - np.exp(-self.tc * time) / (self.tc**2))
        return y

    def angular_velocity_value(self, time: float):
        y = np.array([0.0, 0.0, 0.0])
        # y[0] = aw_x*(np.exp(-time)-1) * np.sin(ft*time)
        # y[1] = aw_y*(np.exp(-time)-1) * np.cos(ft_y * time)
        # y[2] = aw_z*(np.exp(-time)-1)
        y[0] = self.aw_x * np.sin(self.ft * time)
        y[1] = self.aw_y * np.sin(self.ft_y * time)
        y[2] = self.aw_z * np.sin(self.ft_y * time)
        return y

    def plot_trajectories(self):
        plt.figure()
        plt.ylabel('aceleração')
        plt.plot(self.time_vector,self.acceleration_vector[0,:],'.')
        plt.plot(self.time_vector,self.acceleration_vector[1,:],'.')
        plt.plot(self.time_vector,self.acceleration_vector[2,:],'.')
        plt.show()
        plt.figure()
        plt.ylabel('velocidade')
        plt.plot(self.time_vector,self.velocity_vector[0,:],'.')
        plt.plot(self.time_vector,self.velocity_vector[1,:],'.')
        plt.plot(self.time_vector,self.velocity_vector[2,:],'.')
        plt.show()
        plt.figure()
        plt.ylabel('posição')
        plt.plot(self.time_vector,self.position_vector[0,:],'.')
        plt.plot(self.time_vector,self.position_vector[1,:],'.')
        plt.plot(self.time_vector,self.position_vector[2,:],'.')
        plt.show()
        plt.figure()
        plt.ylabel('posição')
        plt.plot(self.time_vector,self.angular_velocity_vector[0,:],'.')
        plt.plot(self.time_vector,self.angular_velocity_vector[1,:],'.')
        plt.plot(self.time_vector,self.angular_velocity_vector[2,:],'.')