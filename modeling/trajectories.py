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
from common_functions import quaternion_functions as fq

plt.style.use("common_functions/roney3.mplstyle")
FIG_L = 6.29
FIG_A = (90.0) / 25.4


class Trajectory(object):
    """docstring for Trajectory."""

    def __init__(self, time_vector: np.ndarray, **kwargs):
        self.time_vector = time_vector
        self.n = time_vector.size
        # inertial states
        self.acceleration_vector_i = np.zeros((3, self.n))
        """ inertial acceleration """
        self.velocity_vector_i = np.zeros((3, self.n))
        """ inertial velocity """
        self.position_vector_i = np.zeros((3, self.n))
        """ inertial position """
        self.q_b_i = np.zeros((4, self.n))
        self.q_b_i[0, :] = 1.0
        """ quaternion transform inertial to body """
        # body states
        self.acceleration_vector_b = np.zeros((3, self.n))
        """ body acceleration """
        self.velocity_vector_b = np.zeros((3, self.n))
        """ body velocity """
        self.position_vector_b = np.zeros((3, self.n))
        """ body position """
        self.angular_velocity_vector_b = np.zeros((3, self.n))
        """ body velocity """
        self.angular_acceleration_vector_b = np.zeros((3, self.n))
        """ body acceleration """
        # find time indexes
        self.ft = 100.0 * 2.0 * np.pi
        self.ft_y = 50 * 2.0 * np.pi

        self.at = 150.0 / self.ft**2.0
        self.at_y = 150.0 / self.ft_y**2
        self.at_z = 10.0

        self.tc = 300.0

        self.ft_w_x = 10.0 * 2.0 * np.pi
        self.ft_w_y = -10.0 * 2.0 * np.pi
        self.ft_w_z = 20.0 * 2.0 * np.pi
        if kwargs["test"] == "translational_movement":
            print("Linear pure trajectory")
            self.aw_x = self.aw_y = self.aw_z = 0.0
        else:
            print("Linear and angular trajectory")
            self.aw_x = 10.0
            self.aw_y = 20.0
            self.aw_z = 30.0

        try:
            self.idx_y = np.where(self.time_vector > 2.0 * (2.0 * np.pi) / self.ft)[0][
                0
            ]
        except IndexError:
            self.idx_y = self.n
        try:
            self.idx_xy = np.where(self.time_vector > 3.0 * (2.0 * np.pi) / self.ft_y)[
                0
            ][0]
        except IndexError:
            self.idx_xy = self.n
        try:
            self.idx_xyz = np.where(self.time_vector > 5.0 * (2.0 * np.pi) / self.ft_y)[
                0
            ][0]
        except IndexError:
            self.idx_xyz = self.n
        for idx, t in enumerate(self.time_vector):
            self.acceleration_vector_i[:, idx] = self.acceleration_value(t)
            self.angular_velocity_vector_b[:, idx] = self.angular_velocity_value(t)
            self.angular_acceleration_vector_b[:, idx] = self.angular_acceleration_value(t)
        # set null values on first block of y and z
        self.acceleration_vector_i[1, : self.idx_y] *= 0.0
        # set null values on second block of
        self.acceleration_vector_i[0, self.idx_y : self.idx_xy] *= 0.0
        # integration with trap
        dt = self.time_vector[1] - self.time_vector[0]
        for idx in range(self.n - 1):
            self.velocity_vector_i[:, idx + 1] = self.velocity_vector_i[:, idx] + (
                self.acceleration_vector_i[:, idx]
                + self.acceleration_vector_i[:, idx + 1]
            ) * 0.5 * (self.time_vector[idx + 1] - self.time_vector[idx])
            self.position_vector_i[:, idx + 1] = self.position_vector_i[:, idx] + (
                self.velocity_vector_i[:, idx] + self.velocity_vector_i[:, idx + 1]
            ) * 0.5 * (self.time_vector[idx + 1] - self.time_vector[idx])
            self.q_b_i[:, idx + 1] = fq.mult_quat(
                p=self.q_b_i[:, idx],
                q=fq.expMap(self.angular_velocity_vector_b[:, idx], dt=dt),
            )
            self.acceleration_vector_b[:, idx + 1] = (
                fq.rotationMatrix(self.q_b_i[:, idx + 1])
                @ self.acceleration_vector_i[:, idx]
            )

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

    def angular_acceleration_value(self, time: float):
        y = np.array([0.0, 0.0, 0.0])
        # y[0] = aw_x*(np.exp(-time)-1) * np.sin(ft*time)
        # y[1] = aw_y*(np.exp(-time)-1) * np.cos(ft_y * time)
        # y[2] = aw_z*(np.exp(-time)-1)
        y[0] = self.ft_w_x * self.aw_x * np.cos(self.ft_w_x * time)
        y[1] = self.ft_w_y * self.aw_y * np.cos(self.ft_w_y * time)
        y[2] = self.ft_w_z * self.aw_z * np.cos(self.ft_w_z * time)
        return y

    def angular_velocity_value(self, time: float):
        y = np.array([0.0, 0.0, 0.0])
        # y[0] = aw_x*(np.exp(-time)-1) * np.sin(ft*time)
        # y[1] = aw_y*(np.exp(-time)-1) * np.cos(ft_y * time)
        # y[2] = aw_z*(np.exp(-time)-1)
        y[0] = self.aw_x * np.sin(self.ft_w_x * time)
        y[1] = self.aw_y * np.sin(self.ft_w_y * time)
        y[2] = self.aw_z * np.sin(self.ft_w_z * time)
        return y

    def plot_trajectories(self):
        self.plot_acceleration()
        self.plot_position()
        self.plot_velocity()
        self.plot_quaternion()
        self.plot_angular_velocity()

    def plot_acceleration(self):
        fig, ax = plt.subplots(3, 1, num=1, sharex=True, figsize=(FIG_L, FIG_A))

        fig.supylabel("aceleração")
        fig.supxlabel(r"Tempo, [$\unit{\ms}$]")
        ax[0].plot(
            1e3 * self.time_vector,
            self.acceleration_vector_i[0, :],
            label=r"Sistema $\mathcal{I}$",
        )
        ax[1].plot(
            1e3 * self.time_vector,
            self.acceleration_vector_i[1, :],
            label=r"Sistema $\mathcal{I}$",
        )
        ax[2].plot(
            1e3 * self.time_vector,
            self.acceleration_vector_i[2, :],
            label=r"Sistema $\mathcal{I}$",
        )
        ax[0].plot(
            1e3 * self.time_vector,
            self.acceleration_vector_b[0, :],
            label=r"Sistema $\mathcal{B}$",
        )
        ax[1].plot(
            1e3 * self.time_vector,
            self.acceleration_vector_b[1, :],
            label=r"Sistema $\mathcal{B}$",
        )
        ax[2].plot(
            1e3 * self.time_vector,
            self.acceleration_vector_b[2, :],
            label=r"Sistema $\mathcal{B}$",
        )
        ax[2].legend()
        plt.show()

    def plot_velocity(self):
        fig, ax = plt.subplots(3, 1, num=2, sharex=True, figsize=(FIG_L, FIG_A))

        fig.supylabel("velocidade")
        fig.supxlabel(r"Tempo, [$\unit{\ms}$]")
        ax[0].plot(1e3 * self.time_vector, self.velocity_vector_i[0, :])
        ax[1].plot(1e3 * self.time_vector, self.velocity_vector_i[1, :])
        ax[2].plot(1e3 * self.time_vector, self.velocity_vector_i[2, :])

        plt.show()

    def plot_position(self):
        fig, ax = plt.subplots(3, 1, num=3, sharex=True, figsize=(FIG_L, FIG_A))

        fig.supylabel("posição")
        fig.supxlabel(r"Tempo, [$\unit{\ms}$]")
        ax[0].plot(1e3 * self.time_vector, self.position_vector_i[0, :])
        ax[1].plot(1e3 * self.time_vector, self.position_vector_i[1, :])
        ax[2].plot(1e3 * self.time_vector, self.position_vector_i[2, :])

        plt.show()

    def plot_angular_velocity(self):
        fig, ax = plt.subplots(3, 1, num=4, sharex=True, figsize=(FIG_L, FIG_A))

        fig.supylabel("velocidade angular")
        fig.supxlabel(r"Tempo, [$\unit{\ms}$]")
        ax[0].plot(1e3 * self.time_vector, self.angular_velocity_vector_b[0, :])
        ax[1].plot(1e3 * self.time_vector, self.angular_velocity_vector_b[1, :])
        ax[2].plot(1e3 * self.time_vector, self.angular_velocity_vector_b[2, :])

    def plot_angular_acceleration(self):
        fig, ax = plt.subplots(3, 1, num=5, sharex=True, figsize=(FIG_L, FIG_A))

        fig.supylabel(r"$\dot{\mathbf{\omega} [\unit{\radian\per\squared\second}]$")
        fig.supxlabel(r"Tempo, [$\unit{\ms}$]")
        ax[0].plot(1e3 * self.time_vector, self.angular_angular_vector_b[0, :])
        ax[1].plot(1e3 * self.time_vector, self.angular_angular_vector_b[1, :])
        ax[2].plot(1e3 * self.time_vector, self.angular_angular_vector_b[2, :])

    def plot_quaternion(self):
        fig, ax = plt.subplots(4, 1, num=6, sharex=True, figsize=(FIG_L, FIG_A))
        fig.supxlabel(r"Tempo, $[\unit{\ms}]$")
        fig.supylabel("quatérnion")
        ax[0].plot(1e3 * self.time_vector, self.q_b_i[0, :])
        ax[1].plot(1e3 * self.time_vector, self.q_b_i[1, :])
        ax[2].plot(1e3 * self.time_vector, self.q_b_i[2, :])
        ax[3].plot(1e3 * self.time_vector, self.q_b_i[3, :])
