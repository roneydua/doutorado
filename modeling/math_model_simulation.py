#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File    :   math_model_simulation.py
@Time    :   2023/10/24 09:04:05
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""

import numpy as np
from tqdm import tqdm
import h5py
import pandas as pd
from modeling.mathModelAccel import (
    AccelModelInertialFrame,
    InverseProblem,
    SimpleSolution,
)
from modeling.trajectories import Trajectory
from common_functions.RK import RungeKutta
from common_functions import quaternion_functions as fq

import matplotlib.pyplot as plt


class statesOfSimulation_object(object):
    def __init__(
        self, dt=1e-5, tf=1.0, hf=h5py.Group, accel=AccelModelInertialFrame()
    ) -> None:
        """
        __init__ states of simulation collect all states too simulate
        Args:
            dt: Time step for integration. Defaults to 1e-5.NOTE: With dt<1e-5 the simulation show artificia damper.
            tf: Time of simulation. Defaults to 1.0.
        """
        self.hf = hf
        self.hf.attrs["dt"] = dt
        """ Step of time of simulation."""
        self.hf["t"] = np.arange(0.0, tf, dt)
        self.hf["t"].attrs["about"] = "Time of simulation in seconds"

        self.hf["x"] = np.zeros((26, self.hf["t"].size))
        self.hf["x"].attrs[
            "about"
        ] = "States  d_rb[:3] d_rm[3:6] rb[6:9] rm[9:12] qb[12:16] qm[16:20] wb[20:23] wm[23:26]"
        self.hf["true_accel"] = np.zeros(shape=(3, self.hf["t"].size))
        self.hf["true_accel"].attrs["about"] = "Exact acceleration"
        """ Relative deformation with respect of the initial length (l-l0)/l. """
        self.hf["f"] = np.zeros(shape=(12, 3, self.hf["t"].size))
        self.hf["fiber_len"] = np.zeros(shape=(12, self.hf["t"].size))

        self.hf.attrs["b_B"] = accel.b_B
        self.hf.attrs["m_M"] = accel.m_M
        self.hf.attrs["k"] = accel.k
        self.hf.attrs["seismic_mass"] = accel.seismic_mass
        self.hf.attrs["fiber_length"] = accel.fiber_length


if __name__ == "__main__":
    accel = AccelModelInertialFrame(
        damper_for_computation_simulations=0.0, fiber_length=3e-3
    )
    f = h5py.File("modeling_data_4.hdf5", "w")
    ff = f.require_group("pure_translational_movement")
    s = statesOfSimulation_object(tf=0.3, dt=5e-5, hf=ff, accel=accel)
    traj = Trajectory(s.hf["t"][:], test="linear")
    # traj.plot_trajectories()

    RK = RungeKutta(s.hf["x"].shape[0], accel.dd_x_forced_body_state)
    # initial conditions for all quaternions as [1, 0, 0, 0]
    s.hf["x"][12, :] = 1.0
    s.hf["x"][16, :] = 1.0
    # initial misalignment
    s.hf["x"][12:16, 0] = fq.eulerQuaternion(yaw=0, pitch=0, roll=0)  # qb
    s.hf["x"][16:20, 0] = fq.eulerQuaternion(yaw=0, pitch=0, roll=0)  # qm
    # Set body sensor and seismic mass initial conditions
    s.hf["x"][:3, 0] = traj.velocity_vector_i[:, 0]  # body velocity
    s.hf["x"][3:6, 0] = traj.velocity_vector_i[:, 0]  # seismic mass  velocity
    s.hf["x"][6:9, 0] = traj.position_vector_i[:, 0]  # body position
    s.hf["x"][9:12, 0] = traj.position_vector_i[:, 0]  # seismic mass position
    s.hf["x"][20:23, 0] = traj.angular_velocity_vector[:, 0]
    s.hf["x"][23:, 0] = traj.angular_velocity_vector[:, 0]

    accel.update_states(
        rb=s.hf["x"][6:9, 0],
        rm=s.hf["x"][9:12, 0],
        qb=s.hf["x"][12:16, 0],
        qm=s.hf["x"][16:20, 0],
    )
    s.hf["f"][:, :, 0] = accel.f
    s.hf["fiber_len"][:, 0] = np.linalg.norm(accel.f, axis=1)

    fibers_with_length_info = np.array([1, 4, 5, 8, 9, 12])
    # deformation = np.zeros((fibers_with_length_info.size, 1))
    ip = InverseProblem(fibers_with_length_info, fiber_length=accel.fiber_length)
    ss = SimpleSolution(
        np.array([1, 5, 9]), fiber_length=accel.fiber_length, push_pull=True
    )
    for i in tqdm(range(s.hf["t"].size - 1)):
        # integration of states
        s.hf["x"][:, i + 1] = RK.integrates_states(
            s.hf["x"][:, i], u=None, h=s.hf.attrs["dt"]
        )
        # Put the body in specific trajectories
        s.hf["x"][:3, i + 1] = traj.velocity_vector_i[:, i + 1]
        s.hf["x"][6:9, i + 1] = traj.position_vector_i[:, i + 1]
        s.hf["x"][20:23, i + 1] = traj.angular_velocity_vector[:, i + 1]

        s.hf["x"][12:16, i + 1] = fq.mult_quat(
            p=s.hf["x"][12:16, i],
            q=fq.expMap(s.hf["x"][20:23, i + 1], dt=s.hf.attrs["dt"]),
        )

        s.hf["true_accel"][:, i + 1] = traj.acceleration_vector_i[:, i + 1]

        # # update class structs to compute f vectors
        accel.update_states(
            rb=s.hf["x"][6:9, i + 1],
            rm=s.hf["x"][9:12, i + 1],
            qb=s.hf["x"][12:16, i + 1],
            qm=s.hf["x"][16:20, i + 1],
        )
        ## take the f vector of integration
        s.hf["f"][:, :, i + 1] = accel.f
        s.hf["fiber_len"][:, i + 1] = np.linalg.norm(accel.f, axis=1)

        
    f.close()
