import locale

import dill
import numpy as np
from tqdm import tqdm

import common_functions.funcoesQuaternion as fq
from mathModelAccel import AccelModel
from RK import RungeKutta


class statesOfSimulation_object(object):

    def __init__(self, dt=1e-5, tf=1.0) -> None:
        '''
        __init__ states of simulation collect all states too simulate 
        Args:
            dt: Time step for integration. Defaults to 1e-5.NOTE: With dt<1e-5 the simulation show artificia damper.
            tf: Time of simulation. Defaults to 1.0.
        '''
        self.dt = dt
        ''' Step of time of simulation.
        NOTE: For step minnor that 1e-5  '''
        self.t = np.arange(0.0, tf, dt)
        ''' Time of simulation'''
        self.x = np.zeros((26, self.t.size))
        ''' States  d_rb[:3] d_rm[3:6] rb[6:9] rm[9:12] qb[12:16] qm[16:20] wb[20:23] wm[23:26]'''
        # if 0:
        self.f_norm = np.zeros(shape=(12, self.t.size - 1))
        ''' Vector of length of fibers. This vector differs of initial length.'''
        self.u = np.zeros(shape=(6, self.t.size))
        '''translation and rotation generalized forces.'''
        self.true_accel = np.zeros(shape=(3, self.t.size))
        """ Exact acceleration"""
        self.recover_accel_simple = np.zeros((3, self.t.size))
        ''' Recover acceleration with simple model with no without consideration of cross effects (As in the doctoral dissertation of the Cazo) '''
        self.recover_accel_with_no_relatives_rotations = np.zeros(
            (3, self.t.size))
        """ Recover acceleration without considering relative rotations."""
        self.deformation = 0.0 * self.f_norm
        """ Relative deformation with respect of the initial length (l-l0)/l. """

        self.accel_instance = AccelModel()
        """ Instance used on simulation"""

    def saveData(self, filename):
        """Save data with Dill Package"""
        with open(filename, "wb") as fp:
            dill.dump(self, fp)
        print("data save as " + filename)


if __name__ == "__main__":
    ft = 2.0 * np.pi
    ft_y = 2.0 * np.pi
    at = 0
    at_y = 0
    at_z = 9.89
    tc = 3.0
    ft_w_x = 1
    ft_w_y = -1
    ft_w_z = -1

    def acceleration_value(y, time):
        y[0] = -at * (ft**2) * np.sin(ft * time)
        y[1] = -at_y * (ft_y**2) * np.sin(ft_y * time)
        y[2] = at_z * (1 - np.exp(-tc * time))

    def velocity_value(y, time):
        y[0] = at * ft * np.cos(ft * time)
        y[1] = at_y * ft_y * np.cos(ft_y * time)
        y[2] = at_z * (time + np.exp(-tc * time) / tc)

    def position_value(y, time):
        y[0] = at * np.sin(ft * time)
        y[1] = at_y * np.sin(ft_y * time)
        y[2] = at_z * (0.5 * time**2 - np.exp(-tc * time) / (tc**2))

    def angular_velocity_value(y, time):
        y[0] = 0.0  #time+ np.exp(-time)-1#* np.sin(ft*time)
        y[1] = 0.0  #time+np.exp(-time)-1#* np.cos(ft_y * time)
        y[2] = 0.0  #time+ np.exp(-time)-1

    accel = AccelModel(fibers_with_info=np.array([1, 5, 9, 12], dtype=np.int8),
                       inverse_problem_full=False,
                       damper_for_computation_simulations=1e-6)
    s = statesOfSimulation_object(tf=2.0, dt=1e-5)
    # save instance accel used on model
    s.accel_instance = accel
    RK = RungeKutta(s.x.shape[0], accel.dd_x_forced_body_state)
    # initial conditions for all quaternions as [1, 0, 0, 0]
    s.x[12, :] = 1.0
    s.x[16, :] = 1.0
    # initial misalignment
    s.x[12:16, 0] = fq.eulerQuaternion(yaw=0, pitch=0, roll=0)  # qb
    s.x[16:20, 0] = fq.eulerQuaternion(yaw=0, pitch=0, roll=0)  # qm
    # Set body sensor and seismic mass initial conditions
    velocity_value(s.x[:3, 0], s.t[0])  # body velocity
    velocity_value(s.x[3:6, 0], s.t[0])  # seismic mass velocity
    position_value(s.x[6:9, 0], s.t[0])  # body position
    position_value(s.x[9:12, 0], s.t[0])  # seismic mass position

    angular_velocity_value(s.x[20:23, 0], s.t[0])
    angular_velocity_value(s.x[23:, 0], s.t[0])

    accel.update_states(rb=s.x[6:9, 0],
                        rm=s.x[9:12, 0],
                        qb=s.x[12:16, 0],
                        qm=s.x[16:20, 0])
    for i in tqdm(range(s.t.size - 1)):
        # integration of states
        s.x[:, i + 1] = RK.integrates_states(s.x[:, i], s.u[:, i], s.dt)
        # Put the body in specific trajectories.
        velocity_value(s.x[:3, i + 1], s.t[i + 1])
        position_value(s.x[6:9, i + 1], s.t[i + 1])
        angular_velocity_value(s.x[20:23, i + 1], s.t[i + 1])

        fq.MultQuat(s.x[12:16, i + 1], s.x[12:16, i],
                    fq.expMap(s.x[20:23, i + 1], dt=s.dt))
        acceleration_value(s.true_accel[:, i + 1], s.t[i + 1])
        # update class structs to compute f vectors
        accel.update_states(rb=s.x[6:9, i + 1],
                            rm=s.x[9:12, i + 1],
                            qb=s.x[12:16, i + 1],
                            qm=s.x[16:20, i + 1])
        s.f_norm[:, i] = np.linalg.norm(accel.f, axis=1) - accel.fiber_length
        s.deformation[:, i] = s.f_norm[:, i] / accel.fiber_length
        s.recover_accel_simple[
            0, i + 1] = 4.0 * accel.k * s.f_norm[0, i] / accel.seismic_mass
        s.recover_accel_simple[
            1, i + 1] = 4.0 * accel.k * s.f_norm[4, i] / accel.seismic_mass
        s.recover_accel_simple[
            2, i +
            1] = 4.0 * accel.k * s.f_norm[9, i] / accel.seismic_mass + accel.G

    s.saveData('data/vertical_acceleration.pickle')


