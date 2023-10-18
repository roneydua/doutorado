import numpy as np
from numpy.linalg import inv
import common_functions.funcoesQuaternion as fq


# Constates coordinates of connections


class states:
    rot = np.eye(3)
    ''' Rotation matrix'''
    euler = np.zeros(3)
    '''Euler angules. psi, theta, phi [deg]'''

    def __init__(self):

        self.q = np.array([1., 0., 0., 0.])
        '''attitude quaternion'''
        self.w = np.array([0, 0, 0])
        '''# angular body velcity'''
        self.r = np.array([0, 0, 0])
        '''inertial position'''
        self.dr = np.array([0, 0, 0])
        '''inertial velocity'''
        self.ddr = np.array([0, 0, 0])
        '''inertial acceleration'''

    def updates_attitude(self, q):
        '''
        updates_attitude Update states of attitude quaterion q, Euler angles and matrix rotation
        Args:
            q: Attitude quaterion.
        '''
        self.q = q
        self.euler = fq.quat2Euler(q, deg=1)
        self.rot = fq.rotationMatrix(self.q)

    def psi(self):
        return self.euler[0]

    def theta(self):
        return self.euler[1]

    def phi(self):
        return self.euler[2]


class AccelModel(object):

    density = 2.6989e3
    '''material density  (aluminiun)'''
    E = 70e9
    '''Young's Module in  GPa'''
    fiber_diameter = 125e-6
    '''fiber diameter in meters'''
    sms = states()
    '''Seismic mass states'''
    bss = states()
    '''Sensor base states'''
    G = 0.0  # -9.89
    '''Gravity in meters per squared second'''
    damper_for_computation_simulations = 0.0
    '''Artificial damper for numérical simulations'''
    
    rm_B =  np.array([0., 0., 0.])
    q_M_B = np.array([1., 0., 0., 0.])
    '''Attitude quaternion of sensor base with respect of seismic mass `{}^{B}_{M}q`'''
    rot_M_B = np.eye(3)
    '''Attitude matrix rotation of sensor base with respect of seismic mass `{}^{B}_{M}R`'''

    def __init__(self, seismic_edge=16.4e-3,
                 fiber_diameter=125e-6, fiber_length=3e-3,  damper_for_computation_simulations=0.0):
        '''
        __init__ Contructor of mathModel

        _extended_summary_

        Args:
            seismic_edge: Seismic mass cube edge. Defaults to 16.4e-3.
            fiber_diameter: diameter of optical fiber in meters. Defaults to 125e-6.
            fiber_length: fiber size in meters. Defaults to 3e-3.
            fibers_with_info: Number of fibers with information of the size. Defaults to np.arange(1, 13).
            damper_for_computation_simulations: Artificial damper coeficient. Defaults to 0.0.
        '''
        self.damper_for_computation_simulations = damper_for_computation_simulations
        self.fibers_with_info = fibers_with_info
        '''index of fiber with deformation informartion. This is important to recover aceleration'''
        self.fiber_diameter = fiber_diameter
        '''Fiber diameter in meters'''
        self.fiber_length = fiber_length
        '''initial fiber length'''
        self.k = (self.E * 0.25 * np.pi *
                  (self.fiber_diameter**2)) / self.fiber_length
        '''stifness of optical fiber'''
        self.seismic_edge = seismic_edge
        '''Length of seismic edge'''
        self.seismic_mass = (self.seismic_edge ** 3.0) * self.density
        '''Sismic mass'''
        self.external_base_sensor_edge = self.seismic_edge + 2 * self.fiber_length + 4e-3
        ''' Approximation of the base of the sensor as a cube.The value 6e-3 refers double the length of the fibers that support the cube.4e-3 is twice the thickness.'''
        self.base_sensor_edge = self.seismic_edge + 2 * self.fiber_length
        self.base_sensor_mass = (
            self.external_base_sensor_edge**3 - self.base_sensor_edge ** 3) * self.density
        self.inertial_seismic_mass = np.eye(
            3) * self.seismic_mass / 6.0 * self.seismic_edge**2
        self.inertial_base_sensor = np.eye(
            3) * self.base_sensor_mass / 6.0 * self.seismic_edge**2
        # NOTE: Only for computational economy
        self.inertial_base_sensor_inv = inv(self.inertial_base_sensor)
        self.inertial_seismic_mass_inv = inv(self.inertial_seismic_mass)
        # perpendicular distance of center
        _d = self.seismic_edge / 2 - 1e-3
        _e = self.seismic_edge / 2
        self.m_M = np.array([[_e, _d, 0.0],  # X+ grade
                             [_e, -_d, 0.0],
                             [-_e, _d, 0.0],
                             [-_e, -_d, 0.0],  # X grade

                             [0.0, _e, _d],  # Y
                             [0.0, _e, -_d],  # Y+ grade
                             [0.0, -_e, _d],  # Y- grade
                             [0.0, -_e, -_d],

                             [_d, 0, _e],
                             [-_d, 0, _e],  # Z+ grade
                             [_d, 0, -_e],  # Z- grade
                             [-_d, 0, -_e],])
        _f = self.base_sensor_edge / 2
        self.b_B = np.array([[_f, _d, 0.0],  # X
                             [_f, -_d, 0.0],
                             [-_f, _d, 0.0],
                             [-_f, -_d, 0.0],

                             [0.0, _f, _d],  # Y
                             [0.0, _f, -_d],
                             [0.0, -_f, _d],
                             [0.0, -_f, -_d],
                             [_d, 0, _f],
                             [-_d, 0, _f],
                             [_d, 0, -_f],
                             [-_d, 0, -_f],])
        # self.b_I = 1.0 * self.b_B
        # ''' Inertial vectors of connections of coils on sensor base'''
        # self.m_I = 1.0 * self.m_M
        # ''' Inertial vectors of connections of coils on seismic mass'''
        self.f_B = 0.0 * self.m_M
        ''' Vector f on base sensor coordinate system. This vector contains the length of the optical fiber at the current instant'''
        self.leg = ['xz', 'x-z', '-xz', '-x-z', 'yz', 'y-z',
                    '-yz', '-y-z', 'zy', 'z-y', '-zy', '-z-y']

    def update_f_vector(self,rm_B):
        ''' update deformation vector f '''
        for i in range(12):
            self.f_B[i,:] = rm_B + self.rot_M_B@self.m_M[i,:] - self.b_B[i,:]

    def update_states(self, rb_I, qb_I, rm_B, q_M_B):
        '''
        update_states Update the states of translation and rotation
        Args:
            rb_I: translation vector of body system with respect of inertial
            qb_I: quaternion of attitude of body sensor
            rm_B: translation vector of seismic mass with respect of base sensor system
            q_M_B: attitude quaternion of seismic mass w.r.t base sensor.
        '''

        self.q_M_B = q_M_B
        self.rot_M_B = fq.rotationMatrix(q_M_B)
        # The base of the sensor has coordinates in the inertial system
        self.bss.r = rb_I
        self.bss.updates_attitude(qb_I)
        _q = np.array([1., 0., 0., 0.])
        self.sms.r = self.bss.r+rm_B
        fq.MultQuat(p=qb_I, q=q_M_B,r=_q)
        self.sms.updates_attitude(_q)

        self.update_f_vector(rm_B)

    def dd_x_forced_body_state(self, d_x, u):
        '''
        dd_x_forced_body_state calc second order of model for numerical integration.
        Args:
            dd_x: second order give a first order states
            d_x: firts order [rb,rm,qb,qm,drb,drm,wb,wm]
        '''
        # d_x = np.arange(1, 21)
        dd_x = 1.0 * d_x
        d_rb = d_x[:3]
        '''Inertial velocity of body sensor'''
        d_rm = d_x[3:6]
        '''velocity of seismic mass wrt body sensor'''
        rb = d_x[6:9]
        '''inertial position of body sensor'''
        rm_B = d_x[9:12]
        '''position of seismic mass wrt body sensor'''
        qb = d_x[12:16]
        '''Atitude quaternion of body sensor'''
        qm_B = d_x[16:20]
        '''Atitude quaternion of seismic mass wrt body sensor'''
        wb = d_x[20:23]
        '''Angular velocity of body sensor'''
        wm = d_x[23:26]
        '''Angular velocity of seismic mass'''
        # calculate deformation vector
        # NOTE: can be optimized, but, for best understanding, we
        f_hat_Del_l = np.zeros((3, 12))
        sum_f_hat_Del_l = np.zeros(3)
        sum_f_hat_Del_l_dfdq_M = np.zeros(4)
        # sum_f_hat_Del_l_dfdq_B = np.zeros(4)
        for i in range(12):
            # compute f
            f_hat_Del_l[:, i] = rm_B + fq.rotationMatrix(
                qm_B) @ self.m_M[i, :]  - self.b_B[i, :]
            # calculate a norm of v to get deformation and versor of f
            f_norm = np.linalg.norm(f_hat_Del_l[:, i])
            # compute f_hat
            f_hat_Del_l[:, i] /= f_norm
            # compute f_hat_Del_l
            f_hat_Del_l[:, i] *= (f_norm - self.fiber_length)
            # sum for compute translational moviments
            sum_f_hat_Del_l += f_hat_Del_l[:, i]
            sum_f_hat_Del_l_dfdq_M += fq.calc_dfdq(qm_B,
                                               self.m_M[i, :]) @ f_hat_Del_l[:, i]
            # NOTE: important signal of dfdq
            # sum_f_hat_Del_l_dfdq_B = -calc_dfdq(qb, self.b_B[i, :]) @ f_hat_Del_l[:, i]
        # calculate dd_rm
        dd_x[3:6] = -self.k * sum_f_hat_Del_l / self.seismic_mass
        dd_x[5] += self.G
        # Artifictial damper
        dd_x[3:6] -= (self.damper_for_computation_simulations /
                      self.seismic_mass) * (d_rm-d_rb)

        # calculate d_rb
        dd_x[6:9] = d_rb
        # calculate d_rm
        dd_x[9:12] = d_rm
        # Left quaternion matrix
        Qb = fq.matrixQ(qb)
        Qm_B = fq.matrixQ(qm_B)
        # calculate attitude quaternion of body sensor
        dd_x[12:16] = 0.5 * Qb @ wb
        # calculate attitude quaternion of seismic mass wrt B
        dd_x[16:20] = 0.5 * Qm_B @ wm
        # calculate a angular acceleration of body sensor
        # NOTE: due to symmetry, the product fq.screwMatrix(wm) @ self.inertial_seismic_mass @ wm it is always zero!
        dd_x[23:26] = - self.inertial_seismic_mass_inv @  (
            0.5*self.k * Qm_B.T @ sum_f_hat_Del_l_dfdq_M)
        return dd_x
