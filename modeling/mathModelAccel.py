import numpy as np
from numpy.linalg import inv
import funcoesQuaternion as fq


# Constates coordinates of connections



class states:
    # rotation matrix
    rot = np.eye(3)
    # psi, theta, phi
    euler = np.zeros(3)

    def __init__(self):

        self.q = np.array([1., 0., 0., 0.])
        """attitude quaternion"""
        self.w = np.array([0, 0, 0])
        """# angular body velcity"""
        self.r = np.array([0, 0, 0])
        """inertial position"""
        self.dr = np.array([0, 0, 0])
        """inertial velocity"""
        self.ddr = np.array([0, 0, 0])
        """inertial acceleration"""

    def updates_attitude(self, q):
        '''
        updates_attitude Update states of attitude
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
    """material density  (aluminiun)"""
    E = 70e9  # GPa
    '''Young's Module'''
    fiber_diameter = 125e-6
    '''fiber diameter'''
    sms = states()
    '''seismic mass states'''
    bss = states()
    '''base sensor states'''
    G = 0.0 #-9.89
    '''Gravity'''
    damper_for_computation_simulations = 0.0
    """Artificial damper for numérical simulations"""


    def __init__(self, seismic_edge=16.4e-3,
                 fiber_diameter=125e-6, fiber_length=3e-3, fibers_with_info=np.arange(1, 13), inverse_problem_full=True,damper_for_computation_simulations=0.0):
        self.damper_for_computation_simulations = damper_for_computation_simulations
        self.fibers_with_info = fibers_with_info
        # index of fiber with deformation informartion
        self.fiber_diameter = fiber_diameter
        # Fiber diameter
        self.fiber_length = fiber_length
        '''initial fiber length'''
        self.k = (self.E * 0.25 * np.pi * (self.fiber_diameter**2)) / self.fiber_length
        """stifness of optical fiber"""
        self.seismic_edge = seismic_edge
        '''Length of seismic edge'''
        self.seismic_mass = (self.seismic_edge ** 3.0) * self.density
        '''Sismic mass'''
        self.external_base_sensor_edge = self.seismic_edge + 2 * self.fiber_length + 4e-3
        """# approximation of the base of the sensor as a cube.The value 6e-3 refers double the length of the fibers that support the cube.4e-3 is twice the thickness."""
        self.base_sensor_edge = self.seismic_edge + 2 * self.fiber_length
        self.base_sensor_mass = (
            self.external_base_sensor_edge**3 - self.base_sensor_edge ** 3) * self.density
        self.inertial_seismic_mass = np.eye(3) * self.seismic_mass / 6.0 * self.seismic_edge**2
        self.inertial_base_sensor = np.eye(3) * self.base_sensor_mass / 6.0 * self.seismic_edge**2
        # NOTE: Only for computational economy
        self.inertial_base_sensor_inv = inv(self.inertial_base_sensor)
        self.inertial_seismic_mass_inv = inv(self.inertial_seismic_mass)
        # perpendicular distance of center
        _d = self.seismic_edge / 2 - 1e-3
        _e = self.seismic_edge / 2
        self.m_M = np.array([[_e, _d, 0.0],  # X+ GRADE
                             [_e, -_d, 0.0],
                             [-_e, _d, 0.0],
                             [-_e, -_d, 0.0],# X GRADE

                             [0.0, _e,_d], # Y 
                             [0.0, _e,-_d], # Y+ GRADE
                             [0.0, -_e,_d],# Y- GRADE
                             [0.0, -_e, -_d],

                             [ _d,0, _e],
                             [ -_d,0, _e], #Z+ grade
                             [ _d,0, -_e],#Z- grade
                             [ -_d,0, -_e],])
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
        self.b_I = 1.0 * self.b_B
        ''' Inertial vectors of connections of coils on sensor base'''
        self.m_I = 1.0 * self.m_M
        ''' Inertial vectors of connections of coils on seismic mass'''
        self.f = 0.0 * self.m_M
        ''' Vector f. This vector contains the length of the optical fiber at the current instant.'''
        self.leg = ['xz', 'x-z', '-xz', '-x-z', 'yz', 'y-z', '-yz', '-y-z', 'zy', 'z-y', '-zy', '-z-y']
        # legend of numerical point

    def update_inertial_coil_connections(self):
        '''
        update_inertial_coil_connections update the inertial coordinates of mass and body connections
        The deformation vector f is also updated.

        '''
        for i in range(12):
            # inertial vectors of connections of coils on sensor base
            self.m_I[i, :] = self.sms.r + self.sms.rot @ self.m_M[i, :]
            # inertial vectors of connections of coils on seismic mass
            self.b_I[i, :] = self.bss.r + self.bss.rot @ self.b_B[i, :]


    def update_f_vector(self):
        """ update deformation vector f """
        self.f = self.m_I - self.b_I


    def update_states(self, rb, qb, rm, qm):
        '''
        update_states Update the states of translation and rotation
        Args:
            rb: translation vector of body system with respect of inertial
            qb: quaternion of attitude of body sensor
            rm: translation vector of seismic mass with respect of inertial system
            qm: attitude quaternion of seismic mass.
        '''
        self.bss.r = rb
        self.bss.updates_attitude(qb)
        self.sms.r = rm
        self.sms.updates_attitude(qm)
        self.update_inertial_coil_connections()
        self.update_f_vector()

    # def dd_x(self, d_x, u):
    #     '''
    #     dd_x calc second order of model for numerical integration.
    #     Args:
    #         dd_x: second order give a first order states
    #         d_x: firts order [rb,rm,qb,qm,drb,drm,wb,wm]
    #     '''
    #     # d_x = np.arange(1, 21)
    #     dd_x = 0.0 * d_x

    #     def calc_dfdq(q, v):
    #         dfdq = np.zeros((4, 3))
    #         dfdq[0, :] = q[0] * v.T - v.T @ fq.screwMatrix(q[1:])
    #         dfdq[1:, :] = q[1:].T @ v * np.eye(3)
    #         dfdq[1:, :] += v.reshape((3, 1)) @ q[1:].reshape((1, 3))
    #         dfdq[1:, :] -= q[1:].reshape((3, 1)) @ v.reshape((1, 3))
    #         dfdq[1:, :] += q[0] * fq.screwMatrix(v)
    #         return 2.0 * dfdq
    #     # @param: inertial velocity of body sensor
    #     d_rb = d_x[:3]
    #     d_rm = d_x[3:6]
    #     rb = d_x[6:9]
    #     rm = d_x[9:12]
    #     qb = d_x[12:16]
    #     qm = d_x[16:20]
    #     wb = d_x[20:23]
    #     wm = d_x[23:26]
    #     # calculate deformation vector
    #     # NOTE: can be optimized, but, for best understanding, we
    #     f_hat_dell = np.zeros((3, 12))
    #     sum_f_hat_dell = np.zeros(3)
    #     sum_f_hat_dell_dfdq_M = np.zeros(4)
    #     sum_f_hat_dell_dfdq_B = np.zeros(4)
    #     for i in range(12):
    #         # compute f
    #         f_hat_dell[:, i] = rm + fq.rotationMatrix(qm) @ self.m_M[i, :] - rb - fq.rotationMatrix(qb) @ self.b_B[i, :]
    #         # calculate a norm of v to get deformation and versor of f
    #         f_norm = np.linalg.norm(f_hat_dell[:, i])
    #         # compute f_hat
    #         f_hat_dell[:, i] /= f_norm
    #         # compute f_hat_dell
    #         f_hat_dell[:, i] *= (f_norm - self.fiber_length)
    #         # sum for compute translational moviments
    #         sum_f_hat_dell += f_hat_dell[:, i]
    #         sum_f_hat_dell_dfdq_M += calc_dfdq(qm, self.m_M[i, :]) @ f_hat_dell[:, i]
    #         # NOTE: important signal of dfdq
    #         sum_f_hat_dell_dfdq_B = -calc_dfdq(qb, self.b_B[i, :]) @ f_hat_dell[:, i]
    #     # calculate dd_rb
    #     dd_x[:3] = (self.k * sum_f_hat_dell + u[:3]) / self.base_sensor_mass
    #     dd_x[2] += self.G
    #     # TEST: Artificial damper
    #     # dd_x[:3] += .2 * d_rb
    #     # HACK! Null acceleration to simulate accelerometer stopped!
    #     # dd_x[2] *=0
    #     # calculate dd_rm
    #     dd_x[3:6] = -self.k * sum_f_hat_dell / self.seismic_mass
    #     dd_x[5] += self.G
    #     # calculate d_rb
    #     dd_x[6:9] = d_rb
    #     # calculate d_rm
    #     dd_x[9:12] = d_rm
    #     # Forcing the sensor's body to stand still
    #     # dd_x[:3] *= 0
    #     # dd_x[6:9] *=0
    #     # Left quaternion matrix
    #     Qb = fq.matrixQ(qb)
    #     Qm = fq.matrixQ(qm)
    #     # calculate attitude quaternion of body sensor
    #     dd_x[12:16] = 0.5 * Qb @ wb
    #     # calculate attitude quaternion of seismic mass
    #     dd_x[16:20] = 0.5 * Qm @ wm
    #     # calculate a angular acceleration of body sensor
    #     dd_x[20:23] = -self.inertial_base_sensor_inv @ (fq.screwMatrix(wb) @
    #                                                     self.inertial_base_sensor @ wb + 0.5 * self.k * Qb.T @ sum_f_hat_dell_dfdq_B - u[3:])
    #     # calculate a angular acceleration of body sensor
    #     dd_x[23:26] = -self.inertial_seismic_mass_inv @ (fq.screwMatrix(wm) @
    #                                                     self.inertial_seismic_mass @ wm + 0.5 * self.k * Qm.T @ sum_f_hat_dell_dfdq_M)
    #     return dd_x

    def dd_x_forced_body_state(self, d_x, u):
        '''
        dd_x_forced_body_state calc second order of model for numerical integration.
        Args:
            dd_x: second order give a first order states
            d_x: firts order [rb,rm,qb,qm,drb,drm,wb,wm]
        '''
        # d_x = np.arange(1, 21)
        dd_x = 1.0 * d_x

        def calc_dfdq(q, v):
            dfdq = np.zeros((4, 3))
            dfdq[0, :] = q[0] * v.T - v.T @ fq.screwMatrix(q[1:])
            dfdq[1:, :] = q[1:].T @ v * np.eye(3)
            dfdq[1:, :] += v.reshape((3, 1)) @ q[1:].reshape((1, 3))
            dfdq[1:, :] -= q[1:].reshape((3, 1)) @ v.reshape((1, 3))
            dfdq[1:, :] += q[0] * fq.screwMatrix(v)
            return 2.0 * dfdq
        d_rb = d_x[:3]
        '''inertial velocity of body sensor'''
        d_rm = d_x[3:6]
        '''inertial velocity of seismic mass'''
        rb = d_x[6:9]
        '''inertial position of body sensor'''
        rm = d_x[9:12]
        '''inertial position of seismic mass'''
        qb = d_x[12:16]
        '''Atitude quaternion of body sensor'''
        qm = d_x[16:20]
        '''Atitude quaternion of seismic mass'''
        wb = d_x[20:23]
        '''Angular velocity of body sensor'''
        wm = d_x[23:26]
        '''Angular velocity of seismic mass'''
        # calculate deformation vector
        # NOTE: can be optimized, but, for best understanding, we
        f_hat_dell = np.zeros((3, 12))
        sum_f_hat_dell = np.zeros(3)
        sum_f_hat_dell_dfdq_M = np.zeros(4)
        # sum_f_hat_dell_dfdq_B = np.zeros(4)
        for i in range(12):
            # compute f
            f_hat_dell[:, i] = rm + fq.rotationMatrix(qm) @ self.m_M[i, :] - rb - fq.rotationMatrix(qb) @ self.b_B[i, :]
            # calculate a norm of v to get deformation and versor of f
            f_norm = np.linalg.norm(f_hat_dell[:, i])
            # compute f_hat
            f_hat_dell[:, i] /= f_norm
            # compute f_hat_dell
            f_hat_dell[:, i] *= (f_norm - self.fiber_length)
            # sum for compute translational moviments
            sum_f_hat_dell += f_hat_dell[:, i]
            sum_f_hat_dell_dfdq_M += calc_dfdq(qm, self.m_M[i, :]) @ f_hat_dell[:, i]
            # NOTE: important signal of dfdq
            # sum_f_hat_dell_dfdq_B = -calc_dfdq(qb, self.b_B[i, :]) @ f_hat_dell[:, i]
        # calculate dd_rm
        dd_x[3:6] = -self.k * sum_f_hat_dell / self.seismic_mass
        dd_x[5] += self.G
        ## Artifictial damper
        dd_x[3:6] -= (self.damper_for_computation_simulations /
                      self.seismic_mass) * (d_rm-d_rb)

        # calculate d_rb
        dd_x[6:9] = d_rb
        # calculate d_rm
        dd_x[9:12] = d_rm
        # Left quaternion matrix
        Qb = fq.matrixQ(qb)
        Qm = fq.matrixQ(qm)
        # calculate attitude quaternion of body sensor
        dd_x[12:16] = 0.5 * Qb @ wb
        # calculate attitude quaternion of seismic mass
        dd_x[16:20] = 0.5 * Qm @ wm
        # calculate a angular acceleration of body sensor
        # dd_x[20:23] = -self.inertial_base_sensor_inv @ (fq.screwMatrix(wb) @
        # self.inertial_base_sensor @ wb + 0.5 * self.k * Qb.T @ sum_f_hat_dell_dfdq_B - u[3:])
        # calculate a angular acceleration of body sensor
        ## ATTENTION! due to symmetry, the product fq.screwMatrix(wm) @ self.inertial_seismic_mass @ wm it is always zero!
        dd_x[23:26] = - self.inertial_seismic_mass_inv @  (0.5*self.k  *Qm.T@ sum_f_hat_dell_dfdq_M)
        return dd_x
