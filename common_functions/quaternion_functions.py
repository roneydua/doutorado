#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 09:36:05 2020
@author: roney
"""
import numpy as np
from numpy import cos, sin


def screw_matrix(v):
    """Computes the screw matrix from the v vector.
    Parameters
    ----------
    @param vector or quaternion.

    Returns
    -------
    q_x : 3x3 screw-Matrix
    """
    q_x = np.zeros((3, 3))
    q_x[0, 1] = -v[-1]
    q_x[1, 0] = v[-1]
    q_x[0, 2] = v[-2]
    q_x[2, 0] = -v[-2]
    q_x[1, 2] = -v[-3]
    q_x[2, 1] = v[-3]
    return q_x


def matrixQ(quat, right=False):
    """Compute Q matrix.

    Parameters
    ----------
    quat : quaternion orientation 4x1
    right : boolean. The default is False.

    Returns
    -------
    Q : Q matrix.
    """
    Q = np.zeros((4, 3))
    Q[0, 0] = -quat[1]
    Q[0, 1] = -quat[2]
    Q[0, 2] = -quat[3]
    if right:
        Q[1:, :] = quat[0] * np.eye(3) - screw_matrix(quat)
    else:
        Q[1:, :] = quat[0] * np.eye(3) + screw_matrix(quat)
    return Q


def matrixS(quat, right=False):
    """Compute S matrix.

    Parameters
    ----------
    quat : quaternion 4x1.
    right : boolean, optional. The default is False.

    Returns
    -------
    S : Matrix S.
    """
    S = np.zeros((4, 4))
    if right:
        S[:, 1:] = matrixQ(quat, right=True)
    else:
        S[:, 1:] = matrixQ(quat)

    S[:, 0] = quat
    return S


def mult_quat(p: np.ndarray, q: np.ndarray):
    """MultQuat Parameters.
    Args:
    -----
        q: quaternion of attitude
        p: quaternion of attitude
    Returns:
    -------
        r: result of multiplication. r = p x q
    """
    r = np.array([0.0, 0.0, 0.0, 0.0])
    r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3]
    r[1] = p[1] * q[0] + p[0] * q[1] - p[3] * q[2] + p[2] * q[3]
    r[2] = p[2] * q[0] + p[3] * q[1] + p[0] * q[2] - p[1] * q[3]
    r[3] = p[3] * q[0] - p[2] * q[1] + p[1] * q[2] + p[0] * q[3]
    return r


def conj(q):
    """Retorna o conjugado do quaternion q.

    Parameters
    ----------
    q : quaternion de atitude.

    Returns
    -------
    p : conjugado do quaternion de atitude q.
    """
    p = 1 * q
    p[:1] *= -1.0
    return p


def quat2Euler(q, deg=0):
    """Parameters.

    q : quaternion attitude.
    deg : TYPE, optional
        DESCRIPTION. 1 to return in degrees. The default is 0.

    Returns
    -------
    TYPE
        DESCRIPTION.
    e : 3x1 euler angles
    """
    e = np.array([0.0, 0.0, 0.0])
    e[2] = np.arctan2(q[0] * q[1] + q[2] * q[3], q[0] ** 2 + q[3] ** 2 - 0.5)
    e[1] = np.arcsin(2 * (q[0] * q[2] - q[1] * q[3]))
    e[0] = np.arctan2(q[0] * q[3] + q[1] * q[2], q[0] ** 2 + q[1] ** 2 - 0.5)
    if deg:
        return 180.0 * e / np.pi
    else:
        return e


def rotationMatrix(q):
    """Computa a matriz de rotacao a partir do quaternio de attitude.

    Parameters
    ----------
    q : quaternion de atitude.

    Returns
    -------
    Matrix de rotacao
    """
    return (
        np.eye(3) * (q[0] * q[0] - np.dot(q[1:], q[1:]))
        + 2.0 * q[0] * screw_matrix(q)
        + 2.0 * q[1:].reshape(3, 1) @ q[1:].reshape(1, 3)
    )


def eulerQuaternion(yaw, pitch, roll, deg=True):
    """Calcula o quaternio a partir dos Angulos de Euler.

    Parameters
    ----------
    yaw : angulo de guinad [rad]
    pitch : angulo de guinad [rad]
    roll : angulo de guinad [rad]
    deg : se 1 informa que os angulos estao em graus
    q : quaternio atitude
    """
    if deg:
        yaw *= np.pi / 180.0
        pitch *= np.pi / 180.0
        roll *= np.pi / 180.0
    cy = cos(yaw * 0.5)
    sy = sin(yaw * 0.5)
    cp = cos(pitch * 0.5)
    sp = sin(pitch * 0.5)
    cr = cos(roll * 0.5)
    sr = sin(roll * 0.5)
    q = np.zeros(4)
    q[0] = cr * cp * cy + sr * sp * sy
    q[1] = sr * cp * cy - cr * sp * sy
    q[2] = cr * sp * cy + sr * cp * sy
    q[3] = cr * cp * sy - sr * sp * cy
    return q


def expMap(w: np.ndarray, dt: float, tolerance=1e-15):
    # Compute norm of angular vector
    expMapArray = np.array([1.0, 0.0, 0.0, 0.0])
    w_norm = np.linalg.norm(w)
    if w_norm > tolerance:
        expMapArray[0] = cos(w_norm * 0.5 * dt)
        expMapArray[1:] = w / w_norm * sin(w_norm * 0.5 * dt)
    return expMapArray


def calc_dfdq(q, a):
    """
    calc_dfdq Calculate jacobian of rotated vector
    Args:
        q: quaternion rotation
        v: vector 3x1
    Returns:
        Jacobian of rotation vector w.r.t quaterion.
        NOTE: This function return transpose of Jacobian
    """
    # dfdq = np.zeros((4, 3))
    # dfdq[0, :] = q[0] * v.T - v.T @ screw_matrix(q[1:])
    # dfdq[1:, :] = q[1:].T @ v * np.eye(3)
    # dfdq[1:, :] += v.reshape((3, 1)) @ q[1:].reshape((1, 3))
    # dfdq[1:, :] -= q[1:].reshape((3, 1)) @ v.reshape((1, 3))
    # dfdq[1:, :] += q[0] * screw_matrix(v)
    # return 2.0 * dfdq
    dfdq = np.zeros((3, 4))
    dfdq[:, 0] = q[0] * a + screw_matrix(q[1:]) @ a
    dfdq[:, 1:] = q[1:] .T @ a * np.eye(3)
    dfdq[:, 1:] += q[1:].reshape((3, 1)) @ a.reshape((1, 3))
    dfdq[:, 1:] -= a.reshape((3, 1)) @ q[1:].reshape((1, 3))
    dfdq[:, 1:] -= q[0] * screw_matrix(a)
    return 2.0 * dfdq
