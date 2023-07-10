#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 15:03:46 2020

@author: roney
"""


import sympy as sp


class QuaternionSymbolic(object):
    """Documentation for quaternion
    """

    def __init__(self, diff=False, sub="", coord=[]):
        if diff:
            q0 = sp.Symbol(r'\dot{q}_{0}'+sub)
            q1 = sp.Symbol(r'\dot{q}_{1}'+sub)
            q2 = sp.Symbol(r'\dot{q}_{2}'+sub)
            q3 = sp.Symbol(r'\dot{q}_{3}'+sub)
            self.quat = sp.Matrix([q0, q1, q2, q3])
        elif len(coord) == 0:
            q0 = sp.Symbol(r'q_'+'0_'+sub)
            q1 = sp.Symbol(r'q_'+'1_'+sub)
            q2 = sp.Symbol(r'q_'+'2_'+sub)
            q3 = sp.Symbol(r'q_'+'3_'+sub)
            self.quat = sp.Matrix([q0, q1, q2, q3])
        else:

            self.quat = sp.Matrix([coord[0], coord[1], coord[2], coord[3]])

    def conj(self):
        # temp1 =[1:, 0] *= -1
        return QuaternionSymbolic(coord=[self.quat[0], -self.quat[1], -self.quat[2], -self.quat[3]])

    def norm(self, sqrt=False):
        if sqrt:
            return sp.sqrt((self.quat.T * self.quat)[0, 0])
        else:
            return (self.quat.T * self.quat)[0, 0]

    def screw(self):
        q_x = sp.Matrix(sp.zeros(3, 3))
        q_x[0, 1] = -self.quat[3]
        q_x[1, 0] = self.quat[3]
        q_x[0, 2] = self.quat[2]
        q_x[2, 0] = -self.quat[2]
        q_x[1, 2] = -self.quat[1]
        q_x[2, 1] = self.quat[1]
        return q_x
    

    def quatRot(self, normalized=0):
        """
        @brief      Function that covert a quaternion on matrix rotation

        @details    detailed description

        @param      normalized = 1 return a matrix with normalized quaternion

        @return     sympy matrix 3 x 3
        """
        if normalized < 1:
            iden = sp.Matrix(sp.Identity(3))
            iden *= self.quat[0, 0]**2 - self.quat[1:, 0].dot(self.quat[1:, 0])
            iden += 2*self.quat[0, 0] * self.screw()
            iden += 2 * self.quat[1:, 0]*self.quat[1:, 0].T
        else:
            iden = sp.Matrix(sp.Identity(3))
            iden *= 2*self.quat[0, 0]**2 - 1
            iden += 2*self.quat[0, 0] * self.screw()
            iden += 2 * self.quat[1:, 0]*self.quat[1:, 0].T
        return iden

    def multQuat(self, other):
        a = QuaternionSymbolic()
        aQuatTemp = sp.Matrix(sp.zeros(4, 1))
        aQuatTemp[0, 0] = self.quat[0, 0] * other.quat[0, 0]-self.quat[1, 0] * other.quat[1, 0] - \
            self.quat[2, 0] * other.quat[2, 0] - \
            self.quat[3, 0] * other.quat[3, 0]
        aQuatTemp[1, 0] = self.quat[0, 0] * other.quat[1, 0]+self.quat[1, 0] * other.quat[0, 0] + \
            self.quat[2, 0] * other.quat[3, 0] - \
            self.quat[3, 0] * other.quat[2, 0]
        aQuatTemp[2, 0] = self.quat[0, 0] * other.quat[2, 0]-self.quat[1, 0] * other.quat[3, 0] + \
            self.quat[2, 0] * other.quat[0, 0] + \
            self.quat[3, 0] * other.quat[1, 0]
        aQuatTemp[3, 0] = self.quat[0, 0] * other.quat[3, 0]+self.quat[1, 0] * other.quat[2, 0] - \
            self.quat[2, 0] * other.quat[1, 0] + \
            self.quat[3, 0] * other.quat[0, 0]
        a.quat = aQuatTemp
        return a

    # def RotVetQuat(p_A_B, a_B):

    #     # =============================================================================
    #     #     para um quaternion A->B  tal que leve transorme uma coordenada do referencial A
    #     #     para o referencial B temos
    #     #     ap : queternion vetor 4x1 de A para B
    #     #     av : vetor a ser rotacionado no sistema A
    #     # =============================================================================
    #     p_A_Bconju = 1 * p_A_B
    #     p_A_Bconju[1:, 0] *= -1
    #     if len(a_B) < 4:
    #         vquat = 0 * p_A_B
    #         vquat[1:, 0] = a_B
    #         return multQuat(p_A_B, multQuat(vquat, p_A_Bconju))
    #     else:
    #         return multQuat(p_A_B, multQuat(a_B, p_A_Bconju))


    def Q(self, right=False):
        if right:
            S = sp.Matrix(sp.zeros(4, 3))
            S[0, 0] = -self.quat[1]
            S[0, 1] = -self.quat[2]
            S[0, 2] = -self.quat[3]
            S[1:, :] = self.quat[0]*sp.Matrix(sp.Identity(3)) - self.screw()
            return S
        else:
            S = sp.Matrix(sp.zeros(4, 3))
            S[0, 0] = -self.quat[1]
            S[0, 1] = -self.quat[2]
            S[0, 2] = -self.quat[3]
            S[1:, :] = self.quat[0]*sp.Matrix(sp.Identity(3)) + self.screw()
            return S

    def S(self, right=False):
        if right:
            S = sp.Matrix(sp.zeros(4, 4))
            S[:, 1:] = self.Q(right=True)
            S[:, 0] = self.quat
            return S
        else:
            S = sp.Matrix(sp.zeros(4, 4))
            S[:, 1:] = self.Q()
            S[:, 0] = self.quat
            return S

    def vec(self):
        return self.quat[1:, 0]
