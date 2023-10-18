#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Runge Kutta para sistemas na represetacao de espaco de estados
"""
import numpy as np


class RungeKutta(object):
    """Numerical integrator of Runge Kutta."""

    def __init__(self, num_states, function_of_integration):
        
        '''
        __init__ constructor of class

        Args:
            num_states: Number of states to be integrated.
            f: function of the problem.
        '''
        self.k = np.zeros(shape=(num_states, 4))
        self.function_of_integration = function_of_integration
        


    def integrates_states(self, q, u, h):
        
        """
        integrates_states
        Integration Equations Differentia with Runge Kutta Method of Order 4.
        Args:
        _extended_summary_
        Args:
            q: states
            u: forced vector
            h: step
        Returns:
            qprop: propagated states q_(k+1)
        """
        
        self.k[:, 0] = h * self.function_of_integration(q, u)
        self.k[:, 1] = h * self.function_of_integration(q + 0.5 * self.k[:, 0], u)
        self.k[:, 2] = h * self.function_of_integration(q + 0.5 * self.k[:, 1], u)
        self.k[:, 3] = h * self.function_of_integration(q + self.k[:, 2], u)
        return q + (self.k[:, 0] + 2.0 * (self.k[:, 1] + self.k[:, 2]) + self.k[:, 3]) / 6.0 # type: ignore
