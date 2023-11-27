#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Runge Kutta for state space representation
"""
import numpy as np


class RungeKutta(object):
    """Numerical integrator of Runge Kutta."""

    def __init__(self, num_states, function_of_integration):
        """
        __init__ constructor of class

        Args:
            num_states: Number of states to be integrated.
            f: function of the problem.
        """
        self.k = np.zeros(shape=(num_states, 4))
        self.function_of_integration = function_of_integration

    def integrates_states(self, q: np.ndarray, h):
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
            propagated states q_(k+1)
        """

        self.k[:, 0] = h * self.function_of_integration(0, q)
        self.k[:, 1] = h * self.function_of_integration(0, q + 0.5 * self.k[:, 0])
        self.k[:, 2] = h * self.function_of_integration(0, q + 0.5 * self.k[:, 1])
        self.k[:, 3] = h * self.function_of_integration(0, q + self.k[:, 2])
        return (
            q
            + (self.k[:, 0] + 2.0 * (self.k[:, 1] + self.k[:, 2]) + self.k[:, 3]) / 6.0
        )


class RungeKuttaFehlberg45(object):
    """docstring for RungeKuttaGJ5."""

    def __init__(self, num_states, function_of_integration):
        """
        __init__ constructor of class

        Args:
            num_states: Number of states to be integrated.
            f: function of the problem.
        """
        qn = np.zeros(shape=(num_states))
        self.k = np.zeros(shape=(num_states, 6))
        self.function_of_integration = function_of_integration
        self.b21 = 1.0 / 4.0

        self.b31 = 3.0 / 32.0
        self.b32 = 9.0 / 32.0

        self.b41 = 1932.0 / 2197.0
        self.b42 = -7200.0 / 2197.0
        self.b43 = 7296.0 / 2197.0

        self.b51 = 439.0 / 216.0
        self.b52 = -8.0
        self.b53 = 3680.0 / 513.0
        self.b54 = -845.0 / 4104.0

        self.b61 = -8.0 / 27.0
        self.b62 = 2.0
        self.b63 = -3544.0 / 2565.0
        self.b64 = 1859.0 / 4104.0
        self.b65 = -11.0 / 40.0

        self.c1 = 16.0 / 135.0
        self.c2 = 0.0
        self.c3 = 6656.0 / 12825.0
        self.c4 = 18561.0 / 56430.0
        self.c5 = -9.0 / 50.0
        self.c6 = 2.0 / 55.0

    def integrates_states(self, q: np.ndarray, h):
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
            propagated states q_(k+1)
        """

        self.k[:, 0] = h * self.function_of_integration(q)
        self.k[:, 1] = h * self.function_of_integration(q + self.b21 * self.k[:, 0])
        self.k[:, 2] = h * self.function_of_integration(
            q + self.b31 * self.k[:, 0] + self.b32 * self.k[:, 1]
        )
        self.k[:, 3] = h * self.function_of_integration(
            q
            + self.b41 * self.k[:, 0]
            + self.b42 * self.k[:, 1]
            + self.b43 * self.k[:, 2]
        )
        self.k[:, 4] = h * self.function_of_integration(
            q
            + self.b51 * self.k[:, 0]
            + self.b52 * self.k[:, 1]
            + self.b53 * self.k[:, 2]
            + self.b54 * self.k[:, 3]
        )
        self.k[:, 5] = h * self.function_of_integration(
            q
            + self.b61 * self.k[:, 0]
            + self.b62 * self.k[:, 1]
            + self.b63 * self.k[:, 2]
            + self.b64 * self.k[:, 3]
            + self.b65 * self.k[:, 4]
        )
        return (
            q
            + self.c1 * self.k[:, 0]
            + self.c2 * self.k[:, 1]
            + self.c3 * self.k[:, 2]
            + self.c4 * self.k[:, 3]
            + self.c5 * self.k[:, 4]
            + self.c6 * self.k[:, 5]
        )
