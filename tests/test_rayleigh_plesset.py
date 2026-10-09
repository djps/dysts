"""
Test the Rayleigh-Plesset bubble models
> python -m unittest tests.test_rayleigh_plesset
"""
#!/usr/bin/env python
import os
import numpy as np
import unittest

WORKING_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import sys

sys.path.insert(1, os.path.join(WORKING_DIR, "dysts"))
from dysts.flows import RayleighPlesset, RayleighPlessetLog


def jac_central(model, x, h=1e-6):
    """Central-difference Jacobian with steps relative to each coordinate,
    since R is of order 1e-6 m"""
    jac = np.zeros((len(x), len(x)))
    for j in range(len(x)):
        step = h * abs(x[j]) if x[j] != 0 else h
        xp, xm = np.copy(x), np.copy(x)
        xp[j] += step
        xm[j] -= step
        jac[:, j] = (np.array(model.rhs(xp, 0.0)) - np.array(model.rhs(xm, 0.0))) / (2 * step)
    return jac


def sample_states(n=10, seed=0):
    """Random (R, U, Omega) states spanning compression and expansion"""
    rng = np.random.default_rng(seed)
    R = 2e-6 * np.exp(rng.uniform(-1.2, 1.0, n))
    U = rng.uniform(-200, 200, n)
    Omega = rng.uniform(0, 2 * np.pi, n)
    return zip(R, U, Omega)


class TestRayleighPlesset(unittest.TestCase):
    """
    Tests the bubble right hand side, Jacobians and log-radius variant
    """
    def test_rhs(self):
        """
        Compare the wall acceleration with the model written out in full
        """
        m = RayleighPlesset()
        R, U, Omega = 1.3e-6, 40.0, 1.1
        Pg = (m.Pa + 2 * m.sigma / m.R0) * ((m.R0**3 - m.a**3) / (R**3 - m.a**3))**m.gamma
        dPg = -3 * m.gamma * Pg * U * R**2 / (R**3 - m.a**3)
        p_inf = m.Pa - m.P * np.sin(Omega)
        dp_inf = -m.P * m.omega * np.cos(Omega)
        B = Pg - p_inf - 2 * m.sigma / R - 4 * m.mu * m.rho * U / R + (R / m.c) * (dPg - dp_inf)
        expected = (B / m.rho - 1.5 * U**2) / R
        out = m.rhs(np.array([R, U, Omega]), 0.0)
        assert np.isclose(out[1], expected, rtol=1e-12), "Wall acceleration does not match the model"
        assert out[0] == U and out[2] == m.omega

    def test_jacobian(self):
        """
        Compare the analytic Jacobian with finite differences
        """
        m = RayleighPlesset()
        for R, U, Omega in sample_states():
            x = np.array([R, U, Omega])
            assert np.allclose(m.jac(x, 0.0), jac_central(m, x), rtol=1e-6, atol=0), \
                "Analytic Jacobian does not match finite differences at {}".format(x)

    def test_jacobian_log(self):
        """
        Compare the log-radius analytic Jacobian with finite differences
        """
        m = RayleighPlessetLog()
        for R, U, Omega in sample_states():
            x = np.array([np.log(R), U, Omega])
            assert np.allclose(m.jac(x, 0.0), jac_central(m, x), rtol=1e-5, atol=0), \
                "Analytic Jacobian does not match finite differences at {}".format(x)

    def test_log_variant(self):
        """
        The log-radius variant should reproduce the base model's trajectory
        """
        sol = RayleighPlesset().make_trajectory(400, resample=False)
        sol_log = RayleighPlessetLog().make_trajectory(400, resample=False)
        assert sol.shape == sol_log.shape == (400, 3)
        assert np.allclose(sol_log[:, 0], sol[:, 0], rtol=0, atol=1e-6 * RayleighPlesset().R0), \
            "Log-radius variant does not reproduce the bubble radius"


if __name__ == "__main__":
    unittest.main()
