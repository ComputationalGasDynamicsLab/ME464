"""Numerical regression checks; run: python test_heat1d.py."""
import unittest
import numpy as np
from heat1d import make_stepper, amplification, stability_limits


class ThetaMethodTests(unittest.TestCase):
    def test_named_endpoints_and_boundary_forcing(self):
        for bc in ('dd', 'dn', 'nd', 'nn'):
            for theta, name in ((0, 'ftcs'), (0.5, 'cn'), (1, 'btcs')):
                step = make_stepper(10, 0.1, 0.4, bc, 0.3, -0.2, theta)
                a = np.cos(np.arange(11)) + 0.2 * np.arange(11)
                b = a.copy()
                for _ in range(3):
                    step(a, 'theta')
                    step(b, name)
                np.testing.assert_allclose(a, b, atol=2e-12)

    def test_grid_eigenmodes_and_stability_thresholds(self):
        N = 8
        for bc in ('dd', 'dn', 'nd', 'nn'):
            phase = np.pi if bc == 'nn' else np.pi * (N - (1 if bc == 'dd' else 0.5)) / N
            mode = (np.sin if bc[0] == 'd' else np.cos)(np.arange(N + 1) * phase)
            if bc[0] == 'd':
                mode[0] = 0
            if bc[1] == 'd':
                mode[-1] = 0
            for theta in (0, 0.25, 0.49, 0.5, 0.75, 1):
                limit, _ = stability_limits(theta, N, bc)
                rs = (limit * 0.99, limit, limit * 1.01) if np.isfinite(limit) else (0.4, 5, 50)
                for r in rs:
                    T = mode.copy()
                    make_stepper(N, 1 / N, r, bc, 0, 0, theta)(T, 'theta')
                    g = amplification(theta, r, np.sin(phase / 2)**2)
                    np.testing.assert_allclose(T, mode * g, atol=2e-12)
                    if np.isfinite(limit):
                        self.assertEqual(abs(g) > 1 + 1e-12, r > limit)

    def test_neumann_heat_balance(self):
        N, dx, r = 20, 0.05, 0.4
        for theta in (0, 0.25, 0.5, 0.75, 1):
            for g0, g1 in ((0, 0), (0.3, 0.8)):
                T = np.exp(-80 * (np.linspace(0, 1, N + 1) - 0.35)**2)
                weights = np.ones(N + 1) * dx
                weights[[0, -1]] *= 0.5
                old = T @ weights
                step = make_stepper(N, dx, r, 'nn', g0, g1, theta)
                for _ in range(20):
                    step(T, 'theta')
                self.assertAlmostEqual(T @ weights, old + 20 * r * dx**2 * (g1 - g0), places=12)

    def test_temporal_order_against_semidiscrete_exact(self):
        # Fixed spatial grid isolates time error, avoiding spatial-error masking.
        N, final = 20, 10
        mode = np.sin(np.arange(N + 1) * np.pi / N)
        mode[[0, -1]] = 0
        decay = np.exp(-4 * 0.01 * final * N**2 * np.sin(np.pi / (2 * N))**2)
        for theta, expected_ratio in ((0.25, 2), (0.5, 4), (0.75, 2), (1, 2)):
            errors = []
            for steps in (80, 160, 320):
                r = 0.01 * final * N**2 / steps
                step = make_stepper(N, 1 / N, r, 'dd', 0, 0, theta)
                T = mode.copy()
                for _ in range(steps):
                    step(T, 'theta')
                errors.append(np.max(np.abs(T - mode * decay)))
            for coarse, fine in zip(errors, errors[1:]):
                self.assertAlmostEqual(coarse / fine, expected_ratio, delta=0.04)

    def test_coarse_dirichlet_exception(self):
        limit, _ = stability_limits(0, 4, 'dd')
        self.assertAlmostEqual(limit, 0.585786437626905)
        self.assertLess(0.52, limit)


if __name__ == '__main__':
    unittest.main()
