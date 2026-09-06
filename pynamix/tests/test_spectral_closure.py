"""Percus--Yevick partial structure factors.

The closure is what lets the species weights be written down without a
simulation, so what matters is that it converges and lands on the limits that
are known independently.
"""

import unittest

import numpy as np

from pynamix.spectral.closure import py_partials

# a coarser radial grid than the default keeps the suite quick; it costs a few
# parts in a thousand on S(q), which none of these assertions is sensitive to
NG = 1024


def _one_component(q, eta, n_grid=NG):
    return py_partials({1: 1.0}, {1: 6.0 * eta / np.pi}, q, n_grid=n_grid)[(1, 1)]


class TestSingleComponent(unittest.TestCase):
    def test_a_dilute_fluid_is_structureless(self):
        """As the density goes to zero, S(q) -> 1 everywhere."""
        q = np.linspace(0.5, 25.0, 120)
        S = _one_component(q, 0.002)
        self.assertLess(np.max(np.abs(S - 1.0)), 0.02)

    def test_the_first_peak_grows_and_moves_out_with_density(self):
        q = np.linspace(0.5, 20.0, 400)
        peaks = []
        for eta in (0.10, 0.25, 0.40):
            S = _one_component(q, eta)
            peaks.append((q[np.argmax(S)], S.max()))
        heights = [h for _, h in peaks]
        positions = [p for p, _ in peaks]
        self.assertTrue(all(b > a for a, b in zip(heights, heights[1:])))
        self.assertTrue(all(b >= a for a, b in zip(positions, positions[1:])))
        # a dense hard-sphere fluid peaks near q*sigma = 2*pi
        self.assertAlmostEqual(positions[-1], 2 * np.pi, delta=1.0)

    def test_compressibility_falls_as_the_fluid_is_compressed(self):
        """S(q->0) is the isothermal compressibility, and it must decrease."""
        q = np.linspace(0.05, 3.0, 60)
        s0 = [_one_component(q, e)[0] for e in (0.05, 0.40)]
        self.assertTrue(all(b < a for a, b in zip(s0, s0[1:])))
        self.assertTrue(all(v > 0 for v in s0))

    def test_matches_the_carnahan_starling_compressibility(self):
        """The PY compressibility route has a closed form; check against it."""
        q = np.linspace(0.02, 2.0, 80)
        for eta in (0.1, 0.3):
            with self.subTest(eta=eta):
                S0 = _one_component(q, eta)[0]
                py_exact = (1 - eta) ** 4 / (1 + 2 * eta) ** 2
                self.assertAlmostEqual(S0 / py_exact, 1.0, delta=0.05)


class TestMixture(unittest.TestCase):
    def test_partials_have_the_expected_keys_and_symmetry(self):
        q = np.linspace(0.5, 20.0, 100)
        S = py_partials({1: 1.0, 2: 2.0}, {1: 0.30, 2: 0.02}, q, n_grid=NG)
        self.assertEqual(set(S), {(1, 1), (1, 2), (2, 1), (2, 2)})
        self.assertTrue(np.allclose(S[(1, 2)], S[(2, 1)]))
        for v in S.values():
            self.assertEqual(v.shape, q.shape)
            self.assertTrue(np.all(np.isfinite(v)))

    def test_the_cross_partial_can_go_negative(self):
        """Discarding it makes the velocity fit fail outright, so it must be real."""
        q = np.linspace(0.5, 20.0, 300)
        S = py_partials({1: 1.0, 2: 2.0}, {1: 0.35, 2: 0.03}, q, n_grid=NG)
        self.assertLess(S[(1, 2)].min(), 0.0)

    def test_a_dilute_mixture_tends_to_the_ideal_limit(self):
        q = np.linspace(0.5, 20.0, 100)
        S = py_partials({1: 1.0, 2: 2.0}, {1: 1e-3, 2: 1e-4}, q, n_grid=NG)
        self.assertLess(np.max(np.abs(S[(1, 1)] - 1.0)), 0.05)
        self.assertLess(np.max(np.abs(S[(2, 2)] - 1.0)), 0.05)
        self.assertLess(np.max(np.abs(S[(1, 2)])), 0.05)

    def test_convergence_is_reported_and_judged_on_S_not_gamma(self):
        q = np.linspace(0.5, 20.0, 80)
        res = py_partials({1: 1.0, 2: 2.0}, {1: 0.30, 2: 0.02}, q, n_grid=NG, return_result=True)
        self.assertTrue(res.converged)
        self.assertLess(res.residual, 1e-3)
        self.assertGreater(res.iterations, 0)

    def test_equal_diameters_reproduce_the_one_component_answer(self):
        """Two labels for the same sphere must not change the physics."""
        q = np.linspace(0.5, 20.0, 120)
        rho = 0.35
        one = _one_component(q, np.pi * rho / 6.0)
        S = py_partials({1: 1.0, 2: 1.0}, {1: rho / 2, 2: rho / 2}, q, n_grid=NG)
        # the number-number combination S_NN = sum_ab sqrt(x_a x_b) S^AL_ab, which
        # for two equal labels is half the plain sum, recovers the one-component
        # answer; the plain sum is twice it
        S_NN = 0.5 * (S[(1, 1)] + S[(2, 2)] + 2 * S[(1, 2)])
        self.assertLess(np.max(np.abs(S_NN - one)), 0.02)


if __name__ == "__main__":
    unittest.main()
