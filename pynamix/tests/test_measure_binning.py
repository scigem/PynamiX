"""Regression tests for the window and the radial/angular binning.

Each of these pins down something that was wrong before v0.36 and would not
have announced itself: an off-centre window, Monte Carlo coefficients that
differed between users, and radial bins labelled by their upper edge.
"""

import unittest
import warnings

import numpy as np

from pynamix.measure import (
    _check_wavelength_window,
    _parabolic_peak,
    angular_binning,
    hanning_window,
    radial_grid,
)


class TestHanningWindowCentring(unittest.TestCase):
    """The window must be centred on the patch, not a pixel off it."""

    def test_window_is_centred_and_symmetric(self):
        w = hanning_window(32)
        self.assertEqual(w.shape, (64, 64))
        # a window centred on the patch is invariant under a 180 degree rotation
        self.assertLess(np.abs(w - w[::-1, ::-1]).max(), 1e-12)
        # peak at the centre, zero outside the inscribed circle
        self.assertAlmostEqual(w[31, 31], w.max(), places=12)
        self.assertEqual(w[0, 0], 0.0)
        self.assertEqual(w[-1, -1], 0.0)

    def test_window_is_a_raised_cosine_in_radius(self):
        patchw = 16
        w = hanning_window(patchw)
        i, j = np.indices(w.shape)
        d = np.sqrt((i - (patchw - 0.5)) ** 2 + (j - (patchw - 0.5)) ** 2)
        inside = d <= patchw
        expect = 0.5 * np.cos(np.pi * d[inside] / patchw) + 0.5
        self.assertTrue(np.allclose(w[inside], expect))


class TestAngularBinning(unittest.TestCase):
    """Deterministic quadrature, not unseeded Monte Carlo."""

    def test_is_deterministic_and_exact(self):
        a, b = angular_binning(16), angular_binning(16)
        self.assertTrue(np.array_equal(a, b))
        # k-hat is a unit vector, so the tensor has unit trace at every pixel
        self.assertLess(np.abs(a[:, :, 0, 0] + a[:, :, 1, 1] - 1.0).max(), 1e-12)
        # and it is symmetric, exactly
        self.assertTrue(np.array_equal(a[:, :, 0, 1], a[:, :, 1, 0]))

    def test_matches_the_analytic_limit_far_from_the_origin(self):
        # a pixel many pixels from the origin subtends a tiny angle, so the
        # average of k-hat k-hat over it is just k-hat k-hat at its centre
        patchw = 16
        Q = angular_binning(patchw)
        row, col = 2, 3
        y, x = row - patchw + 0.5, col - patchw + 0.5
        n = np.hypot(x, y)
        self.assertAlmostEqual(Q[row, col, 0, 0], (x / n) ** 2, delta=2e-3)
        self.assertAlmostEqual(Q[row, col, 1, 1], (y / n) ** 2, delta=2e-3)

    def test_the_retired_monte_carlo_argument_warns(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            angular_binning(8, N=1000)
        self.assertTrue(any(issubclass(c.category, DeprecationWarning) for c in caught))


class TestRadialGrid(unittest.TestCase):
    def test_partitions_every_pixel_exactly(self):
        r_grid, nr_pxr = radial_grid(rnb=64, patchw=16)
        # every pixel's area is fully accounted for across the bins
        self.assertLess(np.abs(nr_pxr.sum(axis=2) - 1.0).max(), 1e-12)
        self.assertTrue(np.all(np.diff(r_grid) > 0))
        self.assertGreaterEqual(r_grid[0], 0.0)

    def test_labels_bins_by_their_mean_radius_not_their_edge(self):
        """The old code shifted r_grid up by half a bin after binning."""
        rnb, patchw = 200, 32
        r_grid, nr_pxr = radial_grid(rnb=rnb, patchw=patchw)
        nodes = np.linspace(0, patchw * 1.5, rnb)
        delta = nodes[1] - nodes[0]
        occupied = nr_pxr.sum(axis=(0, 1)) > 0
        offset = (r_grid - nodes)[occupied]
        # each label sits within half a bin of its node, not systematically above it
        self.assertLess(np.abs(offset).max(), 0.75 * delta)
        self.assertLess(abs(np.median(offset)), 0.25 * delta)

    def test_is_deterministic(self):
        a, _ = radial_grid(rnb=32, patchw=8)
        b, _ = radial_grid(rnb=32, patchw=8)
        self.assertTrue(np.array_equal(a, b))


class TestParabolicPeak(unittest.TestCase):
    def test_finds_a_known_vertex(self):
        x = np.linspace(0.0, 4.0, 5)
        y = (-((x - 2.3) ** 2) + 5.0)[None, :]  # vertex at 2.3
        i = np.array([np.argmax(y[0])])
        self.assertAlmostEqual(_parabolic_peak(x, y, i)[0], 2.3, places=9)

    def test_falls_back_when_there_is_no_maximum(self):
        x = np.linspace(0.0, 4.0, 5)
        y = np.array([[0.0, 1.0, 2.0, 3.0, 4.0]])  # monotonic: no interior vertex
        self.assertAlmostEqual(_parabolic_peak(x, y, np.array([2]))[0], x[2], places=12)


class TestWavelengthWindowGuard(unittest.TestCase):
    def test_fires_near_dc(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _check_wavelength_window(20.0, 32, 1.0, "test")  # 3.2 bins from DC
        self.assertTrue(any("bins from DC" in str(c.message) for c in caught))

    def test_stays_quiet_when_there_is_room(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _check_wavelength_window(6.0, 32, 1.0, "test")  # 10.7 bins from DC
        self.assertEqual(len(caught), 0)


if __name__ == "__main__":
    unittest.main()
