"""The bridge between pynamix.measure and pynamix.spectral conventions.

These tests exist because the two halves of the package disagree in ways that
produce plausible numbers rather than errors.  Each one pins down a
disagreement that has bitten us.
"""

import unittest

import numpy as np

from pynamix.exposure import mean_std
from pynamix.spectral import compat
from pynamix.spectral.psd import radial_average, welch_psd2d


class TestWidthsAndAxes(unittest.TestCase):
    def test_half_width_round_trip(self):
        self.assertEqual(compat.roi_from_patchw(32), 64)
        self.assertEqual(compat.patchw_from_roi(64), 32)
        with self.assertRaises(ValueError):
            compat.patchw_from_roi(65)

    def test_wavelength_and_wavevector_carry_the_two_pi(self):
        # a wavelength axis is cycles per unit length; q is radians per unit length
        q = compat.q_from_wavelength(2.0)
        self.assertAlmostEqual(float(q), np.pi, places=12)
        self.assertAlmostEqual(float(compat.wavelength_from_q(q)), 2.0, places=12)
        lam = np.array([np.inf, 4.0, 2.0, 1.0])
        self.assertEqual(compat.q_from_wavelength(lam)[0], 0.0)


class TestAnnularSum(unittest.TestCase):
    def test_sum_to_mean_divides_by_the_mode_count(self):
        # three bins holding 4, 9 and 16 modes, each of true density 2.0
        n_j = np.array([4.0, 9.0, 16.0])
        nr_pxr = np.zeros((8, 8, 3))
        for k in range(3):
            nr_pxr[:, :, k] = n_j[k] / 64
        radialspec = 2.0 * n_j  # what an annular *sum* reports
        mean, got_n = compat.annular_sum_to_mean(radialspec, nr_pxr)
        self.assertTrue(np.allclose(got_n, n_j))
        self.assertTrue(np.allclose(mean, np.full(3, 2.0)))

    def test_rejects_a_shape_mismatch(self):
        with self.assertRaises(ValueError):
            compat.annular_sum_to_mean(np.zeros(5), np.zeros((4, 4, 7)))

    def test_the_jacobian_moves_the_peak(self):
        """An annular sum peaks at a shorter wavelength than the density it sums.

        This is the whole reason the two halves of the package report different
        'sizes' for the same image, so it is worth having it pinned down.
        """
        q = np.linspace(0.5, 20.0, 400)
        P = np.exp(-0.5 * ((q - 7.0) / 1.5) ** 2)  # density peaks at q = 7
        q_sum = q[np.argmax(q * P)]  # annular sum carries a factor q
        self.assertGreater(q_sum, 7.0)  # shifted to higher q, shorter lambda
        self.assertAlmostEqual(q[np.argmax(P)], 7.0, delta=0.05)


class TestNormalisationGuard(unittest.TestCase):
    def test_standardised_patches_are_refused(self):
        rng = np.random.default_rng(0)
        raw = 3.0 + 0.4 * rng.normal(size=(64, 64))
        compat.assert_unnormalised(raw)  # physical units: fine
        with self.assertRaisesRegex(ValueError, "standardised"):
            compat.assert_unnormalised(mean_std(raw))

    def test_standardisation_really_does_destroy_the_amplitude(self):
        """The reason assert_unnormalised exists: nothing downstream can catch it."""
        rng = np.random.default_rng(1)
        field = 3.0 + 0.4 * rng.normal(size=(128, 128))
        p = 0.05
        raw = radial_average(welch_psd2d(field, p, roi=128, window="boxcar"))
        std = radial_average(welch_psd2d(mean_std(field), p, roi=128, window="boxcar"))
        ratio = np.nanmedian(std.P / raw.P)
        # the shapes agree, so nothing downstream looks wrong ...
        self.assertLess(np.nanstd(std.P / raw.P) / ratio, 1e-6)
        # ... but the absolute amplitude, which sets phi, is off by 1/var
        self.assertAlmostEqual(ratio / (1.0 / np.var(field)), 1.0, places=6)
        self.assertFalse(np.isclose(ratio, 1.0, rtol=0.5))


if __name__ == "__main__":
    unittest.main()
