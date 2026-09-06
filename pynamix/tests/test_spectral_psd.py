"""Spectral estimation: normalisation, windows, and the one case with no free parameters."""

import unittest

import numpy as np

from pynamix.spectral import psd, units
from pynamix.tests.spectral_fixtures import poisson_spheres, project


class TestNormalisation(unittest.TestCase):
    """The absolute amplitude is the whole point, so pin it down."""

    def test_parseval_holds_for_every_window(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(64, 64))
        for kind in ("boxcar", "hann", "tukey"):
            with self.subTest(window=kind):
                self.assertLess(abs(psd.parseval_residual(x, 0.1, window=kind)), 1e-10)

    def test_white_noise_has_a_flat_spectrum_of_the_right_height(self):
        """For white noise of variance s^2 on pitch p, P = s^2 p^2."""
        rng = np.random.default_rng(1)
        sigma, p = 0.7, 0.05
        x = rng.normal(0.0, sigma, size=(512, 512))
        r = psd.radial_average(psd.welch_psd2d(x, p, roi=512, window="boxcar"))
        band = r.q > 0.2 * r.q.max()
        self.assertAlmostEqual(np.median(r.P[band]) / (sigma**2 * p**2), 1.0, delta=0.02)

    def test_psd_has_units_of_length_to_the_fourth(self):
        """Halving the pitch at fixed field must not change the spectral density."""
        rng = np.random.default_rng(2)
        x = rng.normal(size=(256, 256))
        a = psd.radial_average(psd.welch_psd2d(x, 0.02, roi=256, window="boxcar"))
        b = psd.radial_average(psd.welch_psd2d(x, 0.04, roi=256, window="boxcar"))
        # same field, pitch doubled: q halves and P scales as p^2
        self.assertAlmostEqual(np.median(b.P) / np.median(a.P), 4.0, delta=0.02)
        self.assertAlmostEqual(b.q.max() / a.q.max(), 0.5, delta=1e-6)


class TestWindowsAndBins(unittest.TestCase):
    def test_taper_shapes_and_limits(self):
        self.assertTrue(np.allclose(psd.taper(32, "boxcar"), 1.0))
        h = psd.taper(32, "hann")
        self.assertAlmostEqual(h[0], 0.0, places=12)
        self.assertGreater(h[16], 0.99)
        t0 = psd.taper(32, "tukey", alpha=0.0)
        self.assertTrue(np.allclose(t0, 1.0))  # alpha = 0 is a boxcar
        with self.assertRaises(ValueError):
            psd.taper(32, "not-a-window")

    def test_radial_abscissa_is_the_mean_of_the_annulus(self):
        rng = np.random.default_rng(3)
        r = psd.radial_average(psd.welch_psd2d(rng.normal(size=(128, 128)), 0.1, roi=128))
        self.assertTrue(np.all(np.diff(r.q) > 0))
        self.assertTrue(np.all(r.n >= 1))
        # mode counts grow with the circumference
        self.assertGreater(r.n[len(r.n) // 2], r.n[2])

    def test_pixel_mtf_is_anisotropic_and_unity_at_dc(self):
        q = np.linspace(0.0, 30.0, 64)
        m2 = psd.pixel_mtf2(q, q, 0.1)
        self.assertAlmostEqual(m2[0, 0], 1.0, places=12)
        # on-axis and diagonal differ near Nyquist: that is why it must be
        # divided out in 2-D, before any radial averaging
        n = len(q) - 1
        self.assertGreater(abs(m2[0, n] - m2[n, n]), 0.05)

    def test_deconvolving_the_pixel_mtf_raises_high_q_power(self):
        rng = np.random.default_rng(4)
        raw = psd.welch_psd2d(rng.normal(size=(128, 128)), 0.1, roi=128, window="boxcar")
        fixed = psd.deconvolve_pixel_mtf(raw)
        a = psd.radial_average(raw)
        b = psd.radial_average(fixed)
        self.assertGreater(b.P[-1], a.P[-1])
        self.assertAlmostEqual(b.P[0] / a.P[0], 1.0, delta=0.05)

    def test_gaussian_mtf_deconvolution_inverts_a_known_blur(self):
        rng = np.random.default_rng(5)
        p, sigma = 0.05, 0.12
        x = rng.normal(size=(256, 256))
        raw = psd.welch_psd2d(x, p, roi=256, window="boxcar")
        q2 = raw.qy[:, None] ** 2 + raw.qx[None, :] ** 2
        blurred = psd.Psd2D(raw.qx, raw.qy, raw.P * np.exp(-q2 * sigma**2), raw.p, raw.n_avg)
        back = psd.radial_average(psd.deconvolve_gaussian_mtf(blurred, sigma))
        ref = psd.radial_average(raw)
        band = ref.q < 0.5 * ref.q.max()
        self.assertLess(np.max(np.abs(back.P[band] / ref.P[band] - 1.0)), 1e-9)


class TestPenetrableSpheres(unittest.TestCase):
    """The primary standard: S(q) = 1 exactly, so P_X has no free parameters."""

    def test_projected_spectrum_matches_the_closed_form(self):
        L, d, phi = 16.0, 1.0, 0.20
        centres, radii, n = poisson_spheres(phi=phi, L=L, d=d, seed=7)
        X, p = project(centres, radii, L, 256)

        self.assertAlmostEqual(X.mean() / (phi * L), 1.0, delta=0.01)

        r = psd.radial_average(psd.deconvolve_pixel_mtf(psd.welch_psd2d(X, p, roi=256, window="boxcar")))
        v_s = units.sphere_volume(d)
        analytic = (n / L**3) * L * v_s**2 * units.form_factor(r.q * d / 2) ** 2
        band = (r.q > 1.0) & (r.q < 12.0)
        self.assertAlmostEqual(np.median((r.P / analytic)[band]), 1.0, delta=0.06)


if __name__ == "__main__":
    unittest.main()
