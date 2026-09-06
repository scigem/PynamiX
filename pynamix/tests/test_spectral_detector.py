"""The detector chain: attenuation, noise, hardening, scatter and their guards."""

import unittest

import numpy as np

from pynamix.spectral import detector


class TestSpectrum(unittest.TestCase):
    def test_monochromatic_has_no_hardness(self):
        s = detector.Spectrum.monochromatic(0.3)
        self.assertAlmostEqual(s.mu_mean, 0.3, places=12)
        self.assertAlmostEqual(s.hardness, 0.0, places=12)

    def test_three_point_hits_the_requested_mean_and_spread(self):
        for h in (0.05, 0.2, 0.4):
            with self.subTest(hardness=h):
                s = detector.Spectrum.three_point(0.3, h)
                self.assertAlmostEqual(s.mu_mean, 0.3, places=12)
                self.assertAlmostEqual(s.hardness, h, places=10)
                self.assertAlmostEqual(float(s.weight.sum()), 1.0, places=12)


class TestAttenuation(unittest.TestCase):
    def test_monochromatic_transmission_is_beer_lambert(self):
        X = np.linspace(0.0, 8.0, 17)
        s = detector.Spectrum.monochromatic(0.25)
        self.assertTrue(np.allclose(detector.polychromatic_transmission(X, s), np.exp(-0.25 * X)))

    def test_a_hard_beam_transmits_more_than_beer_lambert(self):
        """Beam hardening: the soft end burns off, so the beam gets more penetrating."""
        X = np.linspace(0.5, 8.0, 16)
        soft = detector.polychromatic_transmission(X, detector.Spectrum.monochromatic(0.25))
        hard = detector.polychromatic_transmission(X, detector.Spectrum.three_point(0.25, 0.35))
        self.assertTrue(np.all(hard >= soft - 1e-12))
        self.assertGreater(hard[-1] / soft[-1], 1.05)

    def test_hardening_calibration_inverts_the_curve(self):
        s = detector.Spectrum.three_point(0.25, 0.3)
        curve = detector.hardening_curve(s, X_max=10.0)
        X = np.linspace(0.2, 9.0, 40)
        recovered = detector.linearise_calibrated(detector.polychromatic_transmission(X, s), curve)
        self.assertLess(np.max(np.abs(recovered - X)), 0.02)

    def test_expected_counts_and_linearise_round_trip(self):
        X = np.linspace(0.0, 6.0, 25)
        counts = detector.expected_counts(X, 0.3, 1e5)
        self.assertLess(np.max(np.abs(detector.linearise(counts, 0.3, 1e5) - X)), 1e-9)


class TestBlurAndScatter(unittest.TestCase):
    def test_blur_conserves_the_mean_and_zero_width_is_a_no_op(self):
        rng = np.random.default_rng(0)
        img = 1.0 + 0.1 * rng.normal(size=(64, 64))
        blurred = detector.apply_mtf(img, 1.5)
        self.assertAlmostEqual(blurred.mean(), img.mean(), places=10)
        self.assertLess(blurred.std(), img.std())
        self.assertTrue(np.allclose(detector.apply_mtf(img, 0.0), img))

    def test_scatter_is_an_additive_pedestal(self):
        img = np.full((16, 16), 2.0)
        out = detector.add_scatter(img, 0.1)
        self.assertAlmostEqual(out.mean(), 2.2, places=12)
        self.assertTrue(np.allclose(detector.add_scatter(img, 0.0), img))

    def test_gaussian_mtf2_is_unity_at_dc_and_falls_with_q(self):
        q = np.linspace(0.0, 20.0, 32)
        m2 = detector.gaussian_mtf2(q, q, 0.1)
        self.assertAlmostEqual(m2[0, 0], 1.0, places=12)
        self.assertLess(m2[-1, -1], m2[0, 0])


class TestNoiseModel(unittest.TestCase):
    def test_analytic_variance_matches_a_simulation(self):
        mu, N0, X_mean = 0.3, 1e4, 6.0
        rng = np.random.default_rng(2)
        X = np.full((512, 512), X_mean)
        Y = detector.simulate(X, mu, N0, rng=rng)
        self.assertAlmostEqual(np.var(Y) / detector.noise_variance_X(mu, N0, X_mean), 1.0, delta=0.05)

    def test_noise_psd_is_the_variance_times_the_pixel_area(self):
        mu, N0, X_mean, p = 0.3, 1e4, 6.0, 0.05
        self.assertAlmostEqual(
            detector.noise_psd(mu, N0, X_mean, p) / (detector.noise_variance_X(mu, N0, X_mean) * p**2),
            1.0,
            places=12,
        )

    def test_read_noise_only_ever_adds_variance(self):
        a = detector.noise_variance_X(0.3, 1e4, 6.0, read_noise=0.0)
        b = detector.noise_variance_X(0.3, 1e4, 6.0, read_noise=5.0)
        self.assertGreater(b, a)

    def test_jensen_factor_is_unity_for_a_flat_field(self):
        X = np.full((32, 32), 6.0)
        self.assertAlmostEqual(detector.jensen_factor(X, 0.3, 1e4), 1.0, delta=1e-9)
        rng = np.random.default_rng(3)
        rough = 6.0 + 1.5 * rng.normal(size=(64, 64))
        # a spread in X raises the mean noise above its value at the mean X
        self.assertGreater(detector.jensen_factor(rough, 0.3, 1e4), 1.0)

    def test_linearisation_guard_refuses_a_starved_exposure(self):
        self.assertTrue(detector.linearisation_valid(0.3, 1e5, 6.0))
        self.assertFalse(detector.linearisation_valid(0.3, 1e1, 6.0))
        self.assertGreater(detector.linearisation_snr(0.3, 1e5, 6.0), detector.linearisation_snr(0.3, 1e2, 6.0))

    def test_log_bias_falls_as_the_exposure_rises(self):
        self.assertGreater(detector.log_bias(0.3, 1e2, 6.0), detector.log_bias(0.3, 1e5, 6.0))

    def test_optimal_operating_point_is_near_two(self):
        """Photon-limited SNR peaks at mu*X = 2; read noise pulls it down."""
        self.assertAlmostEqual(detector.optimal_mu_X(), 2.0, delta=0.05)
        self.assertLess(detector.optimal_mu_X(read_noise=20.0, N0=1e3), 2.01)


if __name__ == "__main__":
    unittest.main()
