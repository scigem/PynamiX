"""Uncertainty machinery: coverage, budgets, bounds, and the calibration fit."""

import unittest

import numpy as np

from pynamix.spectral.calibrate import fit_calibration, infer_H
from pynamix.spectral.invert import SizeDistribution, fit_first_peak_ml, schulz_pdf
from pynamix.spectral.uncertainty import (
    Budget,
    coverage,
    crb_displacement,
    crb_first_peak,
    efficiency,
    modes_in_annulus,
)


class TestCoverage(unittest.TestCase):
    def test_on_synthetic_gaussian_errors(self):
        """The coverage test must itself be calibrated before it can test anything."""
        rng = np.random.default_rng(0)
        n = 20000
        s = np.full(n, 0.02)
        c1 = coverage(rng.normal(0, 0.02, n), np.zeros(n), s, k=1.0)
        c2 = coverage(rng.normal(0, 0.02, n), np.zeros(n), s, k=2.0)
        self.assertAlmostEqual(c1["covered"], 0.683, delta=0.02)
        self.assertAlmostEqual(c2["covered"], 0.954, delta=0.01)

    def test_detects_an_understated_uncertainty(self):
        rng = np.random.default_rng(1)
        n = 20000
        c = coverage(rng.normal(0, 0.04, n), np.zeros(n), np.full(n, 0.02))
        self.assertLess(c["covered"], 0.45)
        self.assertLess(c["z"], -10)


class TestBoundsAndBudgets(unittest.TestCase):
    def test_mode_count_scales_with_field_size_not_pixels(self):
        self.assertAlmostEqual(modes_in_annulus(7.4, 32.0), 7.4 * 32.0, places=9)
        self.assertAlmostEqual(modes_in_annulus(7.4, 64.0), 2 * modes_in_annulus(7.4, 32.0), places=9)

    def test_budget_combines_and_correlation_matters(self):
        b = Budget("H", 24.0).add("X", 14.4, 0.02, 1 / 0.6, "A").add("phi", 0.6, 0.005, -24 / 0.6, "A")
        indep = b.combined
        b.correlations[("X", "phi")] = 0.7
        self.assertGreater(b.combined, indep)
        self.assertAlmostEqual(b.expanded(2.0), 2 * b.combined, places=12)
        self.assertEqual(b.by_kind("B"), 0.0)
        self.assertTrue(b.table().contribution.is_monotonic_decreasing)

    def test_crb_is_a_lower_bound_and_scales_with_modes(self):
        q = np.linspace(3.0, 11.0, 120)
        S = 1 + 3.0 * np.exp(-((q - 7.4) ** 2) / (2 * 0.55**2))
        n1 = np.full_like(q, 200.0)
        u1 = crb_first_peak(q, S, n1)
        u4 = crb_first_peak(q, S, 4 * n1)
        self.assertAlmostEqual(u4 / (u1 / 2), 1.0, places=6)
        self.assertAlmostEqual(efficiency(2 * u1, u1), 0.25, places=12)

    def test_crb_displacement_falls_with_wavevector(self):
        qx = np.full(500, 5.0)
        qy = np.zeros(500)
        g = np.full(500, 0.9)
        self.assertAlmostEqual(crb_displacement(2 * qx, qy, g) / (crb_displacement(qx, qy, g) / 2), 1.0, places=9)


class TestCalibration(unittest.TestCase):
    def test_is_a_fit_not_an_interpolant(self):
        """A fit must not pass exactly through noisy points."""
        rng = np.random.default_rng(3)
        phi = np.repeat([0.50, 0.55, 0.60, 0.635], 5)
        Q1 = 4.0 * phi + 5.1 + rng.normal(0, 0.02, len(phi))
        cal = fit_calibration(phi, Q1, deg=1)
        self.assertAlmostEqual(cal.sensitivity(0.55), 4.0, delta=0.3)
        self.assertTrue(cal.is_monotonic)
        self.assertFalse(np.allclose(np.atleast_1d(cal(phi)), Q1))
        # the curve is better determined than any single point on it
        self.assertLess(cal.curve_uncertainty(0.575), cal.u_scatter)

    def test_inversion_propagates_all_three_terms(self):
        phi = np.repeat([0.50, 0.55, 0.60, 0.635], 4)
        cal = fit_calibration(phi, 4.0 * phi + 5.1, deg=1)
        _, u_bare = cal.invert(7.3, 0.0, 0.0)
        _, u_meas = cal.invert(7.3, 0.02, 0.0)
        _, u_both = cal.invert(7.3, 0.02, 0.02)
        self.assertLess(u_bare, u_meas)
        self.assertLess(u_meas, u_both)

    def test_infer_H_error_propagation(self):
        H, u = infer_H(14.4, 0.6, 0.02, 0.005)
        self.assertAlmostEqual(H, 24.0, places=9)
        self.assertAlmostEqual(u / H, np.hypot(0.02 / 14.4, 0.005 / 0.6), places=9)


class TestDistributionAndPeakFit(unittest.TestCase):
    def test_schulz_moments_match_the_closed_form(self):
        for cv in (0.05, 0.15, 0.30):
            with self.subTest(cv=cv):
                d = np.linspace(1e-4, 6.0, 400000)
                p = schulz_pdf(d, 1.0, cv)
                p = p / np.trapezoid(p, d)
                mean = np.trapezoid(p * d, d)
                sd = np.sqrt(np.trapezoid(p * (d - mean) ** 2, d))
                d32 = np.trapezoid(p * d**3, d) / np.trapezoid(p * d**2, d)
                self.assertAlmostEqual(mean, 1.0, delta=1e-4)
                self.assertAlmostEqual(sd / mean, cv, delta=2e-3 * cv)
                self.assertAlmostEqual(SizeDistribution(0.6, mean, cv, 0.0).d32, d32, delta=2e-3 * d32)

    def test_ml_peak_fit_is_unbiased_on_a_known_peak(self):
        q = np.arange(0.2, 12.0, 0.196)
        for true in (7.05, 7.30, 7.55):
            with self.subTest(q1=true):
                S = 1 + 3.0 * np.exp(-((q - true) ** 2) / (2 * 0.55**2))
                fit = fit_first_peak_ml(q, S, np.full_like(q, 200.0), search=(4.0, 11.0))
                self.assertAlmostEqual(fit.q1, true, delta=0.02)
                self.assertTrue(np.isfinite(fit.sigma_q1))
                self.assertGreater(fit.sigma_q1, 0)


if __name__ == "__main__":
    unittest.main()
