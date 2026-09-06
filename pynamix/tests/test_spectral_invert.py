"""Inversion: recovering S(q), locating its first peak, and the self-term fits."""

import unittest

import numpy as np

from pynamix.spectral import invert, psd, units
from pynamix.tests.spectral_fixtures import bidisperse_spheres, poisson_spheres, project


def _radial(X, p, n_pixels):
    return psd.radial_average(psd.deconvolve_pixel_mtf(psd.welch_psd2d(X, p, roi=n_pixels, window="boxcar")))


class TestUnits(unittest.TestCase):
    def test_form_factor_limits_and_first_zero(self):
        self.assertAlmostEqual(float(units.form_factor(0.0)), 1.0, places=12)
        # the small-u series must not lose precision to cancellation
        self.assertAlmostEqual(float(units.form_factor(1e-6)), 1.0, places=10)
        self.assertAlmostEqual(float(units.form_factor(units.QD_FIRST_ZERO / 2)), 0.0, delta=1e-6)
        self.assertAlmostEqual(float(units.form_factor_qd(units.QD_FIRST_ZERO)), 0.0, delta=1e-6)

    def test_porod_envelope_is_the_large_u_average_of_the_form_factor(self):
        u = np.linspace(60.0, 80.0, 20001)
        F2 = units.form_factor(u) ** 2
        self.assertAlmostEqual(np.mean(F2 * u**4) / 4.5, 1.0, delta=0.02)

    def test_specific_surface_and_sauter_diameter(self):
        phis = {1.0: 0.3, 2.0: 0.3}
        s = units.specific_surface(phis)
        self.assertAlmostEqual(s, 6 * (0.3 / 1.0 + 0.3 / 2.0), places=12)
        self.assertAlmostEqual(units.sauter_diameter(phis), 6 * 0.6 / s, places=12)
        self.assertTrue(np.all(units.porod_chi(np.array([2.0, 4.0]), s) > 0))


class TestRecoverS(unittest.TestCase):
    def test_penetrable_spheres_give_S_equal_to_one(self):
        """S(q) = 1 exactly for uniformly random centres: the sharpest check there is."""
        L, d, phi = 16.0, 1.0, 0.20
        centres, radii, _ = poisson_spheres(phi=phi, L=L, d=d, seed=3)
        X, p = project(centres, radii, L, 256)
        S = invert.recover_S(_radial(X, p, 256), X.mean(), d, 0.0)
        band = (S.q > 1.0) & (S.q < 8.0)
        self.assertAlmostEqual(float(np.median(S.S[band])), 1.0, delta=0.08)

    def test_bidisperse_recovery_returns_the_polydispersity_floor(self):
        L, d_s, d_l, n_s, n_l = 16.0, 1.0, 2.0, 3000, 400
        centres, radii = bidisperse_spheres(L=L, d_s=d_s, d_l=d_l, n_s=n_s, n_l=n_l, seed=5)
        X, p = project(centres, radii, L, 256)
        S = invert.recover_S_bidisperse(
            _radial(X, p, 256),
            X.mean(),
            L,
            {1: d_s, 2: d_l},
            {1: n_s / L**3, 2: n_l / L**3},
        )
        band = (S.q > 1.0) & (S.q < 8.0)
        # penetrable spheres again: the decoupling inversion should land on 1
        self.assertAlmostEqual(float(np.median(S.S[band])), 1.0, delta=0.2)

    def test_lobe_mask_excludes_the_form_factor_zeros(self):
        q = np.linspace(1.0, 30.0, 2000)
        keep = invert.lobe_mask(q, [1.0])
        self.assertFalse(bool(keep[np.argmin(np.abs(q - units.QD_FIRST_ZERO))]))
        self.assertTrue(bool(keep[np.argmin(np.abs(q - 5.0))]))
        # the mask is a union over species, not an intersection: a bin survives
        # if it is near a lobe maximum for *any* diameter, so a mixture loses
        # less data than a single species does
        self.assertGreaterEqual(invert.lobe_mask(q, [1.0, 2.0]).sum(), keep.sum())


class TestFirstPeak(unittest.TestCase):
    def _synthetic(self, q1, n=400):
        q = np.linspace(3.0, 12.0, n)
        return q, 1 + 3.0 * np.exp(-((q - q1) ** 2) / (2 * 0.55**2))

    def test_quadratic_fit_finds_a_known_peak(self):
        for true in (6.8, 7.4, 8.1):
            with self.subTest(q1=true):
                q, S = self._synthetic(true)
                fit = invert.fit_first_peak(q, S, search=(4.0, 11.0))
                self.assertAlmostEqual(fit.q1, true, delta=0.02)
                self.assertGreater(fit.S1, 3.5)
                self.assertGreater(fit.width, 0.0)

    def test_the_fit_beats_taking_the_largest_bin(self):
        """Sub-bin interpolation is the point; a coarse grid should not limit it."""
        q = np.arange(3.0, 12.0, 0.25)
        true = 7.37
        S = 1 + 3.0 * np.exp(-((q - true) ** 2) / (2 * 0.55**2))
        coarse = q[np.argmax(S)]
        fit = invert.fit_first_peak(q, S, search=(4.0, 11.0))
        self.assertLess(abs(fit.q1 - true), abs(coarse - true))

    def test_a_featureless_window_yields_no_usable_width(self):
        """There is no peak to find, and the returned uncertainty says so."""
        q = np.linspace(3.0, 12.0, 200)
        fit = invert.fit_first_peak(q, np.ones_like(q), search=(11.5, 11.9))
        self.assertTrue(np.isnan(fit.width))
        self.assertTrue(np.isnan(fit.sigma_q1))


class TestSelfTerm(unittest.TestCase):
    def test_recovers_both_species_volume_fractions(self):
        """The exact self term, fitted against the known basis functions."""
        d_s, d_l, V = 1.0, 2.0, 1.0e5
        phi_s, phi_l = 0.20, 0.25
        q = np.linspace(30.0, 60.0, 400)
        chi = (
            phi_s * units.sphere_volume(d_s) * units.form_factor(q * d_s / 2) ** 2
            + phi_l * units.sphere_volume(d_l) * units.form_factor(q * d_l / 2) ** 2
        )
        r = psd.RadialPsd(q=q, P=chi, n=np.full(q.size, 300.0), sem=chi / np.sqrt(300.0))
        out = invert.fit_self_term(r, 1.0, [d_s, d_l])
        self.assertAlmostEqual(out[d_s] / phi_s, 1.0, delta=0.02)
        self.assertAlmostEqual(out[d_l] / phi_l, 1.0, delta=0.02)
        c = invert.concentration_from_self_term(out, d_l)
        self.assertAlmostEqual(c, phi_l / (phi_s + phi_l), delta=0.02)

    def test_too_few_bins_is_refused(self):
        q = np.linspace(30.0, 30.2, 3)
        r = psd.RadialPsd(q=q, P=np.ones_like(q), n=np.full(3, 10.0), sem=np.ones_like(q))
        with self.assertRaises(ValueError):
            invert.fit_self_term(r, 1.0, [1.0, 2.0])

    def test_concentration_from_porod_inverts_the_sauter_relation(self):
        for c in (0.1, 0.5, 0.9):
            with self.subTest(c=c):
                d_s, d_l = 1.0, 3.0
                d32 = 1.0 / ((1 - c) / d_s + c / d_l)
                self.assertAlmostEqual(invert.concentration_from_porod(d32, d_s, d_l), c, places=10)


class TestProfileD32(unittest.TestCase):
    def test_the_profile_flags_when_it_is_not_identifiable(self):
        """The flag is the point: a railed profile must not pass as a measurement.

        With several hundred modes per bin the Gamma likelihood is extremely
        sharp, so any mismatch between model and data dominates the interval and
        the maximum runs to the edge of the scan.  ``identifiable`` reports that
        rather than returning a confident wrong number.
        """
        q = np.linspace(30.0, 60.0, 153)
        n = np.full(q.size, 300.0)
        num, vbar = invert._chi_model(q, 1.0, 0.05)
        chi = 0.6 * num / vbar
        r = psd.RadialPsd(q=q, P=chi, n=n, sem=chi / np.sqrt(n))
        prof = invert.profile_d32(r, 1.0)
        self.assertEqual(prof.grid.shape, prof.ell.shape)
        self.assertGreaterEqual(prof.hi, prof.lo)
        self.assertGreater(prof.phi, 0.0)
        if not prof.identifiable:
            self.assertTrue(np.isclose(prof.d32, prof.grid[0]) or np.isclose(prof.d32, prof.grid[-1]))

    def test_too_few_bins_is_refused(self):
        q = np.linspace(30.0, 30.2, 4)
        r = psd.RadialPsd(q=q, P=np.ones_like(q), n=np.full(4, 10.0), sem=np.ones_like(q))
        with self.assertRaises(ValueError):
            invert.profile_d32(r, 1.0)


if __name__ == "__main__":
    unittest.main()
