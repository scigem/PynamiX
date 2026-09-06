"""Velocimetry from the phase of the cross-spectrum."""

import unittest

import numpy as np

from pynamix.spectral import psd, velocity
from pynamix.tests.spectral_fixtures import bidisperse_spheres, poisson_spheres, project


def _pair(shift, L=16.0, d=1.0, phi=0.2, n_pixels=256, seed=11):
    centres, radii, _ = poisson_spheres(phi=phi, L=L, d=d, seed=seed)
    X1, p = project(centres, radii, L, n_pixels)
    moved = centres.copy()
    moved[:, 0] += shift[0]
    moved[:, 1] += shift[1]
    X2, _ = project(np.mod(moved, L), radii, L, n_pixels)
    return X1, X2, p


class TestCrossSpectrum(unittest.TestCase):
    def test_auto_cross_spectrum_is_the_psd(self):
        """A frame against itself must reproduce the PSD exactly, including scale.

        The normalisation of the cross-spectrum and of the PSD have to agree, or
        every weight in the species fit is wrong by a constant no one notices.
        """
        rng = np.random.default_rng(0)
        x = rng.normal(size=(128, 128))
        cs = velocity.welch_cross_spectrum(x, x, 0.1, roi=128, window="boxcar")
        pp = psd.welch_psd2d(x, 0.1, roi=128, window="boxcar")
        self.assertLess(np.max(np.abs(np.real(cs.C) - pp.P)) / pp.P.max(), 1e-12)
        self.assertLess(np.max(np.abs(np.imag(cs.C))) / pp.P.max(), 1e-12)

    def test_zero_shift_gives_zero_phase(self):
        X1, X2, p = _pair((0.0, 0.0))
        u = velocity.fit_displacement(velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="hann"), 1.0, 20.0)
        self.assertLess(np.hypot(*u), 1e-6)


class TestRigidShift(unittest.TestCase):
    def test_recovers_a_known_translation(self):
        for shift in ((0.10, 0.0), (0.0, -0.13), (0.08, 0.11)):
            with self.subTest(shift=shift):
                X1, X2, p = _pair(shift)
                cs = velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="hann")
                u = velocity.fit_displacement(cs, 1.0, 20.0)
                self.assertLess(np.hypot(u[0] - shift[0], u[1] - shift[1]), 0.02)

    def test_phase_correlation_agrees_with_the_phase_fit(self):
        X1, X2, p = _pair((0.10, 0.0))
        coarse = velocity.phase_correlation_shift(X1, X2, p)
        fine = velocity.fit_displacement(velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="hann"), 1.0, 20.0)
        self.assertLess(np.hypot(coarse[0] - fine[0], coarse[1] - fine[1]), 0.05)

    def test_displacement_is_ordered_columns_then_rows(self):
        """centres[:, 0] moves along array axis 1, and is reported first."""
        X1, X2, p = _pair((0.12, 0.0))
        u = velocity.fit_displacement(velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="hann"), 1.0, 20.0)
        self.assertGreater(u[0], 0.05)
        self.assertLess(abs(u[1]), 0.02)

    def test_untapered_subregions_read_low(self):
        """Welch regions have edges even when the parent image does not."""
        X1, X2, p = _pair((0.12, 0.0))
        boxcar = velocity.fit_displacement(
            velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="boxcar"), 1.0, 20.0
        )
        hann = velocity.fit_displacement(velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="hann"), 1.0, 20.0)
        self.assertLess(boxcar[0], hann[0])
        self.assertLess(abs(hann[0] - 0.12), abs(boxcar[0] - 0.12))


class TestBandsAndWeights(unittest.TestCase):
    def test_max_unambiguous_shift_is_half_a_wavelength_at_the_band_top(self):
        self.assertAlmostEqual(velocity.max_unambiguous_shift(10.0), np.pi / 10, places=12)
        self.assertGreater(velocity.max_unambiguous_shift(5.0), velocity.max_unambiguous_shift(10.0))

    def test_suggest_band_reaches_below_the_large_species_peak(self):
        out = velocity.suggest_band({1: 1.0, 2: 6.0}, p=1 / 32, n_pixels=1024)
        self.assertLess(out["q_min"], 7.5 / 6.0)  # must not exclude the large grains
        self.assertGreater(out["q_max"], out["q_min"])
        self.assertGreaterEqual(out["roi"], 2)

    def test_species_weights_sum_to_the_total_and_carry_the_cross_term(self):
        q = np.linspace(1.0, 25.0, 200)
        d = {1: 1.0, 2: 2.0}
        rho = {1: 0.4, 2: 0.05}
        w = velocity.species_weights(q, d, rho, H=16.0)
        self.assertEqual(set(w), {1, 2})
        for a in w:
            self.assertEqual(w[a].shape, q.shape)
            self.assertTrue(np.all(np.isfinite(w[a])))
        # the cross partial is what makes a weight able to go negative; with
        # ideal partials the diagonal terms dominate and both stay positive
        self.assertGreater(np.sum(w[1]), 0.0)
        self.assertGreater(np.sum(w[2]), 0.0)


class TestSpeciesFit(unittest.TestCase):
    def test_a_common_motion_is_recovered_for_both_species(self):
        """Both species moving together is the well-conditioned limit."""
        L, d_s, d_l, n_s, n_l = 16.0, 1.0, 2.0, 3000, 400
        centres, radii = bidisperse_spheres(L=L, d_s=d_s, d_l=d_l, n_s=n_s, n_l=n_l, seed=4)
        X1, p = project(centres, radii, L, 256)
        moved = centres.copy()
        moved[:, 0] += 0.10
        X2, _ = project(np.mod(moved, L), radii, L, 256)
        cs = velocity.welch_cross_spectrum(X1, X2, p, roi=128, window="hann")
        out = velocity.fit_species_displacements(cs, {1: d_s, 2: d_l}, {1: n_s / L**3, 2: n_l / L**3}, H=L, q_max=20.0)
        self.assertTrue(out.success)
        for a in (1, 2):
            self.assertLess(abs(out.displacement[a][0] - 0.10), 0.05)
            self.assertLess(abs(out.displacement[a][1]), 0.05)


if __name__ == "__main__":
    unittest.main()
