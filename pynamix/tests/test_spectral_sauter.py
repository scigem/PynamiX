"""The Sauter mean diameter from the Porod tail, and what it must not do.

The first test is a regression guard: an earlier version applied the
zero-avoiding mask of ``recover_S`` inside ``porod_fit``.  That mask is correct
where one *divides* by a vanishing form factor and wrong where one *averages*,
because it deletes the troughs of the oscillation and keeps the peaks.  It
biased ``s`` high by about 12 per cent for a bidisperse packing, and so
``d_32`` low by the same amount, with no other symptom.
"""

import unittest

import numpy as np

from pynamix.spectral.invert import _chi_model, fit_size_distribution, porod_fit
from pynamix.spectral.psd import RadialPsd
from pynamix.spectral.units import form_factor, sphere_volume


def _exact_self_term(q, diameters, counts, volume):
    """chi(q) = (1/V) sum_j v_j^2 F^2(q a_j), the exact self term."""
    d = np.repeat(np.asarray(diameters, float), counts)
    v = sphere_volume(d)
    return np.array([float((v**2 * form_factor(qq * d / 2.0) ** 2).sum() / volume) for qq in q])


class TestPorodSauter(unittest.TestCase):
    def test_recovers_sauter_for_a_bidisperse_mixture(self):
        for R, c in ((1.5, 0.5), (2.0, 0.5), (2.0, 0.25), (3.0, 0.5)):
            with self.subTest(R=R, c=c):
                d_s, d_l, phi, V = 1.0, R, 0.60, 1.0e5
                n_l = c * phi * V / sphere_volume(d_l)
                n_s = (1.0 - c) * phi * V / sphere_volume(d_s)
                counts = [int(round(n_s)), int(round(n_l))]
                q = np.linspace(30.0, 60.0, 400)
                chi = _exact_self_term(q, [d_s, d_l], counts, V)
                psd = RadialPsd(q=q, P=chi, n=np.full(q.size, 300.0), sem=chi / np.sqrt(300.0))

                d = np.repeat([d_s, d_l], counts)
                d32_true = float((d**3).sum() / (d**2).sum())
                d32 = 6.0 * (sphere_volume(d).sum() / V) / porod_fit(psd, 1.0, (30.0, 60.0))
                # 1-2% is the honest finite-band residual; masked it gave 12%
                self.assertLess(abs(d32 / d32_true - 1.0), 0.04)

    def test_porod_beats_the_shape_fit_for_a_broad_distribution(self):
        """Once the oscillations are damped, the two shape parameters are
        unidentifiable but the Porod amplitude is not."""
        for cv in (0.10, 0.20, 0.30, 0.40):
            with self.subTest(cv=cv):
                q = np.linspace(30.0, 60.0, 153)
                n = np.full(q.size, 300.0)
                d_mean, phi = 1.0, 0.60
                d32_true = d_mean * (1.0 + 2.0 * cv**2)
                num, vbar = _chi_model(q, d_mean, cv)
                chi = phi * num / vbar
                psd = RadialPsd(q=q, P=chi, n=n, sem=chi / np.sqrt(n))

                d32_porod = 6.0 * phi / porod_fit(psd, 1.0, (30.0, 60.0))
                self.assertLess(abs(d32_porod / d32_true - 1.0), 0.02)

                shape = fit_size_distribution(psd, 1.0, multistart=True)
                # documents the limitation rather than asserting the fit works
                self.assertGreater(abs(shape.d32 / d32_true - 1.0), abs(d32_porod / d32_true - 1.0))

    def test_multistart_is_the_default_and_does_not_need_a_seed(self):
        q = np.linspace(30.0, 60.0, 153)
        n = np.full(q.size, 300.0)
        num, vbar = _chi_model(q, 1.0, 0.05)
        chi = 0.60 * num / vbar
        psd = RadialPsd(q=q, P=chi, n=n, sem=chi / np.sqrt(n))
        # narrow distribution: the oscillations survive, so the shape fit is sound
        out = fit_size_distribution(psd, 1.0)
        self.assertAlmostEqual(out.d32 / (1.0 + 2.0 * 0.05**2), 1.0, delta=0.02)


if __name__ == "__main__":
    unittest.main()
