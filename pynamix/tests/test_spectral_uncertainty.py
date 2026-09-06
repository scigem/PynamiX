"""Uncertainty machinery: coverage, budgets, bounds, and the calibration fit."""

import numpy as np
import pytest

from pynamix.spectral.calibrate import Calibration, fit_calibration, infer_H
from pynamix.spectral.invert import fit_first_peak_ml, schulz_pdf, SizeDistribution
from pynamix.spectral.uncertainty import (Budget, Contribution, coverage, crb_displacement,
                                 crb_first_peak, efficiency, modes_in_annulus)


def test_coverage_on_synthetic_gaussian_errors():
    """The coverage test must itself be calibrated before it can test anything."""
    rng = np.random.default_rng(0)
    n = 20000
    s = np.full(n, 0.02)
    c1 = coverage(rng.normal(0, 0.02, n), np.zeros(n), s, k=1.0)
    c2 = coverage(rng.normal(0, 0.02, n), np.zeros(n), s, k=2.0)
    assert c1["covered"] == pytest.approx(0.683, abs=0.02)
    assert c2["covered"] == pytest.approx(0.954, abs=0.01)


def test_coverage_detects_an_understated_uncertainty():
    rng = np.random.default_rng(1)
    n = 20000
    true_sigma, stated = 0.04, 0.02
    c = coverage(rng.normal(0, true_sigma, n), np.zeros(n), np.full(n, stated))
    assert c["covered"] < 0.45
    assert c["z"] < -10


def test_mode_count_scales_with_field_size_not_pixels():
    assert modes_in_annulus(7.4, 32.0) == pytest.approx(7.4 * 32.0)
    assert modes_in_annulus(7.4, 64.0) == 2 * modes_in_annulus(7.4, 32.0)


def test_budget_combines_and_correlation_matters():
    b = (Budget("H", 24.0)
         .add("X", 14.4, 0.02, 1 / 0.6, "A")
         .add("phi", 0.6, 0.005, -24 / 0.6, "A"))
    indep = b.combined
    b.correlations[("X", "phi")] = 0.7
    assert b.combined > indep
    assert b.expanded(2.0) == pytest.approx(2 * b.combined)
    assert b.by_kind("B") == 0.0
    assert b.table().contribution.is_monotonic_decreasing


def test_crb_is_a_lower_bound_and_scales_with_modes():
    q = np.linspace(3.0, 11.0, 120)
    S = 1 + 3.0 * np.exp(-((q - 7.4) ** 2) / (2 * 0.55**2))
    n1 = np.full_like(q, 200.0)
    u1 = crb_first_peak(q, S, n1)
    u4 = crb_first_peak(q, S, 4 * n1)
    assert u4 == pytest.approx(u1 / 2, rel=1e-6)
    assert efficiency(2 * u1, u1) == pytest.approx(0.25)


def test_crb_displacement_falls_with_wavevector():
    qx = np.full(500, 5.0)
    qy = np.zeros(500)
    g = np.full(500, 0.9)
    assert crb_displacement(2 * qx, qy, g) == pytest.approx(
        crb_displacement(qx, qy, g) / 2, rel=1e-9)


def test_calibration_is_a_fit_not_an_interpolant():
    """A fit must not pass exactly through noisy points."""
    rng = np.random.default_rng(3)
    phi = np.repeat([0.50, 0.55, 0.60, 0.635], 5)
    Q1 = 4.0 * phi + 5.1 + rng.normal(0, 0.02, len(phi))
    cal = fit_calibration(phi, Q1, deg=1)
    assert cal.sensitivity(0.55) == pytest.approx(4.0, abs=0.3)
    assert cal.is_monotonic
    assert not np.allclose(np.atleast_1d(cal(phi)), Q1)
    # the curve is better determined than any single point on it
    assert cal.curve_uncertainty(0.575) < cal.u_scatter


def test_calibration_inversion_propagates_all_three_terms():
    phi = np.repeat([0.50, 0.55, 0.60, 0.635], 4)
    Q1 = 4.0 * phi + 5.1
    cal = fit_calibration(phi, Q1, deg=1)
    _, u_bare = cal.invert(7.3, 0.0, 0.0)
    _, u_meas = cal.invert(7.3, 0.02, 0.0)
    _, u_both = cal.invert(7.3, 0.02, 0.02)
    assert u_bare < u_meas < u_both


def test_infer_H_error_propagation():
    H, u = infer_H(14.4, 0.6, 0.02, 0.005)
    assert H == pytest.approx(24.0)
    assert u / H == pytest.approx(np.hypot(0.02 / 14.4, 0.005 / 0.6), rel=1e-9)


@pytest.mark.parametrize("cv", [0.05, 0.15, 0.30])
def test_schulz_moments_match_the_closed_form(cv):
    d = np.linspace(1e-4, 6.0, 400000)
    p = schulz_pdf(d, 1.0, cv)
    p = p / np.trapezoid(p, d)
    mean = np.trapezoid(p * d, d)
    sd = np.sqrt(np.trapezoid(p * (d - mean) ** 2, d))
    d32 = np.trapezoid(p * d**3, d) / np.trapezoid(p * d**2, d)
    assert mean == pytest.approx(1.0, rel=1e-4)
    assert sd / mean == pytest.approx(cv, rel=2e-3)
    assert SizeDistribution(0.6, mean, cv, 0.0).d32 == pytest.approx(d32, rel=2e-3)


def test_ml_peak_fit_is_unbiased_on_a_known_peak():
    q = np.arange(0.2, 12.0, 0.196)
    for true in (7.05, 7.30, 7.55):
        S = 1 + 3.0 * np.exp(-((q - true) ** 2) / (2 * 0.55**2))
        fit = fit_first_peak_ml(q, S, np.full_like(q, 200.0), search=(4.0, 11.0))
        assert fit.q1 == pytest.approx(true, abs=0.02)
        assert np.isfinite(fit.sigma_q1) and fit.sigma_q1 > 0
