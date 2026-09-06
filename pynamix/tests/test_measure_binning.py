"""Regression tests for the window and the radial/angular binning.

Each of these pins down something that was wrong before v0.36 and would not
have announced itself: an off-centre window, Monte Carlo coefficients that
differed between users, and radial bins labelled by their upper edge.
"""
import warnings

import numpy as np
import pytest

from pynamix.measure import (
    _check_wavelength_window,
    _parabolic_peak,
    angular_binning,
    hanning_window,
    radial_grid,
)


def test_window_is_centred_and_symmetric():
    w = hanning_window(32)
    assert w.shape == (64, 64)
    # a window centred on the patch is invariant under a 180 degree rotation
    assert np.abs(w - w[::-1, ::-1]).max() == pytest.approx(0.0, abs=1e-12)
    # peak at the centre, zero outside the inscribed circle
    assert w[31, 31] == pytest.approx(w.max())
    assert w[0, 0] == 0.0 and w[-1, -1] == 0.0


def test_window_is_a_raised_cosine_in_radius():
    patchw = 16
    w = hanning_window(patchw)
    i, j = np.indices(w.shape)
    d = np.sqrt((i - (patchw - 0.5)) ** 2 + (j - (patchw - 0.5)) ** 2)
    inside = d <= patchw
    expect = 0.5 * np.cos(np.pi * d[inside] / patchw) + 0.5
    assert np.allclose(w[inside], expect)


def test_angular_binning_is_deterministic_and_exact():
    a, b = angular_binning(16), angular_binning(16)
    assert np.array_equal(a, b)  # no unseeded Monte Carlo
    # k-hat is a unit vector, so the tensor has unit trace at every pixel
    assert np.abs(a[:, :, 0, 0] + a[:, :, 1, 1] - 1.0).max() < 1e-12
    # and it is symmetric, exactly
    assert np.array_equal(a[:, :, 0, 1], a[:, :, 1, 0])


def test_angular_binning_matches_the_analytic_limit_far_from_the_origin():
    patchw = 16
    Q = angular_binning(patchw)
    # a pixel many pixels from the origin subtends a tiny angle, so the average
    # of k-hat k-hat over it is just k-hat k-hat at its centre
    row, col = 2, 3
    y, x = row - patchw + 0.5, col - patchw + 0.5
    n = np.hypot(x, y)
    assert Q[row, col, 0, 0] == pytest.approx((x / n) ** 2, abs=2e-3)
    assert Q[row, col, 1, 1] == pytest.approx((y / n) ** 2, abs=2e-3)


def test_radial_grid_partitions_every_pixel_exactly():
    r_grid, nr_pxr = radial_grid(rnb=64, patchw=16)
    # every pixel's area is fully accounted for across the bins
    assert np.abs(nr_pxr.sum(axis=2) - 1.0).max() < 1e-12
    assert np.all(np.diff(r_grid) > 0)
    assert r_grid[0] >= 0.0


def test_radial_grid_labels_bins_by_their_mean_radius_not_their_edge():
    """The old code shifted r_grid up by half a bin after binning."""
    rnb, patchw = 200, 32
    r_grid, nr_pxr = radial_grid(rnb=rnb, patchw=patchw)
    nodes = np.linspace(0, patchw * 1.5, rnb)
    delta = nodes[1] - nodes[0]
    occupied = nr_pxr.sum(axis=(0, 1)) > 0
    # each label sits within half a bin of its node, not systematically above it
    offset = (r_grid - nodes)[occupied]
    assert np.abs(offset).max() < 0.75 * delta
    assert abs(np.median(offset)) < 0.25 * delta


def test_radial_grid_is_deterministic():
    a, _ = radial_grid(rnb=32, patchw=8)
    b, _ = radial_grid(rnb=32, patchw=8)
    assert np.array_equal(a, b)


def test_parabolic_peak_finds_a_known_vertex():
    x = np.linspace(0.0, 4.0, 5)
    # y = -(x - 2.3)^2 + 5, sampled on the grid
    y = (-((x - 2.3) ** 2) + 5.0)[None, :]
    i = np.array([np.argmax(y[0])])
    assert _parabolic_peak(x, y, i)[0] == pytest.approx(2.3, abs=1e-9)


def test_parabolic_peak_falls_back_when_there_is_no_maximum():
    x = np.linspace(0.0, 4.0, 5)
    y = np.array([[0.0, 1.0, 2.0, 3.0, 4.0]])  # monotonic: no interior vertex
    i = np.array([2])
    assert _parabolic_peak(x, y, i)[0] == pytest.approx(x[2])


def test_wavelength_window_guard_fires_near_dc():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _check_wavelength_window(20.0, 32, 1.0, "test")   # 3.2 bins from DC
    assert any("bins from DC" in str(c.message) for c in caught)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _check_wavelength_window(6.0, 32, 1.0, "test")    # 10.7 bins from DC
    assert not caught
