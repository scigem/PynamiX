"""Bridge between :mod:`pynamix.measure` conventions and :mod:`pynamix.spectral`.

The two halves of the package read the same Fourier transform of the same
radiograph but disagree on three conventions, and none of the disagreements
announces itself.  Everything here exists to make a crossing between them
explicit.
"""

import warnings

import numpy as np

__all__ = [
    "annular_sum_to_mean",
    "q_from_wavelength",
    "wavelength_from_q",
    "assert_unnormalised",
    "roi_from_patchw",
    "patchw_from_roi",
]


def roi_from_patchw(patchw):
    """Full ROI width in pixels from :mod:`pynamix.measure`'s half width."""
    return 2 * int(patchw)


def patchw_from_roi(roi):
    """:mod:`pynamix.measure`'s half width from a full ROI width in pixels."""
    if int(roi) % 2:
        raise ValueError(f"roi must be even to express as a half width, got {roi}")
    return int(roi) // 2


def q_from_wavelength(wavelength):
    r"""Angular wavevector :math:`q=2\pi/\lambda` from a wavelength axis.

    :func:`pynamix.measure.radial_FFT` returns a *wavelength* axis in physical
    units, sorted from long to short, with the first entry infinite (the DC
    bin).  Every estimator in :mod:`pynamix.spectral` is written in terms of the
    angular wavevector instead.  Note the :math:`2\pi`: a wavelength axis
    converts to *cycles* per unit length, not radians.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    with np.errstate(divide="ignore"):
        return 2.0 * np.pi / wavelength


def wavelength_from_q(q):
    r"""Wavelength :math:`\lambda=2\pi/q` from an angular wavevector axis."""
    q = np.asarray(q, dtype=float)
    with np.errstate(divide="ignore"):
        return 2.0 * np.pi / q


def annular_sum_to_mean(radialspec, nr_pxr):
    r"""Undo the annular sum of :func:`pynamix.measure.radial_FFT`.

    ``radial_FFT`` reports ``sum(S * nr_pxr[:, :, k])`` for each radial bin
    ``k``: an annular *sum*, which carries the :math:`2\pi q` Jacobian and so
    moves the apparent peak towards shorter wavelengths relative to the peak of
    the spectral density itself.  The forward model of
    :mod:`pynamix.spectral` predicts the density, so the sum has to be divided
    by the number of Fourier modes in the annulus, which is exactly
    ``nr_pxr[:, :, k].sum()``.

    Args:
        radialspec: the annular sums, with the radial bin as the last axis.
        nr_pxr: the binning weights from :func:`pynamix.measure.radial_grid`.

    Returns:
        The same array converted to an annular mean, and the per-bin mode count
        ``n_j`` that the likelihood-based estimators need.
    """
    radialspec = np.asarray(radialspec, dtype=float)
    n_j = np.asarray(nr_pxr, dtype=float).sum(axis=(0, 1))
    if radialspec.shape[-1] != n_j.shape[0]:
        raise ValueError(
            f"last axis of radialspec ({radialspec.shape[-1]}) must match the number "
            f"of radial bins in nr_pxr ({n_j.shape[0]})"
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(n_j > 0, radialspec / n_j, np.nan)
    return mean, n_j


def assert_unnormalised(patch, tol=0.05):
    """Refuse a patch whose absolute amplitude has been standardised away.

    :func:`pynamix.exposure.mean_std` divides a patch by its own standard
    deviation.  For the peak-position and orientation measures of
    :mod:`pynamix.measure` that is harmless, and it is the default there.  For
    every estimator in :mod:`pynamix.spectral` it is fatal: the solid fraction,
    the species volume fractions and the interfacial area are all read off the
    *absolute* spectral amplitude, and a standardised patch still produces a
    perfectly plausible number with no way to detect the error downstream.

    This checks the cheap signature -- unit variance about a zero mean -- and
    raises rather than warning, because a silent wrong answer is the failure
    mode being guarded against.  It cannot catch every case: a patch that
    happens to have unit variance in physical units is indistinguishable, hence
    the deliberately tight tolerance.

    Args:
        patch: the projected path length, in physical units of length.
        tol: how close to (mean 0, std 1) counts as standardised.
    """
    patch = np.asarray(patch, dtype=float)
    if patch.size == 0:
        raise ValueError("empty patch")
    mean, std = float(np.mean(patch)), float(np.std(patch))
    if abs(mean) < tol and abs(std - 1.0) < tol:
        raise ValueError(
            "patch looks standardised (mean %.3g, std %.3g): pynamix.spectral needs "
            "an unnormalised projected path length in physical units, because the "
            "absolute spectral amplitude is what fixes phi, the species volume "
            "fractions and the interfacial area. If this patch really is "
            "unnormalised and happens to have unit variance, pass it directly to "
            "the estimator instead of through this check." % (mean, std)
        )
    if std == 0.0:
        warnings.warn("patch has zero variance; the spectrum will be identically zero", RuntimeWarning, stacklevel=2)
