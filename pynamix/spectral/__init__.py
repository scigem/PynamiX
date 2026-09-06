"""Quantitative spectral analysis of radiographs of dense granular media.

This subpackage is the model-based counterpart to :mod:`pynamix.measure`.
Both read the same Fourier transform of the same radiograph; they differ in
what they do with it.

:mod:`pynamix.measure` extracts an empirical descriptor -- the wavelength of a
spectral peak, a ratio of two peak heights -- and calibrates it against known
samples.  :mod:`pynamix.spectral` inverts a forward model,

.. math:: P_X(q) = \\langle X \\rangle\\, v_s\\, F^2(qd/2)\\, S(q),

so the same spectrum yields the solid fraction, the sample thickness, the
grain size and a velocity per species, two of them with no calibration curve at
all.

**The conventions differ from** :mod:`pynamix.measure` **and are not
interchangeable.** Three in particular:

``patchw`` versus ``roi``
    :mod:`pynamix.measure` takes ``patchw``, a *half* width.  Everything here
    takes a full width in pixels.

annular sum versus radial mean
    :func:`pynamix.measure.radial_FFT` *sums* over each annulus, which carries a
    factor of :math:`2\\pi q` -- the Jacobian -- and shifts the apparent peak
    towards shorter wavelengths.  :func:`~pynamix.spectral.psd.radial_average`
    takes the arithmetic *mean*, because each bin is an independent sample of
    the same spectral density and the mean is what the forward model predicts.

normalisation
    :func:`pynamix.exposure.mean_std` divides each patch by its own standard
    deviation.  That discards the absolute amplitude, which is precisely what
    fixes the solid fraction and the species volume fractions without a
    calibration.  Estimators here need an *unnormalised* projected path length
    in physical units.  Passing them a ``mean_std``-normalised patch returns a
    confident, meaningless number; there is no way to detect it after the fact.

See :mod:`pynamix.spectral.compat` for the bridge between the two.
"""

from . import (calibrate, closure, compat, detector, invert, psd, uncertainty,
               units, velocity)
from .calibrate import Calibration, fit_calibration, infer_H
from .invert import (
    concentration_from_porod,
    concentration_from_self_term,
    fit_first_peak,
    fit_first_peak_ml,
    fit_self_term,
    fit_size_distribution,
    porod_fit,
    profile_d32,
    recover_S,
)
from .psd import (
    deconvolve_gaussian_mtf,
    deconvolve_pixel_mtf,
    periodogram,
    radial_average,
    welch_psd2d,
)
from .units import form_factor, sauter_diameter, specific_surface, sphere_volume
from .velocity import (
    cross_spectrum,
    fit_species_displacements,
    species_weights,
    suggest_band,
    welch_cross_spectrum,
)

__all__ = [
    "calibrate", "closure", "detector", "invert", "psd", "uncertainty", "units",
    "velocity", "compat",
    "Calibration", "fit_calibration", "infer_H",
    "recover_S", "fit_first_peak", "fit_first_peak_ml", "fit_self_term",
    "fit_size_distribution", "porod_fit", "profile_d32",
    "concentration_from_porod", "concentration_from_self_term",
    "periodogram", "welch_psd2d", "radial_average",
    "deconvolve_pixel_mtf", "deconvolve_gaussian_mtf",
    "form_factor", "sphere_volume", "specific_surface", "sauter_diameter",
    "cross_spectrum", "welch_cross_spectrum", "species_weights",
    "suggest_band", "fit_species_displacements",
]
