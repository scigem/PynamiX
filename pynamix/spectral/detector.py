"""Beer-Lambert attenuation, photon-counting noise, and the noise PSD floor.

The core model is photon (Poisson) statistics plus Gaussian read noise.  The
*pixel-integration* MTF is present whether or not it is wanted, because the ray
tracer area-integrates over each pixel; it is removed analytically in
:func:`biphase.psd.deconvolve_pixel_mtf`.

Three further effects -- beam hardening, additional detector blur, and scattered
radiation -- are the dominant real-world systematics, and are modelled here so
that their cross-sensitivities can be quantified as Type B contributions rather
than assumed negligible.  They are applied in the physically correct order:

    X  ->  polychromatic transmission  ->  detector blur  ->  scatter
       ->  Poisson + read noise  ->  linearisation

Blur precedes noise because photons spread in the scintillator before they are
counted; applying it afterwards would smooth the noise as well and understate it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Spectrum",
    "polychromatic_transmission",
    "hardening_curve",
    "linearise_calibrated",
    "apply_mtf",
    "gaussian_mtf2",
    "add_scatter",
    "expected_counts",
    "add_noise",
    "linearise",
    "simulate",
    "noise_variance_X",
    "noise_variance_field",
    "noise_psd",
    "noise_psd_field",
    "jensen_factor",
    "linearisation_snr",
    "linearisation_valid",
    "log_bias",
    "optimal_mu_X",
]


# --------------------------------------------------------------------------
# beam hardening
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Spectrum:
    """A polychromatic source: photon weights and the solid's attenuation per energy.

    ``weight`` is the incident photon fraction at each energy and is normalised on
    construction.  ``mu`` is the linear attenuation coefficient of the solid at
    that energy, in inverse length units matching the projected path length.
    """

    energy: np.ndarray
    weight: np.ndarray
    mu: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "energy", np.atleast_1d(np.asarray(self.energy, float)))
        w = np.atleast_1d(np.asarray(self.weight, dtype=float))
        object.__setattr__(self, "weight", w / w.sum())
        object.__setattr__(self, "mu", np.atleast_1d(np.asarray(self.mu, dtype=float)))

    @property
    def mu_mean(self) -> float:
        """Photon-weighted mean attenuation, i.e. the zero-thickness limit."""
        return float((self.weight * self.mu).sum())

    @property
    def hardness(self) -> float:
        """Fractional spread of ``mu`` across the spectrum; 0 is monochromatic.

        This is the single number that controls the size of every beam-hardening
        effect, which is why the sensitivity study sweeps it rather than sweeping
        tube voltage or filtration.
        """
        m = self.mu_mean
        var = float((self.weight * (self.mu - m) ** 2).sum())
        return float(np.sqrt(var) / m) if m > 0 else 0.0

    @classmethod
    def monochromatic(cls, mu: float) -> "Spectrum":
        return cls(np.array([1.0]), np.array([1.0]), np.array([float(mu)]))

    @classmethod
    def three_point(cls, mu_mean: float, hardness: float) -> "Spectrum":
        """A minimal spectrum with a prescribed mean and fractional spread.

        Three equally weighted lines at ``mu_mean (1 + h*sqrt(3/2) * {-1,0,1})``
        reproduce any desired ``hardness`` exactly while keeping the mean fixed,
        which makes the sweep a clean one-parameter family.
        """
        d = float(hardness) * np.sqrt(1.5)
        mu = mu_mean * np.array([1.0 - d, 1.0, 1.0 + d])
        return cls(np.array([0.5, 1.0, 1.5]), np.ones(3), np.maximum(mu, 1e-9))


def polychromatic_transmission(X: np.ndarray, spec: Spectrum) -> np.ndarray:
    r"""``I/I_0 = \sum_E w_E \exp(-\mu_E X)`` -- a sum of exponentials, not one."""
    X = np.asarray(X, dtype=float)
    out = np.zeros_like(X)
    for w, mu in zip(spec.weight, spec.mu):
        out += w * np.exp(-mu * X)
    return out


def hardening_curve(spec: Spectrum, X_max: float, n: int = 512):
    """Tabulate the true ``X`` vs ``-ln(I/I_0)`` relation of a step wedge.

    This is what a step-wedge calibration measures, and inverting it removes the
    *mean* beam-hardening nonlinearity exactly.  What such a calibration cannot
    remove is the effect of the nonlinearity on the fluctuations, which is the
    part that matters for a spectral measurement.
    """
    X = np.linspace(0.0, float(X_max), n)
    return X, -np.log(np.maximum(polychromatic_transmission(X, spec), 1e-300))


def linearise_calibrated(I_ratio: np.ndarray, curve) -> np.ndarray:
    """Invert a measured hardening curve instead of assuming Beer--Lambert."""
    X_tab, L_tab = curve
    L = -np.log(np.maximum(np.asarray(I_ratio, dtype=float), 1e-300))
    return np.interp(L, L_tab, X_tab)


# --------------------------------------------------------------------------
# detector blur and scatter
# --------------------------------------------------------------------------
def apply_mtf(I: np.ndarray, sigma_px: float) -> np.ndarray:
    """Gaussian detector blur of standard deviation ``sigma_px`` pixels.

    Applied to the intensity, wrapping periodically to match the periodic image.
    """
    if sigma_px <= 0:
        return np.asarray(I, dtype=float)
    from scipy.ndimage import gaussian_filter

    return gaussian_filter(np.asarray(I, dtype=float), sigma_px, mode="wrap")


def gaussian_mtf2(qx: np.ndarray, qy: np.ndarray, sigma_len: float) -> np.ndarray:
    r"""``|M(q)|^2 = \exp(-q^2\sigma^2)`` for a Gaussian blur of width ``sigma_len``."""
    q2 = np.asarray(qy, dtype=float)[:, None] ** 2 + np.asarray(qx, dtype=float)[None, :] ** 2
    return np.exp(-q2 * float(sigma_len) ** 2)


def add_scatter(I: np.ndarray, fraction: float) -> np.ndarray:
    """Add a flat scattered background at ``fraction`` of the mean primary signal.

    Scatter is smooth on the scale of the microstructure, so to first order it is
    an additive pedestal.  Its effect is to reduce the apparent attenuation, which
    biases the mean -- and therefore the thickness -- directly.
    """
    I = np.asarray(I, dtype=float)
    return I + float(fraction) * float(I.mean()) if fraction > 0 else I


def expected_counts(X: np.ndarray, mu: float, N0: float) -> np.ndarray:
    """Mean transmitted counts ``N0 exp(-mu X)`` per pixel."""
    return N0 * np.exp(-mu * np.asarray(X, dtype=float))


def add_noise(
    counts: np.ndarray, read_noise: float = 0.0, rng: np.random.Generator | None = None
) -> np.ndarray:
    """Poisson photon noise plus Gaussian read noise."""
    rng = np.random.default_rng() if rng is None else rng
    y = rng.poisson(counts).astype(float)
    if read_noise > 0:
        y = y + rng.normal(0.0, read_noise, size=y.shape)
    return y


def linearise(
    Y: np.ndarray, mu: float, N0: float, floor_frac: float = 1e-3
) -> np.ndarray:
    """Invert Beer-Lambert: ``Xhat = -ln(Y/N0)/mu``.

    ``Y`` is clipped from below at ``floor_frac * N0`` because read noise can
    drive counts to zero or negative, where the log diverges.  If the clip
    fires on a non-negligible fraction of pixels the linearisation is invalid --
    check :func:`linearisation_valid` first.
    """
    Y = np.asarray(Y, dtype=float)
    return -np.log(np.maximum(Y, floor_frac * N0) / N0) / mu


def simulate(
    X: np.ndarray,
    mu: float,
    N0: float,
    read_noise: float = 0.0,
    rng: np.random.Generator | None = None,
    spectrum: "Spectrum | None" = None,
    sigma_mtf_px: float = 0.0,
    scatter: float = 0.0,
    calibrate_hardening: bool = True,
    calib_spectrum: "Spectrum | None" = None,
) -> np.ndarray:
    """Full forward-and-back detector chain.

    ``X -> transmission -> blur -> scatter -> Poisson + read -> linearisation``.

    With all optional arguments at their defaults this is the monochromatic,
    blur-free, scatter-free chain used everywhere else in the package.  Supplying
    a ``spectrum``, ``sigma_mtf_px`` or ``scatter`` activates the corresponding
    systematic, which is how their cross-sensitivities are measured.

    ``calibrate_hardening`` inverts a step-wedge relation rather than assuming
    Beer--Lambert, which is what an experiment does.  Note what this implies:
    beam hardening is a *static, monotonic* nonlinearity in ``X``, so an exact
    calibration inverts it exactly -- pixel by pixel, not merely on average, and
    with no residual effect on the fluctuations at all.  The real sensitivity is
    therefore to calibration *error*, which is what ``calib_spectrum`` provides:
    pass a spectrum differing from the true one to invert with the wrong curve.

    ``N0 = inf`` bypasses the counting noise but still applies blur, scatter and
    hardening, so each systematic can be isolated from photon statistics.
    """
    X = np.asarray(X, dtype=float)
    mono = spectrum is None
    spec = Spectrum.monochromatic(mu) if mono else spectrum

    T = np.exp(-mu * X) if mono else polychromatic_transmission(X, spec)
    if sigma_mtf_px > 0:
        T = apply_mtf(T, sigma_mtf_px)
    if scatter > 0:
        T = add_scatter(T, scatter)

    if np.isfinite(N0):
        rng = np.random.default_rng() if rng is None else rng
        counts = rng.poisson(np.maximum(N0 * T, 0.0)).astype(float)
        if read_noise > 0:
            counts = counts + rng.normal(0.0, read_noise, size=counts.shape)
        T = np.maximum(counts, 1e-3 * N0) / N0

    if mono or not calibrate_hardening:
        return -np.log(np.maximum(T, 1e-300)) / spec.mu_mean
    cs = spec if calib_spectrum is None else calib_spectrum
    return linearise_calibrated(T, hardening_curve(cs, float(X.max()) * 1.5 + 1e-6))


def noise_variance_X(
    mu: float, N0: float, X_mean: float, read_noise: float = 0.0
) -> float:
    r"""Per-pixel variance of ``Xhat`` from the log linearisation.

    .. math::
        \sigma_X^2 = \frac{1}{\mu^2}\left[\frac{e^{\mu\langle X\rangle}}{N_0}
                     + \frac{\sigma_r^2 e^{2\mu\langle X\rangle}}{N_0^2}\right]

    obtained from ``Var(Y) = Nbar + sigma_r^2`` and ``dXhat/dY = -1/(mu Y)``.
    """
    if not np.isfinite(N0):
        return 0.0
    e = np.exp(mu * X_mean)
    return (e / N0 + read_noise**2 * e**2 / N0**2) / mu**2


def noise_psd(
    mu: float, N0: float, X_mean: float, p: float, read_noise: float = 0.0
) -> float:
    """Flat noise PSD floor ``P_n = sigma_X^2 p^2``, in length**4.

    Per-pixel white noise has continuum autocovariance ``sigma_X^2 p^2 delta^2(x)``,
    so its PSD is flat at ``sigma_X^2 p^2`` -- consistent with the estimator in
    :mod:`biphase.psd`, whose expectation for white input is ``sigma^2 p^2``.
    """
    return noise_variance_X(mu, N0, X_mean, read_noise) * p**2


def noise_variance_field(
    X: np.ndarray, mu: float, N0: float, read_noise: float = 0.0
) -> float:
    r"""Exact spatially-averaged noise variance for a *given* projected field.

    :func:`noise_variance_X` evaluates the closed form at ``<X>``, but the
    per-pixel variance grows like :math:`e^{\mu X}`, so by Jensen's inequality
    the true spatial average exceeds it:

    .. math::
        \langle\sigma_X^2\rangle = \frac{1}{\mu^2}\left[
            \frac{\langle e^{\mu X}\rangle}{N_0}
            + \frac{\sigma_r^2\langle e^{2\mu X}\rangle}{N_0^2}\right]

    For a roughly Gaussian ``X`` the two terms are inflated by
    :math:`e^{\mu^2\sigma_X^2/2}` and :math:`e^{2\mu^2\sigma_X^2}`
    respectively, which at realistic contrast is tens of per cent on the photon
    term and a factor of two on the read-noise term.

    This is the version to subtract as the noise floor in the inversion: it is
    computable from the measured radiograph itself, since ``X`` is known once
    the image has been linearised.

    The resulting noise PSD is flat to within a per cent above
    ``q ~ 0.15 q_Nyquist``, but carries a real excess at lower ``q``: the noise
    *amplitude* is modulated by the sample structure (thick regions transmit
    fewer photons), so the long-wavelength structure of ``exp(mu X)`` is
    imprinted on the noise power.  Subtracting a flat floor therefore
    over-corrects at low ``q``.  This does not affect the first structural peak,
    which sits near ``q d ~ 7``, far above the affected band.
    """
    if not np.isfinite(N0):
        return 0.0
    X = np.asarray(X, dtype=float)
    e1 = np.exp(mu * X).mean()
    e2 = np.exp(2.0 * mu * X).mean() if read_noise > 0 else 0.0
    return (e1 / N0 + read_noise**2 * e2 / N0**2) / mu**2


def noise_psd_field(
    X: np.ndarray, mu: float, N0: float, p: float, read_noise: float = 0.0
) -> float:
    """Flat noise PSD floor from :func:`noise_variance_field`, in length**4."""
    return noise_variance_field(X, mu, N0, read_noise) * p**2


def jensen_factor(
    X: np.ndarray, mu: float, N0: float, read_noise: float = 0.0
) -> float:
    """Ratio of the exact field-averaged noise variance to the closed form at <X>."""
    X = np.asarray(X, dtype=float)
    approx = noise_variance_X(mu, N0, float(X.mean()), read_noise)
    return noise_variance_field(X, mu, N0, read_noise) / approx if approx else 1.0


def log_bias(mu: float, N0: float, X_mean: float, read_noise: float = 0.0) -> float:
    """Second-order bias ``E[Xhat] - X ~ (Nbar + sigma_r^2)/(2 mu Nbar^2)``."""
    if not np.isfinite(N0):
        return 0.0
    nbar = N0 * np.exp(-mu * X_mean)
    return (nbar + read_noise**2) / (2.0 * mu * nbar**2)


def linearisation_snr(
    mu: float, N0: float, X_mean: float, read_noise: float = 0.0
) -> float:
    """Per-pixel count SNR ``Nbar / sqrt(Nbar + sigma_r^2)``."""
    if not np.isfinite(N0):
        return np.inf
    nbar = N0 * np.exp(-mu * X_mean)
    return nbar / np.sqrt(nbar + read_noise**2)


def linearisation_valid(
    mu: float, N0: float, X_mean: float, read_noise: float = 0.0, threshold: float = 5.0
) -> bool:
    """Whether the log linearisation is trustworthy.

    Below ``Nbar/sqrt(Nbar + sigma_r^2) ~ 5`` the transmitted counts approach
    zero, the log linearisation fails catastrophically, and
    :func:`noise_variance_X` underestimates the true variance by an order of
    magnitude or more.
    """
    return linearisation_snr(mu, N0, X_mean, read_noise) >= threshold


def optimal_mu_X(read_noise: float = 0.0, N0: float = np.inf) -> float:
    """Operating point maximising photon-limited contrast SNR.

    Pure photon noise gives the classic ``mu <X> = 2``.  Read noise shifts the
    optimum down; solved numerically when ``sigma_r > 0``.
    """
    if read_noise <= 0 or not np.isfinite(N0):
        return 2.0
    from scipy.optimize import minimize_scalar

    def neg_snr(t):
        # contrast SNR ~ mu / sigma_X with mu <X> = t held fixed
        return -1.0 / np.sqrt(np.exp(t) / N0 + read_noise**2 * np.exp(2 * t) / N0**2)

    return float(minimize_scalar(neg_snr, bounds=(0.1, 6.0), method="bounded").x)
