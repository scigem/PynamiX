r"""Species-resolved velocimetry from a pair of radiographs.

A rigid translation by :math:`\bm u` between two frames multiplies every Fourier
component by a pure phase, so the cross-spectrum

.. math:: C(\bm q) = \tilde X_2(\bm q)\,\tilde X_1^*(\bm q) = |\tilde X|^2 e^{-\ii\bm q\cdot\bm u}

carries the displacement in its phase and nothing else.  That is standard phase
correlation.  The useful observation for a *mixture* is that the projected field
is a sum over species, and the projection is linear, so

.. math::
    C(\bm q) = \sum_a w_a(q)\,e^{-\ii\bm q\cdot\bm u_a},
    \qquad w_a(q) = H\,\phi_a v_a F^2(q a_a)\,S_{aa}(q),

with the *same* weights that appear in the spectral density.  Because the two
species have different form factors and different structure-factor peaks, the
weights have genuinely different ``q`` dependence, and a fit over the whole
spectrum recovers both displacements at once.

This separates species without tracking a single particle, and without the two
spectral rings needing to be individually resolvable: the peak-fraction
heuristic reads two peak *heights* and fails once the rings merge below
``R ~ 2``, whereas this fit uses the whole spectral shape.

Phase wrapping sets the working range.  The model is written in terms of
``exp(-i q.u)``, so wrapping is represented exactly rather than needing to be
unwrapped, but the optimiser can find a wrong branch when ``q_max |u| > pi``.
:func:`phase_correlation_shift` supplies a starting point good to a fraction of
a pixel, and :func:`max_unambiguous_shift` reports the limit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .units import form_factor, sphere_volume

__all__ = [
    "cross_spectrum",
    "welch_cross_spectrum",
    "phase_correlation_shift",
    "max_unambiguous_shift",
    "suggest_band",
    "species_weights",
    "fit_displacement",
    "fit_species_displacements",
    "SpeciesVelocity",
]


# --------------------------------------------------------------------------
# synthesising a moving packing
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class CrossSpectrum:
    qx: np.ndarray
    qy: np.ndarray
    C: np.ndarray  # complex, unshifted FFT order
    p: float

    def qmag(self) -> np.ndarray:
        return np.hypot(self.qy[:, None], self.qx[None, :])

    @property
    def q_nyquist(self) -> float:
        return np.pi / self.p


def cross_spectrum(
    X1: np.ndarray, X2: np.ndarray, p: float, window: str = "boxcar", alpha: float = 0.25
) -> CrossSpectrum:
    """``C(q) = FFT(X2) conj(FFT(X1))`` with the window-weighted means removed."""
    from .psd import _window2d

    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)
    if X1.shape != X2.shape:
        raise ValueError("frames must have the same shape")
    W = _window2d(X1.shape, window, alpha)
    sw = W.sum()
    f1 = np.fft.fft2(W * (X1 - (W * X1).sum() / sw))
    f2 = np.fft.fft2(W * (X2 - (W * X2).sum() / sw))
    ny, nx = X1.shape
    return CrossSpectrum(
        qx=2.0 * np.pi * np.fft.fftfreq(nx, d=p),
        qy=2.0 * np.pi * np.fft.fftfreq(ny, d=p),
        C=f2 * np.conj(f1) * (p**2) / (W**2).sum(),
        p=p,
    )


def welch_cross_spectrum(
    X1: np.ndarray,
    X2: np.ndarray,
    p: float,
    roi: int = 256,
    overlap: float = 0.5,
    window: str = "boxcar",
    alpha: float = 0.25,
) -> CrossSpectrum:
    """Cross-spectrum averaged over tiled regions of interest.

    This averaging is not optional, it is what makes the species fit possible.
    The single-region cross-spectrum of one realisation is

    ``C = |Xs|^2 e^{-i ts} + |Xl|^2 e^{-i tl} + Xs Xl^* e^{-i ts} + Xl Xs^* e^{-i tl}``,

    in which the per-bin powers are chi^2 speckle rather than the smooth
    ensemble weights, and the two cross terms carry random phase.  A
    deterministic model of the *expectation* cannot fit that bin by bin -- it
    leaves a residual of order unity however good the displacements are.

    Averaging over regions estimates the expectation, in which the cross terms
    vanish and the powers tend to the smooth weights, while leaving the phase
    ``q.u`` untouched: it depends only on ``q``, which is common to every region.
    The cost is a coarser wavevector grid, ``2 pi / (roi * p)``.
    """
    from .psd import _window2d

    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)
    if X1.shape != X2.shape:
        raise ValueError("frames must have the same shape")
    ny, nx = X1.shape
    roi = min(roi, ny, nx)
    step = max(1, int(round(roi * (1.0 - overlap))))
    ys = list(range(0, ny - roi + 1, step)) or [0]
    xs = list(range(0, nx - roi + 1, step)) or [0]

    W = _window2d((roi, roi), window, alpha)
    sw = W.sum()
    # Normalised exactly as biphase.psd.periodogram, so that C and the PSD share
    # units and their ratio is a proper coherence.  Without this the two differ
    # by sum(W^2)/p^2 -- a factor of order 1e8 at these sizes -- which is silent
    # in any fit that carries a free amplitude, and wrong in anything that does
    # not.
    norm = (p**2) / (W**2).sum()
    acc = np.zeros((roi, roi), dtype=complex)
    n = 0
    for y0 in ys:
        for x0 in xs:
            a = X1[y0 : y0 + roi, x0 : x0 + roi]
            b = X2[y0 : y0 + roi, x0 : x0 + roi]
            fa = np.fft.fft2(W * (a - (W * a).sum() / sw))
            fb = np.fft.fft2(W * (b - (W * b).sum() / sw))
            acc += fb * np.conj(fa) * norm
            n += 1
    return CrossSpectrum(
        qx=2.0 * np.pi * np.fft.fftfreq(roi, d=p),
        qy=2.0 * np.pi * np.fft.fftfreq(roi, d=p),
        C=acc / n,
        p=p,
    )


def max_unambiguous_shift(q_max: float) -> float:
    """Largest ``|u|`` for which ``q.u < pi`` everywhere in the fit band.

    The limit is set by the top of the fitted band alone; the pixel pitch enters
    only through which ``q_max`` is reachable in the first place.
    """
    return float(np.pi / q_max)


def phase_correlation_shift(X1: np.ndarray, X2: np.ndarray, p: float) -> np.ndarray:
    """Coarse displacement from the phase-correlation peak, in length units.

    Returns ``u = (ux, uy)`` such that ``X2(x) ~ X1(x - u)``, where ``ux`` is
    along **axis 1** of the array (columns) and ``uy`` along axis 0 (rows).
    Every displacement in this module uses that order, so an image whose first
    array axis is the physical ``x`` will report its motion in ``uy``.

    Sub-pixel position is refined by a parabola through the peak and its two
    neighbours in each direction, which is enough to initialise the phase fit.
    """
    cs = cross_spectrum(X1, X2, p)
    R = cs.C / np.maximum(np.abs(cs.C), 1e-30)
    corr = np.real(np.fft.ifft2(R))
    ny, nx = corr.shape
    k = np.unravel_index(np.argmax(corr), corr.shape)

    def _sub(axis_vals):
        a, b, c = axis_vals
        den = a - 2 * b + c
        return 0.0 if den == 0 else 0.5 * (a - c) / den

    dy = _sub([corr[(k[0] - 1) % ny, k[1]], corr[k], corr[(k[0] + 1) % ny, k[1]]])
    dx = _sub([corr[k[0], (k[1] - 1) % nx], corr[k], corr[k[0], (k[1] + 1) % nx]])
    # a peak at index n > N/2 is a negative shift
    sy = k[0] + dy - (ny if k[0] > ny // 2 else 0)
    sx = k[1] + dx - (nx if k[1] > nx // 2 else 0)
    return np.array([sx * p, sy * p])


def suggest_band(diameters, p: float, n_pixels: int, peak_qd: float = 7.5) -> dict:
    r"""Fit band and region size matched to the *largest* species.

    Each species carries its velocity information mainly around its own
    structure-factor peak, at :math:`q_a \simeq Q_1/d_a`.  For a large size
    ratio the big grains' peak sits at low ``q``, so a band chosen for the small
    grains can exclude it altogether -- at ``R = 6`` the large peak falls at
    ``q d_s = 1.25``, right at the edge of a band starting at 1.0, and the fit
    then has almost no leverage on the large species.  Two things follow: the
    band must start below ``Q_1/d_l``, and the region of interest must be large
    enough that the wavevector grid ``2 pi / (roi p)`` actually resolves it.

    Returns ``{"q_min", "q_max", "roi"}``.
    """
    d = np.atleast_1d([float(x) for x in (diameters.values() if isinstance(diameters, dict) else diameters)])
    d_max, d_min = float(d.max()), float(d.min())
    if d_max == d_min:
        # single species: nothing to separate, so use the whole reliable band
        q_min = 0.5
    else:
        q_min = 0.5 * peak_qd / d_max
    q_max = min(11.0 / d_min, 0.6 * np.pi / p)
    # need several grid steps below the large-species peak: 2 pi / (roi p) < q_min / 2
    roi_needed = 4.0 * np.pi / (max(q_min, 1.0) * p)
    roi = int(2 ** np.ceil(np.log2(max(128.0, roi_needed))))
    return {"q_min": float(q_min), "q_max": float(q_max), "roi": int(min(roi, n_pixels))}


def species_weights(
    q: np.ndarray,
    diameters: dict[int, float],
    number_density: dict[int, float],
    H: float,
    S_ab: dict[tuple[int, int], np.ndarray] | None = None,
) -> dict[int, np.ndarray]:
    r"""Per-species weight in the cross-spectrum, *including* the cross partial.

    Write the frame-1 transform as a sum over species, :math:`\tilde X=A+B`.  A
    differential translation gives
    :math:`\tilde X_2 = Ae^{-\ii\theta_s}+Be^{-\ii\theta_l}`, so

    .. math::
        C = (|A|^2 + AB^*)e^{-\ii\theta_s} + (|B|^2 + BA^*)e^{-\ii\theta_l},

    and in expectation the weight attached to species ``a`` is *not* its own
    spectral density but

    .. math::
        W_a(q) = \sum_b \sqrt{\rho_a\rho_b}\,v_av_b\,F(qa_a)F(qa_b)\,S^{\rm AL}_{ab}(q),

    i.e. its own term plus **half-open share of the cross term**: every ``(a,b)``
    term of the spectral density, Eq.~(bipartial), attaches to species ``a``
    through the factor ``A`` it contains.  By construction
    :math:`\sum_a W_a = \tilde\chi`, the total spectral density.

    With ``S_ab`` omitted the partials default to :math:`\delta_{ab}`, which is
    the uncorrelated-species limit and is exact for a Poisson mixture.

    Including the cross partial is not a refinement.  ``S_sl`` is substantially
    non-zero in a dense mixture and is often negative, so ``W_a`` can be
    negative; a model built from positive weights alone forces the fitted phase
    to lie between the two species' phases, whereas the measured effective
    displacement demonstrably falls outside that interval at low ``q``.  Omitting
    it does not degrade the fit gracefully -- it makes it wrong.
    """
    types = sorted(diameters)
    F = {a: form_factor(q * diameters[a] / 2.0) for a in types}
    v = {a: sphere_volume(diameters[a]) for a in types}
    out = {}
    for a in types:
        w = np.zeros_like(np.asarray(q, dtype=float))
        for b in types:
            # Ashcroft-Langreth partials tend to delta_ab at large k, so the
            # no-information default is the Kronecker delta -- NOT unity for the
            # cross term, which would assert perfect inter-species correlation.
            S = np.ones_like(w) if a == b else np.zeros_like(w)
            if S_ab is not None and (a, b) in S_ab:
                S = S_ab[(a, b)]
            w = w + (np.sqrt(number_density[a] * number_density[b]) * v[a] * v[b] * F[a] * F[b] * S)
        out[a] = H * w
    return out


def _band(cs: CrossSpectrum, q_min: float, q_max: float):
    qm = cs.qmag()
    m = (qm >= q_min) & (qm <= q_max)
    qx = np.broadcast_to(cs.qx[None, :], qm.shape)[m]
    qy = np.broadcast_to(cs.qy[:, None], qm.shape)[m]
    return qx, qy, qm[m], cs.C[m]


def fit_displacement(
    cs: CrossSpectrum, q_min: float = 0.0, q_max: float | None = None, u0: np.ndarray | None = None
) -> np.ndarray:
    """Single-species displacement from the cross-spectrum phase ramp.

    Maximises ``Re sum_q C(q) exp(+i q.u)`` over ``u``: the maximum-likelihood
    estimator for a common phase ramp, weighting each bin by its own ``|C|`` so
    that high-power bins near the structure peak dominate.
    """
    from scipy.optimize import minimize

    q_max = 0.6 * cs.q_nyquist if q_max is None else q_max
    qx, qy, _, C = _band(cs, max(q_min, 1e-12), q_max)
    if u0 is None:
        u0 = np.zeros(2)

    def neg(u):
        return -float(np.real(C * np.exp(1j * (qx * u[0] + qy * u[1]))).sum())

    res = minimize(
        neg, np.asarray(u0, dtype=float), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=4000)
    )
    return res.x


@dataclass(frozen=True)
class SpeciesVelocity:
    """Result of a two-species displacement fit."""

    displacement: dict[int, np.ndarray]
    amplitude: dict[int, complex]
    residual: float
    success: bool

    def velocity(self, dt: float) -> dict[int, np.ndarray]:
        return {t: u / dt for t, u in self.displacement.items()}


def fit_species_displacements(
    cs: CrossSpectrum,
    diameters: dict[int, float],
    number_density: dict[int, float],
    H: float,
    S_ab: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] | None = None,
    q_min: float = 0.5,
    q_max: float | None = None,
    u0: dict[int, np.ndarray] | None = None,
    bound: float | None = None,
) -> "SpeciesVelocity":
    r"""Fit one displacement per species to the mixture cross-spectrum.

    Model: :math:`C(\bm q)=\sum_a A_a W_a(q)\,e^{-\ii\bm q\cdot\bm u_a}` with the
    weights of :func:`species_weights`.

    The amplitudes ``A_a`` enter linearly and the displacements non-linearly, so
    they are separated by variable projection: for each trial set of
    displacements the amplitudes follow in closed form from complex linear least
    squares, and only the displacements are searched over.  That turns a
    six-parameter fit into a four-parameter one and removes the strong
    amplitude--displacement coupling that otherwise stalls the optimiser.

    ``S_ab`` supplies the measured partial structure factors as
    ``{(a, b): (k, S)}``, interpolated onto the fit band.  They matter: see
    :func:`species_weights`.

    ``bound`` limits ``|u|`` per component; it defaults to the phase-wrapping
    limit ``pi / q_max``, beyond which the model has aliased solutions that the
    optimiser will happily fall into.
    """
    from scipy.optimize import least_squares

    q_max = 0.6 * cs.q_nyquist if q_max is None else q_max
    qx, qy, qm, C = _band(cs, max(q_min, 1e-12), q_max)
    types = sorted(diameters)

    S_on_band = None
    if S_ab is not None:
        S_on_band = {ab: np.interp(qm, kS[0], kS[1]) for ab, kS in S_ab.items()}
    Wd = species_weights(qm, diameters, number_density, H, S_on_band)
    W = np.stack([Wd[t] for t in types], axis=1)

    if u0 is None:
        u0 = {t: np.zeros(2) for t in types}
    x0 = np.concatenate([np.asarray(u0[t], dtype=float) for t in types])
    lim = float(np.pi / q_max) if bound is None else float(bound)
    x0 = np.clip(x0, -0.999 * lim, 0.999 * lim)

    def design(x):
        M = np.empty((len(qm), len(types)), dtype=complex)
        for i in range(len(types)):
            M[:, i] = W[:, i] * np.exp(-1j * (qx * x[2 * i] + qy * x[2 * i + 1]))
        return M

    def resid(x):
        M = design(x)
        A, *_ = np.linalg.lstsq(M, C, rcond=None)
        r = C - M @ A
        return np.concatenate([r.real, r.imag])

    out = least_squares(resid, x0, bounds=(-lim, lim), xtol=1e-13, ftol=1e-13, max_nfev=20000)
    M = design(out.x)
    A, *_ = np.linalg.lstsq(M, C, rcond=None)
    r = C - M @ A
    return SpeciesVelocity(
        displacement={t: out.x[2 * i : 2 * i + 2] for i, t in enumerate(types)},
        amplitude={t: complex(A[i]) for i, t in enumerate(types)},
        residual=float(np.linalg.norm(r) / max(np.linalg.norm(C), 1e-30)),
        success=bool(out.success),
    )
