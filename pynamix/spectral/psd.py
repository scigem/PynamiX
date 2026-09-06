"""Absolute-amplitude power spectral density estimation.

The whole method divides a measured spectrum by ``<X> * v_s * F^2(qd/2)`` to get
a dimensionless structure factor, so the *absolute* amplitude of the estimator
matters -- not just its shape.  This module fixes that normalisation once.

Units
-----
``X`` has units of length, so ``P_X`` has units of **length**4**.  Check against
main.tex Eq. (25): ``<X>[L] * v_s[L^3] * F^2[1] * S[1] = L^4``, and
``Var(X) = int P_X d^2q/(2pi)^2`` gives ``L^4 * L^-2 = L^2``.

Estimator
---------
For an ``N x N`` region of interest on pixel pitch ``p`` with separable taper
``W``::

    Phat[k,l] = p^2 * |FFT2(W * dX)[k,l]|^2 / sum(W^2)

The divisor is the window **power** ``sum(W**2)``, not ``(sum(W))**2`` -- the
latter preserves the amplitude of a line spectrum and is wrong for a continuous
PSD.  Derivation: the continuous transform is ``p^2 * fft2``, the PSD is
``|.|^2 / A`` with area ``A = (N p)^2``, giving ``p^4 |fft|^2 / (N^2 p^2)``;
the window then replaces ``N^2`` by ``sum(W^2)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Psd2D",
    "RadialPsd",
    "taper",
    "periodogram",
    "welch_psd2d",
    "deconvolve_pixel_mtf",
    "deconvolve_gaussian_mtf",
    "pixel_mtf2",
    "radial_average",
    "parseval_residual",
]


@dataclass(frozen=True)
class Psd2D:
    """Two-dimensional PSD on the FFT wavevector grid.

    Attributes
    ----------
    qx, qy : 1-D angular wavevector axes (rad / length), from ``fftfreq``.
    P : (ny, nx) PSD in length**4, in unshifted FFT order.
    p : pixel pitch (length).
    n_avg : number of independent periodograms averaged into ``P``.
    """

    qx: np.ndarray
    qy: np.ndarray
    P: np.ndarray
    p: float
    n_avg: int = 1

    @property
    def q_nyquist(self) -> float:
        return np.pi / self.p

    def qmag(self) -> np.ndarray:
        """(ny, nx) array of |q|."""
        return np.hypot(self.qy[:, None], self.qx[None, :])


@dataclass(frozen=True)
class RadialPsd:
    """Radially averaged PSD.

    Attributes
    ----------
    q : bin abscissa, the *mean* |q| within each annulus (not the bin centre).
    P : arithmetic mean of the 2D PSD over the annulus, length**4.
    n : number of 2D bins in the annulus.
    sem : standard error of ``P``; each 2D bin is an independent chi^2_2 sample
        so the relative error is ``1/sqrt(n * n_avg)``.
    """

    q: np.ndarray
    P: np.ndarray
    n: np.ndarray
    sem: np.ndarray


def taper(n: int, kind: str = "tukey", alpha: float = 0.25) -> np.ndarray:
    """Separable 1-D window of length ``n``.

    ``hann`` uses the *periodic* (DFT-even) form, ``np.hanning(n+1)[:n]``.
    ``tukey`` with ``alpha=0.25`` is the default: it suppresses leakage from the
    ROI edge while keeping the equivalent noise bandwidth close to unity, so
    peak *heights* are barely biased.
    """
    if kind == "boxcar":
        return np.ones(n)
    if kind == "hann":
        return np.hanning(n + 1)[:n]
    if kind == "tukey":
        if alpha <= 0:
            return np.ones(n)
        if alpha >= 1:
            return np.hanning(n + 1)[:n]
        w = np.ones(n)
        m = int(np.floor(alpha * (n - 1) / 2.0)) + 1
        k = np.arange(m)
        ramp = 0.5 * (1.0 + np.cos(np.pi * (2.0 * k / (alpha * (n - 1)) - 1.0)))
        w[:m] = ramp
        w[n - m :] = ramp[::-1]
        return w
    raise ValueError(f"unknown window {kind!r}")


def _window2d(shape: tuple[int, int], kind: str, alpha: float) -> np.ndarray:
    ny, nx = shape
    return np.outer(taper(ny, kind, alpha), taper(nx, kind, alpha))


def periodogram(x: np.ndarray, p: float, window: str = "tukey", alpha: float = 0.25) -> Psd2D:
    """Single-ROI absolute-amplitude periodogram of the field ``x``.

    The *window-weighted* mean is removed, so the DC bin is exactly zero.
    Removing the plain mean instead leaves a non-zero DC bin whose leakage
    biases the lowest few annuli.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 2:
        raise ValueError("x must be 2-D")
    W = _window2d(x.shape, window, alpha)
    sw = W.sum()
    mean_w = float((W * x).sum() / sw)
    f = np.fft.fft2(W * (x - mean_w))
    P = (p**2) * np.abs(f) ** 2 / (W**2).sum()
    ny, nx = x.shape
    return Psd2D(
        qx=2.0 * np.pi * np.fft.fftfreq(nx, d=p),
        qy=2.0 * np.pi * np.fft.fftfreq(ny, d=p),
        P=P,
        p=p,
    )


def welch_psd2d(
    X: np.ndarray,
    p: float,
    roi: int = 512,
    overlap: float = 0.5,
    window: str = "tukey",
    alpha: float = 0.25,
) -> Psd2D:
    """Welch-average periodograms over tiled ROIs of size ``roi``.

    ``overlap`` is the fractional overlap between neighbouring ROIs.  If the
    image is smaller than ``roi`` the whole image is used as a single ROI.
    """
    X = np.asarray(X, dtype=float)
    ny, nx = X.shape
    roi = min(roi, ny, nx)
    step = max(1, int(round(roi * (1.0 - overlap))))
    ys = list(range(0, ny - roi + 1, step)) or [0]
    xs = list(range(0, nx - roi + 1, step)) or [0]

    acc = None
    n = 0
    for y0 in ys:
        for x0 in xs:
            pg = periodogram(X[y0 : y0 + roi, x0 : x0 + roi], p, window, alpha)
            acc = pg.P if acc is None else acc + pg.P
            n += 1
    assert acc is not None
    return Psd2D(qx=pg.qx, qy=pg.qy, P=acc / n, p=p, n_avg=n)


def _sinc(u: np.ndarray) -> np.ndarray:
    """``sin(u)/u`` with the removable singularity handled (NOT numpy's sinc)."""
    return np.sinc(np.asarray(u) / np.pi)


def pixel_mtf2(qx: np.ndarray, qy: np.ndarray, p: float) -> np.ndarray:
    """``|M_pix|^2 = [sinc(qx p/2) sinc(qy p/2)]^2`` on the 2-D grid.

    This is the transfer function of *area integration over a square pixel*,
    which the ray tracer reproduces exactly by supersampling.  It is intrinsic
    to any pixelated detector and must be removed before the form-factor
    division, even when no additional detector blur is modelled.
    """
    mx = _sinc(qx[None, :] * p / 2.0)
    my = _sinc(qy[:, None] * p / 2.0)
    return (mx * my) ** 2


def deconvolve_pixel_mtf(psd: Psd2D) -> Psd2D:
    """Divide out the pixel MTF, per 2-D bin.

    Must be done *before* radial averaging: ``|M_pix|^2`` is anisotropic (it
    differs by 10-20% between the axis and the diagonal near Nyquist), so
    averaging first and dividing after injects a systematic, q-dependent error.
    """
    m2 = pixel_mtf2(psd.qx, psd.qy, psd.p)
    return Psd2D(psd.qx, psd.qy, psd.P / m2, psd.p, psd.n_avg)


def deconvolve_gaussian_mtf(psd: Psd2D, sigma_len: float) -> Psd2D:
    r"""Divide out an additional Gaussian detector blur, ``|M|^2=e^{-q^2\sigma^2}``.

    Applied on top of :func:`deconvolve_pixel_mtf`, which handles the pixel
    aperture.  The correction grows exponentially with ``q``, so an error in
    ``sigma`` is amplified at exactly the wavevectors the measurement uses; the
    sensitivity to mis-estimating it is quantified in the paper.
    """
    if sigma_len <= 0:
        return psd
    q2 = psd.qy[:, None] ** 2 + psd.qx[None, :] ** 2
    return Psd2D(psd.qx, psd.qy, psd.P * np.exp(q2 * float(sigma_len) ** 2), psd.p, psd.n_avg)


def radial_average(
    psd: Psd2D,
    drop_low_bins: int = 3,
    q_max_frac_nyquist: float = 0.6,
) -> RadialPsd:
    """Radially average a 2-D PSD.

    Uses the *arithmetic mean* over each annulus -- each 2-D bin is an
    independent chi^2_2 sample, so uniform weights are optimal.  The abscissa is
    the mean ``|q|`` within the annulus rather than the nominal bin centre,
    which removes a real bias at low bin index where annuli are sparse.

    Never *sums* over the annulus.  An annular sum carries a ``2 pi q`` Jacobian
    and produces a spectrum proportional to ``q P(q)``; that is the convention
    used by Dulanjalee et al. (2020) and is handled separately in
    :mod:`biphase.dulanjalee`.
    """
    qmag = psd.qmag()
    dq = 2.0 * np.pi / (psd.P.shape[0] * psd.p)
    idx = np.rint(qmag / dq).astype(np.int64)

    nbin = idx.max() + 1
    counts = np.bincount(idx.ravel(), minlength=nbin)
    qsum = np.bincount(idx.ravel(), weights=qmag.ravel(), minlength=nbin)
    psum = np.bincount(idx.ravel(), weights=psd.P.ravel(), minlength=nbin)

    good = counts > 0
    good[:drop_low_bins] = False  # mean removal high-passes the lowest annuli
    q = qsum[good] / counts[good]
    P = psum[good] / counts[good]
    n = counts[good]

    keep = q <= q_max_frac_nyquist * psd.q_nyquist
    q, P, n = q[keep], P[keep], n[keep]
    sem = P / np.sqrt(n * psd.n_avg)
    return RadialPsd(q=q, P=P, n=n, sem=sem)


def parseval_residual(x: np.ndarray, p: float, window: str = "tukey", alpha: float = 0.25) -> float:
    """Relative Parseval residual for :func:`periodogram`; must be ~1e-15.

    ``sum(Phat) / (N p)^2 == sum((W dX)^2) / sum(W^2)``.  A pure normalisation
    check with no physics in it, so it isolates FFT-convention errors.
    """
    x = np.asarray(x, dtype=float)
    ny, nx = x.shape
    W = _window2d(x.shape, window, alpha)
    mean_w = (W * x).sum() / W.sum()
    lhs = periodogram(x, p, window, alpha).P.sum() / (ny * p) / (nx * p)
    rhs = (W * (x - mean_w)) ** 2
    rhs = rhs.sum() / (W**2).sum()
    return abs(lhs - rhs) / abs(rhs)
