"""Recovering S(q) from a measured radiograph spectrum, and inverting it.

Monodisperse
------------
``P_X(q) = <X> v_s F^2(q d / 2) S(q)`` (thick-sample limit), so

``S(q) = [P_obs(q) - P_n(q)] / [|M(q)|^2 <X> v_s F^2(q d / 2)]``.

The sample thickness has cancelled -- that is the reason to analyse the
spectrum rather than only its variance.

Bidisperse
----------
For a mixture the single-sphere form factor is no longer correct.  Writing the
spectral density as a self (Laue) plus distinct (interference) term,

``chi(k) = sum_a phi_a v_a F^2(k a_a) + sum_ab phi_a phi_b F F h_ab(k)``

the self term is a flat-in-S "polydispersity floor" that the monodisperse
inversion would misread as ``S(q) -> const``.  Under the decoupling
approximation (size and position statistically independent) this gives a clean
inversion, :func:`recover_S_bidisperse`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


from .psd import RadialPsd
from .units import form_factor, sphere_volume

__all__ = [
    "RecoveredS",
    "PeakFit",
    "recover_S",
    "recover_S_bidisperse",
    "fit_first_peak",
    "fit_first_peak_ml",
    "porod_fit",
    "concentration_from_porod",
    "lobe_mask",
    "fit_self_term",
    "concentration_from_self_term",
    "schulz_pdf",
    "fit_size_distribution",
    "SizeDistribution",
]


@dataclass(frozen=True)
class RecoveredS:
    q: np.ndarray
    S: np.ndarray
    sem: np.ndarray
    n: np.ndarray


@dataclass(frozen=True)
class PeakFit:
    """Fitted first-peak descriptors.  ``q1`` is in the same units as ``q``."""

    q1: float
    S1: float
    width: float
    sigma_q1: float
    n_points: int


def recover_S(
    psd: RadialPsd,
    X_mean: float,
    d: float,
    P_noise: float = 0.0,
    f2_floor: float = 1e-3,
) -> RecoveredS:
    """Invert the thick-sample forward model for a monodisperse packing.

    Bins where ``F^2`` has fallen below ``f2_floor`` are discarded rather than
    amplified: division by the form factor is ill-conditioned near its zeros
    (the first is at ``q d = 8.987``), and data there is noise multiplied by a
    large number.

    The default floor of 1e-3 keeps data out to ``q d ~ 8.6``.  It must be kept
    this low: ``F^2`` has already fallen to 0.017 at ``q d = 7.6``, just above
    the first structure peak, so a floor of 0.02 truncates the spectrum *at*
    the peak and leaves nothing to fit on its high-q side.

    The form-factor zero does impose a real limit: the first *minimum* of
    ``S(q)`` lies beyond ``q d = 8.99`` and is therefore not accessible to a
    monodisperse single-species inversion.
    """
    F2 = form_factor(psd.q * d / 2.0) ** 2
    keep = F2 > f2_floor
    denom = X_mean * sphere_volume(d) * F2[keep]
    S = (psd.P[keep] - P_noise) / denom
    return RecoveredS(q=psd.q[keep], S=S, sem=psd.sem[keep] / denom, n=psd.n[keep])


def recover_S_bidisperse(
    psd: RadialPsd,
    X_mean: float,
    H: float,
    diameters: dict[int, float],
    number_density: dict[int, float],
    P_noise: float = 0.0,
    floor: float = 1e-3,
) -> RecoveredS:
    r"""Decoupling-approximation inversion for a two-species mixture.

    With :math:`f_a(q)=v_a F(q a_a)` and :math:`\langle\cdot\rangle=\sum_a x_a(\cdot)`
    over number fractions :math:`x_a`,

    .. math::
        S(q) = 1 + \frac{P_X(q)/H - \rho\langle f^2\rangle}{\rho\langle f\rangle^2}

    The subtracted term is the self/Laue floor; omitting it is what makes the
    monodisperse formula read a flat polydispersity background as structure.
    """
    rho = sum(number_density.values())
    x = {a: n / rho for a, n in number_density.items()}
    f = {a: sphere_volume(dd) * form_factor(psd.q * dd / 2.0) for a, dd in diameters.items()}
    f1 = sum(x[a] * f[a] for a in diameters)
    f2 = sum(x[a] * f[a] ** 2 for a in diameters)

    denom = rho * f1**2
    keep = denom > floor * np.nanmax(denom)
    chi = (psd.P[keep] - P_noise) / H
    S = 1.0 + (chi - rho * f2[keep]) / denom[keep]
    return RecoveredS(q=psd.q[keep], S=S, sem=psd.sem[keep] / H / denom[keep], n=psd.n[keep])


def fit_first_peak(
    q: np.ndarray,
    S: np.ndarray,
    sem: np.ndarray | None = None,
    search: tuple[float, float] = (3.0, 12.0),
    window_frac: float = 0.08,
    model: str = "quadratic",
) -> PeakFit:
    """Locate the first maximum of ``S(q)`` by local fitting.

    Taking the largest discrete FFT bin quantises ``q1`` at the bin spacing,
    which is far coarser than the precision the inversion needs (a 0.3% error
    in ``q1 d`` already costs ~0.005 in phi).  A local quadratic through the
    points within ``window_frac`` of the peak recovers sub-bin position.

    ``search`` brackets the first peak in the *same units as* ``q``; for
    ``q`` in units of 1/d the first peak sits near 7.

    The default ``window_frac = 0.08`` was chosen by fitting synthetic peaks of
    known position: it holds the systematic bias below 0.06% at a bin spacing of
    0.2/d, against 0.3% at ``window_frac = 0.18`` where the parabola is a poor
    model for the peak shape.  Since ``d(q1 d)/d(phi) ~ 4``, a 0.3% bias in
    ``q1 d`` would already cost ~0.005 in the inferred phi.
    """
    m = (q >= search[0]) & (q <= search[1]) & np.isfinite(S)
    if m.sum() < 5:
        raise ValueError("too few points in the search window")
    qs, Ss = q[m], S[m]
    i = int(np.argmax(Ss))
    q_peak = qs[i]

    # Symmetric window.  The usable range is truncated on the right by the
    # form-factor zero (data above q d ~ 8.6 are discarded), so a nominally
    # symmetric fractional window can become one-sided near the peak -- which
    # tilts the parabola and throws the vertex badly, by up to 12% at high
    # density.  Shrink the half-width to whatever is available on both sides.
    half = window_frac * q_peak
    half = min(half, q_peak - qs[0], qs[-1] - q_peak)
    w = (qs >= q_peak - half) & (qs <= q_peak + half)
    if w.sum() < 4:
        w = np.zeros_like(qs, dtype=bool)
        lo = max(0, i - 3)
        w[lo : min(len(qs), lo + 7)] = True
    qq, ss = qs[w], Ss[w]
    weights = None
    if sem is not None:
        e = sem[m][w]
        weights = 1.0 / np.maximum(e, 1e-30) ** 2

    if model == "quadratic":
        y = ss
    elif model == "gaussian":
        y = np.log(np.maximum(ss, 1e-12))
    else:
        raise ValueError(f"unknown model {model!r}")

    # centre the abscissa: fitting a parabola in raw q at q ~ 7 is badly
    # conditioned, and the curvature term is what sets the peak position.
    t = qq - q_peak
    A = np.vstack([t**2, t, np.ones_like(t)]).T
    if weights is not None:
        Wm = np.sqrt(weights)[:, None]
        coef, *_ = np.linalg.lstsq(A * Wm, y * Wm[:, 0], rcond=None)
    else:
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b, c = coef
    if a >= 0:  # not a maximum; fall back to the discrete argmax
        return PeakFit(float(q_peak), float(Ss[i]), np.nan, np.nan, int(w.sum()))

    dq = -b / (2.0 * a)
    # The vertex of a fitted parabola can land far outside the data it was fitted
    # to when the curvature is weak (noisy spectra do this).  Such a fit carries
    # no information about the peak, so fall back to the discrete maximum.
    if not np.isfinite(dq) or abs(dq) > 0.5 * (qq.max() - qq.min()):
        return PeakFit(float(q_peak), float(Ss[i]), np.nan, np.nan, int(w.sum()))
    q1 = float(q_peak + dq)
    S1 = float(c - b**2 / (4.0 * a)) if model == "quadratic" else float(np.exp(c - b**2 / (4 * a)))
    width = float(np.sqrt(-S1 / a)) if model == "quadratic" and a < 0 else float(np.sqrt(-1.0 / (2 * a)))

    resid = y - A @ coef
    dof = max(len(y) - 3, 1)
    s2 = float(resid @ resid / dof)
    try:
        cov = s2 * np.linalg.inv(A.T @ A)
        # var(dq) for dq = -b/2a via the delta method
        var = (
            cov[1, 1] / (4 * a**2)
            + cov[0, 0] * b**2 / (4 * a**4)
            - cov[0, 1] * b / (2 * a**3)
        )
        sigma_q1 = float(np.sqrt(max(var, 0.0)))
    except np.linalg.LinAlgError:
        sigma_q1 = np.nan
    return PeakFit(q1, S1, width, sigma_q1, int(w.sum()))


def lobe_mask(q: np.ndarray, diameters, threshold: float = 0.25) -> np.ndarray:
    r"""Keep only wavevectors near the maxima of the form-factor lobes.

    Between lobes ``F^2`` passes through zero, but the *measured* spectrum cannot:
    a radial bin has finite width and ``F^2`` is convex near a zero, so the bin
    average exceeds ``F^2`` at the bin centre, and statistical scatter is
    strictly positive as well.  Fitting through those bins therefore biases any
    amplitude estimate badly and erratically -- it was the cause of a 20%
    scatter in the recovered solid fraction before this mask was applied.

    A bin is kept when ``F^2`` is at least ``threshold`` times its Porod envelope
    :math:`9/(2u^4)` for *some* species.  For a mixture the species' zeros fall
    at different wavevectors, so little data is lost.
    """
    keep = np.zeros_like(q, dtype=bool)
    for d in np.atleast_1d(diameters):
        u = np.maximum(q * d / 2.0, 1e-9)
        keep |= form_factor(u) ** 2 > threshold * 9.0 / (2.0 * u**4)
    return keep


def fit_self_term(
    psd: RadialPsd,
    H: float,
    diameters,
    q_range: tuple[float, float] = (30.0, 60.0),
    P_noise: float = 0.0,
    lobe_threshold: float = 0.25,
) -> dict[float, float]:
    r"""Recover each species' volume fraction from the high-\ ``q`` self term.

    Beyond the correlation range every interference term ``h_ab`` has died and
    Eq.~(self/distinct) leaves only

    .. math:: 	ilde\chi(q) 	o \sum_a \phi_a v_a F^2(q a_a),

    which is *exact*, not merely asymptotic.  Fitting the amplitudes
    :math:`\phi_a v_a` by linear least squares against the known basis functions
    therefore gives the volume fractions directly, with no calibration and no
    assumption about the packing.

    This is preferable to averaging ``chi q^4`` against the Porod envelope: it
    uses the oscillation structure rather than averaging it away, so the result
    does not depend on the band happening to span whole lobes.  Measured against
    penetrable spheres of known density it is accurate to about 1% and stable
    across bands, where the envelope average swings by +/-10%.

    The band must start high.  In a dense packing ``S(q)`` is still oscillating
    strongly at ``q d = 10-25`` -- that range contains the first minimum and the
    second peak -- so a fit begun there is contaminated by real structure and
    returns wildly band-dependent amplitudes.  It settles by ``q d ~ 30``, hence
    the default.  This is a genuine resolution requirement, and a stiffer one
    than the first-peak method needs: reaching ``q d = 60`` below
    ``0.6 q_Nyquist`` requires roughly 32 pixels per small-particle diameter,
    against about 8 for the peak position.

    A residual few per cent bias remains in the *absolute* volume fractions,
    from the tail of the interference term.  It affects both species almost
    equally, so the concentration ratio is considerably more accurate than
    either amplitude alone.

    Returns ``{diameter: phi_a}``.
    """
    diameters = list(np.atleast_1d(diameters))
    m = (psd.q >= q_range[0]) & (psd.q <= q_range[1])
    m &= lobe_mask(psd.q, diameters, lobe_threshold)
    if m.sum() < 2 * len(diameters):
        raise ValueError("too few usable bins for the self-term fit")

    q = psd.q[m]
    chi = (psd.P[m] - P_noise) / H
    A = np.stack([sphere_volume(d) * form_factor(q * d / 2.0) ** 2 for d in diameters], 1)
    w = np.sqrt(psd.n[m])[:, None]
    coef, *_ = np.linalg.lstsq(A * w, chi * w[:, 0], rcond=None)
    return {d: float(c) for d, c in zip(diameters, coef)}


def concentration_from_self_term(phi_by_d: dict[float, float], d_large: float) -> float:
    """Large-species volume concentration from :func:`fit_self_term` output."""
    tot = sum(phi_by_d.values())
    return float(phi_by_d.get(d_large, 0.0) / tot) if tot > 0 else float("nan")


def porod_fit(
    psd: RadialPsd,
    H: float,
    q_range: tuple[float, float],
    P_noise: float = 0.0,
    free_floor: bool = False,
) -> float:
    r"""Fit the Porod tail to recover the specific surface area.

    At large ``k`` all interference terms vanish and
    :math:`\tilde\chi(k)\to 2\pi s/k^4` exactly, whatever the size
    distribution.  So ``s`` -- and hence the Sauter diameter -- is measurable
    with no model and no calibration.  The form factor oscillates about this
    envelope, so the fit is taken over a band spanning several lobes.

    No zero-avoiding mask is applied, and that is deliberate.  A mask of the
    kind :func:`recover_S` needs -- where dividing by a vanishing :math:`F^2`
    would amplify noise without limit -- is actively harmful here, because this
    estimator *averages* rather than divides.  Dropping the bins near the zeros
    removes the troughs of the oscillation while keeping its peaks, so the mean
    is biased upwards.  For a bidisperse packing the effect is large: masking on
    the small-species zeros inflates ``s`` by about 12 per cent and so reports
    ``d_32`` about 12 per cent low, against 1--2 per cent unmasked.

    The ``q^4`` weight makes this estimator unusually sensitive to the noise
    floor, because it is precisely the top of the band -- where the signal is
    weakest -- that is weighted hardest.  With ``free_floor`` the fit becomes

    .. math:: \tilde\chi(q)\,q^4 = 2\pi s + b\,q^4 ,

    so any residual error in ``P_noise`` is absorbed by ``b`` instead of being
    read as surface area.  At :math:`N_0=10^4` photons per pixel this turns a
    14 per cent bias in ``d_32`` into 2 per cent; on noise-free data it costs
    about one percentage point, the price of the extra free parameter.  At
    :math:`N_0=10^3` neither form works -- the top of the band is then noise.
    """
    m = (psd.q >= q_range[0]) & (psd.q <= q_range[1])
    q = psd.q[m]
    chi = (psd.P[m] - P_noise) / H
    if not free_floor:
        return float(np.mean(chi * q**4) / (2.0 * np.pi))
    A = np.stack([np.ones(q.size), q**4], 1)
    w = np.sqrt(psd.n[m])[:, None]
    coef, *_ = np.linalg.lstsq(A * w, chi * q**4 * w[:, 0], rcond=None)
    return float(coef[0] / (2.0 * np.pi))


def concentration_from_porod(d32: float, d_s: float, d_l: float) -> float:
    r"""Large-species volume concentration from the Sauter diameter.

    Since :math:`1/d_{32} = (1-c)/d_s + c/d_l`,

    .. math:: c = \frac{d_s/d_{32} - 1}{d_s/d_l - 1}

    Calibration-free: it needs only the Sauter diameter from the Porod tail and
    the two grain sizes, both readable from the same spectrum.
    """
    return float((d_s / d32 - 1.0) / (d_s / d_l - 1.0))


# --------------------------------------------------------------------------
# continuous size distributions
# --------------------------------------------------------------------------
def schulz_pdf(d: np.ndarray, d_mean: float, cv: float) -> np.ndarray:
    r"""Schulz (gamma) number density of diameters, mean ``d_mean``, CV ``cv``.

    Chosen over a log-normal because its moments are closed form, so the target
    :math:`d_{32}=d_{\rm mean}(1+cv^2)(1+2cv^2)/(1+cv^2)` and every weighting the
    inversion needs follow analytically rather than by quadrature.
    """
    from scipy.stats import gamma

    k = 1.0 / max(cv, 1e-9) ** 2
    return gamma.pdf(np.asarray(d, dtype=float), a=k, scale=d_mean / k)


@dataclass(frozen=True)
class SizeDistribution:
    """Recovered continuous size distribution."""

    phi: float
    d_mean: float
    cv: float
    residual: float

    @property
    def d32(self) -> float:
        """Sauter mean diameter of the fitted Schulz distribution."""
        c2 = self.cv**2
        return float(self.d_mean * (1.0 + c2) * (1.0 + 2.0 * c2) / (1.0 + c2))

    def moments(self) -> dict[str, float]:
        c2 = self.cv**2
        return {
            "d10": self.d_mean,
            "d32": self.d32,
            "sigma": self.d_mean * self.cv,
        }


def _chi_model(q, d_mean, cv, n_nodes=160):
    r"""``<v^2 F^2>`` per particle for a Schulz distribution, by quadrature."""
    lo = max(d_mean * (1.0 - 6.0 * cv), 1e-3 * d_mean)
    hi = d_mean * (1.0 + 8.0 * cv)
    d = np.linspace(lo, hi, n_nodes)
    w = schulz_pdf(d, d_mean, cv)
    w = w / np.trapezoid(w, d)
    v = sphere_volume(d)
    # integrate over the distribution at every wavevector at once
    F2 = form_factor(np.outer(q, d) / 2.0) ** 2
    return np.trapezoid(w * v**2 * F2, d, axis=1), float(np.trapezoid(w * v, d))


def fit_size_distribution(
    psd: RadialPsd,
    H: float,
    q_range: tuple[float, float] = (30.0, 60.0),
    P_noise: float = 0.0,
    d0: float = 1.0,
    cv0: float | None = None,
    multistart: bool = True,
) -> SizeDistribution:
    r"""Recover a continuous size distribution from the high-\ ``q`` self term.

    The bidisperse self-term fit generalises directly: beyond the correlation
    range the interference terms vanish and

    .. math:: \tilde\chi(q)\to\rho\,\langle v^2F^2(qa)\rangle
              = \phi\,\frac{\langle v^2F^2\rangle}{\langle v\rangle},

    where the averages are over the number distribution of diameters.  Fitting
    the two shape parameters of a Schulz distribution plus one amplitude
    therefore gives the solid fraction and the full distribution -- in
    particular ``d_32``, which is the quantity the Porod law also delivers, here
    obtained with the oscillation structure used rather than averaged away.

    ``phi`` enters linearly and is eliminated by variable projection, leaving a
    two-parameter search over ``(d_mean, cv)``.

    The residual surface has more than one basin, and the spurious one can sit
    *below* the value at the truth.  A fit seeded near the answer will therefore
    look excellent while being unidentifiable, so ``multistart`` (the default)
    runs a grid of starts and keeps the global best.  Pass ``cv0`` with
    ``multistart=False`` only to reproduce a single local solve; never seed from
    a quantity derived from the truth.

    ``d_32`` is much better determined than ``(d_mean, cv)`` separately -- see
    :func:`profile_d32`, which is the honest way to put an interval on it.
    """
    from scipy.optimize import least_squares

    m = (psd.q >= q_range[0]) & (psd.q <= q_range[1])
    m &= lobe_mask(psd.q, [d0])
    if m.sum() < 8:
        raise ValueError("too few usable bins for a distribution fit")
    q = psd.q[m]
    chi = (psd.P[m] - P_noise) / H
    w = np.sqrt(psd.n[m])

    def resid(x):
        d_mean, cv = abs(x[0]), abs(x[1])
        num, vbar = _chi_model(q, d_mean, max(cv, 1e-3))
        basis = num / vbar                      # chi / phi
        amp = float(np.sum(w * chi * basis) / max(np.sum(w * basis**2), 1e-30))
        return w * (chi - amp * basis)

    bounds = ([0.2 * d0, 1e-3], [5 * d0, 1.0])
    if multistart:
        starts = [(dm, cv) for dm in np.linspace(0.5 * d0, 2.2 * d0, 7)
                  for cv in (0.02, 0.10, 0.20, 0.35, 0.55, 0.80)]
    else:
        starts = [(d0, 0.2 if cv0 is None else cv0)]
    out = None
    for s in starts:
        try:
            o = least_squares(resid, list(s), bounds=bounds, xtol=1e-12, ftol=1e-12)
        except Exception:
            continue
        if out is None or o.cost < out.cost:
            out = o
    if out is None:
        raise RuntimeError("size-distribution fit did not converge from any start")
    d_mean, cv = abs(out.x[0]), abs(out.x[1])
    num, vbar = _chi_model(q, d_mean, cv)
    basis = num / vbar
    phi = float(np.sum(w * chi * basis) / max(np.sum(w * basis**2), 1e-30))
    r = w * (chi - phi * basis)
    return SizeDistribution(phi=phi, d_mean=d_mean, cv=cv,
                            residual=float(np.linalg.norm(r) / max(np.linalg.norm(w * chi), 1e-30)))


@dataclass(frozen=True)
class D32Profile:
    """Profile-likelihood interval for the Sauter mean diameter.

    ``d32`` is the maximiser, ``lo``/``hi`` the endpoints where the profile
    log-likelihood has fallen by ``dell`` from its peak (``dell=0.5`` is the
    one-parameter, one-sigma convention), and ``grid``/``ell`` the curve itself.
    ``identifiable`` is False when the interval runs to the edge of the scan,
    which is what an unidentifiable fit looks like.
    """

    d32: float
    lo: float
    hi: float
    phi: float
    cv: float
    grid: np.ndarray
    ell: np.ndarray
    identifiable: bool

    @property
    def half_width(self) -> float:
        return 0.5 * (self.hi - self.lo)


def profile_d32(
    psd: RadialPsd,
    H: float,
    q_range: tuple[float, float] = (30.0, 60.0),
    P_noise: float = 0.0,
    d0: float = 1.0,
    span: tuple[float, float] = (0.6, 2.4),
    n_grid: int = 120,
    dell: float = 0.5,
) -> D32Profile:
    r"""Profile the Gamma log-likelihood over the Sauter mean diameter.

    The two shape parameters of the Schulz distribution are only weakly
    determined by the high-``q`` self term, but the combination
    :math:`d_{32}=d_{\rm mean}(1+2\,{\rm cv}^2)` -- the surface-weighted mean,
    which is the diameter that permeability, drag and interfacial-area
    arguments actually want -- is far better constrained.  Reparameterising in
    :math:`(d_{32},{\rm cv})` and profiling out ``cv`` and the amplitude gives
    an interval for it directly, without ever having to trust the individual
    shape parameters.

    The statistics are the ones stated for every other estimator in this work:
    a radially averaged bin is Gamma distributed with shape ``n_j``, so

    .. math:: \ell=-\sum_j n_j\left[S_j/\mu_j+\ln\mu_j\right],

    which is *not* a sum of squares.
    """
    from scipy.optimize import minimize_scalar

    m = (psd.q >= q_range[0]) & (psd.q <= q_range[1])
    m &= lobe_mask(psd.q, [d0])
    if m.sum() < 8:
        raise ValueError("too few usable bins for a distribution fit")
    q = psd.q[m]
    chi = (psd.P[m] - P_noise) / H
    n = psd.n[m]
    good = chi > 0
    q, chi, n = q[good], chi[good], n[good]
    if q.size < 8:
        raise ValueError("too few positive bins for a distribution fit")

    def ell_at(d32, cv):
        cv = float(np.clip(cv, 1e-3, 0.95))
        num, vbar = _chi_model(q, d32 / (1.0 + 2.0 * cv**2), cv)
        basis = num / vbar
        # amplitude that maximises the Gamma likelihood is the n-weighted mean
        # of the ratio, in closed form
        amp = float(np.sum(n * chi / basis) / max(np.sum(n), 1e-30))
        mu = amp * basis
        return float(-np.sum(n * (chi / mu + np.log(mu)))), amp, cv

    grid = np.linspace(span[0] * d0, span[1] * d0, n_grid)
    prof = np.empty_like(grid)
    cvs = np.empty_like(grid)
    for i, g in enumerate(grid):
        r = minimize_scalar(lambda c: -ell_at(g, c)[0], bounds=(1e-3, 0.95),
                            method="bounded", options={"xatol": 1e-6})
        prof[i] = -r.fun
        cvs[i] = float(np.clip(r.x, 1e-3, 0.95))
    k = int(np.argmax(prof))
    peak = prof[k]
    above = prof >= peak - dell
    idx = np.flatnonzero(above)
    lo, hi = float(grid[idx[0]]), float(grid[idx[-1]])
    identifiable = bool(idx[0] > 0 and idx[-1] < n_grid - 1)
    _, amp, cv = ell_at(grid[k], cvs[k])
    return D32Profile(d32=float(grid[k]), lo=lo, hi=hi, phi=amp, cv=cv,
                      grid=grid, ell=prof, identifiable=identifiable)


def fit_first_peak_ml(
    q: np.ndarray,
    S: np.ndarray,
    n: np.ndarray,
    search: tuple[float, float] = (4.0, 11.0),
    band: tuple[float, float] = (0.72, 1.28),
) -> PeakFit:
    r"""Maximum-likelihood first-peak fit under the correct spectral statistics.

    A radially averaged periodogram value is the mean of ``n_j`` independent
    exponential samples, i.e. Gamma distributed with shape ``n_j`` and mean
    ``mu_j``.  Ordinary least squares assumes Gaussian errors of *equal*
    variance, which is wrong twice over: the errors are Gamma, and their scale is
    proportional to the value itself.  Dropping constants, the log-likelihood is

    .. math:: \ell = -\sum_j n_j\left[S_j/\mu_j + \ln\mu_j\right],

    maximised here over a Gaussian-on-a-sloping-baseline model of the peak.  This
    is the estimator whose variance the Cramer-Rao bound in
    :mod:`biphase.uncertainty` describes, so it is the one that can attain it;
    the local parabola reaches only about a third of that bound.

    ``band`` is the fit window as a fraction of the located peak position.
    """
    from scipy.optimize import minimize

    q = np.asarray(q, dtype=float)
    S = np.asarray(S, dtype=float)
    n = np.maximum(np.asarray(n, dtype=float), 1.0)

    m = (q >= search[0]) & (q <= search[1]) & np.isfinite(S)
    if m.sum() < 6:
        raise ValueError("too few points in the search window")
    qs, Ss, ns = q[m], S[m], n[m]
    q_peak = qs[int(np.argmax(Ss))]

    w = (qs >= band[0] * q_peak) & (qs <= band[1] * q_peak)
    if w.sum() < 6:
        w = np.ones_like(qs, dtype=bool)
    qq, ss, nn = qs[w], Ss[w], ns[w]

    base = float(np.median(ss))
    amp0 = max(float(ss.max()) - base, 1e-6)

    def model(theta):
        q1, lw, la, b0, b1 = theta
        width = np.exp(lw)
        amp = np.exp(la)
        return (b0 + b1 * (qq - q1)
                + amp * np.exp(-0.5 * ((qq - q1) / width) ** 2))

    def neg_ll(theta):
        mu = model(theta)
        if not np.all(np.isfinite(mu)) or np.any(mu <= 1e-9):
            return 1e12
        return float(np.sum(nn * (ss / mu + np.log(mu))))

    x0 = np.array([q_peak, np.log(0.6), np.log(amp0), max(base, 1e-3), 0.0])
    res = minimize(neg_ll, x0, method="Nelder-Mead",
                   options=dict(xatol=1e-8, fatol=1e-10, maxiter=20000, maxfev=20000))
    q1, lw, la, b0, b1 = res.x
    width = float(np.exp(lw))
    if not (qq.min() < q1 < qq.max()) or not res.success:
        return fit_first_peak(q, S, search=search)

    # Marginal, not conditional, uncertainty on the peak position: invert the
    # full Hessian rather than taking 1/sqrt of its (0,0) element.  The latter
    # holds amplitude, width and baseline fixed and understates the uncertainty
    # by more than a factor of two, because those parameters are strongly
    # correlated with the position.
    k = len(res.x)
    steps = np.array([1e-4 * max(abs(q1), 1.0), 1e-3, 1e-3, 1e-3 * max(abs(b0), 1.0), 1e-3])
    Hn = np.zeros((k, k))
    f0 = neg_ll(res.x)
    for i in range(k):
        ei = np.zeros(k); ei[i] = steps[i]
        for j in range(i, k):
            ej = np.zeros(k); ej[j] = steps[j]
            if i == j:
                Hn[i, i] = (neg_ll(res.x + ei) - 2 * f0 + neg_ll(res.x - ei)) / steps[i] ** 2
            else:
                Hn[i, j] = Hn[j, i] = (
                    neg_ll(res.x + ei + ej) - neg_ll(res.x + ei - ej)
                    - neg_ll(res.x - ei + ej) + neg_ll(res.x - ei - ej)
                ) / (4 * steps[i] * steps[j])
    try:
        cov = np.linalg.inv(Hn)
        sigma = float(np.sqrt(cov[0, 0])) if cov[0, 0] > 0 else np.nan
    except np.linalg.LinAlgError:
        sigma = np.nan
    S1 = float(b0 + np.exp(la))
    return PeakFit(float(q1), S1, width, sigma, int(w.sum()))
