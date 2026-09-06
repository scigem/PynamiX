r"""Per-measurement uncertainty, Cramer-Rao bounds, and the uncertainty budget.

An ensemble scatter is not an uncertainty.  A field instrument has to report
``u(y)`` for *the one radiograph that was taken*, and that estimate has to be
validated -- a stated uncertainty that has not been coverage-tested is an
assertion.  This module provides three things:

* **Type A**, estimated per measurement by resampling the spectral modes that
  went into the fit (:func:`bootstrap_first_peak`), with the analytic
  Cramer-Rao floor for comparison (:func:`crb_first_peak`,
  :func:`crb_displacement`);
* **Type B**, propagated through explicit sensitivity coefficients
  (:class:`Budget`);
* a **coverage test** (:func:`coverage`), which asks whether the intervals
  actually contain the truth as often as they claim.

Statistical basis
-----------------
A periodogram bin is exponentially distributed about the true spectral density
(the modulus squared of a complex Gaussian mode), independently bin to bin.  So
the log-likelihood of a spectral model ``P(q; theta)`` is

.. math:: \ell = -\sum_j\left[\hat I_j/P_j + \ln P_j\right],

whose Fisher information is :math:`I_{\theta\theta}=\sum_j n_j(\partial_\theta\ln P_j)^2`
with ``n_j`` the number of independent modes contributing to bin ``j``.  Every
bound below follows from that expression; nothing here assumes Gaussian errors
on the spectrum itself, which would be wrong.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "bootstrap_first_peak",
    "crb_first_peak",
    "crb_displacement",
    "efficiency",
    "Contribution",
    "Budget",
    "coverage",
    "modes_in_annulus",
]


# --------------------------------------------------------------------------
# Type A: resampling and analytic bounds
# --------------------------------------------------------------------------
def modes_in_annulus(q: float, L: float) -> float:
    """Number of independent Fourier modes in the annulus at ``q``.

    The wavevector grid spacing is ``2 pi / L``, so an annulus of that width at
    radius ``q`` holds ``2 pi q / (2 pi / L)^2 * (2 pi / L) = q L`` bins.  The
    count therefore scales with the *linear size of the field of view in
    particle diameters*, and not at all with the pixel count -- which is why
    precision improves as ``(L/d)^{1/2}`` and why buying more pixels does not
    help.
    """
    return float(q * L)


def bootstrap_first_peak(
    q: np.ndarray,
    S: np.ndarray,
    n: np.ndarray,
    n_boot: int = 400,
    search: tuple[float, float] = (4.0, 11.0),
    window_frac: float = 0.20,
    seed: int = 0,
) -> tuple[float, float]:
    """Per-measurement ``(q1, u(q1))`` by parametric bootstrap of the spectrum.

    Each radial value is the mean of ``n_j`` independent exponential samples, so
    it is Gamma-distributed with shape ``n_j`` and mean ``S_j``.  Resampling from
    that distribution and re-fitting propagates exactly the statistics the
    estimator actually sees, without assuming the peak-fit residual is a good
    error estimate (it is not -- it conflates model error with noise).
    """
    from .invert import fit_first_peak

    rng = np.random.default_rng(seed)
    base = fit_first_peak(q, S, search=search, window_frac=window_frac)
    draws = np.empty(n_boot)
    nn = np.maximum(np.asarray(n, dtype=float), 1.0)
    for b in range(n_boot):
        Sb = rng.gamma(shape=nn, scale=np.asarray(S) / nn)
        try:
            draws[b] = fit_first_peak(q, Sb, search=search,
                                      window_frac=window_frac).q1
        except Exception:
            draws[b] = np.nan
    good = np.isfinite(draws)
    return float(base.q1), float(np.nanstd(draws[good])) if good.sum() > 10 else np.nan


def crb_first_peak(
    q: np.ndarray, S: np.ndarray, n: np.ndarray,
    band: tuple[float, float] = (5.0, 9.0),
) -> float:
    r"""Cramer-Rao lower bound on the first-peak position.

    .. math::
        u_{\min}(q_1)=\Big[\textstyle\sum_j n_j\,(\partial\ln S_j/\partial q_1)^2\Big]^{-1/2}

    The derivative with respect to peak *position* is evaluated by translating
    the measured peak shape, i.e. ``d ln S / d q_1 = -d ln S / d q``, which needs
    no parametric model of the peak.
    """
    q = np.asarray(q, dtype=float)
    S = np.maximum(np.asarray(S, dtype=float), 1e-12)
    n = np.asarray(n, dtype=float)
    m = (q >= band[0]) & (q <= band[1])
    if m.sum() < 5:
        return np.nan
    dlnS = np.gradient(np.log(S[m]), q[m])
    info = float(np.sum(n[m] * dlnS**2))
    return float(1.0 / np.sqrt(info)) if info > 0 else np.nan


def crb_displacement(
    qx: np.ndarray, qy: np.ndarray, coherence: np.ndarray, n_avg: np.ndarray | float = 1.0
) -> float:
    r"""Cramer-Rao lower bound on a displacement measured from cross-spectrum phase.

    For a coherent signal of magnitude-squared coherence :math:`\gamma^2` averaged
    over ``n`` segments, the phase variance is
    :math:`\sigma_\theta^2=(1-\gamma^2)/(2n\gamma^2)`, and the displacement enters
    as :math:`\theta=\bm q\cdot\bm u`, so

    .. math:: u_{\min}=\Big[\textstyle\sum_j |\bm q_j|^2/\sigma_{\theta,j}^2\Big]^{-1/2},

    taking the isotropic average over direction.  The bound falls as ``1/q``, so
    it is the *high*-wavevector modes that carry displacement information -- the
    opposite of the peak-position measurement, which is dominated by the modes
    near ``q_1``.
    """
    g2 = np.clip(np.asarray(coherence, dtype=float) ** 2, 1e-9, 1 - 1e-9)
    n = np.asarray(n_avg, dtype=float)
    var_theta = (1 - g2) / (2 * n * g2)
    q2 = np.asarray(qx, dtype=float) ** 2 + np.asarray(qy, dtype=float) ** 2
    info = float(np.sum(q2 / (2.0 * var_theta)))  # /2 for the isotropic direction split
    return float(1.0 / np.sqrt(info)) if info > 0 else np.nan


def efficiency(u_achieved: float, u_bound: float) -> float:
    """Estimator efficiency ``(u_min/u)^2``; 1 means the bound is attained."""
    if not np.isfinite(u_achieved) or u_achieved <= 0:
        return np.nan
    return float((u_bound / u_achieved) ** 2)


# --------------------------------------------------------------------------
# Budget
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Contribution:
    """One row of an uncertainty budget."""

    name: str
    value: float          # x_i
    u: float              # u(x_i)
    sensitivity: float    # dy/dx_i
    kind: str = "B"       # "A" or "B"
    note: str = ""

    @property
    def contribution(self) -> float:
        """|dy/dx_i| u(x_i), in units of the measurand."""
        return abs(self.sensitivity) * self.u


@dataclass
class Budget:
    """A GUM-style uncertainty budget for one measurand.

    ``correlations`` holds any non-zero correlation coefficients between named
    contributions.  They are not a refinement here: ``H = <X>/phi`` takes both
    inputs from the same radiograph, so treating them as independent -- the
    obvious thing to do -- misstates the combined uncertainty.
    """

    measurand: str
    value: float
    contributions: list[Contribution] = field(default_factory=list)
    correlations: dict[tuple[str, str], float] = field(default_factory=dict)

    def add(self, *args, **kw) -> "Budget":
        self.contributions.append(Contribution(*args, **kw))
        return self

    def _rho(self, a: str, b: str) -> float:
        return self.correlations.get((a, b), self.correlations.get((b, a), 0.0))

    @property
    def combined(self) -> float:
        """Combined standard uncertainty, including correlation terms."""
        c = self.contributions
        var = sum(x.contribution**2 for x in c)
        for i, xi in enumerate(c):
            for xj in c[i + 1:]:
                r = self._rho(xi.name, xj.name)
                if r:
                    var += 2 * r * xi.contribution * xj.contribution
        return float(np.sqrt(max(var, 0.0)))

    def expanded(self, k: float = 2.0) -> float:
        return k * self.combined

    def by_kind(self, kind: str) -> float:
        return float(np.sqrt(sum(x.contribution**2
                                 for x in self.contributions if x.kind == kind)))

    def table(self) -> "object":
        import pandas as pd

        rows = [dict(quantity=x.name, value=x.value, u=x.u, kind=x.kind,
                     sensitivity=x.sensitivity, contribution=x.contribution,
                     percent=100 * x.contribution**2 / max(self.combined**2, 1e-30),
                     note=x.note) for x in self.contributions]
        return pd.DataFrame(rows).sort_values("contribution", ascending=False)


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
def coverage(estimate, truth, sigma, k: float = 1.0) -> dict:
    """Fraction of intervals ``estimate +/- k sigma`` that contain ``truth``.

    The decisive check on a claimed uncertainty.  For ``k = 1`` a correctly
    stated Gaussian uncertainty gives 0.683; materially less means the
    uncertainty is understated, materially more that it is inflated.  The
    binomial standard error on the observed fraction is returned so that
    "materially" can be judged rather than eyeballed.
    """
    e = np.asarray(estimate, dtype=float)
    t = np.asarray(truth, dtype=float)
    s = np.asarray(sigma, dtype=float)
    m = np.isfinite(e) & np.isfinite(t) & np.isfinite(s) & (s > 0)
    if m.sum() == 0:
        return dict(n=0, covered=np.nan, expected=np.nan, se=np.nan)
    inside = np.abs(e[m] - t[m]) <= k * s[m]
    p = float(inside.mean())
    n = int(m.sum())
    from scipy.stats import norm
    expected = float(2 * norm.cdf(k) - 1)
    return dict(n=n, covered=p, expected=expected,
                se=float(np.sqrt(max(p * (1 - p), 1e-12) / n)),
                z=float((p - expected) / max(np.sqrt(expected * (1 - expected) / n), 1e-12)))
