"""Calibration of the first-peak position against solid fraction, and inversion.

The method needs a second measurement independent of thickness.  The mean
radiograph gives only ``<X> = phi H``; the first structure-factor peak supplies
``phi`` through a calibration ``q1 d = f1(phi)``, after which ``H = <X>/phi``.

The calibration is *not* a universal equation of state.  It is specific to a
preparation protocol, and the package treats it that way: a ``Calibration``
carries the protocol label it was fitted for.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["Calibration", "fit_calibration", "infer_H"]


@dataclass
class Calibration:
    r"""Calibration ``q1 d = f1(phi)`` as a low-order least-squares fit.

    An interpolant through the per-density means -- the obvious choice -- is a
    poor calibration.  It passes exactly through noisy points, so measurement
    noise propagates directly into the curve *shape*: the local slope wanders,
    and the curve carries no honest uncertainty of its own (the standard error
    of a mean of two or three realisations is both unstable and far too small).
    Coverage testing exposes this immediately.

    A degree-``deg`` polynomial fitted to the individual measurements instead
    gives a stable slope and, through its parameter covariance, a proper
    prediction band.  Inverting a new measurement then propagates three terms:
    the new measurement's own noise, the structural spread of an individual
    packing about the curve, and the uncertainty of the curve itself.
    """

    phi: np.ndarray
    Q1: np.ndarray
    Q1_err: np.ndarray
    protocol: str = "unknown"
    deg: int = 1
    coef: np.ndarray = field(init=False, repr=False)
    cov: np.ndarray = field(init=False, repr=False)
    u_scatter: float = field(init=False, default=0.0)

    def __post_init__(self):
        order = np.argsort(self.phi)
        self.phi = np.asarray(self.phi, dtype=float)[order]
        self.Q1 = np.asarray(self.Q1, dtype=float)[order]
        self.Q1_err = np.asarray(self.Q1_err, dtype=float)[order]

        X = np.vander(self.phi, self.deg + 1, increasing=True)
        w = 1.0 / np.maximum(self.Q1_err, 1e-12) ** 2 if np.all(self.Q1_err > 0) \
            else np.ones_like(self.phi)
        W = np.diag(w)
        XtWX = X.T @ W @ X
        self.coef = np.linalg.solve(XtWX, X.T @ W @ self.Q1)
        resid = self.Q1 - X @ self.coef
        dof = max(len(self.phi) - (self.deg + 1), 1)
        # scale the covariance by the observed residual spread, so the band
        # reflects how well the model actually describes the data
        s2 = float(resid @ (w * resid)) / dof / max(w.mean(), 1e-30)
        self.cov = s2 * np.linalg.inv(XtWX)
        self.u_scatter = float(np.sqrt(float(resid @ resid) / dof))

    def _x(self, phi):
        return np.vander(np.atleast_1d(np.asarray(phi, dtype=float)),
                         self.deg + 1, increasing=True)

    def __call__(self, phi):
        out = self._x(phi) @ self.coef
        return out if np.ndim(phi) else float(out[0])

    def sensitivity(self, phi):
        """``d(q1 d)/d(phi)`` -- this alone controls inversion precision."""
        p = np.atleast_1d(np.asarray(phi, dtype=float))
        d = sum(k * self.coef[k] * p ** (k - 1) for k in range(1, self.deg + 1))
        return d if np.ndim(phi) else float(np.atleast_1d(d)[0])

    def curve_uncertainty(self, phi) -> float:
        """Standard uncertainty of the fitted curve at ``phi`` (not of a point)."""
        x = self._x(phi)
        return float(np.sqrt(np.maximum((x @ self.cov * x).sum(axis=1), 0.0))[0])

    @property
    def is_monotonic(self) -> bool:
        g = np.linspace(self.phi.min(), self.phi.max(), 200)
        s = np.atleast_1d(self.sensitivity(g))
        return bool(np.all(s > 0) or np.all(s < 0))

    def invert(self, Q1_obs: float, Q1_err: float = 0.0,
               u_structural: float = 0.0) -> tuple[float, float]:
        r"""Inverse map with full propagation.

        .. math::
            u(\phi)^2 = \frac{u(q_1)^2_{\rm meas} + u^2_{\rm struct}
                              + u^2_{\rm curve}(\phi)}{(\mathrm{d}f_1/\mathrm{d}\phi)^2}
        """
        if not self.is_monotonic:
            raise ValueError("calibration is not monotonic; cannot invert uniquely")
        grid = np.linspace(self.phi.min(), self.phi.max(), 8001)
        vals = np.atleast_1d(self(grid))
        if not (vals.min() <= Q1_obs <= vals.max()):
            return float("nan"), float("nan")
        phi_hat = float(np.interp(Q1_obs, vals, grid) if vals[0] < vals[-1]
                        else np.interp(Q1_obs, vals[::-1], grid[::-1]))
        slope = float(self.sensitivity(phi_hat))
        u_tot = np.sqrt(Q1_err**2 + u_structural**2 + self.curve_uncertainty(phi_hat) ** 2)
        return phi_hat, (abs(u_tot / slope) if slope != 0 else float("inf"))


def fit_calibration(phi, Q1, Q1_err=None, protocol: str = "unknown",
                    deg: int = 1) -> Calibration:
    """Fit ``f1(phi)`` to the individual measurements (not to per-density means).

    Fitting the individual points keeps the degrees of freedom honest and lets
    the residual spread estimate the structural scatter directly, which is then
    available as ``u_scatter``.
    """
    phi = np.asarray(phi, dtype=float)
    Q1 = np.asarray(Q1, dtype=float)
    Q1_err = (np.zeros_like(Q1) if Q1_err is None
              else np.asarray(Q1_err, dtype=float))
    return Calibration(phi, Q1, Q1_err, protocol, deg)


def infer_H(
    X_mean: float, phi: float, X_mean_err: float = 0.0, phi_err: float = 0.0
) -> tuple[float, float]:
    r"""``H = <X>/phi`` with propagated relative errors.

    .. math::
        \left(\frac{\sigma_H}{H}\right)^2 \simeq
        \left(\frac{\sigma_{\bar X}}{\bar X}\right)^2
        + \left(\frac{\sigma_\phi}{\phi}\right)^2
    """
    H = X_mean / phi
    rel = np.hypot(X_mean_err / X_mean if X_mean else 0.0,
                   phi_err / phi if phi else 0.0)
    return float(H), float(H * rel)
