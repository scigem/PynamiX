r"""Partial structure factors from an integral-equation closure, without simulation.

The species-velocity fit of :mod:`biphase.velocity` needs the Ashcroft--Langreth
partial structure factors as weights.  Obtaining them from a discrete-element
simulation of the material is possible but heavy, so this module solves the
multicomponent Ornstein--Zernike equation under the Percus--Yevick closure
instead, which returns :math:`S_{ab}(q)` for any :math:`(\phi, c, R)` in a second
and needs nothing but those three numbers.

Method
------
Ornstein--Zernike for a mixture, in the symmetrised form
:math:`\hat H = \hat C + \hat C\hat H` with
:math:`\hat C_{ab}=\sqrt{\rho_a\rho_b}\,\tilde c_{ab}`, gives
:math:`\hat H=(I-\hat C)^{-1}\hat C` and hence

.. math:: S^{\rm AL}_{ab}(q)=\delta_{ab}+\sqrt{\rho_a\rho_b}\,\tilde h_{ab}(q)
          \quad\Longleftrightarrow\quad S^{\rm AL}=(I-\hat C)^{-1}.

The Percus--Yevick closure for additive hard spheres is
:math:`c_{ab}(r)=-[1+\gamma_{ab}(r)]` inside the contact distance
:math:`\sigma_{ab}=(d_a+d_b)/2` and zero outside, with
:math:`\gamma=h-c`.  Iterating the two together to self-consistency is a few
lines given a radial Fourier transform.

Caveat
------
Percus--Yevick is an *equilibrium* fluid theory.  A jammed granular packing is
not an equilibrium fluid, and at the densities of interest the two need not
agree.  Whether the resulting partials are accurate enough to be used as velocity
weights is an empirical question, answered in the paper rather than assumed here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["py_partials", "py_structure_factor", "PYResult"]


def _radial_ft(f: np.ndarray, r: np.ndarray, q: np.ndarray) -> np.ndarray:
    r"""3-D isotropic Fourier transform :math:`(4\pi/q)\int r f(r)\sin(qr)\,dr`."""
    from scipy.fft import dst

    dr = float(r[1] - r[0])
    return 4.0 * np.pi * dr * dst(r * f, type=1) / (2.0 * q)


def _inverse_radial_ft(F: np.ndarray, r: np.ndarray, q: np.ndarray) -> np.ndarray:
    from scipy.fft import dst

    dq = float(q[1] - q[0])
    return dq * dst(q * F, type=1) / (4.0 * np.pi**2 * r)


@dataclass(frozen=True)
class PYResult:
    """Partials plus the diagnostics needed to know whether to trust them."""

    S: dict[tuple[int, int], np.ndarray]
    converged: bool
    residual: float
    iterations: int


def _picard(rho, d, sig, r, q, inside, gamma0, mix, tol, max_iter):
    """Damped Picard iteration of OZ+PY at one density, from a given start.

    Convergence is judged on the *output* ``S(q)``, not on the change in the
    indirect correlation ``gamma(r)``.  The latter is dominated by the first
    radial grid point, where the ``1/r`` in the inverse transform amplifies a
    region that carries no physics: it reports failure at a residual of 1e-1
    while the median residual is 1e-16 and ``S(q)`` is stable to 1e-11.  Judging
    convergence on gamma made every dense solve look like a failure when the
    solver was in fact correct to a fraction of a per cent.
    """
    n = len(d)
    sq = np.sqrt(np.outer(rho, rho))
    eye = np.eye(n)
    gamma = gamma0.copy()
    prev, delta, it = np.inf, np.inf, 0
    S_prev = None
    for it in range(1, max_iter + 1):
        c = np.where(inside, -(1.0 + gamma), 0.0)
        ck = np.array([[_radial_ft(c[i, j], r, q) for j in range(n)] for i in range(n)])
        Ck = ck * sq[:, :, None]
        gk = np.empty_like(Ck)
        S_now = np.empty((n, n, len(q)))
        for m in range(len(q)):
            C = Ck[:, :, m]
            try:
                S = np.linalg.solve(eye - C, eye)
            except np.linalg.LinAlgError:
                S = eye
            S_now[:, :, m] = np.real(S)
            gk[:, :, m] = (S - eye) - C
        gnew = np.array([[_inverse_radial_ft(gk[i, j] / sq[i, j], r, q)
                          for j in range(n)] for i in range(n)])
        delta = (np.inf if S_prev is None
                 else float(np.abs(S_now - S_prev).max()))
        S_prev = S_now
        if not np.isfinite(delta) and it > 1:
            break
        if np.isfinite(delta) and delta > prev:
            mix = max(0.5 * mix, 0.02)
        prev = delta
        gamma = (1 - mix) * gamma + mix * gnew
        if delta < tol:
            break
    return gamma, delta, it


def py_partials(
    diameters: dict[int, float],
    number_density: dict[int, float],
    q_out: np.ndarray,
    n_grid: int = 2048,
    r_max: float | None = None,
    mix: float = 0.3,
    tol: float = 1e-4,
    max_iter: int = 3000,
    n_steps: int | None = None,
    return_result: bool = False,
    strict: bool = False,
):
    """Ashcroft--Langreth partials for an additive hard-sphere mixture.

    Returns ``{(a, b): S_ab(q_out)}``, in the same convention and with the same
    keys as :func:`biphase.structure.partial_structure_factors`, so the two are
    directly interchangeable.  ``return_result`` gives a :class:`PYResult`
    carrying convergence diagnostics instead.

    **Density continuation is what makes this converge.** Picard iteration of the
    OZ equation started from ``gamma = 0`` is only marginally stable above
    moderate density and simply fails at the densities of interest -- silently,
    returning a smooth and entirely wrong answer rather than an error.  The
    solution is stepped up in density instead, each step starting from the
    converged result of the previous one, which keeps every step inside the basin
    of attraction.  Always check ``converged``; ``strict`` raises instead.
    """
    types = sorted(diameters)
    n = len(types)
    d = np.array([diameters[t] for t in types], dtype=float)
    rho = np.array([number_density[t] for t in types], dtype=float)
    sig = 0.5 * (d[:, None] + d[None, :])

    r_max = float(r_max if r_max is not None else 24.0 * d.max())
    dr = r_max / n_grid
    r = (np.arange(n_grid) + 1) * dr
    dq = np.pi / (r_max + dr)
    q = (np.arange(n_grid) + 1) * dq
    inside = np.array([[r < sig[i, j] for j in range(n)] for i in range(n)])

    eta = float((np.pi / 6.0) * (rho * d**3).sum())
    if n_steps is None:
        n_steps = max(1, int(np.ceil(eta / 0.10)))

    gamma = np.zeros((n, n, n_grid))
    delta, it_tot = np.inf, 0
    for s in range(1, n_steps + 1):
        gamma, delta, it = _picard(rho * (s / n_steps), d, sig, r, q, inside,
                                   gamma, mix, tol, max_iter)
        it_tot += it

    sq = np.sqrt(np.outer(rho, rho))
    c = np.where(inside, -(1.0 + gamma), 0.0)
    ck = np.array([[_radial_ft(c[i, j], r, q) for j in range(n)] for i in range(n)])
    Ck = ck * sq[:, :, None]
    S_grid = np.empty((n, n, n_grid))
    eye = np.eye(n)
    for m in range(n_grid):
        try:
            S_grid[:, :, m] = np.linalg.solve(eye - Ck[:, :, m], eye)
        except np.linalg.LinAlgError:
            S_grid[:, :, m] = eye

    out = {}
    for i, a in enumerate(types):
        for j, b in enumerate(types):
            out[(int(a), int(b))] = np.interp(np.asarray(q_out, float), q, S_grid[i, j])

    # delta is a change in S(q) itself, so the threshold is directly
    # interpretable: 1e-3 is a thousandth of a typical peak height, far below
    # the accuracy the closure itself can claim against a real packing.
    converged = bool(np.isfinite(delta) and delta < 1e-3)
    if strict and not converged:
        raise RuntimeError(
            f"Percus-Yevick iteration did not converge (residual {delta:.2e}); "
            "the returned partials are not usable")
    if return_result:
        return PYResult(out, converged, float(delta), it_tot)
    return out


def py_structure_factor(q: np.ndarray, eta: float, sigma: float = 1.0) -> np.ndarray:
    """Single-component ``S(q)`` from the same solver, for validation."""
    rho = 6.0 * eta / (np.pi * sigma**3)
    return py_partials({1: sigma}, {1: rho}, q)[(1, 1)]
