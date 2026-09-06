"""Dimensionless conventions and the analytic sphere functions.

Internal convention throughout ``biphase``: the *small* particle diameter is the
unit of length, ``d_s = 1``.  Every length (box edge, slab thickness, pixel
pitch, wavevector) is expressed in that unit, so a wavevector carries units of
1/d_s and the projected path length ``X`` carries units of d_s.

The one place this matters is the absolute amplitude of the power spectrum:
``P_X`` has units of length**4, so numbers quoted in this package are in
d_s**4.  See :mod:`biphase.psd`.
"""

from __future__ import annotations

import numpy as np

#: First zero of the sphere amplitude form factor F(u), i.e. tan(u) = u.
U_FIRST_ZERO = 4.493409457909064

#: ...expressed as q*d, since F is evaluated at u = q*d/2.
QD_FIRST_ZERO = 2.0 * U_FIRST_ZERO  # 8.98682


def sphere_volume(d: np.ndarray | float) -> np.ndarray | float:
    """Volume of a sphere of diameter ``d``."""
    return np.pi * np.asarray(d, dtype=float) ** 3 / 6.0


def form_factor(u: np.ndarray | float) -> np.ndarray:
    r"""Sphere amplitude form factor :math:`F(u)=3(\sin u-u\cos u)/u^3`.

    ``F(0) = 1``.  Evaluate at ``u = q*a = q*d/2``.  The small-``u`` branch uses
    the series :math:`1 - u^2/10 + u^4/280` to avoid catastrophic cancellation:
    the naive expression loses all precision below ``u ~ 1e-2`` because
    ``sin u - u cos u`` is a difference of two nearly equal numbers of size
    ``u`` whose leading terms cancel to order ``u**3``.
    """
    u = np.asarray(u, dtype=float)
    out = np.empty_like(u)
    small = np.abs(u) < 1e-2
    us = u[small]
    out[small] = 1.0 - us**2 / 10.0 + us**4 / 280.0
    ub = u[~small]
    out[~small] = 3.0 * (np.sin(ub) - ub * np.cos(ub)) / ub**3
    return out


def form_factor_qd(qd: np.ndarray | float) -> np.ndarray:
    """``form_factor`` evaluated from ``q*d`` rather than ``q*a``."""
    return form_factor(np.asarray(qd, dtype=float) / 2.0)


def specific_surface(phi_by_diameter: dict[float, float]) -> float:
    r"""Specific surface area :math:`s = 6\sum_a \phi_a/d_a`.

    Interfacial area per unit *total* volume, for a mixture with volume
    fraction ``phi_a`` of spheres of diameter ``d_a``.  This is the quantity in
    the Porod tail :math:`\tilde\chi(k)\to 2\pi s/k^4`.
    """
    return 6.0 * sum(phi / d for d, phi in phi_by_diameter.items())


def sauter_diameter(phi_by_diameter: dict[float, float]) -> float:
    r"""Sauter (surface-weighted) mean diameter :math:`d_{32}=6\phi/s`."""
    phi = sum(phi_by_diameter.values())
    return 6.0 * phi / specific_surface(phi_by_diameter)


def porod_chi(k: np.ndarray, s: float) -> np.ndarray:
    r"""Asymptotic spectral density :math:`\tilde\chi(k)\to 2\pi s/k^4`."""
    return 2.0 * np.pi * s / np.asarray(k, dtype=float) ** 4
