"""Synthetic radiographs for the pynamix.spectral tests.

Deliberately self-contained: fifteen lines of numpy rather than a dependency on
a ray tracer, so the analysis half of the package can be tested without the
forward model that generated the data for the paper.

Penetrable ("Poisson") spheres are the case worth having.  Uniformly random
centres give ``S(q) = 1`` exactly, so the projected spectrum is
``rho H v_s^2 F^2(qa)`` with no free parameters at all -- which makes it the one
case that can check the transform normalisation, the pixel transfer function and
the absolute amplitude simultaneously.
"""

import numpy as np

__all__ = ["project", "poisson_spheres", "bidisperse_spheres"]


def project(centres, radii, L, n_pixels, supersample=4):
    """Projected solid path length, periodic in the two transverse directions.

    ``centres[:, 0]`` is x and lands on axis 1 (columns), matching the
    convention that :mod:`pynamix.spectral.velocity` uses for displacements.
    """
    p = L / n_pixels
    pf = p / supersample
    fine = np.zeros((n_pixels * supersample,) * 2)
    for c, a in zip(centres, radii):
        n = int(np.ceil(a / pf)) + 1
        jc, ic = int(round(c[0] / pf)), int(round(c[1] / pf))
        di = np.arange(-n, n + 1)
        xj = (jc + di) * pf + 0.5 * pf - c[0]
        yi = (ic + di) * pf + 0.5 * pf - c[1]
        b2 = yi[:, None] ** 2 + xj[None, :] ** 2
        fine[np.ix_(np.mod(ic + di, n_pixels * supersample), np.mod(jc + di, n_pixels * supersample))] += 2 * np.sqrt(
            np.maximum(a * a - b2, 0.0)
        )
    return fine.reshape(n_pixels, supersample, n_pixels, supersample).mean(axis=(1, 3)), p


def poisson_spheres(phi=0.2, L=16.0, d=1.0, seed=0):
    """Centres of penetrable spheres at solid fraction ``phi``. S(q) = 1 exactly."""
    rng = np.random.default_rng(seed)
    n = int(phi * L**3 / (np.pi * d**3 / 6))
    return rng.random((n, 3)) * L, np.full(n, d / 2), n


def bidisperse_spheres(L=16.0, d_s=1.0, d_l=2.0, n_s=3000, n_l=400, seed=0):
    """Two penetrable species, for the self-term and velocity fits."""
    rng = np.random.default_rng(seed)
    centres = rng.random((n_s + n_l, 3)) * L
    radii = np.concatenate([np.full(n_s, d_s / 2), np.full(n_l, d_l / 2)])
    return centres, radii
