"""Coarse-to-fine P-DVC.

Gauss-Newton converges only when the starting displacement is within about
one feature size of the solution. Each level here blurs the reference volume
by a Gaussian of standard deviation ``blur`` voxels and the measured
projections by the same blur mapped to the detector (``blur`` times the
magnification, divided by the pixel pitch), then bins the detector. Blurring
both by the same 3D kernel keeps them consistent for a parallel beam and to
first order for a cone beam; the box average of binning is not modelled, so
coarse levels are kept for initialisation only. Each level's mesh may be
coarser; its result is carried to the next level by trilinear prolongation.
"""

from dataclasses import dataclass

import numpy as np
import torch
from scipy import ndimage

from .geometry import bin_projections
from .mesh import RegularHexMesh
from .model import DynamicProjectionModel
from .solver import GaussNewton


@dataclass
class Level:
    bin: int = 1
    blur: float = 0.0  # voxels
    elem: int | tuple = 4  # elements per axis
    max_iter: int = 20
    tol: float = 1e-3


def default_levels(elem=4):
    """Three levels: (bin 4, blur 2, half mesh), (bin 2, blur 1), (bin 1, no blur)."""
    coarse = max(1, elem // 2)
    return [
        Level(bin=4, blur=2.0, elem=coarse, tol=1e-2),
        Level(bin=2, blur=1.0, elem=elem, tol=5e-3),
        Level(bin=1, blur=0.0, elem=elem),
    ]


def blur_projections(p, geometry, blur):
    if blur <= 0:
        return p
    sig = blur * geometry.magnification() / geometry.pixel_pitch()
    out = np.empty_like(p)
    for i in range(p.shape[0]):
        out[i] = ndimage.gaussian_filter(p[i], sig[i], mode="nearest")
    return out


def multiscale_solve(
    f0,
    geometry,
    measured,
    time_basis,
    times,
    levels=None,
    reference=None,
    samples_per_voxel=2.0,
    model_kw=None,
    solver_kw=None,
    verbose=True,
):
    """Run Gauss-Newton level by level. Returns (results, models), one per level.

    ``reference``: optional projections of the undeformed specimen at the same
    poses, blurred and binned like ``measured`` and passed to the solver.
    """
    levels = levels or default_levels()
    model_kw = dict(model_kw or {})
    solver_kw = dict(solver_kw or {})

    def as_np(x):
        return np.asarray(x.cpu() if torch.is_tensor(x) else x)

    p = as_np(measured)
    ref = as_np(reference) if reference is not None else None
    N = float(np.linalg.norm(f0.shape))
    a, prev_mesh = None, None
    results, models = [], []
    for k, lev in enumerate(levels):
        f0_l = ndimage.gaussian_filter(f0, lev.blur) if lev.blur > 0 else f0
        g_l = geometry.binned(lev.bin)
        p_l = bin_projections(blur_projections(p, geometry, lev.blur), lev.bin)
        ref_l = None
        if ref is not None:
            ref_l = bin_projections(blur_projections(ref, geometry, lev.blur), lev.bin)
        mesh = RegularHexMesh(f0.shape, lev.elem)
        n_samples = int(np.ceil(samples_per_voxel * N / max(1.0, lev.blur)))
        model = DynamicProjectionModel(f0_l, g_l, mesh, time_basis, times, n_samples=n_samples, **model_kw)
        a0 = None
        if a is not None:
            a0 = prev_mesh.prolongate(a, mesh) if prev_mesh.node_shape != mesh.node_shape else a
        if verbose:
            print(
                f"-- level {k}: bin {lev.bin}, blur {lev.blur} vox, mesh {mesh.n_elem}, "
                f"detector {g_l.det_shape}, {n_samples} samples/ray"
            )
        res = GaussNewton(model, p_l, reference=ref_l, verbose=verbose, **solver_kw).solve(
            a0=a0, max_iter=lev.max_iter, tol=lev.tol
        )
        a, prev_mesh = res.a, mesh
        results.append(res)
        models.append(model)
    return results, models
