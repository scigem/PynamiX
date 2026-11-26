"""Synthetic reference volumes and imposed space-time displacement fields."""

import numpy as np
from scipy import ndimage


def granular_phantom(shape, grain_size=3.0, porosity=0.4, blur=1.5, seed=0):
    """Textured two-phase volume inside a finite vertical cylinder, values in [0, 1].

    The cylinder has radius and half-height 0.4 of the volume edge, so the
    specimen sits wholly inside the field of view.

    Smoothed white noise thresholded at the ``porosity`` quantile gives
    grain-like features of size ~``grain_size`` voxels; a final blur of
    ``blur`` voxels removes the staircase so the gradient volume is usable.
    """
    rng = np.random.default_rng(seed)
    n = ndimage.gaussian_filter(rng.standard_normal(shape), grain_size / 2)
    solid = (n > np.quantile(n, porosity)).astype(np.float64)
    z, y, x = np.meshgrid(*[np.arange(s) - (s - 1) / 2 for s in shape], indexing="ij")
    r = 0.4 * min(shape[1], shape[2])
    solid *= ((x**2 + y**2) < r**2) & (np.abs(z) < 0.4 * shape[0])
    return ndimage.gaussian_filter(solid, blur)


def compression_shear(points, tau, strain=0.08, shear=0.06, half_height=1.0):
    """Eulerian displacement of a column under growing compression and simple shear.

    ``points`` (..., 3) in voxel coordinates; ``tau`` in [0, 1]. Vertical strain
    and shear strain both grow as tau**1.5, with a lateral Poisson-type bulge.
    Returns (..., 3).
    """
    s = tau**1.5
    x, y, z = points[..., 0], points[..., 1], points[..., 2]
    zn = z / half_height
    ux = 0.3 * strain * s * x + shear * s * z
    uy = 0.3 * strain * s * y
    uz = -strain * s * z * (1 + 0.2 * zn**2)
    return np.stack([ux, uy, uz], axis=-1)


def sphere_packing(radius, half_height, diameters, counts, weights, rng, batch=2000, max_batches=400, contained=False):
    """Random sequential addition of non-overlapping spheres in a vertical cylinder.

    Species are placed in the order given, ``counts[a]`` spheres of diameter
    ``diameters[a]``. A candidate centre of species a is drawn uniformly in
    the cylinder and kept with probability ``weights[a](x)`` (values in
    [0, 1]), so the local number density follows that function until the
    packing crowds. With ``contained`` False the centres fill the whole
    cylinder and grains cross its surface, so the mean solid fraction is
    uniform up to the surface with no depleted layer at the wall; with True
    every grain lies inside. Returns centres (n, 3) in (x, y, z) and
    diameters (n,).
    """
    from scipy.spatial import cKDTree

    centres = np.empty((0, 3))
    diam = np.empty(0)
    dmax = max(diameters)
    for d, n_target, wfun in zip(diameters, counts, weights):
        a = d / 2
        m = a if contained else 0.0
        placed = 0
        for _ in range(max_batches):
            if placed >= n_target:
                break
            r = (radius - m) * np.sqrt(rng.random(batch))
            th = 2 * np.pi * rng.random(batch)
            z = (half_height - m) * (2 * rng.random(batch) - 1)
            c = np.stack([r * np.cos(th), r * np.sin(th), z], 1)
            c = c[rng.random(batch) < wfun(c)]
            if centres.shape[0] and c.shape[0]:
                tree = cKDTree(centres)
                dist, j = tree.query(c, k=min(8, centres.shape[0]), distance_upper_bound=a + dmax / 2)
                dist = np.atleast_2d(dist.T).T if dist.ndim == 1 else dist
                j = np.atleast_2d(j.T).T if j.ndim == 1 else j
                ok = np.ones(c.shape[0], dtype=bool)
                hit = np.isfinite(dist)
                jj = np.where(hit, j, 0)
                ok &= ~np.any(hit & (dist < a + diam[jj] / 2), axis=1)
                c = c[ok]
            keep = []
            for k in range(c.shape[0]):  # greedy within the batch
                if placed + len(keep) >= n_target:
                    break
                if not keep or np.min(np.sum((c[keep] - c[k]) ** 2, 1)) >= d * d:
                    keep.append(k)
            centres = np.vstack([centres, c[keep]])
            diam = np.concatenate([diam, np.full(len(keep), d)])
            placed += len(keep)
    return centres, diam


def project_spheres(centres, diameters, geometry, i, supersample=2):
    """Exact cone-beam path length through spheres for projection ``i``.

    Each detector pixel is the mean of ``supersample``^2 rays across its area,
    which reproduces area integration over the pixel up to aliasing. Returns
    (n_rows, n_cols) in voxel units.
    """
    import torch

    nr, nc = geometry.det_shape
    s = torch.as_tensor(geometry.src[i])
    det, u, v = (torch.as_tensor(x[i]) for x in (geometry.det, geometry.u, geometry.v))
    nrm = torch.linalg.cross(u, v)
    nrm = nrm / torch.linalg.norm(nrm)
    out = torch.zeros(nr * nc, dtype=torch.float64)
    sub = (torch.arange(supersample, dtype=torch.float64) + 0.5) / supersample - 0.5
    for d in np.unique(diameters):
        c = torch.as_tensor(centres[diameters == d])
        a = float(d) / 2
        cs = c - s
        t = ((det - s) @ nrm) / (cs @ nrm)
        hitp = s + t[:, None] * cs - det
        col = hitp @ u / (u @ u)
        row = hitp @ v / (v @ v)
        M = float(t.max()) / float(torch.linalg.norm(u))  # largest magnification, in pixels per voxel
        h = int(np.ceil(a * M + 2))
        o = torch.arange(-h, h + 1, dtype=torch.float64)
        # pixel indices around the projected centre, then back to centred coordinates
        ic = torch.round(col + (nc - 1) / 2)[:, None, None] + o[None, None, :]  # (n, 1, w)
        ir = torch.round(row + (nr - 1) / 2)[:, None, None] + o[None, :, None]  # (n, w, 1)
        pc, pr = ic - (nc - 1) / 2, ir - (nr - 1) / 2
        acc = 0
        for du in sub:
            for dv in sub:
                Q = det + (pc + du)[..., None] * u + (pr + dv)[..., None] * v  # (n, w, w, 3)
                dirn = Q - s
                dirn = dirn / torch.linalg.norm(dirn, dim=-1, keepdim=True)
                proj = (cs[:, None, None, :] * dirn).sum(-1)
                b2 = (cs * cs).sum(-1)[:, None, None] - proj**2
                acc = acc + 2 * torch.sqrt(torch.clamp(a * a - b2, min=0))
        acc = acc / supersample**2
        ic = ic.expand_as(acc).long()
        ir = ir.expand_as(acc).long()
        ok = (ic >= 0) & (ic < nc) & (ir >= 0) & (ir < nr)
        out.index_add_(0, (ir * nc + ic)[ok], acc[ok])
    return out.reshape(nr, nc).numpy()
