"""Cone-beam acquisition geometry, one pose per projection.

Every length is in voxel units of the reference volume. The volume is centred
on the origin; voxel index i along an axis of size N sits at i - (N - 1) / 2.
Axis order of 3-vectors is (x, y, z), with z the rotation axis. Volume arrays
are stored (z, y, x).

Each projection is described by four vectors, as in ASTRA's ``cone_vec``:
source position, detector centre, and the detector's u (column) and v (row)
pixel vectors whose lengths are the pixel pitch. A rotating-gantry scanner
such as the DynaTOM and a rotating-stage scanner give the same relative
geometry, so a calibrated per-projection pose list covers both.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class ConeBeamGeometry:
    src: np.ndarray  # (n_proj, 3)
    det: np.ndarray  # (n_proj, 3) detector centre
    u: np.ndarray  # (n_proj, 3) column pixel vector
    v: np.ndarray  # (n_proj, 3) row pixel vector
    det_shape: tuple  # (n_rows, n_cols)

    @property
    def n_proj(self):
        return self.src.shape[0]

    @classmethod
    def circular(cls, angles, sod, sdd, det_shape, pixel_size, det_offset=(0.0, 0.0)):
        """Circular trajectory about z. ``angles`` in radians, one per projection.

        ``sod`` source-to-axis distance, ``sdd`` source-to-detector distance,
        ``pixel_size`` detector pitch, all in voxel units. ``det_offset`` shifts
        the detector centre by (columns, rows) in pixels.
        """
        a = np.asarray(angles, dtype=np.float64)
        e = np.stack([np.cos(a), np.sin(a), np.zeros_like(a)], axis=1)
        eu = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=1)
        ev = np.tile([0.0, 0.0, 1.0], (a.size, 1))
        src = -sod * e
        det = (sdd - sod) * e + det_offset[0] * pixel_size * eu + det_offset[1] * pixel_size * ev
        return cls(src, det, pixel_size * eu, pixel_size * ev, tuple(det_shape))

    def subset(self, idx):
        idx = np.atleast_1d(idx)
        return ConeBeamGeometry(self.src[idx], self.det[idx], self.u[idx], self.v[idx], self.det_shape)

    def binned(self, factor):
        """Same poses with the detector binned by an integer factor.

        Rows and columns that do not fill a bin are cropped as in
        :func:`bin_projections`, and the detector centre moves with the crop.
        """
        nr, nc = self.det_shape
        nrb, ncb = nr // factor, nc // factor
        r0, c0 = (nr - nrb * factor) // 2, (nc - ncb * factor) // 2
        dr = r0 + (nrb * factor - 1) / 2 - (nr - 1) / 2
        dc = c0 + (ncb * factor - 1) / 2 - (nc - 1) / 2
        det = self.det + dc * self.u + dr * self.v
        return ConeBeamGeometry(self.src, det, self.u * factor, self.v * factor, (nrb, ncb))

    def magnification(self):
        """Per-projection magnification at the rotation axis, |det - src| / |src|."""
        return np.linalg.norm(self.det - self.src, axis=1) / np.linalg.norm(self.src, axis=1)

    def pixel_pitch(self):
        """Mean detector pixel pitch (voxel units)."""
        return float(np.mean(np.linalg.norm(self.u, axis=1)))

    def ray_samples(self, idx, vol_shape, n_samples, device="cpu", dtype=None):
        """Sample points along every ray of the projections ``idx``.

        Returns ``points`` of shape (B, n_rows, n_cols, n_samples, 3) in voxel
        coordinates (x, y, z) and ``step`` of shape (B, n_rows, n_cols), the
        spacing between samples. Rays are clipped to the sphere enclosing the
        volume and sampled at segment midpoints, so a line integral is
        ``sum(f(points), -1) * step``.
        """
        import torch

        dtype = dtype or torch.float32
        idx = np.atleast_1d(idx)
        nr, nc = self.det_shape

        def t(x):
            return torch.as_tensor(x[idx], device=device, dtype=torch.float64)

        src, det, u, v = t(self.src), t(self.det), t(self.u), t(self.v)

        cols = torch.arange(nc, device=device, dtype=torch.float64) - (nc - 1) / 2
        rows = torch.arange(nr, device=device, dtype=torch.float64) - (nr - 1) / 2
        pix = (
            det[:, None, None, :]
            + cols[None, None, :, None] * u[:, None, None, :]
            + rows[None, :, None, None] * v[:, None, None, :]
        )  # (B, nr, nc, 3)
        s = src[:, None, None, :]
        d = pix - s
        d = d / torch.linalg.norm(d, dim=-1, keepdim=True)

        half = (torch.tensor(vol_shape[::-1], dtype=torch.float64, device=device) - 1) / 2
        radius = torch.linalg.norm(half) + 1.0
        b = (s * d).sum(-1)
        c = (s * s).sum(-1) - radius**2
        disc = torch.clamp(b * b - c, min=0.0)
        t0 = -b - torch.sqrt(disc)
        t1 = -b + torch.sqrt(disc)
        step = (t1 - t0) / n_samples
        k = torch.arange(n_samples, device=device, dtype=torch.float64) + 0.5
        tk = t0[..., None] + k * step[..., None]
        points = s[..., None, :] + tk[..., None] * d[..., None, :]
        return points.to(dtype), step.to(dtype)


def bin_projections(p, factor):
    """Average ``factor`` x ``factor`` pixel blocks of projections (n, nr, nc).

    Rows and columns that do not fill a bin are cropped symmetrically, matching
    :meth:`ConeBeamGeometry.binned`.
    """
    if factor == 1:
        return p
    n, nr, nc = p.shape
    nrb, ncb = nr // factor, nc // factor
    r0, c0 = (nr - nrb * factor) // 2, (nc - ncb * factor) // 2
    q = p[:, r0 : r0 + nrb * factor, c0 : c0 + ncb * factor]
    return q.reshape(n, nrb, factor, ncb, factor).mean((2, 4))
