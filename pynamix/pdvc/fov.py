"""Material outside the reconstructed region of the reference volume.

A specimen usually extends beyond the scanner's field of view, and in a cone
beam the reliably reconstructed region is smaller than the volume array. The
reference volume f0 is therefore known only in a region K. A ray can pick up
material that is not in K in two ways: it crosses material that lies outside
K at all times (above or below the field of view, or in the poorly
reconstructed cone-beam corners), or deformation carries material across the
boundary of K, so that f0 is needed at a reference position x - u outside K.

Two remedies are provided, and can be combined.

Masking. A detector pixel is dropped from the fit when its ray samples a
point whose reference position lies outside K but inside an envelope where
material may be (by default a vertical cylinder of given radius around the
rotation axis). The mask depends on the displacement and is recomputed at
every Gauss-Newton iteration. A dropped pixel stays dropped until
:meth:`FieldOfView.reset`, so the mask can only grow during a solve and
cannot oscillate with the iterates. Without difference data, every ray that
meets the specimen outside K is dropped.

Reference difference. With projections p_ref of the undeformed specimen at
the same poses, the solver fits P_i[a] - P_i[0] to p_i - p_ref,i. Material
outside K that has not moved cancels from this difference, so a pixel need
only be dropped when its ray samples unknown material that has moved by more
than ``u_threshold``. Outside the mesh the displacement is the value on the
mesh boundary held constant, so this test extrapolates the boundary motion. Reconstruction errors of f0 inside K also cancel to
first order. This needs a reference scan whose poses match the dynamic
sequence, or reference projections interpolated in angle.
"""

import numpy as np
import torch
import torch.nn.functional as F

from .sampling import sample


def cone_beam_known(geometry, vol_shape, margin=1.0, stride=4, every=1):
    """Voxels seen by every projection of ``geometry`` (bool array, (z, y, x)).

    A voxel counts as seen when its centre projects inside the detector by at
    least ``margin`` pixels. For a circular trajectory this is the region a
    full scan reconstructs without truncation. The test runs on a grid of
    ``stride`` voxels and is expanded by nearest-neighbour; ``every`` thins
    the poses checked.
    """
    vol_shape = tuple(vol_shape)
    axes = [np.arange(0, n, stride) - (n - 1) / 2 for n in vol_shape]
    Z, Y, X = np.meshgrid(*axes, indexing="ij")
    P = np.stack([X, Y, Z], -1).reshape(-1, 3)
    ok = np.ones(P.shape[0], dtype=bool)
    nr, nc = geometry.det_shape
    for i in range(0, geometry.n_proj, every):
        s, d, u, v = geometry.src[i], geometry.det[i], geometry.u[i], geometry.v[i]
        n = np.cross(u, v)
        t = ((d - s) @ n) / ((P - s) @ n)
        hit = s + t[:, None] * (P - s) - d
        col = hit @ u / (u @ u)
        row = hit @ v / (v @ v)
        ok &= (np.abs(col) <= (nc - 1) / 2 - margin) & (np.abs(row) <= (nr - 1) / 2 - margin)
    ok = ok.reshape(Z.shape)
    for ax, n in enumerate(vol_shape):
        ok = np.repeat(ok, stride, axis=ax)[(slice(None),) * ax + (slice(0, n),)]
    return ok


class FieldOfView:
    """Detector-pixel masking for material outside the known reference region.

    ``known``: bool array over the reference volume, True where f0 is
    reliable; None means the whole array. Outside the array nothing is known.
    ``envelope_radius``: radius (voxels) of the vertical cylinder around the
    rotation axis where material may lie; outside it the scene is taken as
    empty. None means material may be anywhere.
    ``difference``: True when the solver is given reference projections, so
    only moving unknown material forces a pixel out.
    ``u_threshold``: displacement (voxels) above which unknown material counts
    as moving.
    ``dilate``: pixels added around every masked pixel.
    """

    def __init__(self, known=None, envelope_radius=None, difference=False, u_threshold=0.25, dilate=1):
        self.known = known
        self.envelope_radius = envelope_radius
        self.difference = difference
        self.u_threshold = u_threshold
        self.dilate = int(dilate)
        self._known_t = None
        self._dropped = {}

    def reset(self):
        """Forget pixels dropped so far (called at the start of each solve)."""
        self._dropped = {}

    def _known_tensor(self, model):
        if self._known_t is None:
            k = np.ones(model.f0.shape) if self.known is None else np.asarray(self.known, dtype=float)
            self._known_t = torch.as_tensor(k, device=model.device, dtype=model.dtype)[None, None]
        return self._known_t

    def _unknown_material(self, model, q):
        """True where a reference position ``q`` may hold material not in f0."""
        kappa = sample(self._known_tensor(model), q)[:, 0]
        bad = kappa < 0.999
        if self.envelope_radius is not None:
            bad &= (q[..., 0] ** 2 + q[..., 1] ** 2) < self.envelope_radius**2
        return bad

    @torch.no_grad()
    def ray_weights(self, model, a, c):
        """(B, nr, nc) weights, 1 for pixels kept and 0 for pixels dropped."""
        pts, step = model.rays(c)
        u = model.mesh.interpolate(model.nodal(a, c), pts)
        q = pts - u
        bad = self._unknown_material(model, q)
        if self.difference:
            bad |= self._unknown_material(model, pts)
            bad &= torch.linalg.norm(u, dim=-1) > self.u_threshold
        bad = bad.any(-1).to(model.dtype)  # (B, nr, nc)
        if self.dilate > 0:
            k = 2 * self.dilate + 1
            bad = F.max_pool2d(bad[:, None], k, stride=1, padding=self.dilate)[:, 0]
        if c in self._dropped and self._dropped[c].shape == bad.shape:
            bad = torch.maximum(bad, self._dropped[c])
        self._dropped[c] = bad
        return 1.0 - bad
