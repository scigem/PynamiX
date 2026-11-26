"""Forward model and Roux-type sensitivities for one-projection-per-step P-DVC.

The deformed state at time tau is described by an Eulerian displacement field
u(x, tau) such that its attenuation is f(x, tau) = f0(x - u(x, tau)), with f0
the reference volume. The displacement is discretised as

    u(x, tau) = sum_k sum_n psi_k(tau) N_n(x) a_{k n}

with N_n the C8 shape functions of a regular hexahedral mesh and psi_k
temporal hat functions. Projection i, acquired at time tau_i with pose i,
is modelled as the line integrals of f(., tau_i) along that pose's rays.

The sensitivity of projection i to a perturbation v of the nodal unknowns is

    (J v)_i = - sum_samples grad f0(x - u(x, tau_i)) . (N v)(x, tau_i) * step

with grad f0 a finite-difference gradient volume, which is the Gauss-Newton
linearisation used in projection-based DVC.
"""

import numpy as np
import torch

from .sampling import sample


def gradient_volume(f0):
    """Central-difference gradient of a (z, y, x) volume, returned as (3, z, y, x) in (x, y, z) order."""
    gz, gy, gx = torch.gradient(f0, dim=(0, 1, 2))
    return torch.stack([gx, gy, gz])


class DynamicProjectionModel:
    def __init__(
        self,
        f0,
        geometry,
        mesh,
        time_basis,
        times,
        n_samples=None,
        fix_first_knot=True,
        chunk=8,
        device="cpu",
        dtype=torch.float32,
        cache_rays=True,
        pad=2,
    ):
        self.device, self.dtype = device, dtype
        self.f0 = torch.as_tensor(np.asarray(f0), device=device, dtype=dtype)
        # Zero margin so that the gradient volume sees the drop to the
        # zero padding used outside the reconstructed region; centred padding
        # keeps voxel coordinates unchanged.
        fp = torch.nn.functional.pad(self.f0, (pad,) * 6)
        self.grad_f0 = gradient_volume(fp)[None]  # (1, 3, D, H, W)
        self.f0_b = fp[None, None]
        self.geom, self.mesh, self.basis = geometry, mesh, time_basis
        self.times = np.asarray(times, dtype=np.float64)
        assert self.times.size == geometry.n_proj
        self.n_samples = int(n_samples or int(np.ceil(np.linalg.norm(f0.shape))))
        self.fix_first = fix_first_knot
        self.n_free = time_basis.n_knots - (1 if fix_first_knot else 0)
        self.param_shape = (self.n_free, 3, *mesh.node_shape)
        self.chunks = [np.arange(s, min(s + chunk, geometry.n_proj)) for s in range(0, geometry.n_proj, chunk)]
        self.psi = torch.as_tensor(time_basis(self.times), device=device, dtype=dtype)
        if fix_first_knot:
            self.psi = self.psi[:, 1:]
        self._rays = {}
        self.cache_rays = cache_rays

    def zeros(self):
        return torch.zeros(self.param_shape, device=self.device, dtype=self.dtype)

    def rays(self, c):
        if c in self._rays:
            return self._rays[c]
        r = self.geom.ray_samples(self.chunks[c], self.f0.shape, self.n_samples, device=self.device, dtype=self.dtype)
        if self.cache_rays:
            self._rays[c] = r
        return r

    def nodal(self, a, c):
        """Nodal displacements (B, 3, nz, ny, nx) for the projections of chunk ``c``."""
        psi = self.psi[self.chunks[c]]
        return torch.einsum("bk,kczyx->bczyx", psi, a)

    def warped_points(self, a, c):
        pts, step = self.rays(c)
        u = self.mesh.interpolate(self.nodal(a, c), pts)
        return pts, pts - u, step

    def project(self, a, c):
        """Simulated projections (B, n_rows, n_cols) of chunk ``c`` for parameters ``a``."""
        _, q, step = self.warped_points(a, c)
        f = sample(self.f0_b, q)[:, 0]
        return f.sum(-1) * step

    def project_all(self, a):
        with torch.no_grad():
            return torch.cat([self.project(a, c) for c in range(len(self.chunks))])

    def linearise(self, a, c):
        """Per-sample weights -grad f0(x - u) * step, shape (B, nr, nc, ns, 3)."""
        pts, q, step = self.warped_points(a, c)
        g = torch.movedim(sample(self.grad_f0, q), 1, -1)
        return pts, -g * step[..., None, None]

    def jvp(self, v, c, lin):
        """J v for chunk ``c`` given the linearisation ``lin`` from :meth:`linearise`."""
        pts, w = lin
        dv = self.mesh.interpolate(self.nodal(v, c), pts)
        return (w * dv).sum((-1, -2))
