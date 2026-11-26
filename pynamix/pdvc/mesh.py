"""Regular 8-node hexahedral mesh and temporal basis functions.

Trilinear (C8) shape functions on a regular node grid are exactly trilinear
interpolation of the nodal values, so evaluating a nodal displacement field at
arbitrary points is a single grid_sample.
"""

import numpy as np


class RegularHexMesh:
    """Node grid spanning the reference volume, ``n_elem`` elements per axis (z, y, x)."""

    def __init__(self, vol_shape, n_elem):
        self.vol_shape = tuple(vol_shape)
        self.n_elem = tuple(int(n) for n in np.broadcast_to(n_elem, 3))
        self.node_shape = tuple(n + 1 for n in self.n_elem)

    @property
    def n_nodes(self):
        return int(np.prod(self.node_shape))

    def node_coords(self):
        """(z, y, x, 3) array of node positions in voxel coordinates (x, y, z)."""
        axes = [np.linspace(-(N - 1) / 2, (N - 1) / 2, n) for N, n in zip(self.vol_shape, self.node_shape)]
        z, y, x = np.meshgrid(*axes, indexing="ij")
        return np.stack([x, y, z], axis=-1)

    def interpolate(self, U, points):
        """Displacement at ``points``.

        ``U``: (B, 3, nz, ny, nx) nodal displacements (components x, y, z).
        ``points``: (B, ..., 3). Returns (B, ..., 3). Outside the mesh the
        border value is held.
        """
        import torch

        from .sampling import sample

        out = sample(U, points, extent=self.vol_shape, padding="border")
        return torch.movedim(out, 1, -1)

    def prolongate(self, a, fine):
        """Nodal field ``a`` (K, 3, nz, ny, nx) on this mesh evaluated at the nodes of ``fine``.

        Exact when every coarse element boundary is also a fine one.
        """
        import torch

        from .sampling import sample

        X = torch.as_tensor(fine.node_coords(), dtype=a.dtype, device=a.device)
        X = X[None].expand(a.shape[0], *X.shape)
        return sample(a, X, extent=self.vol_shape, padding="border")


class LinearTimeBasis:
    """Piecewise-linear hat functions on ``n_knots`` equally spaced knots in [0, 1]."""

    def __init__(self, n_knots):
        self.n_knots = int(n_knots)

    def __call__(self, tau):
        tau = np.clip(np.asarray(tau, dtype=np.float64), 0.0, 1.0)
        if self.n_knots == 1:
            return np.ones((tau.size, 1))
        knots = np.linspace(0.0, 1.0, self.n_knots)
        h = knots[1] - knots[0]
        return np.clip(1.0 - np.abs(tau[:, None] - knots[None, :]) / h, 0.0, None)
