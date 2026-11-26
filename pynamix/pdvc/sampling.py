"""Trilinear sampling of voxel and nodal arrays at physical points."""

import torch
import torch.nn.functional as F


def normalise(points, shape):
    """Voxel coordinates (x, y, z) centred on the origin -> grid_sample's [-1, 1].

    ``shape`` is the (z, y, x) shape of the array being sampled; its first and
    last entries along each axis map to -1 and +1 (align_corners=True).
    """
    half = torch.tensor(shape[::-1], dtype=points.dtype, device=points.device)
    half = torch.clamp((half - 1) / 2, min=1e-12)
    return points / half


def sample(array, points, extent=None, padding="zeros"):
    """Trilinear interpolation of ``array`` at ``points``.

    ``array``: (B, C, D, H, W) or broadcastable batch 1.
    ``points``: (B, ..., 3) in voxel coordinates of the reference volume.
    ``extent``: (z, y, x) size of the reference volume over which the array's
    nodes are spread. Defaults to the array's own shape (a voxel array).
    Returns (B, C, ...).
    """
    B = points.shape[0]
    lead = points.shape[1:-1]
    grid = normalise(points, extent if extent is not None else tuple(array.shape[2:]))
    grid = grid.reshape(B, -1, 1, 1, 3)
    if array.shape[0] != B:
        array = array.expand(B, *array.shape[1:])
    out = F.grid_sample(array, grid, mode="bilinear", padding_mode=padding, align_corners=True)
    return out.reshape(B, array.shape[1], *lead)
