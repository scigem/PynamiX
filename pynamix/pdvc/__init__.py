"""Projection-based digital volume correlation (P-DVC) with one projection per time step.

A space-time displacement field is fitted directly to a sequence of
radiographs taken while the gantry turns, against a reference volume from a
full scan. See ``docs/pdvc/README.md``.

The geometry, mesh and phantom modules need only numpy and scipy. The
projector and solver need PyTorch (``pip install pynamix[pdvc]``) and are
imported on first use.
"""

from .geometry import ConeBeamGeometry
from .mesh import LinearTimeBasis, RegularHexMesh

__all__ = [
    "ConeBeamGeometry",
    "RegularHexMesh",
    "LinearTimeBasis",
    "DynamicProjectionModel",
    "GaussNewton",
]


def __getattr__(name):
    if name == "DynamicProjectionModel":
        from .model import DynamicProjectionModel

        return DynamicProjectionModel
    if name == "GaussNewton":
        from .solver import GaussNewton

        return GaussNewton
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
