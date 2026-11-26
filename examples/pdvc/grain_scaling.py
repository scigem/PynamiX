"""Counting-statistics floor for species fields on a C8 mesh.

The realised volume fraction of species a is a sum of grains, so at mesh
scale it fluctuates. If grain centres are independent, the L2 projection of
the grain field onto the mesh has, at a point x, variance

    phi_a v_a N(x)^T M^-1 N(x)

with N(x) the shape functions at x and M the mass matrix over the specimen.
Averaged over the specimen this is phi_a^2 c / n_a, where n_a is the number
of grains of species a in one element and c depends only on the mesh and the
specimen shape. This script computes c for the meshes of the study and
compares the predicted floor with the errors in a study's results.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from pynamix.pdvc import RegularHexMesh
from pynamix.spectral.tomography import trilinear_weights


def mesh_constant(R, H, margin, elem, h=4.0):
    """c = <N^T M^-1 N> times the element volume, over the cylinder of radius R and half-height H."""
    vol = (int(2 * H + margin) + 1, int(2 * R + margin) + 1, int(2 * R + margin) + 1)
    mesh = RegularHexMesh(vol, elem)
    ax = np.arange(-R, R + h / 2, h)
    az = np.arange(-H, H + h / 2, h)
    Z, Y, X = np.meshgrid(az, ax, ax, indexing="ij")
    q = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    q = q[np.hypot(q[:, 0], q[:, 1]) < R]
    N = trilinear_weights(q, mesh.vol_shape, mesh.node_shape)
    M = (N.T @ N).toarray() * h**3
    M += 1e-9 * np.mean(np.diag(M)) * np.eye(len(M))
    Minv = np.linalg.inv(M)
    var = np.einsum("pi,ij,pj->p", N.toarray(), Minv, N.toarray())  # per unit (phi v)
    elem_vol = float(np.prod([(n - 1) / e for n, e in zip(mesh.vol_shape, mesh.n_elem)]))
    return float(var.mean() * elem_vol)


def s0(ds, dl, phi, c=0.5):
    """Percus-Yevick S_ss(0) and S_ll(0) at solid fraction phi, half the solid in each species."""
    from pynamix.spectral.closure import py_partials

    v = np.pi * np.array([ds, dl]) ** 3 / 6
    r = py_partials({0: ds, 1: dl}, {0: phi * (1 - c) / v[0], 1: phi * c / v[1]}, np.array([1e-3]), return_result=True)
    return [float(r.S[(0, 0)][0]), float(r.S[(1, 1)][0])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", type=Path, nargs="+", help="study.json files")
    ap.add_argument("--plot", type=Path, default=None, help="write a figure of error against grains per element")
    a = ap.parse_args()
    rows = []
    for path in a.results:
        res = json.loads(path.read_text())
        for key, r in res.items():
            if not (key.startswith("elem=") or key.startswith("grain=")):
                continue
            args = r["args"]
            R, H, dl = args["radius"], args["half_height"], args["d_large"]
            c = mesh_constant(R, H, dl, args["elem"])
            n = r["grains_per_element"]
            phi_bar = args["phi"] / 2  # each species holds about half the solid on average
            S0 = s0(args["d_small"], dl, args["phi"])
            # Poisson count variance reduced by the partial structure factor at q -> 0
            floor = {s: phi_bar * np.sqrt(S0[j] * c / n[s]) for j, s in enumerate(("small", "large"))}
            e = r["error_vs_realised"]
            rows.append(
                (key, c, n["small"], n["large"], floor["small"], e["phi_small"], floor["large"], e["phi_large"])
            )
    print(f"{'run':12s} {'c':>6s} {'n_small':>8s} {'n_large':>8s}  floor/err small   floor/err large")
    for k, c, ns, nl, fs, es, fl, el in sorted(rows, key=lambda r: r[3]):
        print(f"{k:12s} {c:6.2f} {ns:8.0f} {nl:8.0f}  {fs:.4f}/{es:.4f}    {fl:.4f}/{el:.4f}")
    if a.plot:
        plot(rows, a.plot)


def plot(rows, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    for j, (name, cols) in enumerate([("small grains", (2, 4, 5)), ("large grains", (3, 6, 7))]):
        for prefix, marker, lab in [("elem=", "o", "mesh varied"), ("grain=", "s", "grain size varied")]:
            sel = sorted((r for r in rows if r[0].startswith(prefix)), key=lambda r: r[cols[0]])
            n = [r[cols[0]] for r in sel]
            ax[j].loglog(
                n,
                [r[cols[2]] for r in sel],
                marker,
                color="C0",
                mfc="none" if "grain" in prefix else "C0",
                label=f"measured, {lab}",
            )
            ax[j].loglog(
                n,
                [r[cols[1]] for r in sel],
                marker + "--",
                color="C3",
                mfc="none",
                alpha=0.7,
                label=f"counting estimate, {lab}",
            )
        ax[j].set_xlabel(f"{name} per element")
        ax[j].set_ylabel(f"rms error in volume fraction of {name}")
        ax[j].legend(fontsize=7)
    fig.savefig(path, dpi=130)


if __name__ == "__main__":
    main()
