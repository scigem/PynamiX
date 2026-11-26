"""Specimen taller than the field of view, with material flowing into it.

A textured column extends beyond the reconstructed slab at the top and
bottom. The reference volume f0 is the slab, restricted to the region a full
circular scan sees in every projection (``cone_beam_known``). The column is
compressed about its mid-height, so material above and below the slab moves
into the field of view while the gantry turns once. Three treatments are
compared:

    naive   everything outside f0 is taken as empty
    mask    pixels whose rays meet unknown material are dropped
    diff    reference projections of the undeformed column at the same poses
            are subtracted, and only rays meeting unknown material that moves
            are dropped
"""

import argparse
import time
from pathlib import Path

import numpy as np
import torch

from pynamix.pdvc import ConeBeamGeometry, DynamicProjectionModel, GaussNewton, LinearTimeBasis, RegularHexMesh
from pynamix.pdvc.fov import FieldOfView, cone_beam_known
from pynamix.pdvc.phantom import compression_shear, granular_phantom


def nodal_truth(mesh, times, half_height, **kw):
    X = mesh.node_coords()
    U = np.stack([compression_shear(X, t, half_height=half_height, **kw) for t in times])
    return np.moveaxis(U, -1, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=["naive", "mask", "diff"], default="diff")
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--n-proj", type=int, default=90)
    ap.add_argument("--elem", type=int, default=4)
    ap.add_argument("--knots", type=int, default=5)
    ap.add_argument("--noise", type=float, default=0.01)
    ap.add_argument("--strain", type=float, default=0.08)
    ap.add_argument("--shear", type=float, default=0.04)
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/fov"))
    args = ap.parse_args()
    torch.manual_seed(0)
    N = args.n
    slab = (N, N, N)
    tall_shape = (2 * N, N, N)
    half = (N - 1) / 2

    # column of radius 0.4 N and half-height 0.8 N: about 0.3 N beyond the slab at each end
    tall = granular_phantom(tall_shape, seed=1)
    angles = np.linspace(0, 2 * np.pi, args.n_proj, endpoint=False)
    times = np.linspace(0, 1, args.n_proj)
    det = (int(1.5 * N), 2 * N)  # vertical field of view truncates the column; lateral does not
    geom = ConeBeamGeometry.circular(angles, sod=4 * N, sdd=6 * N, det_shape=det, pixel_size=1.0)

    known = cone_beam_known(geom, slab, stride=2)
    f0 = tall[N // 2 : N // 2 + N] * known
    print(f"slab {slab} of a {tall_shape} column; known region {known.mean():.1%} of the slab")

    kw = dict(strain=args.strain, shear=args.shear)
    fine = RegularHexMesh(tall_shape, (32, 16, 16))
    gen = DynamicProjectionModel(
        tall,
        geom,
        fine,
        LinearTimeBasis(args.n_proj),
        times,
        n_samples=6 * N,
        fix_first_knot=False,
        cache_rays=False,
    )
    with torch.no_grad():
        U = torch.as_tensor(nodal_truth(fine, times, half, **kw), dtype=torch.float32)
        p = gen.project_all(U)
        p_ref = gen.project_all(torch.zeros_like(U))
    noise_std = args.noise * float(p.std())
    p = p + noise_std * torch.randn_like(p)
    p_ref = p_ref + noise_std * torch.randn_like(p_ref)
    print(f"max |u| in column {np.abs(U.numpy()).max():.2f} vox, case {args.case}")

    mesh = RegularHexMesh(slab, args.elem)
    basis = LinearTimeBasis(args.knots)
    model = DynamicProjectionModel(f0, geom, mesh, basis, times, n_samples=2 * N)
    fov, ref = None, None
    if args.case == "mask":
        fov = FieldOfView(known=known, envelope_radius=0.45 * N)
    elif args.case == "diff":
        fov = FieldOfView(known=known, envelope_radius=0.45 * N, difference=True)
        ref = p_ref
    t0 = time.time()
    res = GaussNewton(model, p, fov=fov, reference=ref).solve(max_iter=15, tol=1e-3)
    print(f"solved in {time.time() - t0:.0f} s")

    psi = basis(times)[:, 1:]
    U_rec = np.einsum("tk,kczyx->tczyx", psi, res.a.numpy())
    U_true = nodal_truth(mesh, times, half, **kw)
    X = mesh.node_coords()
    in_col = (X[..., 0] ** 2 + X[..., 1] ** 2) < (0.4 * N) ** 2
    interior = in_col & (np.abs(X[..., 2]) < 0.3 * N)
    err = U_rec - U_true
    for name, m in [("nodes in column", in_col), ("interior nodes", interior)]:
        e, ef, uf = err[..., m], err[-1][..., m], U_true[-1][..., m]
        print(
            f"{name:16s} ({m.sum():3d}) rms error {np.sqrt((e**2).mean()):.3f} vox, "
            f"final state {np.sqrt((ef**2).mean()):.3f} of imposed {np.sqrt((uf**2).mean()):.3f}"
        )
    h = res.history[-1]
    print(
        f"noise std {noise_std:.3e}, final rms residual {h['rms_residual']:.3e}, "
        f"pixels kept {h.get('kept_fraction', 1.0):.1%}"
    )

    args.out.mkdir(parents=True, exist_ok=True)
    np.savez(args.out / f"{args.case}.npz", U_rec=U_rec, U_true=U_true, times=times)


if __name__ == "__main__":
    main()
