"""Synthetic one-projection-per-step P-DVC test.

A textured column deforms under growing compression and shear while the
gantry turns once; one projection is taken per time step. The measured
sequence is generated with a fine mesh, a temporal knot at every projection
and twice the ray sampling, so the coarse space-time model used for the
reconstruction does not contain the imposed field exactly. The script reports
the error of the recovered nodal displacements against the imposed field.
"""

import argparse
import time
from pathlib import Path

import numpy as np
import torch

from pynamix.pdvc import ConeBeamGeometry, DynamicProjectionModel, GaussNewton, LinearTimeBasis, RegularHexMesh
from pynamix.pdvc.multiscale import default_levels, multiscale_solve
from pynamix.pdvc.phantom import compression_shear, granular_phantom


def nodal_truth(mesh, times, half_height, **kw):
    X = mesh.node_coords()  # (z, y, x, 3)
    U = np.stack([compression_shear(X, t, half_height=half_height, **kw) for t in times])
    return np.moveaxis(U, -1, 1)  # (n_t, 3, z, y, x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=32, help="volume edge length (voxels)")
    ap.add_argument("--n-proj", type=int, default=120)
    ap.add_argument("--elem", type=int, default=4, help="elements per axis")
    ap.add_argument("--knots", type=int, default=6, help="temporal knots")
    ap.add_argument("--noise", type=float, default=0.01, help="noise std / projection std")
    ap.add_argument("--strain", type=float, default=0.08)
    ap.add_argument("--shear", type=float, default=0.06)
    ap.add_argument("--reg", type=float, default=1e-2, help="relative gradient penalty")
    ap.add_argument(
        "--drift", type=float, default=0.0, help="amplitude of synthetic grey-level drift (fraction of mean projection)"
    )
    ap.add_argument("--grey", default="", help="grey-level terms to fit, e.g. offset,gain,tilt,hardening")
    ap.add_argument("--multiscale", action="store_true", help="coarse-to-fine levels")
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/synthetic"))
    args = ap.parse_args()
    torch.manual_seed(0)
    N = args.n
    shape = (N, N, N)
    half = (N - 1) / 2

    f0 = granular_phantom(shape, seed=1)
    angles = np.linspace(0, 2 * np.pi, args.n_proj, endpoint=False)
    times = np.linspace(0, 1, args.n_proj)
    det = (int(1.5 * N), int(1.5 * N))
    geom = ConeBeamGeometry.circular(angles, sod=4 * N, sdd=6 * N, det_shape=det, pixel_size=1.0)

    # measured data: fine mesh, a knot per projection, denser ray sampling
    fine = RegularHexMesh(shape, 16)
    U_fine = nodal_truth(fine, times, half, strain=args.strain, shear=args.shear)
    gen = DynamicProjectionModel(
        f0,
        geom,
        fine,
        LinearTimeBasis(args.n_proj),
        times,
        n_samples=4 * N,
        fix_first_knot=False,
        cache_rays=False,
    )
    with torch.no_grad():
        p = gen.project_all(torch.as_tensor(U_fine, dtype=torch.float32))
    if args.drift > 0:
        # flux / gain drift, additive offset drift, a ramp across the detector
        # and a global beam-hardening-like compression of long path lengths
        A, pm, pmax = args.drift, float(p.mean()), float(p.max())
        tt = torch.as_tensor(times, dtype=p.dtype)[:, None, None]
        s_ = torch.linspace(-1, 1, det[1], dtype=p.dtype)[None, None, :]
        gain = 1 + A * (tt + 0.5 * torch.sin(9 * tt))
        offset = A * pm * torch.cos(5 * tt)
        tilt = 0.5 * A * pm * torch.sin(3 * tt) * s_
        p = gain * (p - 0.5 * A * p**2 / pmax) + offset + tilt
    noise_std = args.noise * float(p.std())
    p = p + noise_std * torch.randn_like(p)
    umax = np.abs(U_fine).max()
    print(f"volume {shape}, detector {det}, {args.n_proj} projections, max |u| {umax:.2f} vox")

    # reconstruction model
    mesh = RegularHexMesh(shape, args.elem)
    basis = LinearTimeBasis(args.knots)
    grey = [t for t in args.grey.split(",") if t] or None
    n_unk = mesh.n_nodes * 3 * (args.knots - 1)
    print(
        f"unknowns: {n_unk} ({mesh.n_nodes} nodes x 3 x {args.knots - 1} knots)"
        + (f", grey terms {grey} eliminated per projection" if grey else "")
    )

    t0 = time.time()
    if args.multiscale:
        results, _ = multiscale_solve(
            f0,
            geom,
            p,
            basis,
            times,
            levels=default_levels(args.elem),
            solver_kw=dict(reg=args.reg, grey=grey),
        )
        res = results[-1]
        history = [h for r in results for h in r.history]
    else:
        model = DynamicProjectionModel(f0, geom, mesh, basis, times, n_samples=2 * N)
        res = GaussNewton(model, p, reg=args.reg, grey=grey).solve(max_iter=15, tol=1e-3)
        history = res.history
    print(f"solved in {time.time() - t0:.0f} s")

    # compare recovered and imposed nodal fields at every projection time
    a = res.a.numpy()
    psi = basis(times)[:, 1:]
    U_rec = np.einsum("tk,kczyx->tczyx", psi, a)
    U_true = nodal_truth(mesh, times, half, strain=args.strain, shear=args.shear)
    X = mesh.node_coords()
    # nodes within the reference specimen (cylinder of radius and half-height 0.4 N)
    inside = ((X[..., 0] ** 2 + X[..., 1] ** 2) < (0.4 * N) ** 2) & (np.abs(X[..., 2]) < 0.4 * N)
    err = U_rec - U_true
    for name, m in [("all nodes", np.ones_like(inside)), ("nodes in sample", inside)]:
        e = err[..., m]
        print(f"{name:16s} rms error {np.sqrt((e**2).mean()):.3f} vox, max {np.abs(e).max():.3f} vox")
    e_end = err[-1][..., inside]
    print(
        f"final state, nodes in sample: rms error {np.sqrt((e_end**2).mean()):.3f} vox "
        f"(imposed rms {np.sqrt((U_true[-1][..., inside]**2).mean()):.3f} vox)"
    )
    print(f"noise std {noise_std:.3e}, final rms residual {history[-1]['rms_residual']:.3e}")

    args.out.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out / "result.npz",
        a=a,
        U_rec=U_rec,
        U_true=U_true,
        times=times,
        history=np.array([[h["rms_residual"], h["max_update"]] for h in history]),
        grey=res.grey.numpy() if res.grey is not None else np.zeros(0),
    )
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        iz, iy, ix = [s - 1 for s in mesh.node_shape]  # top corner of the mesh
        jz, jy, jx = iz, mesh.node_shape[1] // 2, mesh.node_shape[2] // 2  # top centre
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
        for c, lab in enumerate("xyz"):
            ax[0].plot(times, U_true[:, c, jz, jy, jx], "k-", lw=0.8)
            ax[0].plot(times, U_rec[:, c, jz, jy, jx], "--", label=f"$u_{lab}$")
        ax[0].set_xlabel("normalised time")
        ax[0].set_ylabel("displacement (voxels)")
        ax[0].set_title("top-centre node: imposed (solid), recovered (dashed)")
        ax[0].legend()
        ax[1].semilogy([h["rms_residual"] for h in history], "o-")
        ax[1].axhline(noise_std, color="k", ls=":", label="noise std")
        ax[1].set_xlabel("Gauss-Newton iteration")
        ax[1].set_ylabel("rms projection residual")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(args.out / "synthetic.png", dpi=150)
        print(f"wrote {args.out / 'synthetic.png'}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
