"""Species volume fractions of a bidisperse packing from a rotating projection sequence.

A cylinder of spheres of two sizes is imaged while the gantry turns once. The
large grains sit on top, and the front between large-rich and small-rich
material descends during the sequence; the front is curved and tilted, so
the fields are not axisymmetric. Each projection sees a packing drawn from
the target fields at its own time (``--packings`` sets how many distinct
packings are used across the sequence; fewer means the same realisation is
seen from several angles).

The projections are tiled into patches, the patch spectra are taken, and
phi_small(x, tau) and phi_large(x, tau) are reconstructed on a coarse mesh
by :class:`pynamix.spectral.tomography.SpeciesTomography`.
"""

import argparse
import time
from pathlib import Path

import numpy as np

from pynamix.pdvc import ConeBeamGeometry, LinearTimeBasis, RegularHexMesh
from pynamix.pdvc.phantom import project_spheres, sphere_packing
from pynamix.spectral.tomography import SpeciesTomography, patch_spectra, trilinear_weights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=160)
    ap.add_argument("--n-proj", type=int, default=180)
    ap.add_argument("--packings", type=int, default=0, help="distinct packings (0: one per projection)")
    ap.add_argument("--d-small", type=float, default=6.0)
    ap.add_argument("--d-large", type=float, default=12.0)
    ap.add_argument("--phi", type=float, default=0.3)
    ap.add_argument("--elem", type=int, default=4)
    ap.add_argument("--knots", type=int, default=4)
    ap.add_argument("--patch", type=int, default=48)
    ap.add_argument("--stride", type=int, default=24)
    ap.add_argument("--qd-min", type=float, default=10.0, help="lower band edge, as q d_small in the object")
    ap.add_argument("--noise", type=float, default=0.3, help="std of white noise in X (voxels)")
    ap.add_argument("--reg", type=float, default=1e-3)
    ap.add_argument("--spectra-weight", type=float, default=1.0, help="0 fits the patch means alone")
    ap.add_argument("--mean-error", type=float, default=0.01, help="relative error assigned to patch means")
    ap.add_argument("--static", action="store_true", help="fields constant in time")
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/species"))
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    N = args.n
    R = H = 0.4 * N  # cylinder radius and half-height
    ds, dl = args.d_small, args.d_large
    vs, vl = np.pi * ds**3 / 6, np.pi * dl**3 / 6
    angles = np.linspace(0, 2 * np.pi, args.n_proj, endpoint=False)
    times = np.linspace(0, 1, args.n_proj)
    sod, sdd = 4 * N, 6 * N
    M0 = sdd / sod
    nd = int(np.ceil(2 * H * sdd / (sod - R) / 16) * 16)
    geom = ConeBeamGeometry.circular(angles, sod, sdd, (nd, nd), 1.0)

    def large_fraction(x, tau):
        """Fraction of the solid volume in large grains."""
        r2 = (x[:, 0] ** 2 + x[:, 1] ** 2) / R**2
        z0 = (0.0 if args.static else 0.5 * H * (1 - 2 * tau)) + 0.3 * H * r2 + 0.2 * H * x[:, 0] / R
        return np.clip(0.5 + 0.4 * np.tanh((x[:, 2] - z0) / (0.2 * H)), 0.0, 1.0)

    def inside(x):
        return ((x[:, 0] ** 2 + x[:, 1] ** 2) < R**2) & (np.abs(x[:, 2]) < H)

    # projections
    n_pack = args.packings or args.n_proj
    which = np.minimum((times * n_pack).astype(int), n_pack - 1)
    V = np.pi * R**2 * 2 * H
    mc = rng.random((20000, 3)) * [2 * R, 2 * R, 2 * H] - [R, R, H]
    mc = mc[inside(mc)]
    t0 = time.time()
    X = np.empty((args.n_proj,) + geom.det_shape)
    cache = {}
    samples = []  # (points, tau, small indicator, large indicator) per packing, for the realised fields
    for i in range(args.n_proj):
        k = which[i]
        if k not in cache:
            tau = (k + 0.5) / n_pack
            cl = large_fraction(mc, tau).mean()
            counts = [int(args.phi * cl * V / vl), int(args.phi * (1 - cl) * V / vs)]
            w = [lambda x, t=tau: large_fraction(x, t), lambda x, t=tau: 1 - large_fraction(x, t)]
            cache = {k: sphere_packing(R, H, [dl, ds], counts, w, rng)}
            samples.append(realised_samples(*cache[k], tau, R, H, rng, n=max(3000, 300000 // n_pack)))
        c, d = cache[k]
        X[i] = project_spheres(c, d, geom, i)
    X += args.noise * rng.standard_normal(X.shape)
    print(
        f"{args.n_proj} projections of {n_pack} packings ({len(d)} grains in the last), "
        f"detector {geom.det_shape}, {time.time() - t0:.0f} s"
    )

    # patch spectra over a band where the interference terms are small
    q_min = args.qd_min / (ds * M0)
    spec = patch_spectra(
        X, size=args.patch, stride=args.stride, q_band=(q_min, 0.6 * np.pi), noise_var=args.noise**2, pixel_mtf=2
    )
    keep = spec.mean.mean(0) > 0.25 * spec.mean.mean(0).max()
    spec = spec.select(keep)
    print(
        f"{keep.sum()} patches of {args.patch} px, {spec.q.size} annuli, "
        f"q d_small (object, axis) {spec.q.min() * M0 * ds:.1f}-{spec.q.max() * M0 * ds:.1f}"
    )

    mesh = RegularHexMesh((N, N, N), args.elem)
    basis = LinearTimeBasis(1 if args.static else args.knots)
    tomo = SpeciesTomography(geom, mesh, basis, times, [ds, dl], spec, inside=inside)
    t0 = time.time()
    phi = tomo.solve(reg=args.reg, mean_rel_error=args.mean_error, spectra_weight=args.spectra_weight)
    print(f"solved in {time.time() - t0:.0f} s")

    # references: the target fields and the fields the packings realised, each
    # as its best fit in the same space-time basis
    pts = rng.random((60000, 3)) * [2 * R, 2 * R, 2 * H] - [R, R, H]
    pts = pts[inside(pts)]
    tau = rng.random(pts.shape[0])
    cl = large_fraction(pts, tau)
    target = fit_basis(mesh, basis, [(pts, tau, args.phi * (1 - cl), args.phi * cl)])
    realised = fit_basis(mesh, basis, samples)

    Xn = mesh.node_coords().reshape(-1, 3)
    node_in = inside(Xn)
    print(f"nodes inside specimen: {node_in.sum()} x {basis.n_knots} knots")
    for name, ref in [("target", target), ("realised", realised)]:
        report(name, phi[:, :, node_in], ref[:, :, node_in], ds, dl)
    report("realised vs target", realised[:, :, node_in], target[:, :, node_in], ds, dl)

    args.out.mkdir(parents=True, exist_ok=True)
    np.savez(args.out / "result.npz", phi=phi, target=target, realised=realised, node_in=node_in)
    plot(args, mesh, basis, phi, realised, inside, R, H)


def realised_samples(centres, diam, tau, R, H, rng, n=3000):
    """Species indicators at random points of the cylinder, for one packing."""
    from scipy.spatial import cKDTree

    p = rng.random((3 * n, 3)) * [2 * R, 2 * R, 2 * H] - [R, R, H]
    p = p[(p[:, 0] ** 2 + p[:, 1] ** 2 < R**2)][:n]
    # the nearest centre need not be the grain containing the point, so test several
    dist, j = cKDTree(centres).query(p, k=8)
    hit = dist < diam[j] / 2
    small = (hit & (diam[j] == diam.min())).any(1)
    large = (hit & (diam[j] == diam.max())).any(1)
    return p, np.full(len(p), tau), small.astype(float), large.astype(float)


def fit_basis(mesh, basis, samples):
    """Least-squares nodal fields (n_knots, 2, n_nodes) from sampled values."""
    nu = basis.n_knots * mesh.n_nodes
    H = np.zeros((nu, nu))
    b = np.zeros((2, nu))
    for p, tau, fs, fl in samples:
        Nm = trilinear_weights(p, mesh.vol_shape, mesh.node_shape).toarray()
        A = (basis(tau)[:, :, None] * Nm[:, None, :]).reshape(len(p), -1)
        H += A.T @ A
        b += np.stack([A.T @ fs, A.T @ fl])
    H += 1e-6 * np.mean(np.diag(H)) * np.eye(nu)
    x = np.linalg.solve(H, b.T).T
    return x.reshape(2, basis.n_knots, -1).transpose(1, 0, 2)


def report(name, rec, ref, ds, dl):
    """rms differences of the species fractions, solid fraction, large fraction and d32."""

    def rms(x):
        return np.sqrt(np.mean(x**2))

    tot_r, tot_t = rec.sum(1), ref.sum(1)

    def d32(p):
        return p.sum(1) / (p[:, 0] / ds + p[:, 1] / dl)

    print(
        f"{name:20s} rms: phi_small {rms(rec[:, 0] - ref[:, 0]):.3f}, phi_large {rms(rec[:, 1] - ref[:, 1]):.3f}, "
        f"phi {rms(tot_r - tot_t):.3f}, large fraction {rms(rec[:, 1] / tot_r - ref[:, 1] / tot_t):.3f}, "
        f"d32 {rms(d32(rec) - d32(ref)):.2f} vox"
    )


def plot(args, mesh, basis, phi, true, inside, R, H):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = np.linspace(-R, R, 81)
    zs = np.linspace(-H, H, 81)
    Z, Xg = np.meshgrid(zs, xs, indexing="ij")
    pts = np.stack([Xg.ravel(), np.zeros(Xg.size), Z.ravel()], 1)
    Nm = trilinear_weights(pts, mesh.vol_shape, mesh.node_shape)
    taus = [0.0, 0.5, 1.0] if basis.n_knots > 1 else [0.0]
    fig, ax = plt.subplots(2, len(taus), figsize=(3.2 * len(taus), 6), squeeze=False)
    for j, t in enumerate(taus):
        psi = basis(np.array([t]))[0]
        for r, (name, f) in enumerate([("realised", true), ("reconstructed", phi)]):
            p = np.einsum("k,kan->an", psi, f)
            c = (Nm @ p[1]) / (Nm @ p.sum(0))
            im = ax[r, j].imshow(
                c.reshape(Z.shape), origin="lower", extent=[-R, R, -H, H], vmin=0, vmax=1, cmap="viridis"
            )
            ax[r, j].set_title(f"{name}, τ = {t:.1f}", fontsize=9)
            ax[r, j].set_xlabel("x (vox)")
    ax[0, 0].set_ylabel("z (vox)")
    ax[1, 0].set_ylabel("z (vox)")
    fig.colorbar(im, ax=ax, shrink=0.6, label="large-grain fraction of solid")
    fig.savefig(args.out / "species.png", dpi=130, bbox_inches="tight")


if __name__ == "__main__":
    main()
