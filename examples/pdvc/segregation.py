"""Evolving segregation in a bidisperse packing with spatially varying solid fraction.

A cylinder of spheres of two sizes is imaged while the gantry turns once.
During the turn:

- the solid fraction decreases with height and a dilated column grows off
  the axis, so the porosity varies in 3D and in time;
- the initially mixed packing segregates, with large grains moving up and
  outwards.

Each projection sees a packing drawn from the target fields at its own time
(``--packings`` sets how many distinct packings are used across the
sequence). The fields phi_small(x, tau) and phi_large(x, tau) are
reconstructed from detector-patch means and spectra by
:class:`pynamix.spectral.tomography.SpeciesTomography` and compared with the
fields the packings realised.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from pynamix.pdvc import ConeBeamGeometry, LinearTimeBasis, RegularHexMesh
from pynamix.pdvc.phantom import project_spheres, sphere_packing
from pynamix.spectral.tomography import SpeciesTomography, _gradient_operator, patch_spectra, trilinear_weights


def parse(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=128.0, help="cylinder radius (voxels)")
    ap.add_argument("--half-height", type=float, default=128.0)
    ap.add_argument("--n-proj", type=int, default=180)
    ap.add_argument("--packings", type=int, default=12, help="distinct packings (0: one per projection)")
    ap.add_argument("--d-small", type=float, default=5.0)
    ap.add_argument("--d-large", type=float, default=10.0)
    ap.add_argument("--phi", type=float, default=0.30, help="mean solid fraction")
    ap.add_argument("--pixel", type=float, default=1.0, help="detector pixel pitch (voxels at the detector)")
    ap.add_argument("--noise", type=float, default=0.3, help="std of white noise in X per pixel (voxels)")
    ap.add_argument("--elem", type=int, default=2)
    ap.add_argument("--knots", type=int, default=4)
    ap.add_argument("--patch", type=int, default=48, help="spectral patch width (pixels)")
    ap.add_argument("--mean-patch", type=int, default=16, help="patch width for the mean rows (pixels)")
    ap.add_argument("--qd-min", type=float, default=8.0, help="lower band edge, as q d_small in the object")
    ap.add_argument("--reg", type=float, default=1e-3)
    ap.add_argument("--spectra-weight", type=float, default=1.0, help="0 fits the patch means alone")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/segregation"))
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--cache", type=Path, default=None, help="directory for cached projections and samples")
    return ap.parse_args(argv)


class Fields:
    """Target solid fraction and large-grain fraction of the solid."""

    def __init__(self, R, H, phi):
        self.R, self.H, self.phi0 = R, H, phi

    def phi(self, x, tau):
        R, H = self.R, self.H
        column = np.exp(-((x[:, 0] - 0.4 * R) ** 2 + x[:, 1] ** 2) / (0.35 * R) ** 2)
        return self.phi0 - 0.04 * x[:, 2] / H - 0.06 * tau * column

    def large(self, x, tau):
        R, H = self.R, self.H
        r2 = (x[:, 0] ** 2 + x[:, 1] ** 2) / R**2
        c = 0.5 + tau * (0.35 * np.tanh(x[:, 2] / (0.4 * H)) + 0.2 * (r2 - 0.5))
        return np.clip(c, 0.02, 0.98)

    def species(self, x, tau):
        p, c = self.phi(x, tau), self.large(x, tau)
        return p * (1 - c), p * c


def taper(x, R, H, a):
    """Expected fraction of bulk material at x for grains of radius a whose centres fill the cylinder.

    Each face is treated as a plane: a ball of radius a centred at signed
    distance s inside a face keeps the fraction 1/2 + 3t/4 - t^3/4, t = s/a.
    """

    def g(s):
        t = np.clip(s / a, -1, 1)
        return 0.5 + 0.75 * t - 0.25 * t**3

    return g(R - np.hypot(x[:, 0], x[:, 1])) * g(H - x[:, 2]) * g(H + x[:, 2])


def make_data(args, fields, geom, times, rng):
    """Projections (n_proj, nr, nc) and the packings [(centres, diameters, tau)] they were taken of."""
    R, H = args.radius, args.half_height
    ds, dl = args.d_small, args.d_large
    vs, vl = np.pi * ds**3 / 6, np.pi * dl**3 / 6
    V = np.pi * R**2 * 2 * H
    n_pack = args.packings or args.n_proj
    which = np.minimum((times * n_pack).astype(int), n_pack - 1)
    mc = rng.random((40000, 3)) * [2 * R, 2 * R, 2 * H] - [R, R, H]
    mc = mc[mc[:, 0] ** 2 + mc[:, 1] ** 2 < R**2]
    X = np.empty((args.n_proj,) + geom.det_shape)
    packings, cache = [], {}
    for i in range(args.n_proj):
        k = which[i]
        if k not in cache:
            tau = (k + 0.5) / n_pack
            fs, fl = fields.species(mc, tau)
            counts = [int(fl.mean() * V / vl), int(fs.mean() * V / vs)]
            wl = lambda x, t=tau: fields.species(x, t)[1] / fl.max()  # noqa: E731
            ws = lambda x, t=tau: fields.species(x, t)[0] / fs.max()  # noqa: E731
            # pack a larger cylinder and keep the grains centred inside the specimen, so the
            # packing has no wall layer (sequential addition crowds centres against a wall)
            m = 2 * dl
            grow = (R + m) ** 2 * (H + m) / (R**2 * H)
            c, d = sphere_packing(R + m, H + m, [dl, ds], [int(n * grow) for n in counts], [wl, ws], rng)
            keep = (np.hypot(c[:, 0], c[:, 1]) < R) & (np.abs(c[:, 2]) < H)
            cache = {k: (c[keep], d[keep])}
            packings.append((*cache[k], tau))
        X[i] = project_spheres(*cache[k], geom, i)
    return X, packings


def grid(R, H, margin, h=3.0):
    """Regular quadrature points (n, 3) with spacing h covering the cylinder plus a margin."""
    ax = np.arange(-R - margin, R + margin + h / 2, h)
    az = np.arange(-H - margin, H + margin + h / 2, h)
    Z, Y, X = np.meshgrid(az, ax, ax, indexing="ij")
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1), h**3


def l2_fields(mesh, basis, taus, rhs_fn, weight_fns, q, dv):
    """Nodal fields (n_knots, 2, n_nodes) minimising sum_k int (w_a f_a - g_a(tau_k))^2 over x.

    ``rhs_fn(a, k)`` returns int w_a N_n g_a dx (n_nodes,) at time taus[k];
    ``weight_fns[a]`` is w_a, the taper the field is multiplied by.
    """
    Nq = trilinear_weights(q, mesh.vol_shape, mesh.node_shape)
    nn, K = mesh.n_nodes, basis.n_knots
    D = _gradient_operator(mesh.node_shape)
    L = (D.T @ D).toarray() + 1e-6 * np.eye(nn)
    out = np.empty((K, 2, nn))
    psi = basis(np.asarray(taus))
    for a in range(2):
        Nw = Nq.multiply(weight_fns[a](q)[:, None]).tocsr()
        Ms = (Nw.T @ Nw).toarray() * dv
        Hm = np.kron(psi.T @ psi, Ms)
        b = sum(np.kron(psi[k], rhs_fn(a, k)) for k in range(len(taus)))
        # a light gradient penalty keeps nodes with little material in their support finite
        Hm += 1e-4 * np.mean(np.diag(Hm)) * np.kron(np.eye(K), L)
        out[:, a] = np.linalg.solve(Hm, b).reshape(K, nn)
    return out


def realised_fields(mesh, basis, packings, tapers, R, H, margin):
    """Fields of the packings themselves, each grain counted as its volume at its centre.

    Exact up to the variation of the shape functions across a grain.
    """
    q, dv = grid(R, H, margin)
    diams = np.unique(np.concatenate([p[1] for p in packings]))  # small, large

    def rhs(a, k):
        c, d, _ = packings[k]
        sel = d == diams[a]
        v = np.pi * diams[a] ** 3 / 6
        N = trilinear_weights(c[sel], mesh.vol_shape, mesh.node_shape)
        return v * (N.T @ tapers[a](c[sel]))

    return l2_fields(mesh, basis, [p[2] for p in packings], rhs, tapers, q, dv)


def target_fields(mesh, basis, fields, R, H, n_tau=25):
    """Best fit of the target fields inside the cylinder."""
    q, dv = grid(R, H, 0.0)
    q = q[np.hypot(q[:, 0], q[:, 1]) < R]
    taus = np.linspace(0, 1, n_tau)
    N = trilinear_weights(q, mesh.vol_shape, mesh.node_shape)
    one = [lambda x: np.ones(len(x))] * 2
    return l2_fields(mesh, basis, taus, lambda a, k: dv * (N.T @ fields.species(q, taus[k])[a]), one, q, dv)


def errors(rec, ref, ds, dl):
    """rms differences of species fractions, solid fraction, large fraction and d32 at the given nodes."""

    def rms(x):
        return float(np.sqrt(np.mean(x**2)))

    def d32(p):
        return p.sum(1) / (p[:, 0] / ds + p[:, 1] / dl)

    tot_r, tot_t = rec.sum(1), ref.sum(1)
    return dict(
        phi_small=rms(rec[:, 0] - ref[:, 0]),
        phi_large=rms(rec[:, 1] - ref[:, 1]),
        phi=rms(tot_r - tot_t),
        large_fraction=rms(rec[:, 1] / tot_r - ref[:, 1] / tot_t),
        d32=rms(d32(rec) - d32(ref)),
    )


def run(args):
    rng = np.random.default_rng(args.seed)
    R, H = args.radius, args.half_height
    ds, dl = args.d_small, args.d_large
    fields = Fields(R, H, args.phi)
    angles = np.linspace(0, 2 * np.pi, args.n_proj, endpoint=False)
    times = np.linspace(0, 1, args.n_proj)
    sod, sdd = 10 * R, 15 * R
    M0 = sdd / sod
    span = 2 * max(H, R) * sdd / (sod - R) / args.pixel
    nd = int(np.ceil(span / 16) * 16)
    geom = ConeBeamGeometry.circular(angles, sod, sdd, (nd, nd), args.pixel)

    t0 = time.time()
    keys = ["radius", "half_height", "n_proj", "packings", "d_small", "d_large", "phi", "pixel", "seed"]
    tag = "_".join(f"{getattr(args, k)}" for k in keys)
    cached = args.cache / f"data_{tag}.npz" if args.cache else None
    if cached is not None and cached.exists():
        z = np.load(cached)
        X0 = z["X"]
        packings = [(z[f"c{j}"], z[f"d{j}"], float(z[f"t{j}"])) for j in range(int(z["n"]))]
    else:
        X0, packings = make_data(args, fields, geom, times, rng)
        if cached is not None:
            cached.parent.mkdir(parents=True, exist_ok=True)
            flat = {f"{c}{j}": v for j, pk in enumerate(packings) for c, v in zip("cdt", pk)}
            np.savez(cached, X=X0, n=len(packings), **flat)
    # noise is added here so that cached noise-free projections serve every noise level
    X = X0 + args.noise * np.random.default_rng(args.seed + 1).standard_normal(X0.shape)
    t_data = time.time() - t0

    q_min = args.qd_min / (ds * M0)
    spec = patch_spectra(
        X,
        pitch=args.pixel,
        size=args.patch,
        stride=args.patch // 2,
        q_band=(q_min, 0.6 * np.pi / args.pixel),
        noise_var=args.noise**2,
        pixel_mtf=2,
    )
    spec = spec.select(spec.mean.mean(0) > 0.25 * spec.mean.max())
    means = patch_spectra(X, pitch=args.pixel, size=args.mean_patch, q_band=(0, np.inf), pixel_mtf=None)
    means = means.select(means.mean.mean(0) > 0.1 * means.mean.max())

    tapers = [lambda x: taper(x, R, H, ds / 2), lambda x: taper(x, R, H, dl / 2)]

    def inside(x):
        return np.stack([t(x) for t in tapers])

    vol = (int(2 * H + dl) + 1, int(2 * R + dl) + 1, int(2 * R + dl) + 1)
    mesh = RegularHexMesh(vol, args.elem)
    basis = LinearTimeBasis(args.knots)
    tomo = SpeciesTomography(geom, mesh, basis, times, [ds, dl], spec, means=means, inside=inside)
    t0 = time.time()
    phi = tomo.solve(reg=args.reg, spectra_weight=args.spectra_weight, verbose=False)
    t_solve = time.time() - t0

    realised = realised_fields(mesh, basis, packings, tapers, R, H, dl / 2)
    target = target_fields(mesh, basis, fields, R, H)

    Xn = mesh.node_coords().reshape(-1, 3)
    node_in = (np.hypot(Xn[:, 0], Xn[:, 1]) <= R + dl / 2 + 1e-6) & (np.abs(Xn[:, 2]) <= H + dl / 2 + 1e-6)
    result = dict(
        args={k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        n_patches=int(spec.centres.shape[0]),
        n_mean_patches=int(means.centres.shape[0]),
        n_annuli=int(spec.q.size),
        qd_small=[float(spec.q.min() * M0 * ds), float(spec.q.max() * M0 * ds)],
        chi2=dict(means=float(tomo.history[-1]["chi2_means"]), spectra=tomo.history[-1]["chi2_spectra"].tolist()),
        error_vs_realised=errors(phi[:, :, node_in], realised[:, :, node_in], ds, dl),
        error_vs_target=errors(phi[:, :, node_in], target[:, :, node_in], ds, dl),
        realised_vs_target=errors(realised[:, :, node_in], target[:, :, node_in], ds, dl),
        realised_range=dict(
            phi=[float(realised[:, :, node_in].sum(1).min()), float(realised[:, :, node_in].sum(1).max())],
            large_fraction=[
                float((realised[:, 1, node_in] / realised[:, :, node_in].sum(1)).min()),
                float((realised[:, 1, node_in] / realised[:, :, node_in].sum(1)).max()),
            ],
        ),
        seconds=dict(data=t_data, solve=t_solve),
    )
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez(args.out / "fields.npz", phi=phi, realised=realised, target=target, node_in=node_in)
    (args.out / "result.json").write_text(json.dumps(result, indent=1))
    if not args.no_plot:
        plot(args, mesh, basis, {"realised": realised, "reconstructed": phi}, R, H)
    return result


def plot(args, mesh, basis, fields, R, H):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = np.linspace(-R, R, 81)
    zs = np.linspace(-H, H, 81)
    Z, Xg = np.meshgrid(zs, xs, indexing="ij")
    pts = np.stack([Xg.ravel(), np.zeros(Xg.size), Z.ravel()], 1)
    Nm = trilinear_weights(pts, mesh.vol_shape, mesh.node_shape)
    taus = [0.0, 0.5, 1.0]
    fig, ax = plt.subplots(4, 3, figsize=(8.5, 11), constrained_layout=True)
    for j, t in enumerate(taus):
        psi = basis(np.array([t]))[0]
        for r, (name, f) in enumerate(fields.items()):
            p = np.einsum("k,kan->an", psi, f)
            tot = Nm @ p.sum(0)
            c = (Nm @ p[1]) / tot
            a1 = ax[r, j].imshow(
                tot.reshape(Z.shape), origin="lower", extent=[-R, R, -H, H], vmin=0.15, vmax=0.4, cmap="cividis"
            )
            ax[r, j].set_title(f"φ {name}, τ = {t:.1f}", fontsize=9)
            a2 = ax[r + 2, j].imshow(
                c.reshape(Z.shape), origin="lower", extent=[-R, R, -H, H], vmin=0, vmax=1, cmap="viridis"
            )
            ax[r + 2, j].set_title(f"large fraction {name}, τ = {t:.1f}", fontsize=9)
    for a in ax.ravel():
        a.set_xticks([])
        a.set_yticks([])
    fig.colorbar(a1, ax=ax[:2], shrink=0.6, label="solid fraction")
    fig.colorbar(a2, ax=ax[2:], shrink=0.6, label="large-grain fraction of solid")
    fig.savefig(args.out / "segregation.png", dpi=120)


def main():
    args = parse()
    res = run(args)
    print(json.dumps({k: res[k] for k in res if k != "args"}, indent=1))


if __name__ == "__main__":
    main()
