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
from pynamix.spectral.tomography import (
    Field,
    SpeciesTomography,
    _gradient_operator,
    mixture_fields,
    population_field,
    porod_field,
    species_field,
    patch_spectra,
    trilinear_weights,
)


def parse(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=128.0, help="cylinder radius (voxels)")
    ap.add_argument("--half-height", type=float, default=128.0)
    ap.add_argument("--n-proj", type=int, default=180)
    ap.add_argument("--packings", type=int, default=12, help="distinct packings (0: one per projection)")
    ap.add_argument("--d-small", type=float, default=5.0)
    ap.add_argument("--d-large", type=float, default=10.0)
    ap.add_argument("--phi", type=float, default=0.30, help="mean solid fraction")
    ap.add_argument("--cv", type=float, default=0.0, help="coefficient of variation of each population's diameter")
    ap.add_argument("--bins", type=int, default=4, help="diameter bins per population when --cv > 0")
    ap.add_argument(
        "--model", choices=["species", "porod"], default="species", help="species fractions, or phi and surface only"
    )
    ap.add_argument("--pixel", type=float, default=1.0, help="detector pixel pitch (voxels at the detector)")
    ap.add_argument("--noise", type=float, default=0.3, help="std of white noise in X per pixel (voxels)")
    ap.add_argument("--elem", type=int, default=2)
    ap.add_argument("--knots", type=int, default=4)
    ap.add_argument("--patch", type=int, default=48, help="spectral patch width (pixels)")
    ap.add_argument("--mean-patch", type=int, default=16, help="patch width for the mean rows (pixels)")
    ap.add_argument("--qd-min", type=float, default=8.0, help="lower band edge, as q d_small in the object")
    ap.add_argument("--reg", type=float, default=1e-3)
    ap.add_argument(
        "--no-structure", dest="structure", action="store_false", help="self term only, no grain interference"
    )
    ap.add_argument("--spectra-weight", type=float, default=1.0, help="0 fits the patch means alone")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/segregation"))
    ap.add_argument("--ref-elem", type=int, default=None, help="3D mesh of the realised reference (default: same)")
    ap.add_argument("--ref-knots", type=int, default=None, help="time knots of the realised reference")
    ap.add_argument("--eval-taus", type=int, default=11, help="times at which errors are evaluated")
    ap.add_argument("--layers", action="store_true", help="fields vary with height only (baseline)")
    ap.add_argument("--duration", type=float, default=None, help="fraction of the turn over which the change happens")
    ap.add_argument("--wall", action="store_true", help="pack against the container wall, leaving a wall layer")
    ap.add_argument("--wall-model", action="store_true", help="add a field in a shell at the wall per population")
    ap.add_argument("--eval-margin", type=float, default=0.0, help="evaluate errors this far inside the surface")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--cache", type=Path, default=None, help="directory for cached projections and samples")
    return ap.parse_args(argv)


class Fields:
    """Target solid fraction and large-grain fraction of the solid."""

    def __init__(self, R, H, phi, duration=None):
        self.R, self.H, self.phi0, self.duration = R, H, phi, duration

    def progress(self, tau):
        """Fraction of the change done by time tau: tau itself, or a smooth step lasting ``duration``."""
        if self.duration is None:
            return tau
        t = np.clip((tau - 0.5) / self.duration + 0.5, 0, 1)
        return t * t * (3 - 2 * t)

    def phi(self, x, tau):
        R, H = self.R, self.H
        column = np.exp(-((x[:, 0] - 0.4 * R) ** 2 + x[:, 1] ** 2) / (0.35 * R) ** 2)
        return self.phi0 - 0.04 * x[:, 2] / H - 0.06 * self.progress(tau) * column

    def large(self, x, tau):
        R, H = self.R, self.H
        r2 = (x[:, 0] ** 2 + x[:, 1] ** 2) / R**2
        c = 0.5 + self.progress(tau) * (0.35 * np.tanh(x[:, 2] / (0.4 * H)) + 0.2 * (r2 - 0.5))
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


def populations(args):
    """Bin diameters and volume fractions [(d_j, f_j)] of the small and large populations.

    With ``cv`` = 0 each population is one size. Otherwise diameters follow a
    Schulz (gamma) number distribution with the population's mean diameter
    and coefficient of variation ``cv``, split into ``bins`` bins of equal
    number; each bin is represented by the diameter with its mean volume.
    """
    out = []
    for mean in (args.d_small, args.d_large):
        if args.cv <= 0:
            out.append((np.array([mean]), np.array([1.0])))
            continue
        k = 1.0 / args.cv**2
        d = np.sort(np.random.default_rng(12345).gamma(k, mean / k, 200000))
        dj = np.array([np.mean(g**3) ** (1 / 3) for g in np.array_split(d, args.bins)])
        fj = dj**3 / np.sum(dj**3)
        out.append((dj, fj))
    return out


def population_of(d, pops):
    """Population index (0 small, 1 large) of grains with diameters d."""
    return np.where(np.isin(np.round(d, 6), np.round(pops[1][0], 6)), 1, 0)


def make_data(args, fields, geom, times, rng):
    """Projections (n_proj, nr, nc) and the packings [(centres, diameters, tau)] they were taken of."""
    R, H = args.radius, args.half_height
    pops = populations(args)
    diam = np.concatenate([p[0] for p in pops])
    order = np.argsort(-diam)  # place large grains first
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
            fp = fields.species(mc, tau)
            counts, weights = [], []
            for a, (dj, fj) in enumerate(pops):
                for d_, f_ in zip(dj, fj):
                    counts.append(int(fp[a].mean() * f_ * V / (np.pi * d_**3 / 6)))
                    weights.append(lambda x, t=tau, a=a, mx=fp[a].max(): fields.species(x, t)[a] / mx)
            # pack a larger cylinder and keep the grains centred inside the specimen, so the
            # packing has no wall layer (sequential addition crowds centres against a wall)
            m = 0.0 if args.wall else 2 * diam.max()
            grow = (R + m) ** 2 * (H + m) / (R**2 * H)
            c, d = sphere_packing(
                R + m,
                H + m,
                list(diam[order]),
                [int(counts[j] * grow) for j in order],
                [weights[j] for j in order],
                rng,
                contained=args.wall,
            )
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
    out = np.empty((K, len(weight_fns), nn))
    psi = basis(np.asarray(taus))
    for a in range(len(weight_fns)):
        Nw = Nq.multiply(weight_fns[a](q)[:, None]).tocsr()
        Ms = (Nw.T @ Nw).toarray() * dv
        Hm = np.kron(psi.T @ psi, Ms)
        b = sum(np.kron(psi[k], rhs_fn(a, k)) for k in range(len(taus)))
        # a light gradient penalty keeps nodes with little material in their support finite
        Hm += 1e-4 * np.mean(np.diag(Hm)) * np.kron(np.eye(K), L)
        out[:, a] = np.linalg.solve(Hm, b).reshape(K, nn)
    return out


def realised_fields(mesh, basis, packings, groups, R, H, margin):
    """Fields of the packings themselves, each grain counted at its centre.

    ``groups``: list of (taper, select, amount): ``select(d)`` picks the grains
    of a field by diameter and ``amount(d)`` is what each contributes (its
    volume for a volume fraction, its surface area for a specific surface).
    Exact up to the variation of the shape functions across a grain.
    """
    q, dv = grid(R, H, margin)

    def rhs(a, k):
        c, d, _ = packings[k]
        taper_, select, amount = groups[a]
        sel = select(d)
        N = trilinear_weights(c[sel], mesh.vol_shape, mesh.node_shape)
        return N.T @ (amount(d[sel]) * taper_(c[sel]))

    return l2_fields(mesh, basis, [p[2] for p in packings], rhs, [g[0] for g in groups], q, dv)


def target_fields(mesh, basis, fields, R, H, n_tau=25):
    """Best fit of the target fields inside the cylinder."""
    q, dv = grid(R, H, 0.0)
    q = q[np.hypot(q[:, 0], q[:, 1]) < R]
    taus = np.linspace(0, 1, n_tau)
    N = trilinear_weights(q, mesh.vol_shape, mesh.node_shape)
    one = [lambda x: np.ones(len(x))] * 2
    return l2_fields(mesh, basis, taus, lambda a, k: dv * (N.T @ fields.species(q, taus[k])[a]), one, q, dv)


def field_values(f, mesh, basis, R, H, h=8.0, n_tau=11):
    """Fields (n_knots, 2, n_nodes) evaluated at points filling the specimen, (n_tau, 2, n_points)."""
    q, _ = grid(R, H, 0.0, h)
    q = q[(np.hypot(q[:, 0], q[:, 1]) < R) & (np.abs(q[:, 2]) < H)]
    N = trilinear_weights(q, mesh.vol_shape, mesh.node_shape)
    psi = basis(np.linspace(0, 1, n_tau))
    return np.stack([np.stack([N @ (psi[t] @ f[:, a]) for a in range(2)]) for t in range(n_tau)])


def field_sigma(cov, n_fields, mesh, basis, R, H, h=8.0, n_tau=11):
    """Standard deviations (n_tau, 3, n_points) of field 0, field 1 and their sum at the evaluation points."""
    q, _ = grid(R, H, 0.0, h)
    q = q[(np.hypot(q[:, 0], q[:, 1]) < R) & (np.abs(q[:, 2]) < H)]
    N = trilinear_weights(q, mesh.vol_shape, mesh.node_shape).toarray()
    psi = basis(np.linspace(0, 1, n_tau))
    K, nn = basis.n_knots, mesh.n_nodes
    out = np.empty((n_tau, 3, len(q)))
    for t in range(n_tau):
        rows = []
        for sel in ([1, 0], [0, 1], [1, 1]):
            c = np.zeros((K, n_fields))
            c[:, :2] = psi[t][:, None] * np.array(sel)[None]
            rows.append(np.einsum("ka,pn->pkan", c, N).reshape(len(q), -1))
        for j, W in enumerate(rows):
            out[t, j] = np.sqrt(np.maximum(np.einsum("pi,ij,pj->p", W, cov, W), 0))
    return out


def mixture_s0(diameters, phi_mean):
    """Partial structure factors S_ab at q -> 0 (Percus-Yevick) at the mean composition, as a matrix."""
    from pynamix.spectral.closure import py_partials

    d = np.asarray(diameters, dtype=float)
    rho = np.asarray(phi_mean) / (np.pi * d**3 / 6)
    r = py_partials(
        {a: d[a] for a in range(d.size)}, {a: rho[a] for a in range(d.size)}, np.array([1e-3]), return_result=True
    )
    n = d.size
    return np.array([[float(r.S[(min(a, b), max(a, b))][0]) for b in range(n)] for a in range(n)])


def mass_diagonal(mesh, R, H, h=8.0):
    """N(x)^T M^-1 N(x) at the evaluation points, with M the mass matrix over the cylinder."""
    q, _ = grid(R, H, 0.0, h)
    q = q[(np.hypot(q[:, 0], q[:, 1]) < R) & (np.abs(q[:, 2]) < H)]
    qf, dv = grid(R, H, 0.0, 4.0)
    qf = qf[(np.hypot(qf[:, 0], qf[:, 1]) < R) & (np.abs(qf[:, 2]) < H)]
    Nf = trilinear_weights(qf, mesh.vol_shape, mesh.node_shape)
    M = (Nf.T @ Nf).toarray() * dv
    M += 1e-9 * np.mean(np.diag(M)) * np.eye(len(M))
    N = trilinear_weights(q, mesh.vol_shape, mesh.node_shape).toarray()
    return np.einsum("pi,ij,pj->p", N, np.linalg.inv(M), N)


def coverage(vp, vr, sig):
    """Fraction of evaluation points whose error is within 1 and 2 predicted standard deviations."""
    out = {}
    for j, name in enumerate(["phi_small", "phi_large", "phi"]):
        e = vp[:, 0] - vr[:, 0] if j == 0 else vp[:, 1] - vr[:, 1] if j == 1 else vp.sum(1) - vr.sum(1)
        z = np.abs(e) / np.maximum(sig[:, j], 1e-12)
        out[name] = dict(
            rms_error=float(np.sqrt(np.mean(e**2))),
            rms_sigma=float(np.sqrt(np.mean(sig[:, j] ** 2))),
            within_1sigma=float(np.mean(z < 1)),
            within_2sigma=float(np.mean(z < 2)),
        )
    return out


def errors(rec, ref, h):
    """rms differences of species fractions, solid fraction, large fraction and d32.

    ``rec`` and ``ref`` are field values (n_tau, 2, n_points) at points filling the specimen,
    so the errors do not depend on where the mesh nodes lie.
    """

    def rms(x):
        return float(np.sqrt(np.mean(x**2)))

    def d32(p):
        # 1/d32 = sum_a phi_a h_a / phi, with h_a = sum_j f_j / d_j over population a's bins
        return p.sum(1) / (p[:, 0] * h[0] + p[:, 1] * h[1])

    tot_r, tot_t = rec.sum(1), ref.sum(1)
    return dict(
        phi_small=rms(rec[:, 0] - ref[:, 0]),
        phi_large=rms(rec[:, 1] - ref[:, 1]),
        phi=rms(tot_r - tot_t),
        large_fraction=rms(rec[:, 1] / tot_r - ref[:, 1] / tot_t),
        d32=rms(d32(rec) - d32(ref)),
    )


def errors_porod(rec, ref):
    """rms differences of solid fraction and d32 = 6 phi / s from (n_tau, 2, n) values of (phi, s)."""

    def rms(x):
        return float(np.sqrt(np.mean(x**2)))

    return dict(phi=rms(rec[:, 0] - ref[:, 0]), d32=rms(6 * rec[:, 0] / rec[:, 1] - 6 * ref[:, 0] / ref[:, 1]))


def run(args):
    rng = np.random.default_rng(args.seed)
    R, H = args.radius, args.half_height
    ds, dl = args.d_small, args.d_large
    fields = Fields(R, H, args.phi, args.duration)
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
    if args.cv > 0:
        tag += f"_cv{args.cv}_{args.bins}"
    if args.wall:
        tag += "_wall"
    if args.duration is not None:
        tag += f"_dur{args.duration}"
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

    pops = populations(args)
    # volume-weighted mean radius of each population sets its surface taper
    radii = [np.sum(fj * dj) / 2 for dj, fj in pops]
    h = [float(np.sum(fj / dj)) for dj, fj in pops]  # 1/d32 of each population
    if args.wall:
        # centres confined to a grain radius inside the wall: material tapers from R - 2a to R
        tapers = [lambda x, a=a: taper(x, R - a, H - a, a) for a in radii]
    else:
        tapers = [lambda x, a=a: taper(x, R, H, a) for a in radii]
    dl_max = max(p[0].max() for p in pops)

    def shell(x):
        dist = np.minimum(R - np.hypot(x[:, 0], x[:, 1]), H - np.abs(x[:, 2]))
        return np.clip((1.5 * dl_max - dist) / dl_max, 0, 1) * (dist >= 0)

    extra = []
    if args.wall_model:
        extra = [population_field(dj, fj, weight=shell, name=f"wall_{a}") for a, (dj, fj) in enumerate(pops)]
    dmax = max(p[0].max() for p in pops)
    vol = (int(2 * H + dmax) + 1, int(2 * R + dmax) + 1, int(2 * R + dmax) + 1)
    mesh = RegularHexMesh(vol, (args.elem, 1, 1) if args.layers else args.elem)
    basis = LinearTimeBasis(args.knots)
    t0 = time.time()
    if args.model == "porod":
        # solid fraction from the means, specific surface from the Porod tail; taper at the mean radius
        tm = [lambda x: taper(x, R, H, float(np.mean(radii)))]
        flds = [Field("phi", lambda k: 0.0 * k, 1.0, tm[0]), porod_field(weight=tm[0])]
    elif args.cv > 0:
        # each population's kernel averages v F^2 over its bins (self term only)
        flds = [population_field(dj, fj, weight=tapers[a]) for a, (dj, fj) in enumerate(pops)]
    else:
        flds = [species_field(ds, tapers[0]), species_field(dl, tapers[1])]
    tomo = SpeciesTomography(geom, mesh, basis, times, None, spec, means=means, fields=flds + extra)
    phi = tomo.solve(reg=args.reg, spectra_weight=args.spectra_weight, verbose=False)
    if args.structure and args.model == "species":
        # second pass: interference terms at the mean composition of the first reconstruction;
        # every diameter bin is a species of the mixture, grouped back into the two populations
        mean_phi = np.clip(field_values(phi, mesh, basis, R, H).mean((0, 2))[:2], 1e-3, None)
        dbins = np.concatenate([dj for dj, _ in pops])
        phibins = np.concatenate([mean_phi[a] * fj for a, (_, fj) in enumerate(pops)])
        grp = [list(range(len(pops[0][0]))), list(range(len(pops[0][0]), len(dbins)))]
        mfields = mixture_fields(dbins, phibins, weights=tapers, groups=grp)
        tomo = SpeciesTomography(geom, mesh, basis, times, None, spec, means=means, fields=mfields + extra)
        phi = tomo.solve(reg=args.reg, spectra_weight=args.spectra_weight, verbose=False)
    t_solve = time.time() - t0
    n_fields = phi.shape[1]
    phi = phi[:, :2]  # bulk fields only; wall fields are not evaluated

    # errors are measured in the region at least eval_margin inside the surface
    m = args.eval_margin
    Re, He = R - m, H - m
    if m > 0:
        tapers_eval = [lambda x, a=a: taper(x, Re - a, He - a, a) for a in radii]

        def inner(c):
            return (np.hypot(c[:, 0], c[:, 1]) < Re - 1e-9) & (np.abs(c[:, 2]) < He)

        packings_eval = [(c[inner(c)], d[inner(c)], t) for c, d, t in packings]
    else:
        tapers_eval, packings_eval = tapers, packings

    def volume(d):
        return np.pi * d**3 / 6

    groups = [(tapers_eval[a], lambda d, a=a: population_of(d, pops) == a, volume) for a in range(2)]
    # the reference may use its own mesh and time basis, so different reconstructions share one reference
    rmesh = RegularHexMesh(vol, args.ref_elem) if args.ref_elem else mesh
    rbasis = LinearTimeBasis(args.ref_knots) if args.ref_knots else basis
    realised = realised_fields(rmesh, rbasis, packings_eval, groups, Re, He, dmax / 2)
    target = target_fields(rmesh, rbasis, fields, Re, He)

    Xn = mesh.node_coords().reshape(-1, 3)
    node_in = (np.hypot(Xn[:, 0], Xn[:, 1]) <= R + dmax / 2 + 1e-6) & (np.abs(Xn[:, 2]) <= H + dmax / 2 + 1e-6)
    # with a margin, evaluate a large diameter inside the evaluation surface, where its taper is one
    Rv, Hv = (Re - dmax, He - dmax) if m > 0 else (R, H)
    nt = args.eval_taus
    vr, vt = (field_values(f, rmesh, rbasis, Rv, Hv, n_tau=nt) for f in (realised, target))
    vp = field_values(phi, mesh, basis, Rv, Hv, n_tau=nt)
    # grains of each population per element, at the mean packing, for an element wholly inside the specimen
    elem_vol = float(np.prod([(n - 1) / e for n, e in zip(mesh.vol_shape, mesh.n_elem)]))
    spec_vol = np.pi * R**2 * 2 * H
    counts = np.mean([np.bincount(population_of(d, pops), minlength=2) for _, d, _ in packings], 0)
    counts = counts * elem_vol / spec_vol
    result = dict(
        args={k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        grains_per_element=dict(small=float(counts[0]), large=float(counts[1]), total=float(counts.sum())),
        n_patches=int(spec.centres.shape[0]),
        n_mean_patches=int(means.centres.shape[0]),
        n_annuli=int(spec.q.size),
        qd_small=[float(spec.q.min() * M0 * ds), float(spec.q.max() * M0 * ds)],
        chi2=dict(means=float(tomo.history[-1]["chi2_means"]), spectra=tomo.history[-1]["chi2_spectra"].tolist()),
        error_vs_realised=errors(vp, vr, h) if args.model == "species" else None,
        error_vs_target=errors(vp, vt, h) if args.model == "species" else None,
        realised_vs_target=errors(vr, vt, h),
        realised_range=dict(
            phi=[float(vr.sum(1).min()), float(vr.sum(1).max())],
            large_fraction=[
                float((vr[:, 1] / vr.sum(1)).min()),
                float((vr[:, 1] / vr.sum(1)).max()),
            ],
        ),
        seconds=dict(data=t_data, solve=t_solve),
    )
    if args.model == "species":
        sig = field_sigma(tomo.covariance, n_fields, mesh, basis, Rv, Hv, n_tau=nt)
        result["uncertainty"] = coverage(vp, vr, sig)
        # add the grain-count floor of the realised field: S_aa(0) phi_a v_a N^T M^-1 N
        S0 = mixture_s0([ds, dl], np.clip(vp.mean((0, 2)), 1e-3, None))
        gx = mass_diagonal(rmesh, Rv, Hv)
        pv = np.sqrt(np.maximum(vp, 0) * (np.pi * np.array([ds, dl]) ** 3 / 6)[None, :, None])  # sqrt(phi_a v_a)
        var = np.einsum("ab,tap,tbp->tp", S0, pv, pv) * gx  # total phi, with cross-correlation
        fl = np.stack([np.sqrt(S0[0, 0] * gx) * pv[:, 0], np.sqrt(S0[1, 1] * gx) * pv[:, 1], np.sqrt(var)], 1)
        result["uncertainty_with_floor"] = coverage(vp, vr, np.sqrt(sig**2 + fl**2))
    if args.model == "porod":
        # realised specific surface: each grain counts its surface area 6 v / d
        am = float(np.mean(radii))
        tm = flds[0].weight if m == 0 else (lambda x: taper(x, Re - am, He - am, am))
        sgroups = [
            (tm, lambda d: np.ones(d.size, bool), volume),
            (tm, lambda d: np.ones(d.size, bool), lambda d: 6 * volume(d) / d),
        ]
        rs = realised_fields(mesh, basis, packings_eval, sgroups, Re, He, dmax / 2)
        rs = field_values(rs, mesh, basis, Rv, Hv, n_tau=nt)
        result["error_vs_realised"] = errors_porod(vp, rs)
        result["realised_d32_range"] = [float((6 * rs[:, 0] / rs[:, 1]).min()), float((6 * rs[:, 0] / rs[:, 1]).max())]
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez(args.out / "fields.npz", phi=phi, realised=realised, target=target, node_in=node_in)
    (args.out / "result.json").write_text(json.dumps(result, indent=1))
    if not args.no_plot and args.model == "species":
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
