"""Measured patch spectra against the forward model, from the realised fields of a synthetic run.

Uses the cached radiographs of the default segregation case. The model is
evaluated with the self term only and with Percus-Yevick interference at the
mean composition, from the fields the packings realised (not from a
reconstruction), so the comparison tests the forward model alone.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import segregation as S  # noqa: E402

from pynamix.pdvc import ConeBeamGeometry, LinearTimeBasis, RegularHexMesh  # noqa: E402
from pynamix.spectral.tomography import (  # noqa: E402
    SpeciesTomography,
    mixture_fields,
    patch_spectra,
    species_field,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=Path("output/pdvc/segregation_study/cache"))
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/spectra_check.png"))
    ap.add_argument("--qd-range", type=float, nargs=2, default=[3.0, 14.0])
    a = ap.parse_args()
    args = S.parse(["--cache", str(a.cache)])
    R, H, ds, dl = args.radius, args.half_height, args.d_small, args.d_large
    keys = ["radius", "half_height", "n_proj", "packings", "d_small", "d_large", "phi", "pixel", "seed"]
    z = np.load(a.cache / ("data_" + "_".join(f"{getattr(args, k)}" for k in keys) + ".npz"))
    X = z["X"] + args.noise * np.random.default_rng(args.seed + 1).standard_normal(z["X"].shape)
    packings = [(z[f"c{j}"], z[f"d{j}"], float(z[f"t{j}"])) for j in range(int(z["n"]))]
    n = args.n_proj
    times = np.linspace(0, 1, n)
    sod, sdd = 10 * R, 15 * R
    M0 = sdd / sod
    geom = ConeBeamGeometry.circular(np.linspace(0, 2 * np.pi, n, endpoint=False), sod, sdd, X.shape[1:], 1.0)
    band = (a.qd_range[0] / (ds * M0), a.qd_range[1] / (ds * M0))
    spec = patch_spectra(X, size=args.patch, stride=args.patch // 2, q_band=band, noise_var=args.noise**2, pixel_mtf=2)
    spec = spec.select(spec.mean.mean(0) > 0.25 * spec.mean.max())
    tapers = [lambda x: S.taper(x, R, H, ds / 2), lambda x: S.taper(x, R, H, dl / 2)]
    vol = (int(2 * H + dl) + 1,) * 3
    mesh, basis = RegularHexMesh(vol, args.elem), LinearTimeBasis(args.knots)
    pops = S.populations(args)
    groups = [(tapers[k], lambda d, k=k: S.population_of(d, pops) == k, lambda d: np.pi * d**3 / 6) for k in range(2)]
    real = S.realised_fields(mesh, basis, packings, groups, R, H, dl / 2)
    mean_phi = S.field_values(real, mesh, basis, R, H).mean((0, 2))
    models = {
        "self term": [species_field(ds, tapers[0]), species_field(dl, tapers[1])],
        "with interference": mixture_fields([ds, dl], mean_phi, weights=tapers),
    }
    qd = spec.q * M0 * ds
    ratios = {}
    for name, flds in models.items():
        t = SpeciesTomography(geom, mesh, basis, times, None, spec, fields=flds)
        t.build(verbose=False)
        _, P = t.predict(real)
        ratios[name] = spec.P / P
        print(name, "median ratio per annulus:", np.round(np.median(ratios[name], (0, 1)), 3))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(4.6, 3.2), constrained_layout=True)
    for (name, r), c in zip(ratios.items(), ["C1", "C0"]):
        med = np.median(r, (0, 1))
        lo, hi = np.percentile(r, [25, 75], axis=(0, 1))
        ax.fill_between(qd, lo, hi, color=c, alpha=0.15, lw=0)
        ax.plot(qd, med, "o-", color=c, ms=3, label=name)
    ax.axhline(1, color="k", lw=0.6)
    ax.axvspan(8, 14, color="0.9", zorder=-1, lw=0)
    ax.set_xlabel(r"$q\,d_{\rm small}$ (object, on the axis)")
    ax.set_ylabel("measured / modelled annulus power")
    ax.set_ylim(0.5, 2.0)
    ax.legend(fontsize=8)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=150)


if __name__ == "__main__":
    main()
