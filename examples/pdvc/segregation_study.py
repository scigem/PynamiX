"""Sensitivity of the species-fraction reconstruction to noise, resolution and number of projections.

Runs :mod:`segregation` once per setting, varying one parameter at a time
from the default case, and writes ``study.json`` and ``study.png``.
Each run takes several minutes.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import segregation  # noqa: E402

STUDIES = {
    "noise": ("--noise", [0.0, 0.3, 1.0, 3.0, 10.0]),
    "pixel": ("--pixel", [0.75, 1.0, 1.5]),
    "n_proj": ("--n-proj", [90, 180, 360]),
    "seed": ("--seed", [0, 1]),
    # grains per element: a finer mesh on the same data, or smaller grains at a fixed mesh
    "elem": ("--elem", [1, 2, 3, 4]),
    "grain": ("--d-small", [4.0, 5.0, 6.5, 8.0]),
    # size ratio, small grains fixed at 5 voxels
    "ratio": ("--d-large", [7.5, 10.0, 15.0]),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", default="noise,pixel,n_proj")
    ap.add_argument("--out", type=Path, default=Path("output/pdvc/segregation_study"))
    ap.add_argument("--cache", type=Path, default=None, help="data cache (default: <out>/cache)")
    ap.add_argument("extra", nargs="*", help="further arguments passed to every run")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    path = a.out / "study.json"
    results = json.loads(path.read_text()) if path.exists() else {}
    for name in a.studies.split(","):
        flag, values = STUDIES[name]
        for v in values:
            key = f"{name}={v}"
            if key in results:
                continue
            argv = [
                flag,
                str(v),
                "--out",
                str(a.out / key),
                "--no-plot",
                "--cache",
                str(a.cache or a.out / "cache"),
                *a.extra,
            ]
            if name == "pixel":
                # keep the patches the same physical size; the band must stay below 0.6 of Nyquist
                qd_max = 0.6 * np.pi * 1.5 * 5.0 / v
                argv += ["--qd-min", f"{min(8.0, 0.7 * qd_max):.2f}"]
                argv += ["--patch", str(int(round(48 / v))), "--mean-patch", str(int(round(16 / v)))]
            if name == "grain":
                # large grains twice the size; pixel and patches scale with the grains, so the
                # resolution per grain is that of the default case and only the count changes
                px = v / 5.0
                argv += ["--d-large", str(2 * v), "--pixel", f"{px:g}"]
                argv += ["--patch", str(int(round(48 / px))), "--mean-patch", str(int(round(16 / px)))]
            print(f"-- {key}", flush=True)
            res = segregation.run(segregation.parse(argv))
            results[key] = res
            path.write_text(json.dumps(results, indent=1))
            e = res["error_vs_realised"]
            print(
                f"   phi {e['phi']:.4f}  phi_small {e['phi_small']:.4f}  phi_large {e['phi_large']:.4f}  "
                f"large fraction {e['large_fraction']:.3f}  d32 {e['d32']:.2f}",
                flush=True,
            )
    plot(results, a.out)
    plot_grains(results, a.out)


def plot(results, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [n for n in STUDIES if n not in ("elem", "grain") and any(k.startswith(n + "=") for k in results)]
    if not names:
        return
    fig, ax = plt.subplots(1, len(names), figsize=(4 * len(names), 3.4), squeeze=False, constrained_layout=True)
    labels = {
        "noise": "noise std in X (voxels)",
        "pixel": "detector pixel (voxels)",
        "n_proj": "projections per turn",
        "seed": "random seed",
        "elem": "elements per axis",
        "grain": "small-grain diameter (voxels)",
        "ratio": "large-grain diameter (small = 5 voxels)",
    }
    for j, n in enumerate(names):
        keys = sorted((k for k in results if k.startswith(n + "=")), key=lambda k: float(k.split("=")[1]))
        x = [float(k.split("=")[1]) for k in keys]
        for q, style in [("phi", "o-"), ("phi_small", "s--"), ("phi_large", "^--"), ("large_fraction", "d:")]:
            ax[0, j].plot(x, [results[k]["error_vs_realised"][q] for k in keys], style, label=q)
        ax[0, j].set_xlabel(labels[n])
        ax[0, j].set_ylabel("rms error over the specimen")
        ax[0, j].set_ylim(bottom=0)
    ax[0, 0].legend(fontsize=8)
    fig.savefig(out / "study.png", dpi=130)


if __name__ == "__main__":
    main()


def plot_grains(results, out):
    """Errors against the number of grains per element, from the mesh and grain-size studies."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    for name, marker in [("elem", "o"), ("grain", "s")]:
        keys = [k for k in results if k.startswith(name + "=") and "grains_per_element" in results[k]]
        if not keys:
            continue
        keys.sort(key=lambda k: results[k]["grains_per_element"]["large"])
        n = [results[k]["grains_per_element"]["large"] for k in keys]
        lab = "mesh varied" if name == "elem" else "grain size varied"
        for q, c in [("phi", "C0"), ("phi_small", "C1"), ("phi_large", "C2")]:
            e = [results[k]["error_vs_realised"][q] for k in keys]
            ax[0].loglog(n, e, marker + "-", color=c, label=f"{q}, {lab}", mfc="none" if name == "grain" else c)
        e = [results[k]["error_vs_realised"]["large_fraction"] for k in keys]
        ax[1].loglog(n, e, marker + "-", color="C3", label=lab, mfc="none" if name == "grain" else "C3")
    for a in ax:
        a.set_xlabel("large grains per element")
    ax[0].set_ylabel("rms error in volume fraction")
    ax[1].set_ylabel("rms error in large-grain fraction")
    ax[0].legend(fontsize=7)
    ax[1].legend(fontsize=7)
    fig.savefig(out / "grains_per_element.png", dpi=130)
