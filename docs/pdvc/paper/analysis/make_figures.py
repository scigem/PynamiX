"""Figures for the paper, with panel letters only; descriptions are in the captions.

Run from the repository root after metrics.py:

    python docs/pdvc/paper/analysis/make_figures.py
"""

import json
import string
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "output" / "pdvc"
FIG = Path(__file__).resolve().parents[1] / "figures"
sys.path.insert(0, str(ROOT / "examples" / "pdvc"))

from pynamix.pdvc import LinearTimeBasis, RegularHexMesh  # noqa: E402
from pynamix.spectral.tomography import trilinear_weights  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7})


def letter(ax, i, color="k", outside=False):
    x, y, va = (-0.02, 1.02, "bottom") if outside else (0.03, 0.95, "top")
    ax.text(x, y, f"({string.ascii_lowercase[i]})", transform=ax.transAxes, va=va, ha="left", color=color)


def sections(run, rows, path, taus=(0.0, 0.5, 1.0)):
    """Vertical sections (y = 0) of realised and reconstructed fields; rows chosen from ('phi', 'c')."""
    a = json.loads((OUT / run / "result.json").read_text())["args"]
    f = np.load(OUT / run / "fields.npz")
    R, H, dl = a["radius"], a["half_height"], a["d_large"]
    vol = (int(2 * H + dl) + 1,) * 3
    mesh_r = RegularHexMesh(vol, a.get("ref_elem") or a["elem"])
    basis_r = LinearTimeBasis(a.get("ref_knots") or a["knots"])
    mesh_p, basis_p = RegularHexMesh(vol, a["elem"]), LinearTimeBasis(a["knots"])
    xs, zs = np.linspace(-R, R, 81), np.linspace(-H, H, 81)
    Z, Xg = np.meshgrid(zs, xs, indexing="ij")
    pts = np.stack([Xg.ravel(), np.zeros(Xg.size), Z.ravel()], 1)
    grid = []
    for q in rows:
        for name, fld, mesh, basis in (
            ("realised", f["realised"], mesh_r, basis_r),
            ("reconstructed", f["phi"], mesh_p, basis_p),
        ):
            N = trilinear_weights(pts, mesh.vol_shape, mesh.node_shape)
            row = []
            for t in taus:
                p = np.einsum("k,kan->an", basis(np.array([t]))[0], fld[:, :2])
                tot = N @ p.sum(0)
                row.append((tot if q == "phi" else (N @ p[1]) / tot).reshape(Z.shape))
            grid.append((q, row))
    fig, ax = plt.subplots(len(grid), len(taus), figsize=(6.0, 2.0 * len(grid)), constrained_layout=True)
    ax = np.atleast_2d(ax)
    i = 0
    ims = {}
    for r, (q, row) in enumerate(grid):
        for c, img in enumerate(row):
            kw = dict(vmin=0.2, vmax=0.36, cmap="cividis") if q == "phi" else dict(vmin=0, vmax=1, cmap="viridis")
            ims[q] = ax[r, c].imshow(img, origin="lower", extent=[-R, R, -H, H], **kw)
            ax[r, c].set_xticks([])
            ax[r, c].set_yticks([])
            letter(ax[r, c], i, "w")
            i += 1
    labels = {"phi": r"$\phi$", "c": r"$c=\phi_{\rm large}/\phi$"}
    for q, im in ims.items():
        rws = [r for r, (qq, _) in enumerate(grid) if qq == q]
        fig.colorbar(im, ax=ax[rws].ravel().tolist(), shrink=0.8, label=labels[q])
    fig.savefig(path)
    plt.close(fig)


def grains(path):
    m = json.loads((FIG / "metrics.json").read_text())["grains_per_element"]["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(6.4, 3.1), constrained_layout=True)
    for j, sp in enumerate(("small", "large")):
        for prefix, marker, lab in (("elem=", "o", "mesh refined"), ("grain=", "s", "grains smaller")):
            rows = sorted((v for k, v in m.items() if k.startswith(prefix)), key=lambda v: v[f"n_{sp}"])
            n = [v[f"n_{sp}"] for v in rows]
            fill = "C0" if prefix == "elem=" else "none"
            ax[j].loglog(n, [v[f"error_{sp}"] for v in rows], marker, color="C0", mfc=fill, label=f"error, {lab}")
            ax[j].loglog(
                n, [v[f"estimate_{sp}"] for v in rows], marker + "--", color="C3", mfc="none", label=f"estimate, {lab}"
            )
        ax[j].set_xlabel(f"{sp} grains per element")
        ax[j].set_ylabel(rf"rms error in $\phi_{{\rm {sp}}}$")
        letter(ax[j], j, outside=True)
    handles, labels = ax[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False)
    fig.savefig(path)
    plt.close(fig)


def spectra(path):
    z = np.load(OUT / "spectra_check.npz")
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.5), constrained_layout=True)
    for k, c, lab in (("self_term", "C1", "self term"), ("with_interference", "C0", "with interference")):
        r = z[k]
        med = np.median(r, (0, 1))
        lo, hi = np.percentile(r, [25, 75], axis=(0, 1))
        ax.fill_between(z["qd"], lo, hi, color=c, alpha=0.15, lw=0)
        ax.plot(z["qd"], med, "o-", color=c, ms=3, label=lab)
    ax.axhline(1, color="k", lw=0.6)
    ax.axvspan(8.5, 14, color="0.92", zorder=-1, lw=0)
    ax.set_xlabel(r"$q\,d_{\rm small}$")
    ax.set_ylabel("measured / modelled power")
    ax.set_ylim(0.5, 1.6)
    ax.legend(loc="upper right")
    fig.savefig(path)
    plt.close(fig)


def main():
    spectra(FIG / "spectra.pdf")
    sections("segregation", ("phi", "phi", "c", "c")[::2], FIG / "fields.pdf")
    sections("study4/rotating/3d", ("c",), FIG / "rotating.pdf")
    grains(FIG / "grains.pdf")
    print("figures in", FIG)


if __name__ == "__main__":
    main()
