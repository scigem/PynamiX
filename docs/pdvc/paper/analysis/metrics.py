"""Collect every number quoted in the paper into figures/metrics.json.

Reads the outputs of the synthetic runs under output/pdvc (made by the
scripts in examples/pdvc) and computes the few derived quantities the text
uses. Each entry records its source directory. Run from the repository root:

    python docs/pdvc/paper/analysis/metrics.py
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "output" / "pdvc"
sys.path.insert(0, str(ROOT / "examples" / "pdvc"))
import grain_scaling  # noqa: E402
import segregation as S  # noqa: E402

from pynamix.pdvc import LinearTimeBasis, RegularHexMesh  # noqa: E402


def result(path):
    p = OUT / path
    p = p / "result.json" if p.is_dir() else p
    try:
        return json.loads(p.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def errors(path, source=None):
    r = result(path)
    if r is None:
        return None
    e = dict(r["error_vs_realised"])
    e["source"] = source or str(path)
    return e


def study(path):
    return json.loads((OUT / path).read_text())


def fields_at_points(run, n_tau=11):
    """Reconstructed and realised field values (n_tau, 2, n) of a run directory."""
    r = result(run)
    a = r["args"]
    f = np.load(OUT / run / "fields.npz")
    R, H, dl = a["radius"], a["half_height"], a["d_large"]
    vol = (int(2 * H + dl) + 1,) * 3
    ref_elem = a.get("ref_elem") or a["elem"]
    mesh_r, basis = RegularHexMesh(vol, ref_elem), LinearTimeBasis(a.get("ref_knots") or a["knots"])
    mesh_p = RegularHexMesh(vol, (a["elem"], 1, 1) if a.get("layers") else a["elem"])
    vp = S.field_values(f["phi"], mesh_p, LinearTimeBasis(a["knots"]), R, H, n_tau=n_tau)
    vr = S.field_values(f["realised"], mesh_r, basis, R, H, n_tau=n_tau)
    return vp, vr, (R, H, vol)


def main():
    m = {}
    d = result("segregation")
    a = d["args"]
    m["setup"] = dict(
        radius=a["radius"],
        half_height=a["half_height"],
        d_small=a["d_small"],
        d_large=a["d_large"],
        phi=a["phi"],
        n_proj=a["n_proj"],
        packings=a["packings"],
        pixel=a["pixel"],
        noise=a["noise"],
        patch=a["patch"],
        mean_patch=a["mean_patch"],
        elem=a["elem"],
        knots=a["knots"],
        reg=a["reg"],
        source_distance=10 * a["radius"],
        detector_distance=15 * a["radius"],
        magnification=1.5,
        n_patches=d["n_patches"],
        n_mean_patches=d["n_mean_patches"],
        n_annuli=d["n_annuli"],
        qd_band=d["qd_small"],
        grains_per_element=d["grains_per_element"],
        mean_grains=float(
            np.mean(
                [
                    len(z)
                    for k, z in np.load(
                        OUT / "segregation_study" / "cache" / "data_128.0_128.0_180_12_5.0_10.0_0.3_1.0_0.npz"
                    ).items()
                    if k.startswith("d")
                ]
            )
        ),
    )

    # forward model
    z = np.load(OUT / "spectra_check.npz")
    fm = {"qd": z["qd"].tolist(), "source": "spectra_check.npz"}
    for k in ("self_term", "with_interference"):
        med = np.median(z[k], (0, 1))
        fm[k] = dict(
            median=med.tolist(),
            min=float(med.min()),
            max=float(med.max()),
            band_min=float(med[z["qd"] >= 8.5].min()),
            band_max=float(med[z["qd"] >= 8.5].max()),
        )
    m["forward_model"] = fm

    # default case, with and without interference
    m["default"] = errors("segregation")
    m["default_self_term"] = errors("study5/selfterm_default")
    s2 = study("study2/study.json")
    m["pixel_1.5"] = dict(s2["pixel=1.5"]["error_vs_realised"], source="study2 pixel=1.5")
    m["pixel_1.5_self_term"] = errors("study5/selfterm_px15")

    # correlation of the two species' errors
    corr = {}
    for key in ("elem=2", "elem=3", "elem=4", "grain=8.0"):
        vp, vr, _ = fields_at_points(f"study2/elem_grain/{key}")
        e = vp - vr
        corr[key] = float(np.corrcoef(e[:, 0].ravel(), e[:, 1].ravel())[0, 1])
    m["species_error_correlation"] = corr

    # grains per element against the counting estimate
    eg = study("study2/elem_grain/study.json")
    rows = {}
    S0 = grain_scaling.s0(a["d_small"], a["d_large"], a["phi"])
    for key, r in eg.items():
        ra = r["args"]
        C = grain_scaling.mesh_constant(ra["radius"], ra["half_height"], ra["d_large"], ra["elem"])
        n = r["grains_per_element"]
        S0k = grain_scaling.s0(ra["d_small"], ra["d_large"], ra["phi"])
        est = {s: ra["phi"] / 2 * np.sqrt(S0k[j] * C / n[s]) for j, s in enumerate(("small", "large"))}
        e = r["error_vs_realised"]
        rows[key] = dict(
            C=C,
            n_small=n["small"],
            n_large=n["large"],
            estimate_small=float(est["small"]),
            estimate_large=float(est["large"]),
            error_small=e["phi_small"],
            error_large=e["phi_large"],
            ratio_small=e["phi_small"] / est["small"],
        )
    m["grains_per_element"] = dict(rows=rows, S0_small=S0[0], S0_large=S0[1], source="study2/elem_grain")

    # sensitivity
    sens = {k: v["error_vs_realised"]["large_fraction"] for k, v in s2.items()}
    sens.update({k: v["error_vs_realised"]["large_fraction"] for k, v in study("study2/ratio/study.json").items()})
    seeds = [result(f"study4/coverage/seed{i}")["error_vs_realised"]["large_fraction"] for i in range(4)]
    m["sensitivity"] = dict(large_fraction=sens, seeds=seeds, source="study2, study4/coverage")
    m["ratio_large_per_element"] = {
        k: v["grains_per_element"]["large"] for k, v in study("study2/ratio/study.json").items()
    }

    # polydispersity
    m["polydisperse"] = {
        "cv0.25_pixel1": errors("study3/poly/cv025_species_structure"),
        "cv0.25_pixel0.75": errors("study4/poly/cv0.25_species"),
        "cv0.5_pixel0.75": errors("study4/poly/cv0.5_species"),
        "porod_single_pixel0.75": errors("study3/poly/cv0_porod_px075"),
        "porod_cv0.25_pixel0.75": errors("study4/poly/cv0.25_porod"),
        "porod_cv0.5_pixel0.75": errors("study4/poly/cv0.5_porod"),
    }

    # time resolution
    m["time"] = {
        p.stem: result(f"study3/time/{p.stem}")["error_vs_realised"]["large_fraction"]
        for p in sorted((OUT / "study3/time").glob("*.json"))
        if result(f"study3/time/{p.stem}")
    }
    ends = {}
    for name, run in (("rotating", "study4/rotating/3d"), ("default", "study3/layers/3d")):
        vp, vr, _ = fields_at_points(run)
        c = lambda v: v[:, 1] / v.sum(1)  # noqa: E731
        ends[name] = np.sqrt(((c(vp) - c(vr)) ** 2).mean(1)).tolist()
    m["error_against_time"] = ends

    # 3D against layers
    m["layers"] = {
        "default_3d": errors("study3/layers/3d"),
        "default_layers_z2": errors("study4/layers/layered_z2"),
        "default_layers_z4": errors("study4/layers/layered_z4"),
        "lateral_3d": errors("study4/lateral/3d"),
        "lateral_layers_z4": errors("study4/layers/lateral_z4"),
        "rotating_3d": errors("study4/rotating/3d"),
        "rotating_layers_z4": errors("study4/layers/rotating_z4"),
    }
    vp, vr, (R, H, vol) = fields_at_points("study3/layers/3d")
    q, _ = S.grid(R, H, 0.0, 8.0)
    q = q[(np.hypot(q[:, 0], q[:, 1]) < R) & (np.abs(q[:, 2]) < H)]
    # constant across each layer: project onto functions of height only
    zl = np.unique(np.round(q[:, 2], 6))
    lay = np.empty_like(vr)
    for t in range(vr.shape[0]):
        for s_ in range(2):
            v = vr[t, s_]
            means = {z_: v[np.round(q[:, 2], 6) == z_].mean() for z_ in zl}
            lay[t, s_] = [means[z_] for z_ in np.round(q[:, 2], 6)]
    m["default_lateral_signal_phi_large"] = float(np.sqrt(np.mean((vr[:, 1] - lay[:, 1]) ** 2)))

    # local interference, wall, coverage
    m["local_interference"] = {
        "default": errors("study4/local/default_local"),
        "lateral": errors("study4/lateral/3d_local"),
    }
    m["wall"] = {k: errors(f"study2/wall/{k}") for k in ("nowall", "wall_nomodel", "wall_model", "nowall_model")}
    cov = []
    for i in range(4):
        r = result(f"study4/coverage/seed{i}")
        cov.append(dict(posterior=r["uncertainty"], with_floor=r["uncertainty_with_floor"]))
    m["coverage"] = cov

    # windows
    win = {}
    for p in sorted((OUT / "study5").glob("e*.json")):
        r = result(f"study5/{p.stem}")
        if r:
            win[p.stem] = dict(r["error_vs_realised"], chi2=r["chi2"])
    m["windows"] = win
    eg = study("study2/elem_grain/study.json")
    m["window_meshes"] = {
        f"elem={E}": dict(
            hann=eg[f"elem={E}"]["error_vs_realised"], tukey05=result(f"study5/e{E}_tukey05")["error_vs_realised"]
        )
        for E in (1, 2, 3, 4)
    }
    m["window_seeds"] = {
        f"seed={i}": dict(
            hann=result(f"study4/coverage/seed{i}")["error_vs_realised"],
            tukey05=(result(f"study5/tukey05_seed{i}") or result("study5/e2_tukey05"))["error_vs_realised"],
        )
        for i in range(4)
    }

    path = Path(__file__).resolve().parents[1] / "figures" / "metrics.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(m, indent=1, sort_keys=True, default=float))
    print("wrote", path)


if __name__ == "__main__":
    main()
