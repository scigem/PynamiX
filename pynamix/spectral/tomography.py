"""Species volume fractions from the spectra of detector patches.

A radiograph of a granular packing is the sum of the projections of its
grains. Grains far apart along a ray are uncorrelated, so the power spectral
density of a detector patch is the sum of contributions from every depth the
patch's rays pass through. Beyond the range of the interference (structure
factor) terms, the contribution of species a at depth s is

    phi_a(x) v_a M(x)^2 F^2(q M(x) d_a / 2) ds

with v_a the grain volume, F the sphere form factor and M the magnification
at x (an object wavenumber k appears on the detector at q = k / M). The patch
mean of the projected solid path length X is the line integral of
phi = sum_a phi_a. Both are linear in the fields phi_a.

Writing phi_a(x, tau) on a C8 mesh times a temporal hat basis turns a
sequence of projections taken at different angles and times into a linear
least-squares problem for the nodal values. The conventions for the
periodogram (absolute amplitude, window-power normalisation, radial mean,
pixel MTF removed per 2D bin) follow :mod:`pynamix.spectral.psd`. The
geometry, mesh and time basis are those of :mod:`pynamix.pdvc`; nothing here
needs PyTorch.
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .psd import taper
from .units import form_factor


@dataclass
class PatchSpectra:
    """Radially averaged spectra of detector patches.

    ``centres``: (n_patch, 2) patch centres in pixels (row, col), relative to
    the detector centre. ``size``: patch width in pixels. ``q``: (n_q,) mean
    |q| of each annulus (rad per unit length). ``n``: (n_q,) 2D bins per
    annulus. ``P``: (n_proj, n_patch, n_q) PSD in length^4, noise removed.
    ``mean``: (n_proj, n_patch) patch mean of X. ``q_nodes`` (n_r,) and
    ``q_weights`` (n_q, n_r): the expected annulus value of an isotropic
    spectrum f(|q|) is ``q_weights @ f(q_nodes)``. The weights include the
    width of the annulus and the convolution of the spectrum with the
    window's own spectrum, which fills in the zeros of the form factor.
    """

    centres: np.ndarray
    size: int
    pitch: float
    q: np.ndarray
    n: np.ndarray
    P: np.ndarray
    mean: np.ndarray
    q_nodes: np.ndarray = None
    q_weights: np.ndarray = None
    noise: np.ndarray = None

    def select(self, keep):
        """Keep a subset of patches (bool or index array)."""
        return PatchSpectra(
            self.centres[keep],
            self.size,
            self.pitch,
            self.q,
            self.n,
            self.P[:, keep],
            self.mean[:, keep],
            self.q_nodes,
            self.q_weights,
            self.noise,
        )


def patch_spectra(
    X, pitch=1.0, size=48, stride=None, q_band=None, noise_var=0.0, pixel_mtf="area", alpha=1.0, upsample=4, reach=6
):
    """Tile projections ``X`` (n_proj, nr, nc) into patches and take their spectra.

    ``X`` is the linearised projected path length through solid. ``noise_var``
    is the variance of white detector noise in X; its PSD (pitch^2 noise_var)
    is subtracted before the pixel MTF is divided out. ``q_band`` = (q_min,
    q_max) in rad per unit length selects the annuli kept. ``pixel_mtf`` is
    "area" for a detector that integrates over square pixels, an integer s
    for pixels that average s x s point samples (a supersampled simulation),
    or None. ``alpha`` is the Tukey taper fraction; the default 1 (Hann)
    keeps leakage from the strong low-q spectrum out of the band. The window
    spectrum is tabulated at ``upsample`` times the bin spacing out to
    ``reach`` bins.
    """
    X = np.asarray(X, dtype=np.float64)
    n_proj, nr, nc = X.shape
    stride = stride or size
    r0s = np.arange(0, nr - size + 1, stride)
    c0s = np.arange(0, nc - size + 1, stride)
    R0, C0 = np.meshgrid(r0s, c0s, indexing="ij")
    R0, C0 = R0.ravel(), C0.ravel()
    centres = np.stack([R0 + (size - 1) / 2 - (nr - 1) / 2, C0 + (size - 1) / 2 - (nc - 1) / 2], 1)

    W = np.outer(taper(size, "tukey", alpha), taper(size, "tukey", alpha))
    q1 = 2 * np.pi * np.fft.fftfreq(size, d=pitch)
    qmag = np.sqrt(q1[:, None] ** 2 + q1[None, :] ** 2)
    if pixel_mtf == "area":
        t1 = np.sinc(q1 * pitch / 2 / np.pi)
    elif pixel_mtf:
        s_ = int(pixel_mtf)
        h = q1 * pitch / 2
        t1 = np.where(np.abs(h) < 1e-12, 1.0, np.sin(h) / (s_ * np.sin(h / s_) + 1e-300))
    else:
        t1 = np.ones_like(q1)
    m2 = (t1[:, None] * t1[None, :]) ** 2
    dq = 2 * np.pi / (size * pitch)
    idx = np.rint(qmag / dq).astype(int)
    q_min, q_max = q_band if q_band is not None else (3 * dq, 0.6 * np.pi / pitch)
    nb = idx.max() + 1
    counts = np.bincount(idx.ravel(), minlength=nb)
    qmean = np.bincount(idx.ravel(), weights=qmag.ravel(), minlength=nb) / np.maximum(counts, 1)
    bins = np.where((counts > 0) & (qmean >= q_min) & (qmean <= q_max))[0]
    if bins.size == 0:
        raise ValueError(f"no annulus of a {size}-pixel patch lies in q_band {q_min:.3g}-{q_max:.3g}")

    P = np.empty((n_proj, R0.size, bins.size))
    mean = np.empty((n_proj, R0.size))
    for j, (r0, c0) in enumerate(zip(R0, C0)):
        x = X[:, r0 : r0 + size, c0 : c0 + size]
        mean[:, j] = x.mean((1, 2))
        mw = (W * x).sum((1, 2)) / W.sum()
        f = np.fft.fft2(W * (x - mw[:, None, None]))
        p2 = pitch**2 * np.abs(f) ** 2 / (W**2).sum() - pitch**2 * noise_var
        p2 = p2 / m2
        for k, b in enumerate(bins):
            P[:, j, k] = p2[:, idx == b].mean(-1)
    q_nodes, q_weights = _annulus_weights(W, pitch, qmag, idx, bins, dq, upsample, reach)
    # noise PSD left in each annulus after the MTF division, before it was subtracted
    noise = np.array([(pitch**2 * noise_var / m2[idx == b]).mean() for b in bins])
    return PatchSpectra(centres, size, pitch, qmean[bins], counts[bins], P, mean, q_nodes, q_weights, noise)


def _annulus_weights(W, pitch, qmag, idx, bins, dq, upsample, reach):
    """Radial weights mapping an isotropic spectrum to expected annulus means."""
    n = W.shape[0] * upsample
    K2 = np.abs(np.fft.fft2(W, s=(n, n))) ** 2
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=pitch)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    near = (np.abs(KX) <= reach * dq) & (np.abs(KY) <= reach * dq)
    kx, ky, wk = KX[near], KY[near], K2[near]
    wk = wk / wk.sum()
    h = dq / upsample
    q1 = 2 * np.pi * np.fft.fftfreq(W.shape[0], d=pitch)
    QX, QY = np.meshgrid(q1, q1, indexing="ij")
    lo = max(0.0, qmag[np.isin(idx, bins)].min() - (reach + 1) * dq)
    hi = qmag[np.isin(idx, bins)].max() + (reach + 1) * dq
    q_nodes = np.arange(lo, hi + h, h)
    weights = np.zeros((bins.size, q_nodes.size))
    for k, b in enumerate(bins):
        sel = idx == b
        r = np.sqrt((QX[sel][:, None] - kx[None]) ** 2 + (QY[sel][:, None] - ky[None]) ** 2).ravel()
        w = np.broadcast_to(wk, (sel.sum(), wk.size)).ravel()
        g = (r - lo) / h
        i0 = np.clip(np.floor(g).astype(int), 0, q_nodes.size - 2)
        f = np.clip(g - i0, 0, 1)
        weights[k] = np.bincount(i0, w * (1 - f), q_nodes.size) + np.bincount(i0 + 1, w * f, q_nodes.size)
        weights[k] /= weights[k].sum()
    return q_nodes, weights


def trilinear_weights(points, vol_shape, node_shape):
    """Sparse (n_points, n_nodes) trilinear weights; points (n, 3) in (x, y, z)."""
    n = points.shape[0]
    rows, cols, vals = [], [], []
    f, i0 = [], []
    inside = np.ones(n, dtype=bool)
    for ax in range(3):  # x, y, z
        N = vol_shape[2 - ax]
        ne = node_shape[2 - ax] - 1
        g = (points[:, ax] + (N - 1) / 2) / ((N - 1) / ne)
        inside &= (g >= 0) & (g <= ne)
        i = np.clip(np.floor(g), 0, ne - 1).astype(int)
        f.append(g - i)
        i0.append(i)
    nz, ny, nx = node_shape
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                w = (f[0] if dx else 1 - f[0]) * (f[1] if dy else 1 - f[1]) * (f[2] if dz else 1 - f[2])
                node = ((i0[2] + dz) * ny + (i0[1] + dy)) * nx + (i0[0] + dx)
                rows.append(np.arange(n))
                cols.append(node)
                vals.append(w * inside)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, nz * ny * nx))


def _gradient_operator(node_shape):
    """Sparse forward-difference gradient on the node grid, stacked over the three axes."""

    def d1(n):
        return sp.diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n))

    nz, ny, nx = node_shape
    return sp.vstack(
        [
            sp.kron(d1(nz), sp.kron(sp.identity(ny), sp.identity(nx))),
            sp.kron(sp.identity(nz), sp.kron(d1(ny), sp.identity(nx))),
            sp.kron(sp.identity(nz), sp.kron(sp.identity(ny), d1(nx))),
        ]
    ).tocsr()


class SpeciesTomography:
    """Linear reconstruction of species volume fractions phi_a(x, tau).

    ``diameters``: grain diameter of each species (voxel units). ``spectra``:
    :class:`PatchSpectra` whose annuli give the spectral rows. ``means``:
    optional second :class:`PatchSpectra`, usually of smaller patches, whose
    patch means give the rows for the total solid fraction; by default the
    means of ``spectra`` are used. Small patches resolve the total solid
    fraction better, while the spectra need patches holding many grains.
    ``inside``: optional function of points (n, 3) returning the fraction of
    material that the specimen boundary lets through at each point, either
    one array (n,) for all species or (n_species, n). Boolean values cut the
    fields off at a sharp boundary; a taper over one grain radius describes
    grains that straddle the surface. Unknowns are nodal values on ``mesh`` for every knot
    of ``time_basis`` and every species.
    """

    def __init__(
        self,
        geometry,
        mesh,
        time_basis,
        times,
        diameters,
        spectra,
        means=None,
        inside=None,
        rays_per_side=4,
        mean_rays_per_side=2,
        n_samples=None,
    ):
        self.geometry, self.mesh, self.basis = geometry, mesh, time_basis
        self.d = np.asarray(diameters, dtype=float)
        self.v = np.pi * self.d**3 / 6
        self.spectra = spectra
        self.means = means if means is not None else spectra
        self.psi = time_basis(times)
        self.inside = inside
        self.m = rays_per_side
        self.m_mean = mean_rays_per_side if means is not None else rays_per_side
        self.n_samples = n_samples or int(np.linalg.norm(mesh.vol_shape))
        self.n_q = spectra.q.size
        self.n_nodes = mesh.n_nodes
        self.n_species = self.d.size
        self.G = None
        self.Gm = None

    def _rays(self, i, centres, size, m):
        """Sample points (n_patch * m * m * ns, 3), weights and patch index for projection i."""
        g = self.geometry
        off = ((np.arange(m) + 0.5) / m - 0.5) * size
        orr, occ = np.meshgrid(off, off, indexing="ij")
        r = centres[:, 0, None] + orr.ravel()[None]
        c = centres[:, 1, None] + occ.ravel()[None]
        pix = g.det[i] + c[..., None] * g.u[i] + r[..., None] * g.v[i]
        src = g.src[i]
        d = pix - src
        d /= np.linalg.norm(d, axis=-1, keepdims=True)
        radius = np.linalg.norm((np.array(self.mesh.vol_shape[::-1]) - 1) / 2) + 1
        b = d @ src
        cc = src @ src - radius**2
        disc = np.sqrt(np.clip(b * b - cc, 0, None))
        t0, t1 = -b - disc, -b + disc
        ns = self.n_samples
        step = (t1 - t0) / ns
        tk = t0[..., None] + (np.arange(ns) + 0.5) * step[..., None]
        pts = (src + tk[..., None] * d[..., None, :]).reshape(-1, 3)
        w = np.repeat(step.reshape(-1), ns) / (m * m)
        w = np.broadcast_to(w, (self.n_species, w.size))
        if self.inside is not None:
            w = w * np.broadcast_to(np.asarray(self.inside(pts), dtype=float), w.shape)
        patch_of = np.repeat(np.arange(centres.shape[0]), m * m * ns)
        keep = w.max(0) > 0
        return pts[keep], w[:, keep], patch_of[keep]

    def _magnification(self, i, pts):
        g = self.geometry
        nrm = np.cross(g.u[i], g.v[i])
        nrm /= np.linalg.norm(nrm)
        return ((g.det[i] - g.src[i]) @ nrm) / ((pts - g.src[i]) @ nrm)

    def _patch_sum(self, n_patch, patch_of, pts, K):
        """(n_patch, n_rows, n_nodes): sum over each patch's samples of K (n_rows, n) times N."""
        N = trilinear_weights(pts, self.mesh.vol_shape, self.mesh.node_shape)
        out = np.empty((n_patch, K.shape[0], self.n_nodes))
        S = sp.csr_matrix((np.ones(pts.shape[0]), (patch_of, np.arange(pts.shape[0]))), shape=(n_patch, pts.shape[0]))
        for r in range(K.shape[0]):
            out[:, r] = (S.multiply(K[r][None]) @ N).toarray()
        return out

    def spectral_table(self, M):
        """(n_species, n_q, len(M)) kernel v_a M^2 <F^2(q M d_a / 2)> per unit volume fraction and length."""
        sw = self.spectra
        if sw.q_nodes is not None:
            used = sw.q_weights.max(0) > 1e-6 * sw.q_weights.max()
            wq, qn = sw.q_weights[:, used], sw.q_nodes[used]
        else:
            wq, qn = np.eye(sw.q.size), sw.q
        T = np.empty((self.n_species, self.n_q, M.size))
        for a in range(self.n_species):
            T[a] = self.v[a] * M[None] ** 2 * (wq @ form_factor(qn[:, None] * M[None] * self.d[a] / 2) ** 2)
        return T

    def kernels(self, i):
        """Spectral sensitivity (n_patch, n_q, n_species, n_nodes) and mean sensitivity (n_mean, n_species, n_nodes).

        The spectral kernel depends on position only through the
        magnification, so it is tabulated on a fine grid of M and interpolated.
        """
        s = self.spectra
        pts, w, pof = self._rays(i, s.centres, s.size, self.m)
        M = self._magnification(i, pts)
        Mg = np.linspace(M.min(), M.max() + 1e-9, 256)
        T = self.spectral_table(Mg)
        K = np.empty((self.n_q * self.n_species, pts.shape[0]))
        for a in range(self.n_species):
            for b in range(self.n_q):
                K[b * self.n_species + a] = np.interp(M, Mg, T[a, b]) * w[a]
        G = self._patch_sum(s.centres.shape[0], pof, pts, K).reshape(-1, self.n_q, self.n_species, self.n_nodes)
        mm = self.means
        pts, w, pof = self._rays(i, mm.centres, mm.size, self.m_mean)
        Gm = self._patch_sum(mm.centres.shape[0], pof, pts, w)
        return G, Gm

    def build(self, verbose=True):
        n = self.geometry.n_proj
        self.G = np.empty((n, self.spectra.centres.shape[0], self.n_q, self.n_species, self.n_nodes), dtype=np.float32)
        self.Gm = np.empty((n, self.means.centres.shape[0], self.n_species, self.n_nodes), dtype=np.float32)
        for i in range(n):
            self.G[i], self.Gm[i] = self.kernels(i)
            if verbose and (i + 1) % max(1, n // 5) == 0:
                print(f"  kernels {i + 1}/{n}")

    def predict(self, phi):
        """Model patch means (n_proj, n_mean) and PSD (n_proj, n_patch, n_q) for nodal fields (n_knots, n_species, n_nodes)."""
        phi_t = np.einsum("tk,kan->tan", self.psi, phi)
        mean = np.einsum("tpan,tan->tp", self.Gm, phi_t)
        P = np.einsum("tpran,tan->tpr", self.G, phi_t)
        return mean, P

    def solve(
        self, reg=1e-3, mean_rel_error=0.02, spectra_weight=1.0, n_irls=4, variance_components=True, verbose=True
    ):
        """Weighted least squares for the nodal fields.

        Each annulus value is the mean of ``n`` exponentially distributed
        periodogram bins, so its standard deviation is at least its
        expectation, plus the detector noise PSD that was subtracted, over
        sqrt(n). Expectations come from the previous model
        prediction (the data on the first pass), iterated ``n_irls`` times.
        The number of grains in a patch fluctuates, so a patch mean has an
        absolute error growing as the square root of its value, and an
        annulus has a relative grain-count error growing as one over the
        square root of the patch mean, on top of the periodogram scatter.
        ``mean_rel_error`` sets both at the thickest patch. With
        ``variance_components`` the grain-count terms of the means and of each
        annulus are rescaled after every pass so that their reduced chi^2
        approaches one. ``reg`` weights a
        gradient penalty on every field, relative to the mean diagonal of the
        data term. ``spectra_weight`` scales the weight of every spectral row
        (0 fits the patch means alone).
        """
        if self.G is None:
            self.build(verbose)
        K, A, Nn = self.psi.shape[1], self.n_species, self.n_nodes
        nu = K * A * Nn
        ym, yP = self.means.mean, self.spectra.P
        em, eP = np.abs(ym).copy(), np.abs(yP).copy()
        cm, cP = 1.0, np.ones(self.n_q)  # variance-component scales
        xs = np.maximum(self.spectra.mean, 1e-3 * self.spectra.mean.max())
        xmax = xs.max()
        L = (lambda D: (D.T @ D).toarray())(_gradient_operator(self.mesh.node_shape))
        R = np.kron(np.eye(K * A), L)
        phi = None
        self.history = []
        for it in range(n_irls):
            # grain-count scatter: absolute error of a patch mean grows as sqrt(X)
            sm = cm * mean_rel_error * np.sqrt(np.maximum(em, 1e-3 * em.max()) * em.max())
            # periodogram scatter 1/n plus grain-count scatter, relative variance ~ 1/X
            rel2 = 1.0 / self.spectra.n + (cP * mean_rel_error) ** 2 * (xmax / xs)[..., None]
            # the periodogram scatter is that of signal plus detector noise
            pos = np.maximum(eP, 1e-3 * eP.max())
            nP = self.spectra.noise if self.spectra.noise is not None else 0.0
            sP = np.sqrt((pos + nP) ** 2 / self.spectra.n + pos**2 * (rel2 - 1.0 / self.spectra.n))
            sP = sP / max(spectra_weight, 1e-12)
            H = np.zeros((nu, nu))
            b = np.zeros(nu)
            for i in range(self.geometry.n_proj):
                Gs = self.G[i].astype(np.float64).reshape(-1, A * Nn) / sP[i].reshape(-1, 1)
                Gm = self.Gm[i].astype(np.float64).reshape(-1, A * Nn) / sm[i][:, None]
                Gi = np.vstack([Gm, Gs])
                yi = np.concatenate([ym[i] / sm[i], (yP[i] / sP[i]).reshape(-1)])
                Hi, bi = Gi.T @ Gi, Gi.T @ yi
                nz = np.nonzero(self.psi[i])[0]
                for k in nz:
                    sk = slice(k * A * Nn, (k + 1) * A * Nn)
                    b[sk] += self.psi[i, k] * bi
                    for m in nz:
                        sl = slice(m * A * Nn, (m + 1) * A * Nn)
                        H[sk, sl] += self.psi[i, k] * self.psi[i, m] * Hi
            scale = np.mean(np.diag(H))
            x = np.linalg.solve(H + reg * scale * R + 1e-9 * scale * np.eye(nu), b)
            phi = x.reshape(K, A, Nn)
            em, eP = self.predict(phi)
            zm = (((ym - em) / sm) ** 2).mean()
            zP = (((yP - eP) / sP) ** 2).mean((0, 1))
            self.history.append(dict(chi2_means=zm, chi2_spectra=zP))
            if verbose:
                print(f"  pass {it}: reduced chi^2 means {zm:.2f}, spectra " + " ".join(f"{z:.2f}" for z in zP))
            if variance_components:
                cm *= np.sqrt(zm)
                if spectra_weight > 0:
                    # rescale the grain-count term so the annulus chi^2 approaches one
                    cP *= np.sqrt(np.maximum(zP, 1e-2))
        return phi
