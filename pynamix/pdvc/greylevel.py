"""Per-projection grey-level corrections, eliminated by variable projection.

Each measured projection is modelled as

    p_i = P_i + sum_k theta_ik B_ik

where P_i is the projected deformed reference volume and the columns B_ik are

    gain       P_i            (flux or detector-gain drift)
    offset     1              (additive drift of the log-converted signal)
    tilt       s, t           (linear ramps across the detector, s and t the
                               column and row coordinates scaled to [-1, 1])
    hardening  P_i**2         (first-order beam-hardening / non-linearity)

For fixed displacement unknowns the theta_i are a small linear least-squares
problem per projection, solved exactly. The displacement problem then
sees the residual (I - Pi_i)(P_i - p_i), with Pi_i the orthogonal projector
onto the span of the B_ik, and the sensitivity (I - Pi_i) diag(d_i) J_i with
d_i = 1 + theta_gain + 2 theta_hardening P_i. This drops the derivative of
Pi_i with respect to the displacement (Kaufman's approximation to variable
projection), which is small when the grey corrections are small.
"""

import torch

TERMS = ("gain", "offset", "tilt", "hardening")


class GreyModel:
    def __init__(self, terms=("offset", "gain")):
        terms = tuple(terms)
        for t in terms:
            if t not in TERMS:
                raise ValueError(f"unknown grey-level term {t!r}; choose from {TERMS}")
        self.terms = terms
        self.names = []
        for t in terms:
            self.names += ["tilt_s", "tilt_t"] if t == "tilt" else [t]

    @property
    def n_params(self):
        return len(self.names)

    def columns(self, P):
        """(B, n_pix, n_params) column matrix for simulated projections ``P`` (B, nr, nc)."""
        B, nr, nc = P.shape
        cols = []
        for t in self.terms:
            if t == "gain":
                cols.append(P)
            elif t == "offset":
                cols.append(torch.ones_like(P))
            elif t == "tilt":
                s = torch.linspace(-1, 1, nc, device=P.device, dtype=P.dtype)
                r = torch.linspace(-1, 1, nr, device=P.device, dtype=P.dtype)
                cols.append(s.expand(B, nr, nc))
                cols.append(r[:, None].expand(B, nr, nc))
            elif t == "hardening":
                cols.append(P**2)
        return torch.stack([c.reshape(B, -1) for c in cols], dim=-1)

    def fit(self, P, p, w=None):
        """Fit theta for each projection, using only pixels with weight ``w`` = 1.

        Returns ``r`` (model - measured with grey terms fitted, zero on dropped
        pixels, shape of P), a projector ``proj`` for :meth:`project_out`,
        ``d`` (B, nr, nc) sensitivity scaling and ``theta`` (B, n_params).
        The fit uses normal equations with a ridge of 1e-9 of their mean
        diagonal, so a projection with every pixel dropped gives theta = 0.
        """
        B, dt = P.shape[0], P.dtype
        # normal equations square the conditioning, so work in double precision
        P, p = P.double(), p.double()
        w = w.double() if w is not None else None
        A = self.columns(P)
        e = (p - P).reshape(B, -1, 1)
        if w is not None:
            wf = w.reshape(B, -1, 1)
            A, e = A * wf, e * wf
        M = A.transpose(1, 2) @ A
        n = M.shape[-1]
        ridge = 1e-9 * M.diagonal(dim1=1, dim2=2).mean(-1) + 1e-30
        M = M + ridge[:, None, None] * torch.eye(n, dtype=M.dtype, device=M.device)
        Minv = torch.linalg.inv(M)
        theta = (Minv @ (A.transpose(1, 2) @ e))[..., 0]
        r = -(e - A @ theta[..., None]).reshape(P.shape)
        d = torch.ones_like(P)
        if "gain" in self.names:
            d = d + theta[:, self.names.index("gain"), None, None]
        if "hardening" in self.names:
            d = d + 2 * theta[:, self.names.index("hardening"), None, None] * P
        return r.to(dt), (A.to(dt), Minv.to(dt)), d.to(dt), theta.to(dt)

    @staticmethod
    def project_out(proj, x):
        """(I - A M^-1 A^T) x for x of shape (B, nr, nc)."""
        A, Minv = proj
        B = x.shape[0]
        xf = x.reshape(B, -1, 1)
        return (xf - A @ (Minv @ (A.transpose(1, 2) @ xf))).reshape(x.shape)
