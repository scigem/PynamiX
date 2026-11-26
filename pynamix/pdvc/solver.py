"""Matrix-free Gauss-Newton solver for one-projection-per-step P-DVC.

Minimises  1/2 sum_i || G_i(P_i[a]) - p_i ||^2 + lambda/2 || D a ||^2
over the space-time nodal unknowns a, where G_i is an optional per-projection
grey-level model (see :mod:`pdvc.greylevel`, eliminated by variable
projection) and D is the finite-difference gradient of each knot's nodal
field over the mesh. Each Gauss-Newton step solves
(J^T J + lambda D^T D + mu I) da = -g by conjugate gradients, with J v from
the model's sensitivity and J^T w obtained from it by reverse-mode autodiff.
mu is a Levenberg-Marquardt damping: a step that raises the objective is
rejected and retried with mu increased tenfold; an accepted step lowers mu
by a factor of three, down to its initial value.

``reg`` and ``damping`` are relative: both are multiplied by the mean
diagonal of J^T J, estimated once by random probing, so they do not depend
on the grey-level scale or the number of projections. Nodes whose support
holds no material (mesh corners outside a cylindrical specimen) have no
data sensitivity and are fixed only by these terms.
"""

from dataclasses import dataclass, field

import torch

from .greylevel import GreyModel


def _grad_penalty(a):
    """1/2 sum of squared nodal differences along each mesh axis."""
    s = 0.0
    for dim in (2, 3, 4):
        if a.shape[dim] > 1:
            s = s + 0.5 * (torch.diff(a, dim=dim) ** 2).sum()
    return s


@dataclass
class GNResult:
    a: torch.Tensor
    grey: torch.Tensor | None = None  # (n_proj, n_params) fitted grey-level parameters
    grey_names: list = field(default_factory=list)
    history: list = field(default_factory=list)


class GaussNewton:
    def __init__(
        self,
        model,
        measured,
        reg=1e-2,
        damping=1e-3,
        grey=None,
        fov=None,
        reference=None,
        cg_iter=60,
        cg_tol=1e-3,
        verbose=True,
    ):
        """``grey``: None, a :class:`GreyModel`, or an iterable of its term names.

        ``fov``: optional :class:`pdvc.fov.FieldOfView` for pixel masking.
        ``reference``: optional projections of the undeformed specimen at the
        same poses as ``measured``. The fit then targets
        measured - reference + P[0], from which static material outside the
        known region cancels.
        """
        self.m = model
        self.p = torch.as_tensor(measured, device=model.device, dtype=model.dtype)
        if reference is not None:
            ref = torch.as_tensor(reference, device=model.device, dtype=model.dtype)
            self.p = self.p - ref + model.project_all(model.zeros())
        self.fov = fov
        self.reg_rel, self.damp_rel = reg, damping
        self.reg, self.damp, self.damp0 = 0.0, 0.0, 0.0
        if grey is not None and not isinstance(grey, GreyModel):
            grey = GreyModel(grey)
        self.grey = grey
        self.cg_iter, self.cg_tol = cg_iter, cg_tol
        self.verbose = verbose

    # --- residual and grey-level fit -------------------------------------

    def _weights(self, a):
        """Per-chunk pixel weights from the field-of-view model, or None."""
        if self.fov is None:
            return [None] * len(self.m.chunks)
        return [self.fov.ray_weights(self.m, a, c) for c in range(len(self.m.chunks))]

    def _evaluate(self, a, weights):
        """Residual per chunk with grey terms fitted; returns (states, ssq, n_kept).

        Each state is (r, proj, d, theta, w); proj, d and theta are None
        without a grey model and w is None without a field-of-view model.
        """
        states, ssq, n_kept = [], 0.0, 0
        with torch.no_grad():
            for c, idx in enumerate(self.m.chunks):
                P = self.m.project(a, c)
                p = self.p[idx]
                w = weights[c]
                if self.grey is None:
                    r = P - p if w is None else w * (P - p)
                    st = (r, None, None, None, w)
                else:
                    st = (*self.grey.fit(P, p, w), w)
                ssq += float((st[0] ** 2).sum())
                n_kept += P.numel() if w is None else int(w.sum())
                states.append(st)
        return states, ssq, n_kept

    def _objective(self, a, ssq):
        reg = self.reg * float(_grad_penalty(a)) if self.reg else 0.0
        return 0.5 * ssq + reg

    # --- linear operators ------------------------------------------------

    def _jv(self, x, c, lin, st):
        jv = self.m.jvp(x, c, lin)
        _, proj, d, _, w = st
        if w is not None:
            jv = w * jv
        if proj is not None:
            jv = GreyModel.project_out(proj, d * jv)
        return jv

    def _reg_hvp(self, v):
        if self.reg == 0:
            return torch.zeros_like(v)
        x = v.detach().requires_grad_(True)
        (g,) = torch.autograd.grad(_grad_penalty(x), x)
        return self.reg * g

    def _data_matvec(self, v, lins, states):
        """J^T J v."""
        out = torch.zeros_like(v)
        x = v.detach().requires_grad_(True)
        for c, (lin, st) in enumerate(zip(lins, states)):
            jv = self._jv(x, c, lin, st)
            (g,) = torch.autograd.grad(0.5 * (jv**2).sum(), x)
            out += g
        return out

    def _diag_scale(self, lins, states, n_probe=4):
        """Mean diagonal of J^T J from Rademacher probes."""
        s = 0.0
        for _ in range(n_probe):
            v = torch.randint(0, 2, self.m.param_shape, device=self.m.device).to(self.m.dtype) * 2 - 1
            s += float((v * self._data_matvec(v, lins, states)).sum()) / v.numel()
        return s / n_probe

    def _gradient(self, a, lins, states):
        """J^T r + lambda D^T D a."""
        out = torch.zeros_like(a)
        x = torch.zeros_like(a).requires_grad_(True)
        for c, (lin, st) in enumerate(zip(lins, states)):
            (g,) = torch.autograd.grad((self._jv(x, c, lin, st) * st[0]).sum(), x)
            out += g
        return out + self._reg_hvp(a)

    def _cg(self, b, lins, states, damp):
        x = torch.zeros_like(b)
        r = b.clone()
        p = r.clone()
        rr = float((r * r).sum())
        rr0 = rr
        for k in range(self.cg_iter):
            Ap = self._data_matvec(p, lins, states) + self._reg_hvp(p) + damp * p
            alpha = rr / float((p * Ap).sum())
            x += alpha * p
            r -= alpha * Ap
            rr_new = float((r * r).sum())
            if rr_new <= (self.cg_tol**2) * rr0:
                break
            p = r + (rr_new / rr) * p
            rr = rr_new
        return x, k + 1

    # --- driver ----------------------------------------------------------

    def solve(self, a0=None, max_iter=20, tol=1e-3, max_reject=2, callback=None):
        """Iterate until the largest accepted nodal update is below ``tol`` voxels."""
        a = self.m.zeros() if a0 is None else a0.clone()
        hist = []
        if self.fov is not None:
            self.fov.reset()
        for it in range(max_iter):
            weights = self._weights(a)
            states, ssq, n_kept = self._evaluate(a, weights)
            with torch.no_grad():
                lins = [self.m.linearise(a, c) for c in range(len(self.m.chunks))]
            if it == 0:
                scale = self._diag_scale(lins, states)
                self.reg = self.reg_rel * scale
                self.damp = self.damp0 = self.damp_rel * scale
            obj = self._objective(a, ssq)
            g = self._gradient(a, lins, states)
            for n_rej in range(max_reject + 1):
                da, n_cg = self._cg(-g, lins, states, self.damp)
                new_states, new_ssq, _ = self._evaluate(a + da, weights)
                if self._objective(a + da, new_ssq) <= obj:
                    break
                self.damp *= 10
            else:
                if self.verbose:
                    print(f"GN {it:2d}  no descent after {max_reject} damping increases; stopping")
                hist.append(
                    dict(
                        iter=it,
                        rms_residual=(ssq / max(n_kept, 1)) ** 0.5,
                        max_update=0.0,
                        cg_iter=n_cg,
                        rejected=max_reject + 1,
                        kept_fraction=n_kept / self.p.numel(),
                    )
                )
                break
            a = a + da
            states, ssq = new_states, new_ssq
            self.damp = max(self.damp / 3, self.damp0)
            step = float(da.abs().max())
            rec = dict(
                iter=it,
                rms_residual=(ssq / max(n_kept, 1)) ** 0.5,
                max_update=step,
                cg_iter=n_cg,
                rejected=n_rej,
                kept_fraction=n_kept / self.p.numel(),
            )
            hist.append(rec)
            if self.verbose:
                print(
                    f"GN {it:2d}  rms residual {rec['rms_residual']:.4e}  "
                    f"max |da| {step:.3e} vox  CG {n_cg}"
                    + (f"  kept {rec['kept_fraction']:.1%}" if self.fov is not None else "")
                    + (f"  ({n_rej} rejected)" if n_rej else "")
                )
            if callback is not None:
                callback(a, rec)
            if step < tol:
                break
        grey = None
        if self.grey is not None:
            grey = torch.cat([st[3] for st in states])
        return GNResult(a, grey, list(self.grey.names) if self.grey else [], hist)
