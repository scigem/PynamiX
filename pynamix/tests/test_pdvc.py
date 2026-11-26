"""Projection-based DVC (:mod:`pynamix.pdvc`). Needs PyTorch; skipped without it."""

import unittest

import numpy as np

try:
    import torch

    from pynamix.pdvc import DynamicProjectionModel, GaussNewton
except ImportError:  # optional dependency: pip install pynamix[pdvc]
    torch = None
from pynamix.pdvc import ConeBeamGeometry, LinearTimeBasis, RegularHexMesh
from pynamix.pdvc.phantom import granular_phantom


@unittest.skipIf(torch is None, "PyTorch is not installed")
class TestPDVC(unittest.TestCase):
    """Projector, sensitivity, solver, grey levels, multiscale and field of view."""

    @classmethod
    def setUpClass(cls):
        cls._dtype = torch.get_default_dtype()
        torch.set_default_dtype(torch.float64)

    @classmethod
    def tearDownClass(cls):
        torch.set_default_dtype(cls._dtype)

    def test_projection_conserves_mass(self):
        """In a near-parallel beam the detector-summed projection equals the volume's mass."""
        N = 25
        f0 = granular_phantom((N, N, N), seed=3)
        g = ConeBeamGeometry.circular([0.0, 0.7, 2.1], 1e6, 1e6 + 50, (60, 60), 0.5)
        m = DynamicProjectionModel(
            f0,
            g,
            RegularHexMesh(f0.shape, 1),
            LinearTimeBasis(1),
            [0.0, 0.0, 0.0],
            n_samples=200,
            fix_first_knot=False,
            dtype=torch.float64,
        )
        p = m.project_all(m.zeros()).numpy()
        mass = p.sum((1, 2)) * 0.5**2
        self.assertTrue(np.allclose(mass, f0.sum(), rtol=5e-3))

    def test_time_basis_partition_of_unity(self):
        b = LinearTimeBasis(5)(np.linspace(0, 1, 37))
        self.assertTrue(np.allclose(b.sum(1), 1.0))

    def test_sensitivity_matches_finite_difference(self):
        N = 32
        f0 = granular_phantom((N, N, N), seed=1)
        g = ConeBeamGeometry.circular(np.linspace(0, np.pi, 4), 4 * N, 6 * N, (40, 40), 1.0)
        m = DynamicProjectionModel(
            f0,
            g,
            RegularHexMesh(f0.shape, 2),
            LinearTimeBasis(3),
            np.linspace(0, 1, 4),
            n_samples=2 * N,
            dtype=torch.float64,
        )
        torch.manual_seed(0)
        a = 0.3 * torch.randn(m.param_shape)
        v = torch.randn(m.param_shape)
        jv = m.jvp(v, 0, m.linearise(a, 0))
        eps = 1e-2
        fd = (m.project(a + eps * v, 0) - m.project(a - eps * v, 0)) / (2 * eps)
        self.assertTrue(float((jv - fd).norm() / fd.norm()) < 0.15)

    def test_recovers_rigid_translation(self):
        """A time-linear rigid translation is recovered from one projection per step."""
        N, n_proj = 24, 36
        f0 = granular_phantom((N, N, N), seed=2)
        times = np.linspace(0, 1, n_proj)
        g = ConeBeamGeometry.circular(np.linspace(0, 2 * np.pi, n_proj, endpoint=False), 4 * N, 6 * N, (36, 36), 1.0)
        mesh, basis = RegularHexMesh(f0.shape, 1), LinearTimeBasis(2)
        m = DynamicProjectionModel(f0, g, mesh, basis, times, n_samples=2 * N, dtype=torch.float64)
        truth = m.zeros()
        truth[0, 0], truth[0, 1], truth[0, 2] = 0.8, -0.5, 0.6
        p = m.project_all(truth)
        res = GaussNewton(m, p, reg=0.0, verbose=False).solve(max_iter=10, tol=1e-4)
        self.assertTrue(float((res.a - truth).abs().max()) < 0.02)

    def test_grey_fit_removes_affine_drift(self):
        """Projected-out residual is zero when data differ from the model only by grey terms."""
        from pynamix.pdvc.greylevel import GreyModel

        torch.manual_seed(1)
        P = torch.rand(3, 10, 12) + 1.0
        s = torch.linspace(-1, 1, 12)[None, None, :]
        p = 1.07 * P - 0.03 * P**2 + 0.2 + 0.05 * s
        r, Q, d, theta = GreyModel(["gain", "offset", "tilt", "hardening"]).fit(P, p)
        self.assertTrue(float(r.abs().max()) < 1e-6)
        self.assertTrue(np.allclose(theta[:, 0], 0.07) and np.allclose(theta[:, -1], -0.03))

    def test_binned_geometry_matches_binned_projections(self):
        """Projecting with a binned detector agrees with binning full-resolution projections."""
        from pynamix.pdvc.geometry import bin_projections
        from scipy import ndimage

        N = 32
        f0 = ndimage.gaussian_filter(granular_phantom((N, N, N), seed=4), 2.0)
        g = ConeBeamGeometry.circular([0.4, 1.9], 4 * N, 6 * N, (49, 50), 0.75)

        def mk(geom):
            return DynamicProjectionModel(
                f0,
                geom,
                RegularHexMesh(f0.shape, 1),
                LinearTimeBasis(1),
                [0.0, 0.0],
                n_samples=2 * N,
                fix_first_knot=False,
                dtype=torch.float64,
            )

        full = mk(g).project_all(torch.zeros(1, 3, 2, 2, 2)).numpy()
        coarse = mk(g.binned(4)).project_all(torch.zeros(1, 3, 2, 2, 2)).numpy()
        ref = bin_projections(full, 4)
        self.assertTrue(coarse.shape == ref.shape == (2, 12, 12))
        self.assertTrue(np.abs(coarse - ref).max() < 0.03 * np.abs(ref).max())

    def test_prolongation_is_exact_for_nested_meshes(self):
        coarse, fine = RegularHexMesh((20, 20, 20), 2), RegularHexMesh((20, 20, 20), 4)
        torch.manual_seed(2)
        a = torch.randn(2, 3, *coarse.node_shape)
        af = coarse.prolongate(a, fine)
        pts = (torch.rand(2, 50, 3) - 0.5) * 19
        self.assertTrue(torch.allclose(coarse.interpolate(a, pts), fine.interpolate(af, pts), atol=1e-10))

    def test_multiscale_recovers_large_translation(self):
        """A rigid translation of several voxels is recovered by coarse-to-fine levels."""
        from pynamix.pdvc.multiscale import Level, multiscale_solve

        N, n_proj = 24, 24
        f0 = granular_phantom((N, N, N), seed=5)
        times = np.linspace(0, 1, n_proj)
        g = ConeBeamGeometry.circular(np.linspace(0, 2 * np.pi, n_proj, endpoint=False), 4 * N, 6 * N, (36, 36), 1.0)
        basis = LinearTimeBasis(2)
        gen = DynamicProjectionModel(
            f0, g, RegularHexMesh(f0.shape, 1), basis, times, n_samples=2 * N, dtype=torch.float64
        )
        truth = gen.zeros()
        truth[0, 0], truth[0, 1], truth[0, 2] = 3.0, -2.5, 2.0
        p = gen.project_all(truth)
        levels = [Level(bin=4, blur=2.0, elem=1, tol=1e-2), Level(bin=1, blur=0.0, elem=1, tol=1e-3)]
        results, _ = multiscale_solve(
            f0,
            g,
            p,
            basis,
            times,
            levels=levels,
            model_kw=dict(dtype=torch.float64),
            solver_kw=dict(reg=0.0),
            verbose=False,
        )
        self.assertTrue(float((results[-1].a - truth).abs().max()) < 0.05)

    def test_cone_beam_known_region_shrinks_towards_top(self):
        """The fully seen region of a circular cone-beam scan narrows away from the mid-plane."""
        from pynamix.pdvc.fov import cone_beam_known

        N = 32
        g = ConeBeamGeometry.circular(np.linspace(0, 2 * np.pi, 36, endpoint=False), 2 * N, 3 * N, (40, 60), 1.0)
        k = cone_beam_known(g, (N, N, N), stride=1, every=1)
        width = k.sum((1, 2))
        self.assertTrue(width[N // 2] > width[2] and width[N // 2] > width[-3])
        self.assertTrue(not k[0].any() and k[N // 2, N // 2, N // 2])

    def test_mask_drops_rays_through_unknown_material(self):
        """Rays meeting material above the known region are dropped; the rest are kept."""
        from pynamix.pdvc.fov import FieldOfView

        N = 24
        f0 = granular_phantom((N, N, N), seed=6)
        g = ConeBeamGeometry.circular([0.0], 4 * N, 6 * N, (36, 36), 1.0)
        m = DynamicProjectionModel(
            f0, g, RegularHexMesh(f0.shape, 1), LinearTimeBasis(1), [0.0], n_samples=2 * N, fix_first_knot=False
        )
        known = np.ones(f0.shape, bool)
        known[-4:] = False  # top slab unknown
        w = FieldOfView(known=known, envelope_radius=0.4 * N, dilate=0).ray_weights(m, m.zeros(), 0)[0]
        # rows run along +z: the last rows see the unknown top slab; the first two
        # rows leave the bottom of the volume array, where nothing is known either
        self.assertTrue(w[-8:].sum() < w[-8:].numel() and w[2:20].min() == 1 and w[0].sum() < w[0].numel())
        self.assertTrue(w[:, :3].min() == 1)  # rays outside the envelope are kept

    def test_reference_difference_cancels_static_outside_material(self):
        """With reference projections, unmoved material outside the volume does not bias the fit."""
        N, n_proj = 24, 24
        tall = granular_phantom((2 * N, N, N), seed=7)  # specimen extends beyond the known slab
        f0 = tall[N // 2 : N // 2 + N]
        times = np.linspace(0, 1, n_proj)
        g = ConeBeamGeometry.circular(np.linspace(0, 2 * np.pi, n_proj, endpoint=False), 4 * N, 6 * N, (36, 36), 1.0)
        basis = LinearTimeBasis(2)
        gen = DynamicProjectionModel(tall, g, RegularHexMesh(tall.shape, 1), basis, times, n_samples=4 * N)
        model = DynamicProjectionModel(f0, g, RegularHexMesh(f0.shape, 1), basis, times, n_samples=2 * N)
        p_ref = gen.project_all(gen.zeros())
        outside = p_ref - model.project_all(model.zeros())  # static material outside f0
        truth = model.zeros()
        truth[0, 0] = 0.7  # the known slab translates, the rest does not move
        p = model.project_all(truth) + outside
        naive = GaussNewton(model, p, reg=0.0, verbose=False).solve(max_iter=8, tol=1e-4)
        diff = GaussNewton(model, p, reg=0.0, reference=p_ref, verbose=False).solve(max_iter=8, tol=1e-4)
        err_naive = float((naive.a - truth).abs().max())
        err_diff = float((diff.a - truth).abs().max())
        self.assertTrue(err_diff < 0.02 and err_naive > 5 * err_diff)


if __name__ == "__main__":
    unittest.main()
