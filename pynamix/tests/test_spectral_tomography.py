"""Species-fraction tomography (:mod:`pynamix.spectral.tomography`)."""

import unittest

import numpy as np

try:
    import torch

    from pynamix.pdvc import DynamicProjectionModel
except ImportError:
    torch = None
from pynamix.pdvc import ConeBeamGeometry, LinearTimeBasis, RegularHexMesh
from pynamix.pdvc.phantom import project_spheres, sphere_packing
from pynamix.spectral.tomography import (
    SpeciesTomography,
    composition_states,
    mixture_fields_by_composition,
    mixture_fields,
    patch_spectra,
    population_field,
    porod_field,
    species_field,
    trilinear_weights,
)


class TestFields(unittest.TestCase):
    """Kernels of the field types and their limits."""

    def test_population_of_one_size_is_a_species(self):
        k = np.linspace(0.1, 5, 50)
        self.assertTrue(np.allclose(population_field([6.0], [1.0]).chi(k), species_field(6.0).chi(k)))

    def test_mixture_reduces_to_self_term_when_dilute(self):
        k = np.linspace(2.0, 4.0, 50)  # q d = 10-20 for d = 5
        dilute = mixture_fields([5.0, 10.0], [1e-4, 1e-4])
        for f, d in zip(dilute, [5.0, 10.0]):
            self.assertTrue(np.allclose(f.chi(k), species_field(d).chi(k), rtol=1e-2, atol=1e-6 * d**3))

    def test_mixture_groups_average_their_members(self):
        k = np.linspace(2.0, 4.0, 20)
        single = mixture_fields([5.0, 6.0], [0.1, 0.1])
        grouped = mixture_fields([5.0, 6.0], [0.1, 0.1], groups=[[0, 1]])
        self.assertTrue(np.allclose(grouped[0].chi(k), 0.5 * (single[0].chi(k) + single[1].chi(k))))

    def test_porod_matches_the_lobe_average_of_a_sphere(self):
        # the self term v F^2 oscillates about the Porod envelope 2 pi (3 / a) / k^4 per unit volume fraction
        d = 4.0
        k = np.linspace(20, 40, 20000) / d
        sphere = species_field(d).chi(k) * k**4
        porod = porod_field().chi(k) * (6 / d) * k**4
        self.assertAlmostEqual(sphere.mean() / porod.mean(), 1.0, delta=0.02)


class TestCompositionStates(unittest.TestCase):
    """Kernels that follow the local composition."""

    def test_state_weights_are_hat_functions(self):
        w = composition_states([0.0, 0.5, 1.0])(np.array([0.0, 0.25, 0.5, 0.9, 1.2]))
        self.assertTrue(np.allclose(w.sum(0), 1))
        self.assertTrue(np.allclose(w[:, 1], [0.5, 0.5, 0]))
        self.assertTrue(np.allclose(w[:, 4], [0, 0, 1]))

    def test_states_match_mixture_at_their_composition(self):
        k = np.linspace(2.0, 4.0, 30)
        f = mixture_fields_by_composition([5.0, 10.0], 0.3, [0.25, 0.75])
        ref = mixture_fields([5.0, 10.0], [0.3 * 0.25, 0.3 * 0.75])
        self.assertEqual(len(f[0].chi), 2)
        self.assertTrue(np.allclose(f[1].chi[1](k), ref[1].chi(k)))


class TestSpeciesTomography(unittest.TestCase):
    """Patch spectra, sphere projector and the species-fraction kernels."""

    def test_white_noise_spectrum_is_removed(self):
        rng = np.random.default_rng(0)
        X = 0.5 * rng.standard_normal((4, 96, 96))
        s = patch_spectra(X, size=32, noise_var=0.25, pixel_mtf=None)
        self.assertTrue(np.allclose(s.q_weights.sum(1), 1))
        # noise PSD is 0.25 per unit area; after subtraction the residual mean is small
        self.assertTrue(abs(s.P.mean()) < 0.02)

    @unittest.skipIf(torch is None, "PyTorch is not installed")
    def test_single_sphere_path_length(self):
        g = ConeBeamGeometry.circular([0.0], 400, 600, (64, 64), 1.0)
        X = project_spheres(np.zeros((1, 3)), np.array([10.0]), g, 0, supersample=4)
        M = 600 / 400
        self.assertTrue(abs(X.max() - 10) < 0.1)
        self.assertTrue(abs(X.sum() / M**2 - np.pi * 10**3 / 6) / (np.pi * 10**3 / 6) < 0.01)

    @staticmethod
    def _setup(n_proj=4, N=48):
        g = ConeBeamGeometry.circular(np.linspace(0, 2 * np.pi, n_proj, endpoint=False), 4 * N, 6 * N, (80, 80), 1.0)
        return g, RegularHexMesh((N, N, N), 2)

    @unittest.skipIf(torch is None, "PyTorch is not installed")
    def test_patch_mean_kernel_matches_projector(self):
        N = 48
        g, mesh = self._setup(N=N)
        R = 0.4 * N

        def inside(p):
            return ((p[:, 0] ** 2 + p[:, 1] ** 2) < R**2) & (np.abs(p[:, 2]) < R)

        f = 0.3 + 0.05 * np.random.default_rng(0).standard_normal(mesh.n_nodes)
        z, y, x = np.meshgrid(*[np.arange(N) - (N - 1) / 2] * 3, indexing="ij")
        pts = np.stack([x, y, z], -1).reshape(-1, 3)
        vol = ((trilinear_weights(pts, mesh.vol_shape, mesh.node_shape) @ f) * inside(pts)).reshape(N, N, N)
        m = DynamicProjectionModel(
            vol, g, mesh, LinearTimeBasis(2), np.zeros(g.n_proj), n_samples=4 * N, cache_rays=False
        )
        with torch.no_grad():
            X = m.project_all(m.zeros()).numpy()
        spec = patch_spectra(X, size=16, q_band=(1.0, 1.5))
        t = SpeciesTomography(
            g, mesh, LinearTimeBasis(1), np.zeros(g.n_proj), [4.0], spec, inside=inside, rays_per_side=6
        )
        t.build(verbose=False)
        mean, _ = t.predict(f[None, None])
        ok = mean > 0.5 * mean.max()
        self.assertTrue(np.median(np.abs(spec.mean[ok] / mean[ok] - 1)) < 0.02)

        # separate, smaller mean patches and a per-species boundary give the same means
        small = patch_spectra(X, size=8, q_band=(0, np.inf))

        def inside2(p):
            return np.stack([inside(p), inside(p)])

        t2 = SpeciesTomography(
            g,
            mesh,
            LinearTimeBasis(1),
            np.zeros(g.n_proj),
            [4.0, 8.0],
            spec,
            means=small,
            inside=inside2,
            mean_rays_per_side=4,
        )
        t2.build(verbose=False)
        mean2, P2 = t2.predict(np.stack([0.5 * f, 0.5 * f])[None])
        ok = mean2 > 0.5 * mean2.max()
        self.assertEqual(P2.shape, (g.n_proj, spec.centres.shape[0], spec.q.size))
        self.assertTrue(np.median(np.abs(small.mean[ok] / mean2[ok] - 1)) < 0.03)

    @unittest.skipIf(torch is None, "PyTorch is not installed")
    def test_dilute_packing_spectrum_matches_self_term(self):
        N, d = 64, 6.0
        R = 0.4 * N
        g = ConeBeamGeometry.circular(np.linspace(0, np.pi, 4, endpoint=False), 4 * N, 6 * N, (96, 96), 1.0)
        rng = np.random.default_rng(1)
        V = np.pi * R**2 * 2 * R
        X = []
        for i in range(g.n_proj):
            c, dd = sphere_packing(R, R, [d], [int(0.1 * V / (np.pi * d**3 / 6))], [lambda x: np.ones(len(x))], rng)
            X.append(project_spheres(c, dd, g, i))
        phi = len(dd) * np.pi * d**3 / 6 / V
        spec = patch_spectra(np.array(X), size=32, stride=16, q_band=(8 / (d * 1.5), 0.6 * np.pi), pixel_mtf=2)
        keep = (np.abs(spec.centres) < 16).all(1)
        spec = spec.select(keep)

        def inside(p):
            return ((p[:, 0] ** 2 + p[:, 1] ** 2) < R**2) & (np.abs(p[:, 2]) < R)

        mesh = RegularHexMesh((N, N, N), 1)
        t = SpeciesTomography(g, mesh, LinearTimeBasis(1), np.zeros(g.n_proj), [d], spec, inside=inside)
        t.build(verbose=False)
        mean, P = t.predict(np.full((1, 1, mesh.n_nodes), phi))
        ratio = spec.P.mean((0, 1)) / P.mean((0, 1))
        self.assertTrue(np.all(np.abs(ratio - 1) < 0.25), ratio)
        self.assertTrue(abs(spec.mean.mean() / mean.mean() - 1) < 0.1)


if __name__ == "__main__":
    unittest.main()
