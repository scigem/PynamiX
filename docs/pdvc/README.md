# pynamix.pdvc

Projection-based digital volume correlation (P-DVC) with one projection per
time step, in the style developed by Roux, Hild and co-workers. A reference
volume is reconstructed from a full scan of the undeformed specimen. During
loading the scanner keeps rotating and each projection is taken at a
different time. A space-time displacement field is fitted directly to that
projection sequence, so no volume is reconstructed during loading.

The projector and solver use PyTorch, so they run on CPU or Apple MPS without a
GPU projector library. PyTorch is an optional dependency:
`pip install "pynamix[pdvc]"`. The geometry, mesh, phantoms and the
species-fraction reconstruction need only numpy and scipy.

## Model

The deformed attenuation at time τ is f(x, τ) = f0(x − u(x, τ)), with f0 the
reference volume and u an Eulerian displacement. The displacement is

    u(x, τ) = Σ_k Σ_n ψ_k(τ) N_n(x) a_kn

with N_n the trilinear (C8) shape functions of a regular hexahedral mesh and
ψ_k piecewise-linear hat functions on equally spaced knots in time. The first
knot is fixed at zero, which ties τ = 0 to the reference state.

Projection i, taken at time τ_i with pose i, is modelled as line integrals of
f(·, τ_i) through the volume. A line integral is evaluated by sampling each
ray at regular steps inside the sphere bounding the volume. At each sample
two trilinear look-ups give u and then f0(x − u). For a regular mesh the C8
interpolation of nodal values is the same trilinear look-up, so the forward
model reduces to two `grid_sample` calls per ray sample.

The unknowns a minimise

    ½ Σ_i ‖P_i[a] − p_i‖² + ½ λ ‖D a‖²

with D the nodal finite-difference gradient. Gauss–Newton is used with the
P-DVC sensitivity

    (J v)_i = − Σ_samples ∇f0(x − u) · (N v)(x) Δs

where ∇f0 is a central-difference gradient volume. J^T w comes from
reverse-mode autodiff through this linear map. Each step solves the normal
equations by conjugate gradients without forming the matrix, with a
Levenberg–Marquardt damping that is raised when a step increases the
objective. λ and the damping are given relative to the mean diagonal of
J^T J. Without them, mesh nodes whose support contains no material, such as
the corners outside a cylindrical specimen, leave the system singular.

## Grey-level corrections

Each measured projection may differ from the projected reference volume P_i
by per-projection terms:

    p_i = P_i + θ_gain P_i + θ_offset + θ_s s + θ_t t + θ_hard P_i²

covering flux or gain drift, additive drift of the log signal, linear ramps
across the detector (s, t scaled to [−1, 1]) and a first-order beam-hardening
term. Any subset can be chosen. For a given displacement these θ_i are a
linear least-squares problem per projection, solved by normal equations in
double precision with a small ridge, and the
displacement problem sees only the part of the residual and sensitivity
orthogonal to them (variable projection with Kaufman's approximation). The
CG system therefore keeps only the displacement unknowns.

## Multiscale

Gauss–Newton needs a starting displacement within about one feature size of
the solution. `multiscale_solve` runs levels from coarse to fine. Each level
blurs the reference volume by a Gaussian of σ voxels and the measured
projections by the same σ mapped to the detector, bins the detector, and may
use a coarser mesh. The result is carried to the next mesh by trilinear
prolongation, which is exact for nested meshes. The box average of binning is
not modelled, so the coarse levels serve only as initialisation. The default
is bin 4 with σ = 2 on a half-resolution mesh, then bin 2 with σ = 1, then
full resolution.

## Material outside the field of view

A specimen usually extends beyond the vertical field of view, and a circular
cone-beam scan reconstructs reliably only the region every projection sees.
The reference volume is therefore known only in a region K
(`fov.cone_beam_known` computes it for any pose list). Rays can meet material
outside K, and deformation carries material across the boundary of K.
`FieldOfView` handles this in two ways.

- Masking: a detector pixel is dropped when its ray samples a reference
  position outside K but inside an envelope where material may be (a
  vertical cylinder about the rotation axis). The mask is recomputed each
  Gauss–Newton iteration and only grows during a solve, so it cannot
  oscillate with the iterates.
- Reference difference: given projections of the undeformed specimen at the
  same poses (`reference=`), the solver fits P_i[a] − P_i[0] to
  p_i − p_ref,i. Unmoved material outside K cancels, so a pixel is dropped
  only when its ray meets unknown material displaced by more than
  `u_threshold`. Outside the mesh the displacement is held at its boundary
  value, so this test extrapolates the boundary motion.

`examples/pdvc/fov.py` compares no treatment, masking alone and the reference
difference on a column that runs past both ends of the reconstructed slab
and is compressed so that material flows into the field of view.

## Species volume fractions

`pynamix.spectral.tomography` reconstructs the volume fraction of each grain species,
φ_a(x, τ), from the same rotating projection sequence, without a reference
volume. The projections are tiled into detector patches and the power
spectrum of each patch is taken, with the conventions of
`pynamix.spectral.psd`. Grains far apart along a ray are uncorrelated, so at
high q, where interference terms are small, a patch spectrum is the sum over
the patch's rays of

    φ_a(x) v_a M(x)² F²(q M(x) d_a / 2) ds

with v_a and d_a the grain volume and diameter, F the sphere form factor and
M the magnification at x. The patch mean of the path length is the line
integral of φ = Σ φ_a. The model for each annulus includes its width and the
convolution with the spectrum of the window, which fills in the zeros of F².
Both are linear in φ_a, so with φ_a written on a C8 mesh times temporal hat
functions the reconstruction is a weighted linear least-squares problem.
Solid fraction and Sauter mean diameter follow from the φ_a.

The patch means come from smaller patches than the spectra, since they
resolve the total solid fraction better while the spectra need many grains
per patch. Grains that straddle the specimen surface are described by a
taper over one grain radius. The weights follow grain-count statistics: a
patch mean has an absolute error growing as the square root of its value,
and an annulus has the periodogram scatter of signal plus detector noise
together with a relative grain-count term. The grain-count terms are scaled
after each pass so that their reduced χ² approaches one.

`examples/pdvc/segregation.py` images a bidisperse sphere packing (34 000
grains, size ratio 2) whose solid fraction varies in 3D and whose species
segregate during one turn, with exact cone-beam path lengths through every
sphere. With grain interference in the kernel (Percus-Yevick partial
structure factors at the mean composition), on a 2 x 2 x 2 mesh, the rms
errors over the specimen against the fields the packings realised are 0.0045
in solid fraction, about 0.007 in each species fraction and 0.023 in the
large-grain fraction of the solid. `examples/pdvc/segregation_study.py`
varies noise, pixel size, number of projections, size ratio, mesh and grain
size; `examples/pdvc/grain_scaling.py` compares the errors with a counting
estimate. The fields can also be population volume fractions for
polydisperse grains, solid fraction and specific surface (Porod tail), or
extra fields confined to a shell at the container wall
(`pynamix.spectral.tomography.Field`).
`docs/pdvc/research-note.md` gives the results and what they depend on.

## Layout

Modules of `pynamix/pdvc/`:


- `geometry.py`: per-projection cone-beam poses (source, detector centre, u
  and v pixel vectors, as in ASTRA `cone_vec`) and ray sampling. A
  rotating-gantry scanner such as the DynaTOM and a rotating stage give the
  same relative geometry, so a calibrated pose list covers both.
- `mesh.py`: regular C8 mesh and temporal hat basis.
- `model.py`: forward projection of the deformed volume and its sensitivity.
- `solver.py`: matrix-free Gauss–Newton.
- `greylevel.py`: per-projection grey-level model and its elimination.
- `multiscale.py`: coarse-to-fine levels.
- `fov.py`: known region of the reference volume and pixel masking.
- `phantom.py`: textured phantom, an imposed compression and shear field,
  random sphere packings and an exact cone-beam projector for spheres.

`pynamix/spectral/tomography.py` holds the patch spectra and the linear
reconstruction of species volume fractions; it uses the window and form factor
of `pynamix.spectral`.

## Synthetic test

```bash
pip install -e ".[dev,pdvc]"
python examples/pdvc/synthetic.py
```

The measured sequence is generated with a 16-element mesh, a temporal knot at
every projection and twice the ray sampling, plus Gaussian noise. The coarser
reconstruction model therefore cannot represent the imposed field exactly.
The script prints nodal displacement errors and writes
`output/pdvc/synthetic/synthetic.png`. `--drift A` adds grey-level drift of
relative amplitude A to the data, `--grey offset,gain,tilt,hardening` fits
those terms and `--multiscale` uses the default levels.

## Not yet done

- Real DynaTOM data: projection reader, flat/dark correction and log
  conversion, and the calibrated per-projection geometry.
- Grey-level terms that vary smoothly in time instead of independently per
  projection, which would reduce their number and their correlation with
  the displacement.
- Beam hardening beyond a per-projection quadratic in P.
- The reference difference assumes the reference projections share the
  dynamic poses. Interpolating reference projections in angle is not
  implemented. With a gain term, gain drift that multiplies the contribution
  of material outside K is not captured.
- Species fractions: errors double at the first and last time knots; the
  Porod model needs a band well into the Porod regime of the smallest
  grains; the wall-layer term has been tested on synthetic packings only.
  Species velocities from the cross-spectra of consecutive projections are
  not implemented. `docs/pdvc/research-note.md` describes the plan.
- Mechanical (equilibrium-gap) regularisation in place of the gradient penalty.
- The sensitivity is evaluated at the current state with ∇f0 sampled at the
  warped points. For sharp features its mismatch with the exact derivative of
  the trilinear forward model is about 20 %; with a CT-like blur of 1.5
  voxels it is 5–9 %.
