# Paper plan: species fields from rotating radiography

Working plan for a paper on reconstructing 3D, time-resolved solid fraction,
species volume fractions and Sauter mean diameter of a granular mixture from
the spectra of radiographs taken while a cone-beam gantry rotates. Status as
of the synthetic studies in `docs/pdvc/research-note.md`. Experiments are to
be decided.

## Scope

One paper on composition fields: φ(x, t), φ_a(x, t) per species and d32(x, t).
Species velocities from cross-spectra and coupling to P-DVC are left for a
second paper. They are mentioned in the outlook only. Reasons:

- the composition problem is linear and can be validated against a static CT
  scan, so its claims can be proved cleanly;
- velocities need the rotation-corrected cross-spectrum, which does not yet
  exist, and a velocity ground truth, which is harder to obtain;
- P-DVC is an established method and is not the contribution here.

Velocities are harder for geometric reasons. Consecutive radiographs differ
by a small rotation Δθ, so a static grain at distance s from the axis along
the beam shifts on the detector by about M s Δθ. In the synthetic setup
(radius 128 voxels, M = 1.5, 180 radiographs per turn) static grains at the
front and back of the specimen shift by about ±7 pixels relative to each
other, more than a small grain's diameter. The cross-spectrum is then a sum
of phase ramps that disagree along each ray, and it averages away at the
high q where species are told apart. The shift is known and could be
modelled, but the information survives only if q times the spread of shifts
along a ray stays below about π, which needs angular steps of order 0.1°.
Velocity varying along the ray spreads the phases further. Each view sees
only the two components perpendicular to the ray, and the view turns while
the field evolves. A velocity ground truth is also harder to obtain than a
composition one.

Two routes avoid most of this. Once φ_a(x, t) is reliable, species mass
conservation, ∂φ_a/∂t + ∇·(φ_a v_a) = 0, constrains the species fluxes from
the evolving composition alone (at least along composition gradients), which
gives the segregation flux without cross-spectra. For slow flows, frame pairs
one turn apart at the same angle remove the parallax, and P-DVC applies
while the texture persists. The first route could appear in the outlook of
this paper.

## Central claim

The power spectrum of a detector patch of a radiograph of a dense granular
packing is, at high wavenumber, a sum of independent contributions from every
depth the patch's rays cross, each linear in the local species volume
fractions. Radiographs from many angles, as a rotating gantry provides,
therefore determine 3D fields of species volume fraction through a linear
inverse problem. The fields are statistical: no grain is resolved or tracked
and no reference scan is needed. A space-time basis lets the fields evolve
during a single rotation.

The accuracy is limited by how many grains an element of the reconstruction
mesh holds, and hardly at all by detector noise or the number of radiographs.

## Narrative arc

1. Problem. Segregation, dilation and compaction in granular mixtures are 3D
   and time-dependent. Tomography resolves them only when the flow is
   stopped. Radiography is fast but averages through the depth. Image
   correlation methods, including projection-based DVC, need the grain
   texture to persist, which it does not in flow, and give displacement, not
   composition.
2. Existing single-view spectral method. The spectrum of a radiograph gives
   solid fraction, Sauter mean diameter and species fractions without
   calibration (pynamix.spectral), but only as averages along the beam.
3. Observation. Beyond the interference range the spectrum is additive along
   the ray. A rotating scanner sees every point from many directions.
4. Method. Forward model per patch (mean path length and annulus powers),
   space-time C8 basis, weighted linear least squares with an error model
   from grain-count and periodogram statistics.
5. Synthetic evidence. Exact cone-beam projections of evolving bidisperse
   packings: the fields are recovered, and the accuracy depends on grains
   per element and on the spectral band (pixel size), not on noise or the
   number of radiographs. A scaling law for design.
6. Experimental evidence. Static validation against full CT, then a dynamic
   case.
7. Limits and outlook. Wall layers, interference terms, polydispersity,
   velocities, coupling with mass conservation.

## Claims and the evidence each needs

| Claim | Synthetic evidence | Experimental evidence | Status |
|---|---|---|---|
| C1. Patch spectra are additive along rays and the kernel (self term plus Percus-Yevick interference, magnification, annulus width, window) predicts them | Uniform packings: self term within 3%, interference accounts for the rest | Static homogeneous specimen: measured vs predicted spectra at all angles, for two or three solid fractions | Synthetic done; experiment needed |
| C2. The linear reconstruction recovers 3D species fields | Evolving bidisperse packings: 0.0045 in φ, 0.007 in each species fraction, 0.023 in large-grain fraction, 0.10 voxels in d32 | Static layered or graded bidisperse specimen scanned by DynaTOM and by full CT; CT segmented by grain size gives φ_a(x) | Synthetic done; experiment needed |
| C3. The fields can evolve within one rotation | Changes over a whole, half and quarter turn resolved with four knots per turn | A dynamic test with an independent check (CT before and after, conservation of each species' total volume) | Synthetic done; experiment needed |
| C4. Accuracy is set by grains per element, predictably | Mesh and grain-size studies collapse onto a counting estimate with no fitted parameter (large grains); small grains 1.6 times above it | Same static specimen on several meshes against CT | Synthetic done |
| C5. Insensitive to noise and number of radiographs; needs a pixel reaching q d_small of about 10-14 | Noise, projection, pixel and size-ratio studies | Binned and unbinned versions of the same data | Synthetic done |
| C6. Predicted uncertainties are usable | Posterior covariance plus counting floor: near-nominal coverage for φ and φ_large; conservative for φ_small, whose error inherits the large grains' floor | Same against CT | Synthetic done |
| C7. 3D adds to height-only fields | Large-grain fraction error 0.027 against 0.109 for side-to-side segregation, 0.024 against 0.050 for lateral segregation, 0.023 against 0.028-0.032 for a mostly layered phantom | Static specimen with lateral structure | Synthetic done |

## Science still to do

Model and method:

- Ends of the sequence: errors double at the first and last knots.
- Wide size distributions without known populations: the Porod model
  needs a finer pixel than tested so far.
- Wall layer on real packings: the shell field was tested on packings from
  random sequential addition only.

Settled since the first plan: interference in the kernel (needed; the mean
composition suffices), the small-grain errors (inherited from the large
grains through the total), size ratio, polydisperse populations, time
resolution, uncertainty, and the benefit of 3D on phantoms with lateral
structure.

Experimental pipeline:

- Flat and dark correction, linearisation including beam hardening
  (`pynamix.spectral.detector`), detector blur (MTF) measured and divided out.
- Calibrated pose of every radiograph on the DynaTOM.
- Container: its wall adds strong low-q edges and its own material; patches
  that see the wall need masking or a container term.
- Scatter, which adds a smooth background, mostly at low q.
- Data volume and processing time for a full sequence.

## Candidate experiments

Static, for validation:

- Layered bidisperse glass beads in a cylinder (for example large over small
  over mixed), scanned on the DynaTOM and by full CT. Ground truth φ_a(x)
  from CT segmented by grain size.
- A graded solid fraction (loose over dense) of one size, for the mean and
  Porod terms alone.

Dynamic, where the gantry turns around a stationary apparatus:

- Percolation of fines through a coarse bed: composition changes with time
  and depth, slowly enough for a few knots per turn.
- Discharge of a bidisperse silo: segregation in a converging flow.
- A vibrated or tapped column: compaction and size segregation.

Each dynamic test should start and end with a full CT, so at least the first
and last states have ground truth.

## Figures

1. Concept: rays through a packing, a patch, its spectrum as a sum over
   depth, and the rotating gantry.
2. Forward model check: measured vs modelled patch spectra (synthetic, then
   experimental).
3. Synthetic reconstruction: slices of φ and large-grain fraction, realised
   vs reconstructed, at three times.
4. Error against grains per element, mesh and grain-size studies together,
   with the counting estimate.
5. Sensitivity: noise, pixel size, number of radiographs, size ratio, time
   resolution, layers against 3D.
6. Static experiment vs CT.
7. Dynamic experiment.

## Open decisions

- Which dynamic experiment.
- Whether to include the monodisperse case (solid fraction and d32 only) as
  a simpler first demonstration.
- Whether to report φ and the large-grain fraction as the primary outputs,
  with species volume fractions derived, since that is how their errors
  separate.
- Venue: a granular-physics journal if the dynamic experiment carries
  physics, a measurement journal if the method is the main result.

## Literature to read

No references are listed yet; each needs reading before it is cited.
Topics:

- projection-based DVC and space-time regularisation (Roux, Hild and
  co-workers);
- small-angle scattering tomography and dark-field tomography, where
  scattering signals are also reconstructed as line integrals;
- spectral and dynamic-radiography measurements of granular segregation;
- integral-equation structure factors for hard-sphere mixtures;
- wall layering in granular packings;
- the papers underlying `pynamix.spectral`.
