# Granular fields from rotating radiography

Research note on combining projection-based digital volume correlation with
the spectral analysis of radiographs, for experiments on the DynaTOM.

## Purpose

The DynaTOM rotates its source and detector around a stationary specimen and
records radiographs continuously. Two lines of work are relevant to what it
can measure.

- Projection-based digital volume correlation (P-DVC) measures a
  displacement field in 3D and in time from that stream of radiographs. It
  is implemented in `pynamix.pdvc`.
- The spectral analysis in `pynamix.spectral` measures solid fraction, grain
  size (Sauter mean diameter) and per-species velocities from the Fourier
  spectrum of a radiograph. It works on a single view, so what it returns is
  an average through the thickness of the specimen.

This note describes how the two could be combined so that the rotation of
the DynaTOM turns the spectral quantities into 3D fields that evolve in time.
It is a plan with a first prototype, not a finished method.

## From DIC to P-DVC

### DIC

Digital image correlation compares a reference image with a deformed one.
The deformed image is assumed to be the reference with its texture moved:
f(x) = g(x + u(x)). The displacement field u is described by a few unknowns,
for example values at the nodes of a finite-element mesh, and they are found
by minimising the grey-level difference between the two images. Global
(finite-element) DIC solves for all nodes at once with Gauss–Newton
iterations. Local (subset) DIC solves for each subset separately.

### DVC

Digital volume correlation is the same idea in 3D. The images are volumes
reconstructed by X-ray tomography and the mesh is made of hexahedra. The
difficulty is time. A tomographic scan needs hundreds to thousands of
radiographs over at least half a turn, which takes from seconds to hours, and
the specimen must stay still for all of it. Loading is therefore stopped for
each scan, and anything that happens between scans is not seen.

### P-DVC

Projection-based DVC, developed mainly by Roux, Hild and co-workers, skips
the reconstruction during loading. Only one full scan is made, of the
undeformed specimen, and it gives the reference volume f0. During loading
the gantry keeps turning and each radiograph is taken at a different angle
and a different time. The model for radiograph i is:

1. deform the reference volume with the current displacement estimate,
   f(x, t_i) = f0(x − u(x, t_i));
2. project it along the rays of pose i to get a predicted radiograph;
3. compare it with the measured radiograph.

The displacement is written as a coarse mesh in space times a few basis
functions in time, so the whole loading history has a few hundred to a few
thousand unknowns. A radiograph has 10^5 to 10^6 pixels, so a single
radiograph per time step already gives far more equations than unknowns.
The unknowns are found by Gauss–Newton, as in global DIC, with the image
gradient replaced by the gradient of the reference volume summed along each
ray.

The relation to DIC is close. P-DVC is global DIC in which the "image" is a
line integral of the deformed volume instead of the deformed image itself.
It inherits the same requirements: the texture must still match the
reference, and the motion between the current estimate and the solution must
be smaller than about one feature size.

`pynamix.pdvc` implements this on synthetic data, together with grey-level drift
corrections, a coarse-to-fine scheme and material moving in and out of the
field of view (see `docs/pdvc/README.md`).

### Where P-DVC stops

P-DVC tracks texture. In a granular flow the grains rearrange, and after a
shear strain of order one the texture no longer resembles the reference scan.
It also says nothing directly about solid fraction or grain size. Those are
statistical properties of the packing, not of a particular arrangement of
grains.

## Spectral measurements already available

`pynamix.spectral` treats a radiograph of a dense packing as a random field
and models its power spectrum. In the thick-sample limit the spectrum of the
projected solid path length is

    P(q) = <X> v F²(qd/2) S(q)

with <X> the mean path length through solid, v and d the grain volume and
diameter, F the form factor of a sphere and S the structure factor of the
packing. From this, with no calibration:

- the high-q tail (Porod regime) gives the specific surface and hence the
  Sauter mean diameter d32;
- the high-q self term of a mixture gives each species' volume fraction;
- the cross-spectrum between two radiographs gives a displacement per
  species, because each species weights the spectrum differently.

These methods need no reference scan and no grain-by-grain tracking. They
apply to flowing material where P-DVC fails. Their limitation is geometric:
each radiograph gives one number per patch, averaged through the whole
specimen along the beam. A segregated layer that varies with depth cannot be
distinguished from a uniform mixture.

## Combining the two

### Spectral fields

Grains far apart along a ray are uncorrelated. The spectrum of a patch of
the detector is therefore the sum of contributions from every depth that
the patch's rays pass through. At high q, where the interference terms are
small, the contribution from depth s is proportional to the local volume
fraction of each species, with a known kernel that depends on grain size and
magnification. The patch mean of the path length is the line integral of the
total solid fraction.

The same ray machinery P-DVC uses can then be reused. Write each species'
volume fraction as a coarse mesh in space times basis functions in time.
Every patch of every radiograph gives a few linear equations in those nodal
values, and each radiograph is taken from a different angle. Solving them is
a linear tomography problem for slowly varying fields. Small-angle-scattering
tomography is built on the same principle. Solid fraction and Sauter mean
diameter follow from the species fractions at each point.

No reference scan is needed for this. Only the grain sizes and the
calibrated geometry are needed.

### Species velocities

The cross-spectrum method compares two radiographs at the same pose. On the
DynaTOM consecutive radiographs differ by a small rotation, which shifts
grains at different depths by different amounts. That shift follows from the
geometry and can be removed. The phase of the cross-spectrum then gives a
velocity per species per patch, weighted through the depth by that species'
local volume fraction. This is again a projection of a 3D field, here a
vector field, and the rotating views allow it to be reconstructed in 3D.
Whether all three velocity components are recoverable has to be checked on
synthetic data.

### Mass conservation

Each species obeys ∂φ_a/∂t + ∇·(φ_a v_a) = 0. Fractions and velocities
measured separately must satisfy this, which gives either a consistency check
or a constraint that couples the two fits. The segregation flux, the motion
of one species relative to the mixture, comes out of the same equation.

### Where P-DVC fits

In slow or quasi-static tests the texture survives, and P-DVC gives the full
displacement field. Its volumetric strain and the spectral solid fraction
must then agree, which tests both methods. Spectral velocities could also
provide the initial guess that P-DVC needs when the motion is large.

## First prototype

`pynamix/spectral/tomography.py` implements the species volume fraction part
and `examples/pdvc/segregation.py` tests it on synthetic data.

### Test case

A cylinder of spheres of two sizes (diameters 5 and 10 voxels, radius and
half-height 128 voxels, about 34 000 grains) is imaged while the gantry
turns once, with 180 radiographs. During the turn the solid fraction
decreases with height and a dilated column grows off the axis, and the
initially mixed packing segregates, with large grains moving up and
outwards. Twelve packings, drawn from the target fields at successive times,
stand in for the evolving specimen; each is seen by 15 consecutive
radiographs. The radiographs are exact cone-beam path lengths through every
sphere, with Gaussian detector noise.

Each radiograph is cut into 48-pixel patches for the spectra and 16-pixel
patches for the mean path length. The two species' volume fractions are
reconstructed on a 2 x 2 x 2 element mesh with four time knots, and compared
with the fields the packings actually contain (computed from the grains
themselves, not from the targets the packings were drawn from).

### Results

![Reconstructed and realised fields](figures/segregation.png)

In the default case (pixel = 1 voxel at the detector, noise 0.3 voxels of
path length per pixel) the rms nodal errors are 0.015 in solid fraction,
0.024 and 0.023 in the small and large volume fractions, 0.075 in the
large-grain fraction of the solid and 0.32 voxels in the Sauter mean
diameter (which ranges from 5 to 10). The segregation pattern at the end of
the turn and the dilated column are both recovered.

![Sensitivity to noise, pixel size and projections](figures/segregation_study.png)

Varying one setting at a time:

- Detector noise up to 3 voxels of path length per pixel (about 4% of the
  path length through the centre) changes the errors by less than 0.005; at
  10 voxels (13%) the large-grain fraction error rises from 0.075 to 0.096.
- 90, 180 or 360 radiographs per turn give the same errors.
- A finer pixel (0.75 voxels) gives slightly better results, because the
  usable band reaches higher q. A coarser pixel (1.5 voxels) leaves only
  q d_small < 9, where the interference between grains is not negligible,
  and the species errors roughly triple.
- A second random seed gives similar errors.

Since the errors do not fall with more data or less noise, they are set by
the model and by the test itself, not by counting statistics. Two
contributions have been identified. The spectral model is biased by up to
3% in some annuli, probably from the structure factor, which the model
leaves out; with a size ratio of 2 a 3% bias moves the species fractions by
about 0.013. The realised fields themselves are only defined to about
0.015-0.02 at a node, because an element holds a few hundred large grains.
Which of the two dominates has not been established.

### Lessons from building the test

These took most of the effort and apply equally to real data.

- Grains per element. With fewer than a few hundred grains per element the
  mesh-scale field is not well defined, and no method can recover it. On a
  4 x 4 x 4 mesh the same test gives errors about three times larger, for
  this reason.
- Wall layers. A packing generated against a wall crowds grains into the
  outermost layer, by 15-45% here. Rays that graze the surface run through
  it for long distances, and a smooth field cannot represent it, so the fit
  distorts the whole field (the solid fraction error was 0.06 with the
  layer and 0.015 without). The test now avoids the layer by cutting the
  packings from a larger one. Real specimens have wall layers, so the model
  will need a surface-layer term.
- Error model. A patch mean scatters with the number of grains it holds, so
  its absolute error grows as the square root of the path length, and the
  scatter of a spectrum includes the detector noise that was subtracted
  from it. Weighting by a relative error, or leaving out the noise term,
  made the fit fail at the specimen edge and at high noise.
- Forward model details. The annulus width and the window spectrum must be
  included in the model, otherwise the zeros of the form factor give errors
  of up to 40%.

## Open questions

- Structure factor. Adding it to the kernel (from
  `pynamix.spectral.closure`) should show whether the 3% spectral bias is
  what limits the species fractions.
- Wall layers. A term for the excess or deficit of material within a grain
  diameter of the container wall.
- Specimen size. Real specimens hold far more grains per element than the
  test, so finer meshes should be possible; a larger synthetic case would
  show how the error scales.
- Size ratio and polydispersity. Only a ratio of 2 with two sizes has been
  tested.
- Detector. Real data need flat and dark correction, linearisation for beam
  hardening, the detector blur and the calibrated pose of every radiograph.
  `pynamix.spectral.detector` has much of this for a single view.
- Velocities. The rotation-corrected cross-spectrum and the 3D velocity
  reconstruction are not implemented.
- Validation. A packing whose composition is known in 3D, for example layers
  of two sizes, scanned both statically (full CT) and on the DynaTOM.

## Suggested order of work

1. Read the P-DVC literature, starting from the reviews by Roux, Hild and
   co-workers, and run `examples/pdvc/synthetic.py` to see the method work.
2. Run `examples/pdvc/segregation.py`, then extend it to a monodisperse case with a
   solid-fraction gradient, which tests the Porod and mean terms alone.
3. Replace the synthetic radiographs with DynaTOM radiographs of a static,
   layered bidisperse specimen and compare with its full CT scan.
4. Add the rotation-corrected cross-spectrum and reconstruct velocities on
   synthetic flows.
5. Couple fractions and velocities through mass conservation.

## Glossary

- Cone beam: X-rays from a point source onto a flat detector, so features
  closer to the source appear larger (magnification).
- Pose: the source position and detector position and orientation of one
  radiograph.
- Reference volume: the tomographic reconstruction of the specimen before
  loading.
- Self term: the part of a mixture's spectrum from each grain with itself.
  It does not depend on how the grains are arranged.
- Structure factor: the part of the spectrum from correlations between
  grain positions.
- Sauter mean diameter d32: the diameter of a sphere with the same ratio of
  volume to surface area as the whole mixture.
- Time basis: a small set of functions of time (here piecewise linear) whose
  combination describes how each field evolves during the test.
