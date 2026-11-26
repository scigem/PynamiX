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

`pynamix/spectral/tomography.py` and `examples/pdvc/species.py` implement the species volume
fraction part on synthetic data.

It images a cylinder of spheres of two sizes (diameter ratio 2) while the
gantry turns once. The fraction of large grains varies with height and the
front between large-rich and small-rich material moves during the turn. The
radiographs are exact cone-beam path lengths through every sphere. Each
radiograph is cut into patches, the spectrum of each patch is taken and the
two species' volume fractions are reconstructed on a 4 x 4 x 4 element mesh
with four time knots. The result is compared with the fields the packings
actually realised.

What works:

- The forward model. For a uniform packing the measured patch spectra agree
  with the model to within about 5%. This needed the model to include the
  width of each frequency annulus and the blurring of the spectrum by the
  analysis window. Without them the zeros of the form factor are filled in
  and the error reaches 40%.
- The pipeline from radiographs to fields runs end to end, and the
  reconstructed fields show the imposed trend: large grains on top, spreading
  downward over time.

What does not yet work:

- Accuracy. Nodal errors are about 0.06 in each species' volume fraction and
  about 0.2 in the large-grain fraction of the solid. That is too large to be
  useful.
- Sensitivity. With a size ratio of 2 the two species are told apart by a
  small difference in spectral amplitude, so a 10% error in the spectrum
  moves the large-grain volume fraction by about 0.045. Larger size ratios
  will be easier.
- Patch statistics. A patch holds few large grains, so its spectrum scatters
  more than the standard periodogram statistics predict. The weighting in
  the fit does not yet account for this.
- Solid fraction alone. Using only the patch means, the coarse solid
  fraction of a single packing comes back with errors of about 0.05, against
  0.009 for a smooth field without grains. The reason is not yet understood.
  It is the first thing to investigate, since the species fractions cannot
  be better than this.

Two lessons came from building the test. A packing whose grains are all
inside the container has a depleted layer at the wall, which in a small
specimen changes the interior solid fraction by tens of percent. Real
specimens have this layer too. Separately, the realised packing differs from
the intended fields at mesh scale, so results must be compared with what
was actually generated.

## Open questions

- Interference terms. The model keeps only the self term. Its error in the
  chosen q band at realistic solid fractions should be measured, or the
  structure factor of `pynamix.spectral.closure` added to the kernel.
- Statistics. A patch holds few large grains, so its spectrum fluctuates
  more than the periodogram statistics alone suggest. Weights that include
  grain-number fluctuations are needed for honest uncertainties.
- Resolution. The patch size sets the lateral resolution and the number of
  grains per patch sets the noise. The trade-off between the two, and against
  the number of radiographs per unit time, has not been explored.
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
2. Run `examples/pdvc/species.py`, then extend it to a monodisperse case with a
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
