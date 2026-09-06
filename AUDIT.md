# Audit of `pynamix` at `bff6378`, 2026-09-07

Done alongside the `pynamix.spectral` migration. Findings are grouped by whether
they are fixed on this branch or still open. Nine of the ten are fixed.

## Fixed here

**1. `import pynamix` fails on a clean install.** `pynamix/io.py` imports
`requests` and `pynamix/plotting.py` imports `IPython`, neither of which was
declared in `pyproject.toml`; `pynamix/__init__.py` star-imports both modules,
so a correct `pip install pynamix` followed by `import pynamix` raised
`ModuleNotFoundError`. `pandas` was also imported and undeclared. All three are
now declared.

**2. Six dev/docs packages were declared as runtime dependencies.** `sphinx`,
`nbsphinx`, `pandoc`, `pre-commit`, `ipykernel` and `ipywidgets` were installed
for every user. They are now in `[project.optional-dependencies]` as
`docs`, `dev` and `interactive`, and `plotting.py` imports `ipywidgets` lazily
inside `hist_GUI` so the package imports headless without them.

**3. Four failing tests in `main_direction`.** `np.linalg.eig` is the general
non-symmetric routine and returns complex results for matrices it cannot certify
as symmetric; those propagated into `np.arctan2`, which raised `TypeError`. A
nematic order tensor is symmetric by construction, so `np.linalg.eigh` is both
correct and guaranteed real. The input is now symmetrised first, which also
makes the function total for the arbitrary 2x2 arrays the existing test feeds
it. All 84 original tests pass.

## Fixed here, and these change published numbers

**4. `hanning_window` is centred one pixel off.** For `patchw=32` the window is
built about `(32.5, 32.5)` while the 64x64 patch's centre is `(31.5, 31.5)`;
row 0 and column 0 are forced to zero while rows and columns 63 are not, so the
support is lopsided. Measured effect on isotropy: **none** -- shifting a
radially symmetric window does not break the isotropy of the power spectrum, and
an isotropic field gives the same orientation magnitude either way (0.0374 in
both cases, which is the statistical floor of the estimator, not a bias from the
window). Fixed anyway: the window is now built about `patchw - 0.5` and is
exactly invariant under a 180 degree rotation.

**5. `angular_binning` and `radial_grid` use unseeded Monte Carlo.** Both
estimate a quantity with a closed form -- the average of `k-hat (x) k-hat`, and
the annular occupancy, over each pixel -- by throwing ~10^7 random points and
caching the result to `.npy`. Two consequences: the coefficients are not
reproducible between users (no seed), and at ~2400 samples per pixel they carry
~2% noise which is then *fixed* for that user. Deterministic subpixel quadrature
would be exact and faster. **Fixed**: both are now deterministic sub-pixel
quadrature, exact to machine precision, and 0.6 s rather than minutes. Against
the old Monte Carlo at its shipped settings the maximum difference is 0.0067,
which is within its own noise -- two unseeded runs differed from each other by
up to 0.0116. The disk cache is gone; `N` is accepted and ignored.

**6. `radial_grid` labels each bin by its upper edge.** Binning uses
`argmin(|r_grid - dst|)`, i.e. nearest node, so bin `i` is centred on node `i`;
`r_grid` is then shifted up by half a bin *after* the binning loop. The reported
radius is therefore high by about half a bin, and the reported wavelength low by
roughly `Delta/(2n)` -- about 1% at a peak sitting at `n = 7`. **Fixed**: each
bin is now labelled with the area-weighted mean radius of the region actually
assigned to it. Bin membership is unchanged.

**7. The reported wavelength is biased low, and the bias is size-dependent.**
Feeding pure sinusoids of known wavelength through `radial_FFT` at `patchw=32`:

| true lambda (px) | reported | error |
|---|---|---|
| 4.0 | 3.87 | -3.2% |
| 6.0 | 5.96 | -0.6% |
| 8.0 | 7.47 | -6.6% |
| 12.0 | 10.83 | -9.8% |
| 16.0 | 15.16 | -5.2% |

Most of this is the coarse wavelength grid: with a 64-pixel patch the available
wavelengths near `lambda = 12` are `64/5 = 12.8` and `64/6 = 10.7`, so 12.0
cannot be reported. The `2 pi q` Jacobian of the annular sum and the half-bin
offset of (6) add to it. A calibration against beads of one size therefore does
not transfer to another size. Larger patches shrink all three effects.

**Partly fixed.** `average_size_map` now fits a parabola through the peak bin
and its two neighbours instead of reporting the bin. On realistic broad peaks
(band-limited isotropic noise, 12 realisations per wavelength) the median error
falls from 3.6% to 3.2% and the worst from 13.7% to 12.7%; on pure sinusoids,
which are the hardest case for interpolation because the peak is delta-like and
the parabola ends up modelling the window rather than the signal, it is a wash.
Set `interpolate=False` for the old behaviour.

The residual few per cent is the `2 pi q` Jacobian, and that is not a defect to
be fixed: the measurand *is* the peak of the annular sum. Removing the Jacobian
would measure the peak of the spectral density instead, which is a different
quantity -- it is the one `pynamix.spectral` uses, via
`radial_average`, and it is why the two halves of the package report different
sizes for the same image.

What is now surfaced rather than silently returned: below about five Fourier
bins from DC the peak sits inside the window's main lobe and on the steep part
of the Jacobian, and the answer runs 10--17% low however it is read. Both
`average_size_map` and `bidisperse_concentration_map` warn when the longest
wavelength sought crosses that line -- roughly `lambda > 0.4 * patchw` in
pixels.

## Also fixed

**8. `orientation_map(padding_mode=...)` was a silent no-op.** The padded frame
was built and then never used -- the patch was taken from `data[ti, ...]`
directly, so the documented padding never happened. `frame` is now the array the
patch is cut from, with the grid offset by `patchw` to keep patch centres where
they were.

**9. `bidisperse_concentration_map` had no guard on a collapsed search band.**
If the wavelength grid was coarse enough that `min_val_a == max_val_a`, the
`np.argmax` ran on an empty slice and raised an unhelpful error; one bin wider
and it returned that bin with no peak-finding at all. Both now raise with a
message naming the fix. This is the regime the accompanying paper found for
`R = 6` at `patchw = 32`.

## Open

**10. Off-by-one loop bounds.** `hanning_window` and `angular_binning` start
their loops at index 1 while `radial_grid` starts at 0. Both functions have since
been rewritten as vectorised quadrature, so the loops are gone; this is recorded
only because the same pattern may exist elsewhere in the package.
