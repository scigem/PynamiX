# Audit of `pynamix` at `bff6378`, 2026-09-07

Done alongside the `pynamix.spectral` migration. Findings are grouped by whether
they are fixed on this branch, left alone deliberately, or still open.

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

## Left alone, deliberately

These change published numbers, so they are the maintainer's call rather than
ours.

**4. `hanning_window` is centred one pixel off.** For `patchw=32` the window is
built about `(32.5, 32.5)` while the 64x64 patch's centre is `(31.5, 31.5)`;
row 0 and column 0 are forced to zero while rows and columns 63 are not, so the
support is lopsided. Measured effect on isotropy: **none** -- shifting a
radially symmetric window does not break the isotropy of the power spectrum, and
an isotropic field gives the same orientation magnitude either way (0.0374 in
both cases, which is the statistical floor of the estimator, not a bias from the
window). Worth fixing for clarity; not a correctness bug.

**5. `angular_binning` and `radial_grid` use unseeded Monte Carlo.** Both
estimate a quantity with a closed form -- the average of `k-hat (x) k-hat`, and
the annular occupancy, over each pixel -- by throwing ~10^7 random points and
caching the result to `.npy`. Two consequences: the coefficients are not
reproducible between users (no seed), and at ~2400 samples per pixel they carry
~2% noise which is then *fixed* for that user. Deterministic subpixel quadrature
would be exact and faster. Replacing them changes every orientation number
slightly.

**6. `radial_grid` labels each bin by its upper edge.** Binning uses
`argmin(|r_grid - dst|)`, i.e. nearest node, so bin `i` is centred on node `i`;
`r_grid` is then shifted up by half a bin *after* the binning loop. The reported
radius is therefore high by about half a bin, and the reported wavelength low by
roughly `Delta/(2n)` -- about 1% at a peak sitting at `n = 7`.

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

## Open

**8. `orientation_map(padding_mode=...)` is a silent no-op.** The padded frame is
built and then never used -- the patch is taken from `data[ti, ...]` directly, so
the documented padding behaviour does not happen. Either wire `frame` through or
drop the parameter.

**9. `bidisperse_concentration_map` has no guard on a collapsed search band.**
If the wavelength grid is coarse enough that `min_val_a == max_val_a`, the
`np.argmax` runs on an empty slice and raises; one bin wider and it returns that
bin with no peak-finding at all. This is the regime the accompanying paper found
for `R = 6` at `patchw = 32`.

**10. Off-by-one loop bounds.** `hanning_window` and `angular_binning` start
their loops at index 1 while `radial_grid` starts at 0. In `hanning_window` the
skipped row and column are zeroed by the radial cutoff anyway, and in
`angular_binning` one Monte Carlo sample in 10^4 is dropped, so neither changes a
result -- but the inconsistency is worth removing.
