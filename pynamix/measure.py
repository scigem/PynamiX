import sys
import os
import pynamix
import numpy as np
import warnings
from astropy.convolution import convolve
from scipy.signal import correlate2d
from scipy.stats import linregress
from scipy.ndimage import zoom, gaussian_filter
from pynamix.exposure import *

module_loc = pynamix.__file__[:-11]


def main_direction(tensor):
    """Calculate the principal orientation and orientation magnitude of a nematic order tensor.

    A nematic order tensor is symmetric by construction, so its eigenvalues are
    real.  ``numpy.linalg.eigh`` is used rather than ``numpy.linalg.eig``: the
    general routine returns complex results for matrices it cannot certify as
    symmetric, which then propagates into ``arctan2`` and raises.  The input is
    symmetrised first so that the function is total for any 2 by 2 array.

    Args:
        tensor: 2 by 2 array representing the nematic order tensor.

    Returns:
        Two values, one for the principal orientation in radians (from zero to pi) and one for the magnitude of the orientation on a scale of zero to one.
    """
    tensor = np.asarray(tensor, dtype=float)
    tensor = 0.5 * (tensor + tensor.T)  # nematic tensors are symmetric
    v, V = np.linalg.eigh(tensor)  # eigenvalues (v) and eigenvectors (V), both real
    idx = np.argmax(v)  # principal eigenvalue
    angle = np.arctan2(V[1, idx], V[0, idx])
    if angle < 0:
        angle += np.pi
    elif angle > np.pi:
        angle -= np.pi
    dzeta = np.sqrt(np.sum(tensor**2))  # ||T|| (after eq. 5)
    return angle, dzeta


def hanning_window(patchw=32):
    """Compute a radial Hanning window.

    The window is a raised cosine in radius, unity at the centre and zero at
    ``dst = patchw``, and it is centred on the geometric centre of the
    ``2*patchw`` square, which lies at index ``patchw - 0.5``.

    .. note::
       Before v0.36 the window was built about ``patchw + 0.5``, one pixel off
       in each direction, which left row 0 and column 0 identically zero while
       rows and columns ``2*patchw - 1`` were not.  A radially symmetric window
       stays radially symmetric when it is shifted, so this did not bias the
       orientation measurement, but it did make the support lopsided.

    Args:
        patchw (int): The half width of the patch.

    Returns:
        The radial Hanning window, shape ``(2*patchw, 2*patchw)``.
    """
    i, j = np.indices((2 * patchw, 2 * patchw))
    dst = np.sqrt((i - (patchw - 0.5)) ** 2 + (j - (patchw - 0.5)) ** 2)
    w = 0.5 * np.cos(np.pi * dst / patchw) + 0.5
    w[dst > patchw] = 0.0
    return w


def grid(data, logfile, xstep, ystep, patchw, mode="bottom-left"):
    """Generate two 1D vectors that grid an image

    Creates a grid of patch centers for image analysis. The grid respects
    patch boundaries and optionally ROI constraints.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile: The logfile containing detector information and optional ROI
        xstep (int): The spacing between patches in x direction
        ystep (int): The spacing between patches in y direction
        patchw (int): The half width of the patch
        mode (str): Grid alignment mode:
            - "bottom-left": Align grid to bottom-left corner with patchw buffer
            - "center" or "centre": Center the grid in the available space
            - "full": Cover full domain (no patch buffer, may go outside boundaries)

    Returns:
        tuple: (gridx, gridy) - Two 1D numpy arrays with patch center coordinates

    Examples:
        >>> data = np.zeros((10, 100, 80))
        >>> logfile = {"detector": {}}
        >>> gridx, gridy = grid(data, logfile, 16, 16, 8)
        >>> # gridx starts at 8 (patchw) and ends before 92 (100-8)
    """
    nt, nx, ny = data.shape

    # Determine the effective domain considering ROI
    if "detector" in logfile and "ROI" in logfile["detector"]:
        roi = logfile["detector"]["ROI"]
        # ROI defines the region of interest in the full image
        x_min = roi["left"]
        x_max = roi["right"]
        y_min = roi["top"]
        y_max = roi["bottom"]
    else:
        # No ROI defined, use full domain
        x_min = 0
        x_max = nx
        y_min = 0
        y_max = ny

    # Apply patch buffer constraints based on mode
    if mode in ["bottom-left", "bottom_left"]:
        # Start from bottom-left with patch buffer
        x_start = x_min + patchw
        x_end = x_max - patchw
        y_start = y_min + patchw
        y_end = y_max - patchw

        # Validate boundaries
        if x_start < x_end and 0 <= x_start < nx and 0 < x_end <= nx:
            gridx = np.arange(x_start, x_end, xstep)
        else:
            gridx = np.array([])

        if y_start < y_end and 0 <= y_start < ny and 0 < y_end <= ny:
            gridy = np.arange(y_start, y_end, ystep)
        else:
            gridy = np.array([])

    elif mode in ["center", "centre", "centered", "centred"]:
        # Center the grid in the available space
        x_start = x_min + patchw
        x_end = x_max - patchw
        y_start = y_min + patchw
        y_end = y_max - patchw

        # Calculate number of patches that fit
        if x_start < x_end and y_start < y_end:
            nx_patches = int((x_end - x_start) / xstep)
            ny_patches = int((y_end - y_start) / ystep)

            # Calculate total space used
            x_span = nx_patches * xstep
            y_span = ny_patches * ystep

            # Center the grid
            x_offset = (x_end - x_start - x_span) / 2
            y_offset = (y_end - y_start - y_span) / 2

            if nx_patches > 0:
                gridx = np.arange(nx_patches) * xstep + x_start + x_offset
            else:
                gridx = np.array([])

            if ny_patches > 0:
                gridy = np.arange(ny_patches) * ystep + y_start + y_offset
            else:
                gridy = np.array([])
        else:
            gridx = np.array([])
            gridy = np.array([])

    elif mode == "full":
        # Cover full domain without patch buffer (may go outside for edge patches)
        if 0 <= x_min < nx and 0 < x_max <= nx:
            gridx = np.arange(x_min, x_max, xstep)
        else:
            gridx = np.array([])

        if 0 <= y_min < ny and 0 < y_max <= ny:
            gridy = np.arange(y_min, y_max, ystep)
        else:
            gridy = np.array([])
    else:
        raise ValueError(f"Unknown grid mode: {mode}. Use 'bottom-left', 'center', or 'full'")

    return gridx, gridy


def _subpixel_k(patchw, subsample):
    """Sub-pixel Cartesian offsets for pixel-area quadrature.

    Returns ``(X, Y, r)`` broadcast over ``(row, col, sub_y, sub_x)``, with the
    row index carrying ``y`` and the column index carrying ``x``, matching the
    indexing used throughout this module.
    """
    off = (np.arange(subsample) + 0.5) / subsample - 0.5
    idx = np.arange(2 * patchw) - patchw + 0.5  # pixel centres, in k units
    c = idx[:, None] + off[None, :]
    X = c[None, :, None, :]  # column -> x
    Y = c[:, None, :, None]  # row -> y
    return X, Y, np.sqrt(X**2 + Y**2)


def angular_binning(patchw=32, N=None, subsample=16):
    """
    Compute the individual Q(k) coefficients for equation 4 in `Guillard et al. 2017 <https://www.nature.com/articles/s41598-017-08573-y>`_.

    For each pixel of the shifted power spectrum this is the average of the
    outer product of the unit wavevector with itself, taken over the area of
    that pixel.

    .. note::
       Before v0.36 this was estimated by Monte Carlo, with ~10^7 unseeded
       random throws cached to disk.  That made the coefficients irreproducible
       between users -- two unseeded runs differed by up to 0.012 -- and left
       about 1% noise frozen into every orientation measurement.  Deterministic
       sub-pixel quadrature is exact to machine precision, takes milliseconds
       rather than minutes, and needs no cache.  ``N`` is accepted and ignored.

    Args:
        patchw (int): The half width of the patch.
        N: Ignored. Retained so existing calls keep working.
        subsample (int): Quadrature points per pixel edge.

    Returns:
        4D array: An array n_maskQ that stores Q(k) coefficients in a 4D table such that Q(kx, ky)=n_maskQ(kx,ky,:,:)
    """
    if N is not None:
        warnings.warn(
            "angular_binning is now deterministic; the Monte Carlo sample count N "
            "is ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    X, Y, r = _subpixel_k(patchw, subsample)
    with np.errstate(invalid="ignore", divide="ignore"):
        kx = np.where(r > 0, X / r, 0.0)
        ky = np.where(r > 0, Y / r, 0.0)
    n_maskQ = np.empty([patchw * 2, patchw * 2, 2, 2])
    n_maskQ[:, :, 0, 0] = (kx * kx).mean(axis=(2, 3))
    n_maskQ[:, :, 0, 1] = (kx * ky).mean(axis=(2, 3))
    n_maskQ[:, :, 1, 0] = n_maskQ[:, :, 0, 1]
    n_maskQ[:, :, 1, 1] = (ky * ky).mean(axis=(2, 3))
    return n_maskQ


def radial_grid(rnb=200, patchw=32, N=None, subsample=16):
    """Radial binning weights for the orthoradial summation, and the bin radii.

    Bin membership is by nearest node on ``linspace(0, 1.5*patchw, rnb)``,
    unchanged.  What changed in v0.36 is how the weights and the radii are
    obtained.

    .. note::
       The weights were previously estimated by unseeded Monte Carlo; they are
       now computed by deterministic sub-pixel quadrature.  More importantly,
       ``r_grid`` used to be the node array shifted *up* by half a bin after the
       binning loop, so each bin was labelled by its upper edge rather than its
       centre.  That reported the radius high, and so the wavelength low, by
       roughly half a bin -- about 1% at a peak sitting seven bins out.  Each
       bin is now labelled with the area-weighted mean radius of the region
       actually assigned to it, which is the quantity the label is meant to be.

    Args:
        rnb (int): Number of points to discretise in the radial direction
        patchw (int): The half width of the patch.
        N: Ignored. Retained so existing calls keep working.
        subsample (int): Quadrature points per pixel edge.

    Returns:
        ``(r_grid, nr_pxr)``: the mean radius of each bin, and the fraction of
        each pixel's area falling in each bin.
    """
    if N is not None:
        warnings.warn(
            "radial_grid is now deterministic; the Monte Carlo sample count N is "
            "ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    nodes = np.linspace(0, patchw * 1.5, rnb)
    _, _, r = _subpixel_k(patchw, subsample)
    r = np.broadcast_to(r, (2 * patchw, 2 * patchw, subsample, subsample))
    arg = np.abs(r[..., None] - nodes).argmin(axis=-1)  # nearest node, as before

    nr_pxr = np.zeros([patchw * 2, patchw * 2, rnb])
    rsum = np.zeros(rnb)
    rcnt = np.zeros(rnb)
    flat_r = r.reshape(2 * patchw, 2 * patchw, -1)
    flat_a = arg.reshape(2 * patchw, 2 * patchw, -1)
    npts = flat_a.shape[-1]
    for k in range(rnb):
        hit = flat_a == k
        nr_pxr[:, :, k] = hit.sum(axis=-1) / npts
        rsum[k] = flat_r[hit].sum()
        rcnt[k] = hit.sum()
    # label each bin by the area-weighted mean radius of what landed in it;
    # empty bins fall back to the node itself
    r_grid = np.where(rcnt > 0, rsum / np.maximum(rcnt, 1), nodes)
    return r_grid, nr_pxr


def _check_wavelength_window(lam_max, patchw, resolution, what):
    r"""Warn when the longest wavelength sought is too close to DC to measure.

    A patch ``2*patchw`` across resolves a wavelength ``lambda`` as a peak
    ``n = 2*patchw/lambda`` bins from the origin.  Below about five bins the
    peak sits inside the window's main lobe and on the steep part of the
    ``2*pi*q`` Jacobian carried by the orthoradial sum, and the reported
    wavelength runs 10--17% low however it is read off.  Between roughly
    ``patchw/5`` and ``patchw/3`` the bias is a few per cent.
    """
    n_cycles = 2.0 * patchw / max(lam_max * resolution, 1e-12)
    if n_cycles < 5.0:
        warnings.warn(
            f"{what}: the longest wavelength sought sits only {n_cycles:.1f} Fourier "
            f"bins from DC at patchw={patchw}; below about 5 the reported wavelength "
            f"is biased low by 10% or more. Use a larger patchw, or restrict the "
            f"search to shorter wavelengths.",
            RuntimeWarning,
            stacklevel=3,
        )


def _parabolic_peak(x, y, i):
    """Sub-bin peak abscissa from a three-point parabola, vectorised.

    ``x`` is the (possibly non-uniform) abscissa, ``y`` the spectrum with the
    radial bin last, and ``i`` the index of the largest bin.  Where the three
    points do not bracket a maximum, or the vertex falls outside them, the bin
    abscissa itself is returned.

    Reading the peak as the largest bin quantises the answer onto the available
    wavelength grid, which for a 64-pixel patch near ``lambda = 12`` px offers
    only ``64/5`` and ``64/6``.  That quantisation, not the physics, dominates
    the error in the reported size.
    """
    i = np.clip(i, 1, len(x) - 2)
    x0, x1, x2 = x[i - 1], x[i], x[i + 1]
    y0 = np.take_along_axis(y, (i - 1)[..., None], axis=-1)[..., 0]
    y1 = np.take_along_axis(y, i[..., None], axis=-1)[..., 0]
    y2 = np.take_along_axis(y, (i + 1)[..., None], axis=-1)[..., 0]
    d0 = (x0 - x1) * (x0 - x2)
    d1 = (x1 - x0) * (x1 - x2)
    d2 = (x2 - x0) * (x2 - x1)
    a = y0 / d0 + y1 / d1 + y2 / d2
    b_num = y0 * (x1 + x2) / d0 + y1 * (x0 + x2) / d1 + y2 * (x0 + x1) / d2
    with np.errstate(invalid="ignore", divide="ignore"):
        vertex = b_num / (2.0 * a)
    lo, hi = np.minimum(x0, x2), np.maximum(x0, x2)
    ok = (a < 0) & np.isfinite(vertex) & (vertex > lo) & (vertex < hi)
    return np.where(ok, vertex, x1)


def orientation_map(
    data,
    logfile,
    tmin=0,
    tmax=None,
    tstep=1,
    xstep=32,
    ystep=32,
    patchw=32,
    normalisation=mean_std,
    padding_mode=None,
    verbose=0,
):
    """
    Calculate the principal orientation and orientation magnitude at a set of patches in images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile: The logfile.
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        patchw (int): The half width of the patch.
        normalisation: Which normalisation to use.
        padding_mode: What to use for data points outside the experimental domain. Pick one from the list from the documentation of `numpy.pad` or None (default)
        verbose: (int): How noisy this function should be. 0 for nothing (default), 1 for text updates, 2 for figures.

    Returns:
        Four 2D arrays which describe: (1) the x location of the centre of each patch, (2) the y location of the centre of each patch, (3) the principal orientation and (4) the orientation magnitude for each patch. By default, pads with empty numbers to return NaNs.
    """
    w = hanning_window(patchw)
    n_maskQ = angular_binning(patchw)
    nt, nx, ny = data.shape

    if tmax is None:
        tmax = nt  # optionally set end time
    gridx, gridy = grid(data, logfile, xstep, ystep, patchw, mode="bottom-left")
    # Prepare thre result matrices (3D), first 2 indices are the grid, the last index is time
    orient = np.nan * np.zeros([tmax - tmin, len(gridx), len(gridy)])
    dzeta = np.nan * np.zeros([tmax - tmin, len(gridx), len(gridy)])
    Q = np.zeros_like(n_maskQ)
    # Q2 = np.zeros([2, 2, nx, ny, nt])
    Q2 = np.zeros([2, 2])
    # Loop on the movie frames
    for t, ti in enumerate(range(tmin, tmax, tstep)):
        if verbose == 1:
            print("Up to frame " + str(ti) + " out of " + str(tmax), end="\r")
        frame = data[ti]
        if padding_mode is not None:
            # the padded frame was previously built and then never used, so
            # padding_mode was a silent no-op; the offset keeps patch centres
            # aligned with the unpadded grid
            frame = np.pad(frame, patchw, mode=padding_mode)
            shift = patchw
        else:
            shift = 0
        for i, xi in enumerate(gridx):  # Loop over the grid
            for j, yj in enumerate(gridy):
                patch = frame[
                    xi + shift - patchw : xi + shift + patchw,
                    yj + shift - patchw : yj + shift + patchw,
                ]

                patch = normalisation(patch)

                S = np.fft.fftshift(np.abs(np.fft.fft2(patch * w) ** 2))
                if np.sum(S) == 0:
                    Q[:, :, 0, 0] = Q[:, :, 0, 1] = Q[:, :, 1, 0] = Q[:, :, 1, 1] = S
                else:
                    Q[:, :, 0, 0] = n_maskQ[:, :, 0, 0] * S / np.sum(S)
                    Q[:, :, 0, 1] = n_maskQ[:, :, 0, 1] * S / np.sum(S)
                    Q[:, :, 1, 0] = n_maskQ[:, :, 1, 0] * S / np.sum(S)
                    Q[:, :, 1, 1] = n_maskQ[:, :, 1, 1] * S / np.sum(S)
                Q2 = np.sum(np.sum(Q, axis=0), axis=0)
                Q2[0, 0] -= 0.5
                Q2[1, 1] -= 0.5
                Q2 *= np.sqrt(2)
                orient[t, i, j], dzeta[t, i, j] = main_direction(Q2)
    X, Y = np.meshgrid(gridx, gridy, indexing="ij")
    return X, Y, orient, dzeta


def radial_FFT(
    data,
    logfile,
    rnb=200,
    tmin=0,
    tmax=None,
    tstep=1,
    xstep=32,
    ystep=32,
    patchw=32,
    normalisation=mean_std,
):
    """
    Calculate the orthoradial sum of the 2D FFT at a set of patches in images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile (dict): The logfile.
        rnb (int): Number of points to discretise in the radial direction
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        patchw (int): The half width of the patch.

    Returns:
        Four arrays which describe: (1) a 1D vector with the x location of each patch, (2) a 1D vector with the y location of each patch, (3) a 1D vector with the wavelengths corresponding to each spectral dnesity and (4) a 4D array representing the orthoradially summed power spectral density for each patch.
    """
    w = hanning_window(patchw)
    n_maskQ = angular_binning(patchw)

    nt, nx, ny = data.shape

    gridx = np.arange(patchw, nx - patchw, xstep)  # locations of centres of patches in x direction
    gridy = np.arange(patchw, ny - patchw, ystep)  # locations of centres of patches in y direction
    if tmax is None:
        tmax = nt  # optionally set end time

    frequencyconversion = logfile["detector"]["resolution"] / (
        patchw * 2
    )  # #(do 1/(frequencyconversion*peakfreq) to get the spatial caracteristic wavelength)
    radialspec = np.zeros([(tmax - tmin) // tstep, len(gridx), len(gridy), rnb])  #
    # radialspec = np.zeros([rnb,len(gridx)*len(gridy)*(tmax-tmin)//tstep]) # JUST DURING TESTING

    r_grid, nr_pxr = radial_grid(rnb=rnb, patchw=patchw)
    wavelength = 1.0 / (r_grid * frequencyconversion)  # wavelength in mm
    n = 0
    for t, ti in enumerate(range(tmin, tmax, tstep)):  # Loop on the movie frames
        for i, xi in enumerate(gridx):  # Loop over the grid
            for j, yj in enumerate(gridy):
                patch = data[ti, xi - patchw : xi + patchw, yj - patchw : yj + patchw]
                patch = normalisation(patch)
                S = np.fft.fftshift(np.abs(np.fft.fft2(patch * w) ** 2))

                for k in range(rnb):
                    radialspec[t, i, j, k] = np.sum(S * nr_pxr[:, :, k])  # ortho-radially SUMMED power spectral density
                    # radialspec[k,n]=np.sum(S*nr_pxr[:,:,k])  # ortho-radially SUMMED power spectral density - JUST FOR PLOTTING DURING TESTING
                # n += 1

    return gridx, gridy, wavelength, radialspec


def average_size_map(
    data,
    logfile,
    rnb=200,
    tmin=0,
    tmax=None,
    tstep=1,
    xstep=32,
    ystep=32,
    patchw=32,
    normalisation=mean_std,
    wmin=None,
    wmax=10,
    return_FFTs=False,
    interpolate=True,
):
    """
    Calculate the radial average of the 2D FFT at a set of patches in images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile (dict): The logfile.
        rnb (int): Number of points to discretise in the radial direction
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        patchw (int): The half width of the patch.
        normalisation (function): Which function to normalise the patches by
        wmin (float): Minimum wavelength (mm)
        wmax (float): Maximum wavelength (mm)
        return_FFTs (bool): Optionally return the full FFTs corresponding to each patch
        interpolate (bool): Fit a parabola through the peak bin and its two
            neighbours rather than reporting the peak bin itself. Default True.
            Set False to recover the pre-v0.36 behaviour.

    Returns:
        Three 2D arrays which describe: (1) the x location of the centre of each patch, (2) the y location of the centre of each patch and (3) the average size for each patch, defined as the wavelength corresponding to the highest peak in the orthoradially summed power spectral density.
    """
    gridx, gridy, wavelength, radialspec = radial_FFT(
        data,
        logfile,
        rnb=rnb,
        tmin=tmin,
        tmax=tmax,
        tstep=tstep,
        xstep=xstep,
        ystep=ystep,
        patchw=patchw,
        normalisation=normalisation,
    )

    if wmin is None:
        wmin = 2 / logfile["detector"]["resolution"]  # use Nyquist frequency - i.e. 2 pixels per particle
    _check_wavelength_window(wmax, patchw, logfile["detector"]["resolution"],
                             "average_size_map")
    min_val = np.argmin(np.abs(wavelength - wmax))  # this is large wavelength, wavelength is sorted large to small
    max_val = np.argmin(np.abs(wavelength - wmin))  # this is small wavelength, wavelength is sorted large to small
    # print(wavelength[min_val],wavelength[max_val])
    if max_val <= min_val + 1:
        raise ValueError(
            f"the search window {wmin}--{wmax} spans {max_val - min_val} wavelength "
            f"bins; widen it, raise rnb, or use a larger patchw"
        )
    average_size_index = np.argmax(radialspec[:, :, :, min_val:max_val], axis=3)
    average_size_index += min_val

    if interpolate:
        # interpolate in 1/wavelength, which is proportional to the bin radius
        # and so is the variable the spectrum is actually smooth in
        inv = 1.0 / wavelength
        inv_peak = _parabolic_peak(inv, radialspec, average_size_index)
        with np.errstate(divide="ignore"):
            size = 1.0 / inv_peak
    else:
        size = wavelength[average_size_index]

    X, Y = np.meshgrid(gridx, gridy, indexing="ij")

    if return_FFTs:
        return X, Y, size, wavelength, radialspec
    else:
        return X, Y, size


def bidisperse_concentration_map(
    data,
    logfile,
    rnb=200,
    tmin=0,
    tmax=None,
    tstep=1,
    xstep=32,
    ystep=32,
    patchw=32,
    normalisation=mean_std,
    s_a=0.5,
    s_b=2.0,
    pad=1.2,
    return_FFTs=False,
    calib_func=None,
):
    """
    Calculate the concentration of species a of a bidisperse mixture at a set of patches in images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile (dict): The logfile.
        rnb (int): Number of points to discretise in the radial direction
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        patchw (int): The half width of the patch.
        normalisation (function): Which function to normalise the patches by
        s_a (float): Size of one set of particles - will return the concentration of this size. Doesn't matter if it is large or small. (mm)
        s_b (float): Size of other particles (mm)
        pad (float): Range to look for the peak around each size
        calib_func (function): Function to use to calibrate the concentration from the peak fraction. If no function given, the peak fraction is returned.
        return_FFTs (bool): Optionally return the full FFTs corresponding to each patch

    Returns:
        Three 2D arrays which describe: (1) the x location of the centre of each patch, (2) the y location of the centre of each patch and (3) the average size for each patch, defined as the wavelength corresponding to the highest peak in the orthoradially summed power spectral density.
    """
    gridx, gridy, wavelength, radialspec = radial_FFT(
        data,
        logfile,
        rnb=rnb,
        tmin=tmin,
        tmax=tmax,
        tstep=tstep,
        xstep=xstep,
        ystep=ystep,
        patchw=patchw,
        normalisation=normalisation,
    )

    _check_wavelength_window(max(s_a, s_b) * pad, patchw,
                             logfile["detector"]["resolution"],
                             "bidisperse_concentration_map")
    min_val_a = np.argmin(
        np.abs(wavelength - s_a * pad)
    )  # this is large wavelength, wavelength is sorted large to small
    max_val_a = np.argmin(
        np.abs(wavelength - s_a / pad)
    )  # this is small wavelength, wavelength is sorted large to small
    min_val_b = np.argmin(
        np.abs(wavelength - s_b * pad)
    )  # this is large wavelength, wavelength is sorted large to small
    max_val_b = np.argmin(
        np.abs(wavelength - s_b / pad)
    )  # this is small wavelength, wavelength is sorted large to small

    for name, lo, hi in (("a", min_val_a, max_val_a), ("b", min_val_b, max_val_b)):
        if hi <= lo + 1:
            raise ValueError(
                f"the search window for species {name} spans {hi - lo} wavelength "
                f"bins, which is too few to locate a peak in. Widen pad, raise rnb, "
                f"or use a larger patchw: at patchw={patchw} the available "
                f"wavelengths near a peak {2 * patchw // 8} px across are spaced by "
                f"roughly a tenth of the wavelength itself."
            )
    peak_a_index = np.argmax(radialspec[:, :, :, min_val_a:max_val_a], axis=3)
    peak_b_index = np.argmax(radialspec[:, :, :, min_val_b:max_val_b], axis=3)
    peak_a_index += min_val_a
    peak_b_index += min_val_b

    peak_a = np.take_along_axis(radialspec, peak_a_index[..., None], axis=-1)[..., 0]
    peak_b = np.take_along_axis(radialspec, peak_b_index[..., None], axis=-1)[..., 0]

    peak_fraction = peak_a / (peak_a + peak_b)

    if calib_func is not None:
        concentration = calib_func(peak_fraction)
    else:
        concentration = peak_fraction

    X, Y = np.meshgrid(gridx, gridy, indexing="ij")

    if return_FFTs:
        return X, Y, concentration, wavelength, radialspec
    else:
        return X, Y, concentration


def velocity_map(
    data,
    logfile,
    tmin=0,
    tmax=None,
    tstep=1,
    searchw=4,
    xstep=32,
    ystep=32,
    patchw=32,
    zoomvalue=10,
    normalisation=mean_std,
    padding_mode="edge",
    verbose=0,
):
    """
    Calculate the velocity field at a set of patches between two images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile: The logfile.
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series (velocities at this frame will not be calculated)
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        searchw (int): Half width of search area
        patchw (int): The half width of the patch.
        zoom (int): Zoom level for subpixel resolution.
        normalisation: Which normalisation to use.
        padding_mode: What to use for data points outside the experimental domain. Pick one from "fill", "wrap" or "symm"
        verbose (int): 0 for nothing. 1 for text. 2 for graphs.

    .. warning:: THIS WORKS VERY POORLY AND SLOWLY - NOTE SURE WHY!!

    Returns:
        Four 2D arrays which describe: (1) the x location of the centre of each patch, (2) the y location of the centre of each patch, (3) the x velocity (4) the y velocity. By default, pads with empty numbers to return NaNs.
    """
    print("WARNING: Do not use this function it is pretty garbage.")
    nt, nx, ny = data.shape
    px_per_frame_to_mm_per_s = tstep * logfile["detector"]["resolution"] / logfile["detector"]["fps"]

    if tmax is None:
        tmax = nt  # optionally set end time
    # Prepare thre result matrices
    gridx, gridy = grid(data, logfile, xstep, ystep, patchw, mode="bottom-left")
    u = np.nan * np.zeros([tmax - tmin, len(gridx), len(gridy)])
    v = np.nan * np.zeros([tmax - tmin, len(gridx), len(gridy)])

    # Loop over the frames with a given spacing tstep
    for t, ti in enumerate(range(tmin, tmax - 1, tstep)):
        if verbose == 1:
            print("Up to frame " + str(ti) + " out of " + str(tmax), end="\r")
        this_frame = data[ti]  # .copy()
        next_frame = data[ti + tstep]  # .copy()

        patchws = [patchw]  # * 4, patchw * 2, patchw]
        xsteps = [xstep]  # * 4, xstep * 2, xstep]
        ysteps = [ystep]  # * 4, ystep * 2, ystep]
        searchws = [searchw]  # , searchw // 2, searchw // 4]

        for level, pw in enumerate(patchws):
            sw = searchws[level]

            this_frame = np.pad(this_frame, pw, mode=padding_mode)
            next_frame = np.pad(next_frame, pw + sw, mode=padding_mode)

            if zoom != 1:  # how does this work with levels????
                # Linearly interpolate images
                this_frame = zoom(this_frame, zoomvalue, order=1)
                next_frame = zoom(next_frame, zoomvalue, order=1)

            # Loop over the grid
            # xi,yj are the centre of the patch but NOT the search window,
            # which is offset by searchw in both directions
            for i, xi in enumerate(gridx):
                for j, yj in enumerate(gridy):
                    this_patch = this_frame[
                        (xi - pw) * zoomvalue : (xi + pw + 1) * zoomvalue,
                        (yj - pw) * zoomvalue : (yj + pw + 1) * zoomvalue,
                    ]
                    search_window = next_frame[
                        (xi - pw) * zoomvalue : (xi + pw + 2 * sw + 1) * zoomvalue,
                        (yj - pw) * zoomvalue : (yj + pw + 2 * sw + 1) * zoomvalue,
                    ]

                    # this_patch = this_patch - np.mean(this_patch)
                    # search_window = search_window - np.mean(search_window)

                    # this_patch = (this_patch - np.mean(this_patch)) / np.std(this_patch)
                    # search_window = (search_window - np.mean(search_window)) / np.std(
                    # search_window
                    # )

                    corr1 = correlate2d(search_window, this_patch, mode="valid")
                    # corr2 = np.sum(correlate2d(this_patch, this_patch, mode="valid"))
                    # corr3 = np.sum(
                    # correlate2d(search_window, search_window, mode="valid")
                    # )
                    C = corr1
                    # C = corr1 / corr2 / corr3
                    # C = convolve(search_window, this_patch)
                    ybest, xbest = np.unravel_index(np.argmax(C), C.shape)
                    if C[ybest, xbest] > 0:  # pick a random threshold for badness
                        u[t, i, j] = (ybest - searchw * zoomvalue) * px_per_frame_to_mm_per_s
                        v[t, i, j] = -(xbest - searchw * zoomvalue) * px_per_frame_to_mm_per_s
                    if verbose == 2:
                        plt.ion()
                        plt.clf()
                        plt.subplot(221)
                        plt.imshow(
                            search_window,
                            vmin=np.amin(search_window),
                            vmax=np.amax(search_window),
                        )
                        plt.colorbar()
                        plt.plot(xbest + patchw, ybest + patchw, "ro")
                        plt.subplot(222)
                        plt.imshow(
                            this_patch,
                            vmin=np.amin(search_window),
                            vmax=np.amax(search_window),
                        )
                        plt.colorbar()
                        plt.subplot(223)
                        plt.imshow(C)
                        plt.plot(xbest, ybest, "ro")
                        plt.colorbar()
                        plt.pause(0.1)

    X, Y = np.meshgrid(gridx, gridy, indexing="ij")
    return X, Y, u, v


def validate_velocities(u, v, threshold=2, sigma=1):
    """
    Validate a velocity field so that it is within a tolerance of the gaussian smoothed value

    Args:
        u (nd array): horizontal velocity field
        v (nd array): vertical velocity field
        threshold (float): How many times larger/smaller a value should be for it to be masked
        sigma (float): Width of gaussian blurring to define mean field

    Returns:
        Masked u and v arrays.
    """
    u_smoothed = gaussian_filter(u, sigma=sigma)
    v_smoothed = gaussian_filter(v, sigma=sigma)
    mask = (np.abs(u - u_smoothed) / abs(u_smoothed) > threshold) | (
        np.abs(v - v_smoothed) / abs(v_smoothed) > threshold
    )
    u_filtered = np.ma.masked_where(mask, u)
    v_filtered = np.ma.masked_where(mask, v)
    return u_filtered, v_filtered


# ======================= Andres ================================================
def surface_elevation_sideview(data, scale, min_height=-1, threshold=-1):
    """
    Detect the free surface by finding the highest point below the threshold accross the width of the data provided. Gravity should be in the direction of increasing first dimension.

    Args:
        data: The source data. Should be in the shape [ny,nx]
        scale: mm/px conversion.
        min_height: line from which to compute elevation. There should only be material above that line.
        threshold: if value < threshold, we are inside the material.
    Returns:
        Array of height vs. horizontal location.
    """
    res = np.zeros((2, data.shape[1]))
    res[0, :] = [i * scale for i in range(0, data.shape[1])]
    if min_height == -1:
        min_height = data.shape[0]
    for i in range(0, data.shape[1]):
        for j in range(min_height, -1, -1):
            if data[j][i] > threshold:
                res[1, i] = (min_height - j) * scale
                break
    return res


def surface_absorption_calibration(free_surface_profile, top_absorption):
    """
    Fit the absorption for a set of thicknesses and absorption values.

    Args:
        free_surface_profile: thickness of the absorbing layer
        top_absorption: corresponding transmitted xray intensity (not logged)
    """
    res = linregress(free_surface_profile, np.log(top_absorption))
    return {"mu": res.slope, "beta": res.intercept}


def surface_elevation(data, absorption, sigma=0):
    """
    Returns the elevation calculated from the xray provided and the fitted absorption coefficients.

    Args:
        data: the xray radiograph
        absorption: dictionary of absorption values as returned by `surface_absorption_calibration`
        sigma: size of the gaussian filter to post-apply
    """
    ele = (1 / absorption["mu"]) * data - (absorption["beta"] / absorption["mu"])
    if sigma != 0:
        ele = gaussian_filter(ele, 6)
    return ele


# Testing area
if __name__ == "__main__":
    # Try with fake data - WORKS
    logfile = {}
    logfile["detector"] = {}
    logfile["detector"]["resolution"] = 1
    logfile["detector"]["fps"] = 1
    data = np.random.rand(2, 200, 300)
    X, Y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[2]), indexing="ij")
    data[0] += X / 50
    data[1] = data[0]
    data[1] = np.roll(data[1], 1, axis=0)
    data[1] = np.roll(data[1], 3, axis=1)

    # data = np.zeros([2, 200, 300])
    # data[0, 30:40, 40:50] = 1
    # data[1] = data[0]
    # data[1] = np.roll(data[1], 1, axis=0)
    # data[1] = np.roll(data[1], 3, axis=1)

    # Try some real data
    # data, logfile = load_seq("/Volumes/LTS/DynamiX/FRC/Marta/MartaTest1-D0")
    # data = data[1000:1002].copy().astype(np.float64)
    # data[1] = data[0]
    # print(data.shape)

    x, y, u, v = velocity_map(
        data,
        logfile,
        # tmin=1000,
        # tmax=1001,
        tstep=1,
        xstep=32,
        ystep=32,
        patchw=32,
        searchw=5,
        zoomvalue=1,
        normalisation=mean_std,
        padding_mode="constant",
        verbose=2,
    )
    print(u)
    # u2, v2 = validate_velocities(u, v, threshold=2, sigma=1)
    # u2 = u2.filled(0)
    # v2 = v2.filled(0)

    plt.clf()

    plt.pcolormesh(
        X,
        Y,
        data[0],
    )
    plt.quiver(x, y, u[0], v[0], color="g")
    # plt.quiver(x, y, u2[0], v2[0], color="k")
    plt.savefig("tt.png", dpi=200)
