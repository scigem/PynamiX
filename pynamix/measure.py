import os, pynamix
import numpy as np
from astropy.convolution import convolve
from scipy.signal import correlate2d
from scipy.stats import linregress
from scipy.ndimage import zoom, gaussian_filter
from pynamix.exposure import *

module_loc = pynamix.__file__[:-11]


def main_direction(tensor):
    """Calculate the principal orientation and orientation magnitude of a nematic order tensor.

    Args:
        tensor: 2 by 2 array representing the nematic order tensor.

    Returns:
        Two values, one for the principal orientation in radians (from zero to pi) and one for the magnitude of the orientation on a scale of zero to one.
    """
    v, V = np.linalg.eig(tensor)  # eigenvalues (v) and eigenvectors (V)
    idx = np.argmax(v)  # principal eigenvalue
    angle = np.arctan2(V[1, idx], V[0, idx])
    if angle < 0:
        angle += np.pi
    elif angle > np.pi:
        angle -= np.pi
    dzeta = np.sqrt(np.sum(tensor ** 2))  # ||T|| (after eq. 5)
    return angle, dzeta


def hanning_window(patchw=32):
    """Compute a radial hanning window.

    Args:
        patchw (int): The half width of the patch.

    Returns:
        The radial hanning window.
    """
    w = np.zeros([patchw * 2, patchw * 2])
    for i in range(1, patchw * 2):
        for j in range(1, patchw * 2):
            dst = np.sqrt((i - 0.5 - patchw) ** 2 + (j - 0.5 - patchw) ** 2)
            w[i, j] = 0.5 * (np.cos(2 * np.pi * dst / (patchw * 2))) + 0.5
            if dst > patchw:
                w[i, j] = 0
    return w


def grid(data, logfile, xstep, ystep, patchw, mode="bottom-left"):
    """Generate two 1D vectors that grid an image

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        logfile: The logfile.
        xstep (int): The half width of the patch.
        ystep (int): The half width of the patch.
        patchw (int): The half width of the patch.
        mode (str): bottom-left or centre. Either bottom/left aligned or centred in the image.

    .. warning:: centred mode not implemented yet

    Returns:
        The radial hanning window.
    """
    nt, nx, ny = data.shape

    if mode == "bottom-left":
        if "ROI" in logfile["detector"]:
            gridx = np.arange(
                logfile["detector"]["ROI"]["left"] + patchw,
                nx - patchw + logfile["detector"]["ROI"]["right"],
                xstep,
            )
            # locations of centres of patches in y direction
            gridy = np.arange(
                logfile["detector"]["ROI"]["bot"] + patchw,
                ny - patchw + logfile["detector"]["ROI"]["top"],
                ystep,
            )
        else:
            gridx = np.arange(patchw, nx - patchw, xstep)
            gridy = np.arange(patchw, ny - patchw, ystep)
    # else:
    # locations of centres of patches in x direction
    # gridx = np.arange(0, nx, xstep)
    # locations of centres of patches in y direction
    # gridy = np.arange(0, ny, ystep)

    else:
        sys.exit("Sorry, haven't implemeneted centred grid yet")
    return gridx, gridy


def angular_binning(patchw=32, N=10000):
    """
    Use a Monte-Carlo method to compute the individual Q(k) coefficients for equation 4 in `Guillard et al. 2017 <https://www.nature.com/articles/s41598-017-08573-y>`_.

    Angular binning definition

    step_angle=5 / 180 * pi

    bin_angle=range(0 - step_angle / 2,- 180,- step_angle)

    n_angle=floor(pi / step_angle)

    .. warning:: Fix this definition.

    Args:
        patchw (int): The half width of the patch.
        N (int): Number of particle throws per iteration.

    Returns:
        4D array: An array n_maskQ that stores Q(k) coefficients in a 4D table such that Q(kx, ky)=n_maskQ(kx,ky,:,:)
    """
    if os.path.exists(module_loc + "defaults/n_maskQ_" + str(patchw) + "_" + str(N) + ".npy"):
        n_maskQ = np.load(module_loc + "defaults/n_maskQ_" + str(patchw) + "_" + str(N) + ".npy")
    else:
        print("WARNING: Haven't cached these Q coefficients. Run will take longer this time.")
        K = np.zeros([N, 2])
        n_nbmaskQ = np.zeros([patchw * 2, patchw * 2])
        n_maskQ = np.zeros([patchw * 2, patchw * 2, 2, 2])
        for j in range(1, 1000):  # Number of iteration in Monte-Carlo simulation
            r = (np.random.rand(N, 2) - 0.5) * 2 * patchw
            K[:, 0] = r[:, 0] / (np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2))
            K[:, 1] = r[:, 1] / (np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2))
            x = np.floor(r + patchw).astype(np.int)
            for i in range(1, N):
                n_nbmaskQ[x[i, 1], x[i, 0]] += 1
                n_maskQ[x[i, 1], x[i, 0], 0, 0] += K[i, 0] * K[i, 0]
                n_maskQ[x[i, 1], x[i, 0], 0, 1] += K[i, 0] * K[i, 1]
                n_maskQ[x[i, 1], x[i, 0], 1, 0] += K[i, 1] * K[i, 0]
                n_maskQ[x[i, 1], x[i, 0], 1, 1] += K[i, 1] * K[i, 1]

        # Final scaling for the n_maskQ coefficients.
        n_maskQ[:, :, 0, 0] /= n_nbmaskQ
        n_maskQ[:, :, 0, 1] /= n_nbmaskQ
        n_maskQ[:, :, 1, 0] /= n_nbmaskQ
        n_maskQ[:, :, 1, 1] /= n_nbmaskQ
        np.save(
            module_loc + "defaults/n_maskQ_" + str(patchw) + "_" + str(N) + ".npy",
            n_maskQ,
        )
        print("Successfully cached Q coefficients.")
    return n_maskQ


def radial_grid(rnb=200, patchw=32, N=10000):  # validated against MATLAB code
    if os.path.exists(module_loc + "defaults/r_grid_" + str(rnb) + "_" + str(patchw) + "_" + str(N) + ".npy"):
        r_grid = np.load(module_loc + "defaults/r_grid_" + str(rnb) + "_" + str(patchw) + "_" + str(N) + ".npy")
        nr_pxr = np.load(module_loc + "defaults/nr_pxr_" + str(rnb) + "_" + str(patchw) + "_" + str(N) + ".npy")
    else:
        r_grid = np.linspace(0, patchw * 1.5, rnb)
        nr_px = np.zeros([patchw * 2, patchw * 2])
        nr_pxr = np.zeros([patchw * 2, patchw * 2, rnb])

        Niter = patchw * 2 * patchw * 2 * 200 // N
        for j in range(Niter):
            r = (np.random.rand(N, 2) - 0.5) * 2 * patchw
            dst = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
            x = np.floor(r + patchw).astype(int)
            for i in range(N):
                nr_px[x[i, 1], x[i, 0]] += 1
                arg = np.argmin(np.abs(r_grid - dst[i]))
                nr_pxr[x[i, 1], x[i, 0], arg] += 1
        for i in range(rnb):
            nr_pxr[:, :, i] /= nr_px
        r_grid += np.mean(np.diff(r_grid)) * 0.5
        np.save(
            module_loc + "defaults/r_grid_" + str(rnb) + "_" + str(patchw) + "_" + str(N) + ".npy",
            r_grid,
        )
        np.save(
            module_loc + "defaults/nr_pxr_" + str(rnb) + "_" + str(patchw) + "_" + str(N) + ".npy",
            nr_pxr,
        )
    return r_grid, nr_pxr


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
            frame = np.pad(frame, patchw, mode=padding_mode)
        for i, xi in enumerate(gridx):  # Loop over the grid
            for j, yj in enumerate(gridy):
                patch = data[ti, xi - patchw : xi + patchw, yj - patchw : yj + patchw]

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

    if wmin == None:
        wmin = 2 / logfile["detector"]["resolution"]  # use Nyquist frequency - i.e. 2 pixels per particle
    min_val = np.argmin(np.abs(wavelength - wmax))  # this is large wavelength, wavelength is sorted large to small
    max_val = np.argmin(np.abs(wavelength - wmin))  # this is small wavelength, wavelength is sorted large to small
    # print(wavelength[min_val],wavelength[max_val])
    average_size_index = np.argmax(radialspec[:, :, :, min_val:max_val], axis=3)
    average_size_index += min_val
    nt, nx, ny, _ = radialspec.shape

    # There must be a better way to do this...
    size = np.zeros([nt, nx, ny])
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
                size[t, x, y] = wavelength[average_size_index[t, x, y]]

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

    peak_a_index = np.argmax(radialspec[:, :, :, min_val_a:max_val_a], axis=3)
    peak_b_index = np.argmax(radialspec[:, :, :, min_val_b:max_val_b], axis=3)
    peak_a_index += min_val_a
    peak_b_index += min_val_b

    nt, nx, ny, _ = radialspec.shape

    # There must be a better way to do this...
    peak_a = np.zeros([nt, nx, ny])
    peak_b = np.zeros([nt, nx, ny])
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
                #             print(radialspec[:,:,:,peak_a_index[t,x,y]])
                peak_a[t, x, y] = radialspec[t, x, y, peak_a_index[t, x, y]]
                peak_b[t, x, y] = radialspec[t, x, y, peak_b_index[t, x, y]]

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

            if zoom is not 1:  # how does this work with levels????
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

#======================= Andres ================================================
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
    res = np.zeros((2,data.shape[1])) ;
    res[0,:] = [i*scale for i in range(0,data.shape[1])] ; 
    if min_height==-1: 
        min_height = data.shape[0] 
    for i in range(0,data.shape[1]):
        for j in range(min_height, -1, -1):
            if data[j][i]>threshold :
                res[1,i]= (min_height-j) * scale; 
                break ; 
    return res ; 

def surface_absorption_calibration (free_surface_profile, top_absorption):
    """
    Fit the absorption for a set of thicknesses and absorption values. 
    
    Args:
        free_surface_profile: thickness of the absorbing layer
        top_absorption: corresponding transmitted xray intensity (not logged)
    """
    res=linregress(free_surface_profile, np.log(top_absorption)) ; 
    return {'mu':res.slope, 'beta': res.intercept} ;

def surface_elevation (data, absorption, sigma=0):
    """
    Returns the elevation calculated from the xray provided and the fitted absorption coefficients. 
    
    Args:
        data: the xray radiograph
        absorption: dictionary of absorption values as returned by `surface_absorption_calibration`
        sigma: size of the gaussian filter to post-apply
    """
    ele = (1/absorption['mu'])*data- (absorption['beta']/absorption['mu']);
    if sigma!=0:
        ele=gaussian_filter(ele,6)
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
