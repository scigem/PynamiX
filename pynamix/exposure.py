import numpy as np
from progressbar import progressbar

def mean_std(im):
    """
    Normalise an image such that it has zero mean and standard deviation of one.

    Args:
        im (2D array): The image to normalise.

    Returns:
        out (2D array): The normalised image.
    """
    std = np.std(im)
    if std != 0:
        out = (im - np.mean(im)) / std
    else:
        out = im - np.mean(im)
    return out


def no_normalisation(im):
    """
    Do not normalise the input image.

    Args:
        im (2D array): The image to not normalise.

    Returns:
        out (2D array): The unnormalised image.
    """
    return im


def clamp(data, vmin, vmax):
    """
    Clamp an image between two values.

    Args:
        data (ND array): The image to clamp.
        vmin (float): The lowest value to clamp below.
        vmax (float): The highest value to clamp above.

    Returns:
        clamped_data (ND array): The same data as the original, but with no values outside of the defined range.
    """
    # masked_data = np.ma.masked_outside(data, vmin, vmax, copy=True)
    clamped_data = data.copy()
    clamped_data[data < vmin] = vmin
    clamped_data[data > vmax] = vmax
    return clamped_data


def apply_ROI(data, logfile, top=0, left=0, right=None, bottom=None):
    """
    Apply an ROI to an image or a series of images.

    Args:
        data (2D or 3D array): The source data.
        top (int): the index of the top edge of the ROI.
        left (int): the index of the left edge of the ROI.
        right (int): the index of the right edge of the ROI.
        bottom (int): the index of the bottom edge of the ROI.
    """
    N = len(data.shape)  # number of dimensions
    if N == 2:
        nx, ny = data.shape
        if right == None:
            right = nx
        if bottom == None:
            bottom = ny
        data_ROI = data[left:right, top:bottom]
    elif N == 3:
        _, nx, ny = data.shape
        if right == None:
            right = nx
        if bottom == None:
            bottom = ny
        data_ROI = data[:, left:right, top:bottom]
    else:
        raise Exception("ROI only defined for 2D and 3D arrays")
    if not "detector" in logfile:
        logfile["detector"] = {}
    logfile["detector"]["ROI_software"] = {}
    logfile["detector"]["ROI_software"]["top"] = top
    logfile["detector"]["ROI_software"]["left"] = left
    logfile["detector"]["ROI_software"]["right"] = right
    logfile["detector"]["ROI_software"]["bottom"] = bottom

    return data_ROI, logfile

def set_motion_limits(data, logfile, threshold=False, verbose=False):
    """
    Look at the changes in pixel values over time and try to find the start and end frames in the time series. Save these values into the logfile as `start_frame` and `end_frame`.

    Args:
        data (ND array): The image to normalise.
        logfile: The logfile for the source data.
        threshold (optional): The threshold in image pixel values squared to separate motion from no motion.

    Returns:
        logfile: Updated logfile with new motion limits `start_frame` and `end_frame`.
    """

    # rel_diff = np.nan_to_num((data[1:,10:-10,10:-10] - data[:-1,10:-10,10:-10])/data[:-1,10:-10,10:-10])
    # diff = np.sqrt(np.mean(np.mean(np.square(rel_diff),axis=-1),axis=-1))
    diff = np.sqrt(np.mean(np.mean(np.square(data[1:] - data[:-1]),axis=-1),axis=-1))

    if threshold == False:
        alpha = 0.9 # skew towards the lower end of the spectrum
        threshold = (1-alpha)*diff.max() + alpha*diff.min()
    moving = diff > threshold

    logfile['start_frame'] = int(np.nonzero(moving)[0][0] - 1)
    logfile['end_frame']   = int(np.nonzero(moving)[0][-1])

    if verbose:
        import matplotlib.pyplot as plt
        print(f'Total number of frames: {data.shape[0]}')
        print(f"Start frame: {logfile['start_frame']}. End frame: {logfile['end_frame']}")

        plt.plot(diff,'k')
        plt.plot(np.arange(len(moving)),threshold*np.ones_like(moving),'r--')
        plt.plot([logfile['start_frame'],logfile['start_frame']],[diff.min(),diff.max()],'g-')
        plt.plot([logfile['end_frame'],logfile['end_frame']],[diff.min(),diff.max()],'g-')

    return logfile


def set_angles_from_limits(logfile, max_angle=360):
    """
    Sets the angle in degrees for each frame during one (or multiple) rotations

    Args:
        logfile: The logfile for the source data.
        max_angle (optional): The final rotation at the end of the time series. Default is 360.

    Returns:
        logfile: Updated logfile with an `angles` field.
    """

    num_frames = len(logfile['detector']['frames'])
    angles = np.nan*np.ones(num_frames)
    angles[logfile['start_frame']:logfile['end_frame']] = np.linspace(0,max_angle,logfile['end_frame']-logfile['start_frame'])%360
    logfile['detector']['angles'] = angles.tolist()

    return logfile


def normalise_rotation(fg_data, fg_logfile, bg_data, bg_logfile, verbose=False):
    """
    Normalise an image by a reference rotation with an empty container.

    Args:
        fg_data (ND array): The image to normalise.
        fg_logfile: The logfile for the source data.
        bg_data (ND array): The background image to remove from `data`.
        bg_logfile: The logfile for the background image set.

    Returns:
        normalised_data (ND array): The same data as the original, but with the background removed.
    """
    _,nx,ny = fg_data.shape
    nt = fg_logfile['end_frame'] - fg_logfile['start_frame']
    normalised_data = np.zeros([nt,nx,ny])

    num_fg_frames = len(fg_logfile['detector']['frames'])
    num_bg_frames = len(bg_logfile['detector']['frames'])
    bg_angles = np.array(bg_logfile['detector']['angles'])

    with np.errstate(divide='ignore',invalid='ignore'):
        for i,frame in progressbar(enumerate(range(fg_logfile['start_frame'],
                                                   fg_logfile['end_frame']))):
            fg_angle = fg_logfile['detector']['angles'][frame]%360

            j = np.nanargmin(np.abs(fg_angle - bg_angles))

            normalised_data[i] = np.nan_to_num(fg_data[frame]/bg_data[j])

            if verbose:
                plt.subplot(221)
                plt.imshow(fg_data[frame])

                plt.subplot(222)
                plt.imshow(bg_data[j])

                plt.subplot(223)
                plt.imshow(normalised_data[i])

                plt.show()

    return normalised_data
