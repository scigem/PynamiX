import numpy as np

def mean_std(im):
    """
    Normalise an image such that it has zero mean and standard deviation of one.

    Args:
        im (2D array): The image to normalise.

    Returns:
        out (2D array): The normalised image.
    """
    if np.std(im) != 0: out=(im - np.mean(im)) / np.std(im)
    else:               out=(im - np.mean(im))
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

def clamp(data,vmin,vmax):
    """
    Clamp an image between two values, masking it outside of those values.

    Args:
        data (ND array): The image to clamp.
        vmin (float): The lowest value to clamp below.
        vmax (float): The highest value to clamp above.

    Returns:
        masked_data (masked ND array): The same data as the original, but with masked values outside of the defined range.
    """
    masked_data = np.ma.masked_outside(data,vmin,vmax,copy=True)
    return masked_data

def apply_ROI(data,logfile,top=0,left=0,right=None,bottom=None):
    """
    Apply an ROI to an image or a series of images.

    Args:
        data (2D or 3D array): The source data.
        top (int): the index of the top edge of the ROI.
        left (int): the index of the left edge of the ROI.
        right (int): the index of the right edge of the ROI.
        bottom (int): the index of the bottom edge of the ROI.
    """
    N = len(data.shape) # number of dimensions
    if N == 2:
        nx,ny = data.shape
        if right == None: right = nx
        if bottom == None: bottom = ny
        data_ROI = data[left:right,top:bottom]
    elif N == 3:
        _,nx,ny = data.shape
        if right == None: right = nx
        if bottom == None: bottom = ny
        data_ROI = data[:,left:right,top:bottom]
    else:
        raise Exception('ROI only defined for 2D and 3D arrays')

    logfile['detector']['ROI_software'] = {}
    logfile['detector']['ROI_software']["top"] = top
    logfile['detector']['ROI_software']["left"] = left
    logfile['detector']['ROI_software']["right"] = right
    logfile['detector']['ROI_software']["bottom"] = bottom

    return data_ROI, logfile
