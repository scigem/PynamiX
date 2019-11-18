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
        out (2D array): The **un**normalised image.
    """
    return im
