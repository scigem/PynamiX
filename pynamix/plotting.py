import matplotlib.pyplot as plt
import numpy as np

def hist(data):
    """
    Show a histogram of the grey levels in an image.

    Args:
        data (ND array): Your data.
    """
    plt.hist(data.flatten())
    plt.xlabel('Pixel intensity')
    plt.ylabel('Number of pixels')
    plt.show()
