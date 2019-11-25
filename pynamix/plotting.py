import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import interactive, fixed
import ipywidgets as widgets
from IPython.display import display

import matplotlib as mpl
mpl.rc('image', cmap='inferno')

def hist(data,frame,vmin=0,vmax=65535):
    """
    Show a histogram of the grey levels in an image.

    Args:
        data (ND array): Your data.
    """
    im = data[frame]
    im_masked = np.ma.masked_outside(im,vmin,vmax)
    plt.figure(99,figsize=[9,3])
    plt.subplot(131)
    plt.imshow(im)
    plt.subplot(132)
    plt.imshow(im_masked)
    plt.subplot(133)
    n,bins,patches=plt.hist(im.flatten(),bins=256)
    plt.plot([vmin,vmin],[0,np.amax(n)],'k--',lw=2)
    plt.plot([vmax,vmax],[0,np.amax(n)],'k--',lw=2)
    plt.xlabel('Pixel intensity')
    plt.ylabel('Number of pixels')
    plt.show()

def hist_GUI(data,vmin=0,vmax=65535):
    """
    Show a histogram with editable min and max values.

    Args:
        data (ND array): The image to clamp.
        vmin (float): The lowest value to clamp below.
        vmax (float): The highest value to clamp above.

    Returns:
        masked_data (masked ND array): The same data as the original, but with masked values outside of the defined range.
    """
    # frame = 0
    if len(data.shape) == 2:
        nt = 0
    else:
        nt = data.shape[0] - 1

    w = interactive(hist,
                    data  = fixed(data),
                    frame = widgets.IntSlider(min=0, max=nt,    step=1, value=0),
                    vmin  = widgets.IntSlider(min=0, max=65535, step=1, value=vmin),
                    vmax  = widgets.IntSlider(min=0, max=65535, step=1, value=vmax),
                    # continuous_update=False,
                   )

    return w
