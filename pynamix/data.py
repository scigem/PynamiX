import os
import pynamix
import numpy as np
import matplotlib.pyplot as plt
from pynamix.io import load_seq, download_file

_PYNAMIX_ROOT = os.path.abspath(pynamix.__file__)[:-19]  # location where pynamix has been installed to


def pendulum():
    """Load the default dataset of a pendulum swinging, recorded on a DEXELA detector at 4x4 binnging (see the logfile for more information)

    .. warning:: THIS DATA DOES NOT EXIST YET. FRANCOIS THIS IS YOUR JOB.
    """
    if not os.path.exists("pendulum.seq"):
        response = input("Data does not exist locally. Do you want to download it? This may take a while (y/n) ")
        if (response == "y") or (response == "Y"):
            # if not os.path.exists(_PYNAMIX_ROOT + 'data/'): os.mkdir(_PYNAMIX_ROOT + 'data/')
            # download_file('http://www.benjymarks.com/pynamix/data/Test2.log',_PYNAMIX_ROOT + 'data/pendulum.log')
            # download_file('http://www.benjymarks.com/pynamix/data/Test2.log',_PYNAMIX_ROOT + 'data/pendulum.seq')
            download_file("http://www.benjymarks.com/pynamix/data/pendulum.log", "pendulum.log")
            download_file("http://www.benjymarks.com/pynamix/data/pendulum.seq", "pendulum.seq")
            print("Data successfully downloaded")
        else:
            raise Exception("Built in data file does not exist")
    return load_seq("pendulum", varian=True)


def spiral():
    """Generate an image of a spiral to use for calibration purposes"""
    fig = plt.figure(figsize=[4, 4])
    ax = plt.subplot(111)

    N = 100001
    r = np.linspace(0, 1, N)
    theta = np.linspace(0, 2.0 * np.pi * 25, N)
    lw = 4

    x = 0.5 + r * np.cos(theta)
    y = 0.5 + r * np.sin(theta)

    plt.plot(x, y, "k-", lw=lw)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    # plt.show()
    plt.xticks([])
    plt.yticks([])
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.savefig("spiral.png", dpi=200)


def fibres(theta_mean=0.0, kappa=1.0, N=500, dpi=200, lw=4, alpha=0.2, foldername="."):
    """Generate a fake image of fibres following a known distribution.

    Args:
        N (int): Number of fibres to draw. More makes the image darker.
        theta_mean (float): An angle between 0 and 2 pi that is the average fibre orientation.
        kappa (float): A number greater than 0 that defines the alignment of the particles. kappa -> inf is perfectly aligned, kappa -> 0 is heterogeneous. This parameter is used in the vonmises circular distribution to generate random angles.
        dpi (int): Resolution of final image
        lw (int): Width of fibres in image
        alpha (float): Transparency of the fibres. Default is 0.2.
        foldername (path): Where to export the saved frames to. By default the current path

    Returns:
        Saves a file called `fibres_theta_kappa_N.png`
    """
    from scipy.stats import vonmises

    fig = plt.figure(figsize=[4, 4])
    ax = plt.subplot(111)

    for i in range(N):
        theta = np.pi + theta_mean + vonmises.rvs(kappa, size=1)
        x = np.random.rand() - 0.5
        y = np.random.rand() - 0.5
        plt.plot(
            [x - np.cos(theta), x + np.cos(theta)],
            [y - np.sin(theta), y + np.sin(theta)],
            "k-",
            alpha=alpha,
            lw=lw,
        )

    plt.xlim(-0.25, 0.25)
    plt.ylim(-0.25, 0.25)
    # plt.show()
    plt.xticks([])
    plt.yticks([])
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.savefig(
        foldername + "/fibres_" + str(theta_mean) + "_" + str(kappa) + "_" + str(N) + ".png",
        dpi=dpi,
    )

    plt.close()


# Testing area
if __name__ == "__main__":
    data, logfile = pendulum()
