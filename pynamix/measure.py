import numpy as np
from scipy.linalg import eig # NOT SURE WHY NUMPY ONE SEG FAULTS
from pynamix.exposure import *

def main_direction(tensor):
    """Calculate the principal orientation and orientation magnitude of a nematic order tensor.

    Args:
        tensor: 2 by 2 array represnting the nematic order tensor.

    Returns:
        Two values, one for the principal orientation in radians and one for the magnitude of the orientation on a scale of zero to one.
    """
    v,_ = eig(tensor)
    v = np.real(v) # NOTE: NOT SURE WHY THIS IS NECESSARY
    angle=np.arctan2(v[1],v[0]) # HACK: CHECK THIS!!!!!
    if (angle < 0): angle=angle + np.pi
    elif (angle > np.pi): angle=angle - np.pi
    dzeta=np.sqrt(np.sum(tensor ** 2))
    return angle, dzeta

def hanning_window(patchw=32):
    """Compute a radial hanning window.

    Args:
        patchw (int): The half width of the patch.

    Returns:
        The radial hanning window.
    """
    w = np.zeros([patchw*2,patchw*2])
    for i in range(1,patchw * 2):
        for j in range(1,patchw * 2):
            dst=np.sqrt((i - 0.5 - patchw) ** 2 + (j - 0.5 - patchw) ** 2)
            w[i,j]=0.5 * (np.cos(2 * np.pi * dst / (patchw * 2))) + 0.5
            if (dst > patchw):
                w[i,j]=0
    return w

def angular_binning(patchw=32,N=10000):
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
    K = np.zeros([N,2])
    n_nbmaskQ = np.zeros([patchw * 2,patchw * 2])
    n_maskQ   = np.zeros([patchw * 2,patchw * 2,2,2])
    for j in range(1,1000): # Number of iteration in Monte-Carlo simulation
        r = (np.random.rand(N,2) - 0.5) * 2 * patchw
        K[:,0] = r[:,0] / (np.sqrt(r[:,0] ** 2 + r[:,1] ** 2))
        K[:,1] = r[:,1] / (np.sqrt(r[:,0] ** 2 + r[:,1] ** 2))
        x = np.floor(r + patchw).astype(np.int)
        for i in range(1,N):
            n_nbmaskQ[x[i,1],x[i,0]] += 1
            n_maskQ[x[i,1],x[i,0],0,0] += K[i,0] * K[i,0]
            n_maskQ[x[i,1],x[i,0],0,1] += K[i,0] * K[i,1]
            n_maskQ[x[i,1],x[i,0],1,0] += K[i,1] * K[i,0]
            n_maskQ[x[i,1],x[i,0],1,1] += K[i,1] * K[i,1]

    # Final scaling for the n_maskQ coefficients.
    n_maskQ[:,:,0,0] /= n_nbmaskQ
    n_maskQ[:,:,0,1] /= n_nbmaskQ
    n_maskQ[:,:,1,0] /= n_nbmaskQ
    n_maskQ[:,:,1,1] /= n_nbmaskQ

    return n_maskQ

def radial_grid(rnb=200,patchw=32,N=10000): # validated against MATLAB code
    r_grid = np.linspace(0,patchw*1.5,rnb)
    nr_px = np.zeros([patchw*2,patchw*2])
    nr_pxr = np.zeros([patchw*2,patchw*2,rnb])

    Niter = patchw*2*patchw*2*200//N
    for j in range(Niter):
        r = (np.random.rand(N,2)-0.5)*2*patchw
        dst = np.sqrt(r[:,0]**2+r[:,1]**2)
        x = np.floor(r+patchw).astype(int)
        for i in range(N):
            nr_px[ x[i,1],x[i,0]] += 1
            arg = np.argmin(np.abs(r_grid - dst[i]))
            nr_pxr[x[i,1],x[i,0],arg] += 1
    for i in range(rnb):
        nr_pxr[:,:,i] /= nr_px
    r_grid += np.mean(np.diff(r_grid))*0.5
    return r_grid, nr_pxr

def orientation_map(data,tmin=0,tmax=None,tstep=1,xstep=32,ystep=32,patchw=32,normalisation=mean_std):
    """
    Calculate the principal orientation and orientation magnitude at a set of patches in images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        patchw (int): The half width of the patch.

    Returns:
        Two arrays of shape which describe the principal orientation and orientation magnitude for each patch
    """
    w = hanning_window(patchw)
    n_maskQ = angular_binning(patchw,N=100) # NOTE: LOW N FOR TESTING - REMOVE THIS LATER

    nt,nx,ny = data.shape

    gridx = range(patchw,nx-patchw,xstep) # locations of centres of patches in x direction
    gridy = range(patchw,ny-patchw,ystep) # locations of centres of patches in y direction
    if tmax is None: tmax = nt # optionally set end time

    # Prepare thre result matrices (3D), first 2 indices are the grid, the last index is time
    orient = np.nan * np.zeros([tmax-tmin,len(gridx),len(gridy)])
    dzeta  = np.nan * np.zeros([tmax-tmin,len(gridx),len(gridy)])
    Q = np.zeros_like(n_maskQ)
    Q2 = np.zeros([2,2,nx,ny,nt])

    for t,ti in enumerate(range(tmin,tmax,tstep)): # Loop on the movie frames
        for i,xi in enumerate(gridx): # Loop over the grid
            for j,yj in enumerate(gridy):
                patch = data[ti,xi-patchw:xi+patchw,yj-patchw:yj+patchw]

                patch = normalisation(patch)

                S = np.fft.fftshift(np.abs(np.fft.fft2(patch*w) ** 2))

                if np.sum(S) == 0:
                    Q[:,:,0,0] = Q[:,:,0,1] = Q[:,:,1,0] = Q[:,:,1,1] = S
                else:
                    Q[:,:,0,0] = n_maskQ[:,:,0,0].dot(S) / np.sum(S)
                    Q[:,:,0,1] = n_maskQ[:,:,0,1].dot(S) / np.sum(S)
                    Q[:,:,1,0] = n_maskQ[:,:,1,0].dot(S) / np.sum(S)
                    Q[:,:,1,1] = n_maskQ[:,:,1,1].dot(S) / np.sum(S)
                Q2[:,:,i,j,t] = np.sum(np.sum(Q,0),0)
                Q2[0,0,i,j,t] -= 0.5
                Q2[1,1,i,j,t] -= 0.5
                Q2[:,:,i,j,t] *= np.sqrt(2)
                orient[t,i,j],dzeta[t,i,j] = main_direction(Q2[:,:,i,j,t])
    return orient, dzeta

def size_map(data,logfile,tmin=0,tmax=None,tstep=1,xstep=32,ystep=32,patchw=32,normalisation=mean_std):
    """
    Calculate the radial average of the 2D FFT at a set of patches in images in a series.

    Args:
        data: The source data. Should be in the shape [nt,nx,ny]
        tmin (int): First frame to analyse in the series
        tmax (int): Last frame to analyse in the series
        tstep (int): Spacing between frames to analyse
        xstep (int): Spacing between patches in the x direction
        ystep (int): Spacing between patches in the y direction
        patchw (int): The half width of the patch.

    Returns:
        Two arrays of shape which describe the principal orientation and orientation magnitude for each patch
    """
    w = hanning_window(patchw)
    n_maskQ = angular_binning(patchw,N=100) # NOTE: LOW N FOR TESTING - REMOVE THIS LATER

    nt,nx,ny = data.shape

    gridx = range(patchw,nx-patchw,xstep) # locations of centres of patches in x direction
    gridy = range(patchw,ny-patchw,ystep) # locations of centres of patches in y direction
    if tmax is None: tmax = nt # optionally set end time

    frequencyconversion=logfile['resolution']/(patchw*2) # #(do 1/(frequencyconversion*peakfreq) to get the spatial caracteristic wavelength)
    radialspec=zeros([(tmax-tmin)//tstep,length(gridx), length(gridy), rnb]) #

    r_grid, nr_pxr = radial_grid(rnb=rnb,patchw=patchw)

    for t,ti in enumerate(range(tmin,tmax,tstep)): # Loop on the movie frames
        for i,xi in enumerate(gridx): # Loop over the grid
            for j,yj in enumerate(gridy):
                patch = data[ti,xi-patchw:xi+patchw,yj-patchw:yj+patchw]
                patch = normalisation(patch)
                S = np.fft.fftshift(np.abs(np.fft.fft2(patch*w) ** 2))

                for k in range(rnb):
                    radialspec[t,i,j,k]=np.sum(S*nr_pxr[:,:,k])  # radially SUMMED power spectral density
    return size

# Testing area
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # w = hanning_window()
    # plt.imshow(w)
    # plt.colorbar()
    # plt.show()

    # n_maskQ = angular_binning(N=1000)
    # plt.subplot(221)
    # plt.imshow(n_maskQ[:,:,0,0])
    # plt.colorbar()
    # plt.subplot(222)
    # plt.imshow(n_maskQ[:,:,0,1])
    # plt.colorbar()
    # plt.subplot(223)
    # plt.imshow(n_maskQ[:,:,1,0])
    # plt.colorbar()
    # plt.subplot(224)
    # plt.imshow(n_maskQ[:,:,1,1])
    # plt.colorbar()
    # plt.show()

    # m = np.array([[1,0],[1,0]])/np.sqrt(2)
    # angle,dzeta = main_direction(m)
    # print(np.degrees(angle),dzeta)
    # m = np.array([[1,0],[0,1]])/np.sqrt(2)
    # angle,dzeta = main_direction(m)
    # print(np.degrees(angle),dzeta)
    # m = np.array([[0,1],[1,0]])/np.sqrt(2)
    # angle,dzeta = main_direction(m)
    # print(np.degrees(angle),dzeta)

    #
    # r_grid = radial_grid()
    # print(r_grid)

    # from pynamix.data import pendulum
    # data,logfile = pendulum()
    from pynamix.io import load_seq
    from pynamix import color
    virino = color.virino()
    data,logfile = load_seq('/Volumes/LTS/DynamiX/FRC/Marta/MartaTest1-D0.log')
    orient, dzeta = orientation_map(data,tmin=1000,tmax=1002)
    # size = size_map(data,tmin=1000,tmax=1002)
    #
    plt.subplot(131)
    plt.imshow(data[1000])
    plt.subplot(132)
    plt.pcolormesh(orient[0],cmap=virino)
    plt.colorbar()
    plt.subplot(133)
    plt.pcolormesh(size[:,:,0])
    plt.colorbar()
    plt.show()
