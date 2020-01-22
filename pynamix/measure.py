import os, pynamix
import numpy as np
from pynamix.exposure import *
module_loc = pynamix.__file__[:-11]

def main_direction(tensor):
    """Calculate the principal orientation and orientation magnitude of a nematic order tensor.

    Args:
        tensor: 2 by 2 array representing the nematic order tensor.

    Returns:
        Two values, one for the principal orientation in radians (from zero to pi) and one for the magnitude of the orientation on a scale of zero to one.
    """
    v,V = np.linalg.eig(tensor) # eigenvalues (v) and eigenvectors (V)
    idx = np.argmax(v) # principal eigenvalue
    angle=np.arctan2(V[1,idx],V[0,idx])
    if (angle < 0):       angle += np.pi
    elif (angle > np.pi): angle -= np.pi
    dzeta=np.sqrt(np.sum(tensor ** 2))  # ||T|| (after eq. 5)
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
    if os.path.exists(module_loc + 'defaults/n_maskQ_' + str(patchw) + '_' + str(N) + '.npy'):
        n_maskQ = np.load(module_loc + 'defaults/n_maskQ_' + str(patchw) + '_' + str(N) + '.npy')
    else:
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
        np.save(module_loc + 'defaults/n_maskQ_' + str(patchw) + '_' + str(N) + '.npy',n_maskQ)
    return n_maskQ

def radial_grid(rnb=200,patchw=32,N=10000): # validated against MATLAB code
    if os.path.exists(module_loc + 'defaults/r_grid_' + str(rnb) + '_' + str(patchw) + '_' + str(N) + '.npy'):
        r_grid = np.load(module_loc + 'defaults/r_grid_' + str(rnb) + '_' + str(patchw) + '_' + str(N) + '.npy')
        nr_pxr = np.load(module_loc + 'defaults/nr_pxr_' + str(rnb) + '_' + str(patchw) + '_' + str(N) + '.npy')
    else:
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
        np.save(module_loc + 'defaults/r_grid_' + str(rnb) + '_' + str(patchw) + '_' + str(N) + '.npy',r_grid)
        np.save(module_loc + 'defaults/nr_pxr_' + str(rnb) + '_' + str(patchw) + '_' + str(N) + '.npy',nr_pxr)
    return r_grid, nr_pxr

def orientation_map(data,logfile,tmin=0,tmax=None,tstep=1,xstep=32,ystep=32,patchw=32,normalisation=mean_std):
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
        Four 2D arrays which describe: (1) the x location of the centre of each patch, (2) the y location of the centre of each patch, (3) the principal orientation and (4) the orientation magnitude for each patch.
    """
    w = hanning_window(patchw)
    n_maskQ = angular_binning(patchw)

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
                    Q[:,:,0,0] = n_maskQ[:,:,0,0]*S / np.sum(S)
                    Q[:,:,0,1] = n_maskQ[:,:,0,1]*S / np.sum(S)
                    Q[:,:,1,0] = n_maskQ[:,:,1,0]*S / np.sum(S)
                    Q[:,:,1,1] = n_maskQ[:,:,1,1]*S / np.sum(S)
                Q2[:,:,i,j,t] = np.sum(np.sum(Q,0),0)
                Q2[0,0,i,j,t] -= 0.5
                Q2[1,1,i,j,t] -= 0.5
                Q2[:,:,i,j,t] *= np.sqrt(2)
                orient[t,i,j],dzeta[t,i,j] = main_direction(Q2[:,:,i,j,t])

    X,Y = np.meshgrid(gridx,gridy,indexing='ij')
    return X, Y, orient, dzeta

def radial_FFT(data,logfile,rnb=200,tmin=0,tmax=None,tstep=1,xstep=32,ystep=32,patchw=32,normalisation=mean_std):
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

    nt,nx,ny = data.shape

    gridx = np.arange(patchw,nx-patchw,xstep) # locations of centres of patches in x direction
    gridy = np.arange(patchw,ny-patchw,ystep) # locations of centres of patches in y direction
    if tmax is None: tmax = nt # optionally set end time

    frequencyconversion = logfile['detector']['resolution']/(patchw*2) # #(do 1/(frequencyconversion*peakfreq) to get the spatial caracteristic wavelength)
    radialspec = np.zeros([(tmax-tmin)//tstep,len(gridx),len(gridy),rnb]) #
    # radialspec = np.zeros([rnb,len(gridx)*len(gridy)*(tmax-tmin)//tstep]) # JUST DURING TESTING

    r_grid, nr_pxr = radial_grid(rnb=rnb,patchw=patchw)
    wavelength = 1./(r_grid*frequencyconversion) # wavelength in mm
    n = 0
    for t,ti in enumerate(range(tmin,tmax,tstep)): # Loop on the movie frames
        for i,xi in enumerate(gridx): # Loop over the grid
            for j,yj in enumerate(gridy):
                patch = data[ti,xi-patchw:xi+patchw,yj-patchw:yj+patchw]
                patch = normalisation(patch)
                S = np.fft.fftshift(np.abs(np.fft.fft2(patch*w) ** 2))

                for k in range(rnb):
                    radialspec[t,i,j,k]=np.sum(S*nr_pxr[:,:,k])  # ortho-radially SUMMED power spectral density
                    # radialspec[k,n]=np.sum(S*nr_pxr[:,:,k])  # ortho-radially SUMMED power spectral density - JUST FOR PLOTTING DURING TESTING
                # n += 1

    return gridx, gridy, wavelength, radialspec

def average_size_map(data,logfile,rnb=200,tmin=0,tmax=None,tstep=1,xstep=32,ystep=32,patchw=32,normalisation=mean_std,wmin=None,wmax=10,return_FFTs=False):
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
    gridx, gridy, wavelength, radialspec = radial_FFT(data,logfile,rnb=rnb,tmin=tmin,tmax=tmax,tstep=tstep,xstep=xstep,ystep=ystep,patchw=patchw,normalisation=normalisation)

    if wmin == None: wmin = 2/logfile['detector']['resolution'] # use Nyquist frequency - i.e. 2 pixels per particle
    min_val = np.argmin(np.abs(wavelength-wmax)) # this is large wavelength, wavelength is sorted large to small
    max_val = np.argmin(np.abs(wavelength-wmin)) # this is small wavelength, wavelength is sorted large to small
    # print(wavelength[min_val],wavelength[max_val])
    average_size_index = np.argmax(radialspec[:,:,:,min_val:max_val],axis=3)
    average_size_index += min_val
    nt,nx,ny,_ = radialspec.shape

    # There must be a better way to do this...
    size = np.zeros([nt,nx,ny])
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
                size[t,x,y] = wavelength[average_size_index[t,x,y]]

    X,Y = np.meshgrid(gridx,gridy,indexing='ij')

    if return_FFTs:
        return X, Y, size, wavelength, radialspec
    else:
        return X, Y, size

def bidisperse_concentration_map(data,logfile,rnb=200,tmin=0,tmax=None,tstep=1,xstep=32,ystep=32,patchw=32,normalisation=mean_std,s_a=0.5,s_b=2.0,pad=1.2,return_FFTs=False,calib_func=None):
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
    gridx, gridy, wavelength, radialspec = radial_FFT(data,logfile,rnb=rnb,tmin=tmin,tmax=tmax,tstep=tstep,xstep=xstep,ystep=ystep,patchw=patchw,normalisation=normalisation)

    min_val_a = np.argmin(np.abs(wavelength-s_a*pad)) # this is large wavelength, wavelength is sorted large to small
    max_val_a = np.argmin(np.abs(wavelength-s_a/pad)) # this is small wavelength, wavelength is sorted large to small
    min_val_b = np.argmin(np.abs(wavelength-s_b*pad)) # this is large wavelength, wavelength is sorted large to small
    max_val_b = np.argmin(np.abs(wavelength-s_b/pad)) # this is small wavelength, wavelength is sorted large to small

    peak_a_index = np.argmax(radialspec[:,:,:,min_val_a:max_val_a],axis=3)
    peak_b_index = np.argmax(radialspec[:,:,:,min_val_b:max_val_b],axis=3)
    peak_a_index += min_val_a
    peak_b_index += min_val_b

    nt,nx,ny,_ = radialspec.shape

    # There must be a better way to do this...
    peak_a = np.zeros([nt,nx,ny])
    peak_b = np.zeros([nt,nx,ny])
    for t in range(nt):
        for x in range(nx):
            for y in range(ny):
    #             print(radialspec[:,:,:,peak_a_index[t,x,y]])
                peak_a[t,x,y] = radialspec[t,x,y,peak_a_index[t,x,y]]
                peak_b[t,x,y] = radialspec[t,x,y,peak_b_index[t,x,y]]

    peak_fraction =  peak_a / (peak_a + peak_b)

    if calib_func is not None:
        concentration = calib_func(peak_fraction)
    else:
        concentration = peak_fraction

    X,Y = np.meshgrid(gridx,gridy,indexing='ij')

    if return_FFTs:
        return X, Y, concentration, wavelength, radialspec
    else:
        return X, Y, concentration

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

    # r_grid = radial_grid()
    # print(r_grid)

    from pynamix.io import load_seq
    from pynamix.exposure import clamp, apply_ROI
    data,logfile = load_seq('../data/PerpetualAvalanche-Calibration-3-8-60pc')
    data,logfile = apply_ROI(data, logfile, left=600)
    print(np.amax(data))
    data = clamp(data,10000,50000)
    logfile['length'] = {}
    logfile['length']['height'] = 240. # mm
    # print(len(logfile['frames']))
    # plt.imshow(data[10])
    # plt.show()
    # wavelength, radialspec = radial_FFT(data,logfile,tmin=10,tmax=12)
    # plt.semilogx(wavelength,radialspec)
    # plt.show()

    x,y,size = average_size_map(data,logfile,tmin=10,tmax=11)

    nt,nx,ny = data.shape
    X,Y = np.meshgrid(range(nx),range(ny),indexing='ij')
    plt.subplot(121)
    plt.pcolormesh(X,Y,data[10])
    plt.gca().invert_yaxis()
    plt.subplot(122)
    plt.pcolormesh(x,y,size[0])
    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.show()


    # from pynamix.io import load_seq
    # from pynamix import color
    # virino = color.virino()
    # data,logfile = load_seq('/Volumes/LTS/DynamiX/FRC/Marta/MartaTest1-D0.log')
    # orient, dzeta = orientation_map(data,tmin=1000,tmax=1002)
    # plt.subplot(131)
    # plt.imshow(data[1000])
    # plt.subplot(132)
    # plt.pcolormesh(orient[0],cmap=virino)
    # plt.colorbar()
    # plt.subplot(133)
    # plt.pcolormesh(size[:,:,0])
    # plt.colorbar()
    # plt.show()
