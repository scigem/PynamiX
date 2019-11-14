import numpy as np

def main_direction_(X):
    v,V=np.eig(X)
    __,idx=np.max(np.max(V))
    angle=np.atan2(v[2,idx],v[1,idx])
    if (angle < 0):
        angle=angle + pi
    if (angle > pi):
        angle=angle - pi
    dzeta=np.sqrt(np.sum(np.sum(X ** 2)))
    return angle, dzeta

def window_mask(patchw):
    # Window mask definition: compute a radial hanning window
    w = np.zeros([patchw*2,patchw*2])
    for i in range(1,patchw * 2):
        for j in range(1,patchw * 2):
            dst=np.sqrt((i - 0.5 - patchw) ** 2 + (j - 0.5 - patchw) ** 2)
            w[i,j]=0.5 * (np.cos(2 * pi * dst / (patchw * 2))) + 0.5
            if (dst > patchw):
                w[i,j]=0
    return w

def angular_binning(patchw):
    # Angular binning definition
    # step_angle=5 / 180 * pi
    # bin_angle=range(0 - step_angle / 2,- 180,- step_angle)
    # n_angle=floor(pi / step_angle)

    # Monte-Carlo method to compute the individual Q(k) coefficients for equation 4.
    # n_maskQ is a 4D table such as Q(kx, ky)=n_maskQ(kx,ky,:,:)
    # clear('N','n_nbmaskQ','n_maskQ','x','K','r')
    K = np.zeros([1000,2])
    r = np.zeros_like(K)
    n_nbmaskQ=np.zeros(patchw * 2)
    n_maskQ=np.zeros([patchw * 2,patchw * 2,2,2])
    for j in range(1,1000): # Number of iteration in Monte-Carlo simulation
        N=10000 # Number of particle throw per iteration
        r=(np.random.rand(N,2) - 0.5) * 2 * patchw
        K[:,0]=r[:,0] / (np.sqrt(r[:,0] ** 2 + r[:,1] ** 2))
        K[:,1]=r[:,1] / (np.sqrt(r[:,0] ** 2 + r[:,1] ** 2))
        x=floor(r + patchw) + 1
        for i in range(1,N).reshape(-1):
            n_nbmaskQ[x(i,2),x(i,1)]=n_nbmaskQ(x(i,2),x(i,1)) + 1
            n_maskQ[x(i,2),x(i,1),1,1]=n_maskQ(x(i,2),x(i,1),1,1) + K(i,1) * K(i,1)
            n_maskQ[x(i,2),x(i,1),1,2]=n_maskQ(x(i,2),x(i,1),1,2) + K(i,1) * K(i,2)
            n_maskQ[x(i,2),x(i,1),2,1]=n_maskQ(x(i,2),x(i,1),2,1) + K(i,2) * K(i,1)
            n_maskQ[x(i,2),x(i,1),2,2]=n_maskQ(x(i,2),x(i,1),2,2) + K(i,2) * K(i,2)

    # Final scaling for the n_maskQ coefficients.
    n_maskQ[:,:,1,1]=n_maskQ[:,:,1,1] / n_nbmaskQ
    n_maskQ[:,:,1,2]=n_maskQ[:,:,1,2] / n_nbmaskQ
    n_maskQ[:,:,2,1]=n_maskQ[:,:,2,1] / n_nbmaskQ
    n_maskQ[:,:,2,2]=n_maskQ[:,:,2,2] / n_nbmaskQ

    return n_maskQ

def calculate_orientation(f,start,end,xstep,ystep):
    # define a grid where the orientation matrix will be computed
    # gridx=range((ROI[exp](2,1) - ROI[exp](1,1)) / 2,ROI[exp](2,1) - ROI[exp](1,1) - patchw,step)
    # gridx=[gridx,range((ROI[exp](2,1) - ROI[exp](1,1)) / 2,patchw + 1,- step)]
    # gridx=np.unique(gridx)
    # gridy=range(patchw + 1,ROI[exp](2,2) - ROI[exp](1,2) - patchw,step)
    nt,nx,ny = f.shape()
    gridx = range(xstep//2,nx-xstep//2,xstep)
    gridy = range(ystep//2,ny-ystep//2,ystep)

    # Prepare thre result matrices (3D), first 2 indices are the grid, the last index is time
    orient=NaN * np.zeros(len(gridy),len(gridx),end-start)
    dzeta =NaN * np.zeros(len(gridy),len(gridx),end-start)

    for i in range(start,end): # Loop on the movie frames
        for j in range(len(gridx)): # Loop over the grid
            for k in range(len(gridy)):
                patch=I(range(gridx[j] - patchw,gridx[j] + patchw - 1),
                        range(gridy[k] - patchw,gridy[k] + patchw - 1)
                        )
                if std2(patch) != 0:
                    patch=(patch - np.mean(np.mean(patch))) / np.std2(patch)
                else:
                    patch=(patch - np.mean(np.mean(patch)))
                S=np.fftshift(np.abs(np.fft2(patch.dot(w)) ** 2))
                if (np.sum(np.sum(S)) == 0):
                    Q[:,:,1,1]=S
                    Q[:,:,1,2]=S
                    Q[:,:,2,1]=S
                    Q[:,:,2,2]=S
                else:
                    Q[:,:,1,1]=n_maskQ[:,:,1,1].dot(S) / np.sum(np.sum(S))
                    Q[:,:,1,2]=n_maskQ[:,:,1,2].dot(S) / np.sum(np.sum(S))
                    Q[:,:,2,1]=n_maskQ[:,:,2,1].dot(S) / np.sum(np.sum(S))
                    Q[:,:,2,2]=n_maskQ[:,:,2,2].dot(S) / np.sum(np.sum(S))
                Q2[:,:,k,j,i]=np.permute(np.sum(np.sum(Q,1),2),[3,4,1,2])
                Q2[1,1,k,j,i]=Q2[1,1,k,j,i] - 1 / 2
                Q2[2,2,k,j,i]=Q2[2,2,k,j,i] - 1 / 2
                Q2[:,:,k,j,i]=np.sqrt(2) * Q2[:,:,k,j,i]
                orient[k,j,i],dzeta[k,j,i]=main_direction(Q2[:,:,k,j,i])
    return orient, dzeta
