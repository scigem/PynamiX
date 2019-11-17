import os, json, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import inferno
from imageio import imwrite

def load_seq(filename,varian=False):
    """Load an SEQ file and related logfile, reshaping the file into a 3D array if a logfile is present.

    Args:
        filename (str): location of SEQ/log file to load. Can be the full path to the SEQ or log file including or excluding the file format.
        varian (bool): Was this recorded using the Varian software?

    Returns:
        A 3D array with dimensions [nt,nx,ny]
    """
    if filename[-3:] == 'seq': filename = filename[:-4] # strip seq file ending if present
    elif filename[-3:] == 'log': filename = filename[:-4] # strip seq file ending if present

    with open(filename + '.seq','rb') as f:
        if varian:
            data = np.memmap(f,dtype="H",mode='r',offset=2048)
            try: data = data.reshape(-1,960*2,768*2)
            except:
                try: data = data.reshape(-1,960,768)
                except: pass
            logfile = {} # TODO
        else:
            data = np.memmap(f,dtype="H",mode='r')
            if os.path.exists(filename + '.log'):
                with open(filename + '.log') as g:
                    try:
                        logfile = json.load(g)
                        nt = len(logfile['frames'])
                        ny = logfile['resolution']['width']
                        nx = logfile['resolution']['height']
                        try:
                            data = data.reshape(nt,ny,nx)
                        except:
                            data = data.reshape(-1,ny,nx)
                            print('WARNING: REMOVE THIS TRY STATEMENT AFTER GETTING REAL TEST DATA')
                    except:
                        # for line in g:
                            # if line.find('frames') != -1: nt = 0
                        try: data = data.reshape(-1,960,768)
                        except: pass
                        print("WARNING: Haven't implemented old log file checking, assumed no ROI and 2x2 binning")
                        logfile = {} # TODO

            else:
                raise Exception('No log file found!')
    return data, logfile

def save_as_tiffs(foldername,data,vmin=0,vmax=65535,angle=0,colour=False,logscale=False,tmin=0,tmax=None,tstep=1):
    """Convert an appropriately shaped SEQ file into TIFF files. Optionally takes an angle to rotate the images by.

    Args:
        foldername (str): Location to save tiffs into. If this folder does not exist, it will be created for you. including or excluding the file format.
        data (array): The data loaded as a 3D numpy array with dimensions [nt,nx,ny]
        vmin (int): Minimum value to show in tiff
        vmax (int): Maximum value to show in tiff
        angle (float): Angle in degrees to rotate the images
        colour (bool): Flag to save the image in colour instead of grayscale
        logscale (bool): Flag to save the image using a logscale for the colours
        tmin (int): First frame to save
        tmax (int): Last frame to save
    """
    nt = data.shape[0]
    if tmax == None: tmax = nt
    if not os.path.exists(foldername): os.makedirs(foldername)
    if logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    for t in range(tmin,tmax,tstep):
        im = data[t].astype(np.float)

        # Optionally rotate image
        if angle != 0: im = np.rotate(im,angle=angle,mode='edge')

        # Set correct color mapping
        im = (im -vmin)/(vmax - vmin)*255.0
        im[im<0] = 0
        im[im>255] = 255

        if not colour:
            imwrite(foldername + '/' + str(t).zfill(5) + '.tiff',im.astype(np.uint8))
        else:
            plt.clf()
            plt.imshow(im,
    #                    norm = LogNorm(),
                       cmap = inferno,
                      )
            plt.subplots_adjust(left=0,right=1,bottom=0,top=1)
            plt.axis('off')
            plt.savefig(foldername + '/' + str(t).zfill(5) + '.tiff',dpi=72)

def load_PIVLab_txtfiles(foldername):
    '''Load a folder full of PIVLab txtfiles and return the data as a series of arrays

    Args:
        foldername (str): Name of the folder that contains the text files to load.

    Returns:
        4 element tuple containing

        - x (1D array): The `x` location of each patch
        - y (1D array): The `y` location of each patch
        - u (3D array): The horizontal velocity of each patch. Has the shape [nt,nx,y].
        - v (3D array): The vertical velocity of each patch. Has the shape [nt,nx,y].
    '''

    files = glob.glob(foldername + '/*.txt')
    files = sorted(files)
    nt = len(files)
    for i,f in enumerate(files):
        data = np.loadtxt(f,skiprows=3,delimiter=',')
        if i == 0:
            np = data.shape[0] # number of data points in total
            x_all = data[:,0]
            y_all = data[:,1]
            u = np.zeros([nt,np])
            v = np.zeros([nt,np])
        # elif i % 500 == 0: print('Up to file ' + str(i))
        u[i] = data[:,2]
        v[i] = data[:,3]

    x = np.unique(x_all)
    y = np.unique(y_all)
    nx = len(x)
    ny = len(y)
    u = u.reshape([nt,nx,ny])
    v = v.reshape([nt,nx,ny])

    return x,y,u,v

# Testing area
if __name__ == '__main__':
    from pynamix.data import radiographs
    data,logfile = radiographs()
    save_as_tiffs('tt',data,tmin=1000,tmax=1050,tstep=10)
