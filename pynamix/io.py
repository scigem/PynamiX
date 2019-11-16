import os, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import inferno
from imageio import imwrite

def load_seq(filename,varian=False):
    '''Load an SEQ file and related logfile, reshaping the file into a 3D array if a logfile is present.

    :param filename: location of SEQ/log file to load. Can be the full path to the SEQ or log file including or excluding the file format.
    :type filename: string
    :param varian: Was this recorded using the Varian software?
    :type varian: bool
    :returns: A 3D array with dimensions [nt,nx,ny]
    '''
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
    '''Convert an appropriately shaped SEQ file into TIFF files. Optionally takes an angle to rotate the images by.

    :param foldername: Location to save tiffs into. If this folder does not exist, it will be created for you. including or excluding the file format.
    :type foldername: string
    :param data: The data loaded as a 3D numpy array with dimensions [nt,nx,ny]
    :type data: array
    :param vmin: Minimum value to show in tiff
    :param vmax: Maximum value to show in tiff
    :param angle: Angle in degrees to rotate the images
    :param colour: Flag to save the image in colour instead of grayscale
    :type colour: bool
    :param logscale: Flag to save the image using a logscale for the colours
    :type logscale: bool
    :param tmin: First frame to save
    :param tmax: Last frame to save
    '''
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

# Testing area
if __name__ == '__main__':
    from pynamix.data import radiographs
    data,logfile = radiographs()
    save_as_tiffs('tt',data,tmin=1000,tmax=1050,tstep=10)
