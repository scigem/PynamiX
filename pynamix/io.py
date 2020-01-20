import os, json, glob, requests, shutil, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import inferno
from imageio import imwrite
from pynamix.exposure import *

def strip_seq_log(filename):
    """
    Clean up filename by removing trailing .seq or .log

    Args:
        filename (str): filename with possible file type

    Returns:
        filename (str): filename with .seq or .log stripped
    """
    if (filename[-3:] == 'seq') or (filename[-3:] == 'log'):
        filename = filename[:-4] # strip ending if present
    return filename

def load_seq(filename,varian=False):
    """Load an SEQ file and related logfile, reshaping the file into a 3D array if a logfile is present.

    Args:
        filename (str): location of SEQ/log file to load. Can be the full path to the SEQ or log file including or excluding the file format.
        varian (bool): Was this recorded using the Varian software?

    Returns:
        A 3D array with dimensions [nt,nx,ny]
    """
    filename = strip_seq_log(filename)

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

                        nt = len(logfile['detector']['frames'])
                        nx = logfile['detector']['image_size']['height'] # in landscape mode by default
                        ny = logfile['detector']['image_size']['width']
                        try:
                            data = data.reshape(nt,nx,ny)
                        except:
                            data = data.reshape(-1,nx,ny)
                            print('WARNING: WRONG NUMBER OF FRAMES IN LOG FILE. THIS IS VERY BAD.')
                        if logfile['detector']['rotate'] == 0:
                            pass
                        if logfile['detector']['rotate'] >= 1:
                            data = np.rot90(data, axes=(2,1))
                        if logfile['detector']['rotate'] >= 2:
                            data = np.rot90(data, axes=(2,1))
                        if logfile['detector']['rotate'] == 3:
                            data = np.rot90(data, axes=(2,1))
                    except:
                        # for line in g:
                            # if line.find('frames') != -1: nt = 0
                        try:
                            data = data.reshape(-1,960,768)
                            print("""WARNING: Haven't implemented old log file checking, assumed no ROI and 2x2 binning so this may be garbage. Try using pynamix.io.upgrade_logfile() to upgrade your old log file.""")
                        except:
                            raise Exception("Could not find a log file and also could not load your file. Try using pynamix.io.upgrade_logfile() to upgrade your old log file.")
                        logfile = {}

            else:
                raise Exception('No log file found!')
    return data, logfile

def load_radio_txtfiles(foldername,tmin=0,tmax=None):
    """Load a set of radiographs produced with James's Forward Projector from DEM data.

    Args:
        foldername (str): location of data load.
        tmin (int): number of first file to load.
        tmax (int): number of last file to load.
    """
    files = glob.glob(foldername + '/*.txt')
    files = sorted(files)
    if tmax == None: tmax = len(files)

    data = []
    for f in files[tmin:tmax]:
        frame = np.loadtxt(f)
        data.append(frame)

    data = np.array(data)
    logfile = {}

    return data, logfile

def upgrade_logfile(filename):
    """Load an old logfile (non-JSON formatted), and port it to JSON formatted. Deprecates the old logfile.

    Args:
        filename (str): location of log file to load. Can be the full path to the log file including or excluding the file format.
    """
    filename = strip_seq_log(filename)

    log = {'detector': {}, 'geometry': {}, 'X-rays': {}}
    log['detector']['frames'] = []
    with open(filename + '.log') as f:
        for line in f:
            # print(line.split(' '))
            l = line.split(' ')
            if (l[0] == 'Mon') or (l[0] == 'Tue') or (l[0] == 'Wed') or (l[0] == 'Thu') or (l[0] == 'Fri') or (l[0] == 'Sat') or (l[0] == 'Sun'):
                log['timestamp'] = line[:-1] # ignore newline character
            elif l[0] == 'MODE':
                log['detector']['mode'] = int(l[1])
            # elif l[0] == '768x960\n':
            #     log['detector']['image_size'] = {}
            #     log['detector']['image_size']['width']  = 768
            #     log['detector']['image_size']['height'] = 960
            #     log['detector']['image_size']['unit'] = 'px'
            # elif l[0] == '1536x1920\n':
            #     log['detector']['image_size'] = {}
            #     log['detector']['image_size']['width']  = 1536
            #     log['detector']['image_size']['height'] = 1920
            #     log['detector']['image_size']['unit'] = 'px'
            elif re.match(r'[0-9]+x[0-9]+', l[0]) is not None:
                values = l[0].split('x')
                log['detector']['image_size'] = {}
                log['detector']['image_size']['width']  = int(values[0])
                log['detector']['image_size']['height'] = int(values[1][:-1])
                log['detector']['image_size']['unit'] = 'px'

            elif l[0] == 'ROI':
                log['detector']['ROI'] = {}
                log['detector']['ROI']['top']    = int(l[2])
                log['detector']['ROI']['left']   = int(l[3][:-1]) # ignore comma
                log['detector']['ROI']['right']  = int(l[4])
                log['detector']['ROI']['bottom'] = int(l[5])
                log['detector']['ROI']['unit'] = 'px'
            elif l[0] == 'FPS':
                log['detector']['fps'] = int(l[1])
            elif l[0] == '\n':
                pass
            else:
                log['detector']['frames'].append([len(log['detector']['frames']), int(l[0]), 0])

    log['detector']['rotate'] = 0
    log['detector']['length'] = {}
    log['detector']['length']['width'] = 195.0
    log['detector']['length']['height'] = 244.0
    log['detector']['length']['unit'] = 'mm'
    log['detector']['resolution'] = log['detector']['length']['height']/log['detector']['image_size']['height']

    os.rename(filename + '.log', filename + '.log.dep')
    with open(filename + '.log', 'w') as new:
        json.dump(log, new, indent=2)

def write_seq(filename,data):
    """
    Write an SEQ file and corresponding logfile from a 3D numpy array.

    Args:
        filename (str): The name of the SEQ file to write.
        data (array): The data to write.
        logfile: The original logfile of the source data.
    """
    filename = strip_seq_log(filename)

    with open(filename + '.seq', 'wb') as f:
        for t in data[0]:
            f.write(t.tobytes())

    with open(filename + '.log', 'w') as f:
        json.dump(logfile, f, sort_keys=False, indent=1)

def generate_seq(filename, detector, mode, nbframe=10):
    """
    Write an arbitrary SEQ file for testing purposes.

    Args:
        filename (str): The name of the SEQ file to write.
        detector (int): Which detector to simulate (0, 1 or 2).
        mode (int): Which detector mode to simulate (0|1|11|22|44).
        nbframe (int) : # of frames to write
    """
    filename = strip_seq_log(filename)

    if detector == 2 and not (mode == 11 or mode == 22 or mode == 44):
        raise Exception("Mode should be 11, 22 or 44 for detector 2")
    elif detector < 2 and not (mode == 0 or mode == 1):
        raise Exception("Mode should be 0, or 1 for detectors 1 or 2")
    elif detector < 2:
        if mode == 0:
            w=768
            h=960
        elif mode == 1:
            w=1536
            h=1920
    elif detector == 2:
        w=3072
        h=3888
        if mode == 22:
            w /= 2
            h /= 2
        elif mode == 44:
            w /= 4
            h /= 4

    pattern=np.linspace(0, 256*256-1, num=w*h, dtype='<u2') #Little endian 2 bytes unsigned

    with open(filename + '.seq', 'wb') as f:
        delta = int(w*h/nbframe)
        for i in range (0, nbframe):
            f.write(pattern.tobytes())
            pattern = np.hstack((pattern[delta:],pattern[0:delta]))


def save_as_tiffs(foldername,data,vmin=0,vmax=65535,angle=0,colour=False,normalisation=no_normalisation,tmin=0,tmax=None,tstep=1):
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

    for t in range(tmin,tmax,tstep):
        im = data[t].astype(np.float)
        im = normalisation(im)
        if angle != 0: im = np.rotate(im,angle=angle,mode='edge')

        # Set correct color mapping - FIXME: THIS SHOULD BE IN EXPOSURE!!
        # if not logscale:
        #     im = (im -vmin)/(vmax - vmin)*255.0
        #     im[im<0] = 0
        #     im[im>255] = 255
        # else:
        #     im = (np.log(im) - np.log(vmin))/(np.log(vmax) - np.log(vmin))#*255.0
        #     im[im<0] = 0
        #     # im[im>255] = 1
        #     im = np.ma.masked_greater(im,1)
        #     im = im**0.5
        #     im *= 255


        if not colour:
            imwrite(foldername + '/' + str(t).zfill(5) + '.tiff',im.astype(np.uint8))
        else:
            plt.clf()
            plt.imsave(foldername + '/' + str(t).zfill(5) + '.tiff',
                       im,
                       cmap = inferno
                      )

def load_PIVLab_txtfiles(foldername,tmin=0,tmax=None,tstep=1):
    '''Load a folder full of PIVLab txtfiles and return the data as a series of arrays

    Args:
        foldername (str): Name of the folder that contains the text files to load.
        tmin (int): First textfile to load
        tmax (int): Last textfile to load
        tstep (int): Spacing between textfiles to load

    Returns:
        4 element tuple containing

        - X (2D array): The `x` location of each patch
        - Y (2D array): The `y` location of each patch
        - u (3D array): The horizontal velocity of each patch. Has the shape [nt,nx,y].
        - v (3D array): The vertical velocity of each patch. Has the shape [nt,nx,y].
    '''

    files = glob.glob(foldername + '/*.txt')
    files = sorted(files)
    nt = len(files)
    if end == None: end = nt

    for i,f in enumerate(files[tmin:tmax:tstep]):
        data = np.loadtxt(f,skiprows=3,delimiter=',')
        if i == 0:
            npoints = data.shape[0] # number of data points in total
            x_all = data[:,0]
            y_all = data[:,1]
            u = np.zeros([nt,npoints])
            v = np.zeros([nt,npoints])
        u[i] = data[:,2]
        v[i] = data[:,3]

    x = np.unique(x_all)
    y = np.unique(y_all)
    nx = len(x)
    ny = len(y)
    u = u.reshape([nt,nx,ny])
    v = v.reshape([nt,nx,ny])
    X,Y = np.meshgrid(x,y,indexing='ij')

    return X,Y,u,v

def download_file(url,local_filename):
    """
    Download a file from a web site.

    Args:
        url (str): location to get file from.
        local_filename (str): Name of file to save
    """
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

# Testing area
if __name__ == '__main__':

    data,logfile = load_seq('~/Downloads/Test-D0')
    # upgrade_logfile('old.log')

    # from pynamix.data import pendulum
    # data,logfile = pendulum()
    # save_as_tiffs('tt',data,tmin=1000,tmax=1050,tstep=10)

    # x,y,u,v = load_PIVLab_txtfiles('/Volumes/LTS/DynamiX/PerpetualAvalanche/PerpetualAvalanche-3mm-4mm-80-20/PIV/',tmin=1000,tmax=1020,tstep=5)
    # plt.quiver(x,y,u[0,:,:],v[0,:,:])
    # plt.show()

    # generate_seq('test', 1, 1, nbframe=10)
    # pass
