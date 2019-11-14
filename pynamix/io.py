import json
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from imageio import imwrite as imsave

def load_seq(filename,varian=False):
    '''
    Load an SEQ file and related logfile
    '''
    if filename[-3:] == 'seq': filename = filename[:-4] # strip seq file ending if present
    f = open(filename + '.seq','rb')

    if varian:
        a = memmap(f,dtype="H",mode='r',offset=2048)
        logfile = {}
    else:
        a = memmap(f,dtype="H",mode='r')
        if os.path.exists(filename + '.logfile'):
            g = open(filename + '.logfile')
            logfile = json.load(g)
            nt = logfile['']
            ny = logfile['']
            nx = logfile['']
            a.reshape(nt,ny,nx)
        else:
            logfile = {}
    return a, logfile

def save_as_tiffs(pathname,b,vmin=0,vmax=65535,angle=0,gray=True,logscale=False):
    '''
    Convert an appropriately shaped SEQ file into TIFF files. Optionally takes an angle to rotate the images by.
    '''
    nt = b.shape[0]
    if logscale:
        norm = LogNorm()
    else:
        norm = Normalize()
    for i in range(nt):
        im = b[i].astype(float)

        # Optionally rotate image
        im = rotate(im,angle=angle,mode='edge')

        # Set correct color mapping
        im = (im -vmin)/(vmax - vmin)*255.0
        im[im<0] = 0
        im[im>255] = 255

        if gray:
            imsave(pathname + '/' + str(i).zfill(5) + '.tiff',im.astype(uint8))
        else:
            plt.clf()
            plt.imshow(im,
    #                    norm = LogNorm(),
                       cmap = inferno,
                      )
            plt.subplots_adjust(left=0,right=1,bottom=0,top=1)
            plt.axis('off')
            plt.savefig(pathname + '/' + str(i).zfill(5) + '.tiff',dpi=72)
