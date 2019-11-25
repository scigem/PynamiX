import os
import pynamix
from pynamix.io import load_seq, download_file

_PYNAMIX_ROOT = os.path.abspath(pynamix.__file__)[:-19] # location where pynamix has been installed to

def pendulum():
    """Load the default dataset of a pendulum swinging, recorded on a DEXELA detector at 4x4 binnging (see the logfile for more information)

    .. warning:: THIS DATA DOES NOT EXIST YET. FRANCOIS THIS IS YOUR JOB.
    """
    if not os.path.exists('pendulum.seq'):
        response = input('Data does not exist locally. Do you want to download it? This may take a while (y/n) ')
        if (response == 'y') or (response == 'Y'):
            # if not os.path.exists(_PYNAMIX_ROOT + 'data/'): os.mkdir(_PYNAMIX_ROOT + 'data/')
            # download_file('http://www.benjymarks.com/pynamix/data/Test2.log',_PYNAMIX_ROOT + 'data/pendulum.log')
            # download_file('http://www.benjymarks.com/pynamix/data/Test2.log',_PYNAMIX_ROOT + 'data/pendulum.seq')
            download_file('http://www.benjymarks.com/pynamix/data/pendulum.log','pendulum.log')
            download_file('http://www.benjymarks.com/pynamix/data/pendulum.seq','pendulum.seq')
            print('Data successfully downloaded')
        else:
            raise Exception('Built in data file does not exist')
    return load_seq('pendulum',varian=True)


# Testing area
if __name__ == '__main__':
    data, logfile = pendulum()
