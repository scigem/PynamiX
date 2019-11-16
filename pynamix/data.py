from pynamix.io import load_seq

def radiographs():
    """Load the default dataset of a pendulum swinging, recorded on a DEXELA detector at 4x4 binnging (see the logfile for more information)

    .. warning:: THIS DATA DOES NOT EXIST YET
    """
    return load_seq('../data/test_data')

# Testing area
if __name__ == '__main__':
    data, logfile = radiographs()
