from pynamix.io import load_seq

def pendulum():
    """Load the default dataset of a pendulum swinging, recorded on a DEXELA detector at 4x4 binnging (see the logfile for more information)

    .. warning:: THIS DATA DOES NOT EXIST YET. FRANCOIS THIS IS YOUR JOB.
    """
    return load_seq('../data/pendulum')

# Testing area
if __name__ == '__main__':
    data, logfile = pendulum()
