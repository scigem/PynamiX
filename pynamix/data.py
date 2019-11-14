from pynamix.io import load_seq

def radiographs():
    # load the default dataset of a pendulum swinging, recorded on a DEXELA detector at 4x4 binnging (see the logfile for more information)
    return load_seq('test_data')
