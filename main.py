""""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          main file
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  16 / Jan / 2020

"""

import header
reload(header)
from dataHandling import *

# Command Line Arguments
############################################################################################################################################
ap = argparse.ArgumentParser(description='command line inputs')
ap.add_argument('-time', '--time', help='the simulation timestamp', type=int)
args = vars(ap.parse_args())

# Working script
############################################################################################################################################
if __name__ == "__main__":

    testDataDir = "./testData/"
    dens        = load_obj("rho_10")
