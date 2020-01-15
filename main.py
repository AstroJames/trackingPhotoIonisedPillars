""""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          main file
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  16 / Jan / 2020

"""

import py_compile
py_compile.compile("header.py")
from header import *
py_compile.compile("dataHandling.py")
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
    dens        = loadObj(testDataDir + "rho_10")
