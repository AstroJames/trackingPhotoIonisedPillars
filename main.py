#!/usr/bin/env python2

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
ap.add_argument('-time', '--time',default=10,help='the simulation timestamp',type=int)
ap.add_argument('-viz', '--viz',default=None,help='visualisation setting', type=str)
args = vars(ap.parse_args())

# Working script
############################################################################################################################################
if __name__ == "__main__":
    # data test directory
    testDataDir = "./testData/"

    # read in the Data and extract into a np.array
    dens        = loadObj(testDataDir + "rho_{}".format(args['time']))

    # take a slice through (x,y=0,z)
    dens    = dens[:,0,:]
    s       = np.log(dens / dens.mean())

    # just a quick check of the field
    if args['viz'] == "field":
        plt.figure(dpi=200)
        plt.imshow(s,cmap=plt.cm.plasma)
        plt.show()

    s_mask = s > s.max()*0.5

    if args['viz'] == "mask":
        plt.figure(dpi=200)
        plt.imshow(s_mask)
        plt.show()
