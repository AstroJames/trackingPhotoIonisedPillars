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
from skimage.transform import hough_line, hough_line_peaks
import matplotlib.patches as patches

# Command Line Arguments
############################################################################################################################################
ap = argparse.ArgumentParser(description='command line inputs')
ap.add_argument('-time', '--time',default=10,help='the simulation timestamp',type=int)
ap.add_argument('-viz', '--viz',default=None,help='visualisation setting', type=str)
args = vars(ap.parse_args())

# Functions
############################################################################################################################################

def sampleHough(xValues,yValues):
    """
    Extract density field coordinates from the most
    dominant line in the Hough transform

    """

    # x values
    x0 = xValues[0]
    x1 = xValues[1]

    # y values
    y0 = yValues[0]
    y1 = yValues[1]

    # parameters of the line
    m =  ( y1 - y0 ) / ( x1 -  x0 )
    c = y0 - m*x0

    # sample the space between 0 and the dimension of the
    # field
    x = np.linspace(x0,x1,10000)

    # create the coordinates
    y = m*x + c

    # Now downsample the coordinates to only
    # include cooridnates in the density field domain

    # Create some empty arrays
    xIn = []
    yIn = []
    iterCount = 0
    for xCoord, yCoord in zip(x,y):
        # only take those coordinates that are within the
        # x0 and x1 space.
        if yCoord >= x0 and yCoord <= x1:
            xIn.append(xCoord)
            yIn.append(yCoord)

    # make some new arrays and return
    x = np.array(xIn)
    y = np.array(yIn)

    return x, y

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
        f, ax = plt.subplots(dpi=200)
        ax.imshow(s,cmap=plt.cm.plasma)
        ax.set_axis_off()
        plt.show()

    # Create a mask on s for detecting the ionisation front
    s_mask = s > s.max()*0.5

    # run a hough transform on the masked density data
    eps = 0.05 # the d\theta around a vertical line approximation for the ion. front.
    # create a sample of test angles close to a vertical line
    tested_angles = np.linspace(np.pi + eps, np.pi - eps, 90)

    # run the Hough transform
    origin          = np.array((0, s_mask.shape[1]))            # define the origin coordinate
    h, theta, d     = hough_line(s_mask, theta=tested_angles)   # define the H. transform

    # pick the most dominant line from the H. transform
    _, angle, dist  = hough_line_peaks(h, theta, d)             # extract param. values
    y0, y1          = (dist[0] - origin * np.cos(angle[0])) / np.sin(angle[0])

    # pick out the cooridinates in the density field from the line and create a window
    x , y       = sampleHough(origin,[y0, y1])
    windowSize  = 30
    xMin        = x - windowSize
    xMax        = x + windowSize

    # check the boundary has been detected
    if args['viz'] == "mask":
        f, ax   = plt.subplots(dpi=200)
        ax.imshow(s,cmap=plt.cm.plasma)
        rect = patches.Rectangle((xMin[0],y[0]),xMax[0] - xMin[0], y[-1]-y[0],linewidth=1,edgecolor='r',facecolor='b',alpha=0.5);
        ax.add_patch(rect);
        ax.plot(x,y, '-b')
        #ax.plot(xMin,y, '--r')
        #ax.plot(xMax,y, '--r')
        ax.set_ylim((s_mask.shape[0], 0))
        ax.set_axis_off()
        plt.show()

    # Now we have isolated the mixing layer so we can
    sMixingLayer =
