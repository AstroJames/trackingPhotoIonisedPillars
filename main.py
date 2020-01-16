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
from skimage.morphology import square, opening
from skimage import filters, measure
import matplotlib.patches as patches
import matplotlib.patches as mpatches
from skimage.color import label2rgb

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

def labelCreator(image):
    # Some pre-processing of the labels
    openfactor      = 10;
    image_open      = opening(image,square(openfactor))
    allLabels       = measure.label(image)
    return allLabels

# Working script
############################################################################################################################################
if __name__ == "__main__":
    # data test directory
    testDataDir     = "./testData/"
    globaluniqueID  = {}
    times           = np.arange(10,100)
    centerX         = []
    centerY         = []
    tIter           = 0             # the time iteration value
    minDisTol       = 5            # the tolerance in pixel values for tracking centroids across time

    for time in times:
        # read in the Data and extract into a np.array
        dens        = loadObj(testDataDir + "rho_{}".format(time))

        # take a slice through (x,y,z=0)
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
            rect = patches.Rectangle((xMin[0],0),xMax[0] - xMin[0], s.shape[0],linewidth=1,edgecolor='r',facecolor='b',alpha=0.5);
            ax.add_patch(rect);
            ax.plot(x,y, '-b')
            ax.set_ylim((s_mask.shape[0], 0))
            ax.set_axis_off()
            plt.show()


        # Now we have isolated the mixing layer so we can pick out the region in our
        # simulation to extract the pillars
        x0  = int(xMin[0])
        x1  = int(xMax[0])
        y0  = int(y[0])
        y1  = int(y[-1])

        # filter on the densities within the small regions and then create a threshold
        sML      = s[:,x0:x1]
        sML_mask = sML > sML.max()*0.5
        sML_mask = np.pad(sML_mask,((0,0),(x0,s.shape[1]-x1)),mode='constant')
        allLabels = labelCreator(sML_mask)

        # initialise measurement arrays
        regionPerimeter    = [];
        regionArea         = [];

        # initialise plot
        if tIter == 0:
            vMax = s.max()
            vMin = s.min()

        f,ax = plt.subplots(dpi=200)
        plot = ax.imshow(s,cmap=plt.cm.plasma,vmin=vMin,vmax=vMax)
        ax.set_axis_off()
        cb = plt.colorbar(plot)
        cb.set_label(r"$s = \ln(\rho/\rho_0)$",fontsize=16)


        localuniqueID = {}  # initialise a unique ID for each centroid, for this timestep
        regionCounter = 0   # initialise a region counter

        # for each disjoint region (i.e. pillar)
        possibleKeys = globaluniqueID.copy()
        for region in measure.regionprops(allLabels):

            # skip small regions
            if region.area <= 3:
                continue

            # Add to the area vector.
            regionArea.append(region.area)

            # Add to the perimeter vector.
            regionPerimeter.append(region.perimeter)

            # draw rectangle around segmented high-density regions
            minr, minc, maxr, maxc = region.bbox
            rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,fill=False, edgecolor='red', linewidth=2)
            ax.add_patch(rect)

            # calculate the centroids
            centroidX = maxc - (maxc - minc)/2.
            centroidY = maxr - (maxr - minr)/2.

            centerX.append(centroidX)
            centerY.append(centroidY)
            ax.scatter(centerX,centerY,c='b',s=1,marker='.')

            # This is where we need to give it an ID
            # Add region to the dictionary on the first iteration

            euclideanDis = []
            addState = 0    # initialise an on / off state for adding new keys
            if tIter == 0:
                globaluniqueID[str(regionCounter)] = (centroidX,centroidY)
                ax.text(centroidX,centroidY,str(regionCounter),fontsize=16) # add a number annotation
            else:
                for key in possibleKeys.keys():

                    centroidXOld, centroidYOld = globaluniqueID[key]
                    # calculate the Euclidean distance between each centroid pair
                    euclideanDis.append(np.hypot( centroidX - centroidXOld,centroidY - centroidYOld))

                keysAsInts  = map(int,possibleKeys.keys())
                minDis      = min(np.array(euclideanDis))
                minKey      = str(keysAsInts[euclideanDis.index(min(euclideanDis))])

                if minDis < minDisTol:
                    # if it is, then store the value of that centroid in the old key
                    globaluniqueID[minKey] = (centroidX,centroidY)
                    ax.text(centroidX,centroidY,minKey,fontsize=16)
                    possibleKeys.pop(minKey,None)
                    localuniqueID[minKey] = None
                    print possibleKeys.keys()
                else:
                    addState = 1

                # if you get to the end of the keys with nothing satisfying then we
                # will need to add a new key

                # a switch that makes a new ID be incorporated into a set.
                if addState == 1:

                    # create a new key that is one larger than the previous max
                    newKey                      = max(map(int,globaluniqueID.keys())) + 1
                    globaluniqueID[str(newKey)] = (centroidX,centroidY)
                    ax.text(centroidX,centroidY,newKey,fontsize=16)

            regionCounter +=1

        # Now remove keys that are NOT in the local
        difSet = set(map(int,possibleKeys.keys()))
        if difSet != set():
            for element in difSet:
                globaluniqueID.pop(str(element),None)


        plt.tight_layout()
        plt.savefig("Plots/rho_{}.png".format(time))
        plt.close()
        print("Iteration: {} complete".format(tIter))
        tIter += 1
