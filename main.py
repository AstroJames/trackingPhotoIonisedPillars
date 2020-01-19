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

def labelCreator(image,openFactor):
    # Some pre-processing of the labels
    image_open      = opening(image,square(openFactor))
    allLabels       = measure.label(image)
    return allLabels

# Working script
############################################################################################################################################
if __name__ == "__main__":
    # data test directory
    dx = dy         = 0.02  #pc
    dt              = 10    #kyr
    Amin            = 3     #dxdy
    testDataDir     = "./testData/"
    globaluniqueID  = {}
    times           = np.arange(10,100) # the times to run the code on
    centerX         = []            # the x coordiante of the centroid
    centerY         = []            # the y coordinate of the centroid
    tIter           = 0             # the time iteration value
    minDisTol       = 10            # the tolerance in pixel values for tracking centroids across time (in pixel size)
    idsPerTimeStep  = {}            # the dictionary for storing the IDs each time step
    allIDs          = []            # a list of all IDs, for all time
    openFactor      = 5             # the openning factor of the pixels
    windowSize      = 30            # half the size of the window
    statsPerTimeStep = {}

    for time in times:
        # read in the Data and extract into a np.array
        dens        = loadObj(testDataDir + "rho_{}".format(time))

        # turn the time into the proper file dump by adding 100.
        time        += 100

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
        dtheta = 0.05 # the d\theta around a vertical line approximation for the ion. front.
        # create a sample of test angles close to a vertical line
        tested_angles = np.linspace(np.pi + dtheta, np.pi - dtheta, 90)

        # run the Hough transform
        origin          = np.array((0, s_mask.shape[1]))            # define the origin coordinate
        h, theta, d     = hough_line(s_mask, theta=tested_angles)   # define the H. transform

        # pick the most dominant line from the H. transform
        hspace, angle, dist  = hough_line_peaks(h, theta, d)             # extract param. values
        y0, y1          = (dist[0] - origin * np.cos(angle[0])) / np.sin(angle[0])

        # pick out the cooridinates in the density field from the line and create a window
        x , y       = sampleHough(origin,[y0, y1])
        xMin        = x - windowSize
        xMax        = x + windowSize


        labelDy = 0.93
        # check the boundary has been detected
        f, ax   = plt.subplots(1,2,figsize=(9.2,4),dpi=200)
        plt.subplots_adjust(left=0.00, bottom=0.05, right=0.95, top=0.95, wspace=-0.05, hspace=0.05)
        ax[0].imshow(s_mask,cmap=plt.cm.plasma)
        rect = patches.Rectangle((xMin[0],0),xMax[0] - xMin[0], s.shape[0],linewidth=1,edgecolor='r',facecolor='r',alpha=0.2);
        ax[0].plot((xMin[0] + (xMax - xMin)/2.,xMin[0] + (xMax - xMin)/2.),(0,s.shape[0]), '--r',linewidth=0.5)
        ax[0].annotate(r"$\xi = s > s_{\text{max}}/2$", xy=(labelDx,labelDy), xycoords = xyCoords,color="yellow",fontsize=fs-2)
        ax[0].annotate(r"$(\xi \ominus \square ) \oplus \square = $" + " {}".format(openFactor*dx) + r"$\,\text{pc}$", xy=(labelDx,labelDy-0.06), xycoords = 'axes fraction',color="yellow",fontsize=fs-2)
        ax[0].annotate(r"Hough fit", xy=(labelDx-0.45,labelDy), xycoords = xyCoords,color="blue",fontsize=fs-2)
        ax[0].annotate(r"$\mathcal{W}_{\Delta x} = $" + " {}".format(np.round(windowSize*2*dx,1)) + r"$\,\text{pc}$", xy=(labelDx-0.45,labelDy-0.06), xycoords = xyCoords,color="red",fontsize=fs-2)
        ax[0].annotate(r"$\mathcal{A} \geq \mathcal{A}_{\text{min}} = $" + " {}".format(np.round(Amin*dx*dy,3)) + r"$\,\text{pc}^2$", xy=(labelDx,0.03), xycoords = xyCoords,color="yellow",fontsize=fs-2)
        ax[0].add_patch(rect);
        ax[0].plot(x,y, '-b')
        ax[0].set_ylim((s_mask.shape[0], 0))
        ax[0].set_axis_off()

        # Now we have isolated the mixing layer so we can pick out the region in our
        # simulation to extract the pillars
        x0  = int(xMin[0])
        x1  = int(xMax[0])
        y0  = int(y[0])
        y1  = int(y[-1])

        # filter on the densities within the small regions and then create a threshold
        sML         = s[:,x0:x1]
        sML_mask    = sML > sML.max()*0.5
        sML_mask    = np.pad(sML_mask,((0,0),(x0,s.shape[1]-x1)),mode='constant')
        allLabels   = labelCreator(sML_mask,openFactor)

        # initialise plot
        if tIter == 0:
            vMax = s.max()
            vMin = s.min()

        plot    = ax[1].imshow(s,cmap=plt.cm.plasma,vmin=vMin,vmax=vMax)
        timeDt  = time*dt
        ax[1].annotate(r"$t = ${}".format(time) + r"$dt$", xy=(labelDx + 0.275-eps,labelDy-eps), xycoords = xyCoords,color="white",fontsize=fs-2)
        ax[1].annotate(r"$t = ${}".format(time) + r"$dt$", xy=(labelDx + 0.275,labelDy), xycoords = xyCoords,color="black",fontsize=fs-2)
        ax[1].set_axis_off()
        cb = plt.colorbar(plot)
        cb.set_label(r"$s = \ln(\rho/\rho_0)$",fontsize=16)


        localuniqueID   = {}    # initialise a unique ID for each centroid, for this timestep
        statsPerID      = {}    # initialise a dictionary for the stats, for each region
        regionCounter   = 0     # initialise a region counter

        # initialise a dictionary for identfiying what IDs have been used up
        # through the region iterations
        possibleKeys = globaluniqueID.copy()

        # for each disjoint region (i.e. pillar)
        for region in measure.regionprops(allLabels):

            # skip small regions
            if region.area <= Amin:
                continue

            # draw rectangle around segmented high-density regions
            minr, minc, maxr, maxc = region.bbox
            rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,fill=False, edgecolor='red', linewidth=2)

            ########################################################
            # calculate everything from the region
            ########################################################

            # select the region in the density field
            regionDens      = s[minc:maxc,minr:maxr]
            # calculate the dispersion
            densDispersion  = np.var(regionDens)
            # calculate the total mass
            mass            = sum(sum(regionDens))*dx*dy
            # area in parsecs.
            area            = region.area*dx*dy
            # perimeter in parsecs
            per             = region.perimeter*dx
            """
            TODO:
            """
            # calculate the velocity dispersion
            # calculate the forcing parameter, b


            ax[1].add_patch(rect)

            # calculate the centroids
            centroidX = maxc - (maxc - minc)/2.
            centroidY = maxr - (maxr - minr)/2.

            centerX.append(centroidX)
            centerY.append(centroidY)
            ax[1].scatter(centerX,centerY,c='b',s=1,marker='.')

            # This is where we need to give it an ID
            # Add region to the dictionary on the first iteration


            euclideanDis    = []    # initialise an array for storing the euclidena distances between
                                    # successive time-steps
            addState = 0            # initialise an on / off state for adding new keys

            # if we are on the first iteration just store all of the regions into a dictionary
            # and give a unique ID
            if tIter == 0:
                globaluniqueID[regionCounter]  = {'x':centroidX,'y':centroidY}
                statsPerID[regionCounter]      = {'mass':mass,'svar':densDispersion,
                                                       'area':area,'per':per}
                ax[1].text(centroidX,centroidY,str(regionCounter),fontsize=16) # add a number annotation
            else:
                # if not the first iteration (time)

                # for each of the possible keys
                for key in possibleKeys.keys():

                    centroidXOld = globaluniqueID[key]['x']
                    centroidYOld = globaluniqueID[key]['y']
                    # calculate the Euclidean distance between each centroid pair
                    euclideanDis.append(np.hypot( centroidX - centroidXOld,centroidY - centroidYOld))

                keysAsInts  = map(int,possibleKeys.keys())
                # if all of the possible regions are still there then move to the next iteration
                if keysAsInts == []:
                    continue
                minDis      = min(np.array(euclideanDis))
                minKey      = keysAsInts[euclideanDis.index(min(euclideanDis))]
                if minDis < minDisTol:
                    # if it is, then store the value of that centroid in the old key
                    globaluniqueID[minKey]  = {'x':centroidX,'y':centroidY}
                    statsPerID[minKey]      = {'mass':mass,'svar':densDispersion,
                                               'area':area,'per':per}
                    ax[1].text(centroidX,centroidY,minKey,fontsize=16) # annotate plot
                    possibleKeys.pop(minKey,None)
                    localuniqueID[minKey] = None
                    #print possibleKeys.keys()
                else:
                    addState = 1

                # if you get to the end of the keys with nothing satisfying then we
                # will need to add a new key

                # a switch that makes a new ID be incorporated into a set.
                if addState == 1:

                    # create a new key that is one larger than the previous max
                    newKey                 = max(map(int,globaluniqueID.keys())) + 1
                    # add it to the global centroid dictionary
                    globaluniqueID[newKey] = {'x':centroidX,'y':centroidY}
                    statsPerID[newKey]     = {'mass':mass,'svar':densDispersion,
                                               'area':area,'per':per}
                    ax[1].text(centroidX,centroidY,newKey,fontsize=16) # annotate plot

            regionCounter +=1

        # Now remove keys that are NOT in the local
        difSet = set(map(int,possibleKeys.keys()))
        if difSet != set():
            for element in difSet:
                globaluniqueID.pop(element,None)

        # Need to copy the dictionary otherwise it updates per time-step
        updateGlobal            =  globaluniqueID.copy()
        idsPerTimeStep[time]    = updateGlobal

        # Need to copy the dictionary otherwise it updates per time-step
        updateStats             = statsPerID.copy()
        allIDs                  += map(int,globaluniqueID.keys())
        statsPerTimeStep[time]  = updateStats

        #plt.tight_layout()
        plt.savefig("Plots/rho_{}.png".format(time))
        plt.close()
        print("Iteration: {} complete".format(tIter))
        tIter += 1

idsPerTimeStep["all"] = allIDs

saveObj(idsPerTimeStep, "pillarIDs")
saveObj(statsPerTimeStep, "pillarStatistics")
