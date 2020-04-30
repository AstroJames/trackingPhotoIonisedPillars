#!/usr/bin/env python

"""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          this script extracts the shock front and saves it in 3D
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  28 / April / 2020

"""


from header import *
from DataHandling import *
from skimage.transform import hough_line, hough_line_peaks
from skimage.morphology import square, opening

# Functions
############################################################################################################################################

def houghTransform(s, sThreshold):
    """
    DESCRIPTION:
    This function computes the hough transform, for determining the shock front.
    see: https://scikit-image.org/docs/dev/auto_examples/edges/plot_line_hough_transform.html

    INPUT:
    s           - the 2D density field ln(\rho/\rho_0)
    sThreshold  - the threshold for the field. Densities above this threshold
                will be used to identify the shock front

    OUTPUT:
    x       - x coordinate of the Hough transform
    y       - y coordinate of the Hough transform
    xMin    - x coordinate of the window around the shock front
    yMin    - y coordinate of the window around the shock front
    sMask   - the mask, s.max() * sThreshold (for checking)

    """

    print("Calculating the Hough Transform.")
    # Create a mask on s for detecting the ionisation front
    sMask = s > s.max()*sThreshold

    # run a hough transform on the masked density data
    dtheta = 0.05 # the d\theta around a vertical line approximation for the ion. front.
    # create a sample of test angles close to a vertical line
    tested_angles = np.linspace(np.pi + dtheta, np.pi - dtheta, 90)

    # run the Hough transform
    origin          = np.array((0, sMask.shape[1]))            # define the origin coordinate
    h, theta, d     = hough_line(sMask, theta=tested_angles)   # define the H. transform

    # pick the most dominant line from the H. transform
    hspace, angle, dist  = hough_line_peaks(h, theta, d)             # extract param. values
    y0, y1          = (dist[0] - origin * np.cos(angle[0])) / np.sin(angle[0])

    # pick out the cooridinates in the density field from the line and create a window
    x , y       = sampleHough(origin,[y0, y1])
    xMin        = x - windowSize
    xMax        = x + windowSize

    # Just taking a single value at the for the window
    xMin = xMin[0]
    xMax = xMax[0]

    # Make sure the window size isn't larger than the actual data domain
    if xMin < 0:
        xMin = 0
    elif xMin > sMask.shape[0]-1:
        xMin = sMask.shape[0]-1
    # and for xMax too
    if xMax < 0:
        xMax = 0
    elif xMax > sMask.shape[0]-1:
        xMax = sMask.shape[0] - 1

    return x, y, xMin, xMax, sMask


def sampleHough(xValues, yValues):
    """
    DESCRIPTION:
    This function takes the most dominant line from the Hough transform data.
    see: https://scikit-image.org/docs/dev/auto_examples/edges/plot_line_hough_transform.html

    INPUT:
    xValues - the x values from the houghTransform function
    yValues - the y values from the houghTransform function

    OUTPUT:
    x       - the x values for the most dominant line
    y       - the y values for the most dominant line

    """

    print("Extracting the most dominant line from the Hough transform.")
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
    # include coordinates in the density field domain

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

# Main script
############################################################################################################################################

if __name__ == "__main__":
    times               = np.arange(0,100)      # the times to run the code on
    sThreshold          = 0.5                   # density threshold
    sliceIndexes        = np.arange(200)        # the range of slice indexes for the 200^3 dataset
    sliceIndex          = 0                     # the index of which z plane to slice the fields through
    cubeDim             = 200                   # the dimension of the cube
    readDir             = "./testData/"         # read directory
    writeDir            = "./Data/MaskCubes/"   # write directory
    openFactor          = 2                     # the openning factor of the pixels
    windowSize          = 50                    # half the size of the window

    # Initialise a 3D cube for a mask
    mask3D = np.zeros([cubeDim,cubeDim,cubeDim])

    # for each slice through the 3D field
    for time in times:

        print("Starting iteration on file number {}".format(time))
        print("#####################################")

        for sliceIndex in sliceIndexes:

            print("Slice: {}".format(sliceIndex))
            print("#####################################")

            # read in the Data and extract into a np.array
            dens        = loadObj(readDir + "rho_{}".format(time))

            # take a slice through (x,y,z=0)
            dens    = dens[:,sliceIndex,:] # g / cm^3
            s       = np.log(dens / dens.mean())

            x, y, xMin, xMax, sMask = houghTransform(s,sThreshold)

            # Store the 2D mask
            mask3D[:,sliceIndex,:] = sMask

        # Save
        np.save(writeDir + "mask3D_{}".format(time),mask3D)
