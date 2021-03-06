#!/usr/bin/env python

""""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          main file
    Author:         James Beattie & Shyam Menon
    First Created:  16 / Jan / 2020

"""

from header import *
import py_compile
py_compile.compile("header.py")
py_compile.compile("DataHandling.py")
from DataHandling import *
py_compile.compile("ReferenceUnits.py")
from ReferenceUnits import *
from PhysicalConstantsCGS import *

# Command Line Arguments
############################################################################################################################################
ap = argparse.ArgumentParser(description='command line inputs')
ap.add_argument('-tend', '--tend',default=11,help='the tend simulation timestamp',type=int)
ap.add_argument('-viz', '--viz',default=None,help='visualisation setting', type=str)
ap.add_argument('-write','--write',default=False,help='an argument for writing or not writing pickle files',type=bool)
ap.add_argument('-vel','--vel',default=False,help='an argument for loading in the velocity',type=bool)
ap.add_argument('-ionX','--ionX',default=False,help='an argument for loading in the ionised fraction field',type=bool)
ap.add_argument('-clear','--clear',default=False,help='clear all of the old plotting files',type=bool)
args = vars(ap.parse_args())


# Command Examples
############################################################################################################################################
"""

run main -tend 11 -write True -clear True -vel True -ionX True

> ends the detection code on file number 100
> writes the statistics data to some .pickle files
> clears the plotting directory
> computes the velocity statistics, as well as the density statistics
    (i.e. loads in the velocity field)
> loads the ionised fraction and extract neutral coordinates from the region

"""


# Functions
############################################################################################################################################

def labelCreator(dens, openFactor):
    """
    DESCRIPTION:
    This function creates labels for each of the distinct regions in the density field.
    see: https://scikit-image.org/docs/0.12.x/api/skimage.measure.html

    INPUT:
    dens        - the density field for labelling
    openFactor  - the opening factor for the field
                see: https://en.wikipedia.org/wiki/Opening_(morphology)
                which determines how large the square is for the hit-miss operator

    OUTPUT:
    allLabels   - the labelled regions stored in a list

    """

    # Some pre-processing of the labels
    dens_open      = opening(dens,square(openFactor))
    allLabels      = measure.label(dens)

    return allLabels


def getCmap(n, name='hsv'):
    """
    DESCRIPTION:
    Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.

    INPUT:
    n       - the number of unique colours
    name    - the type of color map

    OUTPUT:
    cm      - the n number of colors from the colour map

    """

    cm = plt.cm.get_cmap(name, n)

    return cm


def ionFractionFilter(time, regionCoords, sliceIndex):
    """
    DESCRIPTION:
    This function samples the region for just the neutrals. Ionised contaminents will influence the
    turbulent Mach number calculation.

    INPUT:
    time        - the time index for the data (0-100)
    regionCoords:
    minr        - the minimum x value of the bounding box for the region
    maxr        - the maximum x value of the bounding box for the region
    minc        - the minimum y value of the bounding box for the region
    maxc        - the maximum y value of the bounding box for the region

    OUTPUT:
    maskCoords  - the coordinates of the neutrals in the region

    """

    # read in the region coords
    minr, maxr, minc, maxc = regionCoords

    # Read in the ion fraction and extract the region
    ionX        = loadObj(readDataDir + "ionX_{}".format(time))[:,sliceIndex,:]
    ionXRegion  = ionX[minr:maxr,minc:maxc]

    # Extract ionised region and filter the denisty
    print("Filtering density for  < 1e-4 ionised cells.")
    ionXMask    = ionXRegion < 1e-4
    maskCoords  = np.nonzero( ionXMask )

    # Compute the size of the region and the amount of neutrals in it
    sizeOfRegion    = len(ionXRegion.ravel())
    sizeOfIonised   = len(maskCoords[0])

    # Compute and report the fraction
    ionisationFraction    = float(sizeOfIonised) / float(sizeOfRegion)

    # no need to keep this in memory
    del ionX

    return maskCoords, ionisationFraction


def computeVelocity(dens, time, sliceIndex):
    """
    DESCRIPTION:
    This function computes the turbulent velocity component of the
    vector field by subtracting the average bulk (centre of mass) motion.

    INPUT:
    dens - the density field for computing the centre of momentum frame (g/cm^3)
    time - the time index for the data (0-100)

    OUTPUT:
    v   - the turbulent velocity component of the velocity field (cm / s)

    """

    print("Computing the turbulent component of the velocity field.")
    # construct the magnitude of v field
    vx = loadObj(readDataDir + "vx_{}".format(time))[:,sliceIndex,:] # cm /s
    vy = loadObj(readDataDir + "vy_{}".format(time))[:,sliceIndex,:] # cm /s
    vz = loadObj(readDataDir + "vz_{}".format(time))[:,sliceIndex,:] # cm /s

    # calculate the centre of mass velocity vector
    vx_cm = sum(sum(vx * dens)) / sum(sum(dens))
    vy_cm = sum(sum(vy * dens)) / sum(sum(dens))
    vz_cm = sum(sum(vz * dens)) / sum(sum(dens))

    # construct the velocity, v_turb, which is the bulk
    # velocity minus the center of mass velocity
    vx_turb = vx - vx_cm
    vy_turb = vy - vy_cm
    vz_turb = vz - vz_cm

    # construct the magnitude of v
    v = np.sqrt(vx_turb**2 + vy_turb**2 + vz_turb**2)
    del vx, vy, vz, vx_turb, vy_turb, vz_turb

    return v


def computeVirial(densRegion, regionCoordinates, velocityDispersion, time, sliceIndex):
    """
    DESCRIPTION:


    INPUT:
    dens - the density field for computing the centre of momentum frame (g/cm^3)
    time - the time index for the data (0-100)

    OUTPUT:
    virial - the virial parameter, \alpha  = 2E_{turb} / |E_grav|

    """

    # Read in region coordinates
    minr, maxr, minc, maxc = regionCoordinates


    # read in the gravitational field
    phi     = loadObj(readDataDir + "virial_{}".format(time))[:,sliceIndex,:] # ergs


    # calculate the virial parameter (see Menon et al. 2020; Beattie et al. 2020)
    virial  = sum( sum( densRegion * velocityDispersion**2 ) ) / sum( sum( densRegion * abs(phi) ) )

    return virial


def computeFreeFallTime(densRegion):
    """
    DESCRIPTION:
    Both the mean free-fall time of the region and the multi free-fall time,
    using a spherical approximation.

    INPUT:
    densRegion  - the 2D density field of the region

    OUTPUT:
    tff         -
    tffMean     -

    """

    tff     = np.sqrt( ( 3*np.pi ) / (32 * G_GravityConstant * densRegion) )            # sec
    tffMean = np.sqrt( ( 3*np.pi ) / (32 * G_GravityConstant * np.mean(densRegion) ) )  # sec

    return tff, tffMean


def computeSFR(virial, sRegion, densRegion, regionMach, mass, tff, tffMean):
    """
    DESCRIPTION:
    Compute the star formation rate (SFR) of a pillar region by calculating the critical denisty,
    and approximating the s PDF with a Gaussian fit.

    INPUT:


    OUTPUT:
    SFR     - the predicted star formation rate of the region
    sCrit   - the critical density for stars to collapse.

    """

    phiX = 1

    sCrit   = np.log( phiX * ( np.pi/5. ) * virial * regionMach**2  )

    # compute the s PDF
    # fit a gaussian or Hopkins PDF
    # integrate the gaussian from s_crit to a large vale (infinity)


    # placeholder
    SFR = None


    return SFR, sCrit

class PillarRegion:
    """
    DESCRIPTION:
    The main pillar class where all of the calculations will go.
    """

    def __init__(self,region,densRegion,sRegion,regionCoords):
        """
        DESCRIPTION:

        INPUT:

        OUTPUT:

        """
        self.time               = time
        self.region             = region
        self.sRegion            = sRegion
        self.densRegion         = densRegion
        self.regionCoords       = regionCoords
        self.mass               = sum(sum(self.densRegion))*dx*dy*dz*(Parsec)**3 * (1. / SolarMass)
        self.area               = region.area*dx*dy
        self.perimeter          = region.perimeter*dx
        self.centroids          = region.centroid
        self.neutralCoords      = []
        self.ionisationFraction = []
        self.sigmaS             = []
        self.velRegion          = []
        self.id                 = []
        self.Mach               = []
        self.bDrive             = []
        self.tff                = []
        self.tffMean            = []
        self.virial             = []
        self.SFR                = []


    def ComputeIonFrac(self):
        """
        DESCRIPTION:

        INPUT:

        OUTPUT:

        """
        neutralCoords, ionisationFraction   = ionFractionFilter(self.time,self.regionCoords)
        self.neutralCoords                  = neutralCoords
        self.ionisationFraction             = ionisationFraction

    def ComputeDensVar(sefl):
        """
        DESCRIPTION:

        INPUT:

        OUTPUT:

        """
        if self.neutralCoords == []:
            self.sigmaS = np.var(self.sRegion)
        else:
            self.sigmaS = np.var(self.sRegion[self.neutralCoords])

    def ComputerVelRegion(self,v):
        """
        DESCRIPTION:

        INPUT:

        OUTPUT:

        """
        if self.neutralCoords == []:
            self.regionVel   = v[self.neutralCoords]
        else:
            minr, minc, maxr, maxc = self.region.bbox
            self.regionVel   = v[minr:maxr,minc:maxc]

        velVar      = self.regionVel.std()
        self.Mach   = velVar / cs

    def ComputeB(self):
        self.b    = np.sqrt( ( np.exp(densDispersion) - 1 ) / regionMach**2  )

    def ComputeFreeFallTime(self):
        tff, tffMean    = computeFreeFallTime(self.densRegion)
        self.tff        = tff
        self.tffMean    = tffMean

    def ComputeVirial(self):
        virial = computeVirial(self.densRegion,
                               self.regionCoordinates,
                               self.velocityDispersion,
                               self.time)
    def SetId(self,id):
        self.id = id









# Working script
########################################################################################################################################
if __name__ == "__main__":

    # All code parameters
    ####################################################################################################################################
    dx = dy = dz        = 0.02  # pc
    dt                  = 10    # kyr
    Amin                = 25    # dxdy
    cs                  = 0.28  # km /s, sound-speed for the neutral gas
    readDataDir         = "./testData/"
    readMaskDir         = "./maskData/"
    writeDir            = "./Data/"
    sThreshold          = 0.5   # density threshold
    globaluniqueID      = {}    # initalise the global ID dictionary
    times               = np.arange(50,52)#args['tend'])     # the times to run the code on
    timesArr            = []            # a list for storing the different times for plotting
    centerX             = []            # the x coordiante of the centroid
    centerY             = []            # the y coordinate of the centroid
    tIter               = 0             # the time iteration value
    minDisTol           = 2             # the tolerance in Euc. distance for tracking centroids across time (in pixel size)
    idsPerTimeStep      = {}            # the dictionary for storing the IDs each time step
    allIDs              = []            # a list of all IDs, for all time
    openFactor          = 2             # the openning factor of the pixels
    windowSize          = 50            # half the size of the window
    statsPerTimeStep    = {}            # a dictionary that stores the statistics for each time step
                                        # and then updates the global dictionary
    sliceIndexes        = np.arange(200)# the range of slice indexes for the 200^3 dataset
    sliceIndex          = 0             # the index of which z plane to slice the fields through

    ####################################################################################################################################


    # Clear old plots from the directory
    if args['clear'] == True:
        print("Clearing old plots.")
        os.system("rm ./Plots/*")


    # for each piece of data
    for time in times:



        # for each slice through the 3D field
        #for sliceIndex in sliceIndexes:

        print(f"Starting iteration on file number {time}")
        print("#####################################")


        # read in the Data and extract into a np.array
        dens        = loadObj(readDataDir + f"rho_{time}")
        mask        = np.load(readMaskDir + f"mask3D_{time}.npy")
        mask = mask[:,sliceIndex,:]

        # take a slice through (x,y,z=0)
        dens    = dens[:,sliceIndex,:] # g / cm^3
        s       = np.log(dens / dens.mean())


        # include the velocity information if the velocity argument
        # is true
        if args['vel'] == True:
            v = computeVelocity(dens,time,sliceIndex) * 1e-5 # km / s


        # just a quick check of the field
        if args['viz'] == "field":
            f, ax = plt.subplots(dpi=200)
            ax.imshow(s,cmap=plt.cm.plasma)
            ax.set_axis_off()
            plt.show()


        # filter on the densities within the small regions and then create a threshold
        #sML         = s[:,x0:x1]                    # no longer required
        #sML_mask    = sML > sML.max()*sThreshold    # no longer required

        f, ax   = plt.subplots(1,2,figsize=(9,4),dpi=200)
        plt.subplots_adjust(left=0.00, bottom=0.05, right=0.95, top=0.95, wspace=-0.05, hspace=0.05)
        ax[0].imshow(mask,cmap=plt.cm.plasma)
        #rect = patches.Rectangle((x0,0),x1 - x0, s.shape[0],linewidth=1,edgecolor='r',facecolor='r',alpha=0.2);
        #ax[0].plot((x0 + (xMax - xMin)/2.,x0 + (xMax - xMin)/2.),(0,s.shape[0]), '--r',linewidth=0.5)
        ax[0].annotate(r"$\xi = s > s_{\text{max}}/2$", xy=(labelDx,labelDy), xycoords = xyCoords,color="yellow",fontsize=fs-2)
        ax[0].annotate(r"$(\xi \ominus \square ) \oplus \square = $" + " {}".format(openFactor*dx) + r"$\,\text{pc}$", xy=(labelDx,labelDy-0.06), xycoords = 'axes fraction',color="yellow",fontsize=fs-2)
        ax[0].annotate(r"Hough fit", xy=(labelDx-0.45,labelDy), xycoords = xyCoords,color="blue",fontsize=fs-2)
        ax[0].annotate(r"$\mathcal{W}_{\Delta x} = $" + " {}".format(np.round(windowSize*2*dx,1)) + r"$\,\text{pc}$", xy=(labelDx-0.45,labelDy-0.06), xycoords = xyCoords,color="red",fontsize=fs-2)
        ax[0].annotate(r"$\mathcal{A} \geq \mathcal{A}_{\text{min}} = $" + " {}".format(np.round(Amin*dx*dy,3)) + r"$\,\text{pc}^2$", xy=(labelDx,0.03), xycoords = xyCoords,color="yellow",fontsize=fs-2)
        #ax[0].add_patch(rect)
        #ax[0].plot(x,y, '-b')
        ax[0].set_ylim((s.shape[0], 0))
        ax[0].set_axis_off()

        #sML_mask    = np.pad(sML_mask,((0,0),(x0,s.shape[1]-x1)),mode='constant') #?
        allLabels   = labelCreator(mask,openFactor)

        # create the s (slice) map with time annotations
        plot    = ax[1].imshow(s,cmap=plt.cm.plasma)#,vmin=vMin,vmax=vMax)
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
        print("Labelling pillars.")
        for region in measure.regionprops(allLabels):


            # skip small regions
            if region.area <= Amin:
                continue


            # draw rectangle around segmented high-density regions
            minr, minc, maxr, maxc = region.bbox
            rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,fill=False, edgecolor='red', linewidth=2)
            ax[1].add_patch(rect)

            ############################################################
            # calculate everything from the region (i.e. this is where
            # all of the region statistics are)
            ############################################################

            #TODO: reshape 2D array into a 3D array with dimensions 1 in the z direction,
            """
            > define minz:maxz that is default 0:0 for the 2D problem
            >
            """

            ## DENSITY REGIONS ##
            # define the the s / dens region
            regionCoords    = (minr,maxr,minc,maxc)
            sRegion         = s[minr:maxr,minc:maxc]
            densRegion      = dens[minr:maxr,minc:maxc]


            ## DENSITY DISPERSION ##
            # if args['ionX'] is true then filter the densities by only using the
            # neutrals to calculate the dispersion
            if args['ionX'] == True:
                neutralCoords, ionisationFraction   = ionFractionFilter(time,regionCoords,sliceIndex)
                densDispersion                      = np.var(sRegion[neutralCoords])
            else:
                densDispersion  = np.var(sRegion)

            #BUG TEST#
            if np.isnan(densDispersion) == True:
                print("Density dispersion is nan")
                continue


            ## MASS ##
            # calculate the total mass within the region rho
            mass            = sum(sum(densRegion))*dx*dy*dz*(Parsec)**3 * (1. / SolarMass)

            ## AREA ##
            # area in parsecs^2.
            area            = region.area*dx*dy

            ## PERIMETER ##
            # perimeter in parsecs
            per             = region.perimeter*dx


            # calculate the velocity statistics
            if args['vel'] == True:


                # select the region in the mag v field either
                # using the regular bounding box or the neutral
                # coorindates
                if args['ionX'] == True:
                    regionVel   = v[neutralCoords]
                else:
                    regionVel   = v[minr:maxr,minc:maxc]

                ## VELOCITY DISPERSION ##
                # calculate the velocity dispersion
                velocityDispersion  = regionVel.std()
                regionMach          = velocityDispersion / cs


                ## TURBULENT DRIVING MODE ##
                # calculate the forcing parameter, b
                regionB     = np.sqrt( ( np.exp(densDispersion) - 1 ) / regionMach**2  )

                #BUG TEST#
                if np.isnan(regionB) == True:
                    print("b is nan")
                    continue


            ## CENTROID COORDINATES ##
            # calculate the (mass) centroids
            """
            add a centroid z
            """
            centroidY ,centroidX = region.centroid


            # store the centroids and the time for plotting
            centerX.append(centroidX)
            centerY.append(centroidY)
            ax[1].scatter(centerX,centerY,c='b',s=1,marker='.')

            # initialise an array for storing the Euclidean distances between
            # successive time-steps
            euclideanDis = []


            # initialise an on / off state for adding new keys
            addState = 0


            # if we are on the first iteration just store all of the regions into a dictionary
            # and give a unique ID
            if tIter == 0:
                globaluniqueID[regionCounter]  = {'x':centroidX,'y':centroidY}

                # add to the local dictionary the various interesting pillarStatistics
                # one for if the
                if args['vel'] == True:
                    statsPerID[regionCounter]  = {'mass':mass,
                                                  'svar':densDispersion,
                                                  'area':area,
                                                  'per':per,
                                                  'mach':regionMach,
                                                  'b':regionB,
                                                  'fracX':ionisationFraction}
                else:
                    statsPerID[regionCounter]  = {'mass':mass,
                                                  'svar':densDispersion,
                                                  'area':area,
                                                  'per':per,
                                                  'fracX':ionisationFraction}

                # add a number annotation to each of the centroids
                ax[1].text(centroidX,centroidY,str(regionCounter),fontsize=16)


            # if not the first iteration (time)
            else:


                # for each of the possible keys (reminder we are looking at a single region)
                # in the previous time step
                for key in possibleKeys.keys():


                    # get the centroid coordinates from the previous time-step
                    centroidXOld = globaluniqueID[key]['x']
                    centroidYOld = globaluniqueID[key]['y']


                    # calculate the Euclidean distance between each centroid pair
                    euclideanDis.append(np.hypot( centroidX - centroidXOld,centroidY - centroidYOld))


                # extract all of the keys from the dictionary as integers
                keysAsInts  = list(map(int,possibleKeys.keys()))


                # I am popping keys out of the dictionary as I find them... so
                # if all of the possible regions are still there then move to the next iteration
                if keysAsInts == []:
                    print("All pillars were found from the previous time-step.")
                    continue # to the next pillar


                # Calculate the minimum distance between pillars between two
                # concurrent time-steps

                """
                TODO:
                Add mass constraints ?
                Add mergers and splits

                """

                # if there is just one close (below minDisTol) pillar in the next time-step
                #if np.sum(np.array(euclideanDis) < minDisTol) == 1:
                minDis      = min(np.array(euclideanDis))
                minKey      = keysAsInts[euclideanDis.index(min(euclideanDis))]
                # if there are more lose pillars
                #else:
                    # need the mass of the newly labelled pillars



                # if it is less than some minimum distance tolerance
                if minDis < minDisTol:


                    # if it is, then store the value of that centroid using the old key
                    # which is the "tracking" part of this code
                    globaluniqueID[minKey]  = {'x':centroidX,'y':centroidY}


                    # update the dictionary for single ID
                    if args['vel'] == True:
                        statsPerID[minKey]      = {'mass':mass,
                                                   'svar':densDispersion,
                                                   'area':area,
                                                   'per':per,
                                                   'mach':regionMach,
                                                   'b':regionB,
                                                   'fracX':ionisationFraction}
                    else:
                        statsPerID[minKey]      = {'mass':mass,
                                                   'svar':densDispersion,
                                                   'area':area,
                                                   'per':per,
                                                   'fracX':ionisationFraction}


                    # add a number annotation to each of the centroids
                    ax[1].text(centroidX,centroidY,minKey,fontsize=16)

                    # get rid of the key that we found (i.e.) it is no longer a
                    # possible key for the rest of the regions to be
                    possibleKeys.pop(minKey,None)


                    # get rid of it from the local dictionary too.
                    localuniqueID[minKey] = None


                    # update an "addState", which really means that no candidates have been found
                    # and we need to add a new key
                else:
                    addState = 1


                # if you get to the end of the keys with nothing satisfying then we
                # will need to add a new key
                # a switch that makes a new ID be incorporated into a set.
                if addState == 1:


                    # create a new key that is one larger than the previous max
                    newKey                 = max(list(map(int,globaluniqueID.keys()))) + 1


                    # add it to the global centroid dictionary
                    globaluniqueID[newKey] = {'x':centroidX,'y':centroidY}

                    # update the dictionary for single ID
                    if args['vel'] == True:
                        statsPerID[newKey]      = {'mass':mass,
                                                   'svar':densDispersion,
                                                   'area':area,
                                                   'per':per,
                                                   'mach':regionMach,
                                                   'b':regionB,
                                                   'fracX':ionisationFraction}
                    else:
                        statsPerID[newKey]     = {'mass':mass,
                                                  'svar':densDispersion,
                                                   'area':area,
                                                   'per':per,
                                                   'fracX':ionisationFraction}


                    # add a number annotation to each of the centroids
                    ax[1].text(centroidX,centroidY,newKey,fontsize=16)


            # update the region counter
            regionCounter +=1


        # print how many distinct regions there are
        print("There have been {} distinct structures identified.".format(regionCounter))


        # a counter for recording the elements that are being deleted
        delCounter = 0


        # Now remove keys that did NOT get popped from the possible keys
        # so they never they never satisfied the minimum Euclidean distance condition
        difSet = set(list(map(int,possibleKeys.keys())))
        if difSet != set():


            # update the delete counter for output
            delCounter += 1


            # update the global dictionary, getting rid of keys that correspond to
            # pillars that have died
            for element in difSet:

                globaluniqueID.pop(element,None)

        print("There has been {} killed structures.".format(delCounter))


        # Need to copy the dictionary otherwise it updates per time-step
        updateGlobal            = globaluniqueID.copy()
        idsPerTimeStep[time]    = updateGlobal


        # Need to copy the dictionary otherwise it updates per time-step
        updateStats             = statsPerID.copy()
        allIDs                  += map(int,globaluniqueID.keys())
        statsPerTimeStep[time]  = updateStats


        # final plotting configuration
        plt.tight_layout()
        print("Writing rho_{}_{}.png".format(time,sliceIndex))
        plt.savefig("Plots/rho_{}_{}.png".format(time,sliceIndex))
        plt.close()
        print("Iteration {} complete on slice {}.".format(tIter,sliceIndex))
        print("##################################### \n\n")

        # update the time step
        tIter += 1


# Store ALL of the IDs (for lifetime plotting)
idsPerTimeStep["all"] = allIDs


# Write out all of the statistics
if args['write'] == True:
    print("All iterations done. \n Writing the pillar statistics.")
    saveObj(idsPerTimeStep,writeDir + "pillarIDs")
    saveObj(statsPerTimeStep,writeDir + "pillarStatistics")
