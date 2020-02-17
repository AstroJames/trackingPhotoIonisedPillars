#!/usr/bin/env python2

""""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          this file unloads the dictionary pickles developed in main and computes statistics
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  17 / Jan / 2020

"""

import py_compile
py_compile.compile("header.py")
from header import *
py_compile.compile("dataHandling.py")
from dataHandling import *
import joypy

# Functions
###################################################################################################

def get_cmap(n, name='hsv'):
    '''
    Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.


    '''
    return plt.cm.get_cmap(name, n)



if __name__ == "__main__":

    # pillar ID Statistics
    ###################################################################################################

    pillarIDDic = loadObj("pillarIDs")
    allIDs      = pillarIDDic['all']
    pillarIDDic.pop('all',None)

    # Histogram of lifetimes
    hist, _         = np.histogram(allIDs, bins=max(allIDs))
    pillarLifeTime  = hist * 10e3 #yrs
    plt.figure(dpi=200)
    plt.hist(x=pillarLifeTime,bins=20,normed=True,log=True,color='black',histtype='step')
    plt.xscale('log')
    plt.xlabel(r'$t \, [\text{yrs}]$')
    plt.ylabel(r'$\text{PDF}$')
    plt.tight_layout()
    plt.show()

    # Space-time diagrams for each of the pillars
    colorMap = get_cmap(300, name='flag')
    f, ax = plt.subplots(2,1,dpi=200,sharex=True)
    plt.subplots_adjust(left=0.08, bottom=0.12, right=0.97, top=0.97, wspace=-0.05, hspace=0.05)
    for timeKey in np.arange(110,200):
        for idKey in pillarIDDic[timeKey]:
            ax[0].scatter(timeKey,pillarIDDic[timeKey][idKey]['x']*0.02,
                        c=colorMap(idKey),s=6,edgecolors='k',linewidths=0.1)
            ax[1].scatter(timeKey,pillarIDDic[timeKey][idKey]['y']*0.02,
                        c=colorMap(idKey),s=6,edgecolors='k',linewidths=0.1)
    ax[0].set_ylabel(r"$x \,[\text{pc}]$")
    ax[1].set_ylabel(r"$y \,[\text{pc}]$")
    ax[1].set_xlabel(r"$t \,[\Delta t]$")
    plt.show()


    # pillar observable Statistics
    ###################################################################################################

    pillarStats= loadObj("pillarStatistics")
