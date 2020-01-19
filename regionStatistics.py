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





if __name__ == "__main__":

    # pillar ID Statistics
    ###################################################################################################

    pillarIDDic = loadObj("pillarIDs")
    allIDs      = pillarIDDic['all']
    pillarIDDic.pop('all',None)

    hist, _         = np.histogram(allIDs, bins=max(allIDs))
    pillarLifeTime  = hist * 10e3 #yrs

    plt.figure(dpi=200)
    plt.hist(x=pillarLifeTime,bins=20,normed=True,log=True,color='black',histtype='step')
    plt.xscale('log')
    plt.xlabel(r'$t \, [\text{yrs}]$')
    plt.ylabel(r'$\text{PDF}$')
    plt.tight_layout()
    plt.show()

    # pillar observable Statistics
    ###################################################################################################

    pillarStats= loadObj("pillarStatistics")
