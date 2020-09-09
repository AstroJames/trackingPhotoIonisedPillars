#!/usr/bin/env python

"""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          a script for visualising 3D pillar masks
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  23 / April / 2020

"""

from header import *
from DataHandling import *

# Functions
############################################################################################################################################

def plot(mask,iter):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    p1 = ax.voxels(mask)#,cmap=plt.cm.plasma,norm=colors.LogNorm(vmin=s.min(),vmax=s.max()))
    #cb = plt.colorbar(p1,ax=ax)
    #cb.set_label(r"$\ln\left(\rho/\rho_0\right)$",fontsize=fs)
    saveNum = str(iter).zfill(3)
    plt.savefig("test_{}.png".format(saveNum),dpi=200)
    plt.close()

# Main script
############################################################################################################################################

if __name__ == "__main__":
    for iter in range(50,100):
        print("Currently on: {}".format(iter))
        mask = np.load("Data/MaskCubes/mask3D_{}.npy".format(iter))
        #dens = loadObj("testData/rho_{}".format(iter))
        plot(mask.T,iter)
        del mask
