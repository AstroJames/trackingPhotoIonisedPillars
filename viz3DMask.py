#!/usr/bin/env python2

""""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          a script for visualising 3D pillar masks
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  23 / April / 2020

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mask3D = np.load("mask3D_50.npy")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.voxels(mask3D, edgecolor='k')
plt.show()