#!/usr/bin/env python2

""""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          a script for visualising 3D pillar masks
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  23 / April / 2020

"""

# import numpy as np
# import ipyvolume as ipv
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# 
# mask3D = np.load("Data/mask3D_50.npy")
# 
# ipv.quickvolshow(mask3D,opacity=0.05)

import numpy as np
import ipyvolume as ipv
V = np.zeros((128,128,128)) # our 3d array
# outer box
V[30:-30,30:-30,30:-30] = 0.75
V[35:-35,35:-35,35:-35] = 0.0
# inner box
V[50:-50,50:-50,50:-50] = 0.25
V[55:-55,55:-55,55:-55] = 0.0
ipv.quickvolshow(V, level=[0.25, 0.75], opacity=0.03, level_width=0.1, data_min=0, data_max=1)