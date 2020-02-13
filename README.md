# Tracking High-Density Pillars
An automated system for tracking high-density, photoionised pillars, previously studied in Menon, Federrath & Kuiper, 2020.

The pillars are created using the PLUTO code, http://plutocode.ph.unito.it .

We rely on computer vision techniques, from `skimage`, like the Hough transform, pixel dilation, and simple segmentation algorithms based upon disconnected density regions to identify each distinct pillar (or high-density globular structure).

Currently the code is implemented in 2D, but the intention is to generalise the code to 3D.

## Statistics we intend to extract from each pillar

* lifetime
* mass
* log density dispersion
* area
* perimeter
* rms turbulent Mach number
* turbulent driving parameter
* star formation potential

## A snap-shot of the detction algorithm working at t= 0
![t=0](/Pics/gitPic1.png)
