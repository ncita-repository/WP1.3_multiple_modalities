# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:45:15 2020

@author: ctorti
"""



""" 
Defining Contour1 and Contour2 as the two contours, re-order the 
points in Contour2 that yields the minimum cummulative line length.
Then the points in the Contour1 and ROContour2, when joined point-by-point, 
will create the minimal surface area (i.e. catenoid) of the two contours:

https://en.wikipedia.org/wiki/Minimal_surface_of_revolution
"""


def ReorderContour(Contour, StartInd):
    
    # Initialise the re-ordered contour:
    ROContour = []
    
    # Re-order Contour by starting at the index StartInd
    # then iterating through all N subsequent points, wrapping
    # the indeces past N:
    for i in range(len(Contour)):
        # The index of the point for Contour starts
        # at StartInd and increments to N, then is 
        # "wrapped" back to 0 and so on until N 
        # iterations:
        if StartInd + i < len(Contour):
            Ind = StartInd + i
        else:
            Ind = len(Contour) - (StartInd + i)

        ROContour.append(Contour[Ind])
        
    return ROContour