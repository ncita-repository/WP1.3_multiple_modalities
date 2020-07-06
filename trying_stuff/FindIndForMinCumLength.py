# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:27:50 2020

@author: ctorti
"""



""" 
Using the function AddLineLengths(), loop through each N possible 
combinations of starting index StartInd2, and work out which
starting index yields the minimum cummulative line length.
"""

def FindIndForMinCumLength(Contour1, Contour2):
    # Import function:
    from AddLineLengths import AddLineLengths
    
    # Create list to store N cummulative line lengths for all
    # N combinations of starting index for Contour2:
    CumLineLengths = []

    # Loop through all N possible starting indeces for Contour2:
    for i in range(len(Contour2)):
        CumLineLengths.append(AddLineLengths(Contour1=Contour1, Contour2=Contour2, StartInd2=i))

    
    # Return the index that yields the minimum cummulative line length:
    return CumLineLengths.index(min(CumLineLengths))
