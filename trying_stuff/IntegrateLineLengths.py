# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:14:30 2020

@author: ctorti
"""



""" 
Estimate the area of the surface created by linking points from Contour1 to
points in Contour2, starting with the first point in Contour1 and a point in 
Contour2 given by the index StartInd2, and adding up all line lengths.

Note: There are N points in Contour1 and Contour2.
"""

def IntegrateLineLengths(Contour1, Contour2, StartInd2):
    
    # Initialise the cummulative line lengths:
    LineLength = 0
    
    for i in range(len(Contour1)):
        Point1 = Contour1[i]
        
        for j in range(len(Contour2)):
            # The index of the point for Contour2 starts
            # at StartInd2 and increments to N, then is 
            # "wrapped" back to 0 and so on until N 
            # iterations:
            if StartInd2 + j < len(Contour2):
                Ind2 = StartInd2 + j
            else:
                Ind2 = len(Contour2) - (StartInd2 + j)
                
            Point2 = Contour2[Ind2]
            
            # The vector length between Point1 and Point2:
            VectorL = ((Point1[0] - Point2[0])**2 \
                       + (Point1[1] - Point2[1])**2 \
                       + (Point1[2] - Point2[2])**2)**(1/2)
            
            # Add the line length between Point1 and Point2:
            LineLength = LineLength + VectorL
            
            
    return LineLength