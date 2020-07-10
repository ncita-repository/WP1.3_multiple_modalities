# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 13:34:58 2020

@author: ctorti
"""



""" Generate list of cumulative perimeters of a contour at each node """

def GetCumulativePerimeters(Contour):
    
    # Initialise the list:
    CumPlist = []
    
    # Initialise the cumulative perimeter:
    CumP = 0
    
    # Loop through each point:
    for i in range(len(Contour) - 1):
        Point1 = Contour[i]
        Point2 = Contour[i+1]
        
         # The segment length between Point1 and Point2:
        SegmentL = (\
                   (Point1[0] - Point2[0])**2 \
                   + (Point1[1] - Point2[1])**2 \
                   + (Point1[2] - Point2[2])**2 \
                   )**(1/2)
        
        CumP = CumP + SegmentL
        
        CumPlist.append(CumP)
        
    return CumPlist
