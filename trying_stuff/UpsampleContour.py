# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 10:35:49 2020

@author: ctorti
"""



""" Up-sample a contour """

def UpsampleContour(Contour, N):
    # Import function:
    from UpsamplePoints import UpsamplePoints
    
    #print(f'There are {len(OrigContour)} points in the original contour')
          
    # Initialise list of points for up-sampled contour:
    USContour = []
    
    # Iterate over every two points in Contour:
    for i in range(len(Contour) - 1):
        #print(f'i = {i}')
        
        Point1 = Contour[i]
        
        Point2 = Contour[i+1]
        
        Points = UpsamplePoints(Point1=Point1, Point2=Point2, N=N)
        
        USContour.extend(Points)
        
        
    #print(f'There are {len(USContour)} points in the up-sampled contour, '
    #      f'i.e. {round(len(USContour)/len(Contour), 1)} more points.')
    
    return USContour