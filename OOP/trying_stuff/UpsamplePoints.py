# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 10:35:26 2020

@author: ctorti
"""



""" Up-sample two adjacent points """

def UpsamplePoints(Point1, Point2, N):
    # Import packages:
    import numpy as np
    
    # Create list of up-sampled arrays for x, y and z:
    X = np.linspace(start=Point1[0], stop=Point2[0], num=N)
    Y = np.linspace(start=Point1[1], stop=Point2[1], num=N)
    Z = np.linspace(start=Point1[2], stop=Point2[2], num=N)
    
    # Create a list of [x, y, z] points from each array of X, Y and Z coordinates:
    Points = []
    
    for i in range(N):
        Points.append([X[i], Y[i], Z[i]])
        
        
    #print('Original point pairs:\n\n', Point1, '\n', Point2)
    #print('\nUpsampled points:\n\n', Points)
        
    return Points