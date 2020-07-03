# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:22:39 2020

@author: ctorti
"""



""" Interpolate the points by taking the midway coordinates of each linked point """

def InterpolatePoints(Points1, Points2, Points2IndsNearestToPoint1):

    InterpContour = []
    
    # Loop through every point in Points1:
    for i in range(len(Points1)):
        point_i = Points1[i]
        
        point_j = Points2[Points2IndsNearestToPoint1[i]]
        
        interp_x = (point_i[0] + point_j[0])/2
        
        interp_y = (point_i[1] + point_j[1])/2
        
        interp_z = (point_i[2] + point_j[2])/2
        
        InterpContour.append([interp_x, interp_y, interp_z])
        
    
    #print('\nInterpContour =', InterpContour)
    
    return InterpContour
