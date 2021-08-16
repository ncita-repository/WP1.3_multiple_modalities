# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:22:39 2020

@author: ctorti
"""



""" Interpolate contours by taking the midway coordinates of each linked point """

def InterpolateContours(Contour1, Contour2):
    # Import function:
    from IndsOfNearestPoints import IndsOfNearestPoints
    
    # Find the nearest point pairs between two contours:
    NearestIndsContour1to2 = IndsOfNearestPoints(Contour1=Contour1, 
                                                 Contour2=Contour2)

    InterpContour = []
    
    # Loop through every point in Contour1:
    for i in range(len(Contour1)):
        point_i = Contour1[i]
        
        point_j = Contour2[NearestIndsContour1to2[i]]
        
        interp_x = (point_i[0] + point_j[0])/2
        
        interp_y = (point_i[1] + point_j[1])/2
        
        interp_z = (point_i[2] + point_j[2])/2
        
        InterpContour.append([interp_x, interp_y, interp_z])
        
    
    #print('\nInterpContour =', InterpContour)
    
    return InterpContour
