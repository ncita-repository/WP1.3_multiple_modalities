# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:17:23 2020

@author: ctorti
"""



""" Find the nearest point pairs between two contours """


def IndsOfNearestPoints(Contour1, Contour2):
    # Initialise the indices of the points in Points2 that 
    # are closest to each point in Points1:
    NearestInds = []
        
    # Loop through every point in the Contour1:
    for i in range(len(Contour1)):
        point_i = Contour1[i]
        
        # The distances between point i in Contour1 and every 
        # point in Contour2:
        distances = []
        
        # Loop through every point in Contour2:
        for j in range(len(Contour2)):
            point_j = Contour2[j]
            
            distance = ( (point_i[0] - point_j[0])**2 +  (point_i[1] - point_j[1])**2 + (point_i[2] - point_j[2])**2 )**(1/2)
            
            distances.append(distance)
            
        # The index of the nearest point in Contour2:
        NearestInds.append(distances.index(min(distances)))
        
    
    #print('\nNearestInds =', NearestInds)
    
    return NearestInds