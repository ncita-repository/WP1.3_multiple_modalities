# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:48:40 2020

@author: ctorti
"""



""" Calculate the perimeter of a contour """

def GetContourPerimeter(Contour):
    
    # Initialise the perimeter value:
    Perimeter = 0
    
    # Loop through each point and integrate the vector lengths:
    for i in range(len(Contour) - 1):
        Point1 = Contour[i]
        Point2 = Contour[i+1]
        
         # The vector length between Point1 and Point2:
        VectorL = (\
                   (Point1[0] - Point2[0])**2 \
                   + (Point1[1] - Point2[1])**2 \
                   + (Point1[2] - Point2[2])**2 \
                   )**(1/2)
        
        Perimeter = Perimeter + VectorL
        
    return Perimeter
