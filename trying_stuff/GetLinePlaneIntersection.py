# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:57:01 2020

@author: ctorti
"""




def GetLinePlaneIntersection(PlaneNormal, PlanePoint, 
                            LineDirection, LinePoint):
    """ 
    Get the intersection between a line and a plane.
    
    Inputs:
        PlaneNormal   - A list of the [x, y, z] components of the plane's 
                        normal
                        
        PlanePoint    - A list of the [x, y z] components of a point lying on 
                        the plane

        LineDirection - A list of the [x, y, z] components of the line's
                        vector orientation
                        
        LinePoint     - A list of the [x, y z] components of a point lying on 
                        the point
                        
                        
    Modified from rosettacode.org:
        
    https://rosettacode.org/wiki/Find_the_intersection_of_a_line_with_a_plane
    """
    
    import numpy as np
    
    epsilon=1e-6
    
    PlaneNormal = np.array(PlaneNormal)
    PlanePoint = np.array(PlanePoint)
    LineDirection = np.array(LineDirection)
    LinePoint = np.array(LinePoint)
    
    NormDOTLine = np.dot(PlaneNormal, LineDirection)
    
    if abs(NormDOTLine) < epsilon:
        raise RuntimeError("No intersection or line is within plane")
        
    dPoint = LinePoint - PlanePoint
    
    NormDOTdPoint = - np.dot(PlaneNormal, dPoint) / NormDOTLine
    
    InterPoint = dPoint + NormDOTdPoint * LineDirection + PlanePoint
    
    return list(InterPoint)
    
    
    
    
    