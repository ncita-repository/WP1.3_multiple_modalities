# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:20:30 2020

@author: ctorti
"""



def GetVectorLength(Vector):
    
    if len(Vector) == 3:
        return ( Vector[0]**2 + Vector[1]**2 + Vector[2]**2 )**(1/2)
    
    else:
        print('Vector must have length = 3')
        
        return