# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:29:14 2020

@author: ctorti
"""



def NormaliseVector(Vector):
    if len(Vector) == 3:
        from GetVectorLength import GetVectorLength
        
        VectorL = GetVectorLength(Vector)
        
        return [Vector[i]/VectorL for i in range(len(Vector))]
    
    else:
        print('Vector must have length = 3')
        
        return