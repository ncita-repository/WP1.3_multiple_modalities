# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 07:39:15 2020

@author: ctorti
"""


"""
Function:
    CopyRoi()
    
Purpose:
    Copy a ROI from one DICOM scan (collection) to another.

Input:
    Inds     - list of indeces that may have repeated values

Returns:
    IndsStr  - list of strings with '-1', '-2', '-3', etc. appended to the 
               index for indeces that are repeated
            
"""



def IndecesCount(Inds):
    # Convert Inds to strings:
    IndsStr = [str(Ind) for Ind in Inds]
    
    # Create a list of counts of the number of occurrences of each index in Inds:
    counts = []
    
    # Start the counter at 1:
    n = 1
    
    # For each index in Inds get the number of occurrences of the index:
    for i in range(len(Inds)):
        Ind = Inds[i]
        
        N = Inds.count(Ind)
        
        counts.append(N)
        
        #print(f'\nIndex {Ind} is repeated {N} times.')
        
        # If this index is repeated, append counter to the IndsStr:
        if N > 1:
            IndsStr[i] = IndsStr[i] + '-' + str(n)
            
            n = n + 1 # increment counter
        else:
            n = 1 # reset counter
            
    
    #print('\nIndsStr =', IndsStr)
    
    return IndsStr