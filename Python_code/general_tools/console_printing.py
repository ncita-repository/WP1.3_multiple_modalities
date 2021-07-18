# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:58:18 2021

@author: ctorti
"""

def print_title(Title):
    ul = '*' * len(Title)
    print('\n' + Title)
    print(ul + '\n')
    
    return

def print_inds_by_roi(IndsByRoi):
    """ Print the list of indices in each ROI, e.g. from a list of C2SindsByRoi
    or F2SindsBySeg. """
    
    R = len(IndsByRoi)
    print(f'There are {R} ROIs/segments in IndsByRoi')
    
    C = [len(IndsByRoi[r]) for r in range(R)]
    print(f'There are {C} contours/frames in each ROI/segment')
    
    print('IndsByRoi:')
    
    for r in range(R):
        if r == 0:
            msg = f'   [{IndsByRoi[r]}'
        else:
            msg = f'    {IndsByRoi[r]}'
        
        if r < R - 1:
            msg += ','
        
        if r == R - 1:
            msg += ']'
                
        print(msg)
        
    return

def print_pts_by_cnt(PtsByCnt):
    """ Print the number of points in each contour in PtsByCnt. """
    
    C = len(PtsByCnt)
    print(f'There are {C} contours')
    
    for c in range(C):
        P = len(PtsByCnt[c])
        print(f'   There are {P} points in the {c}^th contour')
    
    return

def print_pts_by_cnt_by_roi(PtsByCntByRoi):
    """ 
    Print the number of points in each contour of each ROI in PtsByCntByRoi. 
    """
    
    R = len(PtsByCntByRoi)
    print('There are two ROIs/segments in IndsByRoi')
    
    #C = [len(PtsByCntByRoi[r]) for r in range(R)]
    #print(f'There are {C} contours/frames in each ROI/segment')
    
    for r in range(R):
        C = len(PtsByCntByRoi[r])
        print(f'There are {C} contours in the {r}^th ROI')
        
        for c in range(C):
            P = len(PtsByCntByRoi[r][c])
            print(f'   There are {P} points in the {c}^th contour')
    
    return

def print_pixarr_by_seg(PixArrBySeg):
    """ 
    If PixArrBySeg is a list of pixel arrays (as the input variable name
    suggests):
        - print the number of pixel arrays in the list
        - print the shape of each pixel array
    
    If PixArrBySeg is a pixel array (not a list as the variable name suggests):
        - print the shape of the pixel array.
    """
    
    if isinstance(PixArrBySeg, list):
        S = len(PixArrBySeg)
        print(f'There are {S} pixel arrays in PixArrBySeg')
        
        for s in range(S):
            print(f'   The {s}^th pixel array has shape {PixArrBySeg[s].shape}')
    else:
        print(f'The pixel array has shape {PixArrBySeg.shape}')
    
    return

def print_pixarr_shape_by_seg(PixArrBySeg):
    #import numpy as np
    
    S = len(PixArrBySeg)
    
    print('\n   Shape of PixArrBySeg:')
    
    for s in range(S):
        if s == 0:
            msg = f'   [{PixArrBySeg[s].shape}'
        else:
            msg = f'    {PixArrBySeg[s].shape}'
        
        if s < S - 1:
            msg += ','
        
        if s == S - 1:
            msg += ']'
                
        print(msg)

    F = [PixArrBySeg[s].shape[0] for s in range(S)]
    
    print(f'   Number of frames in each segment = {F}')
    
    return