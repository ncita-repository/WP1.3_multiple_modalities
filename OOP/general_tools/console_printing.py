# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:58:18 2021

@author: ctorti
"""

def print_title(title):
    ul = '*' * len(title)
    print('\n' + title)
    print(ul + '\n')
    
    return

def print_inds_by_roi(indsByRoi):
    """ 
    Print the list of indices in each ROI, e.g. from a list of f2sIindsByRoi.
    """
    
    R = len(indsByRoi)
    print(f'There are {R} ROIs/segments in indsByRoi')
    
    C = [len(indsByRoi[r]) for r in range(R)]
    print(f'There are {C} contours/frames in each ROI/segment')
    
    print('indsByRoi:')
    
    for r in range(R):
        if r == 0:
            msg = f'   [{indsByRoi[r]}'
        else:
            msg = f'    {indsByRoi[r]}'
        
        if r < R - 1:
            msg += ','
        
        if r == R - 1:
            msg += ']'
                
        print(msg)
        
    return

def print_pts_by_cnt(ptsByCnt):
    """ Print the number of points in each contour in ptsByCnt. """
    
    C = len(ptsByCnt)
    print(f'There are {C} contours')
    
    for c in range(C):
        P = len(ptsByCnt[c])
        print(f'   There are {P} points in the {c}^th contour')
    
    return

def print_pts_by_cnt_by_roi(ptsByCntByRoi):
    """ 
    Print the number of points in each contour of each ROI in ptsByCntByRoi. 
    """
    
    R = len(ptsByCntByRoi)
    print('There are two ROIs/segments in indsByRoi')
    
    #C = [len(ptsByCntByRoi[r]) for r in range(R)]
    #print(f'There are {C} contours/frames in each ROI/segment')
    
    for r in range(R):
        C = len(ptsByCntByRoi[r])
        print(f'There are {C} contours in the {r}^th ROI')
        
        for c in range(C):
            P = len(ptsByCntByRoi[r][c])
            print(f'   There are {P} points in the {c}^th contour')
    
    return

def print_pixarr_by_seg(pixarrBySeg):
    """ 
    If pixarrBySeg is a list of pixel arrays (as the input variable name
    suggests):
        - print the number of pixel arrays in the list
        - print the shape of each pixel array
    
    If pixarrBySeg is a pixel array (not a list as the variable name suggests):
        - print the shape of the pixel array.
    """
    
    if isinstance(pixarrBySeg, list):
        S = len(pixarrBySeg)
        print(f'There are {S} pixel arrays in pixarrBySeg')
        
        for s in range(S):
            print(f'   The {s}^th pixel array has shape {pixarrBySeg[s].shape}')
    else:
        print(f'The pixel array has shape {pixarrBySeg.shape}')
    
    return

def print_pixarr_shape_by_seg(pixarrBySeg):
    #import numpy as np
    
    S = len(pixarrBySeg)
    
    print('\n   Shape of pixarrBySeg:')
    
    for s in range(S):
        if s == 0:
            msg = f'   [{pixarrBySeg[s].shape}'
        else:
            msg = f'    {pixarrBySeg[s].shape}'
        
        if s < S - 1:
            msg += ','
        
        if s == S - 1:
            msg += ']'
                
        print(msg)

    F = [pixarrBySeg[s].shape[0] for s in range(S)]
    
    print(f'   Number of frames in each segment = {F}')
    
    return