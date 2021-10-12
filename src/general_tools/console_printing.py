# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:58:18 2021

@author: ctorti
"""

def print_title(title):
    ul = '*' * len(title)
    print('\n' + title)
    print(ul + '\n')

def print_indsByRoi(indsByRoi):
    """ 
    Print the list of indices in each ROI, e.g. from a list of f2sIndsByRoi.
    """
    
    #print('\n\n', '-'*120)
    #print('Running of print_indsByRoi():')
    
    if indsByRoi:
        R = len(indsByRoi)
        print(f'\nThere are {R} ROIs/segments in indsByRoi')
        
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
                msg += '\n    ]'
            
            print(msg)
    else:
        print('\nNone or empty')
    
    #print('\nEnd of print_indsByRoi\n')
    #print('-'*120)

def print_ptsByCnt(ptsByCnt):
    """ Print the number of points in each contour in ptsByCnt. """
    
    #print('\n\n', '-'*120)
    #print('Running of print_ptsByCnt():')
    
    if ptsByCnt:
        C = len(ptsByCnt)
        print(f'\nThere are {C} contours')
        
        for c in range(C):
            P = len(ptsByCnt[c])
            print(f'   There are {P} points in the {c}^th contour')
    else:
        print('\nNone or empty')
        
    #print('\nEnd of print_ptsByCnt\n')
    #print('-'*120)

def print_ptsByCntByRoi(ptsByCntByRoi):
    """ 
    Print the number of points in each contour of each ROI in ptsByCntByRoi
    or cntdataByCntByRoi. 
    """
    
    #print('\n\n', '-'*120)
    #print('Running of print_ptsByCntByRoi():')
    
    if ptsByCntByRoi:
        R = len(ptsByCntByRoi)
        print(f'\nThere are {R} ROIs/segments in ptsByCntByRoi')
        
        #C = [len(ptsByCntByRoi[r]) for r in range(R)]
        #print(f'There are {C} contours/frames in each ROI/segment')
        
        for r in range(R):
            C = len(ptsByCntByRoi[r])
            print(f'There are {C} contours in the {r}^th ROI')
            
            for c in range(C):
                P = len(ptsByCntByRoi[r][c])
                print(f'   There are {P} points in the {c}^th contour')
    else:
        print('\nNone or empty')
    
    #print('\nEnd of print_ptsByCntByRoi\n')
    #print('-'*120)

def print_pixarrBySeg(pixarrBySeg):
    """ 
    If pixarrBySeg is a list of pixel arrays (as the input variable name
    suggests):
        - print the number of pixel arrays in the list
        - print the shape of each pixel array
    
    If pixarrBySeg is a pixel array (not a list as the variable name suggests):
        - print the shape of the pixel array.
    """
    
    #print('\n\n', '-'*120)
    #print('Running of print_pixarrBySeg():')
    
    if pixarrBySeg:
        if isinstance(pixarrBySeg, list):
            S = len(pixarrBySeg)
            print(f'\nThere are {S} pixel arrays in pixarrBySeg')
            
            for s in range(S):
                print(f'   The {s}^th pixel array has shape {pixarrBySeg[s].shape}')
        else:
            print(f'\nThe pixel array has shape {pixarrBySeg.shape}')
    else:
        print('\nNone or empty')
    
    #print('\nEnd of print_pixarrBySeg\n')
    #print('-'*120)

def print_shape_of_pixarrBySeg(pixarrBySeg):
    #import numpy as np
    
    #print('\n\n', '-'*120)
    #print('Running of print_shape_of_pixarrBySeg():')
    
    if pixarrBySeg:
        S = len(pixarrBySeg)
        
        print('\nShape of pixarrBySeg:')
        
        for s in range(S):
            if s == 0:
                msg = f'[{pixarrBySeg[s].shape}'
            else:
                msg = f' {pixarrBySeg[s].shape}'
            
            if s < S - 1:
                msg += ','
            
            if s == S - 1:
                msg += '\n    ]'
                    
            print(msg)
    
        F = [pixarrBySeg[s].shape[0] for s in range(S)]
        
        print(f'Number of frames in each segment = {F}')
    else:
        print('\nNone or empty')
        
    #print('\nEnd of print_shape_of_pixarrBySeg\n')
    #print('-'*120)

def print_labimBySeg(labimBySeg):
    """ 
    labimBySeg is a list of SimpleITK Images:
        - print the number of images in the list
        - print the size of each image
    """
    
    #print('\n\n', '-'*120)
    #print('Running of print_pixarrBySeg():')
    
    if labimBySeg:
        if isinstance(labimBySeg, list):
            S = len(labimBySeg)
            print(f'\nThere are {S} label images in labimBySeg')
            
            for s in range(S):
                print(f'   The {s}^th label image has size',
                      f'{labimBySeg[s].GetSize()}')
    else:
        print('\nNone or empty')