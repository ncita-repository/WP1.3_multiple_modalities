# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:07:37 2021

@author: ctorti
"""


import SimpleITK as sitk
import matplotlib.pyplot as plt


def plot_two_ims(Im0, Ind0, PlotLabel0, Im1, Ind1, PlotLabel1):
    
    PixArr0 = sitk.GetArrayViewFromImage(Im0)[Ind0,:,:]
    PixArr1 = sitk.GetArrayViewFromImage(Im1)[Ind1,:,:]
    
    z0 = Im0.TransformIndexToPhysicalPoint([0,0,Ind0])[2]
    z1 = Im1.TransformIndexToPhysicalPoint([0,0,Ind1])[2]
    
    PlotLabel0 = f'{PlotLabel0}\nz = {round(z0, 2)} mm'
    PlotLabel1 = f'{PlotLabel1}\nz = {round(z1, 2)} mm'
    
    plt.subplots(1, 2, figsize=(15,8))
    
    plt.subplot(1, 2, 1)
    plt.imshow(PixArr0, cmap=plt.cm.Greys_r);
    plt.title(PlotLabel0)
    plt.axis('off')
    
    plt.subplot(1, 2, 2)
    plt.imshow(PixArr1, cmap=plt.cm.Greys_r);
    plt.title(PlotLabel1)
    plt.axis('off')
    
    plt.show()
    
    return

def plot_pixarrBySeg(PixArrBySeg, F2SindsBySeg, PlotTitle=''):
    """ 
    PixArrBySeg can either be a list of pixel arrays (as the variable name
    suggests), or a pixel array, and F2SindsBySeg can either be a list (for 
    each segment) of a list (for each frame) of frame-to-slice indices, or it
    can be a list (for each frame) of frame-to-slice indices.
    
    The pixel array(s) must have shape (F, R, C) where F is the number of
    frames, which could be any integer number including 0. 
    """
    
    # Set the number of subplot rows and columns:
    if isinstance(PixArrBySeg, list):
        Ncols = len(PixArrBySeg)
        
        NFramesBySeg = []
        
        for i in range(Ncols):
            PixArr = PixArrBySeg[i]
            
            NFramesBySeg.append(PixArr.shape[0])
            
        Nrows = max(NFramesBySeg)
    else:
        Nrows = PixArr.shape[0]
                
        """ Put the pixel array and F2Sinds into a list. """
        PixArrBySeg = [PixArrBySeg]
        F2SindsBySeg = [F2SindsBySeg]
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows))
    
    n = 1 # initialised sub-plot number
        
    # Loop through each pixel array:
    for i in range(Ncols):
        PixArr = PixArrBySeg[i]
        F2Sinds = F2SindsBySeg[i]
        
        print(f'\nPixArr {i} has shape {PixArr.shape}')
        print(f'F2Sinds = {F2Sinds}')
        F = PixArr.shape[0]
        
        if F != len(F2Sinds):
            print(f'The number of frames in PixArr, {F}, does not match the',
                  f'length of F2Sinds, {len(F2Sinds)}')
        
        # Loop through each frame:
        for f in range(F):
            #if f < F:
            ax = plt.subplot(Nrows, Ncols, n)
            ax.imshow(PixArr[f], cmap=plt.cm.Greys_r)
            ax.set_xlabel('pixels'); ax.set_ylabel('pixels')
            ax.set_title(PlotTitle + f'\nFrame {F2Sinds[f]}')
            
            n += 1 # increment sub-plot number
    return