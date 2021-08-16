# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:07:37 2021

@author: ctorti
"""


import SimpleITK as sitk
import matplotlib.pyplot as plt


def plot_two_ims(im0, ind0, plotLabel0, im1, ind1, plotLabel1):
    
    pixarr0 = sitk.GetArrayViewFromImage(im0)[ind0,:,:]
    pixarr1 = sitk.GetArrayViewFromImage(im1)[ind1,:,:]
    
    z0 = im0.TransformIndexToPhysicalPoint([0,0,ind0])[2]
    z1 = im1.TransformIndexToPhysicalPoint([0,0,ind1])[2]
    
    plotLabel0 = f'{plotLabel0}\nz = {round(z0, 2)} mm'
    plotLabel1 = f'{plotLabel1}\nz = {round(z1, 2)} mm'
    
    plt.subplots(1, 2, figsize=(15,8))
    
    plt.subplot(1, 2, 1)
    plt.imshow(pixarr0, cmap=plt.cm.Greys_r);
    plt.title(plotLabel0)
    plt.axis('off')
    
    plt.subplot(1, 2, 2)
    plt.imshow(pixarr1, cmap=plt.cm.Greys_r);
    plt.title(plotLabel1)
    plt.axis('off')
    
    plt.show()
    
    return

def plot_pixarrBySeg(pixarrBySeg, f2sIndsBySeg, plotTitle=''):
    """ 
    pixarrBySeg can either be a list of pixel arrays (as the variable name
    suggests), or a pixel array, and f2sIndsBySeg can either be a list (for 
    each segment) of a list (for each frame) of frame-to-slice indices, or it
    can be a list (for each frame) of frame-to-slice indices.
    
    The pixel array(s) must have shape (F, R, C) where F is the number of
    frames, which could be any integer number including 0. 
    """
    
    # Set the number of subplot rows and columns:
    if isinstance(pixarrBySeg, list):
        Ncols = len(pixarrBySeg)
        
        NFramesBySeg = []
        
        for i in range(Ncols):
            pixarr = pixarrBySeg[i]
            
            NFramesBySeg.append(pixarr.shape[0])
            
        Nrows = max(NFramesBySeg)
    else:
        Nrows = pixarr.shape[0]
                
        """ Put the pixel array and f2sInds into a list. """
        pixarrBySeg = [pixarrBySeg]
        f2sIndsBySeg = [f2sIndsBySeg]
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows))
    
    n = 1 # initialised sub-plot number
        
    # Loop through each pixel array:
    for i in range(Ncols):
        pixarr = pixarrBySeg[i]
        f2sInds = f2sIndsBySeg[i]
        
        print(f'\npixarr {i} has shape {pixarr.shape}')
        print(f'f2sInds = {f2sInds}')
        F = pixarr.shape[0]
        
        if F != len(f2sInds):
            print(f'The number of frames in pixarr, {F}, does not match the',
                  f'length of f2sInds, {len(f2sInds)}')
        
        # Loop through each frame:
        for f in range(F):
            #if f < F:
            ax = plt.subplot(Nrows, Ncols, n)
            ax.imshow(pixarr[f], cmap=plt.cm.Greys_r)
            ax.set_xlabel('pixels'); ax.set_ylabel('pixels')
            ax.set_title(plotTitle + f'\nFrame {f2sInds[f]}')
            
            n += 1 # increment sub-plot number
    return