# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 11:48:03 2021

@author: ctorti
"""

from importlib import reload

import time
import os
from pathlib import Path
import SimpleITK as sitk
#import itk
import matplotlib.pyplot as plt
#from pydicom import dcmread

import seg_tools.seg_data
reload(seg_tools.seg_data)
import general_tools.general
reload(general_tools.general)

#from general_tools.general import get_unique_items
#from seg_tools.seg_data import (
#    get_seg_data_from_list_of_segs, get_seg_data_from_list_of_labimByRoi
#    )
#from conversion_tools.pixarrs_ims import im_to_pixarr


def plot_metricValues_v_iters(
        metricValues, multiresIters, exportPlot=False, 
        exportDir='cwd', fname=''
        ):
    
    plt.plot(metricValues, 'r')
    plt.plot(multiresIters, 
             [metricValues[ind] for ind in multiresIters], 'b*')
    plt.xlabel('Iteration number', fontsize=12)
    plt.ylabel('Metric value', fontsize=12)
    
    if exportPlot:
        if exportDir == 'cwd':
            exportDir = os.getcwd()
        
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        #exportFname = currentDateTime + '_' + fname + '.jpg'
        exportFname =  f'{fname}_{currentDateTime}.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print(f'\nPlot exported to:\n\n{exportFpath}\n')
        
    return

def plot_fix_mov_ims_v2(
        fixIm, movIm, fixInd, movInd, 
        fixTitle='Fixed image', movTitle='Moving image'
        ):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    fixPixArr = sitk.GetArrayViewFromImage(fixIm)[fixInd,:,:]
    movPixArr = sitk.GetArrayViewFromImage(movIm)[movInd,:,:]
    
    fixZ = fixIm.TransformIndexToPhysicalPoint([0,0,fixInd])[2]
    movZ = movIm.TransformIndexToPhysicalPoint([0,0,movInd])[2]
    
    fixTitle = f'{fixTitle}\nz = {round(fixZ, 2)} mm'
    movTitle = f'{movTitle}\nz = {round(movZ, 2)} mm'
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(fixPixArr, cmap=plt.cm.Greys_r);
    plt.title(fixTitle)
    plt.axis('off')
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(movPixArr, cmap=plt.cm.Greys_r);
    plt.title(movTitle)
    plt.axis('off')
    
    plt.show()
    
    return

def compare_res_results(
        resIm0, resIm1, resInd,
        resTitle0='Resampled image 0', resTitle1='Resampled image 1', 
        exportPlot=False, exportDir='cwd', fname='',
        cBarForDiffIm=False
        ):
    """
    Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space.
    
    Can be used to compare two resampled/registered images, or to compare a
    registered/resampled image to the reference/fixed image.
    """
    
    #useCentre = True # 16/04/21 Not sure it's working as expected
    useCentre = False
        
    resPixarr0 = sitk.GetArrayViewFromImage(resIm0)[resInd,:,:]
    resPixarr1 = sitk.GetArrayViewFromImage(resIm1)[resInd,:,:]
    
    if useCentre:
        # Assuming that resIm0 size = resIm1 size:
        R, C, S = resIm0.GetSize()
        
        i = R//2
        j = C//2
    else:
        i = 0
        j = 0
        
    resZ = resIm0.TransformIndexToPhysicalPoint([i,j,resInd])[2]
    
    resTitle0 = f'{resTitle0}\nk = {resInd}\nz = {round(resZ, 2)} mm'
    resTitle1 = f'{resTitle1}\nk = {resInd}\nz = {round(resZ, 2)} mm'
    
    if exportPlot:
        dpi = 300
    else:
        dpi = 80
        
    plt.subplots(2, 2, figsize=(10,10), dpi=dpi)
    
    alpha = 0.5
    
    blendedIm = (1.0 - alpha)*resIm0[:,:,resInd] + alpha*resIm1[:,:,resInd]
    
    diffIm = resIm0[:,:,resInd] - resIm1[:,:,resInd]
    
    plt.subplot(2, 2, 1)
    plt.imshow(resPixarr0, cmap=plt.cm.Greys_r);
    plt.title(resTitle0)
    plt.axis('off')
    
    plt.subplot(2, 2, 2)
    plt.imshow(resPixarr1, cmap=plt.cm.Greys_r);
    plt.title(resTitle1)
    plt.axis('off')
    
    plt.subplot(2, 2, 3)
    plt.imshow(sitk.GetArrayViewFromImage(blendedIm), cmap=plt.cm.Greys_r);
    plt.title('Blended image')
    plt.axis('off')
    
    ax = plt.subplot(2, 2, 4)
    im = ax.imshow(sitk.GetArrayViewFromImage(diffIm), cmap=plt.cm.Greys_r);
    ax.set_title('Difference image')
    ax.axis('off')
    if cBarForDiffIm:
        #cbar = plt.colorbar(im, ax=ax)
        #cbar.mappable.set_clim(0, MaxVal)
        plt.colorbar(im, ax=ax)
    
    #plt.show()
    
    if exportPlot:
        if exportDir == 'cwd':
            exportDir = os.getcwd()
        
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        #exportFname = currentDateTime + '_' + fname + '.jpg'
        exportFname =  f'{fname}_{currentDateTime}.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print(f'\nPlot exported to:\n\n{exportFpath}\n')
    
    return

