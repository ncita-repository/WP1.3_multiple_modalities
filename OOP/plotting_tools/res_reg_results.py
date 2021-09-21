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

#import dicom_tools.seg_data
#reload(dicom_tools.seg_data)
#import general_tools.general
#reload(general_tools.general)

#from general_tools.general import get_unique_items
#from dicom_tools.seg_data import (
#    get_seg_data_from_list_of_segs, get_seg_data_from_list_of_labimByRoi
#    )
#from conversion_tools.pixarrs_ims import im_to_pixarr

from image_tools.operations import normalise_im
from general_tools.pixarr_ops import checkered_frame

def plot_metricValues_v_iters(
        metricValues, multiresIters, exportPlot=False, 
        exportDir='cwd', fname=''
        ):
    
    plt.figure(figsize=(10,8))
    
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
    
    fixIm = normalise_im(fixIm)
    movIm = normalise_im(movIm)
    
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
        im0, im1, k,
        title0='Resampled image 0', title1='Resampled image 1', 
        exportPlot=False, exportDir='cwd', fname='', fontSize=12,
        cBarForDiffIm=False
        ):
    """
    Plot two SimpleITK images, a blended image and a difference image.
    
    Can be used to compare two resampled/registered images, or to compare a
    registered/resampled image to the reference/fixed image.
    
    Parameters
    ----------
    im0 : SimpleITK Image
    im1 : SimpleITK Image
    k : int
        The slice index for im0 and im1.
    cBarForDiffIm : bool, optional
        If True a colorbar will be included for the difference image. The 
        default value is False.
    title0 : str, optional
        The plot title for im0. The default value is 'Resampled image 0'.
    title1 : str, optional
        The plot title for im0. The default value is 'Resampled image 01'.
    exportPlot : bool, optional
        If True, the plot will be exported to disk. The default value is False.
    exportDir : str, optional
        The directory where the plot will be exported if exportPlot is True. 
        The default value is 'cwd' = the current working directory.
    fname : str, optional
        The file name that will be assigned to the exported plot if exportPlot
        is True. The default value is ''.
    fontSize : int, optional
        The font size for the plot title and axes labels. The default value is 
        12.
        
    Returns
    -------
    None
    
    Note
    ----
    It will be assumed that im0 and im1 have the same size/dimensions, as will
    be the case for resampled/registered images.
    """
    
    # Use the origin or central pixel to report the z-position?
    #useCentre = True # 16/04/21 Not sure it's working as expected
    useCentre = False
    
    # Blending of the two images:
    alpha = 0.5
    
    if useCentre:
        # Assuming that im0 size = im1 size:
        R, C, S = im0.GetSize()
        
        i, j = R//2, C//2
    else:
        i, j = 0, 0
    
    im0 = normalise_im(im0)
    im1 = normalise_im(im1)
    
    frame0 = sitk.GetArrayViewFromImage(im0)[k, :, :]
    frame1 = sitk.GetArrayViewFromImage(im1)[k, :, :]
    
    blendedIm = (1.0 - alpha)*im0[:, :, k] + alpha*im1[:, :, k]
    
    diffIm = im0[:, :, k] - im1[:, :, k]
    
    checkeredFrame = checkered_frame(frame0, frame1)
    
    z = im0.TransformIndexToPhysicalPoint([i, j, k])[2]
    
    title0 = f'{title0}\n\nk = {k}\nz = {z:.2f} mm'
    title1 = f'{title1}\n\nk = {k}\nz = {z:.2f} mm'
    blendedTitle = f'Blended frame\n\nk = {k}\nz = {z:.2f} mm'
    diffTitle = f'Difference frame\n\nk = {k}\nz = {z:.2f} mm'
    checkeredTitle = f'Checkered frame\n\nk = {k}\nz = {z:.2f} mm'
    
    if exportPlot:
        dpi = 300
    else:
        dpi = 80
        
    #plt.subplots(2, 2, figsize=(10,11), dpi=dpi)
    plt.subplots(2, 3, figsize=(15,12), dpi=dpi)
    
    plt.subplot(2, 3, 1)
    plt.imshow(frame0, cmap=plt.cm.Greys_r);
    plt.title(title0, fontsize=fontSize)
    plt.axis('off')
    
    plt.subplot(2, 3, 2)
    plt.imshow(frame1, cmap=plt.cm.Greys_r);
    plt.title(title1, fontsize=fontSize)
    plt.axis('off')
    
    plt.subplot(2, 3, 3)
    plt.axis('off')
    
    plt.subplot(2, 3, 4)
    plt.imshow(sitk.GetArrayViewFromImage(blendedIm), cmap=plt.cm.Greys_r);
    plt.title(blendedTitle, fontsize=fontSize)
    plt.axis('off')
    
    ax = plt.subplot(2, 3, 5)
    im = ax.imshow(sitk.GetArrayViewFromImage(diffIm), cmap=plt.cm.Greys_r);
    ax.set_title(diffTitle, fontsize=fontSize)
    ax.axis('off')
    if cBarForDiffIm:
        #cbar = plt.colorbar(im, ax=ax)
        #cbar.mappable.set_clim(0, MaxVal)
        plt.colorbar(im, ax=ax)
    
    plt.subplot(2, 3, 6)
    plt.imshow(checkeredFrame, cmap=plt.cm.Greys_r);
    plt.title(checkeredTitle, fontsize=fontSize)
    plt.axis('off')
    
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

