# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:07:37 2021

@author: ctorti
"""

from importlib import reload

import time
import os
from pathlib import Path
import SimpleITK as sitk
import matplotlib.pyplot as plt
from pydicom import dcmread

import seg_tools.seg_data
reload(seg_tools.seg_data)

from general_tools.general import get_unique_items
from seg_tools.seg_data import get_seg_data_from_list_of_segs
    

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
        
        NframesBySeg = []
        
        for i in range(Ncols):
            pixarr = pixarrBySeg[i]
            
            NframesBySeg.append(pixarr.shape[0])
            
        Nrows = max(NframesBySeg)
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
            ax.set_title(plotTitle + f'\nframe {f2sInds[f]}')
            
            n += 1 # increment sub-plot number
    return

def plot_pixarrs_from_list_of_segs_v3(
        listOfSegs, listOfDicomDirs, listOfPlotTitles,
        exportPlot=False, exportDir='cwd', runID='', useCaseToApply='',
        forceReg=False, useDroForTx=False, regTxName='', initMethod='', 
        resInterp='', txtToAddToFname='', p2c=False
        ):
    """ 
    02/06/2021
    
    Note:
        
    Modification of v1. Rather than plotting 
    allSliceNums = get_unique_items(listOfF2SindsBySeg)
    will just plot frames with segmentations since f2sInds for Src and Trg can 
    relate to very different depths, leading to a confusing plot.
    
    listOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    listOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    listOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    listOfPixarrBySeg, listOfF2SindsBySeg, listOfImSizes, listOfImSpacings,\
        listOfImIPPs, listOfImDirections, listOfDcmFpaths\
            = get_seg_data_from_list_of_segs(listOfSegs, listOfDicomDirs, p2c)
    
    # Get the maximum number of slices containing segmentations in any ROI in 
    # any dataset:
    maxNumSlices = 0
    
    # The list (for each dataset) of the unique set of slice numbers that 
    # correspond to all segmentations in all ROIs:
    listOfUniqueSinds = []
    
    # Loop through each dataset:
    for i in range(len(listOfF2SindsBySeg)):
        f2sIndsBySeg = listOfF2SindsBySeg[i]
        
        uniqueSinds = get_unique_items(f2sIndsBySeg)
        
        listOfUniqueSinds.append(uniqueSinds)
        
        if len(uniqueSinds) > maxNumSlices:
            maxNumSlices = len(uniqueSinds)
    
    
    print(f'listOfF2SindsBySeg = {listOfF2SindsBySeg}')
    print(f'listOfUniqueSinds = {listOfUniqueSinds}')
    print(f'maxNumSlices = {maxNumSlices}')
    
    
    """ Prepare the figure. """
    
    # Set the number of subplot rows and columns:
    Ncols = len(listOfSegs)
    Nrows = maxNumSlices
    
    cMaps = ['Reds', 'Blues', 'Greens', 'Oranges', 'Purples']
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    dcmAlpha = 0.2
    dcmAlpha = 0.5
    #dcmAlpha = 1
    
    segAlpha = 0.5
    #segAlpha = 0.2
    
    if exportPlot:
        dpi = 120
    else:
        dpi = 80
        
    n = 1 # initialised sub-plot number
    
    if Ncols < 3:
        fig, ax = plt.subplots(
            Nrows, Ncols, figsize=(7*Ncols, 9*Nrows), dpi=dpi
            )
    else:
        fig, ax = plt.subplots(
            Nrows, Ncols, figsize=(4*Ncols, 7.5*Nrows), dpi=dpi
            )
    
    # Loop through each slice number (i.e. each row in the plot):
    for rowNum in range(maxNumSlices):
        # Loop through each data set:
        for i in range(len(listOfPixarrBySeg)):
            pixarrBySeg = listOfPixarrBySeg[i]
            f2sIndsBySeg = listOfF2SindsBySeg[i]
            uniqueSinds = listOfUniqueSinds[i]
            
            dicomFpaths = listOfDcmFpaths[i]
            IPPs = listOfImIPPs[i]
            #directions = listOfImDirections[i]
            #spacings = listOfImSpacings[i]
            plotTitle = listOfPlotTitles[i]
            #dicomDir = listOfDicomDirs[i]
            
            if p2c:
                print(f'\n   i = {i}')
                print(f'   plotTitle = {plotTitle}')
                print(f'   f2sIndsBySeg = {f2sIndsBySeg}')
                print(f'   uniqueSinds = {uniqueSinds}')
            
            # Continue if rowNum does not exceed the number of unique slice indeces:
            if rowNum < len(uniqueSinds):
                sInd = uniqueSinds[rowNum]
                
                ax = plt.subplot(Nrows, Ncols, n)
                
                dicomPixarr = dcmread(dicomFpaths[sInd]).pixel_array
            
                #im = ax.imshow(dicomPixarr, cmap=plt.cm.Greys_r)
                ax.imshow(dicomPixarr, cmap=plt.cm.Greys_r, alpha=dcmAlpha)
                
                IPP = IPPs[sInd]
                    
                # Loop through each segment:
                for s in range(len(f2sIndsBySeg)):
                    f2sInds = f2sIndsBySeg[s]
                    pixarr = pixarrBySeg[s]
                    
                    """ There are only len(cMaps) colormaps defined above. Wrap the 
                    segment index s if there are more segments than the number of 
                    defined colormaps. """
                    if s < len(cMaps):
                        cMap = cMaps[s]
                    else:
                        m = s//len(cMaps)
                        
                        cMap = cMaps[s - m*len(cMaps)]
                        
                    if sInd in f2sInds:
                        frameNums = [i for i, e in enumerate(f2sInds) if e==sInd]
                        
                        # Loop through all frame numbers:
                        for f in range(len(frameNums)):
                            frameNum = frameNums[f]
                        
                            frame = pixarr[frameNum]
                            
                            if p2c:
                                print(f'      SliceNum = {sInd} is in f2sInds')
                                print(f'      frameNum = {frameNum}')
                                #print(f'      Contour = {Contour}')
                                print(f'      frame.shape = {frame.shape}')
                        
                            #ax = plt.subplot(Nrows, Ncols, n)
                            #im = ax.imshow(frame, cmap=plt.cm.nipy_spectral, 
                            #               alpha=alpha)
                            ax.imshow(frame, cmap=cMap, alpha=segAlpha)
                            
                            #frameTxt = f'frame {frameNum}'
                        if len(frameNums) > 1:
                            frameTxt = f'frames {frameNums}'
                        else:
                            frameTxt = f'frame {frameNum}'
                            
                    else:
                        frameTxt = 'No frame'
                        
                        if p2c:
                            print(f'      sInd = {sInd} is NOT in f2sInds')
                        
                    sliceTxt = f'Slice {sInd}'
                    zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    plotTitle = plotTitle.replace(' ', '\:')
                    plotTitle = r"$\bf{{{x}}}$".format(x=plotTitle) \
                        + f'\n\n{sliceTxt}\n{frameTxt}\n{zPosTxt}'
                    #plotTitle = r"$\bf{plotTitle}$\," \
                    #            + f'\n\n{sliceTxt}\n{frameTxt}\n{zPosTxt}'
                    
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    ax.set_title(plotTitle)
                    
            n += 1 # increment sub-plot number
        
    
    if exportPlot:
        if exportDir == 'cwd':
            exportDir = os.getcwd()
        
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        #exportFname = currentDateTime + '_' + txtToAddToFname + '.jpg'
        
        exportFname = ''
        if runID:
            exportFname += runID + '_'
        if ('5' in useCaseToApply or forceReg):
            if useDroForTx:
                exportFname += 'useDroForTx_'
            else:
                exportFname += regTxName + '_'
                exportFname += initMethod + '_'
        exportFname += resInterp + '_'
        if txtToAddToFname:
            exportFname += txtToAddToFname + '_'
        exportFname += currentDateTime + '.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {exportFpath}\n')
        
    return