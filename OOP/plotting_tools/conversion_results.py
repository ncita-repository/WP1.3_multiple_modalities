# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 11:46:07 2021

@author: ctorti
"""

from importlib import reload

import time
import os
from pathlib import Path
#import SimpleITK as sitk
#import itk
import matplotlib.pyplot as plt
#from pydicom import dcmread

import seg_tools.seg_data
reload(seg_tools.seg_data)
import general_tools.general
reload(general_tools.general)

from general_tools.general import get_unique_items, unpack
from image_tools.attrs_info import get_im_attrs_from_list_of_dicomDir
#from seg_tools.seg_data import (
#    get_seg_data_from_list_of_segs, get_seg_data_from_list_of_labimByRoi
#    )
#from conversion_tools.pixarrs_ims import im_to_pixarr
from conversion_tools.inds_pts_cntdata import pts_to_inds


def plot_mask_to_cnt_conversion(
        listOfIms, listOfDcmPixarr, listOfF2SindsByRoi, 
        listOfPtsByCntByRoi, listOfPixarrByRoi,
        listOfDcmDirs, listOfPlotTitles,
        exportPlot=False, exportDir=None, exportFname='', p2c=False):
    """ 
    08/06/21: Modelled on PlotPixArrsFromListOfLabImBySeg and 
    PlotContoursFromListOfRtss_v4.
    
    listOfIms is a list (for each dataset, which could be of length 1) of
    SimpleITK Images.
    
    listOfDcmPixarr is a list (for each dataset, which could be of length 1) of
    3D pixel array representations of the DICOM series.
    
    listOfF2SindsByRoi is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI/segment) of a list of the contour-/segmentation-
    to-slice indices in the contour/segmentation data.
    
    listOfPtsByCntByRoi is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI) of a list (for each contour) of points.
    
    listOfPixarrByRoi is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI) of a Numpy pixel array, converted from the
    contour data.
    
    listOfDcmDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each image.
    
    listOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    print(f'len(listOfIms) = {len(listOfIms)}')
    print(f'len(listOfF2SindsByRoi) = {len(listOfF2SindsByRoi)}')
    print(f'len(listOfPtsByCntByRoi) = {len(listOfPtsByCntByRoi)}')
    print(f'len(listOfPixarrByRoi) = {len(listOfPixarrByRoi)}')
    print(f'len(listOfDcmDirs) = {len(listOfDcmDirs)}')
    #print(f'listOfF2SindsByRoi = {listOfF2SindsByRoi}')
    
    listOfDcmFpaths, listOfSizes, listOfSpacings, listOfSlcThick, listOfIPPs,\
        listOfDirections = get_im_attrs_from_list_of_dicomDir(listOfDcmDirs)
        
    # Get the maximum number of slices containing contours/segmentations in any
    # ROI/segment in any dataset:
    maxNumSlices = 0
    
    # The list (for each dataset) of the unique set of slice numbers that 
    # correspond to all contours/segmentations in all ROIs/segments:
    listOfUniqueSinds = []
    
    # Loop through each dataset:
    for i in range(len(listOfIms)):
        f2sIndsByRoi = listOfF2SindsByRoi[i]
        
        uniqueSinds = get_unique_items(f2sIndsByRoi)
        
        listOfUniqueSinds.append(uniqueSinds)
        
        if len(uniqueSinds) > maxNumSlices:
            maxNumSlices = len(uniqueSinds)
    
    print(f'uniqueSinds = {uniqueSinds}')
    print(f'listOfUniqueSinds = {listOfUniqueSinds}')
    print(f'maxNumSlices = {maxNumSlices}')
    
    
    """ Prepare the figure. """
    
    # Set the number of subplot rows and columns:
    Ncols = len(listOfIms)
    Nrows = maxNumSlices
    
    cMaps = ['Reds', 'Blues', 'Greens', 'Oranges', 'Purples']
    
    lineWidth = 0.5
    lineWidth = 2
    colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
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
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(7*Ncols, 9*Nrows), 
                               dpi=dpi)
    else:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(4*Ncols, 8*Nrows), 
                               dpi=dpi)
    
    # Loop through each slice/frame number (i.e. each row in the plot):
    for rowNum in range(maxNumSlices):
        # Loop through each data set:
        for i in range(Ncols):
            dcmIm = listOfIms[i]
            dcmPixarr = listOfDcmPixarr[i]
            ptsByCntByRoi = listOfPtsByCntByRoi[i]
            f2sIndsByRoi = listOfF2SindsByRoi[i]
            pixarrBySeg = listOfPixarrByRoi[i]
            uniqueSinds = listOfUniqueSinds[i]
            
            IPPs = listOfIPPs[i]
            #origin = listOfOrigins[i]
            #directions = listOfDirections[i]
            #spacings = listOfSpacings[i]
            plotTitle = listOfPlotTitles[i]
            #dcmDir = listOfDcmDirs[i]
            
            if p2c:
                print(f'\n   i = {i}')
                print(f'   plotTitle = {plotTitle}')
                print(f'   f2sIndsByRoi = {f2sIndsByRoi}')
                print(f'   uniqueSinds = {uniqueSinds}')
            
            # Continue if rowNum does not exceed the number of unique slice indices:
            if rowNum < len(uniqueSinds):
                sInd = uniqueSinds[rowNum]
                
                dcmFrame = dcmPixarr[sInd]
                IPP = IPPs[sInd]
                
                ax = plt.subplot(Nrows, Ncols, n)
            
                #im = ax.imshow(dcmFrame, cmap=plt.cm.Greys_r)
                ax.imshow(dcmFrame, cmap=plt.cm.Greys_r, alpha=dcmAlpha)
                
                if p2c:
                    print(f'   sInd = {sInd}')
                    print(f'   IPP = {IPP}')
                
                # Loop through each segment:
                for s in range(len(pixarrBySeg)):
                    f2sInds = f2sIndsByRoi[s]
                    segPixarr = pixarrBySeg[s]
                    
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
                        
                            segFrame = segPixarr[frameNum]
                            
                            if p2c:
                                print(f'      sliceNum = {sInd} is in f2sInds')
                                print(f'      frameNum = {frameNum}')
                                print(f'      segFrame.shape = {segFrame.shape}')
                        
                            #ax = plt.subplot(Nrows, Ncols, n)
                            #im = ax.imshow(segFrame, cmap=plt.cm.nipy_spectral, 
                            #               alpha=alpha)
                            ax.imshow(segFrame, cmap=cMap, alpha=segAlpha)
                            
                            #frameTxt = f'frame {frameNum}'
                        #if len(frameNums) > 1:
                        #    frameTxt = f'Conversion to segmentations {frameNums}'
                        #else:
                        #    frameTxt = f'Conversion to segmentations {frameNum}'
                        #frameTxt = f'Contour-to-segmentation # {frameNums}'
                        nobrackets = f'{frameNums}'.replace('[', '').replace(']', '')
                        if i == 0:
                            frameTxt = f'Contour-to-segmentation # {nobrackets}'
                        else:
                            frameTxt = f'Segmentation # {nobrackets}'
                    else:
                        frameTxt = 'No Segmentations'
                        
                        if p2c:
                            print(f'      sInd = {sInd} is NOT in f2sInds')
                        
                    #sliceTxt = f'slice {sInd}'
                    #zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    #plotTitle = plotTitle.replace(' ', '\:')
                    #plotTitle = r"$\bf{{{x}}}$".format(x=plotTitle) \
                    #            + f'\n\n{sliceTxt}\n{frameTxt}\n{zPosTxt}'
                    ##plotTitle = r"$\bf{plotTitle}$\," \
                    ##            + f'\n\n{sliceTxt}\n{frameTxt}\n{zPosTxt}'
                    #
                    #ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    #ax.set_title(plotTitle)
                
                # Loop through each ROI:
                for r in range(len(f2sIndsByRoi)):
                    c2sInds = f2sIndsByRoi[r]
                    ptsByCnt = ptsByCntByRoi[r]
                    
                    """ There are only len(colours) colours defined above. Wrap the 
                    ROI index r if there are more ROIs than the number of defined 
                    colours. """
                    if r < len(colours):
                        colour = colours[r]
                    else:
                        m = r//len(colours)
                        
                        colour = colours[r - m*len(colours)]
                
                    if sInd in c2sInds:
                        #contourNum = c2sInds.index(sliceNum) # there may be more 
                        # than one contour for sliceNum!!
                        contourNums = [i for i, e in enumerate(c2sInds) if e==sInd]
                        
                        #print(f'contourNums = {contourNums}')
                        
                        numPtsByCnt = []
                        
                        # Loop through all contour numbers:
                        for c in range(len(contourNums)):
                            colour2 = colours[c] 
                            
                            contourNum = contourNums[c]
                            
                            pts = ptsByCnt[contourNum]
                            
                            numPtsByCnt.append(len(pts))
                            
                            #print(f'\npts = {pts}')
                            
                            inds = pts_to_inds(points=pts, refIm=dcmIm)
                            
                            if p2c:
                                print(f'      sInd = {sInd} is in c2sInds')
                                print(f'      contourNum = {contourNum}')
                                #print(f'      contour = {contour}')
                                print(f'      len(pts) = {len(pts)}')
                            
                            #ax = plt.subplot(Nrows, Ncols, n)
                            
                            # Unpack the contour points' indices:
                            X, Y, Z = unpack(inds)
                            
                            # Plot the contour points:
                            #ax = plt.plot(X, Y, linewidth=0.5, c='r')
                            #ax.plot(X, Y, linewidth=0.5, c='r')
                            #ax.plot(X, Y, linewidth=lineWidth, c=colour)
                            ax.plot(X, Y, linewidth=lineWidth, c=colour2)
                        
                        #if len(contourNums) > 1:
                        #    contourTxt = f'contours {contourNums}'
                        #else:
                        #    contourTxt = f'contour {contourNum}'
                        #contourTxt = f'contour # {contourNums}'
                        nobrackets = f'{contourNums}'.replace('[', '').replace(']', '')
                        if i == 0:
                            contourTxt = f'contour # {nobrackets}'
                        else:
                            contourTxt = f'Segmentation-to-contour # {nobrackets}'
                        nobrackets = f'{numPtsByCnt}'.replace('[', '').replace(']', '')
                        ptsTxt = f'No. of points: {nobrackets}'
                    else:
                        contourTxt = 'No contour'
                        ptsTxt = 'No. of points: 0'
                        
                        if p2c:
                            print(f'      sInd = {sInd} is NOT in c2sInds')
                    
                    sliceTxt = f'slice # {sInd}'
                    zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    if rowNum == 0:
                        plotTitle = plotTitle.replace(' ', '\:')
                        plotTitle = r"$\bf{{{x}}}$".format(x=plotTitle)
                    else:
                        plotTitle = ''
                    plotTitle += f'\n\n{sliceTxt}'
                    if i == 0:
                        plotTitle += f'\n{contourTxt}\n{ptsTxt}\n{frameTxt}'
                    else:
                        plotTitle += f'\n{frameTxt}\n{contourTxt}\n{ptsTxt}'
                    plotTitle += f'\n{zPosTxt}'
                    
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
            
        exportFname = currentDateTime + '_' + exportFname + '.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {exportFpath}\n')
        
    return