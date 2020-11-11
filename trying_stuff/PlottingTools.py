# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:52:36 2020

@author: ctorti
"""



# Import packages and functions:
import os
import time
import pydicom
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

from DicomTools import ImportDicoms
from DicomTools import GetDicomSOPuids

from ImageTools import InitialiseImage
from ImageTools import AddImages

from RtsTools import GetCIStoDcmInds
from RtsTools import ConvertContourDataToIcsPoints
from RtsTools import UnpackPoints

from SegTools import GetPFFGStoDcmInds



"""
******************************************************************************
******************************************************************************
PLOTTING FUNCTIONS
******************************************************************************
******************************************************************************
"""       



def PlotDicomsAndContours(DicomDir, RtsFpath, ExportPlot, ExportDir, 
                          LogToConsole=False):
    
    
    
    # Import the DICOMs:
    Dicoms = ImportDicoms(DicomDir=DicomDir)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir=DicomDir)
    
    # Import the RTS ROI:
    RtsRoi = pydicom.dcmread(RtsFpath)
    
    # Get the ContourImageSequence-to-DICOM slice indices:
    CIStoDcmInds = GetCIStoDcmInds(RtsRoi, SOPuids)
    
    # Get a list of the Contour Sequences in the ROI:
    ContourSequences = RtsRoi.ROIContourSequence[0].ContourSequence
    
    
    # Prepare the figure.
    
    # All DICOM slices that have contours will be plotted. The number of
    # slices to plot:
    Nslices = len(CIStoDcmInds)
    
    # Set the number of subplot rows and columns:
    Ncols = np.int8(np.ceil(np.sqrt(Nslices)))
    
    # Limit the number of rows to 5 otherwise it's difficult to make sense of
    # the images:
    if Ncols > 5:
        Ncols = 5
    
    Nrows = np.int8(np.ceil(Nslices/Ncols))
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows), dpi=300)
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=300)
    else:
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows))
        
    
    n = 0 # sub-plot number
    
    # Loop through each slice: 
    for i in range(Nslices):
        if LogToConsole:
                print('\ni =', i)
                
        n += 1 # increment sub-plot number
        
        # The DICOM slice number is:
        s = CIStoDcmInds[i]
        
        #plt.subplot(Nrows, Ncols, n)
        #plt.axis('off')
        ax = plt.subplot(Nrows, Ncols, n)

        #ax = plt.subplot(Nrows, Ncols, n, aspect=AR)
          
        # Plot the DICOM pixel array:
        ax.imshow(Dicoms[s].pixel_array, cmap=plt.cm.Greys_r)
        #ax.set_aspect(ar)
        
        # Get the Contour Sequence that corresponds to this slice:
        ContourSequence = ContourSequences[i]
        
        # Get the Contour Data:
        ContourData = [float(item) for item in ContourSequence.ContourData]
                
        # Convert contour data from Patient Coordinate System to Image 
        # Coordinate System:
        PointsICS = ConvertContourDataToIcsPoints(ContourData, DicomDir)
        
        # Unpack the points:
        X, Y, Z = UnpackPoints(PointsICS)
        
        # Plot the contour points:
        plt.plot(X, Y, linewidth=0.5, c='red');
        
        
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        
        ax.set_title(f'Slice {s}')
        #plt.axis('off')
        
            
            
    if ExportPlot:
        RtsFullFname = os.path.split(RtsFpath)[1]
        
        RtsFname = os.path.splitext(RtsFullFname)[0]
        
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
        
        ExportFname = CurrentDateTime + '_' + RtsFname
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
        
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
        
    return





def PlotDicomsAndSegments(DicomDir, SegFpath, ExportPlot, ExportDir, 
                          LogToConsole):
    
    
    
    # Import the DICOMs:
    Dicoms = ImportDicoms(DicomDir=DicomDir)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir=DicomDir)
    
    # Import the SEG ROI:
    SegRoi = pydicom.dcmread(SegFpath)
    
    # Get the Per-frameFunctionalGroupsSequence-to-DICOM slice indices:
    PFFGStoDcmInds = GetPFFGStoDcmInds(SegRoi, SOPuids)
    
    # The 3D labelmap:
    SegNpa = SegRoi.pixel_array
    
    if LogToConsole:
        print(f'Segments exist on slices {PFFGStoDcmInds}')
        print(f'Shape of PixelData = {SegNpa.shape}\n')

    
    
    # Prepare the figure.
    
    # All DICOM slices that have segmentations will be plotted. The number of
    # slices to plot (= number of segments):
    #Nslices = int(SegRoi.NumberOfFrames)
    Nslices = len(PFFGStoDcmInds)
    
    # Set the number of subplot rows and columns:
    Ncols = np.int8(np.ceil(np.sqrt(Nslices)))
    
    # Limit the number of rows to 5 otherwise it's difficult to make sense of
    # the images:
    if Ncols > 5:
        Ncols = 5
    
    Nrows = np.int8(np.ceil(Nslices/Ncols))
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows), dpi=300)
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=300)
    else:
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows))
        
    
    n = 1 # initialised sub-plot number
    
    # Loop through each slice: 
    for i in range(Nslices):
        if LogToConsole:
                print('\ni =', i)
                
        # The DICOM slice number is:
        s = PFFGStoDcmInds[i]
        
        #plt.subplot(Nrows, Ncols, n)
        #plt.axis('off')
        ax = plt.subplot(Nrows, Ncols, n)

        #ax = plt.subplot(Nrows, Ncols, n, aspect=AR)
          
        # Plot the DICOM pixel array:
        ax.imshow(Dicoms[s].pixel_array, cmap=plt.cm.Greys_r)
        #ax.set_aspect(ar)
        
        if LogToConsole:
            print(f'\nShape of Dicoms[{s}] = {Dicoms[s].pixel_array.shape}')
            print(f'\nShape of SegNpa[{i}] = {SegNpa[i].shape}')
        
        # Plot the segment pixel array:
        #ax.imshow(SegNpa[i], cmap='red', alpha=0.1)
        ax.imshow(SegNpa[i], alpha=0.2)
        
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        
        ax.set_title(f'Slice {s}')
        #plt.axis('off')
        
        n += 1 # increment sub-plot number
        
            
            
    if ExportPlot:
        SegFullFname = os.path.split(SegFpath)[1]
        
        SegFname = os.path.splitext(SegFullFname)[0]
        
        ExportFname = 'Plot_' + SegFname
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
        
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
        
    return






def Plot_Src_ResSrc_Trg_Images(SrcIm, ResSrcIm, TrgIm, SrcImLabel, TrgImLabel,
                               ImageType, ExportPlot, ExportDir, 
                               LogToConsole=False, dpi=80):
    
    """ 
    SrcIm/ResSrcIm/TrgIm could either be a SimpleITK image,
    i.e. type(SrcIm).__name__ == 'Image', or a list of SimpleITK images,
    i.e. type(SrcIm).__name__ == 'list'.
    
    e.g. SrcIm can be a list of sitk images of the labelmaps for each of 
    multiple ROIs.
    """
    
    # If ___Im is an Image (i.e. not a list), redefine it as a list with one 
    # item:  
    if type(SrcIm).__name__ == 'Image':
        SrcIm = [SrcIm]
        
    if type(ResSrcIm).__name__ == 'Image':
        ResSrcIm = [ResSrcIm]
        
    if type(TrgIm).__name__ == 'Image':
        TrgIm = [TrgIm]
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSize()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetSize()} \n',
                      f'\n   Target:           {TrgIm[0].GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSpacing()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetSpacing()} \n',
                      f'\n   Target:           {TrgIm[0].GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm[0].GetOrigin()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetOrigin()} \n',
                      f'\n   Target:           {TrgIm[0].GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm[0].GetDirection()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetDirection()} \n',
                      f'\n   Target:           {TrgIm[0].GetDirection()}')

    
    # Convert from SimpleITK image to Numpy arrays, assigning different values
    # for each ROI (i.e. 1, 2, 3, ...):
    SrcNpas = []
    ResSrcNpas = []
    TrgNpas = []
    
    # Loop through all ROIs:
    for r in range(len(SrcIm)):
        SrcNpas.append(sitk.GetArrayFromImage(SrcIm[r]))
        #SrcNpas.append((r + 1)*sitk.GetArrayFromImage(SrcIm[r]))
        
    for r in range(len(ResSrcIm)):
        ResSrcNpas.append(sitk.GetArrayFromImage(ResSrcIm[r]))
        #ResSrcNpas.append((r + 1)*sitk.GetArrayFromImage(ResSrcIm[r]))
        
    for r in range(len(TrgIm)):
        TrgNpas.append(sitk.GetArrayFromImage(TrgIm[r]))
        #TrgNpas.append((r + 1)*sitk.GetArrayFromImage(TrgIm[r]))
            
    
    # Find the maxima for all labelmaps:
    if 'Labelmaps' in ImageType:
        maxima = []
        
        for r in range(len(SrcNpas)):
            maxima.append(np.amax(SrcNpas[r]))
            
        for r in range(len(ResSrcNpas)):
            maxima.append(np.amax(ResSrcNpas[r]))
            
        for r in range(len(TrgNpas)):
            maxima.append(np.amax(TrgNpas[r]))
        
        MaxVal = max(maxima)

    
    # Prepare the figure:
    #Nslices = max(SrcIm.GetSize()[2], TrgIm.GetSize()[2])
    SrcNslices = SrcNpas[0].shape[0]
    ResSrcNslices = ResSrcNpas[0].shape[0]
    TrgNslices = TrgNpas[0].shape[0]
    Nslices = max(SrcNslices, TrgNslices)
    
    # Set the number of subplot rows and columns:
    Ncols = 3
    
    Nrows = Nslices
    
    
    # Create a figure with two subplots and the specified size:
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Define the colourmaps to use for each ROI:
    #cmaps = [plt.cm.spring, plt.cm.cool, plt.cm.autumn, plt.cm.winter, 
    #         plt.cm.hsv, plt.cm.jet]
    
    # Loop through each slice: 
    for s in range(Nslices):
         
        # Plot the Original Source slice:
        if s < SrcNslices:
            ax = plt.subplot(Nrows, Ncols, n)
            
            # Loop through all ROIs:
            #maxima = []
            for r in range(len(SrcNpas)):
                im = ax.imshow(SrcNpas[r][s], cmap=plt.cm.Greys_r)
                #im = ax.imshow(SrcNpas[r][s], cmap=cmaps[r])
                #maxima.append(np.amax(SrcNpas[r][s]))
                
            if 'Labelmaps' in ImageType:
                cbar = fig.colorbar(im, ax=ax)
                cbar.mappable.set_clim(0, MaxVal)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Original Source slice {s}')
            
        n += 1 # increment sub-plot number
            
        # Plot the Resampled Source slice:
        if s < ResSrcNslices:
            ax = plt.subplot(Nrows, Ncols, n)
            
            #maxima = []
            for r in range(len(ResSrcNpas)):
                im = ax.imshow(ResSrcNpas[r][s], cmap=plt.cm.Greys_r)
                #im = ax.imshow(ResSrcNpas[r][s], cmap=cmaps[r])
                #maxima.append(np.amax(ResSrcNpas[r][s]))
                
            if 'Labelmaps' in ImageType:
                cbar = fig.colorbar(im, ax=ax)
                cbar.mappable.set_clim(0, MaxVal)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Resampled Source slice {s}')
            
        n += 1 # increment sub-plot number
        
        # Plot the Target slice:
        if s < TrgNslices:
            ax = plt.subplot(Nrows, Ncols, n)
            
            #maxima = []
            for r in range(len(TrgNpas)):
                im = ax.imshow(TrgNpas[r][s], cmap=plt.cm.Greys_r)
                #im = ax.imshow(TrgNpas[r][s], cmap=cmaps[r])
                #maxima.append(np.amax(TrgNpas[r][s]))
                
            if 'Labelmaps' in ImageType:
                cbar = fig.colorbar(im, ax=ax)
                cbar.mappable.set_clim(0, MaxVal)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Target slice {s}')
            
        n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
                      
        ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
                      + SrcImLabel + '_' + TrgImLabel + '_' + ImageType \
                      + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return





def Plot_Src_ResSrc_Trg_NewTrg_Images(SrcIm, ResSrcIm, TrgIm, 
                                      SrcLMIms, ResSrcLMIms, 
                                      TrgLMIms, NewTrgLMIms,
                                      SrcImLabel, TrgImLabel,
                                      ImageType, ExportPlot, ExportDir, 
                                      LogToConsole=False, dpi=80):
    
    """ 
    ___Im is a 3D SimpleITK image from DICOMs.
    
    ___LMIms is either be a 3D labelmap as a SimpleITK image, or a list
    of 3D labelmaps as a SimpleITK image for each segment for ___Im.
    
    i.e. type(___LMIms).__name__ == 'Image', or
    i.e. type(___LMIms).__name__ == 'list'.
    
    """
    
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm.GetSize()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetSize()} \n',
                      f'\n   Original Target:  {TrgIm.GetSize()}')#,
                      #f'\n   New Target:       {NewTrgIm[0].GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm.GetSpacing()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetSpacing()} \n',
                      f'\n   Original Target:  {TrgIm.GetSpacing()}')#,
                      #f'\n   New Target:       {NewTrgIm[0].GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm.GetOrigin()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetOrigin()} \n',
                      f'\n   Original Target:  {TrgIm.GetOrigin()}')#,
                      #f'\n   New Target:       {NewTrgIm[0].GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm.GetDirection()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetDirection()} \n',
                      f'\n   Original Target:  {TrgIm.GetDirection()}')#,
                      #f'\n   New Target:       {NewTrgIm[0].GetDirection()}')
        
    
    if 'Images' in ImageType:
        # Convert from sitk to numpy array:
        SrcImNpa = sitk.GetArrayFromImage(SrcIm)
        ResSrcImNpa = sitk.GetArrayFromImage(ResSrcIm)
        TrgImNpa = sitk.GetArrayFromImage(TrgIm)
        
        
        
    
    if 'Labelmaps' in ImageType:
        # If ___Im is an Image (i.e. not a list), redefine it as a list with  
        # one item:  
        if SrcLMIms:
            if type(SrcLMIms).__name__ == 'Image':
                SrcLMIms = [SrcLMIms]
            
            if LogToConsole:
                print(f'\nScrIm has {len(SrcLMIms)} segments')
                
            # Add the segments contained in the list of sitk images into a 
            # single sitk image:
            SrcLMImsSum = InitialiseImage(SrcLMIms[0])
            
            for r in range(len(SrcLMIms)):
                SrcLMImsSum = AddImages(SrcLMImsSum, SrcLMIms[r])
        
            # Convert to numpy arrays:
            SrcLMImsSumNpa = sitk.GetArrayFromImage(SrcLMImsSum)
            
            SrcNslices = SrcLMImsSumNpa.shape[0]
    
    
            
        if ResSrcLMIms:
            if type(ResSrcLMIms).__name__ == 'Image':
                ResSrcLMIms = [ResSrcLMIms]
            
            if LogToConsole:
                print(f'ResScrIm has {len(ResSrcLMIms)} segments')
                
            # Add the segments contained in the list of sitk images into a 
            # single sitk image:
            ResSrcLMImsSum = InitialiseImage(ResSrcLMIms[0])
            
            for r in range(len(ResSrcLMIms)):
                ResSrcLMImsSum = AddImages(ResSrcLMImsSum, ResSrcLMIms[r])
                
            # Convert to numpy arrays:
            ResSrcLMImsSumNpa = sitk.GetArrayFromImage(ResSrcLMImsSum)
            
            ResSrcNslices = ResSrcLMImsSumNpa.shape[0]
        
        
            
        if TrgLMIms:
            if type(TrgLMIms).__name__ == 'Image':
                TrgLMIms = [TrgLMIms]
            
            if LogToConsole:
                print(f'TrgIm has {len(TrgLMIms)} segments')
                
            # Add the segments contained in the list of sitk images into a 
            # single sitk image:
            TrgLMImsSum = InitialiseImage(TrgLMIms[0])
            
            for r in range(len(TrgLMIms)):
                TrgLMImsSum = AddImages(TrgLMImsSum, TrgLMIms[r])
                
            # Convert to numpy arrays:
            TrgLMImsSumNpa = sitk.GetArrayFromImage(TrgLMImsSum)
            
            TrgNslices = TrgLMImsSumNpa.shape[0]
        
        
            
        if NewTrgLMIms:
            if type(NewTrgLMIms).__name__ == 'Image':
                NewTrgLMIms = [NewTrgLMIms]
        
            if LogToConsole:
                print(f'NewTrgIm has {len(NewTrgLMIms)} segments')
                
            # Add the segments contained in the list of sitk images into a 
            # single sitk image:
            NewTrgLMImsSum = InitialiseImage(NewTrgLMIms[0])
            
            for r in range(len(NewTrgLMIms)):
                NewTrgLMImsSum = AddImages(NewTrgLMImsSum, NewTrgLMIms[r])
                
            # Convert to numpy arrays:
            NewTrgLMImsSumNpa = sitk.GetArrayFromImage(NewTrgLMImsSum)
            
            NewTrgNslices = NewTrgLMImsSumNpa.shape[0]
            
            
    else:
        SrcNslices = SrcIm.GetSize()[2]
        ResSrcNslices = ResSrcIm.GetSize()[2]
        TrgNslices = TrgIm.GetSize()[2]
        NewTrgNslices = TrgNslices
    
    
    
    Nslices = max(SrcNslices, TrgNslices)

    
    # Get the maxima for all slices for all labelmaps:
    if 'Labelmaps' in ImageType:
        SrcMaximaBySlice = np.amax(SrcLMImsSumNpa, axis=(1,2))
        
        ResSrcMaximaBySlice = np.amax(ResSrcLMImsSumNpa, axis=(1,2))
            
        TrgMaximaBySlice = np.amax(TrgLMImsSumNpa, axis=(1,2))
        
        NewTrgMaximaBySlice = np.amax(NewTrgLMImsSumNpa, axis=(1,2))
           
        if LogToConsole:
            print(f'\nSrcMaximaBySlice   = {SrcMaximaBySlice}')
            print(f'\nResSrcMaximaBySlice = {ResSrcMaximaBySlice}')
            print(f'\nTrgMaximaBySlice    = {TrgMaximaBySlice}')
            print(f'\nNewTrgMaximaBySlice = {NewTrgMaximaBySlice}')
        
        
        #AllMaxima = [max(SrcMaximaBySlice), max(ResSrcMaximaBySlice),
        #             max(TrgMaximaBySlice), max(NewTrgMaximaBySlice)]
            
        
        # Determine which slices have at least one non-zero labelmap across all
        # labelmaps:
        SumOfMaximaBySlice = []
        
        for s in range(Nslices):
            count = 0
            
            if s < SrcNslices:
                count += SrcMaximaBySlice[s]
                
            if s < ResSrcNslices:
                count += ResSrcMaximaBySlice[s]
                
            if s < TrgNslices:
                count += TrgMaximaBySlice[s]
                
            if s < NewTrgNslices:
                count += NewTrgMaximaBySlice[s]
                
            SumOfMaximaBySlice.append(count)
        
        
        #print(f'\nSumOfMaximaBySlice = {SumOfMaximaBySlice}')
        
    
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = 4
    
    if 'Labelmaps' in ImageType:
        Nrows = np.count_nonzero(np.array(SumOfMaximaBySlice))
    else:
        Nrows = Nslices
    
    
    # Create a figure with two subplots and the specified size:
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    alpha = 0.2
    
    # Define the colormap to use for labelmaps overlaid over DICOMs:
    #cmap = plt.cm.Reds
    cmap = plt.cm.cool
    cmap = plt.cm.hsv
    
    
    # Loop through each slice: 
    for s in range(Nslices):
        # Only proceed if labelmaps are not to be displayed or (if they are) at 
        # least one of the labelmaps is non-zero for this slice number:
        if not 'Labelmaps' in ImageType or SumOfMaximaBySlice[s]:
            
            if s < SrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                if 'Images' in ImageType:
                    # Plot the Original Source DICOM image:
                    im = ax.imshow(SrcImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        # Plot the labelmap with tranparency:
                        im = ax.imshow(SrcLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        #m = 1 if SrcMaximaBySlice[s] < 1 else SrcMaximaBySlice[s]
                        #cbar = fig.colorbar(im, ax=ax)
                        #cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    # Plot the labelmap only:
                    im = ax.imshow(SrcLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if SrcMaximaBySlice[s] < 1 else SrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                        
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Original Source slice {s}')
                 
                
            n += 1 # increment sub-plot number
            
            if s < ResSrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                if 'Images' in ImageType:
                    # Plot the Resampled Source DICOM image:
                    im = ax.imshow(ResSrcImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(ResSrcLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        #m = 1 if ResSrcMaximaBySlice[s] < 1 else ResSrcMaximaBySlice[s]
                        #cbar = fig.colorbar(im, ax=ax)
                        #cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(ResSrcLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if ResSrcMaximaBySlice[s] < 1 else ResSrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                        
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Resampled Source slice {s}')
                
                
            n += 1 # increment sub-plot number
                
            if s < TrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
               
                if 'Images' in ImageType:
                    # Plot the Target slice:
                    im = ax.imshow(TrgImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(TrgLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        #m = 1 if TrgMaximaBySlice[s] < 1 else TrgMaximaBySlice[s]
                        #cbar = fig.colorbar(im, ax=ax)
                        #cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(TrgLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if TrgMaximaBySlice[s] < 1 else TrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Target slice {s}')
                
            n += 1 # increment sub-plot number
                
            
            if s < NewTrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                    
                if 'Images' in ImageType:      
                    # Plot the New Target slice:
                    im = ax.imshow(TrgImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(NewTrgLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        #m = 1 if NewTrgMaximaBySlice[s] < 1 else NewTrgMaximaBySlice[s]
                        #cbar = fig.colorbar(im, ax=ax)
                        #cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(NewTrgLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if NewTrgMaximaBySlice[s] < 1 else NewTrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'New Target slice {s}')
                
            n += 1 # increment sub-plot number
        
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
                      
        ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
                      + SrcImLabel + '_' + TrgImLabel + '_New' + TrgImLabel \
                      + '_' + ImageType + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return