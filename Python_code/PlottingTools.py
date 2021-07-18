# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:52:36 2020

@author: ctorti
"""



# Import packages and functions:
#import os
#import time
#from pydicom import dcmread
#import SimpleITK as sitk
#import numpy as np
#import matplotlib.pyplot as plt

##import importlib
##import GeneralTools
##importlib.reload(GeneralTools)
#from GeneralTools import Unpack

#from ConversionTools import GetIndsInRois
#from ConversionTools import ConvertContourDataToIcsPoints

#from DicomTools import ImportDicoms
#from DicomTools import GetDicomSOPuids
#from DicomTools import GetRoiLabels

#from ImageTools import InitialiseImage
#from ImageTools import AddImages
#from ImageTools import SumAllLabmapIms

#from RtsTools import GetCIStoSliceInds
#from RtsTools import GetCStoSliceInds

#from SegTools import GetPFFGStoSliceInds



"""
******************************************************************************
******************************************************************************
PLOTTING FUNCTIONS
******************************************************************************
******************************************************************************
"""       





def PlotImagesAndContours(Image, ImageType, DicomDir, RtsFpath,
                          ExportPlot, ExportDir, Label='', dpi=80,
                          LogToConsole=False):
    """
    Image can either be a 3D DICOM image, or 3D labelmap, or a list of 3D 
    labelmaps.
    """
    
    import SimpleITK as sitk
    from pydicom import dcmread
    import matplotlib.pyplot as plt
    from ConversionTools import GetIndsInRois
    from DicomTools import GetDicomSOPuids
    from RtsTools import GetCStoSliceInds
    from DicomTools import GetRoiLabels
    from ConversionTools import ConvertContourDataToIndices
    from GeneralTools import Unpack
    import time
    import os
    
    # Convert Image to a list of Numpy arrays (even if only one 3D 
    # DICOM/labelmap):
    Ndas = []
    
    if isinstance(Image, list):
        for r in range(len(Image)):
            Ndas.append(sitk.GetArrayFromImage(Image[r]))
            
    elif isinstance(Image, sitk.SimpleITK.Image):
        Ndas.append(sitk.GetArrayFromImage(Image))
        
    else:
        msg = "Image is type " + str(type(Image)) + ". It must either be a " \
              + "SimpleITK Image or a list of SimpleITK Images"
                        
        raise Exception(msg)
    
    # Import the RTS object:
    Rts = dcmread(RtsFpath)
    
    
    IndsInRois = GetIndsInRois(RtsFpath, DicomDir)
    
    Nrois = len(IndsInRois)
    
    
    ContourToSliceNumsByRoi = GetCStoSliceInds(Rts, GetDicomSOPuids(DicomDir))
    
    print(f'ContourToSliceNumsByRoi = {ContourToSliceNumsByRoi}')
    
    # The slices to plot are all slice numbers in ContourToSliceNums and all
    # non-zero masks in Image if ImageType is 'Labelmap':
    #SliceNums = []
    #SliceNums.extend(ContourToSliceNums)
    
    # Get all the slice numbers in ContourToSliceNumsByRoi:
    SliceNums = []
    
    for nums in ContourToSliceNumsByRoi:
        SliceNums.extend(nums)
        
    # If ImageType is 'Labelmap' append the slice numbers of any non-zero 
    # frames:
    
    if ImageType in ['Labelmap', 'labelmap', 'LabelMap']:
        Non0inds = []
        
        for r in range(len(Ndas)):
            # Get the indices of non-zero arrays in Ndas[r]: 
            for i in range(Ndas[r].shape[0]):
                if Ndas[r][i].sum() > 0:
                    Non0inds.append(i)
            
            #Nnon0s = len(Non0inds)
            
            SliceNums.extend(Non0inds)
        
        print(f'Non0inds = {Non0inds}')
        
    
    # Get the list of unique slice numbers:
    SliceNums = list(set(SliceNums))
    
    print(f'SliceNums = {SliceNums}')
    
    #Nslices = len(SliceNums)
    
    # Set the number of subplot rows and columns:
    #Ncols = np.int8(np.ceil(np.sqrt(Nnon0s)))
    Ncols = Nrois
    
    # Limit the number of rows to 5 otherwise it's difficult to make sense of
    # the images:
    #if Ncols > 5:
    #    Ncols = 5
    
    #Nrows = np.int8(np.ceil(Nnon0s/Ncols))
    Nrows = len(SliceNums)
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows), dpi=300)
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=300)
    else:
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows))
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15*Ncols, 6*Nrows), dpi=dpi)
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=dpi)
        
    
    n = 0 # sub-plot number
    
    # List of colours for multiple ROIs per slice (limited to 5 for now):
    colours = ['r', 'c', 'm', 'y', 'g']
    
    # Get ROI labels:
    RoiLabels = GetRoiLabels(Rts)
    
    # Loop through each slice number in SliceNums: 
    for s in SliceNums:
        
        if LogToConsole:
                print('\ns =', s)
                
        
        # Only plot if Nda[s] is not empty (i.e. if it is a DICOM 2D image or
        # a non-empty mask):
        #if Nda[s].sum() > 0:
        
        #n += 1 # increment sub-plot number
        
        ##plt.subplot(Nrows, Ncols, n)
        ##plt.axis('off')
        #ax = plt.subplot(Nrows, Ncols, n)

        #ax = plt.subplot(Nrows, Ncols, n, aspect=AR)
          
        # Plot the image:
        #ax.imshow(Nda[s], cmap=plt.cm.Greys_r)
        ##ax.set_aspect(ar)
        
        
        # Loop through all ROIs:
        for r in range(Nrois):
            n += 1 # increment sub-plot number
            
            ax = plt.subplot(Nrows, Ncols, n)
            
            # Plot the image:
            if len(Ndas) > 1:
                ax.imshow(Ndas[r][s], cmap=plt.cm.Greys_r)
                
            else:
                ax.imshow(Ndas[0][s], cmap=plt.cm.Greys_r)
            
            
            ContourToSliceNums = ContourToSliceNumsByRoi[r]
            
            # Does this slice have a contour?
            if s in ContourToSliceNums:
                #ind = ContourToSliceNums.index(s)
                
                # Find all occurances of slice s in ContourToSliceNums:
                inds = [i for i, val in enumerate(ContourToSliceNums) if val==s]
                
                for i in range(len(inds)):
                    if ImageType in ['DICOM', 'Dicom', 'dicom']:
                        Points = IndsInRois[r][inds[i]]
                        
                        # Convert contour data from Patient Coordinate System to  
                        # Image Coordinate System:
                        Indices = ConvertContourDataToIndices(Points, DicomDir)
                        
                    else:
                        # ImageType in ['Labelmap', 'labelmap']
                        Indices = IndsInRois[r][inds[i]]
                    
                    
                    # Unpack the indices:
                    X, Y, Z = Unpack(Indices)
                    
                    # Plot the contour points:
                    plt.plot(X, Y, linewidth=0.5, c=colours[i]);
                    
                ax.set_title(f'Slice {s} \n' + RoiLabels[r])
                
            else:
                ax.set_title(f'Slice {s}')
                
            
            
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            
            #ax.set_title(f'Slice {s}')
            #plt.axis('off')
        
            
            
    if ExportPlot:
        CDT = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        ExportFname = CDT + '_' + ImageType + '_' + Label
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
        
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
        
    return





def PlotDicomsAndSegments(DicomDir, SegFpath, ExportPlot, ExportDir, 
                          LogToConsole):
    
    from pydicom import dcmread
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    from DicomTools import ImportDicoms
    from DicomTools import GetDicomSOPuids
    from SegTools import GetPFFGStoSliceInds
    
    # Import the DICOMs:
    Dicoms = ImportDicoms(DicomDir=DicomDir)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir=DicomDir)
    
    # Import the SEG ROI:
    SegRoi = dcmread(SegFpath)
    
    # Get the Per-frameFunctionalGroupsSequence-to-slice indices:
    PFFGStoSliceInds = GetPFFGStoSliceInds(SegRoi, SOPuids)
    
    # The 3D labelmap:
    SegNpa = SegRoi.pixel_array
    
    if LogToConsole:
        print(f'Segments exist on slices {PFFGStoSliceInds}')
        print(f'Shape of PixelData = {SegNpa.shape}\n')

    
    
    # Prepare the figure.
    
    # All DICOM slices that have segmentations will be plotted. The number of
    # slices to plot (= number of segments):
    #Nslices = int(SegRoi.NumberOfFrames)
    Nslices = len(PFFGStoSliceInds)
    
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
        s = PFFGStoSliceInds[i]
        
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






def Plot_Src_ResSrc_Trg_Images(SrcIm, ResSrcIm, TrgIm, SrcLabel, TrgLabel, 
                               Case=None, ImageType='', ExportPlot=False, 
                               ExportDir='', LogToConsole=False, dpi=80):
    
    """ 
    SrcIm/ResSrcIm/TrgIm could either be a SimpleITK image,
    i.e. type(SrcIm).__name__ == 'Image', or a list of SimpleITK images,
    i.e. type(SrcIm).__name__ == 'list'.
    
    e.g. SrcIm can be a list of sitk images of the labelmaps for each of 
    multiple ROIs.
    """
    
    import SimpleITK as sitk
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import time
    
    # If ___Im is an Image (i.e. not a list), redefine it as a list with one 
    # item:  
    if type(SrcIm).__name__ == 'Image':
        SrcIm = [SrcIm]
        
    if type(ResSrcIm).__name__ == 'Image':
        ResSrcIm = [ResSrcIm]
        
    if type(TrgIm).__name__ == 'Image':
        TrgIm = [TrgIm]
        
        
    if Case == '5':
        ResOrReg = 'Registered'
        ResOrReg_short = '_Reg'
        #print('\nCase = 5 -->', ResOrReg)
    else:
        ResOrReg = 'Resampled'
        ResOrReg_short = '_Resamp'
        
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSize()} \n',
                      '\n  ', ResOrReg, f'Source: {ResSrcIm[0].GetSize()} \n',
                      f'\n   Target:           {TrgIm[0].GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSpacing()} \n',
                      '\n  ', ResOrReg, f'Source: {ResSrcIm[0].GetSpacing()} \n',
                      f'\n   Target:           {TrgIm[0].GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm[0].GetOrigin()} \n',
                      '\n  ', ResOrReg, f'Source: {ResSrcIm[0].GetOrigin()} \n',
                      f'\n   Target:           {TrgIm[0].GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm[0].GetDirection()} \n',
                      '\n  ', ResOrReg, f'Source: {ResSrcIm[0].GetDirection()} \n',
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
            ax.set_title(f'Original Source ({SrcLabel}) slice {s}')
            
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
            ax.set_title(ResOrReg + f' Source ({SrcLabel}) slice {s}')
            
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
            ax.set_title(f'Target ({TrgLabel}) slice {s}')
            
        n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
        
        if Case:
            ExportFname = CurrentDateTime + '_Case' + Case + '_Orig' \
                          + SrcLabel + ResOrReg_short + SrcLabel + '_' \
                          + TrgLabel + '_' + ImageType + '.jpg'
        else:
            ExportFname = CurrentDateTime + '_Orig' \
                          + SrcLabel + ResOrReg_short + SrcLabel + '_' \
                          + TrgLabel + '_' + ImageType + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return





def Plot_Src_ResSrc_Trg_NewTrg_Images(SrcIm, ResSrcIm, TrgIm, 
                                      SrcLMIms, ResSrcLMIms, 
                                      TrgLMIms, NewTrgLMIms,
                                      SrcLabel, TrgLabel, Case,
                                      ImageType, ExportPlot, ExportDir, 
                                      LogToConsole=False, dpi=80):
    """ 
    Note:
    
    ResSrcIm = Resampled Source image could instead be a Registered Source
    image.
    
    
    ___Im is a 3D SimpleITK image from DICOMs.
    
    ___LMIms is either be a 3D labelmap as a SimpleITK image, or a list
    of 3D labelmaps as a SimpleITK image for each segment for ___Im.
    
    i.e. type(___LMIms).__name__ == 'Image', or
    i.e. type(___LMIms).__name__ == 'list'.
    
    """
    
    if Case == '5':
        ResOrReg = 'Registered'
        ResOrReg_short = '_Reg'
        #print('\nCase = 5 -->', ResOrReg)
    else:
        ResOrReg = 'Resampled'
        ResOrReg_short = '_Resamp'
        #print('\nCase != 5 -->', ResOrReg)
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'   Original Source:  {SrcIm.GetSize()} \n',
                      '  ', ResOrReg, f'Source: {ResSrcIm.GetSize()} \n',
                      f'   Original Target:  {TrgIm.GetSize()}')#,
                      #f'   New Target:       {NewTrgIm[0].GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'   Original Source:  {SrcIm.GetSpacing()} \n',
                      '  ', ResOrReg, f'Source: {ResSrcIm.GetSpacing()} \n',
                      f'   Original Target:  {TrgIm.GetSpacing()}')#,
                      #f'   New Target:       {NewTrgIm[0].GetSpacing()}')
        
        """
        print('\nThe Source and Target Origin:\n',
                      f'   Original Source:  {SrcIm.GetOrigin()} \n',
                      f'   Resampled Source: {ResSrcIm.GetOrigin()} \n',
                      f'   Original Target:  {TrgIm.GetOrigin()}')#,
                      #f'   New Target:       {NewTrgIm[0].GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'   Original Source:  {SrcIm.GetDirection()} \n',
                      f'   Resampled Source: {ResSrcIm.GetDirection()} \n',
                      f'   Original Target:  {TrgIm.GetDirection()}')#,
                      #f'   New Target:       {NewTrgIm[0].GetDirection()}')
        """
    
    
    import SimpleITK as sitk
    import numpy as np
    import matplotlib.pyplot as plt
    from ImageTools import SumAllLabmapIms
    import os
    import time
    
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
        
            # Add all segments pixel-by-pixel:
            SrcLMImsSumNpa = SumAllLabmapIms(SrcLMIms)
            
            SrcNslices = SrcLMImsSumNpa.shape[0]
    
    
            
        if ResSrcLMIms:
            if type(ResSrcLMIms).__name__ == 'Image':
                ResSrcLMIms = [ResSrcLMIms]
            
            if LogToConsole:
                print(f'ResScrIm has {len(ResSrcLMIms)} segments')
                
            # Add all segments pixel-by-pixel:
            ResSrcLMImsSumNpa = SumAllLabmapIms(ResSrcLMIms)
            
            ResSrcNslices = ResSrcLMImsSumNpa.shape[0]
        
        
            
        if TrgLMIms:
            if type(TrgLMIms).__name__ == 'Image':
                TrgLMIms = [TrgLMIms]
            
            if LogToConsole:
                print(f'TrgIm has {len(TrgLMIms)} segments')
                
            # Add all segments pixel-by-pixel:
            TrgLMImsSumNpa = SumAllLabmapIms(TrgLMIms)
            
            TrgNslices = TrgLMImsSumNpa.shape[0]
        
        
            
        if NewTrgLMIms:
            if type(NewTrgLMIms).__name__ == 'Image':
                NewTrgLMIms = [NewTrgLMIms]
        
            if LogToConsole:
                print(f'NewTrgIm has {len(NewTrgLMIms)} segments')
                
            # Add all segments pixel-by-pixel:
            NewTrgLMImsSumNpa = SumAllLabmapIms(NewTrgLMIms)
            
            NewTrgNslices = NewTrgLMImsSumNpa.shape[0]
            
            # Convert only the last segment (last labelmap image):
            NewTrgOnlyNewLM = sitk.GetArrayFromImage(NewTrgLMIms[-1])
            
            
        print('\nThe Source and Target labelmap image size:\n',
              f'   Original Source:  {SrcLMIms[0].GetSize()} \n',
              ResOrReg, f'Source: {ResSrcLMIms[0].GetSize()} \n',
              f'   Original Target:  {TrgLMIms[0].GetSize()}')
            
            
    else:
        SrcNslices = SrcIm.GetSize()[2]
        ResSrcNslices = ResSrcIm.GetSize()[2]
        TrgNslices = TrgIm.GetSize()[2]
        NewTrgNslices = TrgNslices
    
    
    
    Nslices = max(SrcNslices, TrgNslices)

    
    # Get the maxima for all slices for all labelmaps along x (axis=2) and
    # y (axis=1):
    if 'Labelmaps' in ImageType:
        SrcMaximaBySlice = np.amax(SrcLMImsSumNpa, axis=(1,2))
        
        ResSrcMaximaBySlice = np.amax(ResSrcLMImsSumNpa, axis=(1,2))
            
        TrgMaximaBySlice = np.amax(TrgLMImsSumNpa, axis=(1,2))
        
        NewTrgMaximaBySlice = np.amax(NewTrgLMImsSumNpa, axis=(1,2))
        
        NewTrgOnlyNewLmMaximaBySlice = np.amax(NewTrgOnlyNewLM, axis=(1,2))
           
        if LogToConsole:
            print(f'\nSrcMaximaBySlice   = {SrcMaximaBySlice}')
            print(f'ResSrcMaximaBySlice = {ResSrcMaximaBySlice}')
            print(f'TrgMaximaBySlice    = {TrgMaximaBySlice}')
            print(f'NewTrgMaximaBySlice = {NewTrgMaximaBySlice}')
            print(f'NewTrgOnlyNewLmMaximaBySlice = {NewTrgOnlyNewLmMaximaBySlice}')
        
        
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
    #Ncols = 4
    Ncols = 5 # 19/11: Add New Target new segment only in 5th column
    
    if 'Labelmaps' in ImageType:
        NonZeroSlices = np.nonzero(SumOfMaximaBySlice)
        SliceNums = NonZeroSlices[0]
        
        #print(f'\nNonZeroSlices = {NonZeroSlices}')
        
        #Nrows = np.count_nonzero(np.array(SumOfMaximaBySlice))
        Nrows = len(SliceNums)
    else:
        SliceNums = list(range(Nslices))
        
        Nrows = Nslices
        
    
    #if ExportPlot:
    #    dpi = 120
    
    
    # Create a figure with two subplots and the specified size:
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    alpha = 0.2
    alpha = 0.5
    
    # Define the colormap to use for labelmaps overlaid over DICOMs:
    #cmap = plt.cm.Reds
    #cmap = plt.cm.cool
    #cmap = plt.cm.hsv
    cmap = plt.cm.nipy_spectral
    
    
    # Loop through each SliceNums: 
    for s in SliceNums:
        #print(f's = {s}')
        
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
                        
                        m = 1 if SrcMaximaBySlice[s] < 1 else SrcMaximaBySlice[s]
                        cbar = fig.colorbar(im, ax=ax)
                        cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    # Plot the labelmap only:
                    im = ax.imshow(SrcLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if SrcMaximaBySlice[s] < 1 else SrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                        
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                if len(SrcLMIms) == 1:
                    ax.set_title(f'Original Source ({SrcLabel}) slice {s}'\
                                 + '\n(only segmentation to be copied)')
                else:
                    ax.set_title(f'Original Source ({SrcLabel}) slice {s}'\
                                 + '\n(sum of all segments)')
                 
                
            n += 1 # increment sub-plot number
            
            if s < ResSrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                if 'Images' in ImageType:
                    # Plot the Resampled Source DICOM image:
                    im = ax.imshow(ResSrcImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(ResSrcLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        m = 1 if ResSrcMaximaBySlice[s] < 1 else ResSrcMaximaBySlice[s]
                        cbar = fig.colorbar(im, ax=ax)
                        cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(ResSrcLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if ResSrcMaximaBySlice[s] < 1 else ResSrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                        
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                if len(ResSrcLMIms) == 1:
                    ax.set_title(ResOrReg + f' Source ({SrcLabel}) slice {s}'\
                                 + '\n(only segmentation to be copied)')
                else:
                    ax.set_title(ResOrReg + f' Source ({SrcLabel}) slice {s}'\
                                 + '\n(sum of all segments)')
                
                
            n += 1 # increment sub-plot number
                
            if s < TrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
               
                if 'Images' in ImageType:
                    # Plot the Target slice:
                    im = ax.imshow(TrgImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(TrgLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        m = 1 if TrgMaximaBySlice[s] < 1 else TrgMaximaBySlice[s]
                        cbar = fig.colorbar(im, ax=ax)
                        cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(TrgLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if TrgMaximaBySlice[s] < 1 else TrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Original Target ({TrgLabel}) slice {s}'\
                             + '\n(sum of all segments)')
                
            n += 1 # increment sub-plot number
                
            
            if s < NewTrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                    
                if 'Images' in ImageType:      
                    # Plot the New Target slice:
                    im = ax.imshow(TrgImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(NewTrgLMImsSumNpa[s], cmap=cmap,
                                       alpha=alpha)
                        
                        m = 1 if NewTrgMaximaBySlice[s] < 1 else NewTrgMaximaBySlice[s]
                        cbar = fig.colorbar(im, ax=ax)
                        cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(NewTrgLMImsSumNpa[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if NewTrgMaximaBySlice[s] < 1 else NewTrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'New Target ({TrgLabel}) slice {s}'\
                             + '\n(sum of all segments)')
                
            n += 1 # increment sub-plot number
            
                
            
            if s < NewTrgNslices:
                # Plot the new segment of New Target only:
                ax = plt.subplot(Nrows, Ncols, n)
                    
                if 'Images' in ImageType:      
                    # Plot the New Target slice:
                    im = ax.imshow(TrgImNpa[s], cmap=plt.cm.Greys_r)
                    
                    if 'Labelmaps' in ImageType:
                        im = ax.imshow(NewTrgOnlyNewLM[s], cmap=cmap,
                                       alpha=alpha)
                        
                        m = 1 if NewTrgMaximaBySlice[s] < 1 else NewTrgMaximaBySlice[s]
                        cbar = fig.colorbar(im, ax=ax)
                        cbar.mappable.set_clim(0, m)
                    
                elif 'Labelmaps' in ImageType:
                    im = ax.imshow(NewTrgOnlyNewLM[s], cmap=plt.cm.Greys_r)
                    
                    m = 1 if NewTrgOnlyNewLmMaximaBySlice[s] < 1 else NewTrgOnlyNewLmMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'New Target ({TrgLabel}) slice {s}'\
                             + '\n(only new segment)')
                
            n += 1 # increment sub-plot number
        
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
                      
        ExportFname = CurrentDateTime + '_Case' + Case + '_Orig' + SrcLabel \
                      + ResOrReg_short + SrcLabel + '_' + TrgLabel + '_New' \
                      + TrgLabel + '_' + ImageType + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return







def PlotPixelArrays(ListOfPixArrs, ListOfFrameToSliceInds, ListOfPlotTitles, 
                    AddTxt=None, ExportPlot=False, ExportDir=None, dpi=80, 
                    LogToConsole=False):
    """ 
    Note:
    
    ListOfPixArrs is a list (which could be of length 1) of pixel arrays.
    
    ListOfFrameToSliceInds is a list (which could be of length 1) of 
    FrameToSliceInds.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import time
    from copy import deepcopy
    
    for i in range(len(ListOfFrameToSliceInds)):
        print(f'\nListOfPlotTitles[{i}] = {ListOfPlotTitles[i]}')
        
        print(f'\nListOfFrameToSliceInds[{i}] = {ListOfFrameToSliceInds[i]}')
        
        print(f'\nListOfPixArrs[{i}].shape = {ListOfPixArrs[i].shape}')
    
    
    if len(ListOfFrameToSliceInds) > 1:
        # Get the slice numbers that cover all FrameToSliceInds in 
        # ListOfFrameToSliceInds:
        AllSliceNums = []
        
        for FrameToSliceInds in ListOfFrameToSliceInds:
            for SliceInd in FrameToSliceInds:
                AllSliceNums.append(SliceInd)
        
        # Reduce to unique slice numbers:
        AllSliceNums = list(set(AllSliceNums))    
    else:
        AllSliceNums = deepcopy(ListOfFrameToSliceInds[0])
    
        
    # Get the maxima for all frames for all pixel arrays along x (axis=2) and
    # y (axis=1):
    ListOfMaximaByFrame = []
    
    for PixArr in ListOfPixArrs:
    
        ListOfMaximaByFrame.append(np.amax(PixArr, axis=(1,2)))
    
    
    if LogToConsole:
        print(f'\nListOfMaximaByFrame = {ListOfMaximaByFrame}')
    
    
    
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfPixArrs)
    Nrows = len(AllSliceNums)
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Loop through each AllSliceNums: 
    for SliceNum in AllSliceNums:
        #print(f's = {SliceNum}')
        
        # Loop through each pixel array:
        for p in range(len(ListOfPixArrs)):
            PixArr = ListOfPixArrs[p]
            
            Frame2SliceInds = ListOfFrameToSliceInds[p]
            
            PlotTitle = ListOfPlotTitles[p]
            
            #MaximaByFrame = ListOfMaximaByFrame[p]
            
            if SliceNum in Frame2SliceInds:
                FrameNum = Frame2SliceInds.index(SliceNum)
                
                Frame = PixArr[FrameNum]
                
                ax = plt.subplot(Nrows, Ncols, n)
                
                im = ax.imshow(Frame, cmap=plt.cm.Greys_r)
                
                #m = 1 if MaximaByFrame[FrameNum] < 1 else MaximaByFrame[FrameNum]
                #cbar = fig.colorbar(im, ax=ax)
                #cbar.mappable.set_clim(0, m)
                fig.colorbar(im, ax=ax)
                        
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                if AddTxt:
                    plotTitle = PlotTitle\
                                + f' Slice {SliceNum} Frame {FrameNum} \n'\
                                + AddTxt
                else:
                    plotTitle = PlotTitle\
                                + f' Slice {SliceNum} Frame {FrameNum}'
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        if AddTxt:
            AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
            AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
            AddTxt = AddTxt.replace('\n', '')
        else:
            AddTxt = ''
            
        ExportFname = CurrentDateTime + '_PixArrs_' + AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return









def PlotPixArrsFromListOfSegs_v1(ListOfSegs, ListOfDicomDirs, ListOfPlotTitles,
                                 PlotAllSlices=False, ExportPlot=False, 
                                 ExportDir=None, ExportFname=None, dpi=80, 
                                 LogToConsole=False):
    """ 
    Note:
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    from pathlib import Path
    #import time
    #from copy import deepcopy
    from pydicom import dcmread
    #from DicomTools import GetDicomFpaths
    from SegTools import GetSegDataFromListOfSegs
    from GeneralTools import UniqueItems#, Unpack, PrintTitle
    
    ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
    ListOfIPPs, ListOfDirections, ListOfSpacings\
    = GetSegDataFromListOfSegs(ListOfSegs, ListOfDicomDirs, LogToConsole)
    
    # Get the list of slice numbers to be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices)))
    else:
        AllSliceNums = UniqueItems(Items=ListOfF2SindsBySeg)
    
    
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfSegs)
    Nrows = len(AllSliceNums)
    
    Cmaps = ['Reds', 'Blues', 'Greens', 'Oranges', 'Purples']
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    DcmAlpha = 0.2
    DcmAlpha = 0.5
    #DcmAlpha = 1
    
    SegAlpha = 0.5
    #SegAlpha = 0.2
    
    n = 1 # initialised sub-plot number
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=dpi)
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
    
    # Loop through each AllSliceNums: 
    for SliceNum in AllSliceNums:
        if LogToConsole:
            print(f'\nSliceNum = {SliceNum}')
        
        # Loop through each data set:
        for i in range(len(ListOfPixArrBySeg)):
            PixArrBySeg = ListOfPixArrBySeg[i]
            F2SindsBySeg = ListOfF2SindsBySeg[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            #DicomDir = ListOfDicomDirs[i]
            
            if LogToConsole:
                    print(f'\n   i = {i}')
                    print(f'   PlotTitle = {PlotTitle}')
                    print(f'   F2SindsBySeg = {F2SindsBySeg}')
                    
            ax = plt.subplot(Nrows, Ncols, n)
            
            if SliceNum < len(DicomFpaths):
                DicomPixArr = dcmread(DicomFpaths[SliceNum]).pixel_array
            
                #im = ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r, alpha=DcmAlpha)
                
                
            # Loop through each segment:
            for s in range(len(F2SindsBySeg)):
                F2Sinds = F2SindsBySeg[s]
                PixArr = PixArrBySeg[s]
                
                """ There are only len(Cmaps) colormaps defined above. Wrap the 
                segment index s if there are more segments than the number of 
                defined colormaps. """
                if s < len(Cmaps):
                    Cmap = Cmaps[s]
                else:
                    m = s//len(Cmaps)
                    
                    Cmap = Cmaps[s - m*len(Cmaps)]
                    
                if SliceNum in F2Sinds:
                    FrameNums = [i for i, e in enumerate(F2Sinds) if e==SliceNum]
                    
                    # Loop through all frame numbers:
                    for f in range(len(FrameNums)):
                        FrameNum = FrameNums[f]
                    
                        Frame = PixArr[FrameNum]
                        
                        if LogToConsole:
                            print(f'      SliceNum = {SliceNum} is in F2Sinds')
                            print(f'      FrameNum = {FrameNum}')
                            #print(f'      Contour = {Contour}')
                            print(f'      Frame.shape = {Frame.shape}')
                    
                        #ax = plt.subplot(Nrows, Ncols, n)
                        #im = ax.imshow(Frame, cmap=plt.cm.nipy_spectral, 
                        #               alpha=alpha)
                        ax.imshow(Frame, cmap=Cmap, alpha=SegAlpha)
                        
                        FrameTxt = f'Frame {FrameNum}'
                        
                else:
                    FrameTxt = 'No frame'
                    
                    if LogToConsole:
                        print(f'      SliceNum = {SliceNum} is NOT in F2Sinds')
                    
                SliceTxt = f'Slice {SliceNum}'
                
                #if AddTxt:
                #    plotTitle = PlotTitle \
                #                + f'\n' + SliceTxt + '\n'\
                #                + f'\n' + FrameTxt + '\n'\
                #                + AddTxt
                #else:
                #    plotTitle = PlotTitle \
                #                + f'\n' + SliceTxt + '\n'\
                #                + f'\n' + FrameTxt
                plotTitle = f'{PlotTitle}\n\n{SliceTxt}\n{FrameTxt}'
                
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        #CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        #CurrentTime = time.strftime("%H%M%S", time.gmtime())
        #CurrentDateTime = CurrentDate + '_' + CurrentTime
        #
        #if AddTxt:
        #    AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
        #    AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
        #    AddTxt = AddTxt.replace('\n', '')
        #else:
        #    AddTxt = ''
        #    
        ##ExportFname = CurrentDateTime + '_Segs_' + AddTxt + '.jpg'
        #ExportFname = AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
        
    return





def PlotPixArrsFromListOfSegs_v2(ListOfSegs, ListOfDicomDirs, ListOfPlotTitles,
                                 SearchString, PlotAllSlices=False, AddTxt=None, 
                                 ExportPlot=False, ExportDir=None, dpi=80, 
                                 LogToConsole=False):
    """ 
    Note:
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG ROWS and different SLICES ALONG COLUMNS. 
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    from pathlib import Path
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from DicomTools import GetDicomFpaths
    
    ListOfPixArrs, ListOfF2Sinds, ListOfSegNums,\
    ListOfDcmFpaths = GetSegDataFromListOfSegs(ListOfSegs, ListOfDicomDirs, 
                                               SearchString, LogToConsole)
    
    # Get the list of slice numbers to be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices)))
            
        
    else:
        if len(ListOfF2Sinds) > 1:
            # Get the slice numbers that cover all FrameToSliceInds in 
            # ListOfFrameToSliceInds:
            AllSliceNums = []
            
            for F2Sinds in ListOfF2Sinds:
                if F2Sinds:
                    for SliceInd in F2Sinds:
                        AllSliceNums.append(SliceInd)
            
            # Reduce to unique slice numbers:
            AllSliceNums = list(set(AllSliceNums))    
        else:
            AllSliceNums = deepcopy(ListOfF2Sinds[0])
    
    
    if True:#LogToConsole:
        print(f'\nAllSliceNums = {AllSliceNums}')
        
    # Get the maxima for all frames for all pixel arrays along x (axis=2) and
    # y (axis=1):
    #ListOfMaximaByFrame = []
    #
    #for PixArr in ListOfPixArrs:
    #
    #    ListOfMaximaByFrame.append(np.amax(PixArr, axis=(1,2)))
    #
    #
    #if LogToConsole:
    #    print(f'\nListOfMaximaByFrame = {ListOfMaximaByFrame}')
    
    
    
    # Prepare the figure.
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    #alpha = 0.2
    alpha = 0.5
    
    # Set the number of subplot rows and columns:
    Nrows = len(ListOfPixArrs)
    Ncols = len(AllSliceNums)
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=dpi)
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Loop through each data set:
    for i in range(len(ListOfPixArrs)):
        PixArr = ListOfPixArrs[i]
        
        F2Sinds = ListOfF2Sinds[i]
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        PlotTitle = ListOfPlotTitles[i]
        
        #MaximaByFrame = ListOfMaximaByFrame[i]
        
        #print(f'\n\nDataset {i + 1}, len(DcmFpaths) = {len(DcmFpaths)}')
        
        # Loop through each AllSliceNums: 
        for SliceNum in AllSliceNums:
            #print(f's = {SliceNum}')
            
            ax = plt.subplot(Nrows, Ncols, n)
            
            if SliceNum < len(DcmFpaths):
                DcmPixArr = dcmread(DcmFpaths[SliceNum]).pixel_array
            
                #im = ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
                ax.imshow(DcmPixArr, cmap=plt.cm.Greys_r)#, alpha=alpha)
                
                if F2Sinds:
                    if SliceNum in F2Sinds:
                        FrameNum = F2Sinds.index(SliceNum)
                        
                        Frame = PixArr[FrameNum]
                        
                        #ax = plt.subplot(Nrows, Ncols, n)
                        
                        #im = ax.imshow(Frame, cmap=plt.cm.nipy_spectral, 
                        #               alpha=alpha)
                        ax.imshow(Frame, cmap=plt.cm.nipy_spectral, alpha=alpha)
                        
                        FrameTxt = f'Frame {FrameNum}'
                    else:
                        FrameTxt = 'No frame'
                else:
                    FrameTxt = 'No frame'
                    
                 
                SliceTxt = f'Slice {SliceNum + 1}'
                SliceTxt = ''
                
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                
                if AddTxt:
                    plotTitle = PlotTitle\
                                + '\n\n' + SliceTxt + '\n'\
                                + '\n' + FrameTxt + '\n'\
                                + AddTxt
                else:
                    plotTitle = PlotTitle\
                                + '\n\n' + SliceTxt + '\n'\
                                + '\n' + FrameTxt
                                
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        if AddTxt:
            AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
            AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
            AddTxt = AddTxt.replace('\n', '')
        else:
            AddTxt = ''
            
        ExportFname = CurrentDateTime + '_PixArrs_' + AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return






def PlotPixArrsFromListOfSegs_v3(ListOfSegs, ListOfDicomDirs, ListOfPlotTitles,
                                 ExportPlot=False, ExportDir='cwd', 
                                 TxtToAddToFname='', LogToConsole=False):
    """ 
    02/06/2021
    
    Note:
        
    Modification of v1. Rather than plotting 
    AllSliceNums = UniqueItems(Items=ListOfF2SindsBySeg)
    will just plot frames with segmentations since F2Sinds for Src and Trg can 
    relate to very different depths, leading to a confusing plot.
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    from pathlib import Path
    #import time
    #from copy import deepcopy
    from pydicom import dcmread
    #from DicomTools import GetDicomFpaths
    from SegTools import GetSegDataFromListOfSegs
    from GeneralTools import UniqueItems#, Unpack, PrintTitle
    
    ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
    ListOfIPPs, ListOfDirections, ListOfSpacings\
    = GetSegDataFromListOfSegs(ListOfSegs, ListOfDicomDirs, LogToConsole)
    
    # Get the maximum number of slices containing segmentations in any ROI in 
    # any dataset:
    MaxNumSlices = 0
    
    # The list (for each dataset) of the unique set of slice numbers that 
    # correspond to all segmentations in all ROIs:
    ListOfUniqueSinds = []
    
    # Loop through each dataset:
    for i in range(len(ListOfF2SindsBySeg)):
        F2SindsBySeg = ListOfF2SindsBySeg[i]
        
        UniqueSinds = UniqueItems(F2SindsBySeg)
        
        ListOfUniqueSinds.append(UniqueSinds)
        
        if len(UniqueSinds) > MaxNumSlices:
            MaxNumSlices = len(UniqueSinds)
    
    
    print(f'ListOfF2SindsBySeg = {ListOfF2SindsBySeg}')
    print(f'ListOfUniqueSinds = {ListOfUniqueSinds}')
    print(f'MaxNumSlices = {MaxNumSlices}')
    
    
    """ Prepare the figure. """
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfSegs)
    Nrows = MaxNumSlices
    
    Cmaps = ['Reds', 'Blues', 'Greens', 'Oranges', 'Purples']
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    DcmAlpha = 0.2
    DcmAlpha = 0.5
    #DcmAlpha = 1
    
    SegAlpha = 0.5
    #SegAlpha = 0.2
    
    if ExportPlot:
        dpi = 120
    else:
        dpi = 80
        
    n = 1 # initialised sub-plot number
    
    if Ncols < 3:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(7*Ncols, 9*Nrows), 
                               dpi=dpi)
    else:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(4*Ncols, 7.5*Nrows), 
                               dpi=dpi)
    
    # Loop through each slice number (i.e. each row in the plot):
    for rowNum in range(MaxNumSlices):
        # Loop through each data set:
        for i in range(len(ListOfPixArrBySeg)):
            PixArrBySeg = ListOfPixArrBySeg[i]
            F2SindsBySeg = ListOfF2SindsBySeg[i]
            UniqueSinds = ListOfUniqueSinds[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            IPPs = ListOfIPPs[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            #DicomDir = ListOfDicomDirs[i]
            
            if LogToConsole:
                print(f'\n   i = {i}')
                print(f'   PlotTitle = {PlotTitle}')
                print(f'   F2SindsBySeg = {F2SindsBySeg}')
                print(f'   UniqueSinds = {UniqueSinds}')
            
            # Continue if rowNum does not exceed the number of unique slice indeces:
            if rowNum < len(UniqueSinds):
                Sind = UniqueSinds[rowNum]
                
                ax = plt.subplot(Nrows, Ncols, n)
                
                DicomPixArr = dcmread(DicomFpaths[Sind]).pixel_array
            
                #im = ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r, alpha=DcmAlpha)
                
                IPP = IPPs[Sind]
                    
                # Loop through each segment:
                for s in range(len(F2SindsBySeg)):
                    F2Sinds = F2SindsBySeg[s]
                    PixArr = PixArrBySeg[s]
                    
                    """ There are only len(Cmaps) colormaps defined above. Wrap the 
                    segment index s if there are more segments than the number of 
                    defined colormaps. """
                    if s < len(Cmaps):
                        Cmap = Cmaps[s]
                    else:
                        m = s//len(Cmaps)
                        
                        Cmap = Cmaps[s - m*len(Cmaps)]
                        
                    if Sind in F2Sinds:
                        FrameNums = [i for i, e in enumerate(F2Sinds) if e==Sind]
                        
                        # Loop through all frame numbers:
                        for f in range(len(FrameNums)):
                            FrameNum = FrameNums[f]
                        
                            Frame = PixArr[FrameNum]
                            
                            if LogToConsole:
                                print(f'      SliceNum = {Sind} is in F2Sinds')
                                print(f'      FrameNum = {FrameNum}')
                                #print(f'      Contour = {Contour}')
                                print(f'      Frame.shape = {Frame.shape}')
                        
                            #ax = plt.subplot(Nrows, Ncols, n)
                            #im = ax.imshow(Frame, cmap=plt.cm.nipy_spectral, 
                            #               alpha=alpha)
                            ax.imshow(Frame, cmap=Cmap, alpha=SegAlpha)
                            
                            #FrameTxt = f'Frame {FrameNum}'
                        if len(FrameNums) > 1:
                            FrameTxt = f'Frames {FrameNums}'
                        else:
                            FrameTxt = f'Frame {FrameNum}'
                            
                    else:
                        FrameTxt = 'No frame'
                        
                        if LogToConsole:
                            print(f'      Sind = {Sind} is NOT in F2Sinds')
                        
                    SliceTxt = f'Slice {Sind}'
                    zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    PlotTitle = PlotTitle.replace(' ', '\:')
                    plotTitle = r"$\bf{{{x}}}$".format(x=PlotTitle) \
                                + f'\n\n{SliceTxt}\n{FrameTxt}\n{zPosTxt}'
                    #plotTitle = r"$\bf{PlotTitle}$\," \
                    #            + f'\n\n{SliceTxt}\n{FrameTxt}\n{zPosTxt}'
                    
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    ax.set_title(plotTitle)
                    
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        import time
        #import os
        #from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
        
    return





    
def PlotContoursFromListOfRtss_v1_OLD(ListOfRtss, ListOfDicomDirs, ListOfPlotTitles,
                                  SearchString, PlotAllSlices=False, AddTxt=None, 
                                  ExportPlot=False, ExportDir=None, dpi=80, 
                                  LogToConsole=False):
    """ 
    Note:
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from RtsTools import GetRtsDataFromListOfRtss
    from GeneralTools import Unpack
    from ConversionTools import Points2Indices
    from ImageTools import ImportImage
    
    ListOfPtsByCnt, ListOfC2Sinds, ListOfRoiNums,\
    ListOfDcmFpaths, ListOfOrigins, ListOfDirections,\
    ListOfSpacings = GetRtsDataFromListOfRtss(ListOfRtss, ListOfDicomDirs, 
                                              SearchString, LogToConsole)
    
    # Get the list of slice numbers that will be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices)))
            
        
    else:
        if len(ListOfC2Sinds) > 1:
            # Get the slice numbers that cover all CSToSliceInds in 
            # ListOfSliceToSliceInds:
            AllSliceNums = []
            
            for C2Sinds in ListOfC2Sinds:
                if C2Sinds:
                    for SliceInd in C2Sinds:
                        AllSliceNums.append(SliceInd)
            
            # Reduce to unique slice numbers:
            AllSliceNums = list(set(AllSliceNums))    
        else:
            AllSliceNums = list(set(deepcopy(ListOfC2Sinds[0])))
    
    if LogToConsole:
            print(f'\nAllSliceNums = {AllSliceNums}')
            
            
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfRtss)
    Nrows = len(AllSliceNums)
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7.5*Nrows), dpi=dpi)
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
    
    Colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
    n = 1 # initialised sub-plot number
    
    # Loop through each AllSliceNums: 
    for SliceNum in AllSliceNums:
        if LogToConsole:
            print(f'\nSliceNum = {SliceNum}')
        
        # Loop through each list of contours:
        for i in range(len(ListOfPtsByCnt)):
            #Contours = ListOfContours[i]
            PtsByCnt = ListOfPtsByCnt[i]
            #PtsByCntBySlice = ListOfPtsByCntBySlice[i]
            C2Sinds = ListOfC2Sinds[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            Origin = ListOfOrigins[i]
            Directions = ListOfDirections[i]
            Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            DicomDir = ListOfDicomDirs[i]
            
            RefIm = ImportImage(DicomDir)
            
            if LogToConsole:
                    print(f'\n   i = {i}')
                    print(f'   PlotTitle = {PlotTitle}')
                    #print(f'   len(PtsByCntBySlice) = {len(PtsByCntBySlice)}')
                    print(f'   C2Sinds = {C2Sinds}')
                    
            ax = plt.subplot(Nrows, Ncols, n)
            
            if SliceNum < len(DicomFpaths):
                DicomPixArr = dcmread(DicomFpaths[SliceNum]).pixel_array
            
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
            
            
                if SliceNum in C2Sinds:
                    #ContourNum = C2Sinds.index(SliceNum) # there may be more 
                    # than one contour for SliceNum!!
                    ContourNums = [i for i, e in enumerate(C2Sinds) if e==SliceNum]
                    
                    for c in range(len(ContourNums)):
                        ContourNum = ContourNums[c]
                        
                        Pts = PtsByCnt[ContourNum]
                        
                        #print(f'\n\n\nPts = {Pts}')
                        
                        #Indices = Points2Indices(Points=Contour, Origin=Origin, 
                        #                          Directions=Directions, 
                        #                          Spacings=Spacings)
                        Inds = Points2Indices(Points=Pts, RefIm=RefIm)
                        
                        if LogToConsole:
                            print(f'      SliceNum = {SliceNum} is in C2Sinds')
                            print(f'      ContourNum = {ContourNum}')
                            #print(f'      Contour = {Contour}')
                            print(f'      len(Pts) = {len(Pts)}')
                        
                        #ax = plt.subplot(Nrows, Ncols, n)
                        
                        # Unpack the contour points' indices:
                        X, Y, Z = Unpack(Inds)
                        
                        """Since there are only len(Colours) defined above,
                        wrap the index c if there are more contours than the 
                        number of defined colours."""
                        if c < len(Colours):
                            colour = Colours[c]
                        else:
                            m = c//len(Colours)
                            
                            colour = Colours[c - m*len(Colours)]
                            
                        # Plot the contour points:
                        #ax = plt.plot(X, Y, linewidth=0.5, c='r')
                        ax.plot(X, Y, linewidth=0.5, c='r')
                        #ax.plot(X, Y, linewidth=0.5, c=colour)
                        
                        #ContourTxt = f'Contour {ContourNum + 1}'
                        #ContourTxt = f'Contour {ContourNum}'
                        ContourTxt = f'Contour(s) {ContourNums}'
                    
                else:
                    ContourTxt = 'No contour'
                    
                    if LogToConsole:
                        print(f'      SliceNum = {SliceNum} is NOT in C2Sinds')
                
                SliceTxt = f'Slice {SliceNum}'
                
                if AddTxt:
                    plotTitle = PlotTitle \
                                + '\n' + SliceTxt + '\n'\
                                + '\n' + ContourTxt + '\n'\
                                + AddTxt
                else:
                    plotTitle = PlotTitle \
                                + '\n' + SliceTxt + '\n'\
                                + '\n' + ContourTxt
                
                ax.set_xlabel('pixels'); ax.set_ylabel('pixels')
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        if AddTxt:
            AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
            AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
            AddTxt = AddTxt.replace('\n', '')
        else:
            AddTxt = ''
            
        ExportFname = CurrentDateTime + '_Contours_' + AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return








def PlotContoursFromListOfRtss_v1(ListOfRtss, ListOfDicomDirs, ListOfPlotTitles,
                                  PlotAllSlices=False, ExportPlot=False, 
                                  ExportDir=None, ExportFname=None, dpi=80, 
                                  LogToConsole=False):
    """ 
    Note:
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import DicomTools
    importlib.reload(DicomTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    
    import matplotlib.pyplot as plt
    import os
    from pathlib import Path
    #import time
    #from copy import deepcopy
    from pydicom import dcmread
    from RtsTools import GetRtsDataFromListOfRtss
    from GeneralTools import UniqueItems, Unpack, PrintTitle
    from ConversionTools import Points2Indices
    from ImageTools import ImportImage
    
    ListOfPtsByCntByRoi, ListOfC2SindsByRoi, ListOfDcmFpaths,\
    ListOfIPPs, ListOfDirections, ListOfSpacings\
    = GetRtsDataFromListOfRtss(ListOfRtss, ListOfDicomDirs, LogToConsole)
    
    # Get the list of slice numbers that will be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices))) 
    else:
        AllSliceNums = UniqueItems(Items=ListOfC2SindsByRoi)
    
    if LogToConsole:
        PrintTitle('Running PlotContoursFromListOfRtss_v1()')
        print(f'ListOfC2SindsByRoi = {ListOfC2SindsByRoi}')
        print(f'AllSliceNums = {AllSliceNums}')
            
            
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfRtss)
    Nrows = len(AllSliceNums)
    
    LineWidth = 0.5
    LineWidth = 2
    Colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
    # Set the transparency of the DICOM images:
    #alpha = 0.2
    #alpha = 0.5
    alpha = 1
    
    n = 1 # initialised sub-plot number
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7.5*Nrows), dpi=dpi)
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi) # figure size too big
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(4*Ncols, 7.5*Nrows), dpi=dpi)
    
    # Loop through each AllSliceNums: 
    for SliceNum in AllSliceNums:
        if LogToConsole:
            print(f'\nSliceNum = {SliceNum}')
            
        # Loop through each dataset:
        for i in range(len(ListOfPtsByCntByRoi)):
            PtsByCntByRoi = ListOfPtsByCntByRoi[i]
            C2SindsByRoi = ListOfC2SindsByRoi[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            IPPs = ListOfIPPs[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            DicomDir = ListOfDicomDirs[i]
            
            RefIm = ImportImage(DicomDir)
            
            if LogToConsole:
                    print(f'\n   i = {i}')
                    print(f'   PlotTitle = {PlotTitle}')
                    print(f'   C2SindsByRoi = {C2SindsByRoi}')
                    
            ax = plt.subplot(Nrows, Ncols, n)
            
            if SliceNum < len(DicomFpaths):
                DicomPixArr = dcmread(DicomFpaths[SliceNum]).pixel_array
            
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r, alpha=alpha)
                
                IPP = IPPs[SliceNum]
                
            
            # Loop through each ROI:
            for r in range(len(C2SindsByRoi)):
                C2Sinds = C2SindsByRoi[r]
                PtsByCnt = PtsByCntByRoi[r]
                
                """ There are only len(Colours) colours defined above. Wrap the 
                ROI index r if there are more ROIs than the number of defined 
                colours. """
                if r < len(Colours):
                    colour = Colours[r]
                else:
                    m = r//len(Colours)
                    
                    colour = Colours[r - m*len(Colours)]
            
                if SliceNum in C2Sinds:
                    #ContourNum = C2Sinds.index(SliceNum) # there may be more 
                    # than one contour for SliceNum!!
                    ContourNums = [i for i, e in enumerate(C2Sinds) if e==SliceNum]
                    
                    # Loop through all contour numbers:
                    for c in range(len(ContourNums)):
                        ContourNum = ContourNums[c]
                        
                        Pts = PtsByCnt[ContourNum]
                        
                        #print(f'\n\n\nPts = {Pts}')
                        
                        #Indices = Points2Indices(Points=Contour, Origin=Origin, 
                        #                          Directions=Directions, 
                        #                          Spacings=Spacings)
                        Inds = Points2Indices(Points=Pts, RefIm=RefIm)
                        
                        if LogToConsole:
                            print(f'      SliceNum = {SliceNum} is in C2Sinds')
                            print(f'      ContourNum = {ContourNum}')
                            #print(f'      Contour = {Contour}')
                            print(f'      len(Pts) = {len(Pts)}')
                        
                        #ax = plt.subplot(Nrows, Ncols, n)
                        
                        # Unpack the contour points' indices:
                        X, Y, Z = Unpack(Inds)
                        
                        #"""Since there are only len(Colours) defined above,
                        #wrap the index c if there are more contours than the 
                        #number of defined colours."""
                        #if c < len(Colours):
                        #    colour = Colours[c]
                        #else:
                        #    m = c//len(Colours)
                        #    
                        #    colour = Colours[c - m*len(Colours)]
                            
                        # Plot the contour points:
                        #ax = plt.plot(X, Y, linewidth=0.5, c='r')
                        #ax.plot(X, Y, linewidth=0.5, c='r')
                        ax.plot(X, Y, linewidth=LineWidth, c=colour)
                        
                        #ContourTxt = f'Contour {ContourNum + 1}'
                        #ContourTxt = f'Contour {ContourNum}'
                        #ContourTxt = f'Contour(s) {ContourNums}'
                        ContourTxt = f'Contour {ContourNum}'
                    
                else:
                    ContourTxt = 'No contour'
                    
                    if LogToConsole:
                        print(f'      SliceNum = {SliceNum} is NOT in C2Sinds')
                
                SliceTxt = f'Slice {SliceNum}'
                zPosTxt = f'z = {round(IPP[2], 2)} mm'
                
                #if AddTxt:
                #    plotTitle = PlotTitle \
                #                + f'\n' + SliceTxt + '\n'\
                #                + f'\n' + ContourTxt + '\n'\
                #                + AddTxt
                #else:
                #    plotTitle = PlotTitle \
                #                + f'\n' + SliceTxt + '\n'\
                #                + f'\n' + ContourTxt
                
                plotTitle = f'{PlotTitle}\n\n{SliceTxt}\n{ContourTxt}\n{zPosTxt}'
                
                #if AddTxt:
                #    plotTitle += AddTxt
                
                ax.set_xlabel('pixels'); ax.set_ylabel('pixels')
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        #CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        #CurrentTime = time.strftime("%H%M%S", time.gmtime())
        #CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        #if AddTxt:
        #    AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
        #    AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
        #    AddTxt = AddTxt.replace('\n', '')
        #else:
        #    AddTxt = ''
            
        ##ExportFname = CurrentDateTime + '_Contours_' + AddTxt + '.jpg'
        #ExportFname = AddTxt + '.jpg'
        
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
        
    return










def PlotContoursFromListOfRtss_v2(ListOfRtss, ListOfDicomDirs, ListOfPlotTitles,
                                  SearchString, PlotAllSlices=False, AddTxt=None, 
                                  ExportPlot=False, ExportDir=None, dpi=80, 
                                  LogToConsole=False):
    """ 
    Note:
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG ROWS and different SLICES ALONG COLUMNS. 
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from GeneralTools import Unpack
    from RtsTools import GetRtsDataFromListOfRtss
    from ConversionTools import Points2Indices
    from ImageTools import ImportImage
    
    ListOfPtsByCnt, ListOfC2Sinds, ListOfRoiNums,\
    ListOfDcmFpaths, ListOfOrigins, ListOfDirections,\
    ListOfSpacings = GetRtsDataFromListOfRtss(ListOfRtss, ListOfDicomDirs, 
                                               SearchString, LogToConsole)
    
    # Get the list of slice numbers that will be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices)))
            
        
    else:
        if len(ListOfC2Sinds) > 1:
            # Get the slice numbers that cover all CSToSliceInds in 
            # ListOfSliceToSliceInds:
            AllSliceNums = []
            
            for C2Sinds in ListOfC2Sinds:
                if C2Sinds:
                    for SliceInd in C2Sinds:
                        AllSliceNums.append(SliceInd)
            
            # Reduce to unique slice numbers:
            AllSliceNums = list(set(AllSliceNums))    
        else:
            AllSliceNums = list(set(deepcopy(ListOfC2Sinds[0])))
    
    if True:#LogToConsole:
            print(f'\nAllSliceNums = {AllSliceNums}')
            
            
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Nrows = len(ListOfRtss)
    Ncols = len(AllSliceNums)
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7.5*Nrows), dpi=dpi)
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    Colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
    n = 1 # initialised sub-plot number
        
    # Loop through each dataset:
    for i in range(len(ListOfPtsByCnt)):
        #Contours = ListOfContours[i]
        PtsByCnt = ListOfPtsByCnt[i]
        C2Sinds = ListOfC2Sinds[i]
        
        DicomFpaths = ListOfDcmFpaths[i]
        Origin = ListOfOrigins[i]
        Directions = ListOfDirections[i]
        Spacings = ListOfSpacings[i]
        PlotTitle = ListOfPlotTitles[i]
        DicomDir = ListOfDicomDirs[i]
        
        RefIm = ImportImage(DicomDir)
        
        if LogToConsole:
                print(f'\ni = {i}')
                print(f'   PlotTitle = {PlotTitle}')
                print(f'   C2Sinds = {C2Sinds}')
                if C2Sinds:
                    print(f'   len(PtsByCnt) = {len(PtsByCnt)}')
        
        # Loop through each AllSliceNums: 
        for SliceNum in AllSliceNums:
            if LogToConsole:
                print(f'\n   SliceNum = {SliceNum}')
                
            ax = plt.subplot(Nrows, Ncols, n)
            
            if SliceNum < len(DicomFpaths):
                DicomPixArr = dcmread(DicomFpaths[SliceNum]).pixel_array
            
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
            
                if C2Sinds:
                    if SliceNum in C2Sinds:
                        #ContourNum = C2Sinds.index(SliceNum) # there may be more 
                        # than one contour for SliceNum!!
                        ContourNums = [i for i, e in enumerate(C2Sinds) if e==SliceNum]
                        
                        for c in range(len(ContourNums)):
                            ContourNum = ContourNums[c]
                            
                            #Contour = Contours[ContourNum]
                            Pts = PtsByCnt[ContourNum]
                            
                            #Indices = Points2Indices(Points=Contour, Origin=Origin, 
                            #                          Directions=Directions, 
                            #                          Spacings=Spacings)
                            Inds = Points2Indices(Points=Pts, RefIm=RefIm)
                            
                            if LogToConsole:
                                print(f'   SliceNum = {SliceNum} is in C2Sinds')
                                print(f'   ContourNum = {ContourNum}')
                                #print(f'   len(Contour) = {len(Contour)}')
                                #print(f'      Contour = {Contour}')
                                print(f'      len(Pts) = {len(Pts)}')
                            
                            #ax = plt.subplot(Nrows, Ncols, n)
                            
                            # Unpack the contour points' indices:
                            X, Y, Z = Unpack(Inds)
                            
                            """Since there are only len(Colours) defined above,
                            wrap the index c if there are more contours than the 
                            number of defined colours."""
                            if c < len(Colours):
                                colour = Colours[c]
                            else:
                                m = c//len(Colours)
                                
                                colour = Colours[c - m*len(Colours)]
                                
                            # Plot the contour points:
                            #ax = plt.plot(X, Y, linewidth=0.5, c='r');
                            ax.plot(X, Y, linewidth=0.5, c='r');
                            #ax.plot(X, Y, linewidth=0.5, c=colour)
                        
                            #ContourTxt = f'Contour {ContourNum + 1}'
                            #ContourTxt = f'Contour {ContourNum}'
                            ContourTxt = f'Contour(s) {ContourNums}'
                    else:
                        ContourTxt = 'No contour'
                        
                        if LogToConsole:
                            print(f'   SliceNum = {SliceNum} is NOT in C2Sinds')
                else:
                    ContourTxt = 'No contour'
                
                SliceTxt = f'Slice {SliceNum}'
                SliceTxt = ''
                
                if AddTxt:
                    plotTitle = PlotTitle \
                                + '\n' + SliceTxt + '\n'\
                                + '\n' + ContourTxt + '\n'\
                                + AddTxt
                else:
                    plotTitle = PlotTitle \
                                + '\n' + SliceTxt + '\n'\
                                + '\n' + ContourTxt
                                        
                ax.set_xlabel('pixels'); ax.set_ylabel('pixels')
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        if AddTxt:
            AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
            AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
            AddTxt = AddTxt.replace('\n', '')
        else:
            AddTxt = ''
            
        ExportFname = CurrentDateTime + '_Contours_' + AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return





def PlotContoursFromListOfRtss_v3(ListOfRtss, ListOfDicomDirs, ListOfPlotTitles,
                                  PlotAllSlices=False, AddTxt=None, 
                                  ExportPlot=False, ExportDir=None, dpi=80, 
                                  LogToConsole=False):
    """ 
    Note:
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import DicomTools
    importlib.reload(DicomTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    
    import matplotlib.pyplot as plt
    import os
    import time
    #from copy import deepcopy
    from pydicom import dcmread
    from RtsTools import GetRtsDataFromListOfRtss
    from GeneralTools import UniqueItems, Unpack, PrintTitle
    from ConversionTools import Points2Indices
    from ImageTools import ImportImage
    
    ListOfPtsByCntByRoi, ListOfC2SindsByRoi, ListOfDcmFpaths,\
    ListOfIPPs, ListOfDirections, ListOfSpacings\
    = GetRtsDataFromListOfRtss(ListOfRtss, ListOfDicomDirs, LogToConsole)
    
    ListOfOrigins = [IPPs[0] for IPPs in ListOfIPPs]
    
    MinZ = min([Origin[2] for Origin in ListOfOrigins])
    
    MaxZ = max([IPPs[-1][2] for IPPs in ListOfIPPs])
    
    ListOfLengths = [len(items) for items in ListOfIPPs]
    
    UniqueLengths = UniqueItems(ListOfLengths)
    
    
    
    
    # Get the list of slice numbers that will be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices))) 
    else:
        AllSliceNums = UniqueItems(Items=ListOfC2SindsByRoi)
    
    if LogToConsole:
        PrintTitle('Running PlotContoursFromListOfRtss_v1()')
        print(f'ListOfC2SindsByRoi = {ListOfC2SindsByRoi}')
        print(f'AllSliceNums = {AllSliceNums}')
            
            
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfRtss)
    Nrows = len(AllSliceNums)
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7.5*Nrows), dpi=dpi)
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
    
    Colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
    n = 1 # initialised sub-plot number
    
    # Loop through each AllSliceNums: 
    for SliceNum in AllSliceNums:
        if LogToConsole:
            print(f'\nSliceNum = {SliceNum}')
            
        # Loop through each dataset:
        for i in range(len(ListOfPtsByCntByRoi)):
            PtsByCntByRoi = ListOfPtsByCntByRoi[i]
            C2SindsByRoi = ListOfC2SindsByRoi[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            DicomDir = ListOfDicomDirs[i]
            
            RefIm = ImportImage(DicomDir)
            
            if LogToConsole:
                    print(f'\n   i = {i}')
                    print(f'   PlotTitle = {PlotTitle}')
                    print(f'   C2SindsByRoi = {C2SindsByRoi}')
                    
            ax = plt.subplot(Nrows, Ncols, n)
            
            if SliceNum < len(DicomFpaths):
                DicomPixArr = dcmread(DicomFpaths[SliceNum]).pixel_array
            
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
            
            
            # Loop through each ROI:
            for r in range(len(C2SindsByRoi)):
                C2Sinds = C2SindsByRoi[r]
                PtsByCnt = PtsByCntByRoi[r]
                
                # There are only len(Colours) defined above. Wrap the index c
                # if there are more contours than the number of defined colours
                if r < len(Colours):
                    colour = Colours[r]
                else:
                    m = r//len(Colours)
                    
                    colour = Colours[r - m*len(Colours)]
            
                if SliceNum in C2Sinds:
                    #ContourNum = C2Sinds.index(SliceNum) # there may be more 
                    # than one contour for SliceNum!!
                    ContourNums = [i for i, e in enumerate(C2Sinds) if e==SliceNum]
                    
                    for c in range(len(ContourNums)):
                        ContourNum = ContourNums[c]
                        
                        Pts = PtsByCnt[ContourNum]
                        
                        #print(f'\n\n\nPts = {Pts}')
                        
                        #Indices = Points2Indices(Points=Contour, Origin=Origin, 
                        #                          Directions=Directions, 
                        #                          Spacings=Spacings)
                        Inds = Points2Indices(Points=Pts, RefIm=RefIm)
                        
                        if LogToConsole:
                            print(f'      SliceNum = {SliceNum} is in C2Sinds')
                            print(f'      ContourNum = {ContourNum}')
                            #print(f'      Contour = {Contour}')
                            print(f'      len(Pts) = {len(Pts)}')
                        
                        #ax = plt.subplot(Nrows, Ncols, n)
                        
                        # Unpack the contour points' indices:
                        X, Y, Z = Unpack(Inds)
                        
                        #"""Since there are only len(Colours) defined above,
                        #wrap the index c if there are more contours than the 
                        #number of defined colours."""
                        #if c < len(Colours):
                        #    colour = Colours[c]
                        #else:
                        #    m = c//len(Colours)
                        #    
                        #    colour = Colours[c - m*len(Colours)]
                            
                        # Plot the contour points:
                        #ax = plt.plot(X, Y, linewidth=0.5, c='r')
                        #ax.plot(X, Y, linewidth=0.5, c='r')
                        ax.plot(X, Y, linewidth=0.5, c=colour)
                        
                        #ContourTxt = f'Contour {ContourNum + 1}'
                        #ContourTxt = f'Contour {ContourNum}'
                        ContourTxt = f'Contour(s) {ContourNums}'
                    
                else:
                    ContourTxt = 'No contour'
                    
                    if LogToConsole:
                        print(f'      SliceNum = {SliceNum} is NOT in C2Sinds')
                
                SliceTxt = f'Slice {SliceNum}'
                
                if AddTxt:
                    plotTitle = PlotTitle \
                                + '\n' + SliceTxt + '\n'\
                                + '\n' + ContourTxt + '\n'\
                                + AddTxt
                else:
                    plotTitle = PlotTitle \
                                + '\n' + SliceTxt + '\n'\
                                + '\n' + ContourTxt
                
                ax.set_xlabel('pixels'); ax.set_ylabel('pixels')
                ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        if AddTxt:
            AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
            AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
            AddTxt = AddTxt.replace('\n', '')
        else:
            AddTxt = ''
            
        ExportFname = CurrentDateTime + '_Contours_' + AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return








def PlotContoursFromListOfRtss_v4(ListOfRtss, ListOfDicomDirs, ListOfPlotTitles,
                                  ExportPlot=False, ExportDir='cwd',
                                  TxtToAddToFname='', LogToConsole=False):
    """ 
    02/06/2021
    
    Note:
        
    Modification of v1. Rather than plotting 
    AllSliceNums = UniqueItems(Items=ListOfC2SindsByRoi)
    will just plot frames with contours since C2Sinds for Src and Trg can 
    relate to very different depths, leading to a confusing plot.
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import DicomTools
    importlib.reload(DicomTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    
    import matplotlib.pyplot as plt
    import os
    from pathlib import Path
    #import time
    #from copy import deepcopy
    from pydicom import dcmread
    from RtsTools import GetRtsDataFromListOfRtss
    from GeneralTools import UniqueItems, Unpack, PrintTitle
    from ConversionTools import Points2Indices
    from ImageTools import ImportImage
    
    ListOfPtsByCntByRoi, ListOfC2SindsByRoi, ListOfDcmFpaths,\
    ListOfIPPs, ListOfDirections, ListOfSpacings\
    = GetRtsDataFromListOfRtss(ListOfRtss, ListOfDicomDirs, LogToConsole)
    
    if LogToConsole:
        PrintTitle('Running PlotContoursFromListOfRtss_v1()')
        print(f'ListOfC2SindsByRoi = {ListOfC2SindsByRoi}')
            
    
    # Initialise the maximum number of slices containing contours in any ROI in  
    # any dataset:
    MaxNumSlices = 0
    
    # The list (for each dataset) of the unique set of slice numbers that 
    # correspond to all contours in all ROIs:
    ListOfUniqueSinds = []
    
    # Loop through each dataset:
    for i in range(len(ListOfC2SindsByRoi)):
        C2SindsByRoi = ListOfC2SindsByRoi[i]
        
        UniqueSinds = UniqueItems(C2SindsByRoi)
        
        ListOfUniqueSinds.append(UniqueSinds)
        
        if len(UniqueSinds) > MaxNumSlices:
            MaxNumSlices = len(UniqueSinds)
    
    
    print(f'ListOfC2SindsByRoi = {ListOfC2SindsByRoi}')
    print(f'ListOfUniqueSinds = {ListOfUniqueSinds}')
    print(f'MaxNumSlices = {MaxNumSlices}')
    
    """ Prepare the figure. """
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfRtss)
    Nrows = MaxNumSlices
    
    LineWidth = 0.5
    LineWidth = 2
    Colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
    # Set the transparency of the DICOM images:
    #alpha = 0.2
    #alpha = 0.5
    alpha = 1
    
    if ExportPlot:
        dpi = 120
    else:
        dpi = 80
    
    n = 1 # initialised sub-plot number
    
    if Ncols < 3:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(7*Ncols, 9*Nrows), 
                               dpi=dpi)
    else:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(4*Ncols, 7.5*Nrows), 
                               dpi=dpi)
    
    
    # Loop through each slice/frame number (i.e. each row in the plot):
    for rowNum in range(MaxNumSlices):
        # Loop through each dataset:
        for i in range(Ncols):
            PtsByCntByRoi = ListOfPtsByCntByRoi[i]
            C2SindsByRoi = ListOfC2SindsByRoi[i]
            UniqueSinds = ListOfUniqueSinds[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            IPPs = ListOfIPPs[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            DicomDir = ListOfDicomDirs[i]
            
            RefIm = ImportImage(DicomDir)
            
            if LogToConsole:
                print(f'\n   i = {i}')
                print(f'   PlotTitle = {PlotTitle}')
                print(f'   C2SindsByRoi = {C2SindsByRoi}')
                print(f'   UniqueSinds = {UniqueSinds}')
            
            # Continue if rowNum does not exceed the number of unique slice indices:
            if rowNum < len(UniqueSinds):
                Sind = UniqueSinds[rowNum]
                
                ax = plt.subplot(Nrows, Ncols, n)
                
                DicomPixArr = dcmread(DicomFpaths[Sind]).pixel_array
            
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r, alpha=alpha)
                
                IPP = IPPs[Sind]
            
                # Loop through each ROI:
                for r in range(len(C2SindsByRoi)):
                    C2Sinds = C2SindsByRoi[r]
                    PtsByCnt = PtsByCntByRoi[r]
                    
                    """ There are only len(Colours) colours defined above. Wrap the 
                    ROI index r if there are more ROIs than the number of defined 
                    colours. """
                    if r < len(Colours):
                        colour = Colours[r]
                    else:
                        m = r//len(Colours)
                        
                        colour = Colours[r - m*len(Colours)]
                
                    if Sind in C2Sinds:
                        #ContourNum = C2Sinds.index(SliceNum) # there may be more 
                        # than one contour for SliceNum!!
                        ContourNums = [i for i, e in enumerate(C2Sinds) if e==Sind]
                        
                        #print(f'ContourNums = {ContourNums}')
                        
                        NumPtsByCnt = []
                        
                        # Loop through all contour numbers:
                        for c in range(len(ContourNums)):
                            ContourNum = ContourNums[c]
                            
                            Pts = PtsByCnt[ContourNum]
                            
                            NumPtsByCnt.append(len(Pts))
                            
                            #print(f'\n\n\nPts = {Pts}')
                            
                            Inds = Points2Indices(Points=Pts, RefIm=RefIm)
                            
                            if LogToConsole:
                                print(f'      Sind = {Sind} is in C2Sinds')
                                print(f'      ContourNum = {ContourNum}')
                                #print(f'      Contour = {Contour}')
                                print(f'      len(Pts) = {len(Pts)}')
                            
                            #ax = plt.subplot(Nrows, Ncols, n)
                            
                            # Unpack the contour points' indices:
                            X, Y, Z = Unpack(Inds)
                                
                            # Plot the contour points:
                            #ax = plt.plot(X, Y, linewidth=0.5, c='r')
                            #ax.plot(X, Y, linewidth=0.5, c='r')
                            ax.plot(X, Y, linewidth=LineWidth, c=colour)
                        
                        #if len(ContourNums) > 1:
                        #    ContourTxt = f'Contours {ContourNums}'
                        #else:
                        #    ContourTxt = f'Contour {ContourNum}'
                        #ContourTxt = f'Contour # {ContourNums}'
                        nobrackets = f'{ContourNums}'.replace('[', '').replace(']', '')
                        ContourTxt = f'Contour # {nobrackets}'
                        
                        nobrackets = f'{NumPtsByCnt}'.replace('[', '').replace(']', '')
                        PtsTxt = f'No. of points: {nobrackets}'
                    else:
                        ContourTxt = 'No contour'
                        PtsTxt = 'No. of points: 0'
                        
                        if LogToConsole:
                            print(f'      Sind = {Sind} is NOT in C2Sinds')
                    
                    SliceTxt = f'Slice # {Sind}'
                    zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    PlotTitle = PlotTitle.replace(' ', '\:')
                    plotTitle = r"$\bf{{{x}}}$".format(x=PlotTitle) \
                                + f'\n\n{SliceTxt}\n{ContourTxt}\n{PtsTxt}\n{zPosTxt}'
                    
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    ax.set_title(plotTitle)
            
            n += 1 # increment sub-plot number
    
    
    if ExportPlot:
        import time
        #import os
        #from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
        
    return











def GetPixArrsFromListOfImages(ListOfImages, ListOfDicomDirs, 
                               LogToConsole=False):
    """ 
    Note:
    
    ListOfImages is a list (which could be of length 1) of SimpleITK images
    (e.g. labelmaps).
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import numpy as np
    from DicomTools import GetDicomFpaths
    
    from ConversionTools import Image2PixArr
    
    ListOfPixArrs = []
    ListOfF2Sinds = []
    ListOfDcmFpaths = []
    
    for i in range(len(ListOfImages)):
        if LogToConsole:
            print(f'\n\nListOfSegs[{i}]:')
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        ListOfDcmFpaths.append(DcmFpaths)
        
        
        if ListOfImages[i]:
            PixArr, F2Sinds = Image2PixArr(ListOfImages[i])
            
            ListOfF2Sinds.append(F2Sinds)
            
            if LogToConsole:
                print(f'\nListOfF2Sinds[{i}] = {ListOfF2Sinds[i]}')
            
            
            if LogToConsole:      
                print(f'\nPixArr.shape = {PixArr.shape}')
                [print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
                
            """
            If PixArr has a single frame it will have shape (NumOfRows, NumOfCols).
            For compatibility with multi-framed pixel arrays, convert single-framed
            pixel arrays to shape (1, NumOfRows, NumOfCols).
            """
            if len(PixArr.shape) < 3:
                R, C = PixArr.shape
                
                PixArr = np.reshape(PixArr, (1, R, C))
                
                if LogToConsole:      
                    print(f'\nPixArr was a single frame so reshaped to {PixArr.shape}')
                    
        else:
            ListOfF2Sinds.append(None)
            PixArr = None
        
        
        ListOfPixArrs.append(PixArr)
        
        #if LogToConsole: 
        #    print('')
        #    print(f'\nListOfPixArrs[{i}].shape = {ListOfPixArrs[i].shape}')
        
    
    return ListOfPixArrs, ListOfF2Sinds, ListOfDcmFpaths







def PlotPixArrsFromListOfLabIm(ListOfLabIm, ListOfDicomDirs, ListOfPlotTitles,
                               PlotAllSlices=False, AddTxt=None, 
                               ExportPlot=False, ExportDir=None, dpi=80, 
                               LogToConsole=False):
    """ 
    Note:
    
    ListOfLabIm is a list (which could be of length 1) of binary label images
    (SimpleITK Images).
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from DicomTools import GetDicomFpaths
    
    ListOfPixArrs, ListOfF2Sinds,\
    ListOfDcmFpaths = GetPixArrsFromListOfImages(ListOfLabIm, ListOfDicomDirs, 
                                                 LogToConsole)
    
    print(f'\nlen(ListOfPixArrs) = {len(ListOfPixArrs)}')
    print(f'len(ListOfF2Sinds) = {len(ListOfF2Sinds)}')
    print(f'len(ListOfDcmFpaths) = {len(ListOfDcmFpaths)}')
    print(f'ListOfF2Sinds = {ListOfF2Sinds}')
    
    
    # Get the list of slice numbers to be plotted:
    if PlotAllSlices:
        NumOfSlices = []
        
        for i in range(len(ListOfDcmFpaths)):
            NumOfSlices.append(len(ListOfDcmFpaths[i]))
            
        AllSliceNums = list(range(max(NumOfSlices)))
            
        
    else:
        if len(ListOfF2Sinds) > 1:
            # Get the slice numbers that cover all FrameToSliceInds in 
            # ListOfFrameToSliceInds:
            AllSliceNums = []
            
            for F2Sinds in ListOfF2Sinds:
                if F2Sinds:
                    for SliceInd in F2Sinds:
                        AllSliceNums.append(SliceInd)
            
            print(f'\nAllSliceNums = {AllSliceNums}')
            
            # Reduce to unique slice numbers:
            AllSliceNums = list(set(AllSliceNums))
            
            print(f'\nAllSliceNums = {AllSliceNums}')
        else:
            AllSliceNums = deepcopy(ListOfF2Sinds[0])
    
    
    
    print(f'\nAllSliceNums = {AllSliceNums}')
    
    # Get the maxima for all frames for all pixel arrays along x (axis=2) and
    # y (axis=1):
    #ListOfMaximaByFrame = []
    #
    #for PixArr in ListOfPixArrs:
    #
    #    ListOfMaximaByFrame.append(np.amax(PixArr, axis=(1,2)))
    #
    #
    #if LogToConsole:
    #    print(f'\nListOfMaximaByFrame = {ListOfMaximaByFrame}')
    
    
    
    # Prepare the figure.
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    #alpha = 0.2
    alpha = 0.5
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfPixArrs)
    Nrows = len(AllSliceNums)
    
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 7*Nrows), dpi=dpi)
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Loop through each AllSliceNums: 
    for SliceNum in AllSliceNums:
        #print(f's = {SliceNum}')
        
        # Loop through each data set:
        for i in range(len(ListOfPixArrs)):
            PixArr = ListOfPixArrs[i]
            
            F2Sinds = ListOfF2Sinds[i]
            
            DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
            
            PlotTitle = ListOfPlotTitles[i]
            
            #MaximaByFrame = ListOfMaximaByFrame[i]
            
            ax = plt.subplot(Nrows, Ncols, n)
            
            
            if SliceNum < len(DcmFpaths):
                DcmPixArr = dcmread(DcmFpaths[SliceNum]).pixel_array
            
                #im = ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
                ax.imshow(DcmPixArr, cmap=plt.cm.Greys_r, alpha=alpha)
                
                if SliceNum in F2Sinds:
                    FrameNum = F2Sinds.index(SliceNum)
                    
                    Frame = PixArr[FrameNum]
                    
                    #ax = plt.subplot(Nrows, Ncols, n)
                    
                    #im = ax.imshow(Frame, cmap=plt.cm.nipy_spectral, 
                    #               alpha=alpha)
                    ax.imshow(Frame, cmap=plt.cm.nipy_spectral, alpha=alpha)
                    
                    #m = 1 if MaximaByFrame[FrameNum] < 1 else MaximaByFrame[FrameNum]
                    #cbar = fig.colorbar(im, ax=ax)
                    #cbar.mappable.set_clim(0, m)
                    
                    #fig.colorbar(im, ax=ax)
                    #fig.colorbar(ax=ax)
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    if AddTxt:
                        plotTitle = PlotTitle\
                                    + f' Slice {SliceNum}, Frame {FrameNum}\n'\
                                    + AddTxt
                    else:
                        plotTitle = PlotTitle\
                                    + f' Slice {SliceNum}, Frame {FrameNum}'
                    ax.set_title(plotTitle)
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        if AddTxt:
            AddTxt = AddTxt.replace(' = ', '_').replace(' ', '_')
            AddTxt = AddTxt.replace(',', '').replace('(', '').replace(')', '')
            AddTxt = AddTxt.replace('\n', '')
        else:
            AddTxt = ''
            
        ExportFname = CurrentDateTime + '_PixArrs_' + AddTxt + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return




def PlotPixArrsFromListOfLabImBySeg(ListOfLabImBySeg, ListOfDicomDirs, 
                                    ListOfPlotTitles, AddTxt=None,
                                    ExportPlot=False, ExportDir=None, 
                                    TxtToAddToFname='', LogToConsole=False):
    """ 
    08/06/21: Modelled on PlotPixArrsFromListOfSegs_v3.
    
    ListOfLabImBySeg is a list (for each dataset, which could be of length 1) 
    of a list (for each segment) of label images (SimpleITK Images).
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each image.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    #import time
    #from copy import deepcopy
    from pydicom import dcmread
    #from DicomTools import GetDicomFpaths
    from SegTools import GetSegDataFromListOfLabImBySeg
    from GeneralTools import UniqueItems
    
    ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
        ListOfIPPs, ListOfDirections, ListOfSpacings\
        = GetSegDataFromListOfLabImBySeg(ListOfLabImBySeg, 
                                         ListOfDicomDirs, 
                                         LogToConsole)
    
    print(f'\nlen(ListOfPixArrBySeg) = {len(ListOfPixArrBySeg)}')
    print(f'len(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
    print(f'len(ListOfDcmFpaths) = {len(ListOfDcmFpaths)}')
    print(f'ListOfF2SindsBySeg = {ListOfF2SindsBySeg}')
    
    
    # Get the maximum number of slices containing segmentations in any ROI in 
    # any dataset:
    MaxNumSlices = 0
    
    # The list (for each dataset) of the unique set of slice numbers that 
    # correspond to all segmentations in all ROIs:
    ListOfUniqueSinds = []
    
    # Loop through each dataset:
    for i in range(len(ListOfF2SindsBySeg)):
        F2SindsBySeg = ListOfF2SindsBySeg[i]
        
        UniqueSinds = UniqueItems(F2SindsBySeg)
        
        ListOfUniqueSinds.append(UniqueSinds)
        
        if len(UniqueSinds) > MaxNumSlices:
            MaxNumSlices = len(UniqueSinds)
    
    
    print(f'ListOfF2SindsBySeg = {ListOfF2SindsBySeg}')
    print(f'ListOfUniqueSinds = {ListOfUniqueSinds}')
    print(f'MaxNumSlices = {MaxNumSlices}')
    
    
    """ Prepare the figure. """
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfLabImBySeg)
    Nrows = MaxNumSlices
    
    Cmaps = ['Reds', 'Blues', 'Greens', 'Oranges', 'Purples']
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    DcmAlpha = 0.2
    DcmAlpha = 0.5
    #DcmAlpha = 1
    
    SegAlpha = 0.5
    #SegAlpha = 0.2
    
    if ExportPlot:
        dpi = 120
    else:
        dpi = 80
        
    n = 1 # initialised sub-plot number
    
    if Ncols < 3:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(7*Ncols, 9*Nrows), 
                               dpi=dpi)
    else:
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(4*Ncols, 7.5*Nrows), 
                               dpi=dpi)
    
    # Loop through each slice/frame number (i.e. each row in the plot):
    for rowNum in range(MaxNumSlices):
        # Loop through each data set:
        for i in range(len(ListOfPixArrBySeg)):
            PixArrBySeg = ListOfPixArrBySeg[i]
            F2SindsBySeg = ListOfF2SindsBySeg[i]
            UniqueSinds = ListOfUniqueSinds[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            IPPs = ListOfIPPs[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            #DicomDir = ListOfDicomDirs[i]
            
            if LogToConsole:
                print(f'\n   i = {i}')
                print(f'   PlotTitle = {PlotTitle}')
                print(f'   F2SindsBySeg = {F2SindsBySeg}')
                print(f'   UniqueSinds = {UniqueSinds}')
            
            # Continue if rowNum does not exceed the number of unique slice indeces:
            if rowNum < len(UniqueSinds):
                Sind = UniqueSinds[rowNum]
                
                ax = plt.subplot(Nrows, Ncols, n)
                
                DicomPixArr = dcmread(DicomFpaths[Sind]).pixel_array
            
                #im = ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r, alpha=DcmAlpha)
                
                IPP = IPPs[Sind]
                
                if LogToConsole:
                    print(f'   Sind = {Sind}')
                    print(f'   IPP = {IPP}')
                    
                # Loop through each segment:
                for s in range(len(F2SindsBySeg)):
                    F2Sinds = F2SindsBySeg[s]
                    PixArr = PixArrBySeg[s]
                    
                    """ There are only len(Cmaps) colormaps defined above. Wrap the 
                    segment index s if there are more segments than the number of 
                    defined colormaps. """
                    if s < len(Cmaps):
                        Cmap = Cmaps[s]
                    else:
                        m = s//len(Cmaps)
                        
                        Cmap = Cmaps[s - m*len(Cmaps)]
                        
                    if Sind in F2Sinds:
                        FrameNums = [i for i, e in enumerate(F2Sinds) if e==Sind]
                        
                        # Loop through all frame numbers:
                        for f in range(len(FrameNums)):
                            FrameNum = FrameNums[f]
                        
                            Frame = PixArr[FrameNum]
                            
                            if LogToConsole:
                                print(f'      SliceNum = {Sind} is in F2Sinds')
                                print(f'      FrameNum = {FrameNum}')
                                #print(f'      Contour = {Contour}')
                                print(f'      Frame.shape = {Frame.shape}')
                        
                            #ax = plt.subplot(Nrows, Ncols, n)
                            #im = ax.imshow(Frame, cmap=plt.cm.nipy_spectral, 
                            #               alpha=alpha)
                            ax.imshow(Frame, cmap=Cmap, alpha=SegAlpha)
                            
                            #FrameTxt = f'Frame {FrameNum}'
                        if len(FrameNums) > 1:
                            FrameTxt = f'Frames {FrameNums}'
                        else:
                            FrameTxt = f'Frame {FrameNum}'
                            
                    else:
                        FrameTxt = 'No frame'
                        
                        if LogToConsole:
                            print(f'      Sind = {Sind} is NOT in F2Sinds')
                        
                    SliceTxt = f'Slice {Sind}'
                    zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    PlotTitle = PlotTitle.replace(' ', '\:')
                    plotTitle = r"$\bf{{{x}}}$".format(x=PlotTitle) \
                                + f'\n\n{SliceTxt}\n{FrameTxt}\n{zPosTxt}'
                    #plotTitle = r"$\bf{PlotTitle}$\," \
                    #            + f'\n\n{SliceTxt}\n{FrameTxt}\n{zPosTxt}'
                    
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    ax.set_title(plotTitle)
                    
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        import time
        #import os
        from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
        
    return




def PlotResultOfMaskToContourConversion(ListOfIms, 
                                        ListOfPtsByCntByRoi, ListOfC2SindsByRoi,
                                        ListOfPixArrBySeg, ListOfF2SindsBySeg,
                                        ListOfDicomDirs, ListOfPlotTitles, 
                                        AddTxt=None, ExportPlot=False, 
                                        ExportDir=None, TxtToAddToFname='', 
                                        LogToConsole=False):
    """ 
    08/06/21: Modelled on PlotPixArrsFromListOfLabImBySeg and 
    PlotContoursFromListOfRtss_v4.
    
    ListOfIms is a list (for each dataset, which could be of length 1) of
    SimpleITK Images.
    
    ListOfPtsByCntByRoi is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI) of a list (for each contour) of points.
    
    ListOfC2SindsByRoi is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI) of a list of the contour-to-slice indices in the
    contour data.
    
    ListOfPixArrBySeg is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI) of a Numpy pixel array, converted from the
    contour data.
    
    ListOfF2SindsBySeg is a list (for each dataset, which could be of length 1) 
    of a list (for each ROI) of a list of the frame-to-slice indices in the
    pixel arrays.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each image.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    
    
    Plot has different DATASETS ALONG COLUMNS and different SLICES ALONG ROWS. 
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import matplotlib.pyplot as plt
    import os
    #import time
    #from copy import deepcopy
    from pydicom import dcmread
    #from DicomTools import GetDicomFpaths
    from GeneralTools import UniqueItems, Unpack
    from ImageTools import GetImAttributesFromListOfDcmDirs, ImportImage
    from ConversionTools import Points2Indices
    
    ListOfDcmFpaths, ListOfIPPs, ListOfDirections, ListOfSpacings\
        = GetImAttributesFromListOfDcmDirs(ListOfDicomDirs)
    
    print(f'len(ListOfIms) = {len(ListOfIms)}')
    print(f'len(ListOfPtsByCntByRoi) = {len(ListOfPtsByCntByRoi)}')
    print(f'len(ListOfC2SindsByRoi) = {len(ListOfC2SindsByRoi)}')
    print(f'len(ListOfPixArrBySeg) = {len(ListOfPixArrBySeg)}')
    print(f'len(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
    print(f'len(ListOfDcmFpaths) = {len(ListOfDcmFpaths)}')
    print(f'ListOfC2SindsByRoi = {ListOfC2SindsByRoi}')
    print(f'ListOfF2SindsBySeg = {ListOfF2SindsBySeg}')
    
    
    # Get the maximum number of slices containing segmentations in any ROI in 
    # any dataset:
    MaxNumSlices = 0
    
    # The list (for each dataset) of the unique set of slice numbers that 
    # correspond to all segmentations in all ROIs:
    ListOfUniqueSinds = []
    
    # Loop through each dataset:
    for i in range(len(ListOfIms)):
        C2SindsByRoi = ListOfC2SindsByRoi[i]
        F2SindsBySeg = ListOfF2SindsBySeg[i]
        
        UniqueC2Sinds = UniqueItems(C2SindsByRoi)
        UniqueF2Sinds = UniqueItems(F2SindsBySeg)
        
        #UniqueSinds = UniqueItems(F2SindsBySeg)
        UniqueSinds = UniqueItems([UniqueC2Sinds, UniqueF2Sinds])
        
        ListOfUniqueSinds.append(UniqueSinds)
        
        if len(UniqueSinds) > MaxNumSlices:
            MaxNumSlices = len(UniqueSinds)
    
    
    print(f'UniqueC2Sinds = {UniqueC2Sinds}')
    print(f'UniqueF2Sinds = {UniqueF2Sinds}')
    print(f'ListOfUniqueSinds = {ListOfUniqueSinds}')
    print(f'MaxNumSlices = {MaxNumSlices}')
    
    
    """ Prepare the figure. """
    
    # Set the number of subplot rows and columns:
    Ncols = len(ListOfIms)
    Nrows = MaxNumSlices
    
    Cmaps = ['Reds', 'Blues', 'Greens', 'Oranges', 'Purples']
    
    LineWidth = 0.5
    LineWidth = 2
    Colours = ['r', 'b', 'm', 'c', 'k', 'y']
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    DcmAlpha = 0.2
    DcmAlpha = 0.5
    #DcmAlpha = 1
    
    SegAlpha = 0.5
    #SegAlpha = 0.2
    
    if ExportPlot:
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
    for rowNum in range(MaxNumSlices):
        # Loop through each data set:
        for i in range(Ncols):
            PtsByCntByRoi = ListOfPtsByCntByRoi[i]
            C2SindsByRoi = ListOfC2SindsByRoi[i]
            PixArrBySeg = ListOfPixArrBySeg[i]
            F2SindsBySeg = ListOfF2SindsBySeg[i]
            UniqueSinds = ListOfUniqueSinds[i]
            
            DicomFpaths = ListOfDcmFpaths[i]
            IPPs = ListOfIPPs[i]
            #Origin = ListOfOrigins[i]
            #Directions = ListOfDirections[i]
            #Spacings = ListOfSpacings[i]
            PlotTitle = ListOfPlotTitles[i]
            DicomDir = ListOfDicomDirs[i]
            
            RefIm = ImportImage(DicomDir)
            
            if LogToConsole:
                print(f'\n   i = {i}')
                print(f'   PlotTitle = {PlotTitle}')
                print(f'   F2SindsBySeg = {F2SindsBySeg}')
                print(f'   UniqueSinds = {UniqueSinds}')
            
            # Continue if rowNum does not exceed the number of unique slice indices:
            if rowNum < len(UniqueSinds):
                Sind = UniqueSinds[rowNum]
                
                ax = plt.subplot(Nrows, Ncols, n)
                
                DicomPixArr = dcmread(DicomFpaths[Sind]).pixel_array
            
                #im = ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r)
                ax.imshow(DicomPixArr, cmap=plt.cm.Greys_r, alpha=DcmAlpha)
                
                IPP = IPPs[Sind]
                
                if LogToConsole:
                    print(f'   Sind = {Sind}')
                    print(f'   IPP = {IPP}')
                
                                
                # Loop through each segment:
                for s in range(len(F2SindsBySeg)):
                    F2Sinds = F2SindsBySeg[s]
                    PixArr = PixArrBySeg[s]
                    
                    """ There are only len(Cmaps) colormaps defined above. Wrap the 
                    segment index s if there are more segments than the number of 
                    defined colormaps. """
                    if s < len(Cmaps):
                        Cmap = Cmaps[s]
                    else:
                        m = s//len(Cmaps)
                        
                        Cmap = Cmaps[s - m*len(Cmaps)]
                    
                        
                    if Sind in F2Sinds:
                        FrameNums = [i for i, e in enumerate(F2Sinds) if e==Sind]
                        
                        # Loop through all frame numbers:
                        for f in range(len(FrameNums)):
                            FrameNum = FrameNums[f]
                        
                            Frame = PixArr[FrameNum]
                            
                            if LogToConsole:
                                print(f'      SliceNum = {Sind} is in F2Sinds')
                                print(f'      FrameNum = {FrameNum}')
                                print(f'      Frame.shape = {Frame.shape}')
                        
                            #ax = plt.subplot(Nrows, Ncols, n)
                            #im = ax.imshow(Frame, cmap=plt.cm.nipy_spectral, 
                            #               alpha=alpha)
                            ax.imshow(Frame, cmap=Cmap, alpha=SegAlpha)
                            
                            #FrameTxt = f'Frame {FrameNum}'
                        #if len(FrameNums) > 1:
                        #    FrameTxt = f'Conversion to segmentations {FrameNums}'
                        #else:
                        #    FrameTxt = f'Conversion to segmentations {FrameNum}'
                        #FrameTxt = f'Contour-to-segmentation # {FrameNums}'
                        nobrackets = f'{FrameNums}'.replace('[', '').replace(']', '')
                        if i == 0:
                            FrameTxt = f'Contour-to-segmentation # {nobrackets}'
                        else:
                            FrameTxt = f'Segmentation # {nobrackets}'
                    else:
                        FrameTxt = 'No Segmentations'
                        
                        if LogToConsole:
                            print(f'      Sind = {Sind} is NOT in F2Sinds')
                        
                    #SliceTxt = f'Slice {Sind}'
                    #zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    #PlotTitle = PlotTitle.replace(' ', '\:')
                    #plotTitle = r"$\bf{{{x}}}$".format(x=PlotTitle) \
                    #            + f'\n\n{SliceTxt}\n{FrameTxt}\n{zPosTxt}'
                    ##plotTitle = r"$\bf{PlotTitle}$\," \
                    ##            + f'\n\n{SliceTxt}\n{FrameTxt}\n{zPosTxt}'
                    #
                    #ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    #ax.set_title(plotTitle)
                    
                    
                # Loop through each ROI:
                for r in range(len(C2SindsByRoi)):
                    C2Sinds = C2SindsByRoi[r]
                    PtsByCnt = PtsByCntByRoi[r]
                    
                    """ There are only len(Colours) colours defined above. Wrap the 
                    ROI index r if there are more ROIs than the number of defined 
                    colours. """
                    if r < len(Colours):
                        colour = Colours[r]
                    else:
                        m = r//len(Colours)
                        
                        colour = Colours[r - m*len(Colours)]
                
                    if Sind in C2Sinds:
                        #ContourNum = C2Sinds.index(SliceNum) # there may be more 
                        # than one contour for SliceNum!!
                        ContourNums = [i for i, e in enumerate(C2Sinds) if e==Sind]
                        
                        #print(f'ContourNums = {ContourNums}')
                        
                        NumPtsByCnt = []
                        
                        # Loop through all contour numbers:
                        for c in range(len(ContourNums)):
                            colour2 = Colours[c] 
                            
                            ContourNum = ContourNums[c]
                            
                            Pts = PtsByCnt[ContourNum]
                            
                            NumPtsByCnt.append(len(Pts))
                            
                            #print(f'\nPts = {Pts}')
                            
                            Inds = Points2Indices(Points=Pts, RefIm=RefIm)
                            
                            if LogToConsole:
                                print(f'      Sind = {Sind} is in C2Sinds')
                                print(f'      ContourNum = {ContourNum}')
                                #print(f'      Contour = {Contour}')
                                print(f'      len(Pts) = {len(Pts)}')
                            
                            #ax = plt.subplot(Nrows, Ncols, n)
                            
                            # Unpack the contour points' indices:
                            X, Y, Z = Unpack(Inds)
                            
                            # Plot the contour points:
                            #ax = plt.plot(X, Y, linewidth=0.5, c='r')
                            #ax.plot(X, Y, linewidth=0.5, c='r')
                            #ax.plot(X, Y, linewidth=LineWidth, c=colour)
                            ax.plot(X, Y, linewidth=LineWidth, c=colour2)
                        
                        #if len(ContourNums) > 1:
                        #    ContourTxt = f'Contours {ContourNums}'
                        #else:
                        #    ContourTxt = f'Contour {ContourNum}'
                        #ContourTxt = f'Contour # {ContourNums}'
                        nobrackets = f'{ContourNums}'.replace('[', '').replace(']', '')
                        if i == 0:
                            ContourTxt = f'Contour # {nobrackets}'
                        else:
                            ContourTxt = f'Segmentation-to-contour # {nobrackets}'
                        nobrackets = f'{NumPtsByCnt}'.replace('[', '').replace(']', '')
                        PtsTxt = f'No. of points: {nobrackets}'
                    else:
                        ContourTxt = 'No contour'
                        PtsTxt = 'No. of points: 0'
                        
                        if LogToConsole:
                            print(f'      Sind = {Sind} is NOT in C2Sinds')
                    
                    SliceTxt = f'Slice # {Sind}'
                    zPosTxt = f'z = {round(IPP[2], 2)} mm'
                    if rowNum == 0:
                        PlotTitle = PlotTitle.replace(' ', '\:')
                        plotTitle = r"$\bf{{{x}}}$".format(x=PlotTitle)
                    else:
                        plotTitle = ''
                    plotTitle += f'\n\n{SliceTxt}'
                    if i == 0:
                        plotTitle += f'\n{ContourTxt}\n{PtsTxt}\n{FrameTxt}'
                    else:
                        plotTitle += f'\n{FrameTxt}\n{ContourTxt}\n{PtsTxt}'
                    plotTitle += f'\n{zPosTxt}'
                    
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    ax.set_title(plotTitle)
                    
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        import time
        #import os
        from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
        
    return





def PlotSitkImage(Image, FrameNum=None, Colormap='grayscale', 
                  ExportPlot=False, ExportDir='cwd', TxtToAddToFname='',
                  LogToConsole=False):
    """ 
     
    """
    
    import SimpleITK as sitk
    import numpy as np
    import matplotlib.pyplot as plt
    #from ConversionTools import Image2PixArr
    
    #PixArr, F2Sinds = Image2PixArr(Image)
    
    PixArr = sitk.GetArrayFromImage(Image)
    
    if FrameNum != None:
        PixArr = PixArr[FrameNum]
        
        R, C = PixArr.shape
        
        PixArr = np.reshape(PixArr, (1, R, C))

    F, R, C = PixArr.shape
    
    #print(f'\nPixArr.shape = {PixArr.shape}')
    #print(f'\nF2Sinds = {F2Sinds}')
    
    if Colormap == 'color':
        cmap = plt.cm.nipy_spectral
    else:
        cmap = plt.cm.Greys_r
    
    if ExportPlot:
        dpi = 300
    else:
        dpi = 80
    
    # Set the number of subplot rows and columns:
    Ncols = 1
    Nrows = F
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    n = 1 # initialised sub-plot number
    
    if FrameNum == None:
        for f in range(F):
            ax = plt.subplot(Nrows, Ncols, n)
            
            ax.imshow(PixArr[f], cmap=cmap)
                    
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            #ax.set_title(f'Frame {f} = Slice {F2Sinds[f]}')
            ax.set_title(f'Frame {f}')
                    
            n += 1 # increment sub-plot number
    else:
        ax = plt.subplot(Nrows, Ncols, n)
            
        ax.imshow(PixArr[0], cmap=cmap)
                
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        #ax.set_title(f'Frame {f} = Slice {F2Sinds[f]}')
        ax.set_title(f'Frame {FrameNum}')
    
    
    if ExportPlot:
        import time
        import os
        from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
    return




def PlotImageAndLabmapIm(DicomIm, LabmapIm, dpi=80, LogToConsole=False):
    """ 
     
    """
    
    import matplotlib.pyplot as plt
    from ConversionTools import Image2PixArr
    
    DcmPixArr, DcmF2Sinds = Image2PixArr(DicomIm)
    DcmF, DcmR, DcmC = DcmPixArr.shape
    
    LmapPixArr, LmapF2Sinds = Image2PixArr(LabmapIm)
    LmapF, LmapR, LmapC = LmapPixArr.shape
    
    print(f'\nDcmPixArr.shape = {DcmPixArr.shape}')
    print(f'DcmF2Sinds = {DcmF2Sinds}')
    
    print(f'\nLmapPixArr.shape = {LmapPixArr.shape}')
    print(f'LmapF2Sinds = {LmapF2Sinds}')
    
    #if Colormap == 'color':
    #    cmap = plt.cm.nipy_spectral
    #else:
    #    cmap = plt.cm.Greys_r
    
    # Colormap for labelmaps:
    Cmap = plt.cm.nipy_spectral
    #Cmap = plt.cm.hsv
    #Cmap = plt.cm.gist_rainbow
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    #alpha = 0.2
    alpha = 0.5
    
    # Set the number of subplot rows and columns:
    Ncols = 1
    Nrows = LmapF
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    n = 1 # initialised sub-plot number
    
    for f in range(LmapF):
        ax = plt.subplot(Nrows, Ncols, n)
        
        # The DICOM slice number:
        s = LmapF2Sinds[f]
        
        #ax.imshow(DcmPixArr[s], cmap=plt.cm.Greys_r)
        ax.imshow(DcmPixArr[s], cmap=plt.cm.Greys_r, alpha=alpha)
        
        ax.imshow(LmapPixArr[f], cmap=Cmap, alpha=alpha)
                
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        ax.set_title(f'Frame = {f}, Slice {s}')
                
        n += 1 # increment sub-plot number
        
    return




def PlotImageAndLabmapImBySeg_NOT_COMPLETE(DicomIm, LabmapImBySeg, dpi=80, LogToConsole=False):
    """ 
     
    """
    
    import matplotlib.pyplot as plt
    from ConversionTools import Image2PixArr
    
    DcmPixArr, DcmF2Sinds = Image2PixArr(DicomIm)
    DcmF, DcmR, DcmC = DcmPixArr.shape
    
    LmapPixArr, LmapF2Sinds = Image2PixArr(LabmapIm)
    LmapF, LmapR, LmapC = LmapPixArr.shape
    
    print(f'\nDcmPixArr.shape = {DcmPixArr.shape}')
    print(f'DcmF2Sinds = {DcmF2Sinds}')
    
    print(f'\nLmapPixArr.shape = {LmapPixArr.shape}')
    print(f'LmapF2Sinds = {LmapF2Sinds}')
    
    #if Colormap == 'color':
    #    cmap = plt.cm.nipy_spectral
    #else:
    #    cmap = plt.cm.Greys_r
    
    # Colormap for labelmaps:
    Cmap = plt.cm.nipy_spectral
    #Cmap = plt.cm.hsv
    #Cmap = plt.cm.gist_rainbow
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    #alpha = 0.2
    alpha = 0.5
    
    # Set the number of subplot rows and columns:
    Ncols = 1
    Nrows = LmapF
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    n = 1 # initialised sub-plot number
    
    for f in range(LmapF):
        ax = plt.subplot(Nrows, Ncols, n)
        
        # The DICOM slice number:
        s = LmapF2Sinds[f]
        
        #ax.imshow(DcmPixArr[s], cmap=plt.cm.Greys_r)
        ax.imshow(DcmPixArr[s], cmap=plt.cm.Greys_r, alpha=alpha)
        
        ax.imshow(LmapPixArr[f], cmap=Cmap, alpha=alpha)
                
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        ax.set_title(f'Frame = {f}, Slice {s}')
                
        n += 1 # increment sub-plot number
        
    return





def PlotImageAndContours(DicomIm, PtsByCnt, C2Sinds, dpi=80, 
                         LogToConsole=False):
    """ 
     
    """
    
    import matplotlib.pyplot as plt
    from ConversionTools import Image2PixArr, Points2Indices
    from GeneralTools import Unpack
    
    DcmPixArr, DcmF2Sinds = Image2PixArr(DicomIm)
    DcmF, DcmR, DcmC = DcmPixArr.shape
    
    print(f'\nDcmPixArr.shape = {DcmPixArr.shape}')
    print(f'DcmF2Sinds = {DcmF2Sinds}')
    
    C = len(C2Sinds)
    print(f'\nC2Sinds = {C2Sinds}')
    
    #if Colormap == 'color':
    #    cmap = plt.cm.nipy_spectral
    #else:
    #    cmap = plt.cm.Greys_r
    
    # Colormap for labelmaps:
    #Cmap = plt.cm.nipy_spectral
    #Cmap = plt.cm.hsv
    #Cmap = plt.cm.gist_rainbow
    
    # Set the transparency of labelmaps to be overlaid over DICOMs:
    #alpha = 0.2
    alpha = 0.5
    alpha = 1
    
    LineWidth = 0.5
    LineWidth = 1
    LineWidth = 2
    #Colors = ['r', 'b', 'm', 'c', 'k', 'y']
    
    # Set the number of subplot rows and columns:
    Ncols = 1
    Nrows = C
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    n = 1 # initialised sub-plot number
    
    for c in range(C):
        ax = plt.subplot(Nrows, Ncols, n)
        
        # The DICOM slice number:
        s = C2Sinds[c]
        
        #ax.imshow(DcmPixArr[s], cmap=plt.cm.Greys_r)
        ax.imshow(DcmPixArr[s], cmap=plt.cm.Greys_r, alpha=alpha)
        
        Pts = PtsByCnt[c]
        
        Inds = Points2Indices(Points=Pts, RefIm=DicomIm)
        
        # Unpack the contour points' indices:
        X, Y, Z = Unpack(Inds)
        
        ax.plot(X, Y, linewidth=LineWidth, c='r')
                        
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        ax.set_title(f'Contour = {c}, Slice = {s}')
                
        n += 1 # increment sub-plot number
        
    return






def PlotSitkLabelmap(Labelmap, Colormap='grayscale', dpi=80, LogToConsole=False):
    """ 
     
    """
    
    import importlib
    import ConversionTools
    importlib.reload(ConversionTools)
    
    import matplotlib.pyplot as plt
    from ConversionTools import Labelmap2BinaryIm, Image2PixArr
    
    Image = Labelmap2BinaryIm(Labelmap)
    
    PixArr, F2Sinds = Image2PixArr(Image)
    
    F, R, C = PixArr.shape
    
    print(f'\nPixArr.shape = {PixArr.shape}')
    print(f'\nF2Sinds = {F2Sinds}')
    
    if Colormap == 'color':
        cmap = plt.cm.nipy_spectral
    else:
        cmap = plt.cm.Greys_r
    
    # Set the number of subplot rows and columns:
    Ncols = 1
    Nrows = F
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    n = 1 # initialised sub-plot number
    
    for f in range(F):
        ax = plt.subplot(Nrows, Ncols, n)
        
        ax.imshow(PixArr[f], cmap=cmap)
                
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        ax.set_title(f'Frame {f} = Slice {F2Sinds[f]}')
                
        n += 1 # increment sub-plot number
        
    return





def CompareSitkFrames(Image0, Ind0, Image1, Ind1, PlotTxt0='', PlotTxt1='',
                      Colormap='grayscale', dpi=80, ExportPlot=False, 
                      ExportDir=None, TxtToAddToFname='', LogToConsole=False):
    """ 
    Compare two frames from two sitk images at specified slice indices.
    """
    
    import os
    import time
    import matplotlib.pyplot as plt
    from ConversionTools import Image2PixArr
    
    if Image0.GetSize()[1] != Image1.GetSize()[1] or Image0.GetSize()[2] != Image1.GetSize()[2]:
        from ImageTools import ResampleImage
        
        Image0 = ResampleImage(Image0, Image1, 'linear')
        
        """ Since Image0 has been resampled to Image1, use the same Ind0 as 
        Ind1: """
        Ind0 = Ind1
        
        PlotTxt0 = 'Resampled ' + PlotTxt0
    
    PixArr0, F2Sinds0 = Image2PixArr(Image0)
    PixArr1, F2Sinds1 = Image2PixArr(Image1)
    
    Frame0 = PixArr0[Ind0]
    Frame1 = PixArr1[Ind1]
    
    Diff = Frame1 - Frame0
    
    z0 = Image0.TransformIndexToPhysicalPoint((0,0,Ind0))[2]
    z1 = Image1.TransformIndexToPhysicalPoint((0,0,Ind1))[2]
    
    dz = z1 - z0
    
    PlotTitle0 = PlotTxt0 + f'\nFrame {Ind0} \nz = {round(z0, 2)} mm'
    PlotTitle1 = PlotTxt1 + f'\nFrame {Ind1} \nz = {round(z1, 2)} mm'
    PlotTitle2 = f'Difference \ndz = {round(dz, 2)} mm' 
    
    if Colormap == 'color':
        cmap = plt.cm.nipy_spectral
    else:
        cmap = plt.cm.Greys_r
        
    if ExportPlot:
        dpi = 120
    else:
        dpi = 80
    
    # Set the number of subplot rows and columns:
    Ncols = 3
    Nrows = 1
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 8*Nrows), dpi=dpi)
        
    ax = plt.subplot(Nrows, Ncols, 1)
    ax.imshow(Frame0, cmap=cmap)
    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
    ax.set_title(PlotTitle0)
    
    ax = plt.subplot(Nrows, Ncols, 2)
    ax.imshow(Frame1, cmap=cmap)
    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
    ax.set_title(PlotTitle1)
    
    ax = plt.subplot(Nrows, Ncols, 3)
    ax.imshow(Diff, cmap=cmap)
    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
    ax.set_title(PlotTitle2)
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname \
                      + '_ComparisonOfFrames.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return






def PlotPoints(Points, dim=3):
    from GeneralTools import Unpack
    import matplotlib.pyplot as plt
    
    # Unpack Points:
    X, Y, Z = Unpack(Points)
        
    fig = plt.figure(figsize=(13, 13))
    
    if dim==3:
        from mpl_toolkits.mplot3d import Axes3D
        
        ax = fig.add_subplot(111, projection='3d')
        
        ax.scatter(X, Y, Z, marker='o')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
    
    elif dim==2:
        ax = fig.add_subplot(111)
        
        ax.scatter(X, Y, marker='o')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    else:
        msg = 'dim must be 2 or 3'
        
        raise Exception(msg)
    
    return




def PlotTwoSetsOfPoints(Points0, Points1, dim=3):
    from GeneralTools import Unpack
    import matplotlib.pyplot as plt
    
    # Unpack points:
    X0, Y0, Z0 = Unpack(Points0)
    X1, Y1, Z1 = Unpack(Points1)
    
    fig = plt.figure(figsize=(13, 13))
    
    if dim==3:
        from mpl_toolkits.mplot3d import Axes3D
        
        ax = fig.add_subplot(111, projection='3d')
        
        ax.scatter(X0, Y0, Z0, marker='o', color='b')
        ax.scatter(X1, Y1, Z1, marker='o', color='r')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
    
    elif dim==2:
        ax = fig.add_subplot(111)
        
        ax.scatter(X0, Y0, marker='o', color='b')
        ax.scatter(X1, Y1, marker='o', color='r')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    else:
        msg = 'dim must be 2 or 3'
        
        raise Exception(msg)
    
    return


def PlotPixArrBySeg(PixArrBySeg, F2SindsBySeg, PlotTitle=''):
    """ 
    09/07/21:
        See plot_pixarrBySeg in plotting_tools.plotting.py
    
    PixArrBySeg can either be a list of pixel arrays (as the variable name
    suggests), or a pixel array, and F2SindsBySeg can either be a list (for 
    each segment) of a list (for each frame) of frame-to-slice indices, or it
    can be a list (for each frame) of frame-to-slice indices.
    
    The pixel array(s) must have shape (F, R, C) where F is the number of
    frames, which could be any integer number including 0. """
    
    import matplotlib.pyplot as plt
    
    # Prepare the figure.
    
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





def PlotLabmapImBySeg(LabmapImBySeg, F2SindsBySeg, PlotTitle=''):
    """ LabmapImBySeg can either be a list of labelmap images (as the variable 
    name suggests), or a single labelmap image, and F2SindsBySeg can either be 
    a list (for each segment) of a list (for each frame) of frame-to-slice 
    indices, or it can be a list (for each frame) of frame-to-slice indices.
    """
    
    import importlib
    import ConversionTools
    importlib.reload(ConversionTools)
    
    from ConversionTools import ImageBySeg2PixArrBySeg
    
    #PixArrBySeg, F2SindsBySeg = ImageBySeg2PixArrBySeg(LabmapImBySeg, 
    #                                                   F2SindsBySeg)
    PixArrBySeg, F2SindsBySeg = ImageBySeg2PixArrBySeg(LabmapImBySeg)
    
    PlotPixArrBySeg(PixArrBySeg, F2SindsBySeg, PlotTitle)
        
    return






def overlay(img1, img2, title=None, interpolation=None, sizeThreshold=128):
    """
    Description
    -----------
    Displays an overlay of two 2D images using matplotlib.pyplot.imshow().
    
    Source:
    ******
    https://github.com/SuperElastix/SimpleElastix/issues/79

    Args
    ----
    img1 : np.array or sitk.Image
        image to be displayed
    img2 : np.array or sitk.Image
        image to be displayed

    Optional
    --------
    title : string
        string used as image title
        (default None)
    interpolation : string
        desired option for interpolation among matplotlib.pyplot.imshow options
        (default nearest)
    sizeThreshold : integer
        size under which interpolation is automatically set to 'nearest' if all
        dimensions of img1 and img2 are below
        (default 128)
    """
    # Check for type of images and convert to np.array
    if isinstance(img1, sitk.Image):
        img1 = sitk.GetArrayFromImage(img1)
    if isinstance(img2, sitk.Image):
        img2 = sitk.GetArrayFromImage(img2)
    if type(img1) is not type(img2) is not np.ndarray:
        raise NotImplementedError('Please provide images as np.array or '
                                  'sitk.Image.')
    # Check for size of images
    if not img1.ndim == img2.ndim == 2:
        raise NotImplementedError('Only supports 2D images.')

    if interpolation:
        plt.imshow(img1, cmap='summer', interpolation=interpolation)
        plt.imshow(img2, cmap='autumn', alpha=0.5, interpolation=interpolation)
    elif max(max(img1.shape), max(img2.shape)) > sizeThreshold:
        plt.imshow(img1, cmap='summer')
        plt.imshow(img2, cmap='autumn', alpha=0.5)
    else:
        plt.imshow(img1, cmap='summer', interpolation='nearest')
        plt.imshow(img2, cmap='autumn', alpha=0.5, interpolation='nearest')
    plt.title(title)
    plt.axis('off')
    plt.show()
    
    return




""" The next few functions from here:
https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
"""

def PlotFixMovImages(FixInd, MovInd, FixPixArr, MovPixArr,
                     FixTitle='Fixed image', MovTitle='Moving image'):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    import matplotlib.pyplot as plt
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(FixPixArr[FixInd,:,:], cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(MovPixArr[MovInd,:,:], cmap=plt.cm.Greys_r);
    plt.title(MovTitle)
    plt.axis('off')
    
    plt.show()
    
    return




def PlotFixMovImages_v2(FixIm, MovIm, FixInd, MovInd, 
                     FixTitle='Fixed image', MovTitle='Moving image'):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
    FixPixArr = sitk.GetArrayViewFromImage(FixIm)[FixInd,:,:]
    MovPixArr = sitk.GetArrayViewFromImage(MovIm)[MovInd,:,:]
    
    FixZ = FixIm.TransformIndexToPhysicalPoint([0,0,FixInd])[2]
    MovZ = MovIm.TransformIndexToPhysicalPoint([0,0,MovInd])[2]
    
    FixTitle = f'{FixTitle}\nz = {round(FixZ, 2)} mm'
    MovTitle = f'{MovTitle}\nz = {round(MovZ, 2)} mm'
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(FixPixArr, cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(MovPixArr, cmap=plt.cm.Greys_r);
    plt.title(MovTitle)
    plt.axis('off')
    
    plt.show()
    
    return





def PlotFixMovImagesWithMouseInput(FixInd, MovInd, FixPixArr, MovPixArr,
                                   FixTitle='Fixed image', MovTitle='Moving image'):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    import matplotlib.pyplot as plt
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(FixPixArr[FixInd,:,:], cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    FixInput = plt.ginput(0,0)
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(MovPixArr[MovInd,:,:], cmap=plt.cm.Greys_r);
    plt.title(MovTitle)
    plt.axis('off')
    
    
    return FixInput





def PlotFixMovImagesWithMouseInput_v2(FixIm, MovIm, FixInd, MovInd, 
                     FixTitle='Fixed image', MovTitle='Moving image'):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
    FixPixArr = sitk.GetArrayViewFromImage(FixIm)[FixInd,:,:]
    MovPixArr = sitk.GetArrayViewFromImage(MovIm)[MovInd,:,:]
    
    FixZ = FixIm.TransformIndexToPhysicalPoint([0,0,FixInd])[2]
    MovZ = MovIm.TransformIndexToPhysicalPoint([0,0,MovInd])[2]
    
    FixTitle = f'{FixTitle}\nz = {round(FixZ, 2)} mm'
    MovTitle = f'{MovTitle}\nz = {round(MovZ, 2)} mm'
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(FixPixArr, cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    FixInput = plt.ginput(0,0)
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(MovPixArr, cmap=plt.cm.Greys_r);
    plt.title(MovTitle)
    plt.axis('off')
    MovInput = plt.ginput(0,0)
    
    plt.show()
    
    return FixInput, MovInput






def PlotBlendedImage(FixIm, ResIm, Ind, alpha=0.5, PlotTitle='Blended image'):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
    plt.subplots(1, 1, figsize=(6,5))
    
    BlendedIm = (1.0 - alpha)*FixIm[:,:,Ind] + alpha*ResIm[:,:,Ind] 
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.title(PlotTitle)
    plt.axis('off')
    plt.show()
    
    return


def StartPlot():
    """ Callback invoked when the StartEvent happens, sets up our new data. """
    global metric_values, multires_iterations
    
    metric_values = []
    multires_iterations = []
    
    return


def EndPlot():
    """ Callback invoked when the EndEvent happens, do cleanup of data and 
    figure. """
    
    global metric_values, multires_iterations
    import matplotlib.pyplot as plt
    
    del metric_values
    del multires_iterations
    # Close figure, we don't want to get a duplicate of the plot latter on.
    plt.close()
    
    return


def PlotValues(RegMethod):
    """ Callback invoked when the IterationEvent happens, update our data 
    and display new figure. """
    
    global metric_values, multires_iterations
    import matplotlib.pyplot as plt
    from IPython.display import clear_output
    
    metric_values.append(RegMethod.GetMetricValue())                                       
    # Clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True)
    # Plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
    plt.xlabel('Iteration Number',fontsize=12)
    plt.ylabel('Metric Value',fontsize=12)
    plt.show()
    
    return






def PlotSrcResSrcTrgImages(SrcIm, ResSrcIm, TrgIm, SrcInd, ResSrcInd, TrgInd,
                           SrcTitle='Source image', 
                           ResSrcTitle='Resampled Source image', 
                           TrgTitle='Target image'):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
    SrcPixArr = sitk.GetArrayViewFromImage(SrcIm)[SrcInd,:,:]
    ResSrcPixArr = sitk.GetArrayViewFromImage(ResSrcIm)[ResSrcInd,:,:]
    TrgPixArr = sitk.GetArrayViewFromImage(TrgIm)[TrgInd,:,:]
    
    SrcZ = SrcIm.TransformIndexToPhysicalPoint([0,0,SrcInd])[2]
    ResSrcZ = ResSrcIm.TransformIndexToPhysicalPoint([0,0,ResSrcInd])[2]
    TrgZ = SrcIm.TransformIndexToPhysicalPoint([0,0,TrgInd])[2]
    
    SrcTitle = f'{SrcTitle}\nz = {round(SrcZ, 2)} mm'
    ResSrcTitle = f'{ResSrcTitle}\nz = {round(ResSrcZ, 2)} mm'
    TrgTitle = f'{TrgTitle}\nz = {round(TrgZ, 2)} mm'
    
    plt.subplots(1, 3, figsize=(15,8))
    
    # Draw the Source/Moving image in the first subplot.
    plt.subplot(1, 3, 1)
    plt.imshow(SrcPixArr, cmap=plt.cm.Greys_r);
    plt.title(SrcTitle)
    plt.axis('off')
    
    # Draw the Resampled/registered Source/Moving image in the second subplot.
    plt.subplot(1, 3, 2)
    plt.imshow(ResSrcPixArr, cmap=plt.cm.Greys_r);
    plt.title(ResSrcTitle)
    plt.axis('off')
    
    # Draw the Target/Fixed image in the third subplot.
    plt.subplot(1, 3, 3)
    plt.imshow(TrgPixArr, cmap=plt.cm.Greys_r);
    plt.title(TrgTitle)
    plt.axis('off')
    
    plt.show()
    
    return



"""
R0, C0, S0 = Im0.GetSize()
R1, C1, S1 = Im1.GetSize()
R2, C2, S2 = Im2.GetSize()
"""


def PlotTwoImages(Im0, Ind0, PlotLabel0, Im1, Ind1, PlotLabel1):
    """
    08/07/21:
        See plot_two_ims in plotting_tools.plotting.py
    """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
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






def PlotThreeImages(Im0, Ind0, PlotLabel0, 
                    Im1, Ind1, PlotLabel1, 
                    Im2, Ind2, PlotLabel2):
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
    PixArr0 = sitk.GetArrayViewFromImage(Im0)[Ind0,:,:]
    PixArr1 = sitk.GetArrayViewFromImage(Im1)[Ind1,:,:]
    PixArr2 = sitk.GetArrayViewFromImage(Im2)[Ind2,:,:]
    
    z0 = Im0.TransformIndexToPhysicalPoint([0,0,Ind0])[2]
    z1 = Im1.TransformIndexToPhysicalPoint([0,0,Ind1])[2]
    z2 = Im2.TransformIndexToPhysicalPoint([0,0,Ind2])[2]
    
    PlotLabel0 = f'{PlotLabel0}\nz = {round(z0, 2)} mm'
    PlotLabel1 = f'{PlotLabel1}\nz = {round(z1, 2)} mm'
    PlotLabel2 = f'{PlotLabel2}\nz = {round(z2, 2)} mm'
    
    plt.subplots(1, 3, figsize=(15,8))
    
    plt.subplot(1, 3, 1)
    plt.imshow(PixArr0, cmap=plt.cm.Greys_r);
    plt.title(PlotLabel0)
    plt.axis('off')
    
    plt.subplot(1, 3, 2)
    plt.imshow(PixArr1, cmap=plt.cm.Greys_r);
    plt.title(PlotLabel1)
    plt.axis('off')
    
    plt.subplot(1, 3, 3)
    plt.imshow(PixArr2, cmap=plt.cm.Greys_r);
    plt.title(PlotLabel2)
    plt.axis('off')
    
    plt.show()
    
    return






def PlotTwoImagesWithFiducials(Im0, FidsFpath0, PlotLabel0, 
                               Im1, FidsFpath1, PlotLabel1,
                               IndsFromImageJ=False):
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    import importlib
    import ImageTools
    importlib.reload(ImageTools)
    from ImageTools import ImportFiducialsFromTxt
    
    Fids0, FidsType0 = ImportFiducialsFromTxt(FidsFpath0, IndsFromImageJ, Im0)
    Fids1, FidsType1 = ImportFiducialsFromTxt(FidsFpath1, IndsFromImageJ, Im1)
    
    Nrows = min(len(Fids0), len(Fids1))
    
    plt.subplots(Nrows, 2, figsize=(15,7*Nrows))
    
    i = 1 # sub-plot number
    
    for r in range(Nrows):
        i0, j0, k0 = Fids0[r]
        i1, j1, k1 = Fids1[r]
        
        PixArr0 = sitk.GetArrayViewFromImage(Im0)[k0,:,:]
        PixArr1 = sitk.GetArrayViewFromImage(Im1)[k1,:,:]
        
        z0 = Im0.TransformIndexToPhysicalPoint([i0,j0,k0])[2]
        z1 = Im1.TransformIndexToPhysicalPoint([i1,j1,k1])[2]
        
        Lab0 = f'{PlotLabel0}\nz = {round(z0, 2)} mm\nFudicial at {Fids0[r]}'
        Lab1 = f'{PlotLabel1}\nz = {round(z1, 2)} mm\nFudicial at {Fids1[r]}'
    
        plt.subplot(Nrows, 2, i)
        plt.imshow(PixArr0, cmap=plt.cm.Greys_r);
        plt.scatter(i0, j0, marker='+', c='r');
        plt.title(Lab0)
        plt.axis('off')
        
        i += 1
        
        plt.subplot(Nrows, 2, i)
        plt.imshow(PixArr1, cmap=plt.cm.Greys_r);
        plt.scatter(i1, j1, marker='+', c='r');
        plt.title(Lab1)
        plt.axis('off')
        
        i += 1
    
    plt.show()
    
    return






def PlotThreeImagesWithFiducials(Im0, FidsFpath0, PlotLabel0, 
                                 Im1, FidsFpath1, PlotLabel1,
                                 Im2, FidsFpath2, PlotLabel2,
                                 IndsFromImageJ=False, 
                                 ExportPlot=False, 
                                 ExportDir='{cwd}/fiducials',
                                 TxtToAddToFname=''):
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    import importlib
    import ImageTools
    importlib.reload(ImageTools)
    from ImageTools import ImportFiducialsFromTxt
    
    if ExportPlot:
        dpi = 120
    else:
        dpi = 80
    
    Fids0, FidsType0 = ImportFiducialsFromTxt(FidsFpath0, IndsFromImageJ, Im0)
    Fids1, FidsType1 = ImportFiducialsFromTxt(FidsFpath1, IndsFromImageJ, Im1)
    Fids2, FidsType2 = ImportFiducialsFromTxt(FidsFpath2, IndsFromImageJ, Im2)
    
    Nrows = min(len(Fids0), len(Fids1), len(Fids2))
    
    plt.subplots(Nrows, 3, figsize=(20,7*Nrows), dpi=dpi)
    
    marker = '+'
    markersize = 154
    colour = 'r'
    
    i = 1 # sub-plot number
    
    for r in range(Nrows):
        i0, j0, k0 = Fids0[r]
        i1, j1, k1 = Fids1[r]
        i2, j2, k2 = Fids2[r]
        
        PixArr0 = sitk.GetArrayViewFromImage(Im0)[k0,:,:]
        PixArr1 = sitk.GetArrayViewFromImage(Im1)[k1,:,:]
        PixArr2 = sitk.GetArrayViewFromImage(Im2)[k2,:,:]
        
        z0 = Im0.TransformIndexToPhysicalPoint([i0,j0,k0])[2]
        z1 = Im1.TransformIndexToPhysicalPoint([i1,j1,k1])[2]
        z2 = Im2.TransformIndexToPhysicalPoint([i2,j2,k2])[2]
        
        Lab0 = f'{PlotLabel0}\nz = {round(z0, 2)} mm\nFudicial at {Fids0[r]}'
        Lab1 = f'{PlotLabel1}\nz = {round(z1, 2)} mm\nFudicial at {Fids1[r]}'
        Lab2 = f'{PlotLabel2}\nz = {round(z2, 2)} mm\nFudicial at {Fids2[r]}'
    
        plt.subplot(Nrows, 3, i)
        plt.imshow(PixArr0, cmap=plt.cm.Greys_r);
        plt.scatter(i0, j0, marker=marker, s=markersize, c=colour);
        plt.title(Lab0)
        plt.axis('off')
        
        i += 1
        
        plt.subplot(Nrows, 3, i)
        plt.imshow(PixArr1, cmap=plt.cm.Greys_r);
        plt.scatter(i1, j1, marker=marker, s=markersize, c=colour);
        plt.title(Lab1)
        plt.axis('off')
        
        i += 1
        
        plt.subplot(Nrows, 3, i)
        plt.imshow(PixArr2, cmap=plt.cm.Greys_r);
        plt.scatter(i2, j2, marker=marker, s=markersize, c=colour);
        plt.title(Lab2)
        plt.axis('off')
        
        i += 1
    
    #plt.show()
    
    if ExportPlot:
        import time
        import os
        from pathlib import Path
        
        if ExportDir == '{cwd}/fiducials':
            cwd = os.getcwd()
            
            ExportDir = os.path.join(cwd, 'fiducials')
            
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
    
    return







def PlotResResults(FixIm, MovIm, ResIm, FixInd, MovInd=None,
                   FixTitle='Fixed image', MovTitle='Moving image', 
                   ResTitle='Resampled image',
                   ExportPlot=False, ExportDir='cwd', TxtToAddToFname=''):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. 
    """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    from GeneralTools import GetIndOfNearestSlice
    
    #UseCentre = True # 16/04/21 Not sure it's working as expected
    UseCentre = False
    
    if MovInd == None:
        MovInd, zDiff = GetIndOfNearestSlice(FixIm, FixInd, MovIm, UseCentre)
    
    FixPixArr = sitk.GetArrayViewFromImage(FixIm)[FixInd,:,:]
    MovPixArr = sitk.GetArrayViewFromImage(MovIm)[MovInd,:,:]
    ResPixArr = sitk.GetArrayViewFromImage(ResIm)[FixInd,:,:]
    
    if UseCentre:
        R_f, C_f, S_f = FixIm.GetSize()
        R_m, C_m, S_m = MovIm.GetSize()
        
        i_f = R_f//2
        j_f = C_f//2
        
        i_m = R_m//2
        j_m = C_m//2
    else:
        i_f = 0
        j_f = 0
        
        i_m = 0
        j_m = 0
        
    FixZ = FixIm.TransformIndexToPhysicalPoint([i_f,j_f,FixInd])[2]
    MovZ = MovIm.TransformIndexToPhysicalPoint([i_m,j_m,MovInd])[2]
    ResZ = ResIm.TransformIndexToPhysicalPoint([i_f,j_f,FixInd])[2]
    
    FixTitle = f'{FixTitle}\nk = {FixInd}\nz = {round(FixZ, 2)} mm'
    MovTitle = f'{MovTitle}\nk = {MovInd}\nz = {round(MovZ, 2)} mm'
    ResTitle = f'{ResTitle}\nk = {FixInd}\nz = {round(ResZ, 2)} mm'
    
    if ExportPlot:
        dpi = 300
    else:
        dpi = 80
        
    plt.subplots(2, 3, figsize=(18,12), dpi=dpi)
    
    alpha = 0.5
    
    BlendedIm = (1.0 - alpha)*FixIm[:,:,FixInd] + alpha*ResIm[:,:,FixInd]
    
    DiffIm = FixIm[:,:,FixInd] - ResIm[:,:,FixInd]
    
    plt.subplot(2, 3, 1)
    plt.imshow(MovPixArr, cmap=plt.cm.Greys_r);
    plt.title(MovTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 2)
    plt.imshow(FixPixArr, cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 3)
    plt.imshow(ResPixArr, cmap=plt.cm.Greys_r);
    plt.title(ResTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 4)
    plt.axis('off')
    
    plt.subplot(2, 3, 5)
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.title('Blended image')
    plt.axis('off')
    
    plt.subplot(2, 3, 6)
    plt.imshow(sitk.GetArrayViewFromImage(DiffIm), cmap=plt.cm.Greys_r);
    plt.title('Difference image')
    plt.axis('off')
    
    #plt.show()
    
    if ExportPlot:
        import time
        import os
        from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print(f'Plot exported to:\n {ExportFpath}\n')
    
    return






def PlotResResultsWithLabIms(FixIm, MovIm, ResIm, MovInd, MovLabIm, ResLabIm, 
                             FixTitle='Fixed image', MovTitle='Moving image', 
                             ResTitle='Resampled moving image',
                             ExportPlot=False, ExportDir='cwd', FileName=''):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. 
    """
    
    #import numpy as np
    import SimpleITK as sitk
    #import matplotlib as mpl
    import matplotlib.pyplot as plt
    #from matplotlib.colors import colorConverter
    from GeneralTools import GetIndOfNearestSlice
    
    #UseCentre = True # 16/04/21 Not sure it's working as expected
    UseCentre = False
    
    FixInd, zDiff = GetIndOfNearestSlice(MovIm, MovInd, FixIm, UseCentre)
    
    FixPixArr = sitk.GetArrayViewFromImage(FixIm)[FixInd,:,:]
    MovPixArr = sitk.GetArrayViewFromImage(MovIm)[MovInd,:,:]
    ResPixArr = sitk.GetArrayViewFromImage(ResIm)[FixInd,:,:]
    
    FixZ = FixIm.TransformIndexToPhysicalPoint([0,0,FixInd])[2]
    MovZ = MovIm.TransformIndexToPhysicalPoint([0,0,MovInd])[2]
    ResZ = ResIm.TransformIndexToPhysicalPoint([0,0,FixInd])[2]
    
    FixTitle = f'{FixTitle}\nz = {round(FixZ, 2)} mm'
    MovTitle = f'{MovTitle}\nz = {round(MovZ, 2)} mm'
    ResTitle = f'{ResTitle}\nz = {round(ResZ, 2)} mm'
    
    if ExportPlot:
        dpi = 300
    else:
        dpi = 80
    
    plt.subplots(2, 3, figsize=(18,12), dpi=dpi)
    
    # Colormap for labelmaps:
    Cmap = plt.cm.nipy_spectral
    #RedCmap = plt.cm.Reds
    #Cmap = plt.cm.bwr
    Cmap = plt.cm.hsv
    Cmap = plt.cm.Reds
    
    # https://stackoverflow.com/questions/10127284/overlay-imshow-plots-in-matplotlib
    #white = colorConverter.to_rgba('white')
    #black = colorConverter.to_rgba('black')

    #GreyCmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap2', [white, black], 256)
    #ColourCmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap', ['red','blue'], 256)
    
    #ColourCmap._init() # create the _lut array, with rgba values
    
    # create your alpha array and fill the colormap with them.
    # here it is progressive, but you can create whathever you want
    #alphas = np.linspace(0, 0.8, ColourCmap.N + 3)
    #ColourCmap._lut[:,-1] = alphas

    
    alpha = 0.5
    alpha = 0.3
    
    BlendedFixResIm = (1.0 - alpha)*FixIm[:,:,FixInd] + alpha*ResIm[:,:,FixInd]
    
    DiffFixResIm = FixIm[:,:,FixInd] - ResIm[:,:,FixInd]
    
    MovSegArr = sitk.GetArrayViewFromImage(MovLabIm)[MovInd,:,:]
    ResSegArr = sitk.GetArrayViewFromImage(ResLabIm)[FixInd,:,:]
    
    #MovSegMask = MovSegArr > 0
    
    #BlendedMovSegIm = (1.0 - alpha)*MovIm[:,:,MovInd] + alpha*MovLabmapIm[:,:,MovInd]
    #BlendedResSegIm = (1.0 - alpha)*ResIm[:,:,FixInd] + alpha*ResLabmapIm[:,:,FixInd]
    
    #plt.subplot(2, 3, 1)
    #plt.imshow(MovPixArr, cmap=plt.cm.Greys_r);
    ##plt.imshow(MovPixArr, cmap=plt.cm.Greys_r, alpha=alpha);
    #plt.imshow(MovSegArr, cmap=Cmap, alpha=alpha);
    #plt.title(MovTitle)
    #plt.axis('off')
    
    
    plt.subplot(2, 3, 1)
    plt.imshow(MovPixArr, cmap=plt.cm.Greys_r);
    #plt.imshow(MovPixArr, cmap=GreyCmap);
    #plt.imshow(MovPixArr, cmap=GreyCmap, alpha=alpha);
    #plt.imshow(MovSegArr, cmap=ColourCmap);
    #plt.imshow(MovSegArr, cmap=RedCmap);
    plt.imshow(MovSegArr, cmap=Cmap, alpha=alpha);
    #plt.imshow(MovSegArr*0.5, cmap=Cmap, alpha=alpha);
    #plt.imshow(MovSegArr*MovSegMask, cmap=Cmap);
    #plt.imshow(MovSegArr*MovSegMask, cmap=Cmap, alpha=alpha);
    plt.title(MovTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 2)
    plt.imshow(FixPixArr, cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 3)
    plt.imshow(ResPixArr, cmap=plt.cm.Greys_r);
    #plt.imshow(ResPixArr, cmap=plt.cm.Greys_r, alpha=alpha);
    plt.imshow(ResSegArr, cmap=Cmap, alpha=alpha);
    plt.title(ResTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 4)
    plt.axis('off')
    
    plt.subplot(2, 3, 5)
    plt.imshow(sitk.GetArrayViewFromImage(BlendedFixResIm), cmap=plt.cm.Greys_r);
    plt.title('Blended image')
    plt.axis('off')
    
    plt.subplot(2, 3, 6)
    plt.imshow(sitk.GetArrayViewFromImage(DiffFixResIm), cmap=plt.cm.Greys_r);
    plt.title('Difference image')
    plt.axis('off')
    
    #plt.show()
    
    if ExportPlot:
        import os
        
        if FileName == '':
            import time
            
            Fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '.jpg'
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
            
        Fpath = os.path.join(ExportDir, Fname)
        
        plt.savefig(Fpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', Fpath)
        
    
    return






def CompareResResults(ResIm0, ResIm1, ResInd,
                      ResTitle0='Resampled image 0', 
                      ResTitle1='Resampled image 1', 
                      ExportPlot=False, ExportDir='cwd', TxtToAddToFname='',
                      ColorbarForDiffIm=False):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. 
    """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt
    
    #UseCentre = True # 16/04/21 Not sure it's working as expected
    UseCentre = False
    
    ResPixArr0 = sitk.GetArrayViewFromImage(ResIm0)[ResInd,:,:]
    ResPixArr1 = sitk.GetArrayViewFromImage(ResIm1)[ResInd,:,:]
    
    if UseCentre:
        # Assuming that ResIm0 size = ResIm1 size:
        R, C, S = ResIm0.GetSize()
        
        i = R//2
        j = C//2
    else:
        i = 0
        j = 0
        
    ResZ = ResIm0.TransformIndexToPhysicalPoint([i,j,ResInd])[2]
    
    ResTitle0 = f'{ResTitle0}\nk = {ResInd}\nz = {round(ResZ, 2)} mm'
    ResTitle1 = f'{ResTitle1}\nk = {ResInd}\nz = {round(ResZ, 2)} mm'
    
    if ExportPlot:
        dpi = 300
    else:
        dpi = 80
        
    plt.subplots(2, 2, figsize=(10,10), dpi=dpi)
    
    alpha = 0.5
    
    BlendedIm = (1.0 - alpha)*ResIm0[:,:,ResInd] + alpha*ResIm1[:,:,ResInd]
    
    DiffIm = ResIm0[:,:,ResInd] - ResIm1[:,:,ResInd]
    
    plt.subplot(2, 2, 1)
    plt.imshow(ResPixArr0, cmap=plt.cm.Greys_r);
    plt.title(ResTitle0)
    plt.axis('off')
    
    plt.subplot(2, 2, 2)
    plt.imshow(ResPixArr1, cmap=plt.cm.Greys_r);
    plt.title(ResTitle1)
    plt.axis('off')
    
    plt.subplot(2, 2, 3)
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.title('Blended image')
    plt.axis('off')
    
    ax = plt.subplot(2, 2, 4)
    im = ax.imshow(sitk.GetArrayViewFromImage(DiffIm), cmap=plt.cm.Greys_r);
    ax.set_title('Difference image')
    ax.axis('off')
    if ColorbarForDiffIm:
        #cbar = plt.colorbar(im, ax=ax)
        #cbar.mappable.set_clim(0, MaxVal)
        plt.colorbar(im, ax=ax)
    
    #plt.show()
    
    if ExportPlot:
        import time
        import os
        from pathlib import Path
        
        if ExportDir == 'cwd':
            ExportDir = os.getcwd()
        
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        ExportFname = CurrentDateTime + '_' + TxtToAddToFname + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
    
    return





