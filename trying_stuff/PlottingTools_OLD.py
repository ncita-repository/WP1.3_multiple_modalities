# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:34:17 2020

@author: ctorti
"""



def Plot_Src_ResSrc_Trg_Images_OLD(SrcIm, ResSrcIm, TrgIm, SrcImLabel, TrgImLabel,
                               ImageType, ExportPlot, ExportDir, 
                               LogToConsole=False, dpi=80):
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm.GetSize()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetSize()} \n',
                      f'\n   Target:           {TrgIm.GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm.GetSpacing()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetSpacing()} \n',
                      f'\n   Target:           {TrgIm.GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm.GetOrigin()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetOrigin()} \n',
                      f'\n   Target:           {TrgIm.GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm.GetDirection()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetDirection()} \n',
                      f'\n   Target:           {TrgIm.GetDirection()}')

    
    
    # Convert from SimpleITK image to Numpy arrays:
    SrcNpa = sitk.GetArrayFromImage(SrcIm)
    ResSrcNpa = sitk.GetArrayFromImage(ResSrcIm)
    TrgNpa = sitk.GetArrayFromImage(TrgIm)
    
    # Prepare the figure:
    Nslices = max(SrcIm.GetSize()[2], TrgIm.GetSize()[2])
    
    # Set the number of subplot rows and columns:
    Ncols = 3
    
    Nrows = Nslices
    
    
    # Create a figure with two subplots and the specified size:
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Loop through each slice: 
    for i in range(Nslices):
        #if LogToConsole:
        #        print('\ni =', i)
        
        #plt.subplot(Nrows, Ncols, n)
        #plt.axis('off')
        #ax = plt.subplot(Nrows, Ncols, n)

        #ax = plt.subplot(Nrows, Ncols, n, aspect=AR)
         
        # Plot the Original Source slice:
        if i < SrcIm.GetSize()[2]:
            #maximum = SrcNpa[i].max()
            maximum = np.amax(SrcNpa[i])
            if maximum <= 1:
                maximum = 1
            
            ax = plt.subplot(Nrows, Ncols, n)
            im = ax.imshow(SrcNpa[i], cmap=plt.cm.Greys_r)
            cbar = fig.colorbar(im, ax=ax)
            cbar.mappable.set_clim(0, maximum)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Original Source slice {i}')
            
        n += 1 # increment sub-plot number
            
        # Plot the Resampled Source slice:
        if i < ResSrcIm.GetSize()[2]:
            #maximum = ResSrcNpa[i].max()
            maximum = np.amax(ResSrcNpa[i])
            if maximum <= 1:
                maximum = 1
                
            ax = plt.subplot(Nrows, Ncols, n)
            im = ax.imshow(ResSrcNpa[i], cmap=plt.cm.Greys_r)
            cbar = fig.colorbar(im, ax=ax)
            cbar.mappable.set_clim(0, maximum)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Resampled Source slice {i}')
            
        n += 1 # increment sub-plot number
        
        # Plot the Target slice:
        if i < TrgIm.GetSize()[2]:
            #maximum = TrgNpa[i].max()
            maximum = np.amax(TrgNpa[i])
            if maximum <= 1:
                maximum = 1
                
            ax = plt.subplot(Nrows, Ncols, n)
            im = ax.imshow(TrgNpa[i], cmap=plt.cm.Greys_r)
            cbar = fig.colorbar(im, ax=ax)
            cbar.mappable.set_clim(0, maximum)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Target slice {i}')
            
        n += 1 # increment sub-plot number
        
      
    
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
    
    if ExportPlot:
        # Short version of Interpolation for exported filename:
        #if 'inear' in Interpolation:
        #    Interp = 'Lin'
        #if 'pline' in Interpolation:
        #    Interp = 'BSp'
        #if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:    
        #    Interp = 'NN'
            
        # Short version of ResampleImageFilter:
        #if Method == 'ResampleImageFilter':
        #    Method = 'ResampleImFilt'
            
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        #ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
        #              + SrcImLabel + '_' + TrgImLabel + '_images_' \
        #              + Method + '_' + Interp + '.jpg'
                      
        ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
                      + SrcImLabel + '_' + TrgImLabel + '_' + ImageType \
                      + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return







def Plot_Src_ResSrc_Trg_NewTrg_Images(SrcIm, ResSrcIm, TrgIm, NewTrgIm,
                                      SrcImLabel, TrgImLabel,
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
        
    if type(NewTrgIm).__name__ == 'Image':
        NewTrgIm = [NewTrgIm]
    
    
    if LogToConsole:
        print(f'\nScrIm has {len(SrcIm)} segments')
        print(f'ResScrIm has {len(ResSrcIm)} segments')
        print(f'TrgIm has {len(TrgIm)} segments')
        print(f'NewTrgIm has {len(NewTrgIm)} segments')
    
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSize()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetSize()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetSize()}',
                      f'\n   New Target:       {NewTrgIm[0].GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSpacing()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetSpacing()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetSpacing()}',
                      f'\n   New Target:       {NewTrgIm[0].GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm[0].GetOrigin()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetOrigin()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetOrigin()}',
                      f'\n   New Target:       {NewTrgIm[0].GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm[0].GetDirection()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetDirection()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetDirection()}',
                      f'\n   New Target:       {NewTrgIm[0].GetDirection()}')

    
    """
    The code below works but it's less efficient than what follows, i.e.
    instead of converting from sitk images to numpy and adding the segments,
    it's more efficient to add the segments as sitk images the convert to 
    numpy.
    
    # Convert from SimpleITK image to Numpy arrays, assigning different values
    # for each ROI (i.e. 1, 2, 3, ...):
    SrcNpas = []
    ResSrcNpas = []
    TrgNpas = []
    NewTrgNpas = []
    
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
        
    for r in range(len(NewTrgIm)):
        NewTrgNpas.append(sitk.GetArrayFromImage(NewTrgIm[r]))
            
    
    if LogToConsole:
        print(f'\nlen(SrcNpas)    = {len(SrcNpas)}')
        print(f'len(ResSrcNpas) = {len(ResSrcNpas)}')
        print(f'len(TrgNpas)    = {len(TrgNpas)}')
        print(f'len(NewTrgNpas) = {len(NewTrgNpas)}')
        
        print(f'\nSrcNpas[0].shape    = {SrcNpas[0].shape}')
        print(f'ResSrcNpas[0].shape = {ResSrcNpas[0].shape}')
        print(f'TrgNpas[0].shape    = {TrgNpas[0].shape}')
        print(f'NewTrgNpas[0].shape = {NewTrgNpas[0].shape}')
        
    # Take the sum of the labelmaps for all ROIs:
    SrcNpaSum = np.zeros((SrcNpas[0].shape))
    ResSrcNpaSum = np.zeros((ResSrcNpas[0].shape))
    TrgNpaSum = np.zeros((TrgNpas[0].shape))
    NewTrgNpaSum = np.zeros((NewTrgNpas[0].shape))
    
    for r in range(len(SrcNpas)):
        SrcNpaSum += SrcNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(SrcNpaSum) = {np.amax(SrcNpaSum)}')
        
    for r in range(len(ResSrcNpas)):
        ResSrcNpaSum += ResSrcNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(ResSrcNpaSum) = {np.amax(ResSrcNpaSum)}')
        
    for r in range(len(TrgNpas)):
        TrgNpaSum += TrgNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(TrgNpaSum) = {np.amax(TrgNpaSum)}')
        
    for r in range(len(NewTrgNpas)):
        NewTrgNpaSum += NewTrgNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(NewTrgNpaSum) = {np.amax(NewTrgNpaSum)}')
        
        
    print(f'\nSrcNpaSum.shape    = {SrcNpaSum.shape}')
    print(f'ResSrcNpaSum.shape = {ResSrcNpaSum.shape}')
    print(f'TrgNpaSum.shape    = {TrgNpaSum.shape}')
    print(f'NewTrgNpaSum.shape = {NewTrgNpaSum.shape}')
    """
    
    
    # Add the segments contained in the list of sitk images into a single sitk
    # image:
    SrcImSum = InitialiseImage(SrcIm[0])
    for r in range(len(SrcIm)):
        SrcImSum = AddImages(SrcImSum, SrcIm[r])
        
    ResSrcImSum = InitialiseImage(ResSrcIm[0])
    for r in range(len(ResSrcIm)):
        ResSrcImSum = AddImages(ResSrcImSum, ResSrcIm[r])
        
    TrgImSum = InitialiseImage(TrgIm[0])
    for r in range(len(TrgIm)):
        TrgImSum = AddImages(TrgImSum, TrgIm[r])
       
    NewTrgImSum = InitialiseImage(NewTrgIm[0])
    for r in range(len(NewTrgIm)):
        NewTrgImSum = AddImages(NewTrgImSum, NewTrgIm[r])
    
        
    # Convert to numpy arrays:
    SrcNpaSum = sitk.GetArrayFromImage(SrcImSum)
    ResSrcNpaSum = sitk.GetArrayFromImage(ResSrcImSum)
    TrgNpaSum = sitk.GetArrayFromImage(TrgImSum)
    NewTrgNpaSum = sitk.GetArrayFromImage(NewTrgImSum)
    
    
    #Nslices = max(SrcIm.GetSize()[2], TrgIm.GetSize()[2])

    #SrcNslices = SrcNpas[0].shape[0]
    #ResSrcNslices = ResSrcNpas[0].shape[0]
    #TrgNslices = TrgNpas[0].shape[0]
    #NewTrgNslices = NewTrgNpas[0].shape[0]
    
    SrcNslices = SrcNpaSum.shape[0]
    ResSrcNslices = ResSrcNpaSum.shape[0]
    TrgNslices = TrgNpaSum.shape[0]
    NewTrgNslices = NewTrgNpaSum.shape[0]
    
    Nslices = max(SrcNslices, TrgNslices)

    
    # Get the maxima for all slices for all labelmaps:
    if 'Labelmaps' in ImageType:
        SrcMaximaBySlice = np.amax(SrcNpaSum, axis=(1,2))
        
        ResSrcMaximaBySlice = np.amax(ResSrcNpaSum, axis=(1,2))
            
        TrgMaximaBySlice = np.amax(TrgNpaSum, axis=(1,2))
        
        NewTrgMaximaBySlice = np.amax(NewTrgNpaSum, axis=(1,2))
           
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
    
    # Define the colourmaps to use for each ROI:
    #cmaps = [plt.cm.spring, plt.cm.cool, plt.cm.autumn, plt.cm.winter, 
    #         plt.cm.hsv, plt.cm.jet]
    
    # Loop through each slice: 
    for s in range(Nslices):
        # Only proceed if ImageType is not 'Labelmaps' or (if it is) at least 
        # one of the labelmaps is non-zero for this slice number:
        if not 'Labelmaps' in ImageType or SumOfMaximaBySlice[s]:
            # Plot the Original Source slice:
            if s < SrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                # Loop through all ROIs:
                #maxima = []
                #for r in range(len(SrcNpas)):
                #    im = ax.imshow(SrcNpas[r][s], cmap=plt.cm.Greys_r)
                #    #im = ax.imshow(SrcNpas[r][s], cmap=cmaps[r])
                #    #maxima.append(np.amax(SrcNpas[r][s]))
                
                im = ax.imshow(SrcNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if SrcMaximaBySlice[s] < 1 else SrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Original Source slice {s}')
                
            n += 1 # increment sub-plot number
                
            # Plot the Resampled Source slice:
            if s < ResSrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                #maxima = []
                #npa = np.zeros()
                #for r in range(len(ResSrcNpas)):
                #    im = ax.imshow(ResSrcNpas[r][s], cmap=plt.cm.Greys_r)
                #    #im = ax.imshow(ResSrcNpas[r][s], cmap=cmaps[r])
                #    #maxima.append(np.amax(ResSrcNpas[r][s]))
                
                im = ax.imshow(ResSrcNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if ResSrcMaximaBySlice[s] < 1 else ResSrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Resampled Source slice {s}')
                
            n += 1 # increment sub-plot number
            
            # Plot the Target slice:
            if s < TrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                #maxima = []
                #for r in range(len(TrgNpas)):
                #    im = ax.imshow(TrgNpas[r][s], cmap=plt.cm.Greys_r)
                #    #im = ax.imshow(TrgNpas[r][s], cmap=cmaps[r])
                #    #maxima.append(np.amax(TrgNpas[r][s]))
                
                im = ax.imshow(TrgNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if TrgMaximaBySlice[s] < 1 else TrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Target slice {s}')
                
            n += 1 # increment sub-plot number
            
            # Plot the New Target slice:
            if s < NewTrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                #for r in range(len(NewTrgNpas)):
                #    im = ax.imshow(NewTrgNpas[r][s], cmap=plt.cm.Greys_r)
                
                im = ax.imshow(NewTrgNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
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








def PlotDicomsAndContours_OLD(DicomDir, RtsFpath, ExportPlot, ExportDir, 
                          LogToConsole=False):
    
    from pydicom import dcmread
    import numpy as np
    import matplotlib.pyplot as plt
    import time
    import os
    from DicomTools import ImportDicoms
    from DicomTools import GetDicomSOPuids
    from RtsTools import GetCIStoSliceInds
    from ConversionTools import ConvertContourDataToIndices
    from GeneralTools import Unpack
    
    
    # Import the DICOMs:
    Dicoms = ImportDicoms(DicomDir=DicomDir)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir=DicomDir)
    
    # Import the RTS ROI:
    RtsRoi = dcmread(RtsFpath)
    
    # Get the ContourImageSequence-to-DICOM slice indices:
    CIStoSliceInds = GetCIStoSliceInds(RtsRoi, SOPuids)
    
    # Get a list of the Contour Sequences in the ROI:
    ContourSequences = RtsRoi.ROIContourSequence[0].ContourSequence
    
    
    # Prepare the figure.
    
    # All DICOM slices that have contours will be plotted. The number of
    # slices to plot:
    Nslices = len(CIStoSliceInds)
    
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
        s = CIStoSliceInds[i]
        
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
        PointsICS = ConvertContourDataToIndices(ContourData, DicomDir)
        
        # Unpack the points:
        X, Y, Z = Unpack(PointsICS)
        
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




