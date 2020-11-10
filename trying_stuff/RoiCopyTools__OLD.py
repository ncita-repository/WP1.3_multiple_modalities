# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:01:08 2020

@author: ctorti
"""



"""
Old functions from RoiCopyTools.py
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






def GetDIVsGroupedByRoi_OLD1(SegRoi):
    """
    Get the DimensionIndexValues in a SEG ROI grouped by common ROI.
    
    Inputs:
        SegRoi - ROI Object from a SEG file
        
        
    Returns:
        DIVs   - (List of lists of integers) Dimension Index Values grouped by
                 ROI
        
        
    Note:
        DIVs is a nested list of lists of the form:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ],
          [ [2, 6], [2, 7], [2, 8], ... ],
          ... 
         ]
          
         
        where each list of index values that belong to a common ROI are 
        separated into separate lists.  If there's only one ROI, the list of
        index values will still be nested in a list to maintain continuity:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ]
         ]
    """
    
    # Get the flat list of DIVs:
    FlatList = GetDIVs(SegRoi) 
        
    # Group the DIVs by common ROI.
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in FlatList:
        RoiNums.append(DIV[0])
        
        
    UniqueRoiNums = list(set(RoiNums))
    
    
    if len(UniqueRoiNums) == 1:
        DimIndVals = [FlatList]
        
    else:
        DimIndVals = []
        
        RoiNum = 1 # initial value
        
        ListThisRoi = [] # initialise
        
        for i in range(len(FlatList)):
            
            if FlatList[i][0] == RoiNum:
                ListThisRoi.append(FlatList[i])
                
            else:
                DimIndVals.append(ListThisRoi)
                
                RoiNum += 1 # increment
                
                ListThisRoi = [] # re-initialise
                
                ListThisRoi.append(FlatList[i])
        
            
        DimIndVals.append(ListThisRoi)
             
    return DimIndVals




def GetDIVsGroupedByRoi_OLD2(SegRoi):
    """
    Get the DimensionIndexValues in a SEG ROI grouped by common ROI.
    
    Inputs:
        SegRoi - ROI Object from a SEG file
        
        
    Returns:
        DIVs   - (List of lists of integers) Dimension Index Values grouped by
                 ROI
        
        
    Note:
        DIVs is a nested list of lists of the form:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ],
          [ [2, 6], [2, 7], [2, 8], ... ],
          ... 
         ]
          
         
        where each list of index values that belong to a common ROI are 
        separated into separate lists.  If there's only one ROI, the list of
        index values will still be nested in a list to maintain continuity:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ]
         ]
    """
    
    # Get the flat list of DIVs:
    FlatList = GetDIVs(SegRoi) 
        
    # Group the DIVs by common ROI.
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in FlatList:
        RoiNums.append(DIV[0])
        
        
    UniqueRoiNums = list(set(RoiNums))
    
    
    if len(UniqueRoiNums) == 1:
        DimIndVals = [FlatList]
        
    else:
        DimIndVals = []
        
        for n in range(len(UniqueRoiNums)):
            ListThisRoi = [] # initialise
            
            for i in range(len(FlatList)):
                if FlatList[i][0] == n:
                    ListThisRoi.append(FlatList[i])
        
            DimIndVals.append(ListThisRoi)
             
    return DimIndVals





def ResampleImage(Image, RefImage, Method, Interpolation):#, RefImageSpacings):
    
    #NewSrcImage = sitk.Resample(image=SrcImage, 
    #                            size=TrgSitkSize,
    #                            transform=sitk.Transform(),
    #                            interpolator=sitk.sitkLinear,
    #                            outputOrigin=TrgSitkOrigin,
    #                            outputSpacing=TrgSitkSpacing,
    #                            outputDirection=TrgSitkDirection,
    #                            defaultPixelValue=0,
    #                            outputPixelType=TrgImage.GetPixelID())
    
    
    # Define which interpolator to use:
    """ 
    Use the linear (or BSpline) interpolator for intensity images, and
    NearestNeighbour for binary images (e.g. segmentation) so that no new
    labels are introduced.
    """
    
    if 'inear' in Interpolation:
        Interpolator = sitk.sitkLinear
    if 'pline' in Interpolation:
        Interpolator = sitk.sitkBSpline
    if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:    
        Interpolator = sitk.sitkNearestNeighbor
            
        
    #print('\nUsing', Interpolation, 'interpolation\n')
    
    Dimension = Image.GetDimension()
    

    
    if 'ResampleImageFilter' in Method:
        Resampler = sitk.ResampleImageFilter()
        
        Resampler.SetReferenceImage(RefImage)
        
        Resampler.SetInterpolator(Interpolator)
        
        if 'Identity' in Method:
            Resampler.SetTransform(sitk.Transform())
        if 'Affine' in Method:
            Resampler.SetTransform(sitk.AffineTransform(Dimension))
            
            
        #print(f'\nImage.GetPixelIDValue() = {Image.GetPixelIDValue()}')
        
        Resampler.SetOutputSpacing(RefImage.GetSpacing())
        Resampler.SetSize(RefImage.GetSize())
        Resampler.SetOutputDirection(RefImage.GetDirection())
        Resampler.SetOutputOrigin(RefImage.GetOrigin())
        #Resampler.SetDefaultPixelValue(Image.GetPixelIDValue())
        Resampler.SetDefaultPixelValue(0)
        
        #print('\n')
        #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: print("\rProgress: {0:03.1f}%...".format(100*Resampler.GetProgress()),end=''))
        #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: sys.stdout.flush())
        
        ResImage = Resampler.Execute(Image)
    
    
    
    if 'Resample2' in Method:
        """
        Use the 2nd method listed in the SimpleITK Transforms and Resampling 
        tutorial (the new Size, Spacing, Origin and Direction will be 
        determined intrinsically from RefImage):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        
        Also follow the example here:
            
        https://gist.github.com/zivy/79d7ee0490faee1156c1277a78e4a4c4
        """
        
        # Define the transform which maps RefImage to Image with the 
        # translation mapping the image origins to each other:
        Transform = sitk.AffineTransform(Dimension)
        #Transform.SetMatrix(Image.GetDirection())                                      # <-- 9:55 on 30/10
        Transform.SetMatrix(RefImage.GetDirection())                                    # <-- 9:55 on 30/10
        #Transform.SetTranslation(np.array(Image.GetOrigin()) - RefImage.GetOrigin())   # <-- 9:55 on 30/10
        Transform.SetTranslation(np.array(RefImage.GetOrigin()) - Image.GetOrigin())    # <-- 9:55 on 30/10
        
        if 'Centres' in Method:
            # Use the TransformContinuousIndexToPhysicalPoint to compute an 
            # indexed point's physical coordinates as this takes into account  
            # size, spacing and direction cosines. 
            """
            For the vast majority of images the direction cosines are the 
            identity matrix, but when this isn't the case simply multiplying 
            the central index by the spacing will not yield the correct 
            # coordinates resulting in a long debugging session. 
            """
            RefCentre = np.array(RefImage.TransformContinuousIndexToPhysicalPoint(np.array(RefImage.GetSize())/2.0))
            
            # Modify the transformation to align the centres of the original 
            # and reference image instead of their origins:
            CenteringTransform = sitk.TranslationTransform(Dimension)
            
            ImageCentre = np.array(Image.TransformContinuousIndexToPhysicalPoint(np.array(Image.GetSize())/2.0))
            
            CenteringTransform.SetOffset(np.array(Transform.GetInverse().TransformPoint(ImageCentre) - RefCentre))
            
            CenteredTransform = sitk.Transform(Transform)
            
            CenteredTransform.AddTransform(CenteringTransform)

            ResImage = sitk.Resample(Image, RefImage, CenteredTransform, 
                                     Interpolator, Image.GetPixelIDValue())
            
            
        if 'Origins' in Method:
            ResImage = sitk.Resample(Image, RefImage, Transform, Interpolator, 
                                     Image.GetPixelIDValue())
        
        
    
    
    if 'Resample3' in Method:
        """
        Use the 3rd method listed in the SimpleITK Transforms and Resampling 
        tutorial (explicitly define the new Size, Spacing, Origin and 
        Direction):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        """
        
        NewSize = RefImage.GetSize()
        NewOrigin = RefImage.GetOrigin()
        NewSpacing = RefImage.GetSpacing()
        NewDirection = RefImage.GetDirection()
        
        
        # Define the transform which maps RefImage to Image with the 
        # translation mapping the image origins to each other:
        Transform = sitk.AffineTransform(Dimension)
        #Transform.SetMatrix(Image.GetDirection())                                     # <-- 10:03 on 30/10
        Transform.SetMatrix(RefImage.GetDirection())                                   # <-- 10:03 on 30/10
        #Transform.SetTranslation(np.array(Image.GetOrigin()) - RefImage.GetOrigin())  # <-- 10:03 on 30/10
        Transform.SetTranslation(np.array(RefImage.GetOrigin()) - Image.GetOrigin())   # <-- 10:03 on 30/10
        
        if 'Centres' in Method:
            # Use the TransformContinuousIndexToPhysicalPoint to compute an 
            # indexed point's physical coordinates as this takes into account  
            # size, spacing and direction cosines. 
            """
            For the vast majority of images the direction cosines are the 
            identity matrix, but when this isn't the case simply multiplying 
            the central index by the spacing will not yield the correct 
            # coordinates resulting in a long debugging session. 
            """
            RefCentre = np.array(RefImage.TransformContinuousIndexToPhysicalPoint(np.array(RefImage.GetSize())/2.0))
            
            # Modify the transformation to align the centres of the original 
            # and reference image instead of their origins:
            CenteringTransform = sitk.TranslationTransform(Dimension)
            
            ImageCentre = np.array(Image.TransformContinuousIndexToPhysicalPoint(np.array(Image.GetSize())/2.0))
            
            CenteringTransform.SetOffset(np.array(Transform.GetInverse().TransformPoint(ImageCentre) - RefCentre))
            
            CenteredTransform = sitk.Transform(Transform)
            
            CenteredTransform.AddTransform(CenteringTransform)
        
            ResImage = sitk.Resample(Image, NewSize, CenteredTransform, 
                                     Interpolator, NewOrigin, NewSpacing, 
                                     NewDirection, Image.GetPixelIDValue())  
            
            
        if 'Origins' in Method:
            ResImage = sitk.Resample(Image, NewSize, Transform, 
                                     Interpolator, NewOrigin, NewSpacing, 
                                     NewDirection, Image.GetPixelIDValue())  
            
            
            
            
    if Method == 'NotCorrect':
        """ 
        Use the 3rd method listed in the SimpleITK Transforms and Resampling 
        tutorial (explicitly define the new Size, Spacing, Origin and 
        Direction):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        """
        
        ExtremePts = GetExtremePoints(RefImage)
        
        Affine = sitk.Euler3DTransform(RefImage.TransformContinuousIndexToPhysicalPoint([size/2 for size in RefImage.GetSize()]),
                                       RefImage.GetAngleX(), 
                                       RefImage.GetAngleY(), 
                                       RefImage.GetAngleZ(), 
                                       RefImage.GetTranslation())
        
        InvAffine = Affine.GetInverse()
        
        # Transformed ExtremePts:
        TxExtremePts = [InvAffine.TransformPoint(pt) for pt in ExtremePts]
        
        min_x = min(TxExtremePts)[0]
        min_y = min(TxExtremePts, key=lambda p: p[1])[1]
        min_z = min(TxExtremePts, key=lambda p: p[2])[2]
        max_x = max(TxExtremePts)[0]
        max_y = max(TxExtremePts, key=lambda p: p[1])[1]
        max_z = max(TxExtremePts, key=lambda p: p[2])[2]
        
        NewSpacing = RefImage.GetSpacing()
        #NewDirection = np.eye(3).flatten()
        NewDirection = RefImage.GetDirection()
        #NewOrigin = [min_x, min_y, min_z]
        NewOrigin = RefImage.GetOrigin()
        
        #NewSize = [math.ceil((max_x - min_x)/NewSpacing[0]), 
        #           math.ceil((max_y - min_y)/NewSpacing[1]),
        #           math.ceil((max_z - min_z)/NewSpacing[2])]
        NewSize = RefImage.GetSize()
        
        ResImage = sitk.Resample(Image, NewSize, Affine, Interpolator, 
                                 NewOrigin, NewSpacing, NewDirection)
    
    
    #sitk.Show(ResImage)

    return ResImage








def ExportSegRoi(TrgSegRoi, SrcSegFpath, NamePrefix, ExportDir):
    """
    Export contour ROI Object to disc.
    
    
    Inputs:
        TrgSegRoi   - Target SEG ROI object to be exported
        
        SrcSegFpath - (String) Full path of the Source DICOM SEG file (used to
                      generate the filename of the new SEG file)
        
        NamePrefix  - (String) Prefix to be added to the assigned filename and
                      Content Label (e.g. 'MR4_S9_s23_to_MR4_S9_s22') 
                            
        ExportDir   - (String) Directory where the SEG is to be exported
                           
    
                            
    Returns:
        TrgSegFpath - (String) Full path of the exported Target DICOM SEG file
    
    """
    
    """
    The following was moved into the main function ModifySegTagVals:
    # Generate a new SOP Instance UID:
    SegRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    SegRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    """
    
    # Get the filename of the original SEG file:
    SrcSegFname = os.path.split(SrcSegFpath)[1]
    
    # Modify the Content Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    #TrgSegRoi.ContentDate = NewDate
    #TrgSegRoi.ContentTime = NewTime
    
    FnamePrefix = NewDate + '_' + NewTime + '_' + NamePrefix + '_from_' 
    
    ContentLabelPrefix = NamePrefix + '_from_'
    
    # Modify the Content Label (this appears under Name in XNAT):
    #TrgSegRoi.ContentLabel = 'Copy_of_' + TrgSegRoi.ContentLabel
    TrgSegRoi.ContentLabel = ContentLabelPrefix + TrgSegRoi.ContentLabel
    
    # Modify Content Description:
    TrgSegRoi.ContentDescription = ContentLabelPrefix \
                                   + TrgSegRoi.ContentDescription
                                   
    # Modify the Series Description as with the Content Label:
    TrgSegRoi.SeriesDescription = ContentLabelPrefix \
                                  + TrgSegRoi.SeriesDescription
    
    
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = 'New_' + SrcSegFname + '_' + NewDate + '_' + NewTime
    #TrgRoiFname = FnamePrefix + '_SEG_' + NewDate + '_' + NewTime + '_' \
    #              + SrcSegFname
    TrgSegFname = FnamePrefix + SrcSegFname
    
    TrgSegFpath = os.path.join(ExportDir, TrgSegFname)
    
    TrgSegRoi.save_as(TrgSegFpath)
        
    print('\nSEG ROI exported to:\n\n', TrgSegFpath)
    
    return TrgSegFpath






def Compare2dMasksFrom3dMasksOLD(OrigSegRoi, NewSegRoi, OrigSliceNum, NewSliceNum,
                              OrigDicomDir, NewDicomDir):
    """ Compare cropped masks of non-zero elements only. """ 
    
    # Get the DICOM SOP UIDs:
    OrigSOPuids = GetDicomSOPuids(DicomDir=OrigDicomDir)
    NewSOPuids = GetDicomSOPuids(DicomDir=NewDicomDir)
    
    # Get the Per-frameFunctionalGroupsSequence-to-DICOM slice indices:
    OrigPFFGStoDcmInds = GetPFFGStoDcmInds(OrigSegRoi, OrigSOPuids)
    NewPFFGStoDcmInds = GetPFFGStoDcmInds(NewSegRoi, NewSOPuids)
    
    # Get the 3D SEG masks:
    Orig3dMask = OrigSegRoi.pixel_array
    New3dMask = NewSegRoi.pixel_array
    
    OrigShape = Orig3dMask.shape
    NewShape = New3dMask.shape
    
    print(f'Segments exist in OrigSegRoi on slices {OrigPFFGStoDcmInds}')
    print(f'Shape of Orig3dMask = {OrigShape}')
    print(f'\nSegments exist in NewSegRoi on slices {NewPFFGStoDcmInds}')
    print(f'Shape of New3dMask = {NewShape}\n')
    
    
    # Initialise Orig2dMaskCropped and New2dMaskCropped:
    Orig2dMaskCropped = np.zeros((1, 1))
    New2dMaskCropped = np.zeros((1, 1))
    
    
    # Does slice OrigSliceNum in OrigSegRoi have a segment?
    if OrigSliceNum in OrigPFFGStoDcmInds:
        # A segment frame exists for OrigSliceNum:
        OrigFrame = True
        
        # Get the frame number for this slice number:
        OrigFrameNum = OrigPFFGStoDcmInds.index(OrigSliceNum)
        
        Orig2dMask = OrigSegRoi.pixel_array[OrigFrameNum]
        
        print(f'Shape of Orig2dMask = Orig3dMask[{OrigFrameNum}] = {Orig2dMask.shape}')
    
        #print('Cropping Orig2dMask to array of non-zero elements..')
    
        Orig2dMaskCropped = CropNonZerosIn2dMask(Orig2dMask)
        #Orig2dMaskCropped = CropNonZerosIn2dMask(Orig3dMask[OrigFrameNum])
        
    else:
        OrigFrame = False
        
        print(f'A segment does not exist on slice {OrigSliceNum} in OrigSegRoi.')
        
    
    # Does slice NewSliceNum in NewSegRoi have a segment?
    if NewSliceNum in NewPFFGStoDcmInds:
        # A segment frame exists for NewSliceNum:
        NewFrame = True
        
        # Get the frame number for this slice number:
        NewFrameNum = NewPFFGStoDcmInds.index(NewSliceNum)
        
        New2dMask = NewSegRoi.pixel_array[NewFrameNum]
    
        print(f'Shape of New2dMask = New3dMask[{NewFrameNum}]  = {New2dMask.shape}')
    
        #print('Cropping New2dMask to array of non-zero elements..')
        
        New2dMaskCropped = CropNonZerosIn2dMask(New2dMask)
        #New2dMaskCropped = CropNonZerosIn2dMask(New3dMask[FrameNum])
        
        
    else:
        NewFrame = False
        
        print(f'A segment does not exist on slice {NewSliceNum} in NewSegRoi.')
    
    
    # Proceed only if either cropped arrays are not empty (i.e. if at least one
    # of the 2D masks had non-zero elements):
    if Orig2dMaskCropped.any() or New2dMaskCropped.any():
    
        print('')
        
        fig, ax = plt.subplots(1, 2)
        
        if Orig2dMaskCropped.any():
            print('Shape of Orig2dMaskCropped =', Orig2dMaskCropped.shape)
            
            ax = plt.subplot(1, 2, 1, aspect='equal')
            ax.imshow(Orig2dMaskCropped)
        
        else:
            if OrigFrame:
                print(f'Frame {OrigFrameNum} in OrigSegRoi does not have any',
                      'non-zero elements.')
        
        if New2dMaskCropped.any():
            print('Shape of New2dMaskCropped  =', New2dMaskCropped.shape)
        
            ax = plt.subplot(1, 2, 2, aspect='equal')
            ax.imshow(New2dMaskCropped)
        
        else:
            if NewFrame:
                print(f'Frame {NewFrameNum} in NewSegRoi does not have any',
                      'non-zero elements.')
            
    else:
        if OrigFrame and NewFrame:
            print(f'Frame {OrigFrameNum} and {NewFrameNum} does not have any',
                  'non-zero elements in OrigSegRoi and NewSegRoi.')
    
    return







def DELETE_THIS_AddToRIStoDcmInds(RIStoDcmInds, ToSliceNum):
    """
    Append the DICOM slice number (ToSliceNum) representing the segmentation to 
    be added to the Referenced Instance Sequence to RIStoDcmInds (the 
    ReferencedInstanceSequence-to-DICOM slice indeces), and sort the list of 
    indeces.
    
    Inputs:
        RIStoDcmInds    - List of integers of the DICOM slice numbers that 
                          correspond to each Referenced Instance Sequence
        
        ToSliceNum      - (Integer) Slice index in the Target DICOM stack where
                          the segmentation will be copied to (counting from 0)
        
    Returns:
        NewRIStoDcmInds - Modified and sorted list of RIStoDcmInds with the  
                          slice index to be copied
    """
    
    NewRIStoDcmInds = copy.deepcopy(RIStoDcmInds)

    NewRIStoDcmInds.append(ToSliceNum)
    
    NewRIStoDcmInds.sort()
    
    return NewRIStoDcmInds




def DELETE_THIS_AddToPFFGStoDcmInds(PFFGStoDcmInds, ToSliceNum):
    """
    Append the DICOM slice number (ToSliceNum) representing the segmentation to 
    be added to the Per-frame Functional Groups Sequence to PFFGStoDcmInds (the 
    PerFrameFunctionalGroupsSequence-to-DICOM slice indeces), and sort the list 
    of indeces.
    
    Inputs:
        PFFGStoDcmInds    - List of integers of the DICOM slice numbers that 
                            correspond to each Per-frame Functional Groups 
                            Sequence
        
        ToSliceNum        - (Integer) Slice index in the Target DICOM stack 
                            where the segmentation will be copied to (counting 
                            from 0)
        
    Returns:
        NewPFFGStoDcmInds - Modified and sorted list of PFFGStoDcmInds with the 
                            slice index to be copied
    """
    
    NewPFFGStoDcmInds = copy.deepcopy(PFFGStoDcmInds)

    NewPFFGStoDcmInds.append(ToSliceNum)
    
    NewPFFGStoDcmInds.sort()
    
    return NewPFFGStoDcmInds