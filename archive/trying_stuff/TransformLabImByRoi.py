def TransformLabImByRoi(LabImByRoi, F2SindsByRoi, TrgIm, 
                        SelxImFiltOrSitkTx, Interp='NearestNeighbor', 
                        #ApplyPostTxBlur=True, PostTxVariance=(1,1,1),
                        #ApplyPostTxBin=True, #VoxelVolumeRatio=1,
                        ApplyPostTxBlur=False, PostTxVariance=(1,1,1),
                        ApplyPostTxBin=False, #VoxelVolumeRatio=1,
                        LogToConsole=False):
    """
    08/06/21: Default values of ApplyPostTxBlur and ApplyPostTxBin changed
    from True to False due to strang results obtained when testing on CT and MR
    data involving BSpline registration with fiducials alignment.
    
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
    LabImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled. The 
        labelmaps share the same image attributes as MovImage.
    F2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the frame-to-slice 
        indices of each labelmap image in LabImByRoi.
    TrgIm : SimpleITK image
        The 3D image whose grid LabImByRoi is to be transformed to.
    SelxImFiltOrSitkTx : Elastix Image Filter or SimpleITK Transform
        The Elastix Image Filter used to perform image registration (if Elastix
        was used); or the final SimpleITK Transform if registration was 
        performed using SimpleITK.
    Interp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image(s) will be Gaussian blurred.
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image(s) will be binary thresholded.
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Outputs:
    *******
    
    TxLabImByRoi : list of SimpleITK images
        The list (for each ROI) of the 3D labelmap images in LabImByRoi
        transformed to the FixImage grid space.
    TxPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of pixel arrays (converted from each labelmap 
        in TxLabImByRoi).
    TxF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in TxPixArrByRoi.
    
    Notes:
    *****
    If the transformed binary labelmap image is blurred, the resulting labelmap
    image will be non-binary, so will need to be binarised. The threshold used
    will be determined such that the resulting (re-) binarised image will have
    a similar number of truth values (ones) as the pre-transformed (binary)
    one.
    """
    
    import time
    #import SimpleITK as sitk
    from ConversionTools import Image2PixArr
    from GeneralTools import UniqueItems
    from SelxTools import TransformImageSelx
    from copy import deepcopy
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of TransformLabImByRoi():')
        print('\n\n', '-'*120)
    
    """ Start timing. """
    times = []
    times.append(time.time())
    
    """ Transform LabImByRoi to the Target image's grid to get the resampled
    pixel array. """
    TxLabImByRoi = []
    TxPixArrByRoi = []
    TxF2SindsByRoi = [] 
    
    for r in range(len(LabImByRoi)):
        LabIm = LabImByRoi[r]
        
        if LogToConsole:
            print('\n   Image info for LabIm prior to transformation:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(LabIm, LogToConsole)
        
        OrigF2Sinds = deepcopy(F2Sinds)
        
        """ Transform LabIm. """
        if 'ElastixImageFilter' in str(type(SelxImFiltOrSitkTx)):
            TxLabIm = TransformImageSelx(MovIm=LabIm, 
                                         SelxImFilt=SelxImFiltOrSitkTx,
                                         Interp=Interp, 
                                         LogToConsole=LogToConsole)
        elif 'Transform' in str(type(SelxImFiltOrSitkTx)):
            TxLabIm = TransformImageSitk(MovIm=LabIm, FixIm=TrgIm, 
                                         SitkTransform=SelxImFiltOrSitkTx, 
                                         Interp=Interp, 
                                         LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'*Took {Dtime} s to transform the {r}^th labelmap image',
              'using the image registration transformation.\n')
        
        if LogToConsole:
            print('\n   Image info for LabIm after transformation:')
        PixID, PixIDTypeAsStr, UniqueValsPostTx,\
        F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
        if LogToConsole:
            from PlottingTools import PlotTwoImages
            
            PlotTwoImages(Im0=LabIm, Ind0=OrigF2Sinds[0], 
                          PlotLabel0='Original label image', 
                          Im1=TxLabIm, Ind1=F2Sinds[0], 
                          PlotLabel1='Resampled label image')
        
        if ApplyPostTxBlur:
            #""" The binary pixel array: """
            #PixArr_B, F2Sinds = Image2PixArr(TxLabIm)
            
            # 11/06/21:
            """ Ensure TxLabIm is a 32-bit float (PixID = 8). """
            if PixID != 8: 
                if LogToConsole:
                    print(f'\n   TxLabIm has PixelID = {PixID} '\
                          f'({PixIDTypeAsStr})). It will be converted to '\
                          'Float32.')
                    
                """ Convert TxLabIm to Float32: """
                TxLabIm = ChangeImagePixelType(Image=TxLabIm, 
                                               NewPixelType='Float32')
                
                if LogToConsole:
                    print('\n   Image info for TxLabIm after converting to Float32:')
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
            
            """ Gaussian blur TxLabIm (new 14/02/21) to reduce noise that 
            makes selection of an appropriate threshold for binarisation (below)
            difficult. """
            #TxLabIm = GaussianBlurImage(TxLabIm)
            #TxLabIm = GaussianBlurImage(TxLabIm, Variance=PostTxVariance) # 11/06/21
            BlurTxLabIm = GaussianBlurImage(TxLabIm, Variance=PostTxVariance,
                                            LogToConsole=LogToConsole) # 11/06/21
                
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            print(f'*Took {Dtime} s to Gaussian blur the {r}^th labelmap image.\n')
            
            if LogToConsole:
                print('\n   Image info for BlurTxLabIm after Gaussian blur:')
                
            #PixID, PixIDTypeAsStr, UniqueValsPostBlur,\
            #F2Sinds = GetImageInfo(TxLabIm, LogToConsole) # 11/06/21
            PixID, PixIDTypeAsStr, UniqueValsPostBlur,\
            F2Sinds = GetImageInfo(BlurTxLabIm, LogToConsole) # 11/06/21
            
            #""" The non-binary pixel array following transformation: """
            #PixArr_NB, F2Sinds = Image2PixArr(TxLabIm)
            
            if ApplyPostTxBin:
                """Find suitable threshold value that approximately preserves
                the number of pre-transformed truth values: """
                #Thresh = FindSuitableThresh(BinaryIm=LabIm, 
                #                            NonBinaryIm=TxLabIm,
                #                            LogToConsole=LogToConsole) # 11/06/21
                Thresh = FindSuitableThresh(BinaryIm=TxLabIm, 
                                            NonBinaryIm=BlurTxLabIm,
                                            LogToConsole=LogToConsole) # 11/06/21
                
                """ Binary threshold BlurTxLabIm. """
                #TxLabIm = BinaryThresholdImage(Im=TxLabIm, 
                #                               #Thresh=ThreshPostTx)
                #                               Thresh=Thresh)
                BinBlurTxLabIm = BinaryThresholdImage(Im=BlurTxLabIm, 
                                                      #Thresh=ThreshPostTx)
                                                      Thresh=Thresh)
                
                if LogToConsole:
                    print('\n   Image info for BinBlurTxLabIm after binary',
                          #f'thresholding at {ThreshPostTx}:')
                          f'thresholding at {Thresh}:')
                
                #PixID, PixIDTypeAsStr, UniqueVals,\
                #F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(BinBlurTxLabIm, LogToConsole)
                
                # 11/06/21:
                TxLabIm = BinBlurTxLabIm
            else:
                # 11/06/21:
                TxLabIm = BlurTxLabIm
        
        """ Ensure TxLabIm is a 32-bit unsigned integer (PixID = 5). """
        if PixID != 5: 
            if LogToConsole:
                print(f'\n   TxLabIm has PixelID = {PixID} '\
                      f'({PixIDTypeAsStr})).')
            
            """ Convert TxLabIm from float to 32-bit unsigned integer. """
            TxLabIm = ChangeImagePixelType(Image=TxLabIm, 
                                           NewPixelType='UInt32')
            
            if LogToConsole:
                print('\n   Image info for TxLabIm after converting to',
                      '32-bit unsigned int:')
                
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
        if ApplyPostTxBlur and ApplyPostTxBin:
            """ Store the threshold used as metadata: """
            TxLabIm.SetMetaData("PostTxThreshUsed", f"{Thresh}")
            
            if LogToConsole:
                print(f'\nThe key "PostTxThreshUsed" with value "{Thresh}"',
                      'was added to the metadata of TxLabIm.')
            
        TxLabImByRoi.append(TxLabIm)
        
        """ Convert TxLabIm to a pixel array. """
        TxPixArr, TxF2Sinds = Image2PixArr(TxLabIm)
        
        TxPixArrByRoi.append(TxPixArr)
        TxF2SindsByRoi.append(TxF2Sinds)
        
        if LogToConsole:
            print('\n   Info for conversion of TxLabIm to TxPixArr:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
            unique = UniqueItems(Items=TxPixArr, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArr',
                  f'after transforming LabImByRoi[{r}]')
            
            F = TxPixArr.shape[0]
            
            print(f'\nThe segmentation has been registered to {F} frames/',
                  f'contours: {TxF2Sinds}.')
            
        """ Check if there are fewer indices in TxF2Sinds than in 
        F2SindsByRoi[r]. """
        Nin = len(F2SindsByRoi[r])
        Nout = len(TxF2Sinds)
        
        if Nout < Nin and LogToConsole:
            print(f'\nWARNING:  There are {Nin} F2Sinds in the input labelmap',
                  f'but only {Nout} in the transformed labelmap.\n')
            
    if LogToConsole:
        print('-'*120)
        
    return TxLabImByRoi, TxPixArrByRoi, TxF2SindsByRoi