# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:32:13 2020

@author: ctorti
"""



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







def ResampleLabmapImByRoi_OLD(LabmapImByRoi, C2SindsByRoi, SrcImage, TrgImage, 
                          Variance, ThreshLevel, LogToConsole=False, 
                          Interpolation='NearestNeighbor'):
    """
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled.
    
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour.
        
    SrcImage : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of the 
        labelmap images in LabelmapImByRoi.
    
    TrgImage : SimpleITK image
        The Target 3D image whose gridspace the Source labelmap images will be
        resampled to.
    
    Variance : tuple of floats (optional; (1.0, 1.0, 1.0) by default)
        The variance along all dimensions.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if used.  The threshold value will be ThreshLevel*M, where
        M is the maximum intensity value of the image.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Interpolation : string (optional; 'NearestNeighbor' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear' (or 'linear') 
        - 'BSpline' (or 'Bspline' or 'bspline')
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour'; Default value)
        - 'Gaussian' (or 'gaussian')
    
    
    Outputs:
    *******
    
    ResLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of resampled 3D labelmap images.
        
    ResPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of resampled pixel arrays (converted from each
        labelmap in ResLabmapImByRoi).
    
    ResF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in ResPixArrByRoi.

    
    Notes:
    *****
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. F2Sinds = []) even if there
    were contours (i.e. if SrcC2SindsByRoi[r] != []). If this is the case try 
    suggestions from Ziv Yaniv (See Link below):
            
    1. Use the sitkLabelGaussian interpolator.
    
    2. Gaussian blur the segmentation, then resample using a linear 
    interpolator, and then threshold so that values in [a<1<b] are mapped to 1.
    You will need to select appropriate values for a, b.
    
    3. Compute the distance map of the segmentation, then resample using a
    linear interpolator, and then threshold the absolute value of the distance
    map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you need to 
    select an appropriate dist_threshold.


    What has been tried so far:
    
    1. Resulting Gaussian interpolated resampled image had zeros only 
    (the image was converted to float prior to resampling).
    
    2. Works.
    
    
    Comment from Bradley Lowecamp:
        
    The implementation of the sitkLabelGaussian interpolator in the 
    ResampleImageFilter sets the sigma to the pixel spacing of the input image 
    to the resample filter. This is why it does not do a better job in reducing 
    aliasing.
        
    
    Link: 
        
    https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    import numpy as np
    from ConversionTools import Image2PixArr, ConvertImagePixelType
    
    """ Which Workaround to use. """
    #Workaround = 1
    Workaround = 2
    #Workaround = 0 # i.e. no workaround
        
        
    if LogToConsole:
        print(f'\nResults of ResampleLabmapImByRoi():')
        
    ResLabmapImByRoi = []
    ResPixArrByRoi = []
    ResF2SindsByRoi = []
    
    for r in range(len(LabmapImByRoi)):
        """ Resample LabmapImByRoi to the Target image's grid to get the 
        resampled pixel array. """
        ResLabmapIm = ResampleImage(Image=LabmapImByRoi[r], RefImage=TrgImage, 
                                    Interpolation=Interpolation)
        
        ResLabmapImByRoi.append(ResLabmapIm)
        
        """ Convert ResLabelmapIm to a pixel array. """
        ResPixArr, ResCStoSliceInds = Image2PixArr(LabmapIm=ResLabmapIm)
        
        ResPixArrByRoi.append(ResPixArr)
        
        if LogToConsole:
            print(f'   Image info for ResLabmapImByRoi[{r}] after resampling',
                  f'using {Interpolation} interpolation:')
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        ResF2SindsByRoi.append(F2Sinds)
        
        
    """ Check for the presence of any empty labelmaps (i.e. empty frame-to-slice
    indices) where not expected (a list of empty contour-to-slice indices, e.g.
    C2SindsByRoi[1], will be empty if there were no contours of interest for 
    that ROI). """
    
    for r in range(len(ResLabmapImByRoi)):
        if C2SindsByRoi[r] != [] and ResF2SindsByRoi[r] == [] and Workaround:
    
            if LogToConsole:
                print(f'   C2SindsByRoi[{r}] = {C2SindsByRoi[r]}')
                print('   but')
                print(f'   ResF2SindsByRoi[{r}] = {ResF2SindsByRoi[r]}')
                print('   i.e. the resampled image has no non-zero frames',
                      f'(e.g. due to aliasing). Try a workaround.')
            
            """ Convert LabmapImByRoi[r] from 32-bit unsigned integer to float.
            Note:  Result of resampling results in empty labelmap unless
            image is converted to float prior to resampling. """
            LabmapIm = ConvertImagePixelType(Image=LabmapImByRoi[r], 
                                             NewPixelType='sitkFloat32')
            
            if LogToConsole:
                print(f'   Image info for ResLabmapImByRoi[{r}] after'
                      'converting to float:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
        
            if Workaround == 1:
                if LogToConsole:
                    print(f'   Trying workaround #{Workaround}...')
                
                    
                """ Resample LabmapIm to the Target image's grid using a
                Gaussian interpolation. """
                interp = 'Gaussian'
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after'
                          f'resampling using a {interp} interpolation:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                    
                    
                """ Convert back to unsigned int and apply binary threshold. 
                But only proceed if there are non-zero frames in the resampled 
                image. """
                
                if not F2Sinds:
                    Workaround = 2 # try Workaround #2
            
            if Workaround == 2:
                if LogToConsole:
                    print(f'   Trying workaround #{Workaround}...')
                      
                SrcSpacings = SrcImage.GetSpacing()
                TrgSpacings = TrgImage.GetSpacing()
                
                SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, 
                                                                      TrgSpacings))
                
                """ 
                The FWHM of the Gaussian is:
                    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
                    
                Due to the Nyquist-Shannon sampling theorem, the FWHM 
                should be at least 2x the voxel spacings of TrgIm. For 
                symmetry considerations, it makes sense for one Source 
                voxel width to be blurred to three Target voxels widths.
                    
                    sigma = (3*TrgSpacings)/2.355
                """
                
                if LogToConsole:
                    print(f'   SpacingsRatio = {SpacingsRatio}')
                      
                """ Gaussian blur the image. """
                LabmapIm = GaussianBlurImage(Im=LabmapIm, Variance=Variance)
                
                # Use the RecursiveGaussian image filter:
                #LabmapImToCopy = RecursiveGaussianBlurImage(LabmapIm, 
                #                                            Sigma=3,
                #                                            Direction=2)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after',
                          'applying Gaussian blur:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
                
                #if not F2Sinds:
                #    Workaround = 3 # try Workaround #3
                
                if LogToConsole:
                    print(f'   ResLabmapImByRoi[{r}] prior to {interp}',
                          'interpolation:')
                    print(f'      LabmapIm.GetSize() = {LabmapIm.GetSize()}')
                    print(f'      LabmapIm.GetSpacing() = {LabmapIm.GetSpacing()}')
                    print(f'      TrgImage.GetSize() = {TrgImage.GetSize()}')
                    print(f'      TrgImage.GetSpacing() = {TrgImage.GetSpacing()}')
                    
                """ Linearly resample LabmapIm to the Target image's grid. """
                interp = 'Linear'
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after',
                          'linear resampling:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                
                
                """ Binary threshold the image. """
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                ResLabmapIm = BinaryThresholdImage(Im=ResLabmapIm, 
                                                   Thresh=ThreshValue)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after',
                          'binary thresholding at',
                          f'{ThreshLevel}*{max(UniqueVals)} = {ThreshValue}:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
    
        
        
        """ Ensure that ResLabmapIm is a 32-bit unsigned integer. """
        if not PixID == 5:
            if LogToConsole:
                print(f'   ResLabmapImByRoi[{r}] has PixelID = {PixID}',
                      f'({PixIDTypeAsStr})). Converting to unsigned 32-bit',
                      'integer (sitkUInt32)..')
            
            """ Convert ResLabmapIm from float to 32-bit unsigned integer. """
            ResLabmapIm = ConvertImagePixelType(Image=ResLabmapIm, 
                                                NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print(f'   Image info for ResLabmapImByRoi[{r}] after',
                      'converting to 32-bit unsigned int:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        """ Convert ResLabmapIm to a pixel array. """
        ResPixArr, ResF2Sinds = Image2PixArr(LabmapIm=ResLabmapIm)
        
        if LogToConsole:
            print(f'  After converting ResLabmapImByRoi[{r}] to a pixel array:')
            print(f'   ResPixArr.shape = {ResPixArr.shape}')
            print(f'   ResF2Sinds = {ResF2Sinds}')
            
            #PixID, PixIDTypeAsStr, UniqueVals,\
            #F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
            
            #for i in range(ResPixArr.shape[0]):
            #    maxval = np.amax(ResPixArr[i])
            #    minval = np.amin(ResPixArr[i])
            #    
            #    print(f'      Frame {i} has max = {maxval}, min = {minval}')
                
        """ Replace ResLabmapImByRoi[r], ResPixArrByRoi[r] and 
        ResF2SindsByRoi[r] with ResLabelmapIm, ResPixArr and ResC2Sinds. """
        ResLabmapImByRoi[r] = ResLabmapIm
        ResPixArrByRoi[r] = ResPixArr
        #ResF2SindsByRoi[r] = F2Sinds
        ResF2SindsByRoi[r] = ResF2Sinds
        
    return ResLabmapImByRoi, ResPixArrByRoi, ResF2SindsByRoi






def ResampleLabmapImByRoi_OLD2(LabmapImByRoi, F2SindsByRoi, SrcImage, TrgImage, 
                          Interpolation, Variance, ThreshLevel, #PreBlur=False, 
                          LogToConsole=False):
    """
    Resample a list 3D SimpleITK images representing binary labelmaps using
    either a NearestNeighbor or LabelGaussian interpolator, or by Gaussian
    blurring the labelmap before linearly resampling and binary thresholding.  
    If resampling using a NearestNeighbor or LabelGaussian interpolation,
    aliasing effects can result in a empty resampled labelmap.  If so the input 
    labelmap will be Gaussian blurred, linearly resampled and binary 
    thresholded.  See Notes for more details. 
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled.
    
    F2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each frame) of slice numbers that 
        correspond to each frame in the labelmap images.
        
    SrcImage : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of the 
        labelmap images in LabelmapImByRoi.
    
    TrgImage : SimpleITK image
        The Target 3D image whose gridspace the Source labelmap images will be
        resampled to.
    
    Interpolation : string
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding)
        
    Variance : tuple of floats (optional; (1.0, 1.0, 1.0) by default)
        A tuple (for each dimension) of the variance to be applied if a 
        Gaussian blurring + linearly resampling approach is taken.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if a Gaussian blurring + linearly resampling approach is
        taken.  The threshold value will be ThreshLevel*M, where M is the 
        maximum intensity value of the image.
    
    PreBlur : boolean (optional; False by default)
        If True, the labelmap images will be blurred with a Gaussian filter
        prior to resampling.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    
    
    Outputs:
    *******
    
    ResLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of resampled 3D labelmap images.
        
    ResPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of resampled pixel arrays (converted from each
        labelmap in ResLabmapImByRoi).
    
    ResF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in ResPixArrByRoi.

    
    Notes:
    *****
    
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (labelmap) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled labelmap,
    and hence an empty list of frame-to-slice indices, F2Sinds = []) even if 
    there were contours (i.e. if SrcC2SindsByRoi[r] != []). If this is the case 
    try suggestions from Ziv Yaniv (See Link below):
            
    1. Use the sitkLabelGaussian interpolator.
    
    2. Gaussian blur the segmentation, then resample using a linear 
    interpolator, and then threshold so that values in [a<1<b] are mapped to 1.
    You will need to select appropriate values for a, b.
    
    3. Compute the distance map of the segmentation, then resample using a
    linear interpolator, and then threshold the absolute value of the distance
    map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you need to 
    select an appropriate dist_threshold.


    What has been tried so far:
    
    1. Resulting Gaussian interpolated resampled image had zeros only 
    (the image was converted to float prior to resampling).
    
    2. Works.
    
    
    Comment from Bradley Lowecamp:
        
    The implementation of the sitkLabelGaussian interpolator in the 
    ResampleImageFilter sets the sigma to the pixel spacing of the input image 
    to the resample filter. This is why it does not do a better job in reducing 
    aliasing.
        
    
    Link: 
        
    https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    #import numpy as np
    from ConversionTools import Image2PixArr, ConvertImagePixelType
    
    if not 'eighb' in Interpolation or not 'aussian' in Interpolation or not 'inear' in Interpolation:
        msg = f'The chosen interpolation, {Interpolation}, is not one of the '\
              + 'accepted inputs: NearestNeighbor, LabelGaussian, '\
              + 'BlurThenLinear.'
        
        raise Exception(msg)
        
    
    """ Which Workaround to use. """
    #Workaround = 1
    #Workaround = 2
    #Workaround = 0 # i.e. no workaround
    
    """ It only makes sense to try Workaround 1 if Interpolation was not
    LabelGaussian (since Workaround 1 involves trying a LabelGaussian if a
    NearestNeighbor interpolation was first attempted).  Likewise, if the 
    chosen interpolation is 'BlurThenLinear' then need to use Workaround 2. """
    if Workaround == 1 and 'aussian' in Interpolation or 'inear' in Interpolation:
        Workaround = 2
        
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of ResampleLabmapImByRoi():')
        
    ResLabmapImByRoi = []
    ResPixArrByRoi = []
    ResF2SindsByRoi = []
    
    for r in range(len(LabmapImByRoi)):
        if LogToConsole:
            print(f'   Resampling of LabmapImByRoi[{r}]...')
            
        LabmapIm = LabmapImByRoi[r]
        F2Sinds = F2SindsByRoi[r]
        
        #if PreBlur:
        #    """ Gaussian blur the image. """
        #    LabmapIm = GaussianBlurImage(Im=LabmapIm, Variance=Variance)
        #    
        #    if LogToConsole:
        #        print('\n   Image info for LabmapIm after applying Gaussian',
        #              'pre-blur:')
        #        
        #        PixID, PixIDTypeAsStr, UniqueVals,\
        #        F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
        #        
        #    """ THIS ISN'T COMPLETE!  Can't apply nearest neighbor interpolated
        #    resampling (below) after Gaussian blurring without first converting
        #    to binary image!!! """
                    
        
        if not 'inear' in Interpolation:
            """ Attempt to resample LabmapIm to the Target image's grid using 
            the chosen (non-linear) interpolation. """
            ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                        Interpolation=Interpolation)
            
            """ Convert ResLabmapIm to a pixel array. """
            ResPixArr, ResF2Sinds = Image2PixArr(LabmapIm=ResLabmapIm)
            
            if LogToConsole:
                print(f'\n   r = {r}')
                print(f'\n   Image info for ResLabmapImByRoi[{r}] after',
                      f'resampling using {Interpolation} interpolation:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            ResF2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
            
            #print(f'\n   ResF2Sinds = {ResF2Sinds}')
            
            """ Is ResF2Sinds empty and not expected to be? 
            Note:
                If F2Sinds isn't empty, ResF2Sinds shouldn't be empty either.
            
                F2Sinds will be empty if there were no segmentations/contours 
                of interest for the r^th ROI. In this case an empty ResF2Sinds 
                is acceptable. """
        """
        else Skip to the chosen work-around "BlurThenLinear".
            #ResF2Sinds = []
        """
            
        
        #if True:
        if 'inear' in Interpolation or F2Sinds != [] and Workaround and ResF2Sinds == [] or 'inear' in Interpolation:
            WorkaroundWasDone = True
            
            if LogToConsole:
                if ResF2Sinds == []:
                    print(f'   F2Sinds = {F2Sinds} \nbut ResF2Sinds =',
                          f'{ResF2Sinds}. Try a workaround.')
                else:
                    print(f'   ResF2Sinds = {ResF2Sinds} \nbut chosen ',
                          f'interpolation is {Interpolation}.')
            
            """ Convert LabmapIm from 32-bit unsigned integer to float.
            Note:  Result of resampling results in empty labelmap unless
            image is converted to float prior to resampling. """
            LabmapIm = ConvertImagePixelType(Image=LabmapIm, 
                                             NewPixelType='sitkFloat32')
            
            if LogToConsole:
                print(f'\n   Image info for LabmapImByRoi[{r}] after',
                      'converting to float:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
        
            
            if Workaround == 1:
                if LogToConsole:
                    print('\n   Trying sitkLabelGaussian interpolator...')
                
                    
                """ Resample LabmapIm to the Target image's grid using a
                LabelGaussian interpolation. """
                
                interp = 'LabelGaussian'
                
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print('\n   Image info for ResLabmapIm after resampling'
                          f'using a {interp} interpolation:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                    
                    
                """ Convert back to unsigned int and apply binary threshold. 
                But only proceed if there are non-zero frames in the resampled 
                image. """
                
                if not F2Sinds:
                    Workaround = 2 # try Workaround #2
            
            if Workaround == 2:
                if LogToConsole:
                    print('\n   Trying Gaussian blur + resampling with',
                          'linear interpolator + binary thresholding...')
                      
                SrcSpacings = SrcImage.GetSpacing()
                TrgSpacings = TrgImage.GetSpacing()
                
                SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, 
                                                                      TrgSpacings))
                
                """ 
                The FWHM of the Gaussian is:
                    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
                    
                Due to the Nyquist-Shannon sampling theorem, the FWHM should be
                at least 2x the voxel spacings of TrgIm. For symmetry it makes
                sense for one Source voxel width to be blurred to three Target 
                voxels widths:
                    
                    sigma = (3*TrgSpacings)/2.355
                """
                
                if LogToConsole:
                    print(f'\n   SpacingsRatio = {SpacingsRatio}')
                      
                """ Gaussian blur the image. """
                LabmapIm = GaussianBlurImage(Im=LabmapIm, Variance=Variance)
                
                # Use the RecursiveGaussian image filter:
                #ResLabmapIm = RecursiveGaussianBlurImage(ResLabmapIm, 
                #                                         Sigma=3,
                #                                         Direction=2)
                
                if LogToConsole:
                    print(f'\n   Image info for LabmapImByRoi[{r}] after',
                          'applying Gaussian blur:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
                
                #if not F2Sinds:
                #    Workaround = 3 # try Workaround #3
                
                if LogToConsole:
                    print(f'\n   LabmapImByRoi[{r}] prior to resampling:')
                    print(f'      LabmapIm.GetSize() = {LabmapIm.GetSize()}')
                    print(f'      LabmapIm.GetSpacing() = {LabmapIm.GetSpacing()}')
                    print(f'      TrgImage.GetSize() = {TrgImage.GetSize()}')
                    print(f'      TrgImage.GetSpacing() = {TrgImage.GetSpacing()}')
                    
                """ Linearly resample LabmapIm to the Target image's grid. """
                interp = 'Linear'
                
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print(f'\n   Image info for ResLabmapImByRoi[{r}] after',
                          f'{interp} resampling:')
                    
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                
                
                """ Binary threshold the image. """
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                ResLabmapIm = BinaryThresholdImage(Im=ResLabmapIm, 
                                                   Thresh=ThreshValue)
                
                if LogToConsole:
                    print(f'\n   Image info for ResLabmapImByRoi[{r}] after',
                          f'binary thresholding at {ThreshValue}:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        else:
            WorkaroundWasDone = False
    
        
        
        """ Ensure that ResLabmapIm is a 32-bit unsigned integer. """
        if PixID == 5:
            ConversionWasDone = False
        else:
            ConversionWasDone = True
            
            if LogToConsole:
                print(f'\n   ResLabmapImByRoi[{r}] has PixelID = {PixID}',
                      f'({PixIDTypeAsStr})).')
            
            """ Convert ResLabmapIm from float to 32-bit unsigned integer. """
            ResLabmapIm = ConvertImagePixelType(Image=ResLabmapIm, 
                                                NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print(f'\n   Image info for ResLabmapImByRoi[{r}] after',
                      'converting to 32-bit unsigned int:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        
        
        if WorkaroundWasDone or ConversionWasDone:
            """ Convert ResLabmapIm to a pixel array. """
            ResPixArr, ResF2Sinds = Image2PixArr(LabmapIm=ResLabmapIm)
            
            if LogToConsole:
                print(f'\n  After converting ResLabmapIm[{r}] to a pixel array:')
                print(f'   ResPixArr.shape = {ResPixArr.shape}')
                print(f'   ResF2Sinds = {ResF2Sinds}')
                
                #PixID, PixIDTypeAsStr, UniqueVals,\
                #F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                
                #for i in range(ResPixArr.shape[0]):
                #    maxval = np.amax(ResPixArr[i])
                #    minval = np.amin(ResPixArr[i])
                #    
                #    print(f'      Frame {i} has max = {maxval}, min = {minval}')
                
        
        ResLabmapImByRoi.append(ResLabmapIm)
        ResPixArrByRoi.append(ResPixArr)
        ResF2SindsByRoi.append(ResF2Sinds)
        
        #if LogToConsole:
        #    print('\n')
    
    if LogToConsole:
            print('-'*120)
        
    return ResLabmapImByRoi, ResPixArrByRoi, ResF2SindsByRoi






def RegisterImages_OLD(FixIm, MovIm, Tx='affine', LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    RegImFilt : SimpleITK image filter
        Elastix image transformation filter used to transform MovIm to FixIm.
    """
    
    import SimpleITK as sitk
    import time
    
    if not Tx in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
              + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        
        raise Exception(msg)
        
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Initiate RegImFilt: """
    RegImFilt = sitk.ElastixImageFilter()
    #RegImFilt.LogToConsoleOn() # <-- no output in Jupyter
    RegImFilt.LogToConsoleOff() 
    RegImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    """ Define the fixed and moving images: """
    RegImFilt.SetFixedImage(FixIm)
    RegImFilt.SetMovingImage(MovIm)
    
    """ Get the default parameter map template for the chosen transformation: """
    #RegImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    #ParamMap = sitk.GetDefaultParameterMap('affine')
    ParamMap = sitk.GetDefaultParameterMap(Tx)
    
    """ Re-assign some parameters: """
    ParamMap['AutomaticTransformInitialization'] = ['true']
    ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ParamMap['WriteIterationInfo'] = ['true']
    ParamMap['MaximumNumberOfIterations'] = ['512']
    ParamMap['UseDirectionCosines'] = ['true']
    """ 29/05: Trying this instead of trying to change it for Transformix """
    #ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    """ 13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70 """
    #ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    
    # Print the parameters:
    #for keys,values in ElastixParamMap.items():
    #    print(keys, '=', values)
        
    """ Set the parameter map: """
    RegImFilt.SetParameterMap(ParamMap)
    
    if LogToConsole:
        print('\nPerforming registration...')
        
    #RegImFilt.AddCommand(sitk.sitkProgressEvent, 
    #                     lambda: print("\rProgress: {0:03.1f}%...".format(100*RegImFilt.GetProgress()),end=''))
    # 'ElastixImageFilter' object has no attribute 'AddCommand'
    
    """ Register the 3D images: """
    RegImFilt.Execute()
    
    """ Get the registered image: """
    RegIm = RegImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'\n*Took {Dtime} s to register the 3D image stacks.')
        #print('\n', RegImFilt.GetTransform()) # 'ElastixImageFilter' object has 
        ## no attribute 'GetTransform'
        #print('\n', RegImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
        # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
        # at 0x00000220072E8B70> >
    
    
    return RegIm, RegImFilt
    #return RegIm, RegImFilt, ParamMap





def RegisterImages_OLD2(FixIm, MovIm, Tx='affine', LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    RegImFilt : SimpleITK image filter
        Elastix image transformation filter used to transform MovIm to FixIm.
        
        
    Notes:
    *****
    
    Even if the chosen Tx is Bspline, an affine registration will be applied
    first, as advised in the elastix manual (see section 5.3.5 Transform in
    v5.0.0 of the elastix manual).
    """
    
    import SimpleITK as sitk
    import time
    from copy import deepcopy
    
    if not Tx in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
              + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        
        raise Exception(msg)
        
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Initiate RegImFilt: """
    RegImFilt = sitk.ElastixImageFilter()
    RegImFilt.LogToConsoleOn() # <-- no output in Jupyter
    #RegImFilt.LogToConsoleOff() 
    RegImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    """ Define the fixed and moving images: """
    RegImFilt.SetFixedImage(FixIm)
    RegImFilt.SetMovingImage(MovIm)
    
    """ Get the default parameter map template for the chosen transformation: """
    #RegImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    #ParamMap = sitk.GetDefaultParameterMap('affine')
    #ParamMap = sitk.GetDefaultParameterMap(Tx)
    
    if Tx == 'bspline':
        """ First register using an affine transformation: """
        ParamMap = sitk.GetDefaultParameterMap('affine')
    else:
        ParamMap = sitk.GetDefaultParameterMap(Tx)
        
    
    """ Re-assign some parameters: """
    ParamMap['AutomaticTransformInitialization'] = ['true']
    ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
    ParamMap['WriteIterationInfo'] = ['true']
    ParamMap['MaximumNumberOfIterations'] = ['512']
    ParamMap['UseDirectionCosines'] = ['true']
    """ 29/05: Trying this instead of trying to change it for Transformix """
    #ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    """ 13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70 """
    #ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    
    # Print the parameters:
    #for keys,values in ElastixParamMap.items():
    #    print(keys, '=', values)
        
    """ Set the parameter map: """
    RegImFilt.SetParameterMap(ParamMap)
    
    if True:#LogToConsole:
        if Tx == 'bspline':
            ActualTx = 'affine'
        else:
            ActualTx = deepcopy(Tx)
        print(f'\nPerforming registration using {ActualTx} transformation...')
        
    #RegImFilt.AddCommand(sitk.sitkProgressEvent, 
    #                     lambda: print("\rProgress: {0:03.1f}%...".format(100*RegImFilt.GetProgress()),end=''))
    # 'ElastixImageFilter' object has no attribute 'AddCommand'
    
    """ Register the 3D images: """
    RegImFilt.Execute()
    
    """ Get the registered image: """
    RegIm = RegImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n*Took {Dtime} s to register the 3D image stacks using',
              f'{ActualTx} transform.')
        
    
    if Tx == 'bspline':
        """ Re-Initiate RegImFilt: """
        RegImFilt = sitk.ElastixImageFilter()
        RegImFilt.LogToConsoleOn() # <-- no output in Jupyter
        #RegImFilt.LogToConsoleOff() 
        RegImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
        """ Re-define the fixed and moving images: """
        RegImFilt.SetFixedImage(FixIm)
        RegImFilt.SetMovingImage(RegIm)
    
        """ Initialise a new parameter map: """
        #ParamMap = sitk.ParameterMap()
        ParamMap = sitk.GetDefaultParameterMap(Tx)
        
        """ 22/02/21 Try parameters in Appendix A in "Multimodal image registration
        by edge attraction and regularization using a B-spline grid" by Klein et al
        2011. """
        #ParamMap['MaximumNumberOfIterations'] = ['250']
        ParamMap['MaximumNumberOfIterations'] = ['512']
        #ParamMap['FinalGridSpacingInPhysicalUnits'] = ['20']
        ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
        ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent']
        ParamMap['Transform'] = ['BSplineTransform']
        ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        ParamMap['NumberOfResolutions'] = ['4']
        ParamMap['Metric'] = ['AdvancedMattesMutualInformation']
        ParamMap['NumberOfHistogramBins'] = ['32']
        ParamMap['NumberOfSpatialSamples'] = ['2048']
        ParamMap['NewSamplesEveryIteration'] = ['true']
        ParamMap['ImageSampler'] = ['RandomCoordinate']
        ParamMap['Interpolator'] = ['BSplineInterpolator']
        ParamMap['BSplineInterpolatorOrder'] = ['1']
        ParamMap['FinalBSplineInterpolatorOrder'] = ['3']
        ParamMap['Resampler'] = ['DefaultResampler']
        ParamMap['DefaultPixelValue'] = ['0']
        ParamMap['FixedInternalImagePixelType'] = ['float']
        ParamMap['UseDirectionCosines'] = ['true']
        ParamMap['MovingInternalImagePixelType'] = ['float']
        ParamMap['HowToCombineTransforms'] = ['Compose']
    
    
        # Print the parameters:
        #for keys,values in ElastixParamMap.items():
        #    print(keys, '=', values)
            
        """ Set the parameter map: """
        RegImFilt.SetParameterMap(ParamMap)
        
        if True:#LogToConsole:
            print(f'\nPerforming registration using {Tx} transformation...')
            
        #RegImFilt.AddCommand(sitk.sitkProgressEvent, 
        #                     lambda: print("\rProgress: {0:03.1f}%...".format(100*RegImFilt.GetProgress()),end=''))
        # 'ElastixImageFilter' object has no attribute 'AddCommand'
        
        """ Register the 3D images: """
        RegImFilt.Execute()
        
        """ Get the registered image: """
        RegIm = RegImFilt.GetResultImage()
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        if True:#LogToConsole:
            print(f'\n*Took {Dtime} s to register the 3D image stacks using',
                  f'{Tx} transform.')
            #print('\n', RegImFilt.GetTransform()) # 'ElastixImageFilter' object has 
            ## no attribute 'GetTransform'
            #print('\n', RegImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
            # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
            # at 0x00000220072E8B70> >
    
    
    return RegIm, RegImFilt
    #return RegIm, RegImFilt, ParamMap





def TransformLabmapImByRoi_OLD(LabmapImByRoi, F2SindsByRoi, FixImage, MovImage, 
                           Tx='affine', LogToConsole=False):
    """
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled. The 
        labelmaps share the same image attributes as MovImage.
    
    F2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the frame-to-slice 
        indices of each labelmap image in LabmapImByRoi.
        
    FixImage : SimpleITK image
        The 3D image that MovImage will be registered to (e.g. SrcImage).
    
    MovImage : SimpleITK image
        The 3D image that will be registered to FixImage (e.g. TrgImage).
    
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    TxLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of the 3D labelmap images in LabmapImByRoi
        transformed to the FixImage grid space.
        
    TxPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of pixel arrays (converted from each labelmap 
        in TxLabmapImByRoi).
    
    TxF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in TxPixArrByRoi.
    """
    
    import time
    from ConversionTools import Image2PixArr
    from GeneralTools import UniqueItems
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of TransformLabmapImByRoi():')
        
    if not Tx in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
              + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        
        raise Exception(msg)
        
    """ Start timing. """
    times = []
    times.append(time.time())
    
    print('\n*Performing registration...')
    
    """ Register the images. """
    RegIm, RegImFilt = RegisterImages(FixIm=FixImage, MovIm=MovImage, 
                                      Tx=Tx, LogToConsole=False)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\n*Took {Dtime} s to register the 3D image stacks.')
    
        
    """ Transform LabmapImByRoi to the Target image's grid to get the resampled
    pixel array. """
    TxLabmapImByRoi = []
    TxPixArrByRoi = []
    TxF2SindsByRoi = [] 
    
    for r in range(len(LabmapImByRoi)):
        """ Transform SrcLabmapIm using RegImFilt. """
        TxLabmapIm = TransformImage(Im=LabmapImByRoi[r], RegImFilt=RegImFilt,
                                    Interpolation='Nearestneighbor')
        
        TxLabmapImByRoi.append(TxLabmapIm)
        
        """ Convert TxLabmapIm to a pixel array. """
        TxPixArr, TxF2Sinds = Image2PixArr(LabmapIm=TxLabmapIm)
        
        TxPixArrByRoi.append(TxPixArr)
        TxF2SindsByRoi.append(TxF2Sinds)
            
        #PixID, PixIDTypeAsStr, UniqueVals,\
        #F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        if LogToConsole:
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
            unique = UniqueItems(Items=TxPixArr, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArr',
                  f'after transforming LabmapImByRoi[{r}]')
            
            F = TxPixArr.shape[0]
            
            print(f'\nThe segmentation has been registered to {F} frames/',
                  f'contours: {TxF2Sinds}.')
            
        """ Check if there are fewer indices in TxF2Sinds than in 
        F2SindsByRoi[r]. """
        Nin = len(F2SindsByRoi[r])
        Nout = len(TxF2Sinds)
        
        if Nout < Nin:
            print(f'\nWARNING:  There are {Nin} F2Sinds in the input labelmap',
                  f'but only {Nout} in the transformed labelmap.')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\n*Took {Dtime} s to transform the 3D labelmap images using the',
          'image registration transformation.')
    
    if LogToConsole:
        print('-'*120)
        
    return TxLabmapImByRoi, TxPixArrByRoi, TxF2SindsByRoi







def TransformLabmapImByRoi_OLD(LabmapImByRoi, F2SindsByRoi, FixImage, MovImage, 
                           Tx='bspline', ThreshLevel=0.75, LogToConsole=False):
    """
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled. The 
        labelmaps share the same image attributes as MovImage.
    
    F2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the frame-to-slice 
        indices of each labelmap image in LabmapImByRoi.
        
    FixImage : SimpleITK image
        The 3D image that MovImage will be registered to (e.g. SrcImage).
    
    MovImage : SimpleITK image
        The 3D image that will be registered to FixImage (e.g. TrgImage).
    
    Tx : string (optional; 'bspline' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if a Gaussian blurring + linearly resampling approach is
        taken.  The threshold value will be ThreshLevel*M, where M is the 
        maximum intensity value of the image.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    TxLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of the 3D labelmap images in LabmapImByRoi
        transformed to the FixImage grid space.
        
    TxPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of pixel arrays (converted from each labelmap 
        in TxLabmapImByRoi).
    
    TxF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in TxPixArrByRoi.
    """
    
    import time
    from ConversionTools import Image2PixArr, ConvertImagePixelType
    from GeneralTools import UniqueItems
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of TransformLabmapImByRoi():')
        
    if not Tx in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
              + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        
        raise Exception(msg)
    
            
    """ Start timing. """
    times = []
    times.append(time.time())
    
    print('\n*Performing registration...')
    
    """ Register the images. """
    RegIm, RegImFilt = RegisterImages(FixIm=FixImage, MovIm=MovImage, 
                                      Tx=Tx, LogToConsole=False)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\n*Took {Dtime} s to register the 3D image stacks.')
    
        
    """ Transform LabmapImByRoi to the Target image's grid to get the resampled
    pixel array. """
    TxLabmapImByRoi = []
    TxPixArrByRoi = []
    TxF2SindsByRoi = [] 
    
    for r in range(len(LabmapImByRoi)):
        LabmapIm = LabmapImByRoi[r]
        
        if LogToConsole:
            print(f'\n   Image info for LabmapIm prior to transformation:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
            
            
        """ Transform LabmapIm using RegImFilt. """
        TxLabmapIm = TransformImage(Im=LabmapIm, RegImFilt=RegImFilt,
                                    Interpolation='Nearestneighbor')
        
        
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'\n*Took {Dtime} s to transform the {r}^th labelmap image',
              'using the image registration transformation.')
        
        if LogToConsole:
            print(f'\n   Image info for LabmapIm after transformation:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        
        """ Gaussian blur TxLabmapIm (new 14/02/21) to reduce noise that makes
        selection of an appropriate threshold for binarisation (below) very
        difficult. """
        TxLabmapIm = GaussianBlurImage(TxLabmapIm)
            
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'\n*Took {Dtime} s to Gaussian blur the {r}^th labelmap image.\n')
        
        if LogToConsole:
            print(f'\n   Image info for LabmapIm after Gaussian blur:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        
        """ Binary threshold TxLabmapIm. """
        ThreshValue = ThreshLevel*max(UniqueVals)
        ThreshValue = ThreshLevel
        
        TxLabmapIm = BinaryThresholdImage(Im=TxLabmapIm, 
                                          Thresh=ThreshValue)
        
        if LogToConsole:
            print('\n   Image info for TxLabmapIm after binary thresholding',
                  f'at {ThreshValue}:')
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        """ Ensure TxLabmapIm is a 32-bit unsigned integer (PixID = 5). """
        if PixID != 5: 
            if LogToConsole:
                print(f'\n   TxLabmapIm has PixelID = {PixID}'\
                      f'({PixIDTypeAsStr})).')
            
            """ Convert TxLabmapIm from float to 32-bit unsigned integer. """
            TxLabmapIm = ConvertImagePixelType(Image=TxLabmapIm, 
                                               NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print('\n   Image info for TxLabmapIm after converting to',
                      '32-bit unsigned int:')
                
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        
        TxLabmapImByRoi.append(TxLabmapIm)
        
    
        """ Convert TxLabmapIm to a pixel array. """
        TxPixArr, TxF2Sinds = Image2PixArr(TxLabmapIm)
        
        TxPixArrByRoi.append(TxPixArr)
        TxF2SindsByRoi.append(TxF2Sinds)
            
        #PixID, PixIDTypeAsStr, UniqueVals,\
        #F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        if LogToConsole:
            print('\n   Info for conversion of TxLabmapIm to TxPixArr:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
            unique = UniqueItems(Items=TxPixArr, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArr',
                  f'after transforming LabmapImByRoi[{r}]')
            
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
        
    return TxLabmapImByRoi, TxPixArrByRoi, TxF2SindsByRoi








def SitkBsplineRegWithFiducials_OLD(FixIm, MovIm, NumControlPts=8, 
                                FixFidsFpath='fixed_fiducials.txt', 
                                MovFidsFpath='moving_fiducials.txt', 
                                NumIters=100, LogToConsole=False):
    """ Best result from 14/05/21 achieved in 
    20210423_Bspline_registrations.ipynb
    
    NOTE:
        This doesn't work because the initialisation using a deformable
        LandmarkBasedTransformInitializer produces a bizarre result.
    """
    import SimpleITK as sitk
    import time
    import winsound
    
    # Print simple metric values to console or plot of metric value v iteration?
    PlotToConsole = False # simple
    PlotToConsole = True
    
    
    OrigFixFids, FixFidType = ImportFiducialsFromTxt(FixFidsFpath)

    OrigMovFids, MovFidType = ImportFiducialsFromTxt(MovFidsFpath)
    
    print('Fixed fiducials:\n', OrigFixFids, '\n')
    print('Moving fiducials:\n', OrigMovFids, '\n')
    
    """ For BSpline Transforms it's not enough for the points to lie within the
    extent of the image - they also must not lie on the boundary:
    https://itk.org/pipermail/insight-users/2012-June/044975.html
    
    So eliminate all fiducial pairs for which one or both fiducials lie on any
    boundary in the corresponding image.
    """
    xtra = 0 # exclude fiducials on any boundary
    xtra = 3 # extend exclusion to 3 pixels from boundary
    xtra = 1 # extend exclusion to 1 pixel from boundary
    
    FixFids, MovFids,\
    inds = GetFiducialsWithinImExtents(FixIm, OrigFixFids, MovIm, OrigMovFids,
                                       xtra)
    
    if len(FixFids) < len(OrigFixFids):
        import GeneralTools
        import importlib
        importlib.reload(GeneralTools)
        from GeneralTools import xor
        
        allinds = list(range(len(OrigFixFids)))
        
        xorinds = xor(allinds, inds)
        
        print('Fixed fiducials that lie within FixIm:\n', FixFids, '\n')
        print('Moving fiducials that lie within MovIm:\n', MovFids, '\n')
        print('Eliminated Fixed fiducials:\n', 
              [OrigFixFids[i] for i in xorinds], '\n')
        print('Eliminated Moving fiducials:\n', 
              [OrigMovFids[i] for i in xorinds], '\n')
    
    if FixFidType == 'index':
        FixFids = [FixIm.TransformIndexToPhysicalPoint(ind) for ind in FixFids]

    if MovFidType == 'index':
        MovFids = [MovIm.TransformIndexToPhysicalPoint(ind) for ind in MovFids]
    
    #print(f'Fixed fiducials: \n {FixFids}\n')
    #print(f'Moving fiducials: \n {MovFids}\n')
    
    """ Flatten the list of points """
    FixPts = [c for p in FixFids for c in p]
    MovPts = [c for p in MovFids for c in p]
    
    
    if False:
        InitialTx = sitk.LandmarkBasedTransformInitializer(transform=sitk.VersorRigid3DTransform(),
                                                           fixedLandmarks=FixPts, 
                                                           movingLandmarks=MovPts)
    
    if False:
        dim = FixIm.GetDimension()
        
        InitialTx = sitk.LandmarkBasedTransformInitializer(transform=sitk.BSplineTransform(dim),
                                                           fixedLandmarks=FixPts, 
                                                           movingLandmarks=MovPts,
                                                           referenceImage=FixIm,
                                                           numberOfControlPoints=NumControlPts)
                                                           #BSplineNumberOfControlPoints=NumControlPts)
        """ Shouldn't be using BSplineTransform unless objective is to elastically
        deform the image using the fiducials! """
        
    if True:
        # https://github.com/SimpleITK/SimpleITK/issues/197
        MeshSize = [NumControlPts] * MovIm.GetDimension()
    
        InitTx = sitk.BSplineTransformInitializer(image1=FixIm,
                                                  transformDomainMeshSize=MeshSize)
        
        LandmarkTx = sitk.LandmarkBasedTransformInitializerFilter()
        LandmarkTx.SetFixedLandmarks(FixPts)
        LandmarkTx.SetMovingLandmarks(MovPts)
        LandmarkTx.SetBSplineNumberOfControlPoints(NumControlPts)
        LandmarkTx.SetReferenceImage(FixIm)
        InitialTx = LandmarkTx.Execute(InitTx)
    
    
    """ 24/05/21: It seems that despite setting the number of control points
    for LandmarkTx, that info isn't preserved during the call
        InitialTx = LandmarkTx.Execute(InitTx)
        
    and InitialTx ends up with a TransformDomainMeshSize of [20, 20, 20]
    (which is odd given that the default number of control points when
    calling LandmarkBasedTransformInitializerFilter() is 4).
    
    So need to re-set the mesh size for InitialTx, but in order to do that,
    need to convert InitialTx from Transform to BSplineTransform.
    """
    InitialTx = sitk.BSplineTransform(InitialTx)
    InitialTx.SetTransformDomainMeshSize(MeshSize)
    
    print(f'type(InitTx) = {type(InitTx)}')
    print(f'type(LandmarkTx) = {type(LandmarkTx)}')
    print(f'type(InitialTx) = {type(InitialTx)}')
    print('InitialTx:\n')
    print(InitialTx, '\n')
    
    #SamplingPercentage = 0.01
    SamplingPercentage = 0.05
    
    LearningRate = 5.0 # Took 4424.5 s (73.7 min)
    LearningRate = 50.0 # Took 5004.0 s (83.4 min) 
    
    """ Start timing: """
    times = []
    times.append(time.time())
        
    ##MeshSize = [8] * MovIm.GetDimension()
    #MeshSize = [NumControlPts] * MovIm.GetDimension()
    #
    #InitialTx = sitk.BSplineTransformInitializer(image1=FixIm,
    #                                             transformDomainMeshSize=MeshSize)
    
    #print("Initial Parameters:")
    #print(InitialITx.GetParameters())
    
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, True)
    
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    """ Change to original method: """
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
                                                      numberOfIterations=NumIters,
                                                      convergenceMinimumValue=1e-4,
                                                      convergenceWindowSize=5)
    
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    RegMethod.SetShrinkFactorsPerLevel([6, 2, 1])
    RegMethod.SetSmoothingSigmasPerLevel([6, 2, 1])
    
    if True:#LogToConsole:
        if PlotToConsole:
            # Connect all of the observers so that we can perform plotting  
            # during registration.
            RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
            RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    print("-------")
    print('InitialTx:')
    print(InitialTx)
    print("-------")
    print("-------")
    print('FinalTx:')
    print(FinalTx)
    print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print(f"Optimizer stop condition: {RegMethod.GetOptimizerStopConditionDescription()}")
    #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(FixIm)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0.0)
    resampler.SetTransform(FinalTx)
    
    RegIm = resampler.Execute(MovIm)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    Dtime_min = round(Dtime/60, 1)
    print(f'*Took {Dtime} s ({Dtime_min} min) to perform image registration.\n')
    
    # Convert FinalTx from sitk Transform to sitk BSplineTransform:
    #FinalTx = sitk.BSplineTransform(FinalTx)
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    # Took 2283 s = 38 min at iteration 18 (all samples)
    # Took 118.6 s = 2 min (1% samples)
    # Took 214 s = 3.6 min (5% samples)
    
    return RegIm, RegMethod, InitialTx, FinalTx





