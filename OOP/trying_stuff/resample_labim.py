def resample_labim(LabIm, F2Sinds, SrcIm, TrgIm, Interp, 
                   PreResVariance=(1,1,1), ApplyPostResBlur=True, 
                   PostResVariance=(1,1,1), LogToConsole=False):
    """
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------  
    LabIm : SimpleITK image
        The 3D label image to be resampled.
    F2Sinds : list of ints
        List (for each frame) of slice numbers that correspond to each frame in 
        LabIm.
    SrcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of LabIm. 
    TrgIm : SimpleITK image
        The Target 3D image whose gridspace the LabIm will be resampled to.
    Interp : str
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    PreResVariance : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling.
    ApplyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    ResLabIm : SimpleITK images
        The resampled 3D label image.
    ResPixArr : Numpy arrays
        The resampled pixel array (converted from ResLabIm).
    ResF2Sinds : list of ints
        The list (for each frame) of the frame-to-slice indices in ResPixArr.
    
    Note
    ----
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (label) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled label,
    and hence an empty list of frame-to-slice indices, F2Sinds = []) even if 
    there were contours (i.e. if SrcC2Sinds != []). If this is the case try
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
    
    Link: https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of resample_labim():')
        print('\n\n', '-'*120)
        
    if not Interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f'The chosen interpolation, {Interp}, is not one of the '\
              + 'accepted inputs: \'NearestNeighbor\', \'LabelGaussian\', or '\
              + '\'BlurThenLinear\'.'
        
        raise Exception(msg)
    
    # Store the interpolation set as metadata:
    LabIm.SetMetaData("ResInterpSet", Interp)
    
    if Interp in ['NearestNeighbor', 'LabelGaussian']:
        if LogToConsole:
            print('Attempting to resample LabIm using the',
                  f'{Interp} interpolator...\n')
        
        """ Attempt to resample LabIm to the Target image's grid using the
        chosen labelmap interpolator. """
        ResLabIm = resample_im(Im=LabIm, RefIm=TrgIm, Interp=Interp)
        
        if LogToConsole:
            print('Image info for ResLabIm after resampling using',
                  f'{Interp} interpolator:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
        if LogToConsole:
            print('')
        
        if ApplyPostResBlur:
            # Gaussian blur the image.
            ResLabIm = gaussian_blur_im(Im=ResLabIm, 
                                        Variance=PostResVariance)
            
            if LogToConsole:
                print('Image info for ResLabIm after Gaussian blurring:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
            if LogToConsole:
                print('')
            
            # Binarise the resampled labelmap if required.
            if len(UniqueVals) != 2 or sum(UniqueVals) != 1:
                """Find suitable threshold value that approximately preserves
                the number of pre-blurred truth values scaled by VolumeRatio: 
                """
                Thresh = find_thresh(BinaryIm=LabIm, 
                                     NonBinaryIm=ResLabIm,
                                     LogToConsole=LogToConsole)
                
                # Binary threshold the image.
                ResLabIm = binarise_im(Im=ResLabIm, 
                                       #Thresh=PostResThresh)
                                       Thresh=Thresh) 
                
                if LogToConsole:
                    print('\nImage info for ResLabIm after binary',
                          #f'thresholding at {PostResThresh}:')
                          f'thresholding at {Thresh}:')
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = get_im_info(ResLabIm, LogToConsole)
                if LogToConsole:
                    print('')
        
        #print(f'\n   ResF2Sinds = {ResF2Sinds}')
        
        """ Is ResF2Sinds empty and not expected to be? If so, try the 
        "BlurThenLinear" approach.
        
        Note:
            If F2Sinds isn't empty, ResF2Sinds shouldn't be empty either.
        
            F2Sinds will be empty if there were no segmentations/contours 
            of interest for the r^th ROI. In this case an empty ResF2Sinds 
            is acceptable. """
        
        if ResF2Sinds == []:
            print(f'There are {len(F2Sinds)} non-empty masks in the input',
                  f'labelmap image but {len(ResF2Sinds)} non-empty frames in',
                  f'the resampled labelmap image using {Interp}.',
                  'Try using the \'BlurThenLinear\' approach.\n')
            
            Interp = 'BlurThenLinear'

    if Interp == 'BlurThenLinear':
        if LogToConsole:
            print('\nResampling LabIm by applying Gaussian blur,',
                  'resampling using a linear interpolator, then binary',
                  'thresholding...')
        
        #""" The original binary pixel array: """
        #PixArr_B, F2Sinds = im_to_pixarr(LabIm)
        
        """ Convert LabIm from 32-bit unsigned integer to float.
        Note:  Result of resampling results in empty labelmap unless
        image is converted to float prior to resampling. """
        LabIm = change_im_dtype(Image=LabIm, NewPixelType='Float32')
        
        if LogToConsole:
            print('\nImage info for LabIm after converting to float:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = get_im_info(LabIm, LogToConsole)
              
        # Gaussian blur the image.
        LabIm = gaussian_blur_im(Im=LabIm, Variance=PreResVariance)
        
        if LogToConsole:
            print('\nImage info for LabIm after applying Gaussian blur:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = get_im_info(LabIm, LogToConsole)
        
        if LogToConsole:
            print('\nLabIm prior to resampling:')
            print(f'   LabIm.GetSize() = {LabIm.GetSize()}')
            print(f'   LabIm.GetSpacing() = {LabIm.GetSpacing()}')
            print(f'   TrgIm.GetSize() = {TrgIm.GetSize()}')
            print(f'   TrgIm.GetSpacing() = {TrgIm.GetSpacing()}')
        
        # Linearly resample LabIm to the Target image's grid.
        ResLabIm = resample_im(Im=LabIm, RefIm=TrgIm, Interp='Linear')
        
        if LogToConsole:
            print('\nImage info for ResLabIm after resampling using a',
                  'linear interpolator:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = get_im_info(ResLabIm, LogToConsole)
        
        if ApplyPostResBlur:
            # Gaussian blur the image.
            ResLabIm = gaussian_blur_im(Im=ResLabIm, 
                                        Variance=PostResVariance)
        
        #""" The non-binary pixel array following blur + linear resampling: """
        #PixArr_NB, F2Sinds = im_to_pixarr(ResLabIm)
        #
        #if LogToConsole:
        #        print(f'PixArr_B.shape = {PixArr_B.shape}')
        #        print(f'PixArr_NB.shape = {PixArr_NB.shape}\n')
            
        """Find suitable threshold value that approximately preserves the
        number of pre-blurred + linear resampled truth values scaled by
        VolumeRatio: """
        Thresh = find_thresh(BinaryIm=LabIm, NonBinaryIm=ResLabIm,
                             LogToConsole=LogToConsole)
        
        # Binary threshold the image. 
        ResLabIm = binarise_im(Im=ResLabIm, 
                               #Thresh=PostResThresh)
                               Thresh=Thresh) 
        
        if LogToConsole:
            print('\nImage info for ResLabIm after binary thresholding',
                  #f'at {PostResThresh}:')
                  f'at {Thresh}:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = get_im_info(ResLabIm, LogToConsole)
    
    # Ensure that ResLabIm is a 32-bit unsigned integer (PixID = 5).
    if PixID != 5: 
        if LogToConsole:
            print(f'\nResLabIm has PixelID = {PixID} ({PixIDTypeAsStr})).')
        
        # Convert ResLabIm from float to 32-bit unsigned integer. 
        ResLabIm = change_im_dtype(Image=ResLabIm, NewPixelType='UInt32')
        
        if LogToConsole:
            print('\nImage info for ResLabIm after converting to 32-bit',
                  'unsigned int:')
        #print(f'\nThe metadata keys are:', ResLabIm.GetMetaDataKeys())
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = get_im_info(ResLabIm, LogToConsole)
    
    # Convert ResLabIm to a pixel array. 
    ResPixArr, ResF2Sinds = im_to_pixarr(ResLabIm)
    
    # Store the interpolation used as metadata (which may be the same or 
    # different from the interpolation set): 
    ResLabIm.SetMetaData("ResInterpUsed", Interp)
    
    if LogToConsole:
        print(f'\nThe key "ResInterpUsed" with value "{Interp}" was',
              'added to the metadata of ResLabIm. \nThe metadata keys are:',
              ResLabIm.GetMetaDataKeys())
    
    if Interp == 'BlurThenLinear':
        """ Store the threshold used as metadata: """
        ResLabIm.SetMetaData("PostResThreshUsed", f"{Thresh}")
        
        if LogToConsole:
            print(f'\nThe key "PostResThreshUsed" with value "{Thresh}" was',
                  'added to the metadata of ResLabIm. \nThe metadata keys:',
                  ResLabIm.GetMetaDataKeys())
            
    if LogToConsole:
        print('\nAfter converting ResLabIm to a pixel array:')
        print(f'ResPixArr.shape = {ResPixArr.shape}')
        print(f'ResF2Sinds = {ResF2Sinds}')
        print('-'*120)
        
    return ResLabIm, ResPixArr, ResF2Sinds