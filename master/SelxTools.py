# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 11:32:44 2021

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
SimpleElastix Functions
******************************************************************************
******************************************************************************
"""


def RegisterImages_OLD(FixIm, MovIm, Tx='bspline', MaxNumOfIters=512,
                   LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Tx : string (optional; 'bspline' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxNumOfIters : integer
        The maximum number of iterations to be used during optimisation.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    RegImFilt : Elastix image filter
        Elastix image registration transformation filter used to register/
        transform MovIm to FixIm.
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
        """ 22/04/2021: I'm not sure this was the correct approach to 
        performing a BSpline registration... """
    else:
        ParamMap = sitk.GetDefaultParameterMap(Tx)
        
    
    """ Re-assign some parameters: """
    ParamMap['AutomaticTransformInitialization'] = ['true']
    ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
    ParamMap['WriteIterationInfo'] = ['true']
    #ParamMap['MaximumNumberOfIterations'] = ['512']
    ParamMap['MaximumNumberOfIterations'] = [str(MaxNumOfIters)]
    ParamMap['UseDirectionCosines'] = ['true']
    
    """ 29/05: Trying this instead of trying to change it for Transformix """
    #ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    """ 13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70 """
    #ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    
    if Tx == 'bspline':
        """ 22/02/21 Try parameters in Appendix A in "Multimodal image 
        registration by edge attraction and regularization using a B-spline 
        grid" by Klein et al 2011; and from the SimpleElastix docs:
        https://simpleelastix.readthedocs.io/ParameterMaps.html """
        
        #ParamMap['MaximumNumberOfIterations'] = ['250']
        #ParamMap['MaximumNumberOfIterations'] = ['512']
        ParamMap['MaximumNumberOfIterations'] = [str(MaxNumOfIters)]
        #ParamMap['FinalGridSpacingInPhysicalUnits'] = ['20']
        #ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
        ParamMap['Registration'] = ['MultiResolutionRegistration']
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
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # typo? (found 22/04/2021)
        ParamMap['FinalBSplineInterpolationOrder'] = ['3'] # 22/04/2021
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
    
    if LogToConsole:
        print(f'Performing image registration using {Tx} transform...\n')
        
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
        print(f'*Took {Dtime} s to perform image registration.\n')
        #print('\n', RegImFilt.GetTransform()) # 'ElastixImageFilter' object has 
        ## no attribute 'GetTransform'
        #print('\n', RegImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
        # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
        # at 0x00000220072E8B70> >
    
    
    return RegIm, RegImFilt
    #return RegIm, RegImFilt, ParamMap
    
    
    
    
    




def RegisterImagesSelx_OLD(FixIm, MovIm, Transform='affine', 
                       MaxNumOfIters='default', InitMethod='centerofgravity',
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       FlipK=False,
                       LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxNumOfIters : string (optional; 'default' by default)
        If 'default', the maximum number of iterations used during optimisation
        will be whatever is set by default in the parameter map of choice.  
        If != 'default' it must be a string representation of an integer.
    
    InitMethod : string (optional; 'centerofgravity' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'centerofgravity' (or 'moments')
            - 'geometricalcenter' (or 'geometricalcenter')
            - 'origin' (only available for affine transform)
            - 'landmarks'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    SelxImFilt : Elastix image filter
        Elastix image registration transformation filter used to register/
        transform MovIm to FixIm.
    
    ParamMapDict : dictionary
        A dictionary containing the parameter map from SelxImFilt.
    
    TxParamMapDict : dictionary
        A dictionary containing the transform parameter map from SelxImFilt.
    
        
    Notes:
    *****
    
    29/05/20: Trying this instead of trying to change it for Transformix
    ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70
    ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    22/02/21 Try parameters in Appendix A in "Multimodal image registration by 
    edge attraction and regularization using a B-spline grid" by Klein et al 
    2011; and from the SimpleElastix docs:
    https://simpleelastix.readthedocs.io/ParameterMaps.html
    """
    
    import SimpleITK as sitk
    import time
    
    if not Transform in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transform (Transform), {Transform}, is not one of '\
              + 'the accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        raise Exception(msg)
    
    if not InitMethod in ['centerofgravity', 'moments', 'geometricalcenter', 
                          'geometry', 'origins', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'centerofgravity'/'moments',"\
              + " 'geometricalcenter'/'geometry', 'origins' or 'landmarks'."
        raise Exception(msg)
        
    if Transform != 'affine' and InitMethod == 'origins':
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              "is only allowed for the affine transform."
        raise Exception(msg)
    
    
    """ Initiate ElastixImageFilter: """
    SelxImFilt = sitk.ElastixImageFilter()
    if LogToConsole:
        SelxImFilt.LogToConsoleOn() # no output in Jupyter
    else:
        SelxImFilt.LogToConsoleOff() 
    SelxImFilt.LogToFileOn() # output's elastix.log ONLY ONCE per kernel in Jupyter
    
    SelxImFilt.SetFixedImage(FixIm)
    SelxImFilt.SetMovingImage(MovIm)
    
    """ Get the default parameter map for the chosen transform: """
    #SelxImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    #ParamMap = sitk.GetDefaultParameterMap('affine')
    #ParamMap = sitk.GetDefaultParameterMap(Transform) # 10/05/21
    if InitMethod == 'landmarks': # 10/05/21
        #ParamMap = sitk.VectorOfParameterMap() # 14/05/21
        #ParamMap.append(sitk.GetDefaultParameterMap('translation')) # 14/05/21
        SelxImFilt.SetParameterMap(sitk.GetDefaultParameterMap('translation')) # 14/05/21
        if Transform == 'bspline':
            #ParamMap.append(sitk.GetDefaultParameterMap('affine')) # 14/05/21
            #ParamMap.append(sitk.GetDefaultParameterMap('bspline')) # 14/05/21
            SelxImFilt.AddParameterMap(sitk.GetDefaultParameterMap('affine')) # 14/05/21
            SelxImFilt.AddParameterMap(sitk.GetDefaultParameterMap('bspline')) # 14/05/21
        else:
            #ParamMap.append(sitk.GetDefaultParameterMap(Transform)) # 14/05/21
            SelxImFilt.AddParameterMap(sitk.GetDefaultParameterMap(Transform)) # 14/05/21
    else:
        if Transform == 'bspline':
            #ParamMap = sitk.VectorOfParameterMap() # 14/05/21
            #ParamMap.append(sitk.GetDefaultParameterMap('affine')) # 14/05/21
            #ParamMap.append(sitk.GetDefaultParameterMap('bspline')) # 14/05/21
            SelxImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine')) # 14/05/21
            SelxImFilt.AddParameterMap(sitk.GetDefaultParameterMap('bspline')) # 14/05/21
        else:
            #ParamMap = sitk.GetDefaultParameterMap(Transform) # 14/05/21
            SelxImFilt.SetParameterMap(sitk.GetDefaultParameterMap(Transform)) # 14/05/21
    ParamMap = SelxImFilt.GetParameterMap() # 14/05/21
    
    """ 18/03/21:
    
    Kasper said that AutomaticTransformInitialization should be true for rigid,
    and false for affine or bspine: 
    
    https://github.com/SuperElastix/elastix/issues/45
         
    - FixedImagePyramid should be FixedImageSmoothingPyramid for best results 
    (unless memory is really and issue).
    - AutomaticTransformInitialization should be true for rigid and false for 
    affine and bspline.
    - 250 number of iterations is also on the low side. If speed is really an 
    issue, 250 is fine, but if you can afford it use at least 512.
    - Depending on the problem you might want to use at least 3 resolutions. 
    Some places you use 2 which is on the low side.
    - ImageSampler should be RandomSparseMask if users use masks. Otherwise it 
    is fine.
    - If you use a BSplineInterplationOrder of 1 you might as well use 
    LinearInterpolator which produces the same results but is faster.

    """
    
    """ Add/modify parameters: """
    if Transform == 'rigid': # 26/04/21
        ParamMap['AutomaticTransformInitialization'] = ['true']
    else:
        ParamMap['AutomaticTransformInitialization'] = ['false']
    #ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter'] # 19/03/21
    #ParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity'] # 19/03/21
    #ParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
    
    if InitMethod in ['centerofgravity', 'moments']:
        ParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
    elif InitMethod in ['geometricalcenter', 'geometry']:
        ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    elif InitMethod == 'origins':
        ParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
        
    
    ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
    #ParamMap['DefaultPixelValue'] = ['0']
    """ 05/05/21: Not sure why BSplineInterpolatorOrder was set to 1. """
    if False:
        ParamMap['BSplineInterpolatorOrder'] = ['1']
    #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
    #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
    #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
    #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
    #if Transform == 'bspline': # 22/04/2021
    #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
    #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
    ParamMap['FixedInternalImagePixelType'] = ['float']
    ParamMap['MovingInternalImagePixelType'] = ['float']
    ParamMap['HowToCombineTransforms'] = ['Compose']
    #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
    #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
    #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
    #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
    #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
    if MaxNumOfIters != 'default':
        ParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
    #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
    #ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
    #                      'TransformBendingEnergyPenalty'] # default for bspline
    if InitMethod == 'landmarks':
        if Transform in ['rigid', 'affine']:
            ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
                                  'CorrespondingPointsEuclideanDistanceMetric']
            #ParamMap['Metric0Weight'] = ['1.0']
            ParamMap['Metric0Weight'] = ['0.3']
            ParamMap['Metric1Weight'] = ['1.0']
            
        elif Transform == 'bspline':
            ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
                                  'TransformBendingEnergyPenalty',
                                  'CorrespondingPointsEuclideanDistanceMetric']
            ParamMap['Metric0Weight'] = ['1.0']
            ParamMap['Metric1Weight'] = ['1.0']
            ParamMap['Metric2Weight'] = ['1.0']
    #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
    ParamMap['NumberOfHistogramBins'] = ['32']
    #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
    #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
    #if Transform == 'bspline': # 22/04/2021
    #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
    #    ParamMap['Metric0Weight'] = ['0.5']
    #    ParamMap['Metric1Weight'] = ['1.0']
    #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
    if InitMethod == 'landmarks':
        """ MultiResolutionRegistration, which is default for affine, only
        allows for 1 metric. """
        ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration']
    #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
    #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
    #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
    #ParamMap['Transform'] = ['AffineTransform'] # default for affine
    #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
    ParamMap['UseDirectionCosines'] = ['true']
    ParamMap['WriteIterationInfo'] = ['true']
    
    
    # Print the parameters:
    #for keys,values in ParamMap.items():
    #    print(keys, '=', values)
        
    """ Set the parameter map: """
    SelxImFilt.SetParameterMap(ParamMap)
    
    if InitMethod == 'landmarks':
        SelxImFilt.SetFixedPointSetFileName(FixFidsFpath)
        SelxImFilt.SetMovingPointSetFileName(MovFidsFpath)
    
    if LogToConsole:
        print(f'Performing image registration using {Transform} transform...\n')
        
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Register the 3D images: """
    SelxImFilt.Execute()
    
    """ Get the registered image: """
    RegIm = SelxImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
        #print('\n', ElxImFilt.GetTransform()) # 'ElastixImageFilter' object has 
        ## no attribute 'GetTransform'
        #print('\n', ElxImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
        # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
        # at 0x00000220072E8B70> >
        
    
    #""" Get the items in parameter map and in the transform parameter map
    #(both from file and from the Elastix image filter) as a dictionary. """ 
    #
    #""" First get the parameter map from SelxImFilt: """
    #ParamMapDict = ParamMapFromSelxImFilt(SelxImFilt)
    #
    #""" Add the transform parameters from SelxImFilt: """
    #TxParamMapDict = TxParamMapFromSelxImFilt(SelxImFilt, ParamMapDict)
    #
    #""" Add the additional transform parameters (and over-write existing ones*)
    #from the file TransformParameters.0.txt.
    #
    #* There are more decimals in the parameters containted from the file than
    #those stored in ElxImFilt. """
    #TxParamMapDict = TxParamMapFromFile(TxParamMapDict)
    
    """ Get the parameter map and transform parameter map: """
    ParamMapDict = GetParamMapFromSelxImFilt(SelxImFilt)
    TxParamMapDict = GetTxParamMapFromSelxImFilt(SelxImFilt)
    
    
    #return RegIm, SelxImFilt
    return RegIm, SelxImFilt, ParamMapDict, TxParamMapDict






def RegisterImagesSelx(FixIm, MovIm, Transform='affine', 
                       MaxNumOfIters='default', InitMethod='centerofgravity',
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       FlipK=False,
                       LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'translation'
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxNumOfIters : string (optional; 'default' by default)
        If 'default', the maximum number of iterations used during optimisation
        will be whatever is set by default in the parameter map of choice.  
        If != 'default' it must be a string representation of an integer.
    
    InitMethod : string (optional; 'centerofgravity' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'centerofgravity' (or 'moments')
            - 'geometricalcenter' (or 'geometricalcenter')
            - 'origin' (only available for affine transform)
            - 'landmarks'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    SelxImFilt : Elastix image filter
        Elastix image registration transformation filter used to register/
        transform MovIm to FixIm.
    
    ParamMapDict : dictionary
        A dictionary containing the parameter map from SelxImFilt.
    
    TxParamMapDict : dictionary
        A dictionary containing the transform parameter map from SelxImFilt.
    
        
    Notes:
    *****
    
    29/05/20: Trying this instead of trying to change it for Transformix
    ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70
    ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    22/02/21 Try parameters in Appendix A in "Multimodal image registration by 
    edge attraction and regularization using a B-spline grid" by Klein et al 
    2011; and from the SimpleElastix docs:
    https://simpleelastix.readthedocs.io/ParameterMaps.html
    
    18/03/21 Kasper said that AutomaticTransformInitialization should be true 
    for rigid, and false for affine or bspine: 
    
    https://github.com/SuperElastix/elastix/issues/45
         
    - FixedImagePyramid should be FixedImageSmoothingPyramid for best results 
    (unless memory is really and issue).
    - AutomaticTransformInitialization should be true for rigid and false for 
    affine and bspline.
    - 250 number of iterations is also on the low side. If speed is really an 
    issue, 250 is fine, but if you can afford it use at least 512.
    - Depending on the problem you might want to use at least 3 resolutions. 
    Some places you use 2 which is on the low side.
    - ImageSampler should be RandomSparseMask if users use masks. Otherwise it 
    is fine.
    - If you use a BSplineInterplationOrder of 1 you might as well use 
    LinearInterpolator which produces the same results but is faster.
    """
    
    import SimpleITK as sitk
    import time
    
    if not Transform in ['translation', 'rigid', 'affine', 'bspline']:
        msg = f'The chosen transform (Transform), {Transform}, is not one of '\
              + 'the accepted inputs: \'translation\', \'rigid\', \'affine\','\
              + ' or \'bspline\'.'
        raise Exception(msg)
    
    if not InitMethod in ['centerofgravity', 'moments', 'geometricalcenter', 
                          'geometry', 'origins', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'centerofgravity'/'moments',"\
              + " 'geometricalcenter'/'geometry', 'origins' or 'landmarks'."
        raise Exception(msg)
        
    if Transform != 'affine' and InitMethod == 'origins':
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              "is only allowed for the affine transform."
        raise Exception(msg)
    
    
    """ Initiate ElastixImageFilter: """
    SelxImFilt = sitk.ElastixImageFilter()
    if LogToConsole:
        SelxImFilt.LogToConsoleOn() # no output in Jupyter
    else:
        SelxImFilt.LogToConsoleOff() 
    SelxImFilt.LogToFileOn() # output's elastix.log ONLY ONCE per kernel in Jupyter
    
    SelxImFilt.SetFixedImage(FixIm)
    SelxImFilt.SetMovingImage(MovIm)
    
    #print(f'InitMethod = {InitMethod}\n')
    
    """ Define the parameter map. Deal with different cases separately. """
    if InitMethod != 'landmarks' and Transform != 'bspline':
        #print('case 1\n')
        
        #print(f'Transform = {Transform}\n')
        
        """ This is a single-itemed rigid parameter map. """
        ParamMap = sitk.GetDefaultParameterMap(Transform)
        
        #print(f"ParamMap['Transform'] = {ParamMap['Transform']}\n")
        
        """ Set the parameters accordingly: """
        if Transform in ['translation', 'rigid']:
            ParamMap['AutomaticTransformInitialization'] = ['true']
        else:
            ParamMap['AutomaticTransformInitialization'] = ['false']
        
        if InitMethod in ['centerofgravity', 'moments']:
            ParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
        elif InitMethod in ['geometricalcenter', 'geometry']:
            ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
        elif InitMethod == 'origins':
            ParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
        
        ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
        ParamMap['DefaultPixelValue'] = ['0']
        #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
        #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
        #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        ParamMap['FixedInternalImagePixelType'] = ['float']
        ParamMap['MovingInternalImagePixelType'] = ['float']
        #ParamMap['HowToCombineTransforms'] = ['Compose'] # 14/05/21
        #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
        #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
        #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
        #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
        #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
        if MaxNumOfIters != 'default':
            ParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
        #                      'TransformBendingEnergyPenalty'] # default for bspline
        #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
        ParamMap['NumberOfHistogramBins'] = ['32']
        #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
        #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
        #    ParamMap['Metric0Weight'] = ['0.5']
        #    ParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
        #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
        #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
        #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
        #ParamMap['Transform'] = ['AffineTransform'] # default for affine
        #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
        ParamMap['UseDirectionCosines'] = ['true']
        ParamMap['WriteIterationInfo'] = ['true']
        
        """ Set the parameter map: """
        SelxImFilt.SetParameterMap(ParamMap)
        
    elif InitMethod != 'landmarks' and Transform == 'bspline':
        #print('case 2\n')
        
        """ This will require two parameter maps: an affine and a bspline. """
        AffineParamMap = sitk.GetDefaultParameterMap('affine')
        BsplineParamMap = sitk.GetDefaultParameterMap('bspline')
        
        """ Set the parameters for the affine parameter map: """
        AffineParamMap['AutomaticTransformInitialization'] = ['false']
        if InitMethod in ['centerofgravity', 'moments']:
            AffineParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
        elif InitMethod in ['geometricalcenter', 'geometry']:
            AffineParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
        elif InitMethod == 'origins':
            AffineParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
        AffineParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
        AffineParamMap['DefaultPixelValue'] = ['0']
        #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
        #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
        #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        AffineParamMap['FixedInternalImagePixelType'] = ['float']
        AffineParamMap['MovingInternalImagePixelType'] = ['float']
        AffineParamMap['HowToCombineTransforms'] = ['Compose']
        #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
        #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
        #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
        #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
        #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
        if MaxNumOfIters != 'default':
            AffineParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
        #                      'TransformBendingEnergyPenalty'] # default for bspline
        #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
        AffineParamMap['NumberOfHistogramBins'] = ['32']
        #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
        #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
        #    ParamMap['Metric0Weight'] = ['0.5']
        #    ParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
        #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
        #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
        #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
        #ParamMap['Transform'] = ['AffineTransform'] # default for affine
        #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
        AffineParamMap['UseDirectionCosines'] = ['true']
        AffineParamMap['WriteIterationInfo'] = ['true']
        
        """ Set the parameters for the bspline parameter map: """
        BsplineParamMap['AutomaticTransformInitialization'] = ['false']
        BsplineParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
        BsplineParamMap['DefaultPixelValue'] = ['0']
        #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
        #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
        #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        BsplineParamMap['FixedInternalImagePixelType'] = ['float']
        BsplineParamMap['MovingInternalImagePixelType'] = ['float']
        BsplineParamMap['HowToCombineTransforms'] = ['Compose']
        #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
        #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
        #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
        #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
        #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
        if MaxNumOfIters != 'default':
            BsplineParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
        BsplineParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
                                     'TransformBendingEnergyPenalty'] # default for bspline
        BsplineParamMap['Metric0Weight'] = ['1.0']
        BsplineParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
        BsplineParamMap['NumberOfHistogramBins'] = ['32']
        #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
        #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
        #    ParamMap['Metric0Weight'] = ['0.5']
        #    ParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
        #if InitMethod == 'landmarks':
        #    """ MultiResolutionRegistration, which is default for affine, only
        #    allows for 1 metric. """
        #    ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration']
        #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
        #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
        #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
        #ParamMap['Transform'] = ['AffineTransform'] # default for affine
        #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
        BsplineParamMap['UseDirectionCosines'] = ['true']
        BsplineParamMap['WriteIterationInfo'] = ['true']
        
        """ Set the parameter map: """
        SelxImFilt.SetParameterMap(AffineParamMap)
        SelxImFilt.AddParameterMap(BsplineParamMap)
        
    elif InitMethod == 'landmarks' and Transform != 'bspline':
        #print('case 3\n')
        
        """ This will require two parameter maps: A translation
        transformation for landmark alignment and a rigid transformation: """
        TranslationParamMap = sitk.GetDefaultParameterMap('translation')
        ParamMap = sitk.GetDefaultParameterMap(Transform)
        
        """ Set the parameters accordingly: """
        ParamMap['AutomaticTransformInitialization'] = ['false']
        ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
        ParamMap['DefaultPixelValue'] = ['0']
        #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
        #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
        #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        ParamMap['FixedInternalImagePixelType'] = ['float']
        ParamMap['MovingInternalImagePixelType'] = ['float']
        #ParamMap['HowToCombineTransforms'] = ['Compose'] # 14/05/21
        #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
        #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
        #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
        #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
        #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
        if MaxNumOfIters != 'default':
            ParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
        #                      'TransformBendingEnergyPenalty'] # default for bspline
        #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
        ParamMap['NumberOfHistogramBins'] = ['32']
        #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
        #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
        #    ParamMap['Metric0Weight'] = ['0.5']
        #    ParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
        #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
        #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
        #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
        #ParamMap['Transform'] = ['AffineTransform'] # default for affine
        #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
        ParamMap['UseDirectionCosines'] = ['true']
        ParamMap['WriteIterationInfo'] = ['true']
        
        """ Set the parameter map: """
        SelxImFilt.SetParameterMap(TranslationParamMap)
        SelxImFilt.AddParameterMap(ParamMap)
        
    elif InitMethod == 'landmarks' and Transform == 'bspline':
        #print('case 4\n')
        
        """ This will require three parameter maps: A translation
        transformation for landmark alignment, an affine and a bspline: """
        TranslationParamMap = sitk.GetDefaultParameterMap('translation')
        
        """ Set the parameters for the affine parameter map: """
        AffineParamMap['AutomaticTransformInitialization'] = ['false']
        #if InitMethod in ['centerofgravity', 'moments']:
        #    AffineParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
        #elif InitMethod in ['geometricalcenter', 'geometry']:
        #    AffineParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
        #elif InitMethod == 'origins':
        #    AffineParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
        AffineParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
        AffineParamMap['DefaultPixelValue'] = ['0']
        #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
        #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
        #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        AffineParamMap['FixedInternalImagePixelType'] = ['float']
        AffineParamMap['MovingInternalImagePixelType'] = ['float']
        AffineParamMap['HowToCombineTransforms'] = ['Compose']
        #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
        #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
        #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
        #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
        #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
        if MaxNumOfIters != 'default':
            AffineParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
        #                      'TransformBendingEnergyPenalty'] # default for bspline
        #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
        AffineParamMap['NumberOfHistogramBins'] = ['32']
        #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
        #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
        #    ParamMap['Metric0Weight'] = ['0.5']
        #    ParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
        #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
        #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
        #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
        #ParamMap['Transform'] = ['AffineTransform'] # default for affine
        #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
        AffineParamMap['UseDirectionCosines'] = ['true']
        AffineParamMap['WriteIterationInfo'] = ['true']
        
        """ Set the parameters for the bspline parameter map: """
        BsplineParamMap['AutomaticTransformInitialization'] = ['false']
        BsplineParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
        BsplineParamMap['DefaultPixelValue'] = ['0']
        #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
        #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
        #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        BsplineParamMap['FixedInternalImagePixelType'] = ['float']
        BsplineParamMap['MovingInternalImagePixelType'] = ['float']
        BsplineParamMap['HowToCombineTransforms'] = ['Compose']
        #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
        #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
        #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
        #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
        #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
        if MaxNumOfIters != 'default':
            BsplineParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
        #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
        BsplineParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
                                     'TransformBendingEnergyPenalty'] # default for bspline
        BsplineParamMap['Metric0Weight'] = ['1.0']
        BsplineParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
        BsplineParamMap['NumberOfHistogramBins'] = ['32']
        #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
        #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
        #if Transform == 'bspline': # 22/04/2021
        #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
        #    ParamMap['Metric0Weight'] = ['0.5']
        #    ParamMap['Metric1Weight'] = ['1.0']
        #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
        #if InitMethod == 'landmarks':
        #    """ MultiResolutionRegistration, which is default for affine, only
        #    allows for 1 metric. """
        #    ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration']
        #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
        #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
        #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
        #ParamMap['Transform'] = ['AffineTransform'] # default for affine
        #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
        BsplineParamMap['UseDirectionCosines'] = ['true']
        BsplineParamMap['WriteIterationInfo'] = ['true']
        
        """ Set the parameter map: """
        SelxImFilt.SetParameterMap(TranslationParamMap)
        SelxImFilt.AddParameterMap(AffineParamMap)
        SelxImFilt.AddParameterMap(BsplineParamMap)
    
    if InitMethod == 'landmarks':
        SelxImFilt.SetFixedPointSetFileName(FixFidsFpath)
        SelxImFilt.SetMovingPointSetFileName(MovFidsFpath)
    
    
    if LogToConsole:
        print(f'Performing image registration using {Transform} transform...\n')
        
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Register the 3D images: """
    SelxImFilt.Execute()
    
    """ Get the registered image: """
    RegIm = SelxImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    """ Get the parameter map and transform parameter map: """
    ParamMapDict = GetParamMapFromSelxImFilt(SelxImFilt)
    TxParamMapDict = GetTxParamMapFromSelxImFilt(SelxImFilt)
    
    #return RegIm, SelxImFilt
    return RegIm, SelxImFilt, ParamMapDict, TxParamMapDict






def SelxBsplineRegUsingFiducials(FixIm, MovIm, 
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       LogToConsole=False):
    """
    Register two 3D SimpleITK images using a BSpline transformation and
    fiducials for initial alignment in SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file name (if the file is in the current working directory) or
        full file path of the txt file containing fiducials for the fixed image
        in the required format for indices.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file name (if the file is in the current working directory) or
        full file path of the txt file containing fiducials for the moving image
        in the required format for indices.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    SelxImFilt : Elastix image filter
        Elastix image registration transformation filter used to register/
        transform MovIm to FixIm.
    
    ParamMapDict : dictionary
        A dictionary containing the parameter map from SelxImFilt.
    
    TxParamMapDict : dictionary
        A dictionary containing the transform parameter map from SelxImFilt.
    
        
    Notes:
    *****
    
    29/05/20: Trying this instead of trying to change it for Transformix
    ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70
    ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    22/02/21 Try parameters in Appendix A in "Multimodal image registration by 
    edge attraction and regularization using a B-spline grid" by Klein et al 
    2011; and from the SimpleElastix docs:
    https://simpleelastix.readthedocs.io/ParameterMaps.html
    
    18/03/21 Kasper said that AutomaticTransformInitialization should be true 
    for rigid, and false for affine or bspine: 
    
    https://github.com/SuperElastix/elastix/issues/45
         
    - FixedImagePyramid should be FixedImageSmoothingPyramid for best results 
    (unless memory is really and issue).
    - AutomaticTransformInitialization should be true for rigid and false for 
    affine and bspline.
    - 250 number of iterations is also on the low side. If speed is really an 
    issue, 250 is fine, but if you can afford it use at least 512.
    - Depending on the problem you might want to use at least 3 resolutions. 
    Some places you use 2 which is on the low side.
    - ImageSampler should be RandomSparseMask if users use masks. Otherwise it 
    is fine.
    - If you use a BSplineInterpolationOrder of 1 you might as well use 
    LinearInterpolator which produces the same results but is faster.
    
    14/05/21 Use of FinalLinearInterpolator reduced time by 12 s (2.6%).
    """
    
    import SimpleITK as sitk
    import time
    
    SelxImFilt = sitk.ElastixImageFilter()
    SelxImFilt.LogToConsoleOn()
    SelxImFilt.LogToFileOn()
    
    SelxImFilt.SetFixedImage(FixIm)
    SelxImFilt.SetMovingImage(MovIm)
    
    SelxImFilt.SetFixedPointSetFileName(FixFidsFpath)
    SelxImFilt.SetMovingPointSetFileName(MovFidsFpath)
    
    SelxImFilt.SetParameterMap(sitk.GetDefaultParameterMap('translation'))
    #SelxImFilt.SetParameterMap(sitk.GetDefaultParameterMap('rigid'))
    SelxImFilt.AddParameterMap(sitk.GetDefaultParameterMap('affine'))
    SelxImFilt.AddParameterMap(sitk.GetDefaultParameterMap('bspline'))
    
    #SelxImFilt.SetParameter("Transform", "BSplineTransform")
    
    SelxImFilt.SetParameter( "Metric", ("AdvancedMattesMutualInformation", 
                                        "TransformBendingEnergyPenalty",
                                        "CorrespondingPointsEuclideanDistanceMetric",))
    #SelxImFilt.SetParameter("Metric0Weight", "0.3")
    SelxImFilt.SetParameter("Metric0Weight", "1.0")
    SelxImFilt.SetParameter("Metric1Weight", "1.0")
    SelxImFilt.SetParameter("Metric2Weight", "1.0")
    #SelxImFilt.SetParameter("MaximumNumberOfIterations" , '500')
    #SelxImFilt.SetParameter("NumberOfSpatialSamples" , '3000')
    #SelxImFilt.SetParameter("Registration", "AdvancedImageToImageMetric") # error:
    # The following component could not be created: AdvancedImageToImageMetric
    SelxImFilt.SetParameter("Registration", "MultiMetricMultiResolutionRegistration")
    #SelxImFilt.SetParameter("ResampleInterpolator", "FinalBSplineInterpolator") # recommended in elastix manual (5.3.4 Interpolator)
    SelxImFilt.SetParameter("ResampleInterpolator", "FinalLinearInterpolator") # reduces computational time
    
    times = []
    times.append(time.time())
    
    SelxImFilt.Execute()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'*Took {Dtime} s to perform image registration.\n')
    
    RegIm = SelxImFilt.GetResultImage()
    
    ParamMapDict = GetParamMapFromSelxImFilt(SelxImFilt)
    TxParamMapDict = GetTxParamMapFromSelxImFilt(SelxImFilt)
    
    return RegIm, SelxImFilt, ParamMapDict, TxParamMapDict
    
    



def ParamMapFromSelxImFilt_OLD(SelxImFilt, Dict=None):
    """
    Get a dictionary containing the parameter map from an Elastix image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    Dict : dictionary (optional; None by default)
        If a dictionary (e.g. containing parameters from 
        .GetTransformParameterMap()) is provided the parameters (from 
        .GetParameterMap()) will be added to it.  Any keys in Dict (e.g. 
        TransformParameterMap) that are also present in ParameterMap will be 
        overwritten. 
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing Elastix parameters.
    """
    
    if not Dict:
        Dict = {}
    
    for key, value in SelxImFilt.GetParameterMap()[0].items():
        Dict.update({key : value})
    
    return Dict





def GetParamMapFromSelxImFilt(SelxImFilt):
    """
    Get a dictionary containing the parameter map(s) from an Elastix image 
    filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing the parameter values used in the Elastix image 
        filter.  
    
    
    Note:
    ****
    
    The keys of Dict at the first level are integer strings that indicate the 
    parameter map number, e.g. Dict['0'], Dict['1'], Dict['2'], etc.
    
    If there was only one parameter map, there will only be one key,
    e.g. Dict['0'].
    
    Each first level key corresponds to a dictionary containing the parameters 
    for that parameter map, 
    
    e.g. Dict['0'] = {Param0 : value0,
                      Param1 : value1,
                      Param2 : value2,
                      ...}
    """
    
    Dict = {}
    
    for i in range(len(SelxImFilt.GetParameterMap())):
        ParamMap = {}
        
        for key, value in SelxImFilt.GetParameterMap()[i].items():
            ParamMap.update({key : value})
        
        Dict.update({str(i) : ParamMap})
    
    return Dict







def TxParamMapFromSelxImFilt_OLD(SelxImFilt, Dict=None):
    """
    Get a dictionary containing the transform parameter map from an Elastix
    image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    Dict : dictionary (optional; None by default)
        If a dictionary (e.g. containing parameters from .GetParameterMap()) is 
        provided the transform parameters (from .GetTransformParameterMap())
        will be added to it.  Any keys in Dict (e.g. ParameterMap) that are 
        also present in TransformParameterMap will be overwritten. 
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing Elastix parameters.
    
    
    Notes:
    *****
    
    1. The transform parameter map can be obtained as follows:
    
    SelxImFilt.GetTransformParameterMap()[0]
    
    and printed as follows:
        
    sitk.PrintParameterMap(SelxImFilt.GetTransformParameterMap()[0])
    
    however printing to the console does not work in Jupyter.  Although it is
    possible to print the value of a specific key to the console as follows:
        
    SelxImFilt.GetTransformParameterMap()[0]['Direction']
    
    or of the entire dictionary:
        
    for key, value in SelxImFilt.GetTransformParameterMap()[0].items():
        print(f'{key} = {value}')
        
    Hence this function is simply a way of getting around the Jupyter notebook
    limitation.
    
    2. Oddly the transform parameters read from file (TransformParameters.0.txt)
    have more significant digits than that read directly from the Elastix image 
    filter.
    
    e.g. 
    
    list(SelxImFilt.GetTransformParameterMap()[0]['TransformParameters'])
    
    returns:
        
    ['1.0046',
     '-0.00831448',
     '-0.0412239',
     '0.0140337',
     '0.997238',
     '0.0593744',
     '0.0420835',
     '-0.0564159',
     '1.00477',
     '-1.70085',
     '-11.4156',
     '-54.8056']
    
    whereas

    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')

    TxParams = [str(item) for item in TxParams]
    
    TxParams
    
    returns:
    
    ['1.004604',
     '-0.008314',
     '-0.041224',
     '0.014034',
     '0.997238',
     '0.059374',
     '0.042084',
     '-0.056416',
     '1.004766',
     '-1.700848',
     '-11.415621',
     '-54.805552']

    """
    
    if not Dict:
        Dict = {}
    
    for key, value in SelxImFilt.GetTransformParameterMap()[0].items():
        Dict.update({key : value})
    
    return Dict





def GetTxParamMapFromSelxImFilt(SelxImFilt):
    """
    Get a dictionary containing the transform parameter map(s) from an Elastix 
    image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing the transform parameter values used in the 
        Elastix image filter.  
    
    
    Note:
    ****
    
    The keys of Dict at the first level are integer strings that indicate the 
    parameter map number, e.g. Dict['0'], Dict['1'], Dict['2'], etc.
    
    If there was only one transform parameter map, there will only be one key,
    e.g. Dict['0'].
    
    Each first level key corresponds to a dictionary containing the transform
    parameters for that transform parameter map, 
    
    e.g. Dict['0'] = {Param0 : value0,
                      Param1 : value1,
                      Param2 : value2,
                      ...}
    """
    
    Dict = {}
    
    for i in range(len(SelxImFilt.GetTransformParameterMap())):
        ParamMap = {}
        
        for key, value in SelxImFilt.GetTransformParameterMap()[i].items():
            ParamMap.update({key : value})
        
        Dict.update({str(i) : ParamMap})
    
    return Dict






def GetParamOfFinalParamMapFromSelxImFilt(SelxImFilt, Param):
    """
    Get the parameter value(s) for a parameter of interest from the last  
    parameter map of an Elastix image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    Param : string
        The key corresponding to the parameter of interest,
        e.g. 'TransformParameters'; 'GridSize'; 'GridSpacing'.
    
    Outputs:
    *******
    
    ParamVal : string or list of strings
        The parameter value(s) of parameter Param for the final transform 
        parameter map in SelxImFilt.
    """
    
    # Get all transform parameter maps:
    Dict = GetTxParamMapFromSelxImFilt(SelxImFilt)
    
    keys = list(Dict.keys())
    last = keys[-1]
    
    ParamVal = list(Dict[last][Param])
    
    if len(ParamVal) == 1:
        ParamVal = ParamVal[0]
    
    return ParamVal





def GetParamFromSelxImFilt(SelxImFilt, ParamMapInd, Param):
    """
    Get the parameter value for a parameter of interest from a parameter 
    map of interest of an Elastix image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    ParamMapInd : integer
        The index of the parameter map containing the parameter of interest.
    
    Param : string
        The key corresponding to the parameter of interest,
        e.g. 'TransformParameters'; 'GridSize'; 'GridSpacing'.
    
    Outputs:
    *******
    
    ParamVal : string or list of strings or None
        The parameter value(s) of parameter Param for the ParamMapInd^th 
        transform parameter map in SelxImFilt.  The value None will be returned
        if the required parameter does not exist in the chosen parameter map.
    
    
    Note:
    ****
    
    The keys of Dict at the first level are integer strings that indicate the 
    parameter map number, e.g. Dict['0'], Dict['1'], Dict['2'], etc.
    
    If there was only one transform parameter map, there will only be one key,
    e.g. Dict['0'].
    
    Each first level key corresponds to a dictionary containing the transform
    parameters for that transform parameter map, 
    
    e.g. Dict['0'] = {Param0 : value0,
                      Param1 : value1,
                      Param2 : value2,
                      ...}
    
    Hence, to get the parameter value for parameter key 'Param5' in the first
    parameter map:
        GetParamFromSelxImFilt(SelxImFilt, 0, 'Param5')
        
    To get the parameter value for parameter key 'Param5' in the last parameter
    map:
        GetParamFromSelxImFilt(SelxImFilt, -1, 'Param5')
    
    To get the parameter value for parameter key 'Param5' in the second last 
    parameter map:
        GetParamFromSelxImFilt(SelxImFilt, -2, 'Param5')
    """
    
    # Get all transform parameter maps:
    Dict = GetTxParamMapFromSelxImFilt(SelxImFilt)
    
    ParamMapKeys = list(Dict.keys())
    ParamMapKey = ParamMapKeys[ParamMapInd]
    
    try:
        ParamVal = list(Dict[ParamMapKey][Param])
        
        if len(ParamVal) == 1:
            ParamVal = ParamVal[0]
    except KeyError:
        ParamVal = None
    
    return ParamVal






def TxParamMapFromFile(Dict=None):
    """
    Get the Elastix transform parameter map from TransformParameters.0.txt in  
    the current working directory.
    
    Inputs:
    ******
    
    Dict : dictionary (optional; None by default)
        If a dictionary (e.g. containing parameters from .GetParameterMap()) is 
        provided the transform parameters (from TransformParameters.0.txt)
        will be added to it.  Any keys in Dict (e.g. ParameterMap) that are 
        also present in TransformParameters.0.txt will be overwritten.
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing Elastix transform parameters (and other 
        parameters if an input dictionary was provided).
        If the file TransformParameters.0.txt is not found in the current 
        working directory an empty (or the input) dictionary will be returned.
    
    
    Notes:
    *****
    
    The transform parameters read from file (TransformParameters.0.txt)
    have more significant digits than that read directly from the Elastix image 
    filter.
    
    e.g. 
    
    list(SelxImFilt.GetTransformParameterMap()[0]['TransformParameters'])
    
    returns:
        
    ['1.0046',
     '-0.00831448',
     '-0.0412239',
     '0.0140337',
     '0.997238',
     '0.0593744',
     '0.0420835',
     '-0.0564159',
     '1.00477',
     '-1.70085',
     '-11.4156',
     '-54.8056']
    
    whereas

    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')

    TxParams = [str(item) for item in TxParams]
    
    TxParams
    
    returns:
    
    ['1.004604',
     '-0.008314',
     '-0.041224',
     '0.014034',
     '0.997238',
     '0.059374',
     '0.042084',
     '-0.056416',
     '1.004766',
     '-1.700848',
     '-11.415621',
     '-54.805552']

    """
    
    import os
    
    if Dict == None:
        Dict = {}
    
    Fname = 'TransformParameters.0.txt'
    
    if os.path.isfile(Fname):
        # Read in the file:
        with open(Fname, 'r') as FileText:
            lines = FileText.readlines()
            
            
            for line in lines:
                if line[0] == '(':
                    # Get the index of the first white space and ')':
                    a = line.index(' ')
                    b = line.index(')')
                    
                    key = line[1:a]
                    value = line[a+1:b].replace('"', '')
                    
                    # Get values as a list of space-delimited strings:
                    values = value.split()
                    
                    Dict.update({key : tuple(values)})

    
    return Dict






def ParseTransformParameterFile(Fpath, ParamName):
    """
    Parse an Elastix transform parameter map from TransformParameters.0.txt for 
    the value(s) of a specific parameter.
            
    Inputs:
    ******
    
    Fpath : string
        The full filepath of TransformParameters text file.
        
    ParamName : string
        Character string of the parameter name to search and parse.
        
        
    Outputs:
    *******
    
    ParamValues : strings, integers or floats (depends on the parameter of 
    interest)
        The value(s) of the parameter of interest.
    
    
    Notes:
    -----
    
    1.  The text file contains data of the following format:
            
        (FixedImageDimension 2)     <-- single value
        (Size 256 256)              <-- two values
        
    2. See further notes below.
    """

    # Determine whether the parameters need to be converted to integers or
    # float.  
    
    # Create a list of parameter names that belong to int and float data types:
    """
    Notes: 
        1. Many of the parameters are to remain as strings.
        2. The lists below only cover the default parameters, so any additional
           parameters will need to be added.
        3. If the parameter is not found to be in ints or floats it will be 
           assumed that the parameter is a string.
    """
    
    IntParams = ['NumberOfParameters', 'FixedImageDimension', 
                 'MovingImageDimension', 'Size', 'Index', 
                 'FinalBSplineInterpolationOrder']
    
    FloatParams = ['TransformParameters', 'Spacing', 'Origin', 'Direction', 
                   'CenterOfRotationPoint', 'DefaultPixelValue']
    
    if ParamName in IntParams:
        dtype = 'int'
    elif ParamName in FloatParams:
        dtype = 'float'
    else:
        dtype = 'string'
        
    
    # Read in the file:
    FileTxt = open(Fpath, 'r').read()
    
    #print('\nFileTxt:\n', FileTxt)
    
    # Get index of string variable parameter:
    a = FileTxt.index(ParamName)
    
    #print(f'\na = {a}')
    
    """ Note:
    The text file has values encoded in the following format:
    (ParameterName ParameterValue1 ParameterValue2 ...)
    
    so the values begin at a+len(ParameterName)+1 and end at the position 
    before ')'.
    """
    StartInd = a + len(ParamName) + 1 # to account for the space after 
    # ParameterName
    
    #print(f'\nStartInd = {StartInd}')
    
    # Find the position of the first ')' after StartInd:
    b = StartInd + FileTxt[StartInd:].index(')')
    
    #print(f'\nb = {b}')
    
    EndInd = b
    
    #print(f'\nEndInd = {EndInd}')
    
    # Crop the text between StartInd and EndInd:
    CroppedTxt = FileTxt[StartInd:EndInd]
    
    #print(f'\nCroppedTxt = {CroppedTxt}')
    
    """ Note:
    CroppedTxt may be single value or may have multiple values. 
    e.g. (FixedImageDimension 2) <-- single value
         (Size 256 256) <-- two values
         
    If CroppedTxt doesn't contain the space character simply return CroppedTxt.
    If CroppedTxt has a space character, slice and append the elements that 
    preceed it to a new array called ParameterValues, then search for another 
    space, repeat until there are no more spaces.
    """
    # Check if there are any space characters in CroppedTxt:
    AreThereSpaces = ' ' in CroppedTxt
    
    #print('\nAreThereSpaces =', AreThereSpaces)
    
    if AreThereSpaces:
        ParamValues = []
        
        # Continue as long as AreThereSpaces is True:
        while AreThereSpaces:
            # Find the first space character:
            ind = CroppedTxt.index(' ')

            # Slice the characters from 0 to ind:
            ParamValues.append(CroppedTxt[0:ind])
            
            # Remove from CroppedTxt the characters that were appended to 
            # ParamValues but not including the space character:
            CroppedTxt = CroppedTxt[ind+1:]
            
            #print(f'\n   CroppedTxt = {CroppedTxt}')
            
            # Check if there are any space characters in CroppedTxt:
            AreThereSpaces = ' ' in CroppedTxt
            
            #print('\n   AreThereSpaces =', AreThereSpaces)
            
        # There are no remaining spaces so simply append what remains of
        # CroppedTxt to ParamValues:
        if CroppedTxt:
            ParamValues.append(CroppedTxt)
            
        #return ParamValues
        if dtype=='int':
            return [int(item) for item in ParamValues]
        elif dtype=='float':
            return [float(item) for item in ParamValues]
        else:
            return ParamValues
        
    else:
        #return CroppedTxt
        if dtype=='int':
            return [int(item) for item in CroppedTxt]
        elif dtype=='float':
            return [float(item) for item in CroppedTxt]
        else:
            return CroppedTxt
        







def CompareSelxParameterMaps(SelxImFilt0, SelxImFilt1, DiffsOnly=True):
    ParamMap0 = SelxImFilt0.GetParameterMap()[0]
    ParamMap1 = SelxImFilt1.GetParameterMap()[0]
    
    Keys0 = ParamMap0.keys()
    Keys1 = ParamMap1.keys()
    
    print('SelxImFilt0 parameter v SelxImFilt1 parameter')
    print('*********************************************')
    
    for key0 in Keys0:
        if key0 in Keys1:
            if DiffsOnly:
                if ParamMap0[key0] != ParamMap1[key0]:
                    print(f"\n{key0}:")
                    print(f"{ParamMap0[key0]} v {ParamMap1[key0]}")
                else:
                    print(f"\n{key0}:")
                    print(f"{ParamMap0[key0]} v (key not present)")
            else:
                print(f"\n{key0}:")
                print(f"{ParamMap0[key0]} v {ParamMap1[key0]}")
    
    for key1 in Keys1:
        if not key1 in Keys0:
            if DiffsOnly:
                if ParamMap0[key1] != ParamMap1[key1]:
                    print(f"\n{key1}:")
                    print(f"{ParamMap0[key1]} v {ParamMap1[key1]}")
                else:
                    print(f"\n{key1}:")
                    print(f"(key not present) v {ParamMap1[key1]}")
            else:
                print(f"\n{key1}:")
                print(f"{ParamMap0[key1]} v {ParamMap1[key1]}")

    return





def CreateSelxImFiltFromDro(Dro, FixIm, MovIm, TxMatrix, Transform, LogToConsole):
    """
    Create a SimpleElastix image filter from a DRO.
    
    list of elements that make up the
    transform matrix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image
        The 3D image whose grid is to be transformed to.
    
    TxMatrix : list of float strings
        List of float strings representing the transformation that maps the
        moving (Source) image to the fixed (Target) image (e.g. as parsed from
        the FrameOfReferenceTransformationMatrix of a DRO).
    
    SitkTransform : SimpleITK Transform
        The transform (e.g. image registration transform) that maps Im to TrgIm.
    
        
    Outputs:
    *******
    
    SitkTransform : SimpleITK transform
        The 3D transformed image.
    
    
    Note:
    ****
    
    TxMatrix has additional bottom row [0 0 0 1] (see 
    http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
    
    Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. Hence
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
    """
    
    import SimpleITK as sitk
    
    mod = Dro.Modality
    
    if mod != 'REG':
        msg = "The DRO does not have Modality 'REG'. Modality is {mod}."
        raise Exception(msg)
    
    
    """ The indices of the elements in TxMatrix that relate to TxParams: """
    inds = [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]
    
    """ Initiate ElastixImageFilter: """
    SelxImFilt = sitk.ElastixImageFilter()
    if LogToConsole:
        SelxImFilt.LogToConsoleOn() # no output in Jupyter
    else:
        SelxImFilt.LogToConsoleOff() 
    SelxImFilt.LogToFileOn() # output's elastix.log ONLY ONCE per kernel in Jupyter
    
    SelxImFilt.SetFixedImage(FixIm)
    SelxImFilt.SetMovingImage(MovIm)
    
    """ Initialise the parameter map. """
    ParamMap = sitk.VectorOfParameterMap()
    
    if 'Deformable' in Dro.SOPClassUID:
        """ Get the pre-deformation FOR transformation matrix and matrix type. 
        
        Note:
            The DRO contains two sequences in DeformableRegistrationSequence
            - the first relates to the fixed image series, the second to the 
            moving image series. It is the latter sequence that contains 
            PreDeformationMatrixRegistrationSequence, 
            PostDeformationMatrixRegistrationSequence, and
            DeformableRegistrationGridSequence.
        """
        PreDefTxMatrixType = Dro.DeformableRegistrationSequence[1]\
                                .PreDeformationMatrixRegistrationSequence[0]\
                                .FrameOfReferenceTransformationMatrixType\
                                .lower()
        
        PreDefTxMatrix = Dro.DeformableRegistrationSequence[1]\
                            .PreDeformationMatrixRegistrationSequence[0]\
                            .FrameOfReferenceTransformationMatrix
        
        """ If the matrix type is 'rigid' determine if the matrix is not the
        identity matrix. """
        if PreDefTxMatrixType == 'rigid':
            from GeneralTools import AreListsEqualToWithinEpsilon
            
            IdentityMatrix = [1.0, 0.0, 0.0, 0.0, 
                              0.0, 1.0, 0.0, 0.0, 
                              0.0, 0.0, 1.0, 0.0, 
                              0.0, 0.0, 0.0, 1.0]
            
            matrix = [float(item) for item in PreDefTxMatrix]

            if AreListsEqualToWithinEpsilon(matrix, IdentityMatrix):
                """ There's no need to apply the pre-deformation matrix. """
                PreDefTxMatrix = None
            else:
                PreDefTxParams = [PreDefTxMatrix[i] for i in inds]
                
        elif PreDefTxMatrixType == 'rigid_scale':
            """ Need to deal with this case. """
        
        
        """ Get the post-deformation FOR transformation matrix and matrix type. 
        """
        PostDefTxMatrixType = Dro.DeformableRegistrationSequence[1]\
                                 .PostDeformationMatrixRegistrationSequence[0]\
                                 .FrameOfReferenceTransformationMatrixType\
                                 .lower()
        
        PostDefTxMatrix = Dro.DeformableRegistrationSequence[1]\
                             .PostDeformationMatrixRegistrationSequence[0]\
                             .FrameOfReferenceTransformationMatrix
        
        """ If the matrix type is 'rigid' determine if the matrix is not the
        identity matrix. """
        if PostDefTxMatrixType == 'rigid':
            from GeneralTools import AreListsEqualToWithinEpsilon
            
            IdentityMatrix = [1.0, 0.0, 0.0, 0.0, 
                              0.0, 1.0, 0.0, 0.0, 
                              0.0, 0.0, 1.0, 0.0, 
                              0.0, 0.0, 0.0, 1.0]
            
            matrix = [float(item) for item in PreDefTxMatrix]

            if AreListsEqualToWithinEpsilon(matrix, IdentityMatrix):
                """ There's no need to apply the pre-deformation matrix. """
                PostDefTxMatrix = None
            else:
                PostDefTxParams = [PostDefTxMatrix[i] for i in inds]
                
        elif PostDefTxMatrixType == 'rigid_scale':
            """ Need to deal with this case. """
        
        
        """ Get the grid origin, direction, dimensions and resolution, and 
        transform parameters from DeformableRegistrationGridSequence. """
        GridOrig = Dro.DeformableRegistrationSequence[1]\
                      .DeformableRegistrationGridSequence[0]\
                      .ImagePositionPatient
        
        GridDir = Dro.DeformableRegistrationSequence[1]\
                     .DeformableRegistrationGridSequence[0]\
                     .ImageOrientationPatient
        
        GridDims = Dro.DeformableRegistrationSequence[1]\
                      .DeformableRegistrationGridSequence[0]\
                      .GridDimensions
        
        GridRes = Dro.DeformableRegistrationSequence[1]\
                     .DeformableRegistrationGridSequence[0]\
                     .GridResolution
        
        DefTxParams = Dro.DeformableRegistrationSequence[1]\
                         .DeformableRegistrationGridSequence[0]\
                         .VectorGridData
        
        
        """ Add to the parameter map. """
        if PreDefTxMatrix:
            """ Note: 
                Need to deal with case PreDefTxMatrixType = 'rigid_scale'. 
            """
            ParamMap.append(sitk.GetDefaultParameterMap(PreDefTxMatrixType))
        
        ParamMap.append(sitk.GetDefaultParameterMap('bspline'))
        
        if PostDefTxMatrix:
            """ Note: 
                Need to deal with case PostDefTxMatrixType = 'rigid_scale'. 
            """
            ParamMap.append(sitk.GetDefaultParameterMap(PostDefTxMatrixType))
        
            
    

    else:
        """ Get the FOR transformation matrix and matrix type. 
        
        Note:
            The DRO contains two sequences in RegistrationSequence - the first
            relates to the fixed image series, the second to the moving image 
            series. It is the latter sequence that contains the matrix and 
            matrix type of interest.
        """
        TxMatrixType = Dro.RegistrationSequence[1]\
                          .MatrixRegistrationSequence[0]\
                          .MatrixSequence[0]\
                          .FrameOfReferenceTransformationMatrixType\
                          .lower()
        
        TxMatrix = Dro.RegistrationSequence[1]\
                      .MatrixRegistrationSequence[0]\
                      .MatrixSequence[0]\
                      .FrameOfReferenceTransformationMatrix
        
        """ Add to the parameter map.
        Note: 
            Need to deal with case PreDefTxMatrixType = 'rigid_scale'. 
        """
        ParamMap.append(sitk.GetDefaultParameterMap(TxMatrixType))
    
    
    ParamMap['HowToCombineTransforms'] = ['Compose']
    
    
    #if 
    
    
    return SelxImFilt






def TransformImageSelx(MovIm, SelxImFilt, Interp='NearestNeighbor', 
                       LogToConsole=False):
    """
    Transform a 3D SimpleITK image using SimpleElastix's Transformix.
    
    Inputs:
    ******
    
    MovIm : SimpleITK image 
        The 3D image to be transformed.
        
    SelxImFilt : Elastix image filter
        SimpleElastix image transformation filter used to perform image 
        registration.
        
    Interp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
        
    Outputs:
    *******
    
    TxIm : SimpleITK image
        The 3D transformed image.
        
        
    Notes:
    *****
    
    Simple use of Transformix as below:
        
    TxIm = sitk.Transformix(MovIm, ElxImFilt.GetTransformParameterMap())
        
    
    
    https://github.com/SuperElastix/SimpleElastix/issues/409
    """
    
    import SimpleITK as sitk
            
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of TransformImageElx():')
        
    # Initiate TransformixImageFilter:
    TxImFilt = sitk.TransformixImageFilter()
    #TxImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TxImFilt.LogToConsoleOff() # <-- no output in Jupyter
    TxImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    TxImFilt.SetMovingImage(MovIm)
    
    # Get the parameter map used for the registration:
    TxParamMap = SelxImFilt.GetParameterMap()
    
    # Set the parameter map:
    #TxImFilt.SetTransformParameterMap(SelxImFilt.GetTransformParameterMap())
    
    """ Get the transform parameter map that was used for registration: """
    TxMap = SelxImFilt.GetTransformParameterMap()
    
    if LogToConsole:
        print('TxParamMap[0]["Interpolator"] =',
              f'{TxParamMap[0]["Interpolator"]}\n')
        print(f'TxMap[0]["ResampleInterpolator"] =',
              f'{TxMap[0]["ResampleInterpolator"]}\n')
        print(f'TxMap[0]["FinalBSplineInterpolationOrder"] =',
              f'{TxMap[0]["FinalBSplineInterpolationOrder"]}\n')
    
    
    if Interp == 'NearestNeighbor':
        interp = TxMap[0]["ResampleInterpolator"]
        
        TxMap[0]["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]
        
        if LogToConsole:
            print(f"The resample interpolator has been changed from {interp}",
                  "to 'FinalNearestNeighborInterpolator'\n")
    
    
    """ Get the image data type: (new 16/03/21)"""
    dtype = MovIm.GetPixelIDTypeAsString()
    
    if 'integer' in dtype:
        a = TxMap[0]["FinalBSplineInterpolationOrder"]
        
        """ The image is binary so set the FinalBSplineInterpolationOrder to 0
        to avoid "overshoot"-related "garbage" as described in the elastix
        manual section (4.3 The transform parameter file): """
        TxMap[0]["FinalBSplineInterpolationOrder"] = ["0"]
        
        if LogToConsole:
            print("The final BSpline interpolation order has been changed",
                  f"from {a} to '0'\n")
        
    TxImFilt.SetTransformParameterMap(TxMap)
    
    # Set up for computing deformation field:
    #TxImFilt.ComputeDeformationFieldOn()
    
    """ Transform MovIm: """
    TxImFilt.Execute()
    
    TxIm = TxImFilt.GetResultImage()
    
    if LogToConsole:
        print('-'*120)
    
    return TxIm







