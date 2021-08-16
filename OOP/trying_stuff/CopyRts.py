# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 11:15:07 2021

@author: ctorti
"""

def CopyRts(SrcRtsFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgRtsFpath=None, TrgRoiName=None, TrgSliceNum=None, 
            ResInterp='BlurThenLinear', PreResVariance=(1,1,1),
            ApplyPostResBlur=False, PostResVariance=(1,1,1), 
            ForceReg=False, UseDroForTx=True, SelxOrSitk='Sitk',
            Transform='affine', MaxIters='512', InitMethod='centerofgravity', 
            SrcFidsFpath='moving_fiducials.txt', 
            TrgFidsFpath='target_fiducials.txt', TxInterp='NearestNeighbor', 
            ApplyPostTxBin=True, ApplyPostTxBlur=True, PostTxVariance=(1,1,1),
            #TxMatrix=None, GridDims=None, GridRes=None, VectGridData=None,
            Dro=None, TxtToAddToTrgRoiName='', #OverwriteDro=False, 
            LogToConsole=False,
            DictOfInputs={}, ListOfInputs=[], ListOfTimings=[]):
    """
    Note 01/02/2021:
        This function will, depending on the inputs, copy either:
            1. A specific contour
            2. All contours in a specific ROI
            3. All contours in all ROIs
        
     
    Inputs:
    ******
    
    SrcRtsFpath : string
        Filepath of the Source RTS file.
    
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
        
    SrcRoiName : string
        All or part of the ROIName of the ROI containing the contour(s) to be
        copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
        
    TrgRtsFpath : string or None (optional; None by default)
        Filepath to the Target RTS file that the contour(s) is/are to be copied 
        to. TrgRtsFpath != None only if an existing Target RTS exists and is to
        be added to. An existing RTS will be added to only for Direct copy 
        operations.
    
    TrgRoiName : string or None (optional; None by default) 
        All or part of the ROIName of the destination ROI.
        
    TrgSliceNum : integer (optional; None by default)
        The slice index within the Target DICOM stack where the contour will be
        copied to (applies only for the case of direct copies of single 
        contours). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the contour(s) will be copied to will
        not depend on user input.
    
    ResInterp : string
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    PreResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to perform binary thresholding following
    #    resampling using a non-label (e.g. linear) interpolator. Example use is 
    #    on a labelmap image if a Gaussian blurring + linearly resampling 
    #    approach is taken to combat aliasing effects.
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
        
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
        
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, even 
        if the Source and Target images have the same FrameOfReferenceUID.
    
    UseDroForTx : boolean (optional; True by default)
        If True the transformation parameters from any suitable DRO (for the
        required Transform) will be used instead of performing image
        registration.
    
    SelxOrSitk : string (optional; 'Sitk' by default)
        Denotes which package to use for image registration and transformation.
        Acceptable values include:
            - 'Selx' for SimpleElastix
            - 'Sitk' for SimpleITK
     
    Transform : string (optional; 'affine' by default)
        The transformation to used for image registration (if applicable).  
        Acceptable values are:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Transform. If != 'default' it must be a string 
        representation of an integer.
    
    InitMethod : string (optional, 'centerofgravity' by default)
        The initialisation method used for initial alignment prior to image
        registration (if applicable).
    
    SrcFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Source/Moving image.
    
    TrgFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Target/Fixed image.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
        
    #ThreshPostTx : float (optional; 0.05 by default)
    #    The threshold level used to perform binary thresholding after 
    #    registration transformation.  
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
        
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
    
    TxtToAddToTrgRoiName : string (optional, '' by default)
        String of text to add to the ROI Name of the new Target RTS.
    
    #TxMatrix : list of float strings or None (optional; None by default)
    #    List of float strings representing the non-deformable transformation 
    #    that transforms the moving (Source) image to the fixed (Target) image;
    #    or None.  If not None, image registration will be skipped to save 
    #    computational time.
    #
    #GridDims : list of integers or None (optional; None by default)
    #    List of integers representing the dimensions of the BSpline grid used 
    #    to deform the moving (Source) image to the fixed (Target) image; 
    #    or None.  If not None, image registration will be skipped to save 
    #    computational time.
    # 
    #GridRes : list of floats or None (optional; None by default)
    #    List of floats representing the resolution of the BSpline grid used to 
    #    deform the moving (Source) image to the fixed (Target) image; or None.
    #    If not None, image registration will be skipped to save computational
    #    time.
    #
    #VectGridData : list of float strings or None (optional; None by default)
    #    List of floats representing the vector deformations that deform the
    #    moving (Source) image to the fixed (Target) image; or None.
    #    If not None, image registration will be skipped to save computational
    #    time.
    
    Dro : Pydicom Object (optional; None by default)
        A Spatial or Deformable Spatial Registration Object that contains the
        parameters required to register the moving (Source) image to the fixed
        (Target) image.
    
    OverwriteDro : boolean (optional; False by default)
        If True, Dro (if not None) will be overwritten.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    DictOfInputs : dictionary (optional; empty by default)
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings (optional; empty by default)
        A list containing the inputs that were called to CopyRoi().
        
    LogOfTimings : list of strings (optional; empty by default)
        A list of the time to execute certain tasks during the calling of 
        CopyRts().
          
                  
    Outputs:
    *******
        
    TrgRts : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target RTS object.
    
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgRts, None otherwise.
        
    DictOfInputs : dictionary
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
        
    ListOfTimings : list of strings
        A list (for each message) of timings for individual operations.
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours will appear displaced w.r.t. 
    the anatomical features highlighted in the Source image.
    """
    
    import importlib
    #import RtsTools
    #importlib.reload(RtsTools)
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    import DroTools
    importlib.reload(DroTools)
    import RoiCopyTools
    importlib.reload(RoiCopyTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    ##from RtsTools import GetPtsByCntByRoi
    #from RtsTools import GetRtsDataOfInterest, ProportionOfRoisInExtent
    #from RtsTools import CreateRts
    from GeneralTools import UniqueItems, PrintIndsByRoi#, PrintTitle
    from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    from GeneralTools import ShiftPtsByCntByRoi
    #from GeneralTools import ShiftFrame
    #from DicomTools import GetRoiNum
    #from ConversionTools import Points2ContourData
    from ConversionTools import PtsByCntByRoi2CntDataByCntByRoi
    from ImageTools import ImportImage
    from DroTools import CreateSpaDro, CreateDefDro
    from RoiCopyTools import CheckValidityOfInputs
    
    """ Start timing. """
    times = []
    times.append(time.time())
    #""" Log timing messages. """
    #TimingMsgs = []
    
    
    PlotResults = DictOfInputs['PlotResults']
    ExportPlot = DictOfInputs['ExportPlot']
    RtsPlotExportDir = DictOfInputs['RtsPlotExportDir']

    
    
    SrcRts = dcmread(SrcRtsFpath)
    if TrgRtsFpath:
        TrgRts = dcmread(TrgRtsFpath)
    else:
        TrgRts = None
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print('Running of CopyRts():')
        print('\n\n', '-'*120)
        
        SrcRoiLabels = GetRoiLabels(SrcRts)
        
        print(f'\nSrcRoiLabels = {SrcRoiLabels}')
        
        if TrgRtsFpath:
            TrgRoiLabels = GetRoiLabels(TrgRts)
        
            print(f'\nTrgRoiLabels = {TrgRoiLabels}')
    
    
    """ Check the validity of the combination of inputs. """
    CheckValidityOfInputs(SrcRts, SrcRoiName, SrcSliceNum, TrgRts, TrgRoiName,
                          TrgSliceNum, LogToConsole)
    
    #""" Determine which Use Case to apply. """
    #UseCaseThatApplies,\
    #UseCaseToApply = WhichUseCase(SrcSliceNum, SrcRoiName, TrgSliceNum, 
    #                              SrcDcmDir, TrgDcmDir, ForceRegistration,
    #                              LogToConsole=True)
    
    #""" Manipulate UseCase if ForceRegistration = True and UseCase was 3 or 4. """
    #if ForceRegistration and not '5' in UseCase:
    #    msg = f'\nUseCase was {UseCase} but will be treated as UseCase '
    #    
    #    UseCase.replace('3', '5').replace('4', '5')
    #    
    #    msg = f'{UseCase}.'
    #    
    #    print(msg)
    
    """ Get the data of interest from the Source RTS. """
    SrcPtsByCntByRoi,\
    SrcC2SindsByRoi = GetRtsDataOfInterest(SrcRtsFpath, SrcSliceNum, 
                                           SrcRoiName, SrcDcmDir, 
                                           LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source RTS file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Load the Source RTS: """
    SrcRts = dcmread(SrcRtsFpath)
    
    """ Raise exception if there is no data of interest from the Source RTS."""
    if not UniqueItems(SrcC2SindsByRoi):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        #from RtsTools import GetNumOfContoursByRoi, GetCStoSliceIndsByRoi
        
        #SrcRts = dcmread(SrcRtsFpath)
        
        Names = GetRoiLabels(SrcRts)
        N = len(Names)
        
        NumByRoi = GetNumOfContoursByRoi(SrcRts)
        
        SOPuids = GetDicomSOPuids(SrcDcmDir)
        
        C2SindsByRoi = GetCStoSliceIndsByRoi(SrcRts, SOPuids)
        
        msg = f"There are no contours on slice {SrcSliceNum} in any ROI in "\
              + f"the Source RTS. There are {N} ROIs in the Source RTS:"
        
        for i in range(N):
            msg += f"\nROI {i+1} with name '{Names[i]}' has {NumByRoi[i]} "\
                   + f"contours that correspond to slices {C2SindsByRoi[i]}" 
        
        raise Exception(msg)
    
    """
    Determine whether the contours that make up the ROI(s) intersect with the
    Target image extent (14/01/2021). 
    """
    FracProp = ProportionOfRoisInExtent(PtsByCntByRoi=SrcPtsByCntByRoi,
                                        TrgDicomDir=TrgDcmDir,
                                        LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to check whether the Source ROI(s) intersect the '\
          + 'Target grid space.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
        
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the points in the ROI lie ',
              'within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the points in the ROI lie within the physical extent '\
              + 'of the Target image.'
        raise Exception(msg)
        return None
    
    
    """ Import the 3D images. """
    SrcIm = ImportImage(SrcDcmDir)
    TrgIm = ImportImage(TrgDcmDir)
    
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    #if UseCase in ['1', '2a', '3a']: # 08/02
    if UseCaseToApply in ['1', '2a']: # 08/02
        """ Direct copy of a contour without resampling or registration
        transformation. """
        
        from GeneralTools import ReplaceIndInC2SindsByRoi
        
        """ The contour on SrcSliceNum is to be copied to TrgSliceNum. Modify 
        the contour-to-slice index (from SrcSliceNum) to TrgSliceNum. """
        SrcC2SindsByRoi = ReplaceIndInC2SindsByRoi(C2SindsByRoi=SrcC2SindsByRoi,
                                                   IndToReplace=SrcSliceNum, 
                                                   ReplacementInd=TrgSliceNum)
        
        if LogToConsole:
                print(f'\nModfied SrcC2SindsByRoi = {SrcC2SindsByRoi}.')
                
        if TrgRtsFpath:
            """ Direct copy of a contour to an existing RTS. """
            
            """ An existing Target RTS is to be modified.  Since this is a  
            Direct copy, any existing contours in Target with ROI label  
            matching SrcRoiName are to be preserved. """
            
            """ Get the data of interest from the Target RTS. """
            TrgPtsByCntByRoi,\
            TrgC2SindsByRoi = GetRtsDataOfInterest(TrgRtsFpath, TrgSliceNum,
                                                   #SrcRoiName, TrgDcmDir,
                                                   TrgRoiName, TrgDcmDir,
                                                   LogToConsole)
            """ 05/03/21: Previously SrcRoiName was used as an input in 
            GetRtsDataOfInterest for Target, hence requiring that the Target
            ROI have the same name as the Source ROI. Now TrgRoiName has been
            added to the list of inputs allowing for the Target ROI to have a
            different name. """
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to parse the Target RTS file for '\
                  + 'existing contours within the ROI.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
            
            
            if LogToConsole:
                print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}.')
        else:
            TrgPtsByCntByRoi = None
            TrgC2SindsByRoi = None
        

    #if UseCase in ['1', '2a']:
        #""" Direct copy without resampling or registration transformation.
        """
        The contour will be copied to a different slice location, so the
        z-components of the coordinates will need to be modified. """
        
        """ Shift the out-of-plane (z) components of the points in
        SrcPtsByCntByRoi. """
        SrcPtsByCntByRoi = ZshiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi,
                                               NewZind=TrgSliceNum,
                                               DicomDir=SrcDcmDir)
        
        """ Shift the points to account for any difference in the origin of
        SrcIm and TrgIm, and to apply the necessary shift from SrcSliceNum to
        TrgSliceNum. """
        SrcPtsByCntByRoi,\
        SrcC2SindsByRoi = ShiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi, 
                                             C2SindsByRoi=SrcC2SindsByRoi, 
                                             SrcImage=SrcIm, 
                                             SrcSliceNum=SrcSliceNum, 
                                             TrgImage=TrgIm, 
                                             TrgSliceNum=TrgSliceNum, 
                                             RefImage=TrgIm, ShiftInX=False,  
                                             ShiftInY=False, ShiftInZ=True, 
                                             Fractional=False, 
                                             LogToConsole=LogToConsole)
        
        if LogToConsole:
                print(f'\nShifted SrcC2SindsByRoi = {SrcC2SindsByRoi}.')
        
        
        if TrgRtsFpath:
            """ An existing Target RTS is to be modified. 
            Concatenate TrgC2SindsByRoi and SrcC2SindsByRoi, and 
            TrgPtsByCntByRoi and SrcPtsByCntByRoi. """
            
            TrgC2SindsByRoi = [TrgC2SindsByRoi[r] + SrcC2SindsByRoi[r] for r in range(len(TrgC2SindsByRoi))]
            
            TrgPtsByCntByRoi = [TrgPtsByCntByRoi[r] + SrcPtsByCntByRoi[r] for r in range(len(TrgPtsByCntByRoi))]
        
        else:
            """ A Target RTS was not provided. """
            
            TrgC2SindsByRoi = deepcopy(SrcC2SindsByRoi)
            
            TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            
            if LogToConsole:
                print(f'\nSrcC2SindsByRoi = {SrcC2SindsByRoi}')
                print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}')
           
        
        """ Convert points to contour data format. """
        TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
    
    
    if UseCaseToApply == '2b':
        """ Relationship-preserving copy without resampling or registration
        transformation. """
        
        if False:
            """ TrgPtsByCntByRoi is simply SrcPtsByCntByRoi, and 
            TrgC2SindsByRoi is SrcC2SindsByRoi. """
            
            TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            TrgC2SindsByRoi = deepcopy(SrcC2SindsByRoi)
        
        
        
        
        if False:
            """ Shift the points to account for any difference in the origin (i.e.
            slice number 0) of SrcIm and TrgIm, and to apply the necessary shift in 
            the z components of the points, and to the indices in SrcC2SindsByRoi. """
            SrcPtsByCntByRoi,\
            SrcC2SindsByRoi = ShiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi, 
                                                 C2SindsByRoi=SrcC2SindsByRoi, 
                                                 SrcImage=SrcIm, 
                                                 SrcSliceNum=0, 
                                                 TrgImage=TrgIm, 
                                                 TrgSliceNum=0, 
                                                 RefImage=TrgIm, ShiftInX=False,  
                                                 ShiftInY=False, ShiftInZ=True, 
                                                 Fractional=False, 
                                                 LogToConsole=LogToConsole)
            
        
        #from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1
        from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2
        
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        #TrgC2SindsByRoi = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1(SrcC2SindsByRoi, 
        #                                                        SrcIm, TrgIm)
        
        TrgC2SindsByRoi = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcC2SindsByRoi, 
                                                                SrcIm, TrgIm)
        
        """ TrgPtsByCntByRoi is simply SrcPtsByCntByRoi. """
        TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            
        """ Convert points to contour data format. """
        TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
        
        
        if LogToConsole:
            print(f'\nSrcC2SindsByRoi = {SrcC2SindsByRoi}')
            print('TrgC2SindsByRoi (= shifted SrcC2SindsByRoi) =',
                  f'{SrcC2SindsByRoi}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        from ImageTools import GetImageInfo, ResampleLabImByRoi
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PtsByCntByRoi2PixArrByRoi
        from ConversionTools import PixArrByRoi2LabImByRoi
        from ConversionTools import PixArrByRoi2PtsByCntByRoi
        #from RtsTools import GetMaskFromContoursForSliceNum
        #from SegTools import ProportionOfSegsInExtent
        #from GeneralTools import NumOfListsAtDepthTwo
        
        
        #""" Import the 3D images. """
        #SrcIm = ImportImage(SrcDcmDir)
        #TrgIm = ImportImage(TrgDcmDir)
        
        """ Convert the RTS data of interest to a list of 3D pixel arrays for 
        each ROI in SrcPtsByCntByRoi. """
        SrcPixArrByRoi = PtsByCntByRoi2PixArrByRoi(SrcPtsByCntByRoi, SrcIm, 
                                                   LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to convert the Source points to pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        if PlotResults:
            """ Plot result of conversion from contours to segmentations: """
            import importlib
            import PlottingTools
            importlib.reload(PlottingTools)
            from PlottingTools import PlotResultOfMaskToContourConversion
            
            ListOfIms = [SrcIm]
            ListOfPtsByCntByRoi = [SrcPtsByCntByRoi]
            ListOfC2SindsByRoi = [SrcC2SindsByRoi]
            ListOfPixArrBySeg = [SrcPixArrByRoi]
            ListOfF2SindsBySeg = [SrcC2SindsByRoi]
            ListOfDicomDirs = [SrcDcmDir]
            ListOfPlotTitles = ['Src contours and conversion to masks']
            TxtToAddToFname = 'Conversion_of_Src_contours_to_masks'
            
            PlotResultOfMaskToContourConversion(ListOfIms=ListOfIms, 
                                                ListOfPtsByCntByRoi=ListOfPtsByCntByRoi,
                                                ListOfC2SindsByRoi=ListOfC2SindsByRoi,
                                                ListOfPixArrBySeg=ListOfPixArrBySeg,
                                                ListOfF2SindsBySeg=ListOfF2SindsBySeg,
                                                ListOfDicomDirs=ListOfDicomDirs,
                                                ListOfPlotTitles=ListOfPlotTitles, 
                                                AddTxt=None, ExportPlot=ExportPlot, 
                                                ExportDir=RtsPlotExportDir, 
                                                TxtToAddToFname=TxtToAddToFname, 
                                                LogToConsole=False)
        
        
        
        """ Determine whether the pixels that make up SrcPixArrByRoi intersect 
        with the Target image extent. 
        
        Note:
            This is redundant since this was already checked for the points.  
            Only useful as a confirmation that the conversion from points to 
            pixel arrays was correct.
        
        FracProp = ProportionOfSegsInExtent(PixArrBySeg=SrcPixArrByRoi, 
                                            F2SindsBySeg=SrcC2SindsByRoi,
                                            SrcDicomDir=SrcDcmDir,
                                            TrgDicomDir=TrgDcmDir,
                                            LogToConsole=LogToConsole)
        
        if LogToConsole:
            print(f'\n{round(FracProp*100, 1)}% of the voxels in the',
                  'segmentation to be copied lie within the physical extent',
                  'of the Target image.')
        
        if FracProp == 0:
            msg = 'None of the voxels in the segmentation to be copied lie '\
                  + 'within the physical extent of the Target image.'
            raise Exception(msg)
            return None
        """
        
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabImByRoi = PixArrByRoi2LabImByRoi(PixArrByRoi=SrcPixArrByRoi, 
                                               F2SindsByRoi=SrcC2SindsByRoi, 
                                               RefIm=SrcIm,
                                               LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to generate labelmap images from the Source '\
              + 'pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the frame-to-slice indices in the order associated with
        SrcLabImByRoi.
        Comment:
            Not sure the potential ordering difference between SrcC2SindsByRoi
            and SrcF2SindsByRoi is important for further steps.. """
        SrcF2SindsByRoi = []

        for SrcLabIm in SrcLabImByRoi:
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(SrcLabIm, LogToConsole=False)
            
            SrcF2SindsByRoi.append(F2Sinds)
    else:
        SrcPixArrByRoi = None
        SrcLabImByRoi = None
        SrcF2SindsByRoi = None
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        
        """ Resample the source labelmaps.
        Comment:
            Should use SrcF2SindsByRoi rather than SrcC2SindsByRoi?..."""
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        = ResampleLabImByRoi(LabImByRoi=SrcLabImByRoi, 
                             #F2SindsByRoi=SrcC2SindsByRoi,
                             F2SindsByRoi=SrcF2SindsByRoi,
                             SrcIm=SrcIm, 
                             TrgIm=TrgIm,
                             Interp=ResInterp,
                             PreResVariance=PreResVariance,
                             ApplyPostResBlur=ApplyPostResBlur,
                             PostResVariance=PostResVariance,
                             LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the interpolation that was used (the actual interpolation used 
        might have been different from the interpolation set) and the threshold
        used to re-binarise the resampled labelmap image (if applicable): """
        
        #print(f'\n\n\nMetadata keys in ResSrcLabImByRoi[0]:',
        #      ResSrcLabImByRoi[0].GetMetaDataKeys())
        
        for i in range(len(ResSrcLabImByRoi)):
            Interp = ResSrcLabImByRoi[i].GetMetaData("ResInterpUsed")
            
            ListOfInputs.append(f'ResInterpUsedForRoi{i} = {Interp}')
            DictOfInputs[f'ResInterpUsedForRoi{i}'] = Interp
            
            if Interp == 'BlurThenLinear':
                Thresh = ResSrcLabImByRoi[i].GetMetaData("PostResThreshUsed")
                
                ListOfInputs.append(f'PostResThreshUsedForRoi{i} = {Thresh}')
                DictOfInputs[f'PostResThreshUsedForRoi{i}'] = Thresh
    else:
        ResSrcLabImByRoi = None
        ResSrcPixArrByRoi = None
        ResSrcF2SindsByRoi = None
    
        
        
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        if SelxOrSitk == 'Selx':
            from SelxTools import RegisterImagesSelx, GetParamFromSelxImFilt
        else:
            from ImageTools import RegisterImagesSitk
        from ImageTools import TransformImWithListOfSitkTxs
        #from ImageTools import TransformLabImByRoi
        from ImageTools import TransformLabImByRoiWithListOfSitkTxs
        #from DroTools import CreateSitkTxFromTxMatrix_INCOMPLETE
        from DroTools import CreateSitkTxFromSpaDro
        from DroTools import CreatePreDefSitkTxFromDefDro
        from DroTools import CreateSitkTxFromDefDro
        from GeneralTools import AreListsEqualToWithinEpsilon
        
        if not Transform in ['rigid', 'affine', 'bspline']:
            msg = f'The chosen transform, {Transform}, is not one of the '\
                  + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
            
            raise Exception(msg)
        
        
        #print(f"type(Dro) within CopyRts() = {type(Dro)}\n")
        
        #if (Transform in ['rigid', 'affine'] and TxMatrix == None) \
        #or (Transform == 'bspline' and VectGridData == None):
        #if Dro == None:
        if Dro == None or (Dro != None and not UseDroForTx):
            """ The transform parameters required to register SrcIm to TrgIm 
            were not found from any matching subject assessors in XNAT, or
            a suitable DRO was found but UseDroForTx = False. 
            Image registration to be performed. """
            
            if Dro != None and not UseDroForTx:
                print('Image registration will be performed despite a',
                      'suitable DRO found on XNAT.\n')
            
            if LogToConsole == False:
                print(f'Performing image registration using {Transform}',
                      'transform...\n')
            
            times.append(time.time())
            
            """ Register SrcIm to TrgIm. """
            if SelxOrSitk == 'Selx':
                RegIm, SelxImFiltOrSitkTx, TxParamMapDict\
                = RegisterImagesSelx(FixIm=TrgIm, MovIm=SrcIm, 
                                     Transform=Transform,
                                     MaxNumOfIters=MaxIters,
                                     InitMethod=InitMethod,
                                     FixFidsFpath=TrgFidsFpath, 
                                     MovFidsFpath=SrcFidsFpath,
                                     FlipK=False, LogToConsole=LogToConsole)
                
                #""" Get the registration transform parameters: """
                #TxParams = list(TxParamMapDict['TransformParameters'])
                
                """ Get the registration transform parameters for the final
                transform parameter map: (10/05/21) """
                TxParams = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                                  Param='TransformParameters')
                
                
            if SelxOrSitk == 'Sitk':
                """ 08/06/21: Note that SelxImFiltOrSitkTx is a list of sitk 
                Transforms. """
                
                #RegIm, LandmarkTx, SelxImFiltOrSitkTx\
                RegIm, SelxImFiltOrSitkTx\
                = RegisterImagesSitk(FixIm=TrgIm, MovIm=SrcIm, 
                                     Transform=Transform, InitMethod=InitMethod,
                                     FixFidsFpath=TrgFidsFpath, 
                                     MovFidsFpath=SrcFidsFpath,
                                     LogToConsole=LogToConsole)
                
                """ Get the final registration transform parameters: """
                #TxParams = [str(item) for item in SelxImFiltOrSitkTx.GetParameters()]
                TxParams = [str(item) for item in SelxImFiltOrSitkTx[-1].GetParameters()]
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to register the Source image to the Target '\
                  + f'image using {Transform} transform.\n'
            ListOfTimings.append(msg)
            if LogToConsole == False:
                print(f'*{msg}')
            
            #print(f'TxParams = {TxParams}')
            
            #""" Add the row vector [0, 0, 0, 1] to TxParams to complete the 
            #Frame of Reference Transformation Matrix 
            #(http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
            #"""
            ##TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
            #TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
            #            str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
            #            str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
            #            str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
            
            
            from DicomTools import ImportDicom
            from PlottingTools import PlotResResults
            
            SrcDcmMod = ImportDicom(SrcDcmDir).Modality
            TrgDcmMod = ImportDicom(TrgDcmDir).Modality
            
            SrcScanId = DictOfInputs['SrcScanId']
            TrgScanId = DictOfInputs['TrgScanId']
            SrcExpLabel = DictOfInputs['SrcExpLabel']
            TrgExpLabel = DictOfInputs['TrgExpLabel']
            
            FixTitle = f'{TrgDcmMod} {TrgScanId} (fixed)' 
            MovTitle = f'{SrcDcmMod} {SrcScanId} (moving)' 
            
            ResTitle = f'{SrcDcmMod} {SrcScanId} {Transform} registered '\
                          + f'to {TrgDcmMod} {TrgScanId}'
            TxtToAddToFname = f'{SrcExpLabel}_{Transform}_reg_to_{TrgExpLabel}'
            
            ExportPlot = DictOfInputs['ExportPlot']
            RtsPlotExportDir = DictOfInputs['RtsPlotExportDir']
            
            """ This plot is useful but would require prior knowledge of
            which indices to plot. So would need to be moved further down.
            """
            PlotResResults(FixIm=TrgIm, MovIm=SrcIm, ResIm=RegIm, 
                           #FixInd=MrT2Im.GetSize()[2]//2,
                           #FixInd=MrT2Im_S - 25, MovInd=CtIm_S - 242,
                           FixInd=15, MovInd=223,
                           FixTitle=FixTitle, MovTitle=MovTitle, ResTitle=ResTitle,
                           ExportPlot=ExportPlot, ExportDir=RtsPlotExportDir, 
                           TxtToAddToFname=TxtToAddToFname)
            
            """ 
            04/06/2021:
            
            Due to issues with initialising a SimpleITK BSpline Transform
            using fiducials, the workflow for deformable registrations using
            fiducial pre-alignment consists of two resampling steps (as opposed 
            to one for deformable registrations not using fiducial pre-alignment
            or non-deformable registrations using fiducials or not). Hence
            create a list ListOfSitkTxs containing SelxImFiltOrSitkTx so
            that TransformLabImByRoi() can be looped for each item in the list.
            """
            #ListOfSitkTxs = [SelxImFiltOrSitkTx] # 08/06/21
            """ 08/06/21: Now SelxImFiltOrSitkTx is a list of Transforms:
            (this is messy - need to tidy up) """
            ListOfSitkTxs = SelxImFiltOrSitkTx
            
            DictOfInputs['SourceOfTransformation'] = 'ImageRegistration'
        
        else:
            """ The transform parameters required to register SrcIm to TrgIm 
            are obtained from the transform matrix, found from the DRO (subject 
            assessor) in XNAT. Create a SimpleITK transform from Dro:
            """
            
            #print(f'Transform = {Transform}\n')
            
            """ 04/06/21 Workaround described above. """
            ListOfSitkTxs = []
            
            if Transform in ['rigid', 'affine']:
                #SelxImFiltOrSitkTx = CreateSitkTxFromTxMatrix_INCOMPLETE(TrgIm, TxMatrix, 
                #                                              Transform)
                #
                #""" TxMatrix has additional bottom row [0 0 0 1] (see 
                #http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
                #
                #Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. 
                #Hence:
                #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
                #"""
                #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 06/05/21
                #
                #print(f'TxParams from TxMatrix (from DRO) = {TxParams}\n')
                #
                #print('\n')
                
                SitkTx = CreateSitkTxFromSpaDro(Dro, LogToConsole)
                
                TxParams = list(SitkTx.GetParameters())
                
                ListOfSitkTxs.append(SitkTx)
                
            else:
                #print(type(Dro))
                
                #print(Dro)
                
                PreRegTx = CreatePreDefSitkTxFromDefDro(Dro, LogToConsole)
                
                TxParams = list(PreRegTx.GetParameters())
                
                Identity = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
                
                if not AreListsEqualToWithinEpsilon(TxParams, Identity):
                    """ The pre-registration matrix is not the identity matrix
                    so need to resample/transform the label image using 
                    PreRegTx. """
                
                    ListOfSitkTxs.append(PreRegTx)
                
                SitkTx = CreateSitkTxFromDefDro(Dro, RefIm=TrgIm,
                                                LogToConsole=LogToConsole)
                
                ListOfSitkTxs.append(SitkTx)
            
            
            """ Since image registration was skipped, create RegIm by applying
            ListOfSitkTxs iteratively: """
            RegIm = TransformImWithListOfSitkTxs(MovIm=SrcIm, FixIm=TrgIm, 
                                                 ListOfSitkTxs=ListOfSitkTxs, 
                                                 Interp='Linear', 
                                                 LogToConsole=LogToConsole)
            
            if LogToConsole:
                from PlottingTools import PlotResResults
                
                """ This plot is useful but would require prior knowledge of
                which indices to plot. So would need to be moved further down.
                """
                PlotResResults(FixIm=TrgIm, MovIm=SrcIm, ResIm=RegIm, 
                               #FixInd=MrT2Im.GetSize()[2]//2,
                               #FixInd=MrT2Im_S - 25, MovInd=CtIm_S - 242,
                               FixInd=15, MovInd=223,
                               FixTitle='Fixed image', MovTitle='Moving image', 
                               ResTitle='Registered image')
            
            DictOfInputs['SourceOfTransformation'] = 'ParamsFromDRO'
            
            """ Set SelxImFiltOrSitkTx to be the last SimpleITK Transform in
            ListOfSitkTxs: """
            SelxImFiltOrSitkTx = ListOfSitkTxs[-1]
        
        
        """ Transform SrcLabImByRoi using the registration transformation. 
        
        Notes:
            
        Even though the data is being transformed the prefix 'Res' will be
        used since further operations below will need to be applied to
        either resampled or transformed data. 
            
        04/06/21: Workaround to issue with initialising a SimpleITK BSpline
        transform using fiducials involves a list (of potentially length 1)
        for the variable ListOfSitkTxs.  This is not ideal not only due to
        stacking of resampling errors but also due to use of a nearest
        neighbor interpolator.
        """
        
        #ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        #= TransformLabImByRoi(LabImByRoi=SrcLabImByRoi,
        #                      F2SindsByRoi=SrcF2SindsByRoi,
        #                      TrgIm=TrgIm,
        #                      SelxImFiltOrSitkTx=SelxImFiltOrSitkTx,
        #                      Interp=TxInterp,
        #                      ApplyPostTxBlur=ApplyPostTxBlur,
        #                      PostTxVariance=PostTxVariance,
        #                      ApplyPostTxBin=ApplyPostTxBin,
        #                      LogToConsole=LogToConsole)
        
        ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        = TransformLabImByRoiWithListOfSitkTxs(LabImByRoi=SrcLabImByRoi,
                                               F2SindsByRoi=SrcF2SindsByRoi,
                                               TrgIm=TrgIm,
                                               ListOfSitkTxs=ListOfSitkTxs,
                                               Interp=TxInterp,
                                               ApplyPostTxBlur=ApplyPostTxBlur,
                                               PostTxVariance=PostTxVariance,
                                               ApplyPostTxBin=ApplyPostTxBin,
                                               LogToConsole=LogToConsole)
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
        if ApplyPostTxBlur:
            """ Get the threshold used to re-binarise the transformed label 
            image: """
            for i in range(len(ResSrcLabImByRoi)):
                Thresh = ResSrcLabImByRoi[i].GetMetaData("PostTxThreshUsed")
                
                ListOfInputs.append(f'PostTxThreshUsedForRoi{i} = {Thresh}')
                DictOfInputs[f'PostTxThreshUsedForRoi{i}'] = Thresh
    else:
        RegIm = None
        SelxImFiltOrSitkTx = None
        TxParams = None
        ListOfSitkTxs = None
        
        
    
    if UseCaseToApply in ['3a', '4a', '5a']:
        """ Direct copy with resampling/registration and averaging of multi-
        framed pixel arrays.
        
        Note:
        
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled/registered labelmap image).  So if after resampling there are 
        more than one frame, the frames will be averaged or a logical OR 
        operation will be performed.
        """
        
        from GeneralTools import ShiftFramesInPixArrBySeg
        from PlottingTools import PlotPixArrBySeg
        
        R = len(ResSrcPixArrByRoi) # number of ROIs
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrByRoi[r].shape[0] for r in range(R)]
        
        MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
        
        if MaxNumOfFramesBySeg > 1:
            """ Reduce the number of frames in each pixel array to 1. """
            
            """ Perform averaging or OR operation? """
            #ReduceUsing = 'mean'
            ReduceUsing = 'OR'
            
            if LogToConsole:
                print(f'\nThere are {NumOfFramesBySeg} frames in each',
                      f'ROI/segment. Using {ReduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                PlotPixArrBySeg(ResSrcPixArrByRoi, ResSrcF2SindsByRoi, 
                                f'Prior to {ReduceUsing} operation')
            
            
            if ReduceUsing == 'mean':
                from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5
    
                ResSrcPixArrByRoi, ResF2SindsByRoi\
                = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
                                         F2SindsBySeg=ResSrcF2SindsByRoi,
                                         MakeBinary=True, 
                                         BinaryThresh=BinaryThresh,
                                         LogToConsole=LogToConsole)
            
            else:
                from GeneralTools import OrFrameOfPixArrBySeg
                
                #print(f'\n\n\nBefore running OrPixArrByRoi():')
                #print(f'ResSrcPixArrByRoi = {ResSrcPixArrByRoi}')
                #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
                #for r in range(len(ResSrcPixArrByRoi)):
                #    print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
                    
                ResSrcPixArrByRoi, ResF2SindsByRoi\
                = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
                                       F2SindsBySeg=ResSrcF2SindsByRoi,
                                       LogToConsole=LogToConsole)
                
                #print(f'\nAfter running OrPixArrByRoi():')
                #print(f'ResSrcPixArrByRoi = {ResSrcPixArrByRoi}')
                #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
                #for r in range(len(ResSrcPixArrByRoi)):
                #    print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
            
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrByRoi[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment.')
                
                PlotPixArrBySeg(ResSrcPixArrByRoi, ResSrcF2SindsByRoi, 
                                'After {ReduceUsing} operation')
        
        
        
        """ Shift the out-of-plane elements in ResSrcPixArrByRoi to account for
        the shift from SrcSliceNumber to TrgSliceNum. 
        Note:
            The function below was incorrectly applying in-plane shifts as well
            up until the added inputs ShiftInX / Y / Z was added on 10/02. 
        
        14/02: I think SrcSliceNum shouldn't be SrcSliceNum but instead
        ResSrcF2SindsByRoi[0][0], since this is the frame-to-slice number of the
        resampled/transformed labelmap.  The slice number of the original
        (Source) labelmap is not relevant.
        
        Likewise, the inputs SrcImage and TrgImage should both be TrgImage 
        since the resampled labelmap shares the TrgImage grid, and it will be
        shifted to another location, still within the TrgImage grid.
        """
        
        """ SrcSliceNum --> ResSrcSliceNum following resampling or 
        transformation. """
        ResSrcSliceNum = ResSrcF2SindsByRoi[0][0]
        
        if True:#LogToConsole:
            #print(f'\nSrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
            #      f'{ResSrcSliceNum}')
            #print(f'TrgSliceNum = {TrgSliceNum}')
            print(f'SrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
                  f'{ResSrcSliceNum} following resampling/transformation -->',
                  f'TrgSliceNum = {TrgSliceNum} for direct copy\n')
            print('ResSrcF2SindsByRoi prior to shifting of out-of-plane',
                  f'elements: {ResSrcF2SindsByRoi}\n')
        
        ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi, 
                                   F2SindsBySeg=ResSrcF2SindsByRoi, 
                                   ##SrcImage=SrcIm, 
                                   ##SrcSliceNum=SrcSliceNum,
                                   ##TrgImage=TrgIm,
                                   #SrcImage=ResSrcLabImByRoi[0], 
                                   #SrcSliceNum=ResSrcSliceNum,
                                   #TrgImage=ResSrcLabImByRoi[0],
                                   SrcImage=TrgIm, 
                                   SrcSliceNum=ResSrcSliceNum,
                                   TrgImage=TrgIm,
                                   TrgSliceNum=TrgSliceNum, 
                                   RefImage=TrgIm,
                                   ShiftInX=False,
                                   ShiftInY=False,
                                   ShiftInZ=True,
                                   Fractional=False,
                                   LogToConsole=LogToConsole)
        
        if LogToConsole:
            print('After running ShiftFramesInPixArrBySeg():')
            print(f'ResSrcF2SindsByRoi = {ResSrcF2SindsByRoi}')
            print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
            for r in range(len(ResSrcPixArrByRoi)):
                print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}\n')
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copying with resampling or
        registration-transformation. """
        
        """ Convert ResSrcPixArrByRoi to a list (for each ROI) of a list (for
        each contour) of ContourData, and a list (for each ROI) of a list (for
        each contour) of a list (for each point) of coordinates. """
        
        #print(f'\nBefore running PixArrByRoi2PtsByCntByRoi():')
        #print(f'ResSrcF2SindsByRoi = {ResSrcF2SindsByRoi}')
        #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
        #for r in range(len(ResSrcPixArrByRoi)):
        #    #print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
        #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
        
        """ Is Thresh = 0.5 below too high? """
            
        ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi, ResSrcC2SindsByRoi\
        = PixArrByRoi2PtsByCntByRoi(PixArrByRoi=ResSrcPixArrByRoi,
                                    F2SindsByRoi=ResSrcF2SindsByRoi,
                                    DicomDir=TrgDcmDir,
                                    Thresh=0.5,
                                    LogToConsole=LogToConsole)
        
        #print(f'\nAfter running PixArrByRoi2PtsByCntByRoi():')
        #print(f'ResSrcC2SindsByRoi = {ResSrcC2SindsByRoi}')
        #print(f'len(ResSrcPtsByCntByRoi) = {len(ResSrcPtsByCntByRoi)}')
        #for r in range(len(ResSrcPtsByCntByRoi)):
        #    C = len(ResSrcPtsByCntByRoi[r])
        #    print(f'   There are {C} contours in the {r}^th ROI')
        #    for c in range(C):
        #        P = len(ResSrcPtsByCntByRoi[r][c])
        #        print(f'      There are {P} points in the {c}^th contour')
            
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to convert the Source pixel arrays to points.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if LogToConsole:
            #print('\n\n', '-'*120)
            print('After converting ResPixArrToCopy to ContourData and',
                  'points, ResSrcC2SindsByRoi =')
            PrintIndsByRoi(ResSrcC2SindsByRoi)
            print('')
        
        
        if PlotResults:
            """ Plot result of conversion from segmentations to contours: """
            import importlib
            import PlottingTools
            importlib.reload(PlottingTools)
            from PlottingTools import PlotResultOfMaskToContourConversion
            
            ListOfIms = [SrcIm, TrgIm]
            ListOfPtsByCntByRoi = [SrcPtsByCntByRoi, ResSrcPtsByCntByRoi]
            ListOfC2SindsByRoi = [SrcC2SindsByRoi, ResSrcC2SindsByRoi]
            ListOfPixArrBySeg = [SrcPixArrByRoi, ResSrcPixArrByRoi]
            ListOfF2SindsBySeg = [SrcF2SindsByRoi, ResSrcF2SindsByRoi]
            ListOfDicomDirs = [SrcDcmDir, TrgDcmDir]
            ListOfPlotTitles = ['Src contours and conversion to masks',
                                'Tx masks and conversion to contours']
            TxtToAddToFname = 'Conversion_ResSrc_masks_to_contours'
            
            PlotResultOfMaskToContourConversion(ListOfIms=ListOfIms, 
                                                ListOfPtsByCntByRoi=ListOfPtsByCntByRoi,
                                                ListOfC2SindsByRoi=ListOfC2SindsByRoi,
                                                ListOfPixArrBySeg=ListOfPixArrBySeg,
                                                ListOfF2SindsBySeg=ListOfF2SindsBySeg,
                                                ListOfDicomDirs=ListOfDicomDirs,
                                                ListOfPlotTitles=ListOfPlotTitles, 
                                                AddTxt=None, ExportPlot=ExportPlot, 
                                                ExportDir=RtsPlotExportDir, 
                                                TxtToAddToFname=TxtToAddToFname, 
                                                LogToConsole=False)
        
    
        if UseCaseToApply in ['3a', '4a', '5a']:
            """ Direct copying with resampling/transformation. """
            
            if TrgRtsFpath:
                """ An existing Target RTS is to be modified. 
                Concatenate TrgC2SindsByRoi and ResSrcC2SindsByRoi, and
                TrgPtsByCntByRoi and ResSrcPtsByCntByRoi. """
                
                TrgC2SindsByRoi = [TrgC2SindsByRoi[r] + ResSrcC2SindsByRoi[r] for r in range(len(TrgC2SindsByRoi))]
                
                TrgPtsByCntByRoi = [TrgPtsByCntByRoi[r] + ResSrcPtsByCntByRoi[r] for r in range(len(TrgPtsByCntByRoi))]
                
            else:
                """ A Target RTS was not provided. """
                
                TrgC2SindsByRoi = deepcopy(ResSrcC2SindsByRoi)
                
                TrgPtsByCntByRoi = deepcopy(ResSrcPtsByCntByRoi)
                
            
            """ Convert points to contour data format. """
            TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
        
        else:
            """ UseCase is '3b', '4b' or '5b'. """
            
            TrgPtsByCntByRoi = deepcopy(ResSrcPtsByCntByRoi)
            
            TrgCntDataByCntByRoi = deepcopy(ResSrcCntDataByCntByRoi)
            
            TrgC2SindsByRoi = deepcopy(ResSrcC2SindsByRoi)
    else:
        ResSrcPtsByCntByRoi = None
        ResSrcCntDataByCntByRoi = None
        ResSrcC2SindsByRoi = None
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print('Prior to running CreateRts():')
        print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}')
        for r in range(len(TrgPtsByCntByRoi)):
            C = len(TrgPtsByCntByRoi[r])
            print(f'   There are {C} in the {r}^th ROI')
            for c in range(len(TrgPtsByCntByRoi[r])):
                P = len(TrgPtsByCntByRoi[r][c])
                print(f'      There are {P} pts in the {c}^th contour')
        print('')
    
    
    """ Create the Target RTS. """
    
    print('*Creating new Target RTS object...\n')
    
    TrgRts = CreateRts(SrcRtsFpath, TrgRtsFpath, TrgCntDataByCntByRoi,
                       TrgPtsByCntByRoi, TrgC2SindsByRoi, TrgDcmDir,
                       TxtToAddToTrgRoiName, LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to create the RTS object.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Create a DICOM Registration Object (if applicable). """
    #if UseCaseToApply in ['5a', '5b'] and (Dro != None and OverwriteDro
    #                                       or Dro == None):
    if UseCaseToApply in ['5a', '5b']:
        from DroTools import CreateSpaDroFromSitkTx, CreateDefDroFromBsplineTx
        
        # ContentDescription:
        """ Consider adding this as an input variable for the user to specify.
        """
        ContentDesc = ''
        
        if Transform in ['rigid', 'affine']:
            """ Create a Spatial Registration DRO. """
            #Dro = CreateSpaDro(SrcDcmDir, TrgDcmDir, TxParams, ContentDesc, 
            #                   LogToConsole)
            Dro = CreateSpaDroFromSitkTx(SrcDcmDir, TrgDcmDir, 
                                         SelxImFiltOrSitkTx, 
                                         ContentDesc, LogToConsole)
        else:
            #""" Get the bspline grid related registration transform parameters 
            #for the final transform parameter map: (10/05/21) """
            #GridOrig = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1, 
            #                                  'GridOrigin')
            #
            #GridDir = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                 'GridDirection')
            #
            #GridDims = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                  'GridSize')
            #GridDims = [int(item) for item in GridDims]
            #
            #GridRes = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                 'GridSpacing')
            #GridRes = [float(item) for item in GridRes]
            #
            #VectGridData = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                      'TransformParameters')
            #VectGridData = [float(item) for item in VectGridData]
            #
            #""" Get the non-deformable registration transform parameters 
            #for the second-last transform parameter map: (12/05/21) """
            #TxParams = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -2,
            #                                  'TransformParameters')
            #if TxParams != None:
            #    TxParams = [float(item) for item in TxParams]
            
            """ Create a Deformable Spatial Registration DRO. """
            #Dro = CreateDefDro(SrcDcmDir, TrgDcmDir, GridOrig, GridDir, 
            #                   GridDims, GridRes, VectGridData, TxParams,
            #                   ContentDesc, LogToConsole)
            #Dro = CreateDefDroFromBsplineTx(SrcDcmDir, TrgDcmDir, 
            #                                BsplineTx=SelxImFiltOrSitkTx, 
            #                                PreRegTx=LandmarkTx, 
            #                                Description=ContentDesc, 
            #                                LogToConsole=LogToConsole)
            """ Note:
                The first Transform in ListOfSitkTxs is the landmark-based
                Transform, the second/last one is the final BSpline Transform.
            """
            Dro = CreateDefDroFromBsplineTx(SrcDcmDir, TrgDcmDir, 
                                            BsplineTx=ListOfSitkTxs[-1], 
                                            PreRegTx=ListOfSitkTxs[0], 
                                            Description=ContentDesc, 
                                            LogToConsole=LogToConsole)
        
        
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to create the DRO.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
    else:
        Dro = None
    
    
    if LogToConsole:
        print('\nEnd of CopyRts().')
        print('-'*120)
        
    #return TrgRts, Dro, DictOfInputs, ListOfInputs, ListOfTimings
    ##return TrgRts, TrgPtsByCntByRoi, TrgC2SindsByRoi, ListOfTimings
    return SrcDcmDir, TrgDcmDir, SrcIm, TrgIm, SrcRts,\
           SrcPtsByCntByRoi, SrcC2SindsByRoi, TrgPtsByCntByRoi, TrgC2SindsByRoi,\
           SrcPixArrByRoi, SrcLabImByRoi, SrcF2SindsByRoi,\
           RegIm, ListOfSitkTxs, SelxImFiltOrSitkTx, TxParams,\
           ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi,\
           ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi, ResSrcC2SindsByRoi,\
           TrgRts, Dro, DictOfInputs, ListOfInputs, ListOfTimings