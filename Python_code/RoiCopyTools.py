# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:12:23 2020

@author: ctorti
"""


"""
ROI copying tools



Make a Direct or Relationship-preserving copy of a single contour/segmentation 
in a given ROI/segment in a Source RTS/SEG for a given slice to one or more
contours/segmentations in a new (if making a relationship-preserving copy) or 
existing (if making a direct copy) ROI/segment in a new Target RTS/SEG either 
using the original Target RTS/SEG (if one has been provided) or using the 
Source RTS/SEG (if a Target RTS/SEG was not provided) as a template.  

For Direct copies the in-plane coordinates (currently the (x,y) coords) will be
preserved only, whilst for a Relationship-preserving copy, all coordinates are 
spatially preserved.

There are five possible Main Use Cases (indicated by the numbers 1 to 5)
and two Sub-cases (indicated by the letters a or b) that apply for a subset
of the Main Use Cases that are accommodated.  The characters i and ii will 
further split Main Use Case #3 into two: i for a sub-case that considers a
copy from a high to low voxel sampling image, and ii for the sub-case that
considers a copy from a low to a high voxel sampling image.  In all cases a 
single contour/segmentation will be copied from the Source RTS/SEG within the 
ROI/segment collection containing the ROIName/SegmentLabel that matches the
characters in a search string (FromSearchString):

1a. Direct copy of the Source contour/segmentation on slice SrcSliceNum to a 
    new/existing ROI/segment for slice TrgSliceNum in Target. 
    The Source and Target images are the same. 
    
2a. Direct copy of the Source contour/segmentation on slice SrcSliceNum to a 
    new/existing ROI for slice TrgSliceNum in Target. 
    The Source and Target images have the same FOR, IOP, IPP, and VS.

2b. Relationship-preserving copy of the Source contour/segmentation on slice 
    SrcSliceNum to a new ROI/segment for slice TrgSliceNum in Target. 
    The Source and Target images have the same FOR, IOP, IPP, and VS. 

3a. Direct copy of the Source contour/segmentation on slice SrcSliceNum to a 
    new/existing ROI/segment for slice TrgSliceNum in Target. 
    The Source and Target images have the same FOR and IOP but have 
    different VS (they will likely also have different IPP).
    
3b. Relationship-preserving copy of the Source contour/segmentation on slice 
    SrcSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have the same FOR and IOP but have 
    different VS (they will likely also have different IPP).

4.  Relationship-preserving copy of the Source contour/segmentation on slice 
    SrcSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have the same FOR but have different IOP 
    and IPP, and possibly different VS.
   
5.  Relationship-preserving copy of the Source contour/contour on slice 
    SrcSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have different FOR, likely different IPP and 
    IOP, and possibly different VS.
    
FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
IOP = ImageOrientationPatient, VS = Voxel spacings
"""





"""
******************************************************************************
******************************************************************************
CHECK VALIDITY OF INPUTS TO COPYROI
******************************************************************************
******************************************************************************
"""


def CheckValidityOfInputs(SrcRoi, SrcRoiName, SrcSliceNum, TrgRoi, TrgRoiName, 
                          TrgSliceNum, LogToConsole=False):
    """
    Establish whether the inputs SrcRoiName, SrcSliceNum, TrgRoiName and
    TrgSliceNum are valid combinations. If not raise exception.
    
    Inputs:
    ******
    
    SrcRoi : Pydicom object
        The Source RTS/SEG object.
        
    SrcRoiName : string
        The Source ROIName or SegmentLabel.
    
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
    
    TrgRts : Pydicom object or None
        The Target RTS/SEG object if applicable; None otherwise.
        
    TrgRoiName : string or None
        The Target ROIName or SegmentLabel if applicable; None otherwise.
        
    TrgSliceNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.  
                        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.

              
    Outputs:
    ******
        
    IsValid : boolean
        True if the combinations of inputs are valid.
    """
    
    """
    if TrgRoi:
        from DicomTools import IsSameModalities
        
        SrcModality = SrcRoi.Modality
        TrgModality = TrgRoi.Modality
        
        if not IsSameModalities(SrcRoi, TrgRoi):
            msg = f"The Source ({SrcModality}) and Target ({TrgModality}) " \
                  + "modalities are different. This is not valid."
            
            raise Exception(msg)
    """
    
    from DicomTools import GetRoiLabels
    
    #print(f'\nTrgSliceNum = {TrgSliceNum}')
    
    msg = ''
              
    #if SrcRoiName == None:
    if not SrcRoiName:
        #if SrcSliceNum != None:
        if SrcSliceNum:
            msg += f"SrcRoiName = {SrcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut SrcSliceNum (= {SrcSliceNum}) is "\
                   + "specified (i.e. != None). "
        
        #if TrgRoiName != None:
        if TrgRoiName:
            msg += f"SrcRoiName = {SrcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut TrgRoiName (= {TrgRoiName}) is specified "\
                   + "(i.e. != None). "
        
        #if TrgSliceNum != None:
        if TrgSliceNum:
            msg += f"SrcRoiName = {SrcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut TrgSliceNum (= {TrgSliceNum}) is "\
                   + "specified (i.e. != None). "
        
    else:
        #if SrcSliceNum == None and TrgSliceNum != None:
        if not SrcSliceNum and TrgSliceNum:
            msg += f"SrcSliceNum = {SrcSliceNum} => All contours/segmentations"\
                   + f"in the ROI/segment with name/label '{SrcRoiName}' are "\
                   + f"to be copied. \nBut TrgSliceNum (= {TrgSliceNum}) is "\
                   + f"specified. "
        
        #if SrcSliceNum != None and TrgSliceNum == None:
        if SrcSliceNum and not TrgSliceNum:
            msg += f"SrcSliceNum = {SrcSliceNum} => The contour(s)/"\
                   + "segmentation(s) for the specified DICOM slice in the "\
                   + f"ROI/segment with name/label '{SrcRoiName}' are to be "\
                   + f"copied. \nBut TrgSliceNum (= {TrgSliceNum}) is not "\
                   + "specified. "
        
        if TrgRoi != None:
            """ Get the names/labels of all ROIs/segments in TrgRoi: """
            TrgRoiNames = GetRoiLabels(TrgRoi)
            
            T = len(TrgRoiNames)
            
            #if T > 1 and TrgRoiName == None:
            if T > 1 and not TrgRoiName:
                msg += f"SrcRoiName (= {SrcRoiName}) has been specified but "\
                       + f"TrgRoiName (= {TrgRoiName}) has not been specified"\
                       + f" and there are {T} ROIs/segments in TrgRoi. "
        
    if msg:
        msg += "This is not a valid combination of inputs."
            
        raise Exception(msg)
    
    return True





"""
******************************************************************************
******************************************************************************
WHICH USE CASE APPLIES?
******************************************************************************
******************************************************************************
"""

def WhichUseCase(SrcSliceNum, TrgSliceNum, SrcDcmDir, TrgDcmDir, ForceReg, 
                 LogToConsole=False):
    """
    Determine which of 5 Use Cases applies.
    
    Inputs:
    ******
        
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
    
    
    TrgSliceNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.

              
    Outputs:
    ******
        
    UseCaseThatApplies : string
        String that denotes the Use Case that applies.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
    """
    
    from DicomTools import GetDicomUids
    from ImageTools import GetImageAttributes
    from GeneralTools import AreListsEqualToWithinEpsilon#, PrintTitle
    
    
    """ Get the DICOM UIDs. """
    SrcStudyuid, SrcSeriesuid, SrcFORuid,\
    SrcSOPuids = GetDicomUids(DicomDir=SrcDcmDir)
    
    TrgStudyuid, TrgSeriesuid, TrgFORuid,\
    TrgSOPuids = GetDicomUids(DicomDir=TrgDcmDir)
        
    if TrgSliceNum:
        if TrgSliceNum > len(TrgSOPuids) - 1:
            msg = f'The input argument "TrgSliceNum" = {TrgSliceNum} exceeds '\
                  + f'the number of Target DICOMs ({len(TrgSOPuids)}).'
            
            raise Exception(msg)
    
    
    """ Get the Image Attributes for Source and Target using SimpleITK. """
    SrcSize, SrcSpacing, SrcST, SrcIPPs, SrcDirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=SrcDcmDir, Package='sitk')
    
    TrgSize, TrgSpacing, TrgST, TrgIPPs, TrgDirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=TrgDcmDir, Package='sitk')
    
    
    """ Are the Source and Target part of the same study, series and do they 
    have the same SOP UIDs? """
    if SrcStudyuid == TrgStudyuid:
        SameStudy = True
    else:
        SameStudy = False
        
    if SrcSeriesuid == TrgSeriesuid:
        SameSeries = True
    else:
        SameSeries = False
        
    if SrcFORuid == TrgFORuid:
        SameFOR = True
    else:
        SameFOR = False
        
    if SrcSOPuids == TrgSOPuids:
        SameSOPs = True
    else:
        SameSOPs = False
    
    """ Are the Source and Target spacings, directions and origins equal to  
    within 1e-06? """
    SameSpacings = AreListsEqualToWithinEpsilon(SrcSpacing, TrgSpacing)
    SameDirs = AreListsEqualToWithinEpsilon(SrcDirs, TrgDirs)
    SameOrigins = AreListsEqualToWithinEpsilon(SrcIPPs[0], TrgIPPs[0])
    SameST = AreListsEqualToWithinEpsilon([SrcST], [TrgST])
    
    if SrcSize == TrgSize:
        SameSize = True
    else:
        SameSize = False
    
        
    """ Check which case follows. """
    if SameFOR:
        if SameDirs:
            """ 
            Comment 11/02:
                
            It's not entirely necessary to resample the Source image to the 
            Target domain just because their origins differ. It depends on
            whether the slice in Source coincide with slices in Target.
            
            Even if the slices don't coincide, an approximation can be made as
            follows.  When dealing contours, the in-plane coordinates don't 
            need to be changed, and the destination slice can be determined by 
            matching on the out-of-plane (z) coordinates as best as possible. 
            When dealing with segmentations, the destination slice can be 
            determined by matching on the z components of the IPPs as best as
            possible.  Once the destination slice number is known, the 
            necessary shift in pixels for the in-plane components can be 
            applied accordingly. 
            
            Perhaps the approximation can be justified if the difference 
            between the planes in Source and Target are close to an integer.
            If the difference is large, the Source image should be resampled to
            the Target image.  For now I will treat UseCase 2 data sets as 
            UseCase 3 if their origins are not the same. Further complexity
            (taking into account the error in the above approximation) can be
            considered for future work. """
            if SameSpacings and SameOrigins:
            #if SameSpacings:
                if SameSeries:
                    UseCaseThatApplies = '1' # Direct copy
                else:
                    if isinstance(TrgSliceNum, int):
                        UseCaseThatApplies = '2a' # Direct copy
                    else:
                        UseCaseThatApplies = '2b' # Relationship-preserving copy
            else:
                """ Comment 11/02:
                    
                Arguably resampling is only required if the spacings along any
                direction are different when dealing with segmentations.  When
                dealing with contours, resampling is probably unnecessary if
                the spacings are different only in the out-of-plane (z)
                direction. """
                if isinstance(TrgSliceNum, int):
                    UseCaseThatApplies = '3a' # Direct copy
                else:
                    UseCaseThatApplies = '3b' # Relationship-preserving copy
        else:
            if isinstance(TrgSliceNum, int):
                UseCaseThatApplies = '4a' # Direct copy
            else:
                UseCaseThatApplies = '4b' # Relationship-preserving copy
    else:
        if isinstance(TrgSliceNum, int):
            UseCaseThatApplies = '5a' # Direct copy
        else:
            UseCaseThatApplies = '5b' # Relationship-preserving copy
    
        
    """ Print comparisons to the console. """
    if LogToConsole:
        from ImageTools import ImportImage, GetImageExtent
        
        SrcIm = ImportImage(SrcDcmDir)
        TrgIm = ImportImage(TrgDcmDir)
        
        """ Get the image extents. """
        SrcImExtent = GetImageExtent(SrcIm)
        TrgImExtent = GetImageExtent(TrgIm)
        
        """ Are the image extents equal to within 1e-06? """
        SameExtents = AreListsEqualToWithinEpsilon(SrcImExtent, TrgImExtent)
    
        #CompareImageAttributes(SrcDcmDir, TrgDcmDir)
        
        print('\n\n', '-'*120)
        print('Results of running WhichUseCase():')
        print('The Source and Target are in the/have:')
        
        if SameStudy:
            print('   * same study')
        else:
            print('   * different studies')
        
        if SameFOR:
            print('   * same Frame of Reference')
        else:
            print('   * different Frames of Reference')
            
        if SameSeries:
            print('   * same series')
        else:
            print('   * different series')
            
        if SameSOPs:
            print('   * same DICOMs')
        else:
            print('   * different DICOMs')
        
        if SameSize:
            print(f'   * same image size\n      Source/Target: {SrcSize}')
        else:
            print(f'   * different image sizes:\n      Source: {SrcSize}',
                  f'\n      Target: {TrgSize}')
        
        if SameSpacings:
            print(f'   * same voxel sizes\n      Source/Target: {SrcSpacing}')
        else:
            print(f'   * different voxel sizes:\n      Source: {SrcSpacing}',
                  f'\n      Target: {TrgSpacing}')
        
        if SameST:
            print(f'   * same slice thickness\n      Source/Target: {SrcST}')
        else:
            print(f'   * different slice thickness:\n      Source: {SrcST}',
                  f'\n      Target: {TrgST}')
        
        if SameExtents:
            print(f'   * same image extents\n      Source/Target: {SrcImExtent}')
        else:
            print(f'   * different image extents:\n      Source: {SrcImExtent}',
                  f'\n      Target: {TrgImExtent}')
        
        if SameOrigins:
            print(f'   * same origin\n      Source/Target: {SrcIPPs[0]}')
        else:
            print(f'   * different origins:\n      Source: {SrcIPPs[0]}',
                  f'\n      Target: {TrgIPPs[0]}')
            
        if SameDirs:
            print('   * same patient orientation\n',
                  f'      Source/Target: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},\n',
                  f'                      {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},\n',
                  f'                      {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]')
        else:
            print('   * different patient orientations:\n',
                  f'      Source: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},\n',
                  f'               {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},\n',
                  f'               {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]\n',
                  f'      Target: [{TrgDirs[0]}, {TrgDirs[1]}, {TrgDirs[2]},\n',
                  f'               {TrgDirs[3]}, {TrgDirs[4]}, {TrgDirs[5]},\n',
                  f'               {TrgDirs[6]}, {TrgDirs[7]}, {TrgDirs[8]}]')
        
    
    if True:#LogToConsole:
        #PrintTitle(f'Case {UseCase} applies.')
        if ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
            UseCaseToApply = UseCaseThatApplies.replace('3', '5').replace('4', '5')
            
            msg = f'\n*UseCase {UseCaseThatApplies} applies but since '\
                  + f'ForceReg = {ForceReg}, UseCase {UseCaseToApply} will '\
                  + f'be applied.'
            
            print(msg)
        else:
            UseCaseToApply = UseCaseThatApplies
            
            print(f'\n*UseCase {UseCaseThatApplies} applies.')
            
        print('-'*120)
        
    return UseCaseThatApplies, UseCaseToApply






"""
******************************************************************************
******************************************************************************
COPY A CONTOUR / ROI / RTSTRUCT
******************************************************************************
******************************************************************************
"""

def CopyRts(SrcRtsFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgRtsFpath=None, TrgRoiName=None, TrgSliceNum=None, 
            ResInterp='BlurThenLinear', PreResVar=(1,1,1), #PostResThresh=0.75, 
            ForceReg=False, Tx='affine', TxMaxIters='512',
            TxInterp='NearestNeighbor', 
            ApplyPostTxBlur=True, PostTxVar=(2,2,2), 
            ApplyPostTxBin=True, #ThreshPostTx=0.05, 
            TxtToAddToTrgRoiName='', LogToConsole=False,
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
    
    PreResVar : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to perform binary thresholding following
    #    resampling using a non-label (e.g. linear) interpolator. Example use is 
    #    on a labelmap image if a Gaussian blurring + linearly resampling 
    #    approach is taken to combat aliasing effects.
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    Tx : string (optional; 'affine' by default)
        The transformation to used for image registration (if applicable).  
        Acceptable values are:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    TxMaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Tx. If != 'default' it must be a string 
        representation of an integer.
        
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
        
    PostTxVar : tuple of floats (optional; (2,2,2) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
        
    #ThreshPostTx : float (optional; 0.05 by default)
    #    The threshold level used to perform binary thresholding after 
    #    registration transformation.  
    
    TxtToAddToTrgRoiName : string (optional, '' by default)
        String of text to add to the ROI Name of the new Target RTS.
                           
    LogToConsole : boolean (default False)
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
    import RtsTools
    importlib.reload(RtsTools)
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    #from RtsTools import GetPtsByCntByRoi
    from RtsTools import GetRtsDataOfInterest, ProportionOfRoisInExtent
    from RtsTools import CreateRts
    from GeneralTools import UniqueItems, PrintIndsByRoi#, PrintTitle
    from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    from GeneralTools import ShiftPtsByCntByRoi
    #from GeneralTools import ShiftFrame
    #from DicomTools import GetRoiNum
    #from ConversionTools import Points2ContourData
    from ConversionTools import PtsByCntByRoi2CntDataByCntByRoi
    from ImageTools import ImportImage
    
    """ Start timing. """
    times = []
    times.append(time.time())
    #""" Log timing messages. """
    #TimingMsgs = []
    
    
    SrcRts = dcmread(SrcRtsFpath)
    if TrgRtsFpath:
        TrgRts = dcmread(TrgRtsFpath)
    else:
        TrgRts = None
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print(f'Running CopyRts():')
        
        SrcRoiLabels = GetRoiLabels(SrcRts)
        
        print(f'\nSrcRoiLabels = {SrcRoiLabels}')
        
        if TrgRtsFpath:
            TrgRoiLabels = GetRoiLabels(TrgRts)
        
            print(f'\nTrgRoiLabels = {TrgRoiLabels}')
    
    
    """ Check the validity of the combination of inputs. """
    CheckValidityOfInputs(SrcRts, SrcRoiName, SrcSliceNum, TrgRts, TrgSliceNum,
                          LogToConsole)
    
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
    
    
    """ Raise exception if there is no data of interest from the Source RTS."""
    if not UniqueItems(SrcC2SindsByRoi):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        from RtsTools import GetNumOfContoursByRoi, GetCStoSliceIndsByRoi
        
        SrcRts = dcmread(SrcRtsFpath)
        
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
            TrgC2SindsByRoi = GetRtsDataOfInterest(TrgRtsFpath, SrcSliceNum,
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
            print(f'TrgC2SindsByRoi (= shifted SrcC2SindsByRoi) =',
                  f'{SrcC2SindsByRoi}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        from ImageTools import GetImageInfo, ResampleLabmapImByRoi
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PtsByCntByRoi2PixArrByRoi
        from ConversionTools import PixArrByRoi2LabmapImByRoi
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
        SrcLabmapImByRoi = PixArrByRoi2LabmapImByRoi(PixArrByRoi=SrcPixArrByRoi, 
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
        SrcLabmapImByRoi.
        Comment:
            Not sure the potential ordering difference between SrcC2SindsByRoi
            and SrcF2SindsByRoi is important for further steps.. """
        SrcF2SindsByRoi = []

        for SrcLabmapIm in SrcLabmapImByRoi:
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(SrcLabmapIm, LogToConsole=False)
            
            SrcF2SindsByRoi.append(F2Sinds)
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        
        """ Resample the source labelmaps.
        Comment:
            Should use SrcF2SindsByRoi rather than SrcC2SindsByRoi?..."""
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabmapImByRoi, ResSrcPixArrByRoi,\
        ResSrcF2SindsByRoi = ResampleLabmapImByRoi(LabmapImByRoi=SrcLabmapImByRoi, 
                                                   #F2SindsByRoi=SrcC2SindsByRoi,
                                                   F2SindsByRoi=SrcF2SindsByRoi,
                                                   SrcImage=SrcIm, 
                                                   TrgImage=TrgIm,
                                                   Interpolation=ResInterp,
                                                   PreResVar=PreResVar,
                                                   #PostResThresh=PostResThresh,
                                                   LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the interpolation that was used (the actual interpolation used 
        might have been different from the interpolation set) and the threshold
        used to re-binarise the resampled labelmap image (if applicable): """
        
        #print(f'\n\n\nMetadata keys in ResSrcLabmapImByRoi[0]:',
        #      ResSrcLabmapImByRoi[0].GetMetaDataKeys())
        
        for i in range(len(ResSrcLabmapImByRoi)):
            Interp = ResSrcLabmapImByRoi[i].GetMetaData("ResInterpUsed")
            
            ListOfInputs.append(f'ResInterpUsedForRoi{i} = {Interp}')
            DictOfInputs[f'ResInterpUsedForRoi{i}'] = Interp
            
            if Interp == 'BlurThenLinear':
                Thresh = ResSrcLabmapImByRoi[i].GetMetaData("PostResThreshUsed")
                
                ListOfInputs.append(f'PostResThreshUsedForRoi{i} = {Thresh}')
                DictOfInputs[f'PostResThreshUsedForRoi{i}'] = Thresh
    
        
        
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        from ImageTools import RegisterImages
        from ImageTools import TransformLabmapImByRoi
        
        if not Tx in ['rigid', 'affine', 'bspline']:
            msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
                  + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
            
            raise Exception(msg)
        
        if LogToConsole == False:
            print(f'Performing image registration using {Tx} transform...\n')
        
        times.append(time.time())
        
        """ Register SrcIm to TrgIm. """
        RegIm, RegImFilt = RegisterImages(FixIm=TrgIm, MovIm=SrcIm, Tx=Tx,
                                          MaxNumOfIters=TxMaxIters, 
                                          LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to register the Source image to the Target '\
              + f'image using {Tx} transform.\n'
        ListOfTimings.append(msg)
        if LogToConsole == False:
            print(f'*{msg}')
        
        
        """ Transform SrcLabmapImByRoi using the registration transformation. 
        
        Note:
            Even though the data is being transformed the prefix 'Res' will be
            used since further operations below will need to be applied to
            either resampled or transformed data. """
        
        ResSrcLabmapImByRoi, ResSrcPixArrByRoi,\
        ResSrcF2SindsByRoi = TransformLabmapImByRoi(LabmapImByRoi=SrcLabmapImByRoi,
                                                    F2SindsByRoi=SrcF2SindsByRoi,
                                                    RegImFilt=RegImFilt,
                                                    Interpolation=TxInterp,
                                                    ApplyPostTxBlur=ApplyPostTxBlur,
                                                    PostTxVar=PostTxVar,
                                                    ApplyPostTxBin=ApplyPostTxBin,
                                                    #ThreshPostTx=ThreshPostTx,
                                                    LogToConsole=LogToConsole)
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
        """ Get the threshold used to re-binarise the transformed labelmap 
        image: """
        for i in range(len(ResSrcLabmapImByRoi)):
            Thresh = ResSrcLabmapImByRoi[i].GetMetaData("PostTxThreshUsed")
            
            ListOfInputs.append(f'PostTxThreshUsedForRoi{i} = {Thresh}')
            DictOfInputs[f'PostTxThreshUsedForRoi{i}'] = Thresh
        
        
    
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
    
                ResSrcPixArrByRoi,\
                ResF2SindsByRoi = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
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
                    
                ResSrcPixArrByRoi,\
                ResF2SindsByRoi = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
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
        
        ResSrcPixArrByRoi,\
        ResSrcF2SindsByRoi = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi, 
                                                      F2SindsBySeg=ResSrcF2SindsByRoi, 
                                                      ##SrcImage=SrcIm, 
                                                      ##SrcSliceNum=SrcSliceNum,
                                                      ##TrgImage=TrgIm,
                                                      #SrcImage=ResSrcLabmapImByRoi[0], 
                                                      #SrcSliceNum=ResSrcSliceNum,
                                                      #TrgImage=ResSrcLabmapImByRoi[0],
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
            print(f'After running ShiftFramesInPixArrBySeg():')
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
            
        ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi,\
        ResSrcC2SindsByRoi = PixArrByRoi2PtsByCntByRoi(PixArrByRoi=ResSrcPixArrByRoi,
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
            print(f'After converting ResPixArrToCopy to ContourData and',
                  f'points, ResSrcC2SindsByRoi =')
            PrintIndsByRoi(ResSrcC2SindsByRoi)
            print('')
        
    
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
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print(f'Prior to running CreateRts():')
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
    
    if LogToConsole:
        print('\nEnd of CopyRts().')
        print('-'*120)
        
    return TrgRts, DictOfInputs, ListOfInputs, ListOfTimings
    #return TrgRts, TrgPtsByCntByRoi, TrgC2SindsByRoi, ListOfTimings






"""
******************************************************************************
******************************************************************************
COPY A SEGMENTATION / SEGMENT / SEG
******************************************************************************
******************************************************************************
"""

def CopySeg(SrcSegFpath, SrcSliceNum, SrcSegLabel, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgSegFpath=None, TrgSegLabel=None, TrgSliceNum=None, 
            ResInterp='BlurThenLinear', PreResVar=(1,1,1), #PostResThresh=0.75, 
            ForceReg=False, Tx='affine', TxMaxIters='512',
            TxInterp='NearestNeighbor', 
            ApplyPostTxBlur=True, PostTxVar=(2,2,2), 
            ApplyPostTxBin=True, #ThreshPostTx=0.05, 
            TxtToAddToSegLabel='', LogToConsole=False, 
            DictOfInputs={}, ListOfInputs=[], ListOfTimings=[]):
    """
    Note 09/02/2021:
        This function will, depending on the inputs, copy either:
            1. A specific segmentation
            2. All segmentations in a specific segment
            3. All segmentations in all segments
        
     
    Inputs:
    ******
    
    SrcSegFpath : string
        Filepath of the Source SEG file.
    
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a  
        single segmentation). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
        
    SrcSegLabel : string
        All or part of the Source SegmentLabel for the segment containing the
        segmentation(s) to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
        
    TrgRtsFpath : string or None (optional; None by default)
        Filepath to the Target SEG file that the segmentation(s) is/are to be 
        copied to. TrgSegFpath != None only if an existing Target SEG exists 
        and is to added to. An existing SEG will be added to only for Direct 
        copy operations.     
    
    TrgSegLabel : string or None (optional; None by default) 
        All or part of the SegmentLabel of the destination segment.
        
    TrgSliceNum : integer (optional; None by default)
        The slice index within the Target DICOM stack where the segmentation 
        will be copied to (applies only for the case of direct copies of single 
        segmentations). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the segmentation(s) will be copied to 
        will not depend on user input.
        
    ResInterp : string
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    PreResVar : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to perform binary thresholding following 
    #    resampling if a non-label (e.g. linear) interpolator is used. Example 
    #    use is on a labelmap image if a Gaussian blurring + linearly resampling 
    #    approach is taken to combat aliasing effects.
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
        
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    TxMaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Tx. If != 'default' it must be a string 
        representation of an integer.
        
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
    
    PostTxVar : tuple of floats (optional; (2,2,2) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
        
    #ThreshPostTx : float (optional; 0.05 by default)
    #    The threshold level used to perform binary thresholding after 
    #    registration transformation.
    
    TxtToAddToSegLabel : string (optional, '' by default)
        String of text to add to the segment label of the new Target SEG.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
    
    DictOfInputs : dictionary (optional; empty by default)
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings (optional; empty by default)
        A list containing the inputs that were called to CopyRoi().
    
    ListOfTimings : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopySeg().
          
                  
    Outputs:
    *******
        
    TrgSeg : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target SEG object.
      
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
    using resampling techniques the Target segmentations will appear displaced  
    w.r.t. the anatomical features highlighted in the Source image.
    """
    
    import importlib
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from ImageTools import ImportImage#, GetImageInfo
    from SegTools import GetSegDataOfInterest, ProportionOfSegsInExtent
    from SegTools import CreateSeg
    from GeneralTools import UniqueItems#, PrintIndsByRoi#, PrintTitle
    #from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    #from GeneralTools import ShiftFrame
    
    """ Start timing. """
    times = []
    times.append(time.time())
    #""" Log timing messages. """
    #TimingMsgs = []
    
    
    SrcSeg = dcmread(SrcSegFpath)
    if TrgSegFpath:
        TrgSeg = dcmread(TrgSegFpath)
    else:
        TrgSeg = None
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print(f'Running CopySeg():')
        
        SrcSegLabels = GetRoiLabels(SrcSeg)
        
        print(f'\nSrcSegLabels = {SrcSegLabels}')
        
        if TrgSegFpath:
            TrgSegLabels = GetRoiLabels(TrgSeg)
        
            print(f'\nTrgSegLabels = {TrgSegLabels}')
    
    
    """ Check the validity of the combination of inputs. """
    CheckValidityOfInputs(SrcSeg, SrcSegLabel, SrcSliceNum, TrgSeg, TrgSliceNum,
                          LogToConsole)
    
    #""" Determine which Use Case to apply. """
    #UseCaseThatApplies,\
    #UseCaseToApply = WhichUseCase(SrcSliceNum, FromSegLabel, TrgSliceNum, 
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
        
    """ Get the data of interest from the Source SEG. """
    SrcPixArrBySeg,\
    SrcF2SindsBySeg = GetSegDataOfInterest(SrcSegFpath, SrcSliceNum, 
                                           SrcSegLabel, SrcDcmDir, 
                                           LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source SEG file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    """ Raise exception if there is no data of interest from the Source SEG."""
    if not UniqueItems(SrcF2SindsBySeg):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        from SegTools import GetPFFGStoSliceIndsBySeg
        
        SrcSeg = dcmread(SrcSegFpath)
        
        Labels = GetRoiLabels(SrcSeg)
        L = len(Labels)
        
        SOPuids = GetDicomSOPuids(SrcDcmDir)
        
        F2SindsBySeg = GetPFFGStoSliceIndsBySeg(SrcSeg, SOPuids)
        
        msg = f"There are no segmentations on slice {SrcSliceNum} in any "\
              + f"segment in the Source SEG. There are {L} segments in the "\
              + f"Source SEG:"
        
        for i in range(L):
            msg += f"\nSegment {i+1} with label '{Labels[i]}' has "\
                   + f"{len(F2SindsBySeg[i])} segmentations that correspond "\
                   + f"to slices {F2SindsBySeg[i]}" 
                   
        raise Exception(msg)
    
    """
    Determine whether the segmentations that make up the segment(s) intersect 
    with the Target image extent. 
    """
    FracProp = ProportionOfSegsInExtent(PixArrBySeg=SrcPixArrBySeg,
                                        F2SindsBySeg=SrcF2SindsBySeg, 
                                        SrcDicomDir=SrcDcmDir,
                                        TrgDicomDir=TrgDcmDir,
                                        LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to check whether the Source segment(s) intersect '\
          + 'the Target grid space.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the indices in the SEG lie',
              'within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the indices in the SEG lie within the physical extent'\
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
        """ Direct copy of a segmentation without resampling or registration
        transformation. """
        
        from GeneralTools import ReplaceIndInC2SindsByRoi
        from GeneralTools import ShiftFramesInPixArrBySeg
        
        """ The segmentation on SrcSliceNum is to be copied to TrgSliceNum. 
        Modify the frame-to-slice index (from SrcSliceNum) to TrgSliceNum. """
        SrcF2SindsBySeg = ReplaceIndInC2SindsByRoi(C2SindsByRoi=SrcF2SindsBySeg,
                                                   IndToReplace=SrcSliceNum, 
                                                   ReplacementInd=TrgSliceNum)
        
        if LogToConsole:
                print(f'\nModfied SrcF2SindsBySeg = {SrcF2SindsBySeg}.')
                
        if TrgSegFpath:
            """ Direct copy of a segmentation to an existing SEG. """
            
            """ An existing Target SEG is to be modified.  Since this is a  
            Direct copy, any existing segmentations in Target with segment 
            label matching FromSegLabel are to be preserved. """
            
            """ Get the data of interest from the Target RTS. """
            TrgPixArrBySeg,\
            TrgF2SindsBySeg = GetSegDataOfInterest(TrgSegFpath, SrcSliceNum,
                                                   #SrcSegLabel, TrgDcmDir,
                                                   TrgSegLabel, TrgDcmDir,
                                                   LogToConsole)
            """ 05/03/21: Previously SrcSegLabel was used as an input in 
            GetSegDataOfInterest for Target, hence requiring that the Target
            segment have the same label as the Source segment. Now TrgSegLabel 
            has been added to the list of inputs allowing for the Target segment
            to have a different label. """
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to parse the Target SEG file for '\
                  + 'existing segmentations within the ROI.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
            
            if LogToConsole:
                print(f'TrgF2SindsByRoi = {TrgF2SindsBySeg}.')
        

    #if UseCase in ['1', '2a']:
        #""" Direct copy without resampling or registration transformation.
        """
        The segmentation will be copied to a different slice location, so the
        z-components of the indices will need to be modified. """
        
        """ Shift the in-plane (x & y) and out-of-plane (z) components of the 
        indices in SrcPixArrBySeg to account for the shift from SrcSliceNum to
        TrgSliceNum. """
        SrcPixArrBySeg,\
        SrcF2SindsBySeg = ShiftFramesInPixArrBySeg(PixArrBySeg=SrcPixArrBySeg, 
                                                   F2SindsBySeg=SrcF2SindsBySeg, 
                                                   SrcImage=SrcIm, 
                                                   SrcSliceNum=SrcSliceNum,
                                                   TrgImage=TrgIm, 
                                                   TrgSliceNum=TrgSliceNum, 
                                                   RefImage=TrgIm,
                                                   ShiftInX=True,
                                                   ShiftInY=True,
                                                   ShiftInZ=True,
                                                   Fractional=False,
                                                   LogToConsole=LogToConsole)
        
        if LogToConsole:
                print(f'\nShifted SrcF2SindsBySeg = {SrcF2SindsBySeg}.')
                
        
        if TrgSegFpath:
            """ An existing Target SEG is to be modified. 
            Concatenate TrgF2SindsBySeg and SrcF2SindsBySeg, and 
            TrgPixArrBySeg and SrcPixArrBySeg. """
            
            TrgF2SindsBySeg = [TrgF2SindsBySeg[s] + SrcF2SindsBySeg[s] for s in range(len(TrgF2SindsBySeg))]
            
            TrgPixArrBySeg = [TrgPixArrBySeg[s] + SrcPixArrBySeg[s] for s in range(len(TrgPixArrBySeg))]
        
        else:
            """ A Target SEG was not provided. """
            
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            
            if LogToConsole:
                print(f'\nSrcF2SindsBySeg = {SrcF2SindsBySeg}')
                print(f'TrgF2SindsBySeg = {TrgF2SindsBySeg}')
    
    
    
    if UseCaseToApply == '2b':
        """ Relationship-preserving copy without resampling or registration
        transformation. """
        
        if False:
            """ TrgPixArrBySeg is simply SrcPixArrBySeg, and 
            TrgF2SindsBySeg is SrcF2SindsBySeg. """
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
        
        #from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1
        from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2
        
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        TrgF2SindsBySeg = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcF2SindsBySeg, 
                                                                SrcIm, TrgIm)
        
        
        """ TrgPixArrBySeg is simply SrcPixArrBySeg. """
        TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
        
        if LogToConsole:
            print(f'\nSrcF2SindsBySeg = {SrcF2SindsBySeg}')
            print(f'TrgF2SindsBySeg (= shifted SrcF2SindsBySeg) =',
                  f'{SrcF2SindsBySeg}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        #from ImageTools import ResampleLabmapImByRoi#, GetImageInfo
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PixArrByRoi2LabmapImByRoi
        #from GeneralTools import NumOfListsAtDepthTwo
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabmapImBySeg = PixArrByRoi2LabmapImByRoi(PixArrByRoi=SrcPixArrBySeg, 
                                                     F2SindsByRoi=SrcF2SindsBySeg, 
                                                     RefIm=SrcIm,
                                                     LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to generate labelmap images from the Source '\
              'pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'\n*{msg}')
    
    
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        from ImageTools import ResampleLabmapImByRoi
        
        
        """ Resample the source labelmaps. """
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabmapImBySeg, ResSrcPixArrBySeg,\
        ResSrcF2SindsBySeg = ResampleLabmapImByRoi(LabmapImByRoi=SrcLabmapImBySeg, 
                                                   #F2SindsByRoi=SrcC2SindsByRoi,
                                                   F2SindsByRoi=SrcF2SindsBySeg,
                                                   SrcImage=SrcIm, 
                                                   TrgImage=TrgIm,
                                                   Interpolation=ResInterp,
                                                   PreResVar=PreResVar,
                                                   #PostResThresh=PostResThresh,
                                                   LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the interpolation that was used (the actual interpolation used 
        might have been different from the interpolation set) and the threshold
        used to re-binarise the resampled labelmap image (if applicable): """
        
        #print(f'\n\n\nMetadata keys in ResSrcLabmapImBySeg[0]:',
        #      ResSrcLabmapImBySeg[0].GetMetaDataKeys())
        
        for i in range(len(ResSrcLabmapImBySeg)):
            Interp = ResSrcLabmapImBySeg[i].GetMetaData("ResInterpUsed")
            
            ListOfInputs.append(f'ResInterpUsedForSeg{i} = {Interp}')
            DictOfInputs[f'ResInterpUsedForSeg{i}'] = Interp
            
            if Interp == 'BlurThenLinear':
                Thresh = ResSrcLabmapImBySeg[i].GetMetaData("PostResThreshUsed")
                
                ListOfInputs.append(f'PostResThreshUsedForSeg{i} = {Thresh}')
                DictOfInputs[f'PostResThreshUsedForSeg{i}'] = Thresh
                
        
    
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        from ImageTools import RegisterImages
        from ImageTools import TransformLabmapImByRoi
        
        if not Tx in ['rigid', 'affine', 'bspline']:
            msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
                  + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
            
            raise Exception(msg)
        
        if LogToConsole == False:
            print(f'Performing image registration using {Tx} transform...\n')
        
        times.append(time.time())
        
        """ Register SrcIm to TrgIm. """
        RegIm, RegImFilt = RegisterImages(FixIm=TrgIm, MovIm=SrcIm, Tx=Tx,
                                          MaxNumOfIters=TxMaxIters,
                                          LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to register the Source image to the Target '\
              + f'image using {Tx} transform.\n'
        ListOfTimings.append(msg)
        if LogToConsole == False:
            print(f'*{msg}')
        
        
        """ Transform SrcLabmapImBySeg using the registration transformation. 
        
        Note:
            Even though the data is being transformed the prefix 'Res' will be
            used since further operations below will need to be applied to
            either resampled or transformed data. """
        
        ResSrcLabmapImBySeg, ResSrcPixArrBySeg,\
        ResSrcF2SindsBySeg = TransformLabmapImByRoi(LabmapImByRoi=SrcLabmapImBySeg,
                                                    F2SindsByRoi=SrcF2SindsBySeg,
                                                    RegImFilt=RegImFilt,
                                                    Interpolation=TxInterp,
                                                    ApplyPostTxBlur=ApplyPostTxBlur,
                                                    PostTxVar=(2,2,2),
                                                    ApplyPostTxBin=ApplyPostTxBin,
                                                    #ThreshPostTx=ThreshPostTx, 
                                                    LogToConsole=LogToConsole)
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
        """ Get the threshold used to re-binarise the transformed labelmap 
        image: """
        for i in range(len(ResSrcLabmapImBySeg)):
            Thresh = ResSrcLabmapImBySeg[i].GetMetaData("PostTxThreshUsed")
            
            ListOfInputs.append(f'PostTxThreshUsedForSeg{i} = {Thresh}')
            DictOfInputs[f'PostTxThreshUsedForSeg{i}'] = Thresh
              
                
    
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
        
        R = len(ResSrcPixArrBySeg) # number of segments
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
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
                
                PlotPixArrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                f'Prior to {ReduceUsing} operation')
                
            if ReduceUsing == 'mean':
                from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5 # <-- is this too high?
    
                ResSrcPixArrBySeg, ResF2SindsBySeg\
                = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                         F2SindsBySeg=ResSrcF2SindsBySeg,
                                         MakeBinary=True, 
                                         BinaryThresh=BinaryThresh,
                                         LogToConsole=LogToConsole)
            
            else:
                from GeneralTools import OrFrameOfPixArrBySeg
                
                ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
                = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                       F2SindsBySeg=ResSrcF2SindsBySeg,
                                       LogToConsole=LogToConsole)
        
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment, and',
                      f'the new ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}.')
                
                PlotPixArrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                'After {ReduceUsing} operation')
                
                
        
        """ Shift the out-of-plane elements in ResSrcPixArrBySeg to account for
        the shift from SrcSliceNumber to TrgSliceNum. 
        Note:
            The function below was incorrectly applying in-plane shifts as well
            up until the added inputs ShiftInX / Y / Z was added on 10/02. 
        
        14/02: I think SrcSliceNum shouldn't be SrcSliceNum but instead
        ResSrcF2SindsBySeg[0][0], since this is the frame-to-slice number of the
        resampled/transformed labelmap.  The slice number of the original
        (Source) labelmap is not relevant.
        
        Likewise, the inputs SrcImage and TrgImage should both be TrgImage 
        since the resampled labelmap shares the TrgImage grid, and it will be
        shifted to another location, still within the TrgImage grid.
        """
        
        """ SrcSliceNum --> ResSrcSliceNum following resampling or 
        transformation. """
        ResSrcSliceNum = ResSrcF2SindsBySeg[0][0]
        
        if LogToConsole:
            #print(f'\nSrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
            #      f'{ResSrcSliceNum}')
            #print(f'TrgSliceNum = {TrgSliceNum}')
            print(f'SrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
                  f'{ResSrcSliceNum} following resampling/transformation -->',
                  f'TrgSliceNum = {TrgSliceNum} for direct copy')
            print('\nResSrcF2SindsBySeg prior to shifting of out-of-plane',
                  f'elements: {ResSrcF2SindsBySeg}\n')
            
        ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
        = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg, 
                                   F2SindsBySeg=ResSrcF2SindsBySeg, 
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
            print(f'After running ShiftFramesInPixArrBySeg():')
            print(f'ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}')
            print(f'len(ResSrcPixArrBySeg) = {len(ResSrcPixArrBySeg)}')
            for r in range(len(ResSrcPixArrBySeg)):
                print(f'type(ResSrcPixArrBySeg[{r}]) = {type(ResSrcPixArrBySeg[r])}')
                print(f'ResSrcPixArrBySeg[{r}].shape = {ResSrcPixArrBySeg[r].shape}')
            print('')
                
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copying with resampling or
        registration-transformation. """
    
        if UseCaseToApply in ['3a', '4a', '5a']:
            """ Direct copying with resampling/transformation. """
            
            if TrgSegFpath:
                """ An existing Target SEG is to be modified. 
                Concatenate TrgF2SindsBySeg and ResSrcF2SindsBySeg, and
                TrgPixArrBySeg and ResSrcPixArrBySeg. """
                
                TrgF2SindsBySeg = [TrgF2SindsBySeg[s] + ResSrcF2SindsBySeg[s] for s in range(len(TrgF2SindsBySeg))]
                
                TrgPixArrBySeg = [TrgPixArrBySeg[s] + ResSrcPixArrBySeg[s] for s in range(len(TrgPixArrBySeg))]
                
            else:
                """ A Target SEG was not provided. """
                
                TrgF2SindsBySeg = deepcopy(ResSrcF2SindsBySeg)
                
                TrgPixArrBySeg = deepcopy(ResSrcPixArrBySeg)
        
        
        else:
            """ UseCaseToApply is '3b', '4b' or '5b'. """
            
            TrgPixArrBySeg = deepcopy(ResSrcPixArrBySeg)
            
            TrgF2SindsBySeg = deepcopy(ResSrcF2SindsBySeg)
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    #print(f'\n\n\n\n\nTrgPixArrBySeg = {TrgPixArrBySeg}')
    
    
    """ Create/overwrite the Target SEG. """
    
    print('*Creating new Target SEG object...\n')
    
    #print(f'\nTrgDcmDir = {TrgDcmDir}\n')
    
    TrgSeg = CreateSeg(SrcSegFpath, TrgSegFpath, TrgSegLabel, TrgPixArrBySeg, 
                       TrgF2SindsBySeg, TrgDcmDir, SrcSegLabel,
                       TxtToAddToSegLabel, LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to create the SEG object.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print('\nEnd of CopySeg().')
        print('-'*120)
        
    return TrgSeg, DictOfInputs, ListOfInputs, ListOfTimings
    #return TrgSeg, TrgPixArrBySeg, TrgF2SindsBySeg, ListOfTimings







"""
******************************************************************************
******************************************************************************
CREATE A LIST OF INPUTS AND A DICTIONARY OF INPUTS FOR COPYROI()
******************************************************************************
******************************************************************************
"""

def CreateListOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, SrcStudyLabel,
                       SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName,
                       SrcSliceNum, TrgStudyLabel, TrgSeriesLabel, TrgAsrMod, 
                       TrgAsrName, TrgRoiName, TrgSliceNum, TxtToAddToTrgAsrName, 
                       LogToConsole, ExportLogFiles, ResInterp, PreResVar, 
                       ForceReg, Tx, TxMaxIters, TxInterp, ApplyPostTxBlur, 
                       PostTxVar, ApplyPostTxBin):
    
    ListOfInputs = [f"RunDateTime: {RunDateTime}\n",
                    f"XnatUrl = {XnatUrl}\n",
                    f"ProjId = {ProjId}\n",
                    f"SubjLabel = {SubjLabel}\n",
                    f"SrcStudyLabel = {SrcStudyLabel}\n",
                    f"SrcSeriesLabel = {SrcSeriesLabel}\n",
                    f"SrcAsrMod = {SrcAsrMod}\n",
                    f"SrcAsrName = {SrcAsrName}\n",
                    f"SrcRoiName = {SrcRoiName}\n",
                    f"SrcSliceNum = {SrcSliceNum}\n",
                    f"TrgStudyLabel = {TrgStudyLabel}\n",
                    f"TrgSeriesLabel = {TrgSeriesLabel}\n",
                    f"TrgAsrMod = {TrgAsrMod}\n",
                    f"TrgAsrName = {TrgAsrName}\n",
                    f"TrgRoiName = {TrgRoiName}\n",
                    f"TrgSliceNum = {TrgSliceNum}\n",
                    f"TxtToAddToTrgAsrName = {TxtToAddToTrgAsrName}\n",
                    f"LogToConsole = {LogToConsole}",
                    f"ExportLogFiles = {ExportLogFiles}",
                    f"ResInterpSet = {ResInterp}\n",
                    f"PreResVar = {PreResVar}\n",
                    f"ForceReg = {ForceReg}\n",
                    f"Tx = {Tx}\n",
                    f"TxMaxIters = {TxMaxIters}\n",
                    f"TxInterp = {TxInterp}\n",
                    f"ApplyPostTxBlur = {ApplyPostTxBlur}\n",
                    f"PostTxVar = {PostTxVar}\n",
                    f"ApplyPostTxBin = {ApplyPostTxBin}\n"
                    ]
    
    return ListOfInputs





def CreateDictOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, SrcStudyLabel,
                       SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName,
                       SrcSliceNum, TrgStudyLabel, TrgSeriesLabel, TrgAsrMod, 
                       TrgAsrName, TrgRoiName, TrgSliceNum, TxtToAddToTrgAsrName, 
                       LogToConsole, ExportLogFiles, ResInterp, PreResVar, 
                       ForceReg, Tx, TxMaxIters, TxInterp, ApplyPostTxBlur, 
                       PostTxVar, ApplyPostTxBin):
    
    DictOfInputs = {'RunDateTime' : RunDateTime,
                    'XnatUrl' : XnatUrl,
                    'ProjId' : ProjId,
                    'SubjLabel' : SubjLabel,
                    'SrcStudyLabel' : SrcStudyLabel,
                    'SrcSeriesLabel' : SrcSeriesLabel,
                    'SrcAsrMod' : SrcAsrMod,
                    'SrcAsrName' : SrcAsrName,
                    'SrcRoiName' : SrcRoiName,
                    'SrcSliceNum' : SrcSliceNum,
                    'TrgStudyLabel' : TrgStudyLabel,
                    'TrgSeriesLabel' : TrgSeriesLabel,
                    'TrgAsrMod' : TrgAsrMod,
                    'TrgAsrName' : TrgAsrName,
                    'TrgRoiName' : TrgRoiName,
                    'TrgSliceNum' : TrgSliceNum,
                    'TxtToAddToTrgAsrName' : TxtToAddToTrgAsrName,
                    'LogToConsole' : LogToConsole,
                    'ExportLogFiles' : ExportLogFiles,
                    'ResInterpSet' : ResInterp,
                    'PreResVar' : PreResVar,
                    'ForceReg' : ForceReg,
                    'Tx' : Tx,
                    'TxMaxIters' : TxMaxIters,
                    'TxInterp' : TxInterp,
                    'ApplyPostTxBlur' : ApplyPostTxBlur,
                    'PostTxVar' : PostTxVar,
                    'ApplyPostTxBin' : ApplyPostTxBin
                    }
    
    return DictOfInputs

    
    
    

"""
******************************************************************************
******************************************************************************
COPY A CONTOUR / ROI / RTSTRUCT / SEGMENTATION / SEGMENT / SEG
******************************************************************************
******************************************************************************
"""


def CopyLocalRoi(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir,
                 TrgRoiFpath=None, TrgRoiName=None, TrgSliceNum=None, 
                 ResInterp='BlurThenLinear', PreResVar=(1,1,1), ForceReg=False, 
                 Tx='affine', TxMaxIters='512', TxInterp='NearestNeighbor', 
                 ApplyPostTxBlur=True, PostTxVar=(2,2,2), ApplyPostTxBin=True,  
                 TxtToAddTrgRoiName='', LogToConsole=False, ExportLogFiles=False):
    """
    Copy a contour/ROI/RTSTRUCT/segmentation/segment/SEG from locally sourced 
    data.
        
    Inputs:
    ******
    
    SrcRoiFpath : string
        Filepath of Source RTS/SEG file.
    
    SrcSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied (counting from 0).
        
    SrcRoiName : string
        All or part of the Source ROIName/SegmentLabel of the ROI/segment 
        containing the contour/segmentation to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
        
    TrgRoiFpath : string or None (optional; None by default)
        Filepath to the Target RTS/SEG file that the contour(s)/segmentation(s) 
        is/are to be copied to. TrgRoiFpath != None only if an existing Target 
        RTS/SEG exists and is to be added to. An existing RTS/SEG will be added 
        to only for Direct copy operations.        
    
    TrgRoiName : string or None (optional; None by default)
        The ROIName or SegmentLabel of the destination ROI/segment.
    
    TrgSliceNum : integer (optional; None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to (counting from 0).  This only applies 
        for Direct copies, hence the default value None.
        
    ResInterp : string (optional; 'BlurThenLinear' by default)
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding)
    
    PreResVar : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
        
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    TxMaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Tx. If != 'default' it must be a string 
        representation of an integer.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
    
    PostTxVar : tuple of floats (optional; (2,2,2) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
    
    TxtToAddTrgRoiName : string (optional, '' by default)
        String of text to pass to CreateRts()/CreateSeg(). The string will be
        appended to the RTS StructureSetLabel or SEG SeriesDescription, and to 
        the filename of the exported (new) Target RTS/SEG.
            
    LogToConsole : boolean (default False)
        If True, some results will be logged to the console.
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
          
                  
    Outputs:
    *******
        
    TrgRoi : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing RTS/SEG) 
        Target RTS/SEG object.
        
    DictOfInputs : dictionary
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
        
    ListOfTimings : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopyRoi() and subsequent tasks called within CopyRts() or CopySeg().
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours/segmentations will appear 
    displaced w.r.t. the anatomical features highlighted in the Source image.
    """

    import time
    import os
    from pydicom import dcmread
    #from DicomTools import IsSameModalities
    from GeneralTools import ExportDictionaryToJson, ExportListToTxt
    
    RunDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    ListOfInputs = [f'RunDateTime: {RunDateTime}\n',
                    f'SrcRoiFpath = {SrcRoiFpath}\n',
                    f'SrcSliceNum = {SrcSliceNum}\n',
                    f'SrcRoiName = {SrcRoiName}\n',
                    f'SrcDcmDir = {SrcDcmDir}\n',
                    f'TrgDcmDir = {TrgDcmDir}\n',
                    f'TrgRoiFpath = {TrgRoiFpath}\n',
                    f'TrgRoiName = {TrgRoiName}\n',
                    f'TrgSliceNum = {TrgSliceNum}\n',
                    f'ResInterpSet = {ResInterp}\n',
                    f'PreResVar = {PreResVar}\n',
                    #f'PostResThresh = {PostResThresh}\n',
                    f'ForceReg = {ForceReg}\n',
                    f'Tx = {Tx}\n',
                    f'TxMaxIters = {TxMaxIters}\n',
                    f'TxInterp = {TxInterp}\n',
                    f'ApplyPostTxBlur = {ApplyPostTxBlur}\n',
                    f'PostTxVar = {PostTxVar}\n',
                    f'ApplyPostTxBin = {ApplyPostTxBin}\n',
                    #f'ThreshPostTx = {ThreshPostTx}\n',
                    f'TxtToAddTrgRoiName = {TxtToAddTrgRoiName}\n',
                    f'LogToConsole = {LogToConsole}']
    
    DictOfInputs = {'RunDateTime' : RunDateTime,
                    'SrcRoiFpath' : SrcRoiFpath,
                    'SrcSliceNum' : SrcSliceNum,
                    'SrcRoiName' : SrcRoiName,
                    'SrcDcmDir' : SrcDcmDir,
                    'TrgDcmDir' : TrgDcmDir,
                    'TrgRoiFpath' : TrgRoiFpath,
                    'TrgRoiName' : TrgRoiName,
                    'TrgSliceNum' : TrgSliceNum,
                    'ResInterpSet' : ResInterp,
                    'PreResVar' : PreResVar,
                    #'PostResThresh' : PostResThresh,
                    'ForceReg' : ForceReg,
                    'Tx' : Tx,
                    'TxMaxIters' : TxMaxIters,
                    'TxInterp' : TxInterp,
                    'ApplyPostTxBlur' : ApplyPostTxBlur,
                    'PostTxVar' : PostTxVar,
                    'ApplyPostTxBin' : ApplyPostTxBin,
                    #'ThreshPostTx' : ThreshPostTx,
                    'TxtToAddTrgRoiName' : TxtToAddTrgRoiName,
                    'LogToConsole' : LogToConsole}
    
    
    ListOfTimings = [f'RunDateTime: {RunDateTime}\n']
    
    times = []
    
    
    """ Start timing. """
    times.append(time.time())
    
    """ Establish whether the inputs are valid. """
    CheckValidityOfInputs(SrcRoiFpath, SrcRoiName, SrcSliceNum, TrgSliceNum, 
                          TrgRoiFpath, LogToConsole)
    
    times.append(time.time())
    Dtime = round(1000*(times[-1] - times[-2]), 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} ms to check the validity of inputs to CopyRoi().\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
    """ Determine which Use Case to apply. """
    UseCaseThatApplies,\
    UseCaseToApply = WhichUseCase(SrcSliceNum, SrcRoiName, TrgSliceNum, 
                                  SrcDcmDir, TrgDcmDir, ForceReg,
                                  LogToConsole=True)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        #print(f'\n\nDone.  Took {Dtime} s to run.')
        msg = f'Took {Dtime} s to determine which UseCase applies.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
    
    DictOfInputs['UseCaseThatApplies'] = UseCaseThatApplies
    DictOfInputs['UseCaseToApply'] = UseCaseToApply
    
    #print(f'\n\n\nSrcRoiFpath = {SrcRoiFpath}')
    
    Modality = dcmread(SrcRoiFpath).Modality
    
    if Modality == 'RTSTRUCT':        
        TrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyRts(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, 
                  TrgDcmDir, UseCaseToApply, TrgRoiFpath, TrgRoiName, 
                  TrgSliceNum, ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, 
                  TxInterp, ApplyPostTxBlur, PostTxVar, ApplyPostTxBin,    
                  TxtToAddTrgRoiName, LogToConsole,
                  DictOfInputs, ListOfInputs, ListOfTimings)
    
    elif Modality == 'SEG':
        TrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopySeg(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir,
                  TrgDcmDir, UseCaseToApply, TrgRoiFpath, TrgRoiName, 
                  TrgSliceNum, ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, 
                  TxInterp, ApplyPostTxBlur, PostTxVar, ApplyPostTxBin,    
                  TxtToAddTrgRoiName, LogToConsole, 
                  DictOfInputs, ListOfInputs, ListOfTimings)
        
    else:
        msg = f'The Source modality ({Modality}) must be "RTS" or "SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to copy the ROI(s) and create the {Modality}'\
              + ' object.\n'  
        ListOfTimings.append(msg)
        print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export DictOfInputs as a json, and ListOfInputs as a txt, to a sub-
        directory "logs" in the current working directory. """
        JsonFname = RunDateTime + '_DictOfInputs.json'
        TxtFname = RunDateTime + '_ListOfInputs.txt'
        
        CWD = os.getcwd()
        
        ExportDir = os.path.join(CWD, 'logs')
        
        JsonFpath = os.path.join(ExportDir, JsonFname)
        TxtFpath = os.path.join(ExportDir, TxtFname)
        
        ExportDictionaryToJson(DictOfInputs, JsonFname, ExportDir)
        ExportListToTxt(ListOfInputs, TxtFname, ExportDir)
            
        print('Dictionary of inputs to CopyRoi() saved to:\n', JsonFpath, 
              '\nand list of inputs saved to:\n', TxtFpath, '\n')
        
    return TrgRoi, DictOfInputs, ListOfInputs, ListOfTimings






def CopyXnatRoi(XnatUrl, XnatSession, ProjId, SubjLabel, SrcStudyLabel, 
                SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName, SrcSliceNum,  
                TrgStudyLabel, TrgSeriesLabel, TrgAsrMod=None, TrgAsrName=None, 
                TrgRoiName=None, TrgSliceNum=None, TxtToAddToTrgAsrName='', 
                PathsDict=None, XnatDownloadDir='default', 
                LogToConsole=False, ExportLogFiles=False,
                ResInterp='BlurThenLinear', PreResVar=(1,1,1), ForceReg=False,   
                Tx='affine', TxMaxIters='512', TxInterp='NearestNeighbor', 
                ApplyPostTxBlur=True, PostTxVar=(2,2,2), ApplyPostTxBin=True):
    
    """
    Copy a contour/ROI/RTSTRUCT/segmentation/segment/SEG from data downloaded
    from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
        
    XnatSession : requests session
        
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    SrcStudyLabel : string
        The Source study / experiment label. 
    
    SrcSeriesLabel : string
        The Source series label / scan ID.
    
    SrcAsrMod : string
        The Source assessor modality.
        
    SrcAsrName : string
        The Source assessor name (StructureSetLabel or SeriesDescription).
        
    SrcRoiName : string
        The Source ROIName or SegmentLabel.
    
    SrcSliceNum : integer (0-indexed)
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied.
    
    TrgStudyLabel : string
        The Target study / experiment label. 
    
    TrgSeriesLabel : string
        The Target series label / scan ID.
    
    TrgAsrMod : string (optional but required if TrgAsrName != None; 
    None by default)
        The Target assessor modality.
        
    TrgAsrName : string (optional but required if TrgRoiMod != None; 
    None by default)
        The Target assessor name (StructureSetLabel or SeriesDescription). If 
        provided and if a direct copy is to be made, existing contours/
        segmentations for the ROI/segment will be preserved.
    
    TrgRoiName : string (optional; None by default)
        The Target ROIName or SegmentLabel.
        
    TrgSliceNum : integer (optional unless making a direct copy; 0-indexed; 
    None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to.  This only applies for direct copies, 
        hence the default value None.
    
    TxtToAddToTrgAsrName : string (optional, '' by default)
        If provided the string of text will be appended to the assessor name 
        (StructureSetLabel or SeriesDescription), and to the filename of the 
        exported (new) Target RTS/SEG.
    
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded. If provided, paths of
        newly downloaded data will be added to PathsDict. If not provided, a
        new dictionary will be initialised.
    
    XnatDownloadDir : string (optional; 'default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjectLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    
    LogToConsole : boolean (optional; False by default)
        If True, some results will be logged to the console.
    
    ExportLogFiles : boolean (optional; False by default)
        If True, log files will be exported.
        
    ResInterp : string (optional; 'BlurThenLinear' by default)
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding)
    
    PreResVar : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
        
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    TxMaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Tx. If != 'default' it must be a string 
        representation of an integer.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
    
    PostTxVar : tuple of floats (optional; (2,2,2) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
    
                  
    Outputs:
    *******
        
    TrgRoi : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing RTS/SEG) 
        Target RTS/SEG object.
    
    PathsDict : dictionary
        Dictionary containing paths of data downloaded.
        
    XnatSession : requests session
        
    DictOfInputs : dictionary
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
        
    TimingMsgs : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopyRoi() and subsequent tasks called within CopyRts() or CopySeg().
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours/segmentations will appear 
    displaced w.r.t. the anatomical features highlighted in the Source image.
    """

    import time
    import os
    #from pypref import Preferences
    #from pydicom import dcmread
    import importlib
    import XnatTools
    importlib.reload(XnatTools)
    from XnatTools import DownloadScan, DownloadAssessor
    from GeneralTools import ExportDictionaryToJson, ExportListToTxt
    
    RunDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #Defaults = Preferences(directory=os.getcwd(), filename='CopyRoiDefaults.py')
    #
    #ResInterp = Defaults.get('ResInterp')
    #PreResVar = Defaults.get('PreResVar')
    #ForceReg = Defaults.get('ForceReg')
    #Tx = Defaults.get('Tx')
    #TxMaxIters = Defaults.get('TxMaxIters')
    #TxInterp = Defaults.get('TxInterp'),
    #ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
    #PostTxVar = Defaults.get('PostTxVar')
    #ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
    
    ListOfInputs\
    = CreateListOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, SrcStudyLabel, 
                         SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName, 
                         SrcSliceNum, TrgStudyLabel, TrgSeriesLabel, TrgAsrMod, 
                         TrgAsrName, TrgRoiName, TrgSliceNum, 
                         TxtToAddToTrgAsrName, LogToConsole, ExportLogFiles, 
                         ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, 
                         TxInterp, ApplyPostTxBlur, PostTxVar, ApplyPostTxBin)
    
    DictOfInputs\
    = CreateDictOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, SrcStudyLabel, 
                         SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName, 
                         SrcSliceNum, TrgStudyLabel, TrgSeriesLabel, TrgAsrMod, 
                         TrgAsrName, TrgRoiName, TrgSliceNum, 
                         TxtToAddToTrgAsrName, LogToConsole, ExportLogFiles, 
                         ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, 
                         TxInterp, ApplyPostTxBlur, PostTxVar, ApplyPostTxBin)
    
    
    TimingMsgs = [f'RunDateTime: {RunDateTime}\n']
    
    times = []
    
    
    """ Start timing. """
    times.append(time.time())
    
    """ Download data from XNAT. """
    #if PathsDict == None:
    #    PathsDict = {}
    
    #print(f'XnatSession = {XnatSession}')
    
    PathsDict,\
    XnatSession = DownloadScan(XnatUrl, ProjId, SubjLabel, 
                               SrcStudyLabel, SrcSeriesLabel, 
                               XnatSession, XnatDownloadDir, PathsDict)
    
    PathsDict,\
    XnatSession = DownloadScan(XnatUrl, ProjId, SubjLabel, 
                               TrgStudyLabel, TrgSeriesLabel, 
                               XnatSession, XnatDownloadDir, PathsDict)
    
    PathsDict,\
    XnatSession = DownloadAssessor(XnatUrl, ProjId, SubjLabel, SrcStudyLabel, 
                                   SrcSeriesLabel, SrcAsrMod, SrcAsrName, 
                                   XnatSession, XnatDownloadDir, PathsDict)
    
    if TrgRoiName != None:
        PathsDict,\
        XnatSession = DownloadAssessor(XnatUrl, ProjId, SubjLabel, 
                                       TrgStudyLabel, TrgSeriesLabel, TrgAsrMod,
                                       TrgAsrName, XnatSession, XnatDownloadDir, 
                                       PathsDict)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to download data from XNAT.\n'
        TimingMsgs.append(msg)
        print(f'*{msg}')
    
    
    """ Get the DICOM directories and assessor file paths: """
    SrcDcmDir = PathsDict[ProjId][SubjLabel][SrcStudyLabel][SrcSeriesLabel]\
                ['DicomDir']
    
    TrgDcmDir = PathsDict[ProjId][SubjLabel][TrgStudyLabel][TrgSeriesLabel]\
                ['DicomDir']
                
    #print(f'\nTrgDcmDir = {TrgDcmDir}\n')
    
    SrcAsrFpath = PathsDict[ProjId][SubjLabel][SrcStudyLabel][SrcSeriesLabel]\
                  ['assessors'][SrcAsrMod][SrcAsrName]['AsrFpath']
    
    if TrgAsrName == None:
        TrgAsrFpath = None
    else:
        TrgAsrFpath = PathsDict[ProjId][SubjLabel][TrgStudyLabel][TrgSeriesLabel]\
                      ['assessors'][TrgAsrMod][TrgAsrName]['AsrFpath']
    
    
    #""" Establish whether the inputs are valid. """
    #CheckValidityOfInputs(SrcRoiName, SrcSliceNum, TrgRoiName, TrgSliceNum, 
    #                      LogToConsole)
    
    #print(f'\nTrgSliceNum = {TrgSliceNum}')
    
    """ Determine which Use Case to apply. """
    UseCaseThatApplies,\
    UseCaseToApply = WhichUseCase(SrcSliceNum, TrgSliceNum, SrcDcmDir, 
                                  TrgDcmDir, ForceReg, LogToConsole=True)
    
    #print(f'\nTrgSliceNum = {TrgSliceNum}')
    #print(f'\nTrgDcmDir = {TrgDcmDir}\n')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        #print(f'\n\nDone.  Took {Dtime} s to run.')
        msg = f'Took {Dtime} s to determine which UseCase applies.\n'
        TimingMsgs.append(msg)
        print(f'*{msg}')
        
    
    DictOfInputs['UseCaseThatApplies'] = UseCaseThatApplies
    DictOfInputs['UseCaseToApply'] = UseCaseToApply
    
    #print(f'\n\n\nSrcRoiFpath = {SrcRoiFpath}')
    
    #Modality = dcmread(SrcRoiFpath).Modality
    
    if SrcAsrMod == 'RTSTRUCT':        
        TrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyRts(SrcAsrFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir, 
                  UseCaseToApply, TrgAsrFpath, TrgRoiName, TrgSliceNum, 
                  ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, TxInterp, 
                  ApplyPostTxBlur, PostTxVar, ApplyPostTxBin,  
                  TxtToAddToTrgAsrName, LogToConsole,
                  DictOfInputs, ListOfInputs, TimingMsgs)
    
    elif SrcAsrMod == 'SEG':
        #print(f'\nTrgSliceNum = {TrgSliceNum}')
        #print(f'\nTrgDcmDir = {TrgDcmDir}\n')
        
        TrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopySeg(SrcAsrFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir, 
                  UseCaseToApply, TrgAsrFpath, TrgRoiName, TrgSliceNum, 
                  ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, TxInterp,
                  ApplyPostTxBlur, PostTxVar, ApplyPostTxBin,    
                  TxtToAddToTrgAsrName, LogToConsole, 
                  DictOfInputs, ListOfInputs, TimingMsgs)
        
    else:
        msg = f'The modality of the Source assessor ({SrcAsrMod}) must be '\
              + '"RTSTRUCT" or "SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to copy the ROI(s) and create the {SrcAsrMod}'\
              + ' object.\n'  
        TimingMsgs.append(msg)
        print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export DictOfInputs as a json, and ListOfInputs as a txt, to a sub-
        directory "logs" in the current working directory. """
        JsonFname = RunDateTime + '_DictOfInputs.json'
        TxtFname = RunDateTime + '_ListOfInputs.txt'
        
        CWD = os.getcwd()
        
        ExportDir = os.path.join(CWD, 'logs')
        
        JsonFpath = os.path.join(ExportDir, JsonFname)
        TxtFpath = os.path.join(ExportDir, TxtFname)
        
        ExportDictionaryToJson(DictOfInputs, JsonFname, ExportDir)
        ExportListToTxt(ListOfInputs, TxtFname, ExportDir)
            
        print('Dictionary of inputs to CopyRoi() saved to:\n', JsonFpath, 
              '\nand list of inputs saved to:\n', TxtFpath, '\n')
        
    return TrgRoi, PathsDict, XnatSession, DictOfInputs, ListOfInputs,\
           TimingMsgs







"""
******************************************************************************
******************************************************************************
CHECK FOR ERRORS IN NEW TARGET RTS / SEG
******************************************************************************
******************************************************************************
"""

def ErrorCheckRoi(Roi, DicomDir, LogToConsole=False, DictOfInputs=None,
                  ExportLogFiles=False):
    """
    Check a RTS/SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the RTS/SEG.  
    
    Inputs:
    ******
    
    Roi : Pydicom object
        RTS/SEG ROI object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Roi.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
    
    DictOfInputs : dictionary or None (optional; None by default)
        If !=None, it is a dictionary containing the inputs that were called to 
        CopyRoi().
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe info or any errors that are found.
    
    Nerrors : integer
        The number of errors found.
    """
    
    import time
    import os
    from GeneralTools import ExportListToTxt
    
    Modality = Roi.Modality
    
    """ Start timing. """
    times = []
    times.append(time.time())
    
    if Modality == 'RTSTRUCT':
        #import RtsTools
        #import importlib
        #importlib.reload(RtsTools)
        from RtsTools import ErrorCheckRts
        
        LogList, Nerrors = ErrorCheckRts(Roi, DicomDir, LogToConsole)
        
    elif Modality == 'SEG':
        #import importlib
        #import SegTools
        #importlib.reload(SegTools)
        from SegTools import ErrorCheckSeg
        
        LogList, Nerrors = ErrorCheckSeg(Roi, DicomDir, LogToConsole)
        
    else:
        msg = f'The modality ({Modality}) must be either "RTS" or "SEG".'
        
        raise Exception(msg)
        
        LogList = None
        Nerrors = None
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to error check the {Modality} object.\n')
    
    
    if ExportLogFiles:
        """ Export log of error check results. """
        if not DictOfInputs == None:
            DateTime = DictOfInputs['RunDateTime']
        else:
            DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        Fname = DateTime + '_ErrorCheckLog.txt'
        
        CWD = os.getcwd()
        ExportDir = os.path.join(CWD, 'logs')
        
        Fpath = os.path.join(ExportDir, Fname)
        
        ExportListToTxt(LogList, Fname, ExportDir)
            
        print('Log of error checks saved to:\n', Fpath, '\n')
    
    return LogList, Nerrors






"""
******************************************************************************
******************************************************************************
EXPORT NEW TARGET RTS / SEG TO DISK
******************************************************************************
******************************************************************************
"""

def ExportTrgRoi(TrgRoi, SrcRoiFpath, ExportDir, Fname='', DictOfInputs=None):
    """
    Export RTS/SEG to disk.  
    
    Inputs:
    ******
    
    TrgRoi : Pydicom object
        Target RTS/SEG assessor to be exported.
        
    SrcRoiFpath : string
        Full path of the Source RTS/SEG file (used to generate the filename of
        the new RTS/SEG file).
                              
    ExportDir : string
        Directory where the new RTS/SEG is to be exported.
        
    Fname : string
        File name to assign to new assessor.
    
    NamePrefix : string (optional; '' by default)
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
    
    DictOfInputs : dictionary or None (optional; None by default)
        If not None, a dictionary containing the inputs that were called to 
        CopyRoi().
                           
    
    Outputs:
    *******
    
    TrgRoiFpath : string
        Full path of the exported Target RTS/SEG file.
    """
    
    import os
    import time
    from pathlib import Path
    
    # Get the filename of the original RTS/SEG file:
    #SrcRoiFname = os.path.split(SrcRoiFpath)[1]
    
    if not os.path.isdir(ExportDir):
        #os.mkdir(ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    if not DictOfInputs == None:
        DateTime = DictOfInputs['RunDateTime']
    else:
        DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #FnamePrefix = DateTime + '_' + NamePrefix + '_from_'
    #FnamePrefix = DateTime + '_' + TxtToAddToFname.replace(' ', '_')
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = FnamePrefix + SrcRoiFname
    #TrgRoiFname = FnamePrefix
    if Fname == '':
        TrgRoiFname = DateTime + '.dcm'
    else:
        TrgRoiFname = Fname.replace(' ', '_') + '.dcm'
    
    TrgRoiFpath = os.path.join(ExportDir, TrgRoiFname)
    
    TrgRoi.save_as(TrgRoiFpath)
        
    print(f'New Target {TrgRoi.Modality} exported to:\n', TrgRoiFpath, '\n')
    
    return TrgRoiFpath







"""
******************************************************************************
******************************************************************************
RUN COPYROI()
******************************************************************************
******************************************************************************
"""


def RunCopyLocalRoi(TestNums, LogToConsole=False, TxtToAddTrgRoiName='', 
                    ResInterp='BlurThenLinear', PreResVar=(1,1,1), 
                    ForceReg=False, Tx='affine', TxMaxIters='512',
                    TxInterp='NearestNeighbor', ApplyPostTxBlur=True, 
                    PostTxVar=(2,2,2), ApplyPostTxBin=True, ExportRoi=True, 
                    PlotResults=False, PlotAllSlices=False, 
                    ExportPlot=False, ExportLogFiles=True):
    """
    
    Inputs:
    ******
    
    TestNums : list of strings
        A list of strings denoting the tests to run (as defined in 
        CreateDictOfPaths()). Several tests can be run in series by providing
        a list of comma-separated strings (e.g. ['SR1', 'SR2', 'SR3']). To 
        run a single test a single-itemed list must be provided (e.g. ['RD4']).
    
    LogToConsole : boolean (optional; False by default)
        If True, intermediate results will be logged to the console during the
        running of CopyRoi().
    
    TxtToAddTrgRoiName : string (optional, '' by default)
        String of text to pass to CreateRts()/CreateSeg(). The string will be
        appended to the RTS StructureSetLabel or SEG SeriesDescription, and to 
        the filename of the exported (new) Target RTS/SEG.
    
    ResInterp : string
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    PreResVar : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to perform binary thresholding following 
    #    resampling using a non-label (e.g. linear) interpolator. Example use is 
    #    on a labelmap image if a Gaussian blurring + linearly resampling 
    #    approach is taken to combat aliasing effects.
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    TxMaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Tx. If != 'default' it must be a string 
        representation of an integer.
     
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
    
    PostTxVar : tuple of floats (optional; (2,2,2) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
        
    #ThreshPostTx : float (optional; 0.05 by default)
    #    The threshold level used to perform binary thresholding after 
    #    registration transformation.
    
    ExportRoi : boolean (optional; True by default)
        If True, the new Target RTS/SEG will be exported to a directory called
        "new_RTS"/"new_SEG" in the current working directory.
        
    PlotResults : boolean (optional; False by default)
        If True, the new RTS/SEG will be parsed and all contours/segmentations
        will be plotted.
    
    PlotAllSlices : boolean (optional; False by default)
        If True, and if PlotResults is True, all slices in Source and Target
        will be plotted. If False, only slices containing contours/segmentations
        will be plotted.
    
    ExportPlot : boolean (optional; False by default)
        If True, and if PlotResults is True, the plot will be exported to a .jpg
        in a directory called "plots_RTS"/"plots_SEG" in the current working 
        directory.
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
    """

    import time
    import os
    from pydicom import dcmread
    #from DicomTools import IsSameModalities
    from GeneralTools import PrintTitle, ExportListToTxt, ExportDictionaryToJson
    from TestingTools import GetPathInputs
    
    if not isinstance(TestNums, list):
        msg = 'The input argument "TestNums" must be a list of character '\
              + 'strings.'
        print(msg)
    
    CWD = os.getcwd()
    RtsExportDir = os.path.join(CWD, 'new_RTS')
    SegExportDir = os.path.join(CWD, 'new_SEG')
    RtsPlotExportDir = os.path.join(CWD, 'plots_RTS')
    SegPlotExportDir = os.path.join(CWD, 'plots_SEG')
    LogExportDir = os.path.join(CWD, 'logs')
    
    ListOfTimings = [] # list of strings (messages outputted to console)
    
    """ Start timing. """
    Times = []
    Times.append(time.time())
    
    
    for TestNum in TestNums:
        TestRunTimes = [] # this this Test run only
        TestRunTimes.append(time.time())
        
        PrintTitle(f'Running Test {TestNum}:')
    
        SrcRoiFpath, SrcSliceNum, SrcRoiName, TrgSliceNum, SrcDcmDir,\
        TrgDcmDir, TrgRoiFpath, TxtToAddTrgRoiName, SrcLabel,\
        TrgLabel = GetPathInputs(TestNum)
        
        Times.append(time.time())
        Dtime = round(1000*(Times[-1] - Times[-2]), 1)
        msg = f'Took {Dtime} ms to get the inputs to CopyRoi().\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if LogToConsole:
            print(f'Inputs to CopyRoi():\n SrcRoiName = {SrcRoiName}\n',
                      f'SrcSliceNum = {SrcSliceNum}\n',
                      f'TrgSliceNum = {TrgSliceNum}\n')
        
    
        """ Import the RTS or SEGs: """
        SrcRoi = dcmread(SrcRoiFpath)
        SrcModality = SrcRoi.Modality
        
        #print(f'\n\nSrcRoiFpath = {SrcRoiFpath}')
        #print(f'SrcModality = {SrcModality}')
    
        if TrgRoiFpath:
            TrgRoi = dcmread(TrgRoiFpath)
        else:
            TrgRoi = None
        
        """ New ROI label for RTS StructureSetLabel / SEG Series Description 
        and RTS/SEG filename: """
        if 'RR' in TestNum or 'RD' in TestNum:
            NewRoiLabel = SrcRoi.StructureSetLabel + TxtToAddTrgRoiName
        else:
            """ 'SR' in TestNum or 'SD' in TestNum: """
            NewRoiLabel = SrcRoi.SeriesDescription + TxtToAddTrgRoiName
            
        
        """ Text to add to file names: """
        TxtToAddToFname = f'Test_{TestNum}_{NewRoiLabel}'
        
        
        """ Copy Contours/Segmentations """
        
        NewTrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyLocalRoi(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, 
                       TrgDcmDir, TrgRoiFpath, TrgSliceNum, ResInterp, PreResVar, #PostResThresh, 
                       ForceReg, Tx, TxMaxIters, TxInterp, ApplyPostTxBlur, 
                       PostTxVar, ApplyPostTxBin, #ThreshPostTx, 
                       TxtToAddTrgRoiName, LogToConsole, ExportLogFiles=False)
        
        
        if ExportLogFiles:
            """ Export DictOfInputs as a json, and ListOfInputs as a txt, to a 
            sub-directory "logs" in the current working directory. """
            
            if not DictOfInputs == None:
                DateTime = DictOfInputs['RunDateTime']
            else:
                DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
                
            JsonFname = DateTime + f'_TestRun_{TestNum}_DictOfInputs.json'
            TxtFname = DateTime + f'_TestRun_{TestNum}_ListOfInputs.txt'
            
            CWD = os.getcwd()
            
            ExportDir = os.path.join(CWD, 'logs')
            
            JsonFpath = os.path.join(ExportDir, JsonFname)
            TxtFpath = os.path.join(ExportDir, TxtFname)
            
            ExportDictionaryToJson(DictOfInputs, JsonFname, ExportDir)
            ExportListToTxt(ListOfInputs, TxtFname, ExportDir)
                
            print('Dictionary of inputs to CopyRoi() saved to:\n', JsonFpath, 
                  '\nand list of inputs saved to:\n', TxtFpath, '\n')
        
        
        
        """ Error check RTS/SEG """
        Times.append(time.time())
        
        ErrorList, Nerrors = ErrorCheckRoi(NewTrgRoi, TrgDcmDir, LogToConsole,
                                           DictOfInputs, False)
        
        Times.append(time.time())
        Dtime = round(Times[-1] - Times[-2], 1)
        msg = f'Took {Dtime} s to error check the {SrcModality}.  There were '\
              + f'{Nerrors} errors found.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if ExportLogFiles:
            """ Export log of error check results. """
            
            Fname = DateTime + f'_TestRun_{TestNum}_ErrorCheckLog.txt'
            
            CWD = os.getcwd()
            ExportDir = os.path.join(CWD, 'logs')
            
            Fpath = os.path.join(ExportDir, Fname)
            
            ExportListToTxt(ErrorList, Fname, ExportDir)
                
            print('Log of error checks saved to:\n', Fpath, '\n')
        
        
    
        """ Export RTS/SEG """
    
        if ExportRoi:# and not Nerrors:
            UseCaseThatApplies = DictOfInputs['UseCaseThatApplies']
            UseCaseToApply = DictOfInputs['UseCaseToApply']
            
            if ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
                TxtToAddToFname += f'_ForcedReg_{Tx}'
                
            if not ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
                #TxtToAddToFname += f'_{ResInterp}'
                
                """ The actual interpolation used for resampling might not be
                ResInterp. """
                DiffInterpByRoi = []
                
                for key, val in DictOfInputs.items():
                    if 'ResInterpUsedFor' in key and not ResInterp in key:
                        DiffInterpByRoi.append(val)
                
                """ If DiffInterpByRoi is not empty, use the first item. """
                if DiffInterpByRoi:
                    TxtToAddToFname += f'_{DiffInterpByRoi[0]}'
                else:
                    if ResInterp == 'NearestNeighbor':
                        TxtToAddToFname += '_NN'
                    else:
                        TxtToAddToFname += f'_{ResInterp}'
                
                
                
            if UseCaseToApply in ['5a', '5b']:
                if TxInterp == 'NearestNeighbor':
                    TxtToAddToFname += '_NN'
                else:
                    TxtToAddToFname += f'_{TxInterp}'
            
            if NewTrgRoi.Modality == 'RTSTRUCT':
                NewTrgRoiFpath = ExportTrgRoi(NewTrgRoi, SrcRoiFpath, RtsExportDir, 
                                              TxtToAddToFname, DictOfInputs)
            else:
                NewTrgRoiFpath = ExportTrgRoi(NewTrgRoi, SrcRoiFpath, SegExportDir, 
                                              TxtToAddToFname, DictOfInputs)
        
        
        """ Plot Contours """
    
        if PlotResults:
            if ExportPlot:
                dpi = 120
            else:
                dpi = 80
    
    
    
            if TrgRoi:
                ListOfRois = [SrcRoi, TrgRoi, NewTrgRoi]
    
                ListOfDicomDirs = [SrcDcmDir, TrgDcmDir, TrgDcmDir]
    
                ListOfPlotTitles = ['Source', 'Original Target', 'New Target']
            else:
                ListOfRois = [SrcRoi, NewTrgRoi]
    
                ListOfDicomDirs = [SrcDcmDir, TrgDcmDir]
    
                ListOfPlotTitles = ['Source', 'New Target']
    
    
            #TxtToAdd = f'(RunCase = {RunCase}, \nSrcSliceNum = {SrcSliceNum}, '\
            #           + f'\nTrgSliceNum = {TrgSliceNum})'
            
            #print(f'\nListOfDicomDirs = {ListOfDicomDirs}')
            
            Times.append(time.time())
            
            if NewTrgRoi.Modality == 'RTSTRUCT':
                from PlottingTools import PlotContoursFromListOfRtss_v1
                
                PlotContoursFromListOfRtss_v1(ListOfRois, ListOfDicomDirs, 
                                              ListOfPlotTitles,
                                              PlotAllSlices=PlotAllSlices,
                                              AddTxt=TxtToAddToFname, 
                                              ExportPlot=ExportPlot, 
                                              ExportDir=RtsPlotExportDir, 
                                              dpi=dpi, LogToConsole=False)
    
            else:
                from PlottingTools import PlotPixArrsFromListOfSegs_v1
                
                PlotPixArrsFromListOfSegs_v1(ListOfRois, ListOfDicomDirs, 
                                             ListOfPlotTitles,
                                             PlotAllSlices=PlotAllSlices, 
                                             AddTxt=TxtToAddToFname, 
                                             ExportPlot=ExportPlot, 
                                             ExportDir=SegPlotExportDir,  
                                             dpi=dpi, LogToConsole=False)
                
            Times.append(time.time())
            Dtime = round(Times[-1] - Times[-2], 1)
            msg = f'Took {Dtime} s to plot the results.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
        
            TestRunTimes.append(time.time())
            Dtime = round(TestRunTimes[-1] - TestRunTimes[-2], 1)
            msg = f'Took {Dtime} s to run test {TestNum}.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
    
        
    print('All Test run iterations complete.\n')
    
    Times.append(time.time())
    Dtime = round(Times[-1] - Times[0], 1)
    msg = f'Took {Dtime} s to run all tests.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export the ListOfTimings """
        
        Fname = DateTime + f'_TestRun_{TestNum}_ListOfTimings.txt'
        
        Fpath = os.path.join(LogExportDir, Fname)
        
        ExportListToTxt(ListOfTimings, Fname, LogExportDir)
        
        print('Log file saved to:\n', Fpath, '\n')
    
    
    return 
    #return NewTrgRoi
    







def RunCopyXnatRoi(XnatSession, TestNum, LogToConsole=False, PlotResults=False, 
                   ExportPlot=False, XnatUrl=None, ProjId=None, SubjLabel=None, 
                   SrcStudyLabel=None, SrcSeriesLabel=None, SrcAsrMod=None, 
                   SrcAsrName=None, SrcRoiName=None, SrcSliceNum=None, 
                   TrgStudyLabel=None, TrgSeriesLabel=None, TrgAsrMod=None, 
                   TrgAsrName=None, TrgRoiName=None, TrgSliceNum=None, 
                   TxtToAddToTrgAsrName='', 
                   PathsDict=None, XnatDownloadDir='default', 
                   ExportLogFiles=False):
    """
    Wrapper function for CopyXnatRoi().
    
    Inputs:
    ******
    
    XnatSession : requests session (optional; None by default)
        
    TestNum : string or None
        Either a string denoting the test to run (as defined in 
        CopyRoiTestConfig.py), e.g. 'SR1', or None. 
        To run the algorithm on a different combination of data, either update 
        CopyRoiTestConfig.py or set TestNums = None and provide the necessary 
        optional inputs as described below.
    
    LogToConsole : boolean (optional if TestNums != None; False by default)
        If True, intermediate results will be logged to the console during the
        running of CopyRoi().
    
    PlotResults : boolean (optional; False by default)
    
    ExportPlot : boolean (optional; False by default)
    
    XnatUrl : string (optional if TestNums != None; None by default)
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string (optional if TestNums != None; None by default)
        The project ID of interest.
    
    SubjLabel : string (optional if TestNums != None; None by default)
        The subject label of interest. 
    
    SrcStudyLabel : string (optional if TestNums != None; None by default)
        The Source study / experiment label. 
    
    SrcSeriesLabel : string (optional if TestNums != None; None by default)
        The Source series label / scan ID.
    
    SrcAsrMod : string
        The Source assessor modality.
        
    SrcAsrName : string
        The Source assessor name (StructureSetLabel or SeriesDescription).
    
    SrcSliceNum : integer (optional; None by default; 0-indexed)
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied.
    
    TrgStudyLabel : string (optional if TestNums != None; None by default)
        The Target study / experiment label. 
    
    TrgSeriesLabel : string (optional if TestNums != None; None by default)
        The Target series label / scan ID.
    
    TrgAsrMod : string (optional but required if TrgAsrName != None; 
    None by default)
        The Target assessor modality.
        
    TrgAsrName : string (optional but required if TrgRoiMod != None; 
    None by default)
        The Target assessor name (StructureSetLabel or SeriesDescription). If 
        provided and if a direct copy is to be made, existing contours/
        segmentations for the ROI/segment will be preserved.
    
    TrgSliceNum : integer (optional unless making a direct copy; 0-indexed; 
    None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to.  This only applies for direct copies, 
        hence the default value None.
    
    TxtToAddToTrgAsrName : string (optional, '' by default)
        If provided the string of text will be appended to the assessor name 
        (StructureSetLabel or SeriesDescription), and to the filename of the 
        exported (new) Target RTS/SEG.
        
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded. If provided, paths of
        newly downloaded data will be added to PathsDict. If not provided, a
        new dictionary will be initialised.
    
    XnatDownloadDir : string (optional; 'default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjectLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    
    LogToConsole : boolean (optional; False by default)
        If True, some results will be logged to the console.
    
    ExportLogFiles : boolean (optional; False by default)
        If True, log files will be exported.
        
    
    Outputs:
    *******
    
    Times : list of floats
        A list of time stamps at certain task executions.
        
    TimingMsgs : list of strings
        A list of the time to execute certain tasks.
    
    
    Notes:
    *****
    
    If running one or more pre-configured tests, CopyRoiTestConfig.py must be
    copied to the current working directory.
    """

    import time
    import os
    from pypref import Preferences
    #from pydicom import dcmread
    #from DicomTools import IsSameModalities
    #from XnatTools import DownloadScan, DownloadAssessor
    from GeneralTools import PrintTitle, ExportListToTxt#, ExportDictionaryToJson
    
    """ Start timing. """
    Times = []
    Times.append(time.time())
    #RunDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    #TimingMsgs = [f'RunDateTime: {RunDateTime}\n']
    TimingMsgs = []
    
    
    #cwd = os.getcwd()
    #RtsExportDir = os.path.join(cwd, 'new_RTS')
    #SegExportDir = os.path.join(cwd, 'new_SEG')
    #RtsPlotExportDir = os.path.join(cwd, 'plots_RTS')
    #SegPlotExportDir = os.path.join(cwd, 'plots_SEG')
    #LogExportDir = os.path.join(cwd, 'logs')
    
    
    """ The pre-configured inputs: """
    TestCfg = Preferences(directory=os.getcwd(), 
                          filename='CopyRoiTestConfig.py')
    
    if TestNum:
        if not isinstance(TestNum, str):
            msg = 'The input argument "TestNum" must be a character string.'
            print(msg)
            raise Exception(msg)
        
        XnatUrl = TestCfg.get('XnatUrl')
        ProjId = TestCfg.get('ProjId')
        SubjLabel = TestCfg.get('SubjLabel')
        SrcStudyLabel = TestCfg.get(TestNum)['SrcStudyLabel']
        SrcSeriesLabel = TestCfg.get(TestNum)['SrcSeriesLabel']
        SrcAsrMod = TestCfg.get(TestNum)['SrcAsrMod']
        SrcAsrName = TestCfg.get(TestNum)['SrcAsrName']
        SrcRoiName = TestCfg.get(TestNum)['SrcRoiName']
        SrcSliceNum = TestCfg.get(TestNum)['SrcSliceNum']
        TrgStudyLabel = TestCfg.get(TestNum)['TrgStudyLabel']
        TrgSeriesLabel = TestCfg.get(TestNum)['TrgSeriesLabel']
        TrgAsrMod = TestCfg.get(TestNum)['TrgAsrMod']
        TrgAsrName = TestCfg.get(TestNum)['TrgAsrName']
        TrgRoiName = TestCfg.get(TestNum)['TrgRoiName']
        TrgSliceNum = TestCfg.get(TestNum)['TrgSliceNum']
        TxtToAddToTrgAsrName = TestCfg.get(TestNum)['TxtToAddToTrgAsrName']
        
        PrintTitle(f'Running Test {TestNum}:')
    else:
        if XnatUrl == None:
            XnatUrl = TestCfg.get('XnatUrl')
            
        if ProjId == None:
            raise Exception("ProjId is required.")
        
        if SubjLabel == None:
            raise Exception("SubjLabel is required.")
        
        if SrcStudyLabel == None:
            raise Exception("SrcStudyLabel is required.")
        
        if SrcSeriesLabel == None:
            raise Exception("SrcSeriesLabel is required.")
        
        if SrcAsrMod == None:
            raise Exception("SrcAsrMod is required.")
        
        if TrgStudyLabel == None:
            raise Exception("TrgStudyLabel is required.")
        
        if TrgAsrMod == None:
            TrgAsrMod = SrcAsrMod
            print("TrgAsrMod has not been provided so using SrcAsrMod:",
                  f"{SrcAsrMod}.")
    
    
    """ The pre-defined defaults for remaining inputs: """
    Defaults = Preferences(directory=os.getcwd(), 
                           filename='CopyRoiDefaults.py')
    
    ResInterp = Defaults.get('ResInterp')
    PreResVar = Defaults.get('PreResVar')
    ForceReg = Defaults.get('ForceReg')
    Tx = Defaults.get('Tx')
    TxMaxIters = Defaults.get('TxMaxIters')
    TxInterp = Defaults.get('TxInterp')
    ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
    PostTxVar = Defaults.get('PostTxVar')
    ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
    ExportNewTrgAsr = Defaults.get('ExportNewTrgAsr')
    
    print("Default input arguments read from 'CopyRoiDefaults.py':\n",
          f"   ResInterp = {ResInterp}\n   PreResVar = {PreResVar}\n",
          f"   ForceReg = {ForceReg}\n   Tx = {Tx}\n   TxMaxIters =",
          f"{TxMaxIters}\n   TxInterp = {TxInterp}\n   ApplyPostTxBlur =",
          f"{ApplyPostTxBlur}\n   PostTxVar = {PostTxVar}\n   ApplyPostTxBin =",
          f"{ApplyPostTxBin}\n")
    
    
    NewTrgAsr, PathsDict, XnatSession, DictOfInputs, ListOfInputs, TimingMsgs\
    = CopyXnatRoi(XnatUrl, XnatSession, ProjId, SubjLabel, SrcStudyLabel, 
                  SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName, SrcSliceNum, 
                  TrgStudyLabel, TrgSeriesLabel, TrgAsrMod, TrgAsrName, 
                  TrgRoiName, TrgSliceNum, TxtToAddToTrgAsrName, 
                  PathsDict, XnatDownloadDir, LogToConsole, ExportLogFiles,
                  ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, TxInterp, 
                  ApplyPostTxBlur, PostTxVar, ApplyPostTxBin)
    
    
    """ Error check the new Target RTS/SEG """
    Times.append(time.time())
    
    TrgDcmDir = PathsDict[ProjId][SubjLabel][TrgStudyLabel][TrgSeriesLabel]\
                ['DicomDir']
    
    ErrorList, Nerrors = ErrorCheckRoi(NewTrgAsr, TrgDcmDir, LogToConsole,
                                       DictOfInputs, False)
    
    Times.append(time.time())
    Dtime = round(Times[-1] - Times[-2], 1)
    msg = f'Took {Dtime} s to error check the {SrcAsrMod}.  There were '\
          + f'{Nerrors} errors found.\n'
    TimingMsgs.append(msg)
    print(f'*{msg}')
    
    
    DateTime = DictOfInputs['RunDateTime']
    
    
    if ExportLogFiles:
        """ Export log of error check results. """
        if TestNum:
            Fname = f'{DateTime}_TestRun_{TestNum}_ErrorCheckLog.txt'
        else:
            Fname = f'{DateTime}_ErrorCheckLog.txt'
        
        CWD = os.getcwd()
        ExportDir = os.path.join(CWD, 'logs')
        
        Fpath = os.path.join(ExportDir, Fname)
        
        ExportListToTxt(ErrorList, Fname, ExportDir)
            
        print('Log of error checks saved to:\n', Fpath, '\n')
    
    

    """ Export the new Target RTS/SEG """

    if ExportNewTrgAsr:# and not Nerrors:
        UseCaseThatApplies = DictOfInputs['UseCaseThatApplies']
        UseCaseToApply = DictOfInputs['UseCaseToApply']
        if TestNum:
            NewTrgAsrFname = f'TestNum_{TestNum}_{SrcAsrName}'\
                             + f'{TxtToAddToTrgAsrName}'
        else:
            NewTrgAsrFname = f'{SrcAsrName}{TxtToAddToTrgAsrName}'
        
        if ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
            NewTrgAsrFname += f'_ForcedReg_{Tx}'
            
        if not ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
            #NewTrgAsrFname += f'_{ResInterp}'
            
            """ The actual interpolation used for resampling might not be
            ResInterp. """
            DiffInterpByRoi = []
            
            for key, val in DictOfInputs.items():
                if 'ResInterpUsedFor' in key and not ResInterp in key:
                    DiffInterpByRoi.append(val)
            
            """ If DiffInterpByRoi is not empty, use the first item. """
            if DiffInterpByRoi:
                NewTrgAsrFname += f'_{DiffInterpByRoi[0]}'
            else:
                if ResInterp == 'NearestNeighbor':
                    NewTrgAsrFname += '_NN'
                else:
                    NewTrgAsrFname += f'_{ResInterp}'
            
            
            
        if UseCaseToApply in ['5a', '5b']:
            if TxInterp == 'NearestNeighbor':
                NewTrgAsrFname += '_NN'
            else:
                NewTrgAsrFname += f'_{TxInterp}'
        
        SrcAsrFpath = PathsDict[ProjId][SubjLabel][SrcStudyLabel]\
                      [SrcSeriesLabel]['assessors'][SrcAsrMod][SrcAsrName]\
                      ['AsrFpath']
        
        RtsExportDir = Defaults.get('RtsExportDir')
        SegExportDir = Defaults.get('SegExportDir')
        
        if NewTrgAsr.Modality == 'RTSTRUCT':
            NewTrgAsrFpath = ExportTrgRoi(NewTrgAsr, SrcAsrFpath, RtsExportDir, 
                                          NewTrgAsrFname, DictOfInputs)
        else:
            NewTrgAsrFpath = ExportTrgRoi(NewTrgAsr, SrcAsrFpath, SegExportDir, 
                                          NewTrgAsrFname, DictOfInputs)
    
    
    """ Plot Contours """

    if PlotResults:
        from pydicom import dcmread
        
        SrcAsr = dcmread(SrcAsrFpath)
        
        SrcDcmDir = PathsDict[ProjId][SubjLabel][SrcStudyLabel]\
                    [SrcSeriesLabel]['DicomDir']
        
        if TrgAsrName:
            TrgAsrFpath = PathsDict[ProjId][SubjLabel][TrgStudyLabel]\
                  [TrgSeriesLabel]['assessors'][TrgAsrMod][TrgAsrName]\
                  ['AsrFpath']
            
            TrgAsr = dcmread(TrgAsrFpath)
                  
            ListOfRois = [SrcAsr, TrgAsr, NewTrgAsr]

            ListOfDicomDirs = [SrcDcmDir, TrgDcmDir, TrgDcmDir]

            ListOfPlotTitles = ['Source', 'Original Target', 'New Target']
        else:
            ListOfRois = [SrcAsr, NewTrgAsr]

            ListOfDicomDirs = [SrcDcmDir, TrgDcmDir]

            ListOfPlotTitles = ['Source', 'New Target']
        
        
        
        if ExportPlot:
            dpi = 120
        else:
            dpi = 80

        RtsPlotExportDir = Defaults.get('RtsPlotExportDir')
        SegPlotExportDir = Defaults.get('SegPlotExportDir')
        PlotAllSlices = Defaults.get('PlotAllSlices')
        #if TestNum:
        #    PlotFname = f'{DateTime}_TestNum_{TestNum}_{SrcAsrName}'\
        #                + f'{TxtToAddToTrgAsrName}'
        #else:
        #    PlotFname = f'{DateTime}_{SrcAsrName}{TxtToAddToTrgAsrName}'
        
        #TxtToAdd = f'(RunCase = {RunCase}, \nSrcSliceNum = {SrcSliceNum}, '\
        #           + f'\nTrgSliceNum = {TrgSliceNum})'
        
        #print(f'\nListOfDicomDirs = {ListOfDicomDirs}')
        
        Times.append(time.time())
        
        if NewTrgAsr.Modality == 'RTSTRUCT':
            import PlottingTools
            import importlib
            importlib.reload(PlottingTools)
            from PlottingTools import PlotContoursFromListOfRtss_v1
            
            PlotContoursFromListOfRtss_v1(ListOfRois, ListOfDicomDirs, 
                                          ListOfPlotTitles,
                                          PlotAllSlices=PlotAllSlices,
                                          AddTxt=NewTrgAsrFname, 
                                          ExportPlot=ExportPlot, 
                                          ExportDir=RtsPlotExportDir, 
                                          dpi=dpi, LogToConsole=False)

        else:
            from PlottingTools import PlotPixArrsFromListOfSegs_v1
            
            PlotPixArrsFromListOfSegs_v1(ListOfRois, ListOfDicomDirs, 
                                         ListOfPlotTitles,
                                         PlotAllSlices=PlotAllSlices, 
                                         AddTxt=NewTrgAsrFname, 
                                         ExportPlot=ExportPlot, 
                                         ExportDir=SegPlotExportDir,  
                                         dpi=dpi, LogToConsole=False)
            
        Times.append(time.time())
        Dtime = round(Times[-1] - Times[-2], 1)
        msg = f'Took {Dtime} s to plot the results.\n'
        TimingMsgs.append(msg)
        print(f'*{msg}')
    
    
    Dtime = round(Times[-1] - Times[0], 1)
    msg = f'Took {Dtime} s to run test {TestNum}.\n'
    TimingMsgs.append(msg)
    print(f'*{msg}')
        
    """ change stuff below """    
    #print('All Test run iterations complete.\n')
    
    #Times.append(time.time())
    #Dtime = round(Times[-1] - Times[0], 1)
    #msg = f'Took {Dtime} s to run all tests.\n'
    #TimingMsgs.append(msg)
    #print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export the ListOfTimings """
        
        LogExportDir = Defaults.get('LogExportDir')
        
        Fname = DateTime + f'_TestRun_{TestNum}_ListOfTimings.txt'
        
        Fpath = os.path.join(LogExportDir, Fname)
        
        ExportListToTxt(TimingMsgs, Fname, LogExportDir)
        
        print('Log file saved to:\n', Fpath, '\n')
    
    
    #return Times, TimingMsgs # output to send if looping through TestNums
    #return NewTrgRoi
    return
    








"""
******************************************************************************
******************************************************************************
EXPORT A DICOM SPATIAL OR DEFORMABLE SPATIAL REGISTRATION OBJECT (SRO) TO DISK
******************************************************************************
******************************************************************************
"""


def ExportDro_OLD(SrcRoi, TrgRoi, ExportDir, Description=''):
    """
    Export a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO) to disk.  
    
    Inputs:
    ******
    
    SrcRoi : Pydicom object
        Source RTS/SEG ROI object. Any tags associated with the "moving" image
        will be assigned to the corresponding tag values in SrcRoi. 
        
    TrgRoi : Pydicom object
        The Target RTS/SEG ROI object that was created from the copying 
        process.  Any tags associated with the "fixed" image will be assigned
        to the corresponding tag values in TrgRoi.
            
    ExportDir : string
        Directory where the new DRO is to be exported.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
                           
    
    Outputs:
    *******
    
    DroFpath : string
        Full path of the exported DRO.
    """
    
    import os
    import time
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    #from copy import deepcopy
    from GeneralTools import ParseTransformParameterFile
    from DicomTools import ImportDicoms
    
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    
    # The Transform parameter file:
    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    # Parse the Transform parameter file:
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')

    TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0']
    
    
    
    # Sample DICOM Registration objects:
    DroRootDir = r'C:\Users\ctorti\Documents\DICOM Registration Object\IJ_923_dicom-sro.tar\IJ_923_dicom-sro\dicom-sro'
    
    DroDir1 = r"rect-offset-sro-rigid" # Spatial Registration object
    DroDir2 = r"sphere-centered-sro-deformable" # Deformable Spatial 
    # Registration object
    
    DroFname1 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm"
    DroFname2 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm"
    
    DroFpath1 = os.path.join(DroRootDir, DroDir1, DroFname1)
    DroFpath2 = os.path.join(DroRootDir, DroDir2, DroFname2)
    
    #Dro1 = dcmread(DroFpath1) # Spatial Registration Storage
    #Dro2 = dcmread(DroFpath2) # Deformable Spatial Registration Storage
    
    """
    If the registration was non-elastic use Dro1 as a template.
    If the registration was elastic use Dro2.
    """
    
    if True:
        Dro = dcmread(DroFpath1)
    else:
        Dro = dcmread(DroFpath2)
    
    
    CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    CurrentTime = time.strftime("%H%M%S", time.gmtime())
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    Dro.InstanceCreationData = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    if TrgRoi.Modality == 'RTSTRUCT':
        ContentDate = TrgRoi.StructureSetDate
        ContentTime = TrgRoi.StructureSetTime
        SeriesDescription = TrgRoi.StructureSetLabel
        
        TrgRefSOPClassUID = TrgRoi.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPClassUID
        
        SrcRefSOPClassUID = SrcRoi.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPClassUID
        
        
                                          
    else:
        ContentDate = TrgRoi.ContentDate
        ContentTime = TrgRoi.ContentTime
        SeriesDescription = TrgRoi.SeriesDescription
        
        TrgRefSOPClassUID = TrgRoi.ReferencedSeriesSequence[0]\
                                  .ReferencedInstanceSequence[0]\
                                  .ReferencedSOPClassUID
        
        SrcRefSOPClassUID = SrcRoi.ReferencedSeriesSequence[0]\
                                  .ReferencedInstanceSequence[0]\
                                  .ReferencedSOPClassUID
        
        
    
    
    Dro.StudyDate = TrgRoi.StudyDate
    Dro.SeriesDate = TrgRoi.SeriesDate
    Dro.ContentDate = ContentDate
    Dro.StudyTime = TrgRoi.StudyTime
    Dro.SeriesTime = TrgRoi.SeriesTime
    Dro.ContentTime = ContentTime
    Dro.Manufacturer = TrgRoi.Manufacturer
    """Consider modifying StudyDescription (= '' in the template DRO."""
    Dro.SeriesDescription = SeriesDescription
    Dro.ManufacturerModelName = TrgRoi.ManufacturerModelName
    Dro.PatientName = TrgRoi.PatientName
    Dro.PatientID = TrgRoi.PatientID
    Dro.PatientBirthDate = TrgRoi.PatientBirthDate
    Dro.PatientSex = TrgRoi.PatientSex
    Dro.PatientAge = TrgRoi.PatientAge
    Dro.StudyInstanceUID = TrgRoi.StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    Dro.SeriesInstanceUID = generate_uid()
    Dro.StudyID = TrgRoi.StudyID
    Dro.SeriesNumber = TrgRoi.SeriesNumber
    Dro.InstanceNumber = TrgRoi.InstanceNumber
    Dro.FrameOfReferenceUID = TrgRoi.FrameOfReferenceUID
    Dro.PositionReferenceIndicator = TrgRoi.PositionReferenceIndicator
    """Keep ContentLabel as 'REGISTRATION'."""
    Dro.ContentDescription = Description
    if TrgRoi.Modality == 'SEG':
        Dro.ContentCreatorName = TrgRoi.ContentCreatorName
    
    
    
    """ Modify the RegistrationSequence for the fixed domain. """
    
    Dro.RegistrationSequence[0]\
       .FrameOfReferenceUID = TrgRoi.FrameOfReferenceUID
    
    
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' = value in template DRO
    'RIGID_SCALE'
    'AFFINE'    
    
    Need to cover 'RIGID_SCALE' and 'AFFINE' cases...
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = 
    """
    
    """
    Modify the FrameOfReferenceTransformationMatrix.  No change required to
     that of the template DRO is the identity transformation is desired, i.e.
    ['1.0', '0.0', '0.0', '0.0', 
     '0.0', '1.0', '0.0', '0.0', 
     '0.0', '0.0', '1.0', '0.0', 
     '0.0', '0.0', '0.0', '1.0']
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = 
    """
    
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the RegistrationSequence for the moving domain. """
    
    Dro.RegistrationSequence[1]\
       .FrameOfReferenceUID = SrcRoi.FrameOfReferenceUID
    
    
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' = value in template DRO
    'RIGID_SCALE'
    'AFFINE'    
    
    Need to cover 'RIGID_SCALE' and 'AFFINE' cases...
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = 
    """
    
    """ Modify the FrameOfReferenceTransformationMatrix. """
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxParams
    
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ continue with Series and Instance Reference Macro attributes, as 
    covered in 10.3 (page 15) of the DICOM Supplement 73. """
    
    
    """ Modify the ReferencedSeriesSequence for the fixed domain. """
    
    Dro.ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = TrgRefSOPClassUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the first 
    ReferencedSeriesSequence referenced the StudyInstanceUID of the fixed 
    series."""
    Dro.ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = TrgRoi.StudyInstanceUID
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgRoi.SeriesInstanceUID
    
    
    
    
    
    """ Modify the ReferencedSeriesSequence for the moving domain. """
    
    Dro.ReferencedSeriesSequence[1]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = SrcRefSOPClassUID
    
    Dro.ReferencedSeriesSequence[1]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = SrcRoi.StudyInstanceUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the second 
    ReferencedSeriesSequence referenced the SOPInstanceUID of 31st DICOM (in a
    collection of 100 DICOMs) of the moving series.
    Dro.ReferencedSeriesSequence[1]\
       .SeriesInstanceUID = SrcRoi.?????????
    """
    
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    fixed domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = TrgRefSOPClassUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the first 
    StudiesContainingOtherReferencedInstancesSequence referenced the 
    SOPInstanceUID of the last DICOM in the fixed series.
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = ?????
    """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgRoi.SeriesInstanceUID
    
    
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    moving domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = SrcRefSOPClassUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the second 
    StudiesContainingOtherReferencedInstancesSequence referenced the 
    SOPInstanceUID of the last DICOM in the moving series.
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = ?????
    """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = SrcRoi.SeriesInstanceUID
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    """
    Delete all GroupLength tags?
    """
    
    
    DroFname = CurrentDateTime + '_' + Description.replace(' ', '_') + '.dcm'
    
    DroFpath = os.path.join(ExportDir, DroFname)
    
    Dro.save_as(DroFpath)
        
    print('\nNew DICOM Registration Object exported to:\n', DroFpath)
    
    return DroFpath





def CreateDro(SrcDicomDir, TrgDicomDir, Description='', LogToConsole=False):
    """
    Create a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO).  
    
    *COMMENT*:  At the moment only the Spatial Registration Object Storage IOD
    is covered.
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    """
    
    import os
    import time
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    from copy import deepcopy
    from DicomTools import ImportDicoms
    from GeneralTools import ParseTransformParameterFile
    from GeneralTools import GetTxMatrixType
    
    # Start timing:
    times = []
    times.append(time.time())
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
        
    # The Transform parameter file:
    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    # Parse the Transform parameter file:
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n*Took {Dtime} s to parse the transform parameter file.')

    TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0']
    
    
    
    # Sample DICOM Registration objects:
    DroRootDir = r'C:\Users\ctorti\Documents\DICOM Registration Object\IJ_923_dicom-sro.tar\IJ_923_dicom-sro\dicom-sro'
    
    DroDir1 = r"rect-offset-sro-rigid" # Spatial Registration object
    DroDir2 = r"sphere-centered-sro-deformable" # Deformable Spatial 
    # Registration object
    
    DroFname1 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm"
    DroFname2 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm"
    
    DroFpath1 = os.path.join(DroRootDir, DroDir1, DroFname1)
    DroFpath2 = os.path.join(DroRootDir, DroDir2, DroFname2)
    
    #Dro1 = dcmread(DroFpath1) # Spatial Registration Storage
    #Dro2 = dcmread(DroFpath2) # Deformable Spatial Registration Storage
    
    """
    If the registration was non-elastic use Dro1 as a template.
    If the registration was elastic use Dro2.
    """
    
    if True:
        Dro = dcmread(DroFpath1)
    else:
        Dro = dcmread(DroFpath2)
    
    
    CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    CurrentTime = time.strftime("%H%M%S", time.gmtime())
    #CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    Dro.InstanceCreationData = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    
    Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.StudyDate = TrgDicoms[0].StudyDate
    Dro.SeriesDate = TrgDicoms[0].SeriesDate
    Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.StudyTime = TrgDicoms[0].StudyTime
    Dro.SeriesTime = TrgDicoms[0].SeriesTime
    Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.Manufacturer = TrgDicoms[0].Manufacturer
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    Dro.SeriesDescription = TrgDicoms[0].SeriesDescription
    Dro.ManufacturerModelName = TrgDicoms[0].ManufacturerModelName
    Dro.PatientName = TrgDicoms[0].PatientName
    Dro.PatientID = TrgDicoms[0].PatientID
    Dro.PatientBirthDate = TrgDicoms[0].PatientBirthDate
    Dro.PatientSex = TrgDicoms[0].PatientSex
    Dro.PatientAge = TrgDicoms[0].PatientAge
    Dro.StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    Dro.SeriesInstanceUID = generate_uid()
    Dro.StudyID = TrgDicoms[0].StudyID
    Dro.SeriesNumber = TrgDicoms[0].SeriesNumber
    Dro.InstanceNumber = TrgDicoms[0].InstanceNumber
    Dro.FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
    Dro.PositionReferenceIndicator = TrgDicoms[0].PositionReferenceIndicator
    """ Keep ContentLabel as 'REGISTRATION'. """
    """ Keep ContentCreatorName as ''. """
    
    
    
    """ Modify the RegistrationSequence for the fixed domain. """
    
    Dro.RegistrationSequence[0]\
       .FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
    
    """
    Since the first sequence in RegistrationSequence will be the fixed image
    domain, FrameOfReferenceTransformationMatrix will be the identity matrix
    ['1.0', '0.0', '0.0', '0.0', 
     '0.0', '1.0', '0.0', '0.0', 
     '0.0', '0.0', '1.0', '0.0', 
     '0.0', '0.0', '0.0', '1.0'].
    
    Hence no change is required FrameOfReferenceTransformationMatrix for this
    sequence.
    """
    
    
    """
    Since the FrameOfReferenceTransformationMatrix is the identity matrix, the 
    FrameOfReferenceTransformationMatrixType will be 'RIGID'.  Hence no change
    is required to FrameOfReferenceTransformationMatrixType for this sequence.
    """
    
    
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the RegistrationSequence for the moving domain. """
    
    Dro.RegistrationSequence[1]\
       .FrameOfReferenceUID = SrcDicoms[0].FrameOfReferenceUID
    
    """ Modify the FrameOfReferenceTransformationMatrix. """
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxParams
       
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' = value in template DRO
    'RIGID_SCALE'
    'AFFINE'    
    """
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to get the transform matrix type.')
        
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the ReferencedSeriesSequence for the fixed domain. """
    
    print('\nModify the ReferencedSeriesSequence for the fixed domain.')
    
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(TrgDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in ReferencedSeriesSequence for',
              'the fixed series.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = TrgDicoms[i].SOPClassUID
           
        Dro.ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = TrgDicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the fixed series.')
    
    
    
    
    """ Modify the ReferencedSeriesSequence for the moving domain. """
    
    print('\nModify the ReferencedSeriesSequence for the moving domain.')
    
    Dro.ReferencedSeriesSequence[1]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[1]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(SrcDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.ReferencedSeriesSequence[1]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.ReferencedSeriesSequence[1]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.ReferencedSeriesSequence[1]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in ReferencedSeriesSequence for',
              'the moving series.')
        
    # Modify the sequences:
    for i in range(ReqRIS):
        Dro.ReferencedSeriesSequence[1]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = SrcDicoms[i].SOPClassUID
           
        Dro.ReferencedSeriesSequence[1]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = SrcDicoms[i].SOPInstanceUID

    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the moving series.')
        
    
    print('\nModify the StudiesContainingOtherReferencedInstancesSequence',
          'for the fixed domain.')
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    fixed domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
       
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    # The number of sequences:
    OrigRIS = len(Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
                     .ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(TrgDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
                                   .ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'fixed series.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = TrgDicoms[i].SOPClassUID
           
        Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = TrgDicoms[i].SOPInstanceUID
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to modify the',
              'ReferencedInstanceSequences in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'fixed series.')
        
        
    print('\nModify the StudiesContainingOtherReferencedInstancesSequence for',
          'the moving domain.')
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    moving domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .StudyInstanceUID = SrcDicoms[0].StudyInstanceUID
       
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
    
    # The number of sequences:
    OrigRIS = len(Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
                     .ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(SrcDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
                                   .ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'moving series.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = SrcDicoms[i].SOPClassUID
           
        Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = SrcDicoms[i].SOPInstanceUID
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to modify the',
              'ReferencedInstanceSequences in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'moving series.')
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    print('Deleting GroupLength tags.')
    
    """
    The SRO template contains 4 retired tags with name "GroupLength". 
    Delete them.
    """
    
    del Dro[0x08, 0x00]
    del Dro[0x10, 0x00]
    del Dro[0x20, 0x00]
    del Dro[0x70, 0x00]
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to delete all GroupLength tags.')
        
    
    Dtime = round(times[-1] - times[0], 1)
    if True:#LogToConsole:
        print(f'\n   *Took a total of {Dtime} s to create the DRO.')
        
        
    return Dro





def ExportDro(Dro, ExportDir, Description='', LogToConsole=False):
    """
    Export a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO) to disk.  
    
    *COMMENT*:  At the moment only the Spatial Registration Object Storage IOD
    is covered.
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Registration Object to export.
        
    ExportDir : string
        Directory where the new DRO is to be exported.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    DroFpath : string
        Full path of the exported DRO.
    """
    
    import os
    import time
    
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    DroFname = CurrentDateTime + '_' + Description.replace(' ', '_') + '.dcm'
    
    DroFpath = os.path.join(ExportDir, DroFname)
    
    print('Exporting DRO.')
    
    Dro.save_as(DroFpath)
        
    print('\nNew DICOM Registration Object exported to:\n', DroFpath)
    
    return DroFpath