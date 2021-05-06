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
    
    TrgRoi : Pydicom object or None
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
CREATE A LIST OF INPUTS AND A DICTIONARY OF INPUTS FOR COPYROI()
******************************************************************************
******************************************************************************
"""

def CreateListOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, 
                       SrcExpLabel, SrcScanId, SrcSliceNum,
                       SrcRoiColMod, SrcRoiColName, SrcRoiName,  
                       TrgExpLabel, TrgScanId, TrgSliceNum, 
                       TrgRoiColMod, TrgRoiColName, TrgRoiName, 
                       TxtToAddToTrgRoiColName, 
                       LogToConsole, ExportLogFiles, 
                       ResInterp, PreResVariance, ApplyPostResBlur, 
                       PostResVariance, ForceReg, TxMatrix, SelxOrSitk, 
                       Transform, MaxIters, TxInterp, ApplyPostTxBin, 
                       ApplyPostTxBlur, PostTxVariance):
    
    ListOfInputs = [f"RunDateTime: {RunDateTime}\n",
                    f"XnatUrl = {XnatUrl}\n",
                    f"ProjId = {ProjId}\n",
                    f"SubjLabel = {SubjLabel}\n",
                    f"SrcExpLabel = {SrcExpLabel}\n",
                    f"SrcScanId = {SrcScanId}\n",
                    f"SrcSliceNum = {SrcSliceNum}\n",
                    f"SrcRoiColMod = {SrcRoiColMod}\n",
                    f"SrcRoiColName = {SrcRoiColName}\n",
                    f"SrcRoiName = {SrcRoiName}\n",
                    f"TrgExpLabel = {TrgExpLabel}\n",
                    f"TrgScanId = {TrgScanId}\n",
                    f"TrgSliceNum = {TrgSliceNum}\n",
                    f"TrgRoiColMod = {TrgRoiColMod}\n",
                    f"TrgRoiColName = {TrgRoiColName}\n",
                    f"TrgRoiName = {TrgRoiName}\n",
                    f"TxtToAddToTrgRoiColName = {TxtToAddToTrgRoiColName}\n",
                    f"LogToConsole = {LogToConsole}",
                    f"ExportLogFiles = {ExportLogFiles}",
                    f"ResInterpSet = {ResInterp}\n",
                    f"PreResVariance = {PreResVariance}\n",
                    f"ApplyPostResBlur = {ApplyPostResBlur}\n",
                    f"PostResVariance = {PostResVariance}\n",
                    f"ForceReg = {ForceReg}\n",
                    f"TxMatrix = {TxMatrix}\n",
                    f"SelxOrSitk = {SelxOrSitk}\n",
                    f"Transform = {Transform}\n",
                    f"MaxIters = {MaxIters}\n",
                    f"TxInterp = {TxInterp}\n",
                    f"ApplyPostTxBin = {ApplyPostTxBin}\n",
                    f"ApplyPostTxBlur = {ApplyPostTxBlur}\n",
                    f"PostTxVariance = {PostTxVariance}\n"
                    ]
    
    return ListOfInputs





def CreateDictOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, 
                       SrcExpLabel, SrcScanId, SrcSliceNum,
                       SrcRoiColMod, SrcRoiColName, SrcRoiName,  
                       TrgExpLabel, TrgScanId, TrgSliceNum, 
                       TrgRoiColMod, TrgRoiColName, TrgRoiName, 
                       TxtToAddToTrgRoiColName, 
                       LogToConsole, ExportLogFiles, 
                       ResInterp, PreResVariance, ApplyPostResBlur, 
                       PostResVariance, ForceReg, TxMatrix, SelxOrSitk, 
                       Transform, MaxIters, TxInterp, ApplyPostTxBin, 
                       ApplyPostTxBlur, PostTxVariance):
    
    DictOfInputs = {'RunDateTime' : RunDateTime,
                    'XnatUrl' : XnatUrl,
                    'ProjId' : ProjId,
                    'SubjLabel' : SubjLabel,
                    'SrcExpLabel' : SrcExpLabel,
                    'SrcScanId' : SrcScanId,
                    'SrcSliceNum' : SrcSliceNum,
                    'SrcRoiColMod' : SrcRoiColMod,
                    'SrcRoiColName' : SrcRoiColName,
                    'SrcRoiName' : SrcRoiName,
                    'TrgExpLabel' : TrgExpLabel,
                    'TrgScanId' : TrgScanId,
                    'TrgSliceNum' : TrgSliceNum,
                    'TrgRoiColMod' : TrgRoiColMod,
                    'TrgRoiColName' : TrgRoiColName,
                    'TrgRoiName' : TrgRoiName,
                    'TxtToAddToTrgRoiColName' : TxtToAddToTrgRoiColName,
                    'LogToConsole' : LogToConsole,
                    'ExportLogFiles' : ExportLogFiles,
                    'ResInterpSet' : ResInterp,
                    'PreResVariance' : PreResVariance,
                    'ApplyPostResBlur' : ApplyPostResBlur,
                    'PostResVariance' : PostResVariance,
                    'ForceReg' : ForceReg,
                    'TxMatrix' : TxMatrix,
                    'SelxOrSitk' : SelxOrSitk,
                    'Transform' : Transform,
                    'MaxIters' : MaxIters,
                    'TxInterp' : TxInterp,
                    'ApplyPostTxBin' : ApplyPostTxBin,
                    'ApplyPostTxBlur' : ApplyPostTxBlur,
                    'PostTxVariance' : PostTxVariance
                    }
    
    return DictOfInputs

    
    
    

"""
******************************************************************************
******************************************************************************
COPY A CONTOUR / ROI / RTSTRUCT / SEGMENTATION / SEGMENT / SEG
******************************************************************************
******************************************************************************
"""


def CopyLocalRoiCol(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, 
                    TrgDcmDir, TrgRoiFpath=None, TrgRoiName=None, TrgSliceNum=None, 
                    ResInterp='BlurThenLinear', PreResVariance=(1,1,1), 
                    ApplyPostResBlur=False, PostResVariance=(1,1,1),
                    ForceReg=False, SelxOrSitk='Selx', Transform='affine', 
                    MaxIters='512', TxInterp='NearestNeighbor', 
                    ApplyPostTxBin=True, ApplyPostTxBlur=True, 
                    PostTxVariance=(1,1,1), TxtToAddToTrgRoiColName='', 
                    LogToConsole=False, ExportLogFiles=False):
    """
    Copy an ROI Collection (contour/ROI/RTSTRUCT/segmentation/segment/SEG) from 
    locally sourced data.
        
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
    
    PreResVariance : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
        
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    SelxOrSitk : string (optional; 'Selx' by default)
        Denotes which package to use for image registration and transformation.
        Acceptable values include:
            - 'Selx' for SimpleElastix
            - 'Sitk' for SimpleITK
    
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Transform. If != 'default' it must be a string 
        representation of an integer.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
        
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    TxtToAddToTrgRoiColName : string (optional, '' by default)
        String of text to pass to CreateRts()/CreateSeg(). The string will be
        appended to the RTS StructureSetLabel or SEG SeriesDescription, and to 
        the filename of the exported (new) Target RTS/SEG.
            
    LogToConsole : boolean (default False)
        If True, some results will be logged to the console.
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
          
                  
    Outputs:
    *******
        
    TrgRoiCol : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing RTS/SEG) 
        Target RTS/SEG object.
    
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgRoi, None otherwise.
    
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
    from RtsTools import CopyRts
    from SegTools import CopySeg
    
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
                    f'PreResVariance = {PreResVariance}\n',
                    f'ApplyPostResBlur = {ApplyPostResBlur}\n',
                    f'PostResVariance = {PostResVariance}\n',
                    f'ForceReg = {ForceReg}\n',
                    f'SelxOrSitk = {SelxOrSitk}\n',
                    f'Transform = {Transform}\n',
                    f'MaxIters = {MaxIters}\n',
                    f'TxInterp = {TxInterp}\n',
                    f'ApplyPostTxBin = {ApplyPostTxBin}\n',
                    f'ApplyPostTxBlur = {ApplyPostTxBlur}\n',
                    f'PostTxVariance = {PostTxVariance}\n',
                    f'TxtToAddToTrgRoiColName = {TxtToAddToTrgRoiColName}\n',
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
                    'PreResVariance' : PreResVariance,
                    'ApplyPostResBlur' : ApplyPostResBlur,
                    'PostResVariance' : PostResVariance,
                    'ForceReg' : ForceReg,
                    'SelxOrSitk' : SelxOrSitk,
                    'Transform' : Transform,
                    'MaxIters' : MaxIters,
                    'TxInterp' : TxInterp,
                    'ApplyPostTxBin' : ApplyPostTxBin,
                    'ApplyPostTxBlur' : ApplyPostTxBlur,
                    'PostTxVariance' : PostTxVariance,
                    'TxtToAddToTrgRoiColName' : TxtToAddToTrgRoiColName,
                    'LogToConsole' : LogToConsole}
    
    
    ListOfTimings = [f'RunDateTime: {RunDateTime}\n']
    
    times = []
    
    
    """ Start timing. """
    times.append(time.time())
    
    SrcRoiCol = dcmread(SrcRoiFpath)
    if TrgRoiFpath:
        TrgRoiCol = dcmread(TrgRoiFpath)
    else:
        TrgRoiCol = None
    
    """ Establish whether the inputs are valid. """
    CheckValidityOfInputs(SrcRoiCol, SrcRoiName, SrcSliceNum, 
                          TrgRoiCol, TrgRoiName, TrgSliceNum, TrgRoiFpath, 
                          LogToConsole)
    
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
        TrgRoiCol, Dro, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyRts(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, 
                  TrgDcmDir, UseCaseToApply, TrgRoiFpath, TrgRoiName, 
                  TrgSliceNum, ResInterp, PreResVariance, ApplyPostResBlur, 
                  PostResVariance, ForceReg, SelxOrSitk, Transform, MaxIters, 
                  TxInterp, ApplyPostTxBin, ApplyPostTxBlur, PostTxVariance,    
                  TxtToAddToTrgRoiColName, LogToConsole,
                  DictOfInputs, ListOfInputs, ListOfTimings)
    
    elif Modality == 'SEG':
        TrgRoiCol, Dro, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopySeg(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, 
                  TrgDcmDir, UseCaseToApply, TrgRoiFpath, TrgRoiName, 
                  TrgSliceNum, ResInterp, PreResVariance, ApplyPostResBlur, 
                  PostResVariance, ForceReg, SelxOrSitk, Transform, MaxIters, 
                  TxInterp, ApplyPostTxBin, ApplyPostTxBlur, PostTxVariance,    
                  TxtToAddToTrgRoiColName, LogToConsole,
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
        
    return TrgRoiCol, Dro, DictOfInputs, ListOfInputs, ListOfTimings






def CopyXnatRoiCol(XnatUrl, XnatSession, ProjId, SubjLabel, 
                   SrcExpLabel, SrcScanId, SrcSliceNum,
                   SrcRoiColMod, SrcRoiColName, SrcRoiName, 
                   TrgExpLabel, TrgScanId, TrgSliceNum=None, 
                   TrgRoiColMod=None, TrgRoiColName=None, TrgRoiName=None, 
                   TxtToAddToTrgRoiColName='',
                   PathsDict=None, XnatDownloadDir='default', 
                   LogToConsole=False, ExportLogFiles=False,
                   ResInterp='BlurThenLinear', PreResVariance=(1,1,1), 
                   ApplyPostResBlur=False, PostResVariance=(1,1,1),
                   ForceReg=False, TxMatrix=None, SelxOrSitk='Selx',  
                   Transform='affine', MaxIters='512', 
                   TxInterp='NearestNeighbor', ApplyPostTxBin=True, 
                   ApplyPostTxBlur=True, PostTxVariance=(1,1,1)):
    
    """
    Copy a contour/ROI/RTSTRUCT/segmentation/segment/SEG from data downloaded
    from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
        
    XnatSession : requests session or None
        
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    SrcExpLabel : string
        The Source DICOM study / XNAT experiment label. 
    
    SrcScanId : string
        The Source DICOM series label / XNAT scan ID.
    
    SrcSliceNum : integer (0-indexed)
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied.
    
    SrcRoiColMod : string
        The modality of the Source ROI Collection.
        
    SrcRoiColName : string
        The Source ROI Collection name (StructureSetLabel or SeriesDescription).
        
    SrcRoiName : string
        The Source ROIName or SegmentLabel.
    
    TrgExpLabel : string
        The Target DICOM study / XNAT experiment label. 
    
    TrgScanId : string
        The Target DICOM series label / XNAT scan ID.
    
    TrgSliceNum : integer (optional unless making a direct copy; 0-indexed; 
    None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to.  This only applies for direct copies, 
        hence the default value None.
    
    TrgRoiColMod : string (optional but required if TrgRoiColName != None; 
    None by default)
        The Target image assessor modality.
        
    TrgRoiColName : string (optional but required if TrgRoiColMod != None; 
    None by default)
        The Target ROI Collection name (StructureSetLabel or SeriesDescription). 
        If provided and if a direct copy is to be made, existing contours/
        segmentations for the ROI/segment will be preserved.
    
    TrgRoiName : string (optional; None by default)
        The Target ROIName or SegmentLabel.
    
    TxtToAddToTrgRoiColName : string (optional, '' by default)
        If provided the string of text will be appended to the ROI Collection 
        name (StructureSetLabel or SeriesDescription), and to the filename of 
        the exported (new) Target RTS/SEG ROI Collection.
    
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
    
    PreResVariance : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
        
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    TxMatrix : None or list of float strings (optional; None by default)
        List of float strings representing the transformation that would 
        tranform the moving (Source) image to the fixed (Target) image. 
        If not None, image registration will be skipped to save computational
        time.
    
    SelxOrSitk : string (optional; 'Selx' by default)
        Denotes which package to use for image registration and transformation.
        Acceptable values include:
            - 'Selx' for SimpleElastix
            - 'Sitk' for SimpleITK
     
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Transform. If != 'default' it must be a string 
        representation of an integer.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
        
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        tranformed labelmap image(s).
        
    
    
                  
    Outputs:
    *******
        
    TrgRoiCol : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing RTS/SEG) 
        Target ROI Collection (RTS/SEG) object.
    
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgSeg, None otherwise.
        
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
    from XnatTools import CreateXnatSession, DownloadScan, DownloadImAsr
    from XnatTools import GetFnameAndIdFromName
    from GeneralTools import ExportDictionaryToJson, ExportListToTxt
    import RtsTools, SegTools
    importlib.reload(RtsTools)
    importlib.reload(SegTools)
    from RtsTools import CopyRts
    from SegTools import CopySeg
    
    RunDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #Defaults = Preferences(directory=os.getcwd(), filename='CopyRoiDefaults.py')
    #
    #ResInterp = Defaults.get('ResInterp')
    #PreResVariance = Defaults.get('PreResVariance')
    #ForceReg = Defaults.get('ForceReg')
    #Tx = Defaults.get('Tx')
    #TxMaxIters = Defaults.get('TxMaxIters')
    #TxInterp = Defaults.get('TxInterp'),
    #ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
    #PostTxVar = Defaults.get('PostTxVar')
    #ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
    
    ListOfInputs\
    = CreateListOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, 
                         SrcExpLabel, SrcScanId, SrcSliceNum, 
                         SrcRoiColMod, SrcRoiColName, SrcRoiName, 
                         TrgExpLabel, TrgScanId, TrgSliceNum,
                         TrgRoiColMod, TrgRoiColName, TrgRoiName,  
                         TxtToAddToTrgRoiColName, LogToConsole, ExportLogFiles, 
                         ResInterp, PreResVariance, ApplyPostResBlur, 
                         PostResVariance, ForceReg, TxMatrix, SelxOrSitk, 
                         Transform, MaxIters, TxInterp, ApplyPostTxBin, 
                         ApplyPostTxBlur, PostTxVariance)
    
    DictOfInputs\
    = CreateDictOfInputs(RunDateTime, XnatUrl, ProjId, SubjLabel, 
                         SrcExpLabel, SrcScanId, SrcSliceNum, 
                         SrcRoiColMod, SrcRoiColName, SrcRoiName, 
                         TrgExpLabel, TrgScanId, TrgSliceNum,
                         TrgRoiColMod, TrgRoiColName, TrgRoiName,  
                         TxtToAddToTrgRoiColName, LogToConsole, ExportLogFiles, 
                         ResInterp, PreResVariance, ApplyPostResBlur, 
                         PostResVariance, ForceReg, TxMatrix, SelxOrSitk, 
                         Transform, MaxIters, TxInterp, ApplyPostTxBin, 
                         ApplyPostTxBlur, PostTxVariance)
    
    
    TimingMsgs = [f'RunDateTime: {RunDateTime}\n']
    
    times = []
    
    """ Create XNAT session. """
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl)
    
    
    """ Start timing. """
    times.append(time.time())
    
    """ Download data from XNAT. """
    PathsDict, XnatSession\
    = DownloadScan(XnatUrl, ProjId, SubjLabel, SrcExpLabel, SrcScanId, 
                   XnatSession, XnatDownloadDir, PathsDict)
    
    PathsDict, XnatSession\
    = DownloadScan(XnatUrl, ProjId, SubjLabel, TrgExpLabel, TrgScanId, 
                   XnatSession, XnatDownloadDir, PathsDict)

    PathsDict, XnatSession\
    = DownloadImAsr(XnatUrl, ProjId, SubjLabel, SrcExpLabel, 
                    SrcScanId, SrcRoiColMod, SrcRoiColName, 
                    XnatSession, XnatDownloadDir, PathsDict)
    
    if TrgRoiColName != None:
        PathsDict, XnatSession\
        = DownloadImAsr(XnatUrl, ProjId, SubjLabel, TrgExpLabel, 
                        TrgScanId, TrgRoiColMod, TrgRoiColName, 
                        XnatSession, XnatDownloadDir, PathsDict)
        
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to download data from XNAT.\n'
        TimingMsgs.append(msg)
        print(f'*{msg}')
    
    
    #print('\n', PathsDict, '\n')
    
    """ Get the DICOM directories and ROI Collection file paths: """
    SrcDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][SrcExpLabel]['scans'][SrcScanId]\
                ['resources']['DICOM']['files']['Dir']
    
    TrgDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][TrgExpLabel]['scans'][TrgScanId]\
                ['resources']['DICOM']['files']['Dir']
    
    SrcRoiColFname, SrcAsrId\
    = GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, SrcExpLabel,
                            SrcRoiColMod, SrcRoiColName)
    
    SrcRoiColFpath = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                     ['experiments'][SrcExpLabel]['assessors'][SrcAsrId]\
                     ['resources'][SrcRoiColMod]['files'][SrcRoiColFname]\
                     ['Fpath']
    
    if TrgRoiColName == None:
        TrgRoiColFpath = None
    else:
        TrgRoiColFname, TrgAsrId\
        = GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, TrgExpLabel,
                                TrgRoiColMod, TrgRoiColName)
    
        TrgRoiColFpath = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                         ['experiments'][TrgExpLabel]['assessors'][TrgAsrId]\
                         ['resources'][TrgRoiColMod]['files'][TrgRoiColFname]\
                         ['Fpath']
    
    
    #""" Establish whether the inputs are valid. """
    #CheckValidityOfInputs(SrcRoiName, SrcSliceNum, TrgRoiName, TrgSliceNum, 
    #                      LogToConsole)
    
    """ Determine which Use Case to apply. """
    UseCaseThatApplies,\
    UseCaseToApply = WhichUseCase(SrcSliceNum, TrgSliceNum, SrcDcmDir, 
                                  TrgDcmDir, ForceReg, LogToConsole=True)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        #print(f'\n\nDone.  Took {Dtime} s to run.')
        msg = f'Took {Dtime} s to determine which UseCase applies.\n'
        TimingMsgs.append(msg)
        print(f'*{msg}')
        
    
    DictOfInputs['UseCaseThatApplies'] = UseCaseThatApplies
    DictOfInputs['UseCaseToApply'] = UseCaseToApply
    
    if SrcRoiColMod == 'RTSTRUCT':        
        TrgRoiCol, Dro, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyRts(SrcRoiColFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir, 
                  UseCaseToApply, TrgRoiColFpath, TrgRoiName, TrgSliceNum, 
                  ResInterp, PreResVariance, ApplyPostResBlur, PostResVariance, 
                  ForceReg, TxMatrix, SelxOrSitk, Transform, MaxIters, TxInterp, 
                  ApplyPostTxBin, ApplyPostTxBlur, PostTxVariance, 
                  TxtToAddToTrgRoiColName, LogToConsole,
                  DictOfInputs, ListOfInputs, TimingMsgs)
    
    elif SrcRoiColMod == 'SEG':
        TrgRoiCol, Dro, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopySeg(SrcRoiColFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir, 
                  UseCaseToApply, TrgRoiColFpath, TrgRoiName, TrgSliceNum, 
                  ResInterp, PreResVariance, ApplyPostResBlur, PostResVariance, 
                  ForceReg, TxMatrix, SelxOrSitk, Transform, MaxIters, TxInterp, 
                  ApplyPostTxBin, ApplyPostTxBlur, PostTxVariance, 
                  TxtToAddToTrgRoiColName, LogToConsole,
                  DictOfInputs, ListOfInputs, TimingMsgs)
        
    else:
        msg = f'The modality of the Source ROI Collection ({SrcRoiColMod}) '\
              + 'must be "RTSTRUCT" or "SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to copy the ROI(s) and create the '\
              + f'{SrcRoiColMod} object.\n'  
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
        
    return TrgRoiCol, Dro, PathsDict, XnatSession, DictOfInputs,\
           ListOfInputs, TimingMsgs







"""
******************************************************************************
******************************************************************************
CHECK FOR ERRORS IN NEW TARGET RTS / SEG
******************************************************************************
******************************************************************************
"""

def ErrorCheckRoiCol(RoiCol, DicomDir, LogToConsole=False, DictOfInputs=None,
                     ExportLogFiles=False):
    """
    Check a ROI Collection (RTS/SEG) for errors in dependencies based on 
    provided directory of the DICOMs that relate to the RTS/SEG.  
    
    Inputs:
    ******
    
    RoiCol : Pydicom object
        RTS/SEG ROI Collection to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to RoiCol.
        
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
    
    Mod = RoiCol.Modality
    
    """ Start timing. """
    times = []
    times.append(time.time())
    
    if Mod == 'RTSTRUCT':
        #import RtsTools
        #import importlib
        #importlib.reload(RtsTools)
        from RtsTools import ErrorCheckRts
        
        LogList, Nerrors = ErrorCheckRts(RoiCol, DicomDir, LogToConsole)
        
    elif Mod == 'SEG':
        #import importlib
        #import SegTools
        #importlib.reload(SegTools)
        from SegTools import ErrorCheckSeg
        
        LogList, Nerrors = ErrorCheckSeg(RoiCol, DicomDir, LogToConsole)
        
    else:
        msg = f'The modality ({Mod}) must be either "RTS" or "SEG".'
        
        raise Exception(msg)
        
        LogList = None
        Nerrors = None
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to error check the {Mod} object.\n')
    
    
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

def ExportTrgRoiCol(TrgRoiCol, SrcRoiColFpath, ExportDir, Fname='', 
                    DictOfInputs=None):
    """
    Export ROI Collection (RTS/SEG) to disk.  
    
    Inputs:
    ******
    
    TrgRoiCol : Pydicom object
        Target ROI Collection (RTS/SEG) to be exported.
        
    SrcRoiColFpath : string
        Full path of the Source RTS/SEG file (used to generate the filename of
        the new RTS/SEG file).
                              
    ExportDir : string
        Directory where the new RTS/SEG is to be exported.
        
    Fname : string
        File name to assign to new assessor.
    
    DictOfInputs : dictionary or None (optional; None by default)
        If not None, a dictionary containing the inputs that were called to 
        CopyRoi().
                           
    
    Outputs:
    *******
    
    TrgRoiColFpath : string
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
        TrgRoiColFname = DateTime + '.dcm'
    else:
        #TrgRoiFname = Fname.replace(' ', '_') + '.dcm'
        TrgRoiColFname = DateTime + '_' + Fname.replace(' ', '_') + '.dcm'
    
    TrgRoiColFpath = os.path.join(ExportDir, TrgRoiColFname)
    
    TrgRoiCol.save_as(TrgRoiColFpath)
        
    print(f'New Target {TrgRoiCol.Modality} exported to:\n {TrgRoiColFpath}\n')
    
    return TrgRoiColFpath









"""
******************************************************************************
******************************************************************************
WRAP UP CopyLocalRoi(), ErrorCheckRoi(), ExportTrgAsr(), and
PlotContoursFromListOfRtss_v1() or PlotPixArrsFromListOfSegs_v1()

31/03: Need to include CreateDro()..
******************************************************************************
******************************************************************************
"""


def RunCopyLocalRoi_NOT_UP_TO_DATE(TestNums, LogToConsole=False, TxtToAddTrgRoiColName='', 
                    ResInterp='BlurThenLinear', PreResVariance=(1,1,1),
                    ApplyPostResBlur=False, PostResVariance=(1,1,1),
                    ForceReg=False, SelxOrSitk='Selx', Transform='affine', 
                    MaxIters='512', TxInterp='NearestNeighbor', 
                    ApplyPostTxBin=True, ApplyPostTxBlur=True, 
                    PostTxVariance=(1,1,1), ExportRoi=True, PlotResults=False, 
                    PlotAllSlices=False, ExportPlot=False, ExportLogFiles=True):
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
    
    TxtToAddTrgRoiColName : string (optional, '' by default)
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
    
    PreResVariance : tuple of floats (optional; (1, 1, 1) by default)
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
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    SelxOrSitk : string (optional; 'Selx' by default)
        Denotes which package to use for image registration and transformation.
        Acceptable values include:
            - 'Selx' for SimpleElastix
            - 'Sitk' for SimpleITK
    
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Transform. If != 'default' it must be a string 
        representation of an integer.
     
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
        tranformed labelmap image(s).
    
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
    
        SrcRoiColFpath, SrcSliceNum, SrcRoiName, TrgSliceNum, SrcDcmDir,\
        TrgDcmDir, TrgRoiColFpath, TxtToAddToTrgRoiColName, SrcLabel,\
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
        SrcRoiCol = dcmread(SrcRoiColFpath)
        SrcModality = SrcRoiCol.Modality
        
        #print(f'\n\nSrcRoiFpath = {SrcRoiFpath}')
        #print(f'SrcModality = {SrcModality}')
    
        if TrgRoiColFpath:
            TrgRoiCol = dcmread(TrgRoiColFpath)
        else:
            TrgRoiCol = None
        
        """ New ROI label for RTS StructureSetLabel / SEG Series Description 
        and RTS/SEG filename: """
        if 'RR' in TestNum or 'RD' in TestNum:
            NewRoiLabel = SrcRoiCol.StructureSetLabel + TxtToAddToTrgRoiColName
        else:
            """ 'SR' in TestNum or 'SD' in TestNum: """
            NewRoiLabel = SrcRoiCol.SeriesDescription + TxtToAddToTrgRoiColName
            
        
        """ Text to add to file names: """
        TxtToAddToFname = f'Test_{TestNum}_{NewRoiLabel}'
        
        
        """ Copy Contours/Segmentations """
        
        NewTrgRoiCol, TxParams, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyLocalRoiCol(SrcRoiColFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, 
                          TrgDcmDir, TrgRoiColFpath, TrgSliceNum, ResInterp, 
                          PreResVariance, ApplyPostResBlur, PostResVariance, 
                          ForceReg, SelxOrSitk, Transform, MaxIters, TxInterp, 
                          ApplyPostTxBin, ApplyPostTxBlur, PostTxVariance, 
                          TxtToAddToTrgRoiColName, LogToConsole, ExportLogFiles=False)
        
        
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
        
        ErrorList, Nerrors = ErrorCheckRoiCol(NewTrgRoiCol, TrgDcmDir, 
                                              LogToConsole, DictOfInputs, False)
        
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
                TxtToAddToFname += f'_ForcedReg_{Transform}'
                
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
            
            if NewTrgRoiCol.Modality == 'RTSTRUCT':
                NewTrgRoiColFpath = ExportTrgRoiCol(NewTrgRoiCol, SrcRoiColFpath,
                                                    RtsExportDir, TxtToAddToFname, 
                                                    DictOfInputs)
            else:
                NewTrgRoiColFpath = ExportTrgRoiCol(NewTrgRoiCol, SrcRoiColFpath, 
                                                    SegExportDir, TxtToAddToFname,
                                                    DictOfInputs)
        
        
        """ Plot Contours """
    
        if PlotResults:
            if ExportPlot:
                dpi = 120
            else:
                dpi = 80
    
    
    
            if TrgRoiCol:
                ListOfRois = [SrcRoiCol, TrgRoiCol, NewTrgRoiCol]
    
                ListOfDicomDirs = [SrcDcmDir, TrgDcmDir, TrgDcmDir]
    
                ListOfPlotTitles = ['Source', 'Original Target', 'New Target']
            else:
                ListOfRois = [SrcRoiCol, NewTrgRoiCol]
    
                ListOfDicomDirs = [SrcDcmDir, TrgDcmDir]
    
                ListOfPlotTitles = ['Source', 'New Target']
    
    
            #TxtToAdd = f'(RunCase = {RunCase}, \nSrcSliceNum = {SrcSliceNum}, '\
            #           + f'\nTrgSliceNum = {TrgSliceNum})'
            
            #print(f'\nListOfDicomDirs = {ListOfDicomDirs}')
            
            Times.append(time.time())
            
            if NewTrgRoiCol.Modality == 'RTSTRUCT':
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
    





"""
******************************************************************************
******************************************************************************
WRAP UP CopyXnatRoiCol(), ErrorCheckRoiCol(), ExportTrgRoiCol(), and
PlotContoursFromListOfRtss_v1() or PlotPixArrsFromListOfSegs_v1()
******************************************************************************
******************************************************************************
"""

def RunCopyXnatRoiCol(XnatSession, TestNum, LogToConsole=False, 
                      PlotResults=False, ExportPlot=False, DevOutputs=False,
                      XnatUrl=None, ProjId=None, SubjLabel=None, 
                      SrcExpLabel=None, SrcScanId=None, SrcSliceNum=None, 
                      SrcRoiColMod=None, SrcRoiColName=None, SrcRoiName=None, 
                      TrgExpLabel=None, TrgScanId=None, TrgSliceNum=None,  
                      TrgRoiColMod=None, TrgRoiColName=None, TrgRoiName=None, 
                      TxtToAddToTrgRoiColName='', 
                      PathsDict=None, XnatDownloadDir='default', 
                      ExportLogFiles=False):
    """
    Wrapper function for CopyXnatRoiCol(), ErrorCheckRoiCol(), ExportTrgRoiCol(),
    and PlotContoursFromListOfRtss_v1() or PlotPixArrsFromListOfSegs_v1().
    
    Inputs:
    ******
    
    XnatSession : requests session or None
        
    TestNum : string or None
        Either a string denoting the test to run (as defined in 
        CopyRoiTestConfig.py, e.g. 'SR1'), or None. 
        To run the algorithm on a different combination of data, either update 
        CopyRoiTestConfig.py or set TestNums = None and provide the necessary 
        optional inputs as described below.
    
    LogToConsole : boolean (optional; False by default)
        If True, intermediate results will be logged to the console during the
        running of CopyRoiCol().
    
    PlotResults : boolean (optional; False by default)
    
    ExportPlot : boolean (optional; False by default)
    
    DevOutputs : boolean (optional; False by default)
        If True various outputs will be returned.  If False, no outputs will be
        returned.
    
    XnatUrl : string (optional if TestNum != None; None by default)
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string (optional if TestNum != None; None by default)
        The project ID of interest.
    
    SubjLabel : string (optional if TestNum != None; None by default)
        The subject label of interest. 
    
    SrcExpLabel : string (optional if TestNum != None; None by default)
        The Source DICOM study / XNAT experiment label. 
    
    SrcScanId : string (optional if TestNum != None; None by default)
        The Source DICOM series label / XNAT scan ID.
    
    SrcSliceNum : integer (0-indexed)
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied.
    
    SrcRoiColMod : string
        The modality of the Source ROI Collection.
        
    SrcRoiColName : string
        The Source ROI Collection name (StructureSetLabel or SeriesDescription).
        
    SrcRoiName : string
        The Source ROIName or SegmentLabel.
    
    TrgExpLabel : string (optional if TestNum != None; None by default)
        The Target DICOM study / XNAT experiment label. 
    
    TrgScanId : string (optional if TestNum != None; None by default)
        The Target DICOM series label / XNAT scan ID.
    
    TrgSliceNum : integer (optional unless making a direct copy; 0-indexed; 
    None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to.  This only applies for direct copies, 
        hence the default value None.
    
    TrgRoiColMod : string (optional but required if TrgRoiColName != None; 
    None by default)
        The Target image assessor modality.
        
    TrgRoiColName : string (optional but required if TrgRoiColMod != None; 
    None by default)
        The Target ROI Collection name (StructureSetLabel or SeriesDescription). 
        If provided and if a direct copy is to be made, existing contours/
        segmentations for the ROI/segment will be preserved.
    
    TrgRoiName : string (optional; None by default)
        The Target ROIName or SegmentLabel.
    
    TxtToAddToTrgRoiColName : string (optional, '' by default)
        If provided the string of text will be appended to the scan assessor  
        name (StructureSetLabel or SeriesDescription), and to the filename of 
        the exported (new) Target RTS/SEG.
        
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
    ********
    
    XnatSession : requests session
    
    SrcRoiCol : Pydicom object
        The Source ROI Collection (RTS/SEG) object. Returned only if 
        DevOutputs = True.
    
    TrgRoiCol : Pydicom object or None
        The Target ROI Collection (RTS/SEG) object or None. Returned only if 
        DevOutputs = True.
    
    NewTrgRoiCol : Pydicom object
        The new Target ROI Collection (RTS/SEG) object. Returned only if 
        DevOutputs = True.
        
    Dro : Pydicom object or None
        The DICOM Registration Object (if image registration was used) or None
        (if not).
        
    PathsDict : dictionary
        Dictionary of directory and file paths. Returned only if 
        DevOutputs = True.
        
    DictOfInputs : dictionary
        Dictionary of inputs. Returned only if DevOutputs = True.
    
    ListOfInputs : list
        List of inputs. Returned only if DevOutputs = True.
    
    TimingMsgs : list of strings
        A list of the time to execute certain tasks. Returned only if 
        DevOutputs = True.
               
    Times : list of floats
        A list of time stamps at certain task executions. Returned only if 
        DevOutputs = True.
    
    
    
    
    Notes:
    *****
    
    If running one or more pre-configured tests, CopyRoiTestConfig.py must be
    copied to the current working directory.
    """
    
    import importlib
    import XnatTools
    importlib.reload(XnatTools)

    import time
    import os
    from pypref import Preferences
    from pydicom import dcmread
    #from DicomTools import IsSameModalities
    #from XnatTools import DownloadScan, DownloadAssessor
    from XnatTools import GetFnameAndIdFromName, UploadSubjAsr, SearchTxMatrix
    from GeneralTools import PrintTitle, ExportListToTxt#, ExportDictionaryToJson
    
    
    """ Ensure that inputs are valid. """
    if TestNum == None and (XnatUrl == None or ProjId == None \
                            or SubjLabel == None or SrcExpLabel == None \
                            or SrcScanId == None or SrcRoiColMod == None \
                            or SrcRoiColName == None or SrcRoiName == None \
                            or TrgExpLabel == None or TrgScanId == None):
        msg = "If a TestNum is not provided, inputs to XnatUrl, ProjId, "\
              + "SubjLabel, SrcExpLabel, SrcScanId, SrcRoiColMod,\nSrcRoiColName"\
              + "SrcRoiName, TrgExpLabel and TrgScanId must all be provided."
        
        raise Exception(msg)
    
    if (TrgRoiColMod != None and TrgRoiColName == None) \
    or (TrgRoiColMod == None and TrgRoiColName != None):
        msg = "Both a Target ROI Collection name and modality must be provided"\
              + " - not only one or the other."
        
        raise Exception(msg)
    
    
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
    
    if LogToConsole:
        print(f'\nThe pre-configured inputs for TestNum {TestNum}:')
        for key, value in TestCfg.get(TestNum).items():
            print(f'   {key} : {value}')
        
    
    if TestNum:
        if not isinstance(TestNum, str):
            msg = 'The input argument "TestNum" must be a character string.'
            print(msg)
            raise Exception(msg)
        
        XnatUrl = TestCfg.get('XnatUrl')
        ProjId = TestCfg.get('ProjId')
        SubjLabel = TestCfg.get('SubjLabel')
        SrcExpLabel = TestCfg.get(TestNum)['SrcExpLabel']
        SrcScanId = TestCfg.get(TestNum)['SrcScanId']
        SrcSliceNum = TestCfg.get(TestNum)['SrcSliceNum']
        SrcRoiColMod = TestCfg.get(TestNum)['SrcRoiColMod']
        SrcRoiColName = TestCfg.get(TestNum)['SrcRoiColName']
        SrcRoiName = TestCfg.get(TestNum)['SrcRoiName']
        TrgExpLabel = TestCfg.get(TestNum)['TrgExpLabel']
        TrgScanId = TestCfg.get(TestNum)['TrgScanId']
        TrgSliceNum = TestCfg.get(TestNum)['TrgSliceNum']
        TrgRoiColMod = TestCfg.get(TestNum)['TrgRoiColMod']
        TrgRoiColName = TestCfg.get(TestNum)['TrgRoiColName']
        TrgRoiName = TestCfg.get(TestNum)['TrgRoiName']
        TxtToAddToTrgRoiColName = TestCfg.get(TestNum)['TxtToAddToTrgRoiColName']
        
        PrintTitle(f'Running Test {TestNum}:')
    else:
        if XnatUrl == None:
            XnatUrl = TestCfg.get('XnatUrl')
            
        if ProjId == None:
            raise Exception("ProjId is required.")
        
        if SubjLabel == None:
            raise Exception("SubjLabel is required.")
        
        if SrcExpLabel == None:
            raise Exception("SrcExpLabel is required.")
        
        if SrcScanId == None:
            raise Exception("SrcScanId is required.")
        
        if SrcRoiColMod == None:
            raise Exception("SrcRoiColMod is required.")
        
        if TrgExpLabel == None:
            raise Exception("TrgExpLabel is required.")
        
        if TrgRoiColMod == None:
            TrgRoiColMod = SrcRoiColMod
            print("TrgRoiColMod has not been provided so using SrcRoiColMod:",
                  f"{SrcRoiColMod}.")
    
    
    """ The pre-defined defaults for remaining inputs: """
    Defaults = Preferences(directory=os.getcwd(), 
                           filename='CopyRoiDefaults.py')
    
    ResInterp = Defaults.get('ResInterp')
    PreResVar = Defaults.get('PreResVariance')
    ApplyPostResBlur = Defaults.get('ApplyPostResBlur')
    PostResVar = Defaults.get('PostResVariance')
    ForceReg = Defaults.get('ForceReg')
    SelxOrSitk = Defaults.get('SelxOrSitk')
    Transform = Defaults.get('Transform')
    MaxIters = Defaults.get('MaxIters')
    TxInterp = Defaults.get('TxInterp')
    ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
    ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
    PostTxVar = Defaults.get('PostTxVariance')
    ExportNewTrgRoiCol = Defaults.get('ExportNewTrgRoiCol')
    RtsExportDir = Defaults.get('RtsExportDir')
    SegExportDir = Defaults.get('SegExportDir')
    ExportNewDro = Defaults.get('ExportNewDro')
    DroExportDir = Defaults.get('DroExportDir')
    UploadNewDro = Defaults.get('UploadNewDro')
    
    print("Default input arguments read from 'CopyRoiDefaults.py':",
          f"\n   ResInterp = {ResInterp}\n   PreResVar = {PreResVar}",
          f"\n   ApplyPostResBlur = {ApplyPostResBlur}",
          f"\n   PostResVar = {PostResVar}",
          f"\n   ForceReg = {ForceReg}\n   SelxOrSitk = {SelxOrSitk}",
          f"\n   Transform = {Transform}\n   MaxIters = {MaxIters}",
          f"\n   TxInterp = {TxInterp}\n   ApplyPostTxBin = {ApplyPostTxBin}",
          f"\n   ApplyPostTxBlur = {ApplyPostTxBlur}",
          f"\n   PostTxVar = {PostTxVar}",
          f"\n   ExportNewTrgRoiCol = {ExportNewTrgRoiCol}",
          f"\n   ExportNewDro = {ExportNewDro}"
          f"\n   UploadNewDro = {UploadNewDro}\n")
    
    
    """ Check whether there exists a DRO for this set of Source and Target
    images: """
    TxMatrix = SearchTxMatrix(XnatUrl, ProjId, SubjLabel, 
                              SrcExpLabel, SrcScanId, TrgExpLabel, TrgScanId,
                              Transform, XnatSession)
    
    #print(type(TxMatrix))
    
    NewTrgRoiCol, Dro, PathsDict, XnatSession, DictOfInputs, ListOfInputs,\
    TimingMsgs\
    = CopyXnatRoiCol(XnatUrl, XnatSession, ProjId, SubjLabel, 
                     SrcExpLabel, SrcScanId, SrcSliceNum,
                     SrcRoiColMod, SrcRoiColName, SrcRoiName,  
                     TrgExpLabel, TrgScanId, TrgSliceNum, 
                     TrgRoiColMod, TrgRoiColName, TrgRoiName, 
                     TxtToAddToTrgRoiColName, 
                     PathsDict, XnatDownloadDir, LogToConsole, ExportLogFiles,
                     ResInterp, PreResVar, ApplyPostResBlur, PostResVar, 
                     ForceReg, TxMatrix, SelxOrSitk, Transform, MaxIters, 
                     TxInterp, ApplyPostTxBin, ApplyPostTxBlur, PostTxVar)
    
    
    if TestNum:
        # Add TestNum to DictOfInputs:
        DictOfInputs.update({'TestNum' : TestNum})
    
    
    """ Error check the new Target RTS/SEG """
    Times.append(time.time())
    
    TrgDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][TrgExpLabel]['scans'][TrgScanId]\
                ['resources']['DICOM']['files']['Dir']
    
    ErrorList, Nerrors = ErrorCheckRoiCol(NewTrgRoiCol, TrgDcmDir, LogToConsole,
                                          DictOfInputs, False)
    
    Times.append(time.time())
    Dtime = round(Times[-1] - Times[-2], 1)
    msg = f'Took {Dtime} s to error check the {SrcRoiColMod}.  There were '\
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
    
    
    
    """ Info required if exporting the new assessor and/or if DevOutputs = True
    """
    UseCaseThatApplies = DictOfInputs['UseCaseThatApplies']
    UseCaseToApply = DictOfInputs['UseCaseToApply']
    
    if TestNum:
        #NewTrgAsrFname = f'TestNum_{TestNum}_{SrcRoiColName}'\
        #                 + f'{TxtToAddToTrgRoiColName}'
        NewTrgRoiColFname = f'{TestNum}_{SrcRoiColName}{TxtToAddToTrgRoiColName}'
    else:
        NewTrgRoiColFname = f'{SrcRoiColName}{TxtToAddToTrgRoiColName}'
    
    if ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
        NewTrgRoiColFname += f'_ForcedReg_{Transform}'
        
    if not ForceReg and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
        #NewTrgRoiFname += f'_{ResInterp}'
        
        """ The actual interpolation used for resampling might not be
        ResInterp. """
        DiffInterpByRoi = []
        
        for key, val in DictOfInputs.items():
            if 'ResInterpUsedFor' in key and not ResInterp in key:
                DiffInterpByRoi.append(val)
        
        """ If DiffInterpByRoi is not empty, use the first item. """
        if DiffInterpByRoi:
            NewTrgRoiColFname += f'_{DiffInterpByRoi[0]}'
        else:
            if ResInterp == 'NearestNeighbor':
                NewTrgRoiColFname += '_NN'
            else:
                NewTrgRoiColFname += f'_{ResInterp}'
        
        
        
    if UseCaseToApply in ['5a', '5b']:
        if TxInterp == 'NearestNeighbor':
            NewTrgRoiColFname += '_NN'
        else:
            NewTrgRoiColFname += f'_{TxInterp}'
    
    
    
    """ Import the Source and (original) Target (if applicable) ROI Collection(s) 
    for plotting (further below). The filename of the Source ROI Collection will
    be used to generate the filename of the new Target ROI Collection. """
    
    SrcRoiColFname, SrcAsrId\
    = GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, SrcExpLabel,
                            SrcRoiColMod, SrcRoiColName)
    
    SrcRoiColFpath = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                     ['experiments'][SrcExpLabel]['assessors'][SrcAsrId]\
                     ['resources'][SrcRoiColMod]['files'][SrcRoiColFname]\
                     ['Fpath']
                           
    
    SrcRoiCol = dcmread(SrcRoiColFpath)
    
    SrcDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][SrcExpLabel]['scans'][SrcScanId]\
                ['resources']['DICOM']['files']['Dir']
    
    TrgDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][TrgExpLabel]['scans'][TrgScanId]\
                ['resources']['DICOM']['files']['Dir']
                         
    if TrgRoiColName:
        TrgRoiColFname, TrgAsrId\
        = GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, TrgExpLabel,
                                TrgRoiColMod, TrgRoiColName)
        
        TrgRoiColFpath = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                         ['experiments'][TrgExpLabel]['assessors']\
                         [TrgAsrId]['resources'][TrgRoiColMod]['files']\
                         [TrgRoiColFname]['Fpath']
        
        TrgRoiCol = dcmread(TrgRoiColFpath)
    else:
        TrgRoiCol = None
    
    
    
    
    """ Export the new Target RTS/SEG """

    if ExportNewTrgRoiCol:# and not Nerrors:
        if NewTrgRoiCol.Modality == 'RTSTRUCT':
            NewTrgRoiColFpath\
            = ExportTrgRoiCol(NewTrgRoiCol, SrcRoiColFpath, RtsExportDir, 
                              NewTrgRoiColFname, DictOfInputs)
        else:
            NewTrgRoiColFpath\
            = ExportTrgRoiCol(NewTrgRoiCol, SrcRoiColFpath, SegExportDir, 
                              NewTrgRoiColFname, DictOfInputs)
        
        NewTrgRoiColFname = os.path.split(NewTrgRoiColFpath)[1]
        
        """ Update PathsDict: 
            14/04/21: Commenting this out so that PathsDict only relates to
            data that has been downloaded from XNAT.
        keys = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
        ['experiments'][TrgExpLabel]['assessors']\
        .update({NewTrgScanAsrId : {'resources' : {TrgRoiColMod : {'files' : {NewTrgScanAsrFname : {'ScanAsrDate' : NewTrgScanAsr.ContentDate,
                                                                                                     'ScanAsrTime' : NewTrgScanAsr.ContentTime,
                                                                                                     #'AsrLabel' : , # generated by XNAT
                                                                                                     'RoiMod' : TrgRoiColMod,
                                                                                                     #'AsrId' : , # generated by XNAT
                                                                                                     'ScanAsrDir' : NewTrgScanAsrFname,
                                                                                                     'ScanAsrFpath' : NewTrgScanAsrFpath
                                                                                                     }
                                                                               }
                                                                    }
                                                   }
                                    }
               }) 
        """
        
                               
    
    
    """ Export the DRO (if applicable) """
    
    #if TxParams != None:
    #if UseCaseToApply in ['5a', '5b']:
    if Dro != None and ExportNewDro:
        from DroTools import ExportDro
        
        DroLabel = f'Study_{SrcExpLabel}_Series_{SrcScanId}_to_'\
                   + f'Study_{TrgExpLabel}_Series_{TrgScanId}'
        
        DroFname = f'{DateTime}_DRO_'
        
        if TestNum:
            DroFname += f'{TestNum}_'
        
        DroFname += DroLabel + '.dcm'
            
        DroFpath = ExportDro(Dro, DroExportDir, DroFname, DictOfInputs)
        
        """ Update PathsDict: 
            14/04/21: Commenting this out so that PathsDict only relates to
            data that has been downloaded from XNAT.
            Also, this needs to be updated to the current format (as above).
        keys = PathsDict[ProjId][SubjLabel].keys()
        
        if DroLabel in keys:
            PathsDict[ProjId][SubjLabel][DroLabel]\
            .update({'DroDate' : Dro.ContentDate,
                     'DroTime' : Dro.ContentTime,
                     'DroDir' : os.path.split(DroFpath)[0],
                     'DroFpath' : DroFpath
                     })
        else:
            PathsDict[ProjId][SubjLabel]\
            .update({DroLabel : {'DroDate' : Dro.ContentDate,
                                 'DroTime' : Dro.ContentTime,
                                 'DroDir' : os.path.split(DroFpath)[0],
                                 'DroFpath' : DroFpath
                                 }
                     })
        """
            
    
    
    """ Upload the DRO to XNAT (if applicable) """
    
    if Dro != None and UploadNewDro:
        XnatSession = UploadSubjAsr(DroFpath, XnatUrl, ProjId, SubjLabel, 
                                    ContentLabel='SRO-DRO', 
                                    XnatSession=XnatSession)
    
    
    
    
    """ Plot Contours """

    if PlotResults:
        import PlottingTools
        import importlib
        importlib.reload(PlottingTools)
            
        if TrgRoiColName:
            ListOfRois = [SrcRoiCol, TrgRoiCol, NewTrgRoiCol]

            ListOfDicomDirs = [SrcDcmDir, TrgDcmDir, TrgDcmDir]

            ListOfPlotTitles = ['Source', 'Original Target', 'New Target']
        else:
            ListOfRois = [SrcRoiCol, NewTrgRoiCol]

            ListOfDicomDirs = [SrcDcmDir, TrgDcmDir]

            ListOfPlotTitles = ['Source', 'New Target']
        
        if ExportPlot:
            #dpi = 120
            dpi = 80
        else:
            dpi = 80

        RtsPlotExportDir = Defaults.get('RtsPlotExportDir')
        SegPlotExportDir = Defaults.get('SegPlotExportDir')
        PlotAllSlices = Defaults.get('PlotAllSlices')
        
        #if TestNum:
        #    PlotFname = f'{DateTime}_TestNum_{TestNum}_{SrcRoiColName}'\
        #                + f'{TxtToAddToTrgRoiColName}'
        #    
        #else:
        #    PlotFname = f'{DateTime}_{SrcRoiColName}{TxtToAddToTrgRoiColName}'
        
        #TxtToAdd = f'(RunCase = {RunCase}, \nSrcSliceNum = {SrcSliceNum}, '\
        #           + f'\nTrgSliceNum = {TrgSliceNum})'
        
        #print(f'\nListOfDicomDirs = {ListOfDicomDirs}')
        
        NewTrgRoiColFnameNoExt = os.path.splitext(NewTrgRoiColFname)[0]
        
        PlotFname = f'{DateTime}_{NewTrgRoiColFnameNoExt}'
        PlotFname = PlotFname.replace(' ', '_') + '.jpg'
        
        Times.append(time.time())
        
        if NewTrgRoiCol.Modality == 'RTSTRUCT':
            from PlottingTools import PlotContoursFromListOfRtss_v1
            
            PlotContoursFromListOfRtss_v1(ListOfRois, ListOfDicomDirs, 
                                          ListOfPlotTitles,
                                          PlotAllSlices=PlotAllSlices,
                                          ExportPlot=ExportPlot, 
                                          ExportDir=RtsPlotExportDir,
                                          ExportFname=PlotFname,
                                          dpi=dpi, LogToConsole=False)

        else:
            from PlottingTools import PlotPixArrsFromListOfSegs_v1
            
            #print(f'SegPlotExportDir = {SegPlotExportDir}\n')
            #print(f'PlotFname = {PlotFname}\n')
            
            PlotPixArrsFromListOfSegs_v1(ListOfRois, ListOfDicomDirs, 
                                         ListOfPlotTitles,
                                         PlotAllSlices=PlotAllSlices,  
                                         ExportPlot=ExportPlot, 
                                         ExportDir=SegPlotExportDir,
                                         ExportFname=PlotFname,
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
    
    
    
    if DevOutputs:
        return XnatSession, SrcRoiCol, TrgRoiCol, NewTrgRoiCol, Dro,\
               PathsDict, DictOfInputs, ListOfInputs, TimingMsgs, Times
    else:
        return XnatSession