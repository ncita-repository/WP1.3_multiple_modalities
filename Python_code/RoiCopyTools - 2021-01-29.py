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

1a. Direct copy of the Source contour/segmentation on slice FromSliceNum to a 
    new/existing ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images are the same. 

1b. Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images are the same. 
    
2a. Direct copy of the Source contour/segmentation on slice FromSliceNum to a 
    new/existing ROI for slice ToSliceNum in Target. 
    The Source and Target images have the same FOR, IOP, IPP, and VS.

2b. Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images have the same FOR, IOP, IPP, and VS. 

3a. Direct copy of the Source contour/segmentation on slice FromSliceNum to a 
    new/existing ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images have the same FOR and IOP but have 
    different VS (they will likely also have different IPP).
    
3b. Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have the same FOR and IOP but have 
    different VS (they will likely also have different IPP).

4.  Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have the same FOR but have different IOP 
    and IPP, and possibly different VS.
   
5.  Relationship-preserving copy of the Source contour/contour on slice 
    FromSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have different FOR, likely different IPP and 
    IOP, and possibly different VS.
    
FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
IOP = ImageOrientationPatient, VS = Voxel spacings
"""




"""
******************************************************************************
******************************************************************************
WHICH USE CASE APPLIES?
******************************************************************************
******************************************************************************
"""

def WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole=False):
    """
    Determine which of 5 Use Cases applies.
    
    Inputs:
    ******
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    ToSliceNum : integer
        Slice index within the Target DICOM stack where the segmentation is to
        be copied to (counting from 0).  This only applies for Direct copies.
                        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.

              
    Outputs:
    ******
        
    UseCase : string
        String that denotes the Use Case that applies.
    """
    
    from DicomTools import GetDicomUids
    from ImageTools import GetImageAttributes
    from GeneralTools import AreListsEqualToWithinEpsilon
    
    # Get the DICOM UIDs:
    SrcStudyuid, SrcSeriesuid, SrcFORuid,\
    SrcSOPuids = GetDicomUids(DicomDir=SrcDcmDir)
    
    TrgStudyuid, TrgSeriesuid, TrgFORuid,\
    TrgSOPuids = GetDicomUids(DicomDir=TrgDcmDir)
    
    if LogToConsole:
        print(f'\nToSliceNum = {ToSliceNum}')
        print(f'len(TrgSOPuids) = {len(TrgSOPuids)}')
        
    if ToSliceNum:
        if ToSliceNum > len(TrgSOPuids) - 1:
            msg = f'The input argument "ToSliceNum" = {ToSliceNum} (zero-'\
                  + 'indexed) exceeds the number of Target DICOMs '\
                  + f'({len(TrgSOPuids)}).'
            
            raise Exception(msg)
    
    
    # Get the Image Attributes for Source and Target using SimpleITK:
    SrcSize, SrcSpacing, SrcST, SrcIPPs,\
    SrcDirs = GetImageAttributes(DicomDir=SrcDcmDir, Package='sitk')
    
    TrgSize, TrgSpacing, TrgST, TrgIPPs,\
    TrgDirs = GetImageAttributes(DicomDir=TrgDcmDir, Package='sitk')
    
    # Are the Source and Target part of the same study, series and do they have
    # the same SOP UIDs?
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
    
    # Are the Source and Target spacings, directions and origins equal to  
    # within 1e-06?
    SameSpacings = AreListsEqualToWithinEpsilon(SrcSpacing, TrgSpacing)
    SameDirs = AreListsEqualToWithinEpsilon(SrcDirs, TrgDirs)
    SameOrigins = AreListsEqualToWithinEpsilon(SrcIPPs[0], TrgIPPs[0])
    SameST = AreListsEqualToWithinEpsilon([SrcST], [TrgST])
    
    if SrcSize == TrgSize:
        SameSize = True
    else:
        SameSize = False
    
        
    # Check which case follows:
    if SameFOR:
        if SameDirs:
            if SameSpacings and SameOrigins:
                if SameSeries:
                    if isinstance(ToSliceNum, int):
                        UseCase = '1a' # Direct copy
                    else:
                        UseCase = '1b' # Relationship-preserving copy
                else:
                    if isinstance(ToSliceNum, int):
                        UseCase = '2a' # Direct copy
                    else:
                        UseCase = '2b' # Relationship-preserving copy
            else:
                if isinstance(ToSliceNum, int):
                    UseCase = '3a' # Direct copy
                else:
                    UseCase = '3b' # Relationship-preserving copy
        else:
            """ May need to break UseCase 4 into a and b for Direct and
            Relationship-preserving copies. """
            UseCase = '4'
    else:
        """ May need to break UseCase 5 into a and b for Direct and
        Relationship-preserving copies. """
        UseCase = '5'
    
        
    # Print comparisons to the console:    
    if LogToConsole:
        #CompareImageAttributes(SrcDcmDir, TrgDcmDir)
        
        print('\n\nThe Source and Target are in the/have:')
        
        if SameStudy:
            print('\n   same study')
        else:
            print('\n   different studies')
        
        if SameFOR:
            print('\n   same Frame of Reference')
        else:
            print('\n   different Frames of Reference')
            
        if SameSeries:
            print('\n   same series')
        else:
            print('\n   different series')
            
        if SameSOPs:
            print('\n   same DICOMs')
        else:
            print('\n   different DICOMs')
        
        if SameSize:
            print('\n   same image size\n',
                  f'      Source/Target: {SrcSize}')
        else:
            print('\n   different image sizes:\n',
                  f'      Source: {SrcSize} \n',
                  f'      Target: {TrgSize}')
        
        if SameSpacings:
            print('\n   same voxel sizes\n',
                  f'      Source/Target: {SrcSpacing}')
        else:
            print('\n   different voxel sizes:\n',
                  f'      Source: {SrcSpacing} \n',
                  f'      Target: {TrgSpacing}')
        
        if SameST:
            print('\n   same slice thickness\n',
                  f'      Source/Target: {SrcST}')
        else:
            print('\n   different slice thickness:\n',
                  f'      Source: {SrcST} \n',
                  f'      Target: {TrgST}')
        
        if SameOrigins:
            print('\n   same origin\n',
                  f'      Source/Target: {SrcIPPs[0]}')
        else:
            print('\n   different origins:\n',
                  f'      Source: {SrcIPPs[0]} \n',
                  f'      Target: {TrgIPPs[0]}')
            
        if SameDirs:
            print('\n   same patient orientation\n',
                  f'      Source/Target: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},\n',
                  f'                      {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},\n',
                  f'                      {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]')
        else:
            print('\n   different patient orientations:\n',
                  f'      Source: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},\n',
                  f'               {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},\n',
                  f'               {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]\n',
                  f'      Target: [{TrgDirs[0]}, {TrgDirs[1]}, {TrgDirs[2]},\n',
                  f'               {TrgDirs[3]}, {TrgDirs[4]}, {TrgDirs[5]},\n',
                  f'               {TrgDirs[6]}, {TrgDirs[7]}, {TrgDirs[8]}]')
        
        
    return UseCase






"""
******************************************************************************
******************************************************************************
COPY A CONTOUR
******************************************************************************
******************************************************************************
"""

def CopyRts(SrcRts, FromSearchString, SrcDcmDir, FromSliceNum, TrgDcmDir,
            TrgRts=None, ToSliceNum=None, Sigma=(1,1,1), ThreshLevel=0.75,
            AddText='', LogToConsole=False):
    """
    
    Inputs:
    ******
    
    SrcRts : Pydicom object
        Source RTS object.
                             
    FromSearchString : string
        All or part of the Source ROI Name of the ROI containing the contour to
        be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour to
        be copied (counting from 0).
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
        
    TrgRts : Pydicom object (optional; None by default)
        Target RTS object that the contour is to be copied to. TrgRts is only 
        non-None if a Direct copy of the contour is to be made and added to 
        existing Direct-copied contour(s).     
        
    ToSliceNum : integer (optional; None by default)
        Slice index within the Target DICOM stack where the contour is to
        be copied to (counting from 0).  This only applies for Direct copies,
        hence the default value None.
        
    Sigma : tuple (optional; (1,1,1) by default)
        The sigma values that will define the Gaussian if Gaussian image 
        blurring is used.  The tuple will define the sigma along x, y and z
        directions.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if used.  The threshold value will be ThreshLevel*M, where
        M is the maximum intensity value of the image.
    
    AddText : string (optional, '' by default)
        Sting of text to pass to ModifyRts/ModifySeg to be used when generating
        a filename for the new Target RTS/SEG.
                           
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    *******
        
    TrgRts : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour to an existing Direct-copied contour)
        Target RTS object.
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    
    from copy import deepcopy
    from GeneralTools import PrintTitle
    from GeneralTools import UniqueItems
    from DicomTools import GetRoiNum
    #from RtsTools import GetPtsInContour
    from RtsTools import GetPtsByCntForSliceNum
    from GeneralTools import ProportionOfContourInExtent
    from GeneralTools import GetPixelShiftBetweenSlices
    from GeneralTools import ShiftFrame
    from RtsTools import InitialiseRts
    from RtsTools import ModifyRts
    
    
    # Verify that an ROI exists whose name matches FromSearchString:
    #SrcRoiLabels = GetRoiLabels(SrcRts)
    FromRoiNum = GetRoiNum(Roi=SrcRts, SearchString=FromSearchString)
    
    
    #PtsToCopy = GetPtsInContour(Rts=SrcRts, DicomDir=SrcDcmDir, 
    #                            SearchString=FromSearchString, 
    #                            SliceNum=FromSliceNum)
    PtsByCnt = GetPtsByCntForSliceNum(Rts=SrcRts, DicomDir=SrcDcmDir, 
                                      SliceNum=FromSliceNum,
                                      SearchString=FromSearchString,
                                      LogToConsole=LogToConsole)
    
    """
    At present the algorithm can only copy a single contour so irrespective of
    the number of contours on FromSliceNum, only work with the first contour:
    """
    PtsToCopy = PtsByCnt[0]
    
    if LogToConsole:
        print(f'\n\nThe contour to be copied from slice {FromSliceNum} has',
              f'{len(PtsToCopy)} points.')# and CStoSliceIndsToCopy =',
              #f'{CStoSliceIndsToCopy}.')
    
    
    """
    First must determine whether the pixels that make up the mask coincide 
    with the Target image extent (14/01/2021). 
    """
    FracProp = ProportionOfContourInExtent(Points=PtsToCopy,
                                           TrgDicomDir=TrgDcmDir,
                                           LogToConsole=LogToConsole)
    
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the points in the contour to',
              'be copied lie within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the points in the contour to be copied lie '\
              + 'within the physical extent of the Target image.'
              
        raise Exception(msg)
        
        return None
    
    
    # Determine which Use Case applies:
    UseCase = WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole)
    
    if LogToConsole:
        PrintTitle(f'Case {UseCase} applies.')
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    if 'a' in UseCase:
        #from RtsTools import GetPtsInRoi
        from RtsTools import GetPtsByCntForRoi
        from RtsTools import AddCopiedPtsByCnt
        from ImageTools import GetImageAttributes
        from ImageTools import ImportImage
        from ConversionTools import Points2Indices
        from GeneralTools import ChangeZinds
        from ConversionTools import Indices2Points
        
        
    if UseCase in ['1a', '2a']:
        """
        The contour will be copied to a different slice location, so the
        z-components of the coordinates will need to be modified.  For 
        UseCase 3a this operation will be done in index space (i.e. by
        shifting pixels).
        """
        
        SrcSize, SrcSpacings, SrcST,\
        SrcIPPs, SrcDirs = GetImageAttributes(DicomDir=SrcDcmDir)
        
        # Import the source image:
        SrcIm = ImportImage(SrcDcmDir)
        
        # Convert PtsToCopy to indices:
        #IndsToCopy = Points2Indices(Points=PtsToCopy, Origin=SrcIPPs[0], 
        #                            Directions=SrcDirs, 
        #                            Spacings=SrcSpacings, 
        #                            Rounding=False)
        IndsToCopy = Points2Indices(Points=PtsToCopy, 
                                    RefIm=SrcIm, 
                                    Rounding=False)
        
        if LogToConsole:
            print(f'\n\nIndsToCopy prior to changing z-index = {IndsToCopy[:6]}...')
        
        # Modify the z-component (k) of the indices to ToSliceNum:
        IndsToCopy = ChangeZinds(Indices=IndsToCopy, NewZind=ToSliceNum)
        
        if LogToConsole:
            print(f'\n\nIndsToCopy after changing z-index = {IndsToCopy[:6]}...')
        
            print(f'\n\nPtsToCopy prior to changing z-index = {PtsToCopy[:6]}...')
        
        # Convert back to physical points:
        #PtsToCopy = Indices2Points(Indices=IndsToCopy, Origin=SrcIPPs[0], 
        #                           Directions=SrcDirs, Spacings=SrcSpacings)
        PtsToCopy = Indices2Points(Indices=IndsToCopy, RefIm=SrcIm)
        
        if LogToConsole:
            print(f'\n\nPtsToCopy after changing z-index = {PtsToCopy[:6]}...')
        
        
        if TrgRts:
            """
            Since this is a Direct copy, any existing contours in Target with ROI
            label matching FromSearchString are to be preserved.
            """
            
            # Get all points (by contour) in the Target ROI containing the contour
            # to be copied:
            #TrgPtsInRoiByCnt,\
            #TrgCStoSliceIndsInRoi = GetPtsInRoi(Rts=TrgRts, 
            #                                    DicomDir=TrgDcmDir,
            #                                    SearchString=FromSearchString)
            TrgPtsInRoiByCnt,\
            TrgCStoSliceIndsInRoi = GetPtsByCntForRoi(Rts=TrgRts, 
                                                      DicomDir=TrgDcmDir,
                                                      SearchString=FromSearchString)
            
            print('\nThe Target ROI containing the contour to be copied has',
                  f'{len(TrgPtsInRoiByCnt)} contours. \nTrgCStoSliceIndsInRoi =',
                  f'{TrgCStoSliceIndsInRoi}.')
            
            
            """
            PtsToCopy need to be added to previously Direct-copied contour points
            in TrgPtsInRoiByCnt.  The contour-sequence-to-slice-index will be
            [ToSliceNum] rather than [FromSliceNum].
            
            Note:
                Since PtsToCopy is a list of points it needs to be enclosed in []
                so that it is a list (of length 1) of a list of points so that it
                is in the form expected of PtsToAddByCnt.
            """
            
            if LogToConsole:
                print('\n\nInputs to AddCopiedPtsByCnt:')
                print(f'   len(TrgPtsInRoiByCnt) = {len(TrgPtsInRoiByCnt)}')
                [print(f'   len(TrgPtsInRoiByCnt[{i}]) = {len(TrgPtsInRoiByCnt[i])}') for i in range(len(TrgPtsInRoiByCnt))]
                print(f'   TrgCStoSliceIndsInRoi = {TrgCStoSliceIndsInRoi}')
                print(f'   len(PtsToCopy) = {len(PtsToCopy)}')
                print(f'   ToSliceNum = {ToSliceNum}')
            
            
            
            #TrgCntDataByCnt,\
            #TrgPtsByCnt,\
            #TrgCStoSliceInds = AddCopiedPtsByCnt(OrigPtsByCnt=TrgPtsInRoiByCnt, 
            #                                     OrigCStoSliceInds=TrgCStoSliceIndsInRoi, 
            #                                     PtsToAddByCnt=[PtsToCopy], 
            #                                     CStoSliceIndsToAddByCnt=[ToSliceNum],
            #                                     LogToConsole=LogToConsole)
            
            TrgCntDataByObjByFrame,\
            TrgPtsByObjByFrame,\
            TrgCStoSliceInds = AddCopiedPtsByCnt(OrigPtsByCnt=TrgPtsInRoiByCnt, 
                                                 OrigCStoSliceInds=TrgCStoSliceIndsInRoi, 
                                                 PtsToAddByCnt=[PtsToCopy], 
                                                 CStoSliceIndsToAddByCnt=[ToSliceNum],
                                                 LogToConsole=LogToConsole)
        
        else:
            """
            A Target RTS was not provided.
            """
            
            from ConversionTools import Points2ContourData
            
            TrgCStoSliceInds = deepcopy([FromSliceNum])
            
            TrgPtsByObjByFrame = [[deepcopy(PtsToCopy)]]
            
            TrgCntDataByObjByFrame = [[Points2ContourData(PtsToCopy)]]
        
    
    if UseCase in ['1b', '2b']:
        """
        TrgPtsByCnt is simply [PtsToCopy], so that it is a list (of length 1, 
        for the single contour).  Likewise for TrgCntDataByCnt.  The Target
        contour-sequence-to-slice-index will be [FromSliceNum].
        """
        
        from ConversionTools import Points2ContourData
        
        #TrgPtsByCnt = [deepcopy(PtsToCopy)]
        TrgPtsByObjByFrame = [[deepcopy(PtsToCopy)]]
        
        #TrgCntDataByCnt = [Points2ContourData(PtsToCopy)]
        TrgCntDataByObjByFrame = [[Points2ContourData(PtsToCopy)]]
        
        #TrgCStoSliceInds = deepcopy(CStoSliceIndsToCopy)
        TrgCStoSliceInds = deepcopy([FromSliceNum])
        
        #print(f'\n\nlen(PtsToCopy) = {len(PtsToCopy)}')
        #print(f'len(TrgPtsByCnt[0]) = {len(TrgPtsByCnt[0])}')
        #print(f'len(TrgCntDataByCnt[0]) = {len(TrgCntDataByCnt[0])}')
        
        
    
    if UseCase in ['3a', '3b', '4', '5']:
        import numpy as np
        from ImageTools import ImportImage
        from ImageTools import ResampleImage
        from ImageTools import GaussianBlurImage
        from ImageTools import BinaryThresholdImage
        from ImageTools import GetImageInfo
        from ConversionTools import ConvertImagePixelType
        #from RtsTools import GetMaskFromRts
        from RtsTools import GetMaskFromContoursForSliceNum
        from ConversionTools import PixArr2Image
        from ConversionTools import Image2PixArr
        #from ConversionTools import PixArr2Contours
        from ConversionTools import PixArr2PtsByContour
        from GeneralTools import ProportionOfPixArrInExtent
        from ConversionTools import NumOfListsAtDepthTwo
        
        
        # Import the 3D images:
        SrcIm = ImportImage(SrcDcmDir)
        TrgIm = ImportImage(TrgDcmDir)
        
        
        """
        ***********************************************************************
        Convert the contour to be copied to a 3D labelmap image.
        ***********************************************************************
        """    
        # Convert the contour data to be copied to a 2D mask:
        #PixArrToCopy = GetMaskFromRts(Rts=SrcRts, 
        #                              SearchString=FromSearchString, 
        #                              SliceNum=FromSliceNum, 
        #                              DicomDir=SrcDcmDir, 
        #                              RefImage=SrcIm)
        PixArrToCopy = GetMaskFromContoursForSliceNum(Rts=SrcRts,  
                                                      SliceNum=FromSliceNum, 
                                                      DicomDir=SrcDcmDir,
                                                      RefImage=SrcIm,
                                                      SearchString=FromSearchString,
                                                      LogToConsole=LogToConsole)
        
        if LogToConsole:
            print(f'\nPixArrToCopy.shape = {PixArrToCopy.shape}')
        
        
        """
        Determine whether the pixels that make up the mask coincide 
        with the Target image extent (14/01/2021). 
        """
        FracProp = ProportionOfPixArrInExtent(PixArr=PixArrToCopy, 
                                              FrameToSliceInds=[FromSliceNum],
                                              SrcDicomDir=SrcDcmDir,
                                              TrgDicomDir=TrgDcmDir)
        
        if LogToConsole:
            print(f'\n{round(FracProp*100, 1)}% of the voxels in the',
                  'segmentation to be copied lie within the physical extent',
                  'of the Target image.')
        
        if FracProp == 0:
            msg = 'None of the voxels in the segmentation to be copied lie '\
                  + 'within the physical extent of the Target image.'
                  
            raise Exception(msg)
            
            return None
        
        
        # Convert PixArrToCopy to a labelmap image:
        LabmapImToCopy = PixArr2Image(PixArr=PixArrToCopy,
                                      FrameToSliceInds=[FromSliceNum],
                                      RefIm=SrcIm)
        
        if LogToConsole:
            print(f'\nLabmapImToCopy.GetSize() = {LabmapImToCopy.GetSize()}')
        
        
        
        
    
    if UseCase in ['3a', '3b', '4']:
        
        """ First try resampling the image using a NearestNeighbor 
        interpolation. """
        interp = 'NearestNeighbor'
        
        if LogToConsole:
            print('\nUsing', interp, 'interpolation...')
        
        # Resample LabmapImToCopy to the Target image's grid:
        ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                          RefImage=TrgIm, 
                                          Interpolation=interp)
        
        #print(f'\nResLabmapImToCopy.GetSize() = {ResLabmapImToCopy.GetSize()}')
        
        if LogToConsole:
            print('\nAfter resampling:')
            
        
        # PixID required for if statement below.
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
            
        """ 
        If resampling to a smaller grid size aliasing can occur, leading to
        the appearent disappearance of segmentations (i.e. F2Sinds = []). 
        
        Try suggestions from Ziv Yaniv here:
            https://github.com/SimpleITK/SimpleITK/issues/1277
                
        1. Use the sitkLabelGaussian interpolator.
        
        2. Gaussian blur the segmentation, then resample using a linear 
        interpolator, and then threshold so that values in [a<1<b] are mapped
        to 1. You will need to select appropriate values for a, b.
        
        3. Compute the distance map of the segmentation, then resample using a
        linear interpolator, and then threshold the absolute value of the 
        distance map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you 
        need to select an appropriate dist_threshold.

        Tried:
        
        1. Resulting Gaussian interpolated resampled image had zeros only 
        (the image was converted to float prior to resampling).
        
        2. Works for UseCase 3B-i_SEG with FromSliceNum = 26 with binary
        threshold value of 0.04.
        
        
        Comment from Bradley Lowecamp:
            The implementation of the sitkLabelGaussian interpolator in the 
            ResampleImageFilter sets the sigma to the pixel spacing of the 
            input image to the resample filter. This is why it does not do a 
            better job in reducing aliasing.
        """
        #Workaround = 1
        Workaround = 2
        #Workaround = 0 # i.e. no workaround
        
        
        if not F2Sinds and Workaround:
            if LogToConsole:
                print('\nThe resampled image has no non-zero frames (e.g. due',
                      f'to aliasing).  \nTry workaround #{Workaround}.')
            
            # Convert LabmapImToCopy from 32-bit unsigned integer to a float:
            LabmapImToCopy = ConvertImagePixelType(Image=LabmapImToCopy, 
                                                   NewPixelType='sitkFloat32')
            
            if LogToConsole:
                print('\nAfter converting to float:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(LabmapImToCopy, LogToConsole)
        
            if Workaround == 1:
                if LogToConsole:
                    print(f'\nTrying workaround #{Workaround}...')
                
                interp = 'Gaussian' # <-- get zeros only (image was converted 
                # to float prior to resampling)
                
                if LogToConsole:
                    print('\nUsing', interp, 'interpolation...')
                    
                # Resample LabmapImToCopy to the Target image's grid:
                ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                                  RefImage=TrgIm, 
                                                  Interpolation=interp)
                
                if LogToConsole:
                    print('\nAfter resampling:')
                
                
                # PixID required for if statement below.
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
                    
                    
                """ Will need to convert back to unsigned int and may need to
                do some binary thresholding, but only proceed if there are
                non-zero frames in the resampled image. """
                
                if not F2Sinds:
                    Workaround = 2 # try Workaround #2
                    
                    
            if Workaround == 2:
                if LogToConsole:
                    print(f'\nTrying workaround #{Workaround}...')
                      
                SrcSpacings = SrcIm.GetSpacing()
                TrgSpacings = TrgIm.GetSpacing()
                
                SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, TrgSpacings))
                
                """ 
                The FWHM of the Gaussian is:
                    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
                    
                Due to the Nyquist-Shannon sampling theorem, the FWHM should be
                at least 2x the voxel spacings of the Target image. For 
                symmetry considerations, it makes sense for one Source voxel 
                width to be blurred to three Target voxels widths.  :
                    
                    sigma = (3*TrgSpacings)/2.355
                    
                """
                
                if LogToConsole:
                    print(f'\nSpacingsRatio = {SpacingsRatio}')
                      
                # Gaussian blur the image:
                #LabmapImToCopy = GaussianBlurImage(LabmapImToCopy)
                LabmapImToCopy = GaussianBlurImage(Im=LabmapImToCopy,
                                                   Variance=Sigma)
                
                # Use the RecursiveGaussian image filter:
                #LabmapImToCopy = RecursiveGaussianBlurImage(LabmapImToCopy, 
                #                                            Sigma=3,
                #                                            Direction=2)
                
                if LogToConsole:
                    print('\nAfter applying Gaussian blur:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(LabmapImToCopy, LogToConsole)
                
                #if not F2Sinds:
                #    Workaround = 3 # try Workaround #3
                    
                # Linearly resample the blurred image:
                interp = 'Linear'
                
                if LogToConsole:
                    print('\nUsing', interp, 'interpolation...')
                    print(f'   LabmapImToCopy.GetSize() = {LabmapImToCopy.GetSize()}')
                    print(f'   LabmapImToCopy.GetSpacing() = {LabmapImToCopy.GetSpacing()}')
                    print(f'   TrgIm.GetSize() = {TrgIm.GetSize()}')
                    print(f'   TrgIm.GetSpacing() = {TrgIm.GetSpacing()}')
                    
                # Resample LabmapImToCopy to the Target image's grid:
                ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                                  RefImage=TrgIm, 
                                                  Interpolation=interp)
                
                if LogToConsole:
                    print('\nAfter resampling:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
                
                
                """ Prior to converting back to unsigned integer, must first
                binary threshold the image. """
                # Thresh values with results for Case3B-i_SEG and 
                # FromSliceNum = 26 (max value in Pix array prior to binary
                # thresholding = 0.0468482):
                #Thresh = 0.5 # no frames
                #Thresh = 0.005 # two frames (slices 10 & 11)
                #Thresh = 0.04 # one frame (slice 11)
                
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                if LogToConsole:
                    print('\nBinary thresholding at',
                      f'{ThreshLevel}*{max(UniqueVals)} = {ThreshValue}...')
                
                #ResLabmapImToCopy = BinaryThresholdImage(ResLabmapImToCopy, 
                #                                         Thresh)
                ResLabmapImToCopy = BinaryThresholdImage(Im=ResLabmapImToCopy, 
                                                         Thresh=ThreshValue)
                
                if LogToConsole:
                    print(f'\nAfter binary thresholding:')
                
                
                # PixID required for if statement below.
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
            
                
                
            #if Workaround == 3:
            #    print(f'\nTrying workaround #{Workaround}...')
            #          
            #    print('\nStill need to write Workaround 3...')
        
        
        
        """ Ensure that ResLabmapImToCopy is a 32-bit unsigned integer. """
        if not PixID == 5:
            if LogToConsole:
                print(f'\nPixelID = {PixID} ({PixIDTypeAsStr})). Converting',
                      f'to unsigned 32-bit integer (sitkUInt32)..')
            
            # Convert ResLabmapImToCopy from float to 32-bit unsigned integer:
            ResLabmapImToCopy = ConvertImagePixelType(Image=ResLabmapImToCopy, 
                                                      NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print('\nAfter converting to 32-bit unsigned int:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
        
        
        
        # Convert ResLabmapImToCopy to a pixel array:
        ResPixArrToCopy,\
        TrgCStoSliceInds = Image2PixArr(LabmapIm=ResLabmapImToCopy)
        
        #print(f'\nResPixArrToCopy.shape = {ResPixArrToCopy.shape}')
        #
        #unique = UniqueItems(Items=ResPixArrToCopy, IgnoreZero=False)
        #print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy')
        #
        #F = ResPixArrToCopy.shape[0]
        #
        #print(f'\nThe segmentation (from the contour) on slice {FromSliceNum}',
        #      f'has been resampled to {F} frames (contours).\nTrgCStoSliceInds',
        #      f'= {TrgCStoSliceInds}.')
        
        if LogToConsole:
            print('\nAfter converting to a pixel array:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
            
            for i in range(ResPixArrToCopy.shape[0]):
                maxval = np.amax(ResPixArrToCopy[i])
                minval = np.amin(ResPixArrToCopy[i])
                
                print(f'Frame {i} has max = {maxval}, min = {minval}')
        
    
        
    if UseCase == '3a':
        """
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled labelmap-image-to-copy).  So if after resampling there is 
        more than one frame, the frames will be averaged.
        """
        
        from GeneralTools import MeanPixArr
        
        # F required to determine whether ResPixArrToCopy will be averaged. 
        F = ResPixArrToCopy.shape[0]
        
        
        if LogToConsole:
            unique = UniqueItems(Items=ResPixArrToCopy, IgnoreZero=False)
            
            print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy')
            print('\nThe segmentation from the contour on FromSliceNum =',
                  f'{FromSliceNum} has been resampled to {F} frames/contours:',
                  f'{TrgCStoSliceInds}')
        
        if F > 1: 
            if LogToConsole:
                print(f'\nResPixArrToCopy has {F} frames so the pixel arrays',
                      'will be averaged...')
            
            # Take pixel-by-pixel mean across all frames in ResPixArrToCopy:
            ResPixArrToCopy = MeanPixArr(ResPixArrToCopy)
            
            if LogToConsole:
                print(f'\nResPixArrToCopy.shape = {ResPixArrToCopy.shape}')
            
                F = ResPixArrToCopy.shape[0]
            
                print(f'\nResPixArrToCopy now has {F} frames')
        else:
            if LogToConsole:
                print(f'\nResPixArrToCopy has F frames')
        
        
        # Get pixel shift between FromSliceNum in SrcIm and ToSliceNum in 
        # TrgIm in the Target image domain:
        PixShift = GetPixelShiftBetweenSlices(Image0=SrcIm, 
                                              SliceNum0=FromSliceNum, 
                                              Image1=TrgIm,
                                              SliceNum1=ToSliceNum,
                                              RefImage=TrgIm)
        
        if LogToConsole:
            print(f'\nPixShift = {PixShift}')
        
        # Shift the the in-plane elements in ResPixArrToCopy:
        ResPixArrToCopy = ShiftFrame(Frame=ResPixArrToCopy, PixShift=PixShift)
        
        
        if LogToConsole:
            unique = UniqueItems(Items=ResPixArrToCopy, IgnoreZero=False)
            
            print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy',
                  'after shifting the frame.')
        
        # Modify TrgCStoSliceInds to reflect the slice to which ResPixArrToCopy
        # relates to:
        TrgCStoSliceInds = deepcopy([ToSliceNum])
        
        
            
    #if '3' in UseCase: # 18/01/2021
    #if UseCase in ['3', '4']: # 18/01/2021
    if UseCase in ['3a', '3b', '4']: # 20/01/2021 <-- need to check this is ok
        # Convert ResPixArrToCopy to a list of ContourData by contour and a
        # list of points by contour: 
        #TrgCntDataByCnt, TrgPtsByCnt = PixArr2Contours(PixArr=ResPixArrToCopy)
        #TrgPtsByCnt,\
        #TrgCntDataByCnt = PixArr2PtsByContour(PixArr=ResPixArrToCopy,
        #                                      FrameToSliceInds=TrgCStoSliceInds,
        #                                      DicomDir=TrgDcmDir)
        TrgPtsByObjByFrame,\
        TrgCntDataByObjByFrame = PixArr2PtsByContour(PixArr=ResPixArrToCopy,
                                                     FrameToSliceInds=TrgCStoSliceInds,
                                                     DicomDir=TrgDcmDir)
        
        """
        03/12/20: The z-coordinates are suspiciously all the same.  Should look
        into this.
        """
        
        if LogToConsole:
            print(f'\nAfter converting ResPixArrToCopy to ContourData and',
                  f'points TrgCStoSliceInds = {TrgCStoSliceInds}.')
        
            #print(f'\nlen(TrgCntDataByCnt) = {len(TrgCntDataByCnt)}')
            #print(f'\nlen(TrgPtsByCnt) = {len(TrgPtsByCnt)}')
            print(f'\nThere are {len(TrgCntDataByObjByFrame)} frames in',
                  'TrgCntDataByObjByFrame.')
            Ncontours = NumOfListsAtDepthTwo(TrgCntDataByObjByFrame)
            print(f'There are {Ncontours} contours in TrgCntDataByObjByFrame.')
            if len(TrgCntDataByObjByFrame) != Ncontours:
                for f in range(len(TrgCntDataByObjByFrame)):
                    print(f'   Frame {f} has {len(TrgCntDataByObjByFrame[f])}',
                          'contours.')
            print(f'\nThere are {len(TrgPtsByObjByFrame)} frames in',
                  'TrgPtsByObjByFrame.')
            Ncontours = NumOfListsAtDepthTwo(TrgPtsByObjByFrame)
            print(f'There are {Ncontours} contours in TrgPtsByObjByFrame.')
            if len(TrgPtsByObjByFrame) != Ncontours:
                for f in range(len(TrgPtsByObjByFrame)):
                    print(f'   Frame {f} has {len(TrgPtsByObjByFrame[f])}',
                          'contours.')
        
        
    
    if UseCase == '3a':
        if TrgRts:
            """
            Since this is a Direct copy, any existing contours in Target with 
            ROI label matching FromSearchString are to be preserved.
            """
            
            # Get all points (by contour) in the Target ROI containing the 
            # contour to be copied:
            #TrgPtsInRoiByCnt,\
            #TrgCStoSliceIndsInRoi = GetPtsInRoi(Rts=TrgRts, 
            #                                    DicomDir=TrgDcmDir,
            #                                    SearchString=FromSearchString)
            TrgPtsInRoiByCnt,\
            TrgCStoSliceIndsInRoi = GetPtsByCntForRoi(Rts=TrgRts, 
                                                      DicomDir=TrgDcmDir,
                                                      SearchString=FromSearchString)
            
            if LogToConsole:
                print('\nThe Target ROI containing the contour to be copied',
                      f'has {len(TrgPtsInRoiByCnt)} contours.',
                      f'\nTrgCStoSliceIndsInRoi = {TrgCStoSliceIndsInRoi}.')
            
            
            """
            The contour points that arose from ResPixArrToCopy, 
            TrgPtsByObjByFrame (previously TrgPtsByCnt), need to be added to 
            previously Direct-copied contour points TrgPtsInRoiByCnt.
            
            NOTE:
                Need to check behaviour of AddCopiedPtsByCnt following change
                to format of contour data from TrgPtsByCnt to 
                TrgPtsByObjByFrame.
            """
            
            #TrgCntDataByCnt,\
            #TrgPtsByCnt,\
            #TrgCStoSliceInds = AddCopiedPtsByCnt(OrigPtsByCnt=TrgPtsInRoiByCnt, 
            #                                     OrigCStoSliceInds=TrgCStoSliceIndsInRoi, 
            #                                     PtsToAddByCnt=TrgPtsByCnt, 
            #                                     CStoSliceIndsToAddByCnt=TrgCStoSliceInds,
            #                                     LogToConsole=LogToConsole)
            TrgCntDataByCnt,\
            TrgPtsByCnt,\
            TrgCStoSliceInds = AddCopiedPtsByCnt(OrigPtsByCnt=TrgPtsInRoiByCnt, 
                                                 OrigCStoSliceInds=TrgCStoSliceIndsInRoi, 
                                                 PtsToAddByCnt=TrgPtsByObjByFrame[0], 
                                                 CStoSliceIndsToAddByCnt=TrgCStoSliceInds,
                                                 LogToConsole=LogToConsole)
            """ 19/01/21:  Need to verify that indexing TrgPtsByObjByFrame with
            [0] will result in behaviour previously obtained with TrgPtsByCnt.
            """
            
        #else:
            """
            A Target RTS was not provided.
            """
            
            #TrgPtsByCnt = deepcopy()
            #TrgCntDataByCnt = deepcopy()
            #TrgCStoSliceInds = deepcopy()
            
        
    #if UseCase == '4':
    
    
    if UseCase == '5':
        from ImageTools import RegisterImages
        from ImageTools import TransformImage
        from ImageTools import BinaryThresholdImage
    
        # Register the images using an affine transformation:
        RegSrcIm, RegImFilt = RegisterImages(FixIm=TrgIm, MovIm=SrcIm, 
                                             Tx='affine',
                                             LogToConsole=LogToConsole)
        
        # Transform LabmapImToCopy using RegImFilt:
        TxLabmapImToCopy = TransformImage(Im=LabmapImToCopy, 
                                          RegImFilt=RegImFilt,
                                          Interpolation='Nearestneighbor')
        
        """ This isn't required: """
        if False:
            # Store the pre-transformed labelmap image (LabmapImToCopy) and the 
            # pre-binary thresholded labelmap image (TxLabmapImToCopy):
            ListOfLabmapIms = [LabmapImToCopy, TxLabmapImToCopy]
            
            # Convert LabmapImToCopy to a pixel array:
            PixArrToCopy,\
            F2Sinds = Image2PixArr(LabmapIm=LabmapImToCopy)
            
            # Convert TxLabmapImToCopy to a pixel array:
            TxPixArrPreBinary,\
            F2SindsPreBinary = Image2PixArr(LabmapIm=TxLabmapImToCopy)
            
            unique = UniqueItems(Items=TxPixArrPreBinary, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArrToCopy',
                  'prior to binary thresholding TxLabmapImToCopy')
        """ end """
        
        
        
        """
        21/12/20:  
            Following the use of a nearestneighbor interpolator when 
            transforming LabmapImToCopy, i.e.:
            
            TxMap[0]["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]
            
            it's no longer necessary to perform the binary thresholding step.
            
        # Binary threshold TxLabmapImToCopy:
        TxLabmapImToCopy = BinaryThresholdImage(Im=TxLabmapImToCopy,
                                                Thresh=Thresh)
        """
        
        # Convert TxLabmapImToCopy to a pixel array:
        TxPixArrToCopy,\
        TrgCStoSliceInds = Image2PixArr(LabmapIm=TxLabmapImToCopy)
        
        if LogToConsole:
            unique = UniqueItems(Items=TxPixArrToCopy, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArrToCopy',
                  'after binary thresholding TxLabmapImToCopy')
            
            F = TxPixArrToCopy.shape[0]
            
            print('\nThe segmentation (from the contour) on slice',
                  f'{FromSliceNum} has been registered to {F} frames',
                  f'(contours): {TrgCStoSliceInds}.')
        
        
        if False:
            # Append the post-binary thresholded labelmap image (TxLabmapImToCopy):
            ListOfLabmapIms.append(TxLabmapImToCopy)
            
            # Create lists of the pre-transformed, pre-binary thresholded and 
            # post-binary thresholded pixel arrays and frame-to-slice indices:
            ListOfPixArrs = [PixArrToCopy, TxPixArrPreBinary, TxPixArrToCopy]
            ListOfFrameToSliceInds = [F2Sinds, F2SindsPreBinary, TrgCStoSliceInds]
            ListOfPlotTitles = ['PixArrToCopy', 'TxPixArrToCopyPreBinary', 
                                'TxPixArrToCopyPostBinary']
        
            return ListOfLabmapIms, ListOfPixArrs, ListOfFrameToSliceInds,\
                   ListOfPlotTitles
        
        
        if False:
            import PlottingTools
            importlib.reload(PlottingTools)
            from PlottingTools import Plot_PixelArrays
            
            Plot_PixelArrays(ListOfPixArrs=ListOfPixArrs, 
                             ListOfFrameToSliceInds=ListOfFrameToSliceInds, 
                             ListOfPlotTitles=ListOfPlotTitles)
        
        
            
    
        #return LabmapImToCopy, RegImFilt
        
        
        # Convert TxPixArrToCopy to a list of ContourData by contour and a
        # list of points by contour: 
        #TrgCntDataByCnt, TrgPtsByCnt = PixArr2Contours(PixArr=TxPixArrToCopy)
        #TrgPtsByCnt,\
        #TrgCntDataByCnt = PixArr2PtsByContour(PixArr=TxPixArrToCopy,
        #                                      FrameToSliceInds=TrgCStoSliceInds,
        #                                      DicomDir=TrgDcmDir)
        TrgPtsByObjByFrame,\
        TrgCntDataByObjByFrame = PixArr2PtsByContour(PixArr=TxPixArrToCopy,
                                                     FrameToSliceInds=TrgCStoSliceInds,
                                                     DicomDir=TrgDcmDir)
        
        if LogToConsole:
            print(f'\nThere are {len(TrgCntDataByObjByFrame)} frames in',
                  'TrgCntDataByObjByFrame.')
            Ncontours = NumOfListsAtDepthTwo(TrgCntDataByObjByFrame)
            print(f'There are {Ncontours} contours in TrgCntDataByObjByFrame.')
            if len(TrgCntDataByObjByFrame) != Ncontours:
                for f in range(len(TrgCntDataByObjByFrame)):
                    print(f'   Frame {f} has {len(TrgCntDataByObjByFrame[f])}',
                          'contours.')
            print(f'\nThere are {len(TrgPtsByObjByFrame)} frames in',
                  'TrgPtsByObjByFrame.')
            Ncontours = NumOfListsAtDepthTwo(TrgPtsByObjByFrame)
            print(f'There are {Ncontours} contours in TrgPtsByObjByFrame.')
            if len(TrgPtsByObjByFrame) != Ncontours:
                for f in range(len(TrgPtsByObjByFrame)):
                    print(f'   Frame {f} has {len(TrgPtsByObjByFrame[f])}',
                          'contours.')
    
    
    #if 'a' in UseCase:
    #    """ 
    #    For UseCases 1a/b, 2a/b or 3a, there will always be a single 
    #    ContourImageSequence, and the single index will be slice ToSliceNum. 
    #    Over-write previous TrgCStoSliceInds.
    #    """
    #    
    #    TrgCStoSliceInds = [ToSliceNum]
    #    
    #    print(f'\nTrgCStoSliceInds over-written to {TrgCStoSliceInds}')
    
    
    """ 
    If UseCase is 1a/2a/3a (i.e. a Direct copy is to be made) and a RTS has 
    been provided for Target this is not a clean copy, i.e. the copied contour
    is to be added to an existing ROI containing previously copied contour(s). 
    If a Relationship-preserving copy is to be made or if a Target RTS has not 
    been provided, the Source RTS will be used as a template.
    """
    #if not 'a' in UseCase or not TrgRts:
    #    TrgRts = InitialiseRts(RtsTemplate=SrcRts, RoiNum=FromRoiNum,
    #                           CStoSliceInds=TrgCStoSliceInds, 
    #                           DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    if LogToConsole:
        print('\n\nInputs to InitialiseRts:')
        print(f'   FromRoiNum = {FromRoiNum}')
        print(f'   TrgCStoSliceInds = {TrgCStoSliceInds}')
    
    
    if TrgRts:
        """ Use TrgRts to initialise the Target RTS. """
        #TrgRts = InitialiseRts(RtsTemplate=TrgRts, RoiNum=FromRoiNum,
        #                       CStoSliceInds=TrgCStoSliceInds,
        #                       DicomDir=TrgDcmDir, NamePrefix=AddText,
        #                       LogToConsole=LogToConsole)
        TrgRts = InitialiseRts(RtsTemplate=TrgRts, RoiNum=FromRoiNum,
                               CStoSliceInds=TrgCStoSliceInds,
                               CntDataByObjByFrame=TrgCntDataByObjByFrame,
                               DicomDir=TrgDcmDir, NamePrefix=AddText,
                               LogToConsole=LogToConsole)
    else:
        """ Use SrcRts to initialise the Target RTS. """
        #TrgRts = InitialiseRts(RtsTemplate=SrcRts, RoiNum=FromRoiNum,
        #                       CStoSliceInds=TrgCStoSliceInds, 
        #                       DicomDir=TrgDcmDir, NamePrefix=AddText,
        #                       LogToConsole=LogToConsole)
        TrgRts = InitialiseRts(RtsTemplate=SrcRts, RoiNum=FromRoiNum,
                               CStoSliceInds=TrgCStoSliceInds,
                               CntDataByObjByFrame=TrgCntDataByObjByFrame,
                               DicomDir=TrgDcmDir, NamePrefix=AddText,
                               LogToConsole=LogToConsole)
    
    if LogToConsole:
        print('\n\nInputs to ModifyRts:')
        #print(f'   len(TrgCntDataByCnt) = {len(TrgCntDataByCnt)}')
        #[print(f'   len(TrgCntDataByCnt[{i}]) = {len(TrgCntDataByCnt[i])}') for i in range(len(TrgCntDataByCnt))]
        #print(f'   len(TrgPtsByCnt) = {len(TrgPtsByCnt)}')
        #[print(f'   len(TrgPtsByCnt[{i}]) = {len(TrgPtsByCnt[i])}') for i in range(len(TrgPtsByCnt))]
        print(f'   len(TrgCntDataByObjByFrame) = {len(TrgCntDataByObjByFrame)}')
        [print(f'   len(TrgCntDataByObjByFrame[{i}]) = {len(TrgCntDataByObjByFrame[i])}') for i in range(len(TrgCntDataByObjByFrame))]
        print(f'   len(TrgPtsByObjByFrame) = {len(TrgPtsByObjByFrame)}')
        [print(f'   len(TrgPtsByObjByFrame[{i}]) = {len(TrgPtsByObjByFrame[i])}') for i in range(len(TrgPtsByObjByFrame))]
        print(f'   TrgCStoSliceInds = {TrgCStoSliceInds}')
    
    # Modify the tags in TrgRts:
    #TrgRts = ModifyRts(Rts=TrgRts, CntDataByCnt=TrgCntDataByCnt, 
    #                   PtsByCnt=TrgPtsByCnt, CStoSliceInds=TrgCStoSliceInds,
    #                   DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    TrgRts = ModifyRts(Rts=TrgRts, CntDataByObjByFrame=TrgCntDataByObjByFrame, 
                       PtsByObjByFrame=TrgPtsByObjByFrame, 
                       CStoSliceInds=TrgCStoSliceInds,
                       DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    
    return TrgRts






"""
******************************************************************************
******************************************************************************
COPY A SEGMENTATION
******************************************************************************
******************************************************************************
"""



def CopySeg(SrcSeg, FromSearchString, SrcDcmDir, FromSliceNum, TrgDcmDir,
            TrgSeg=None, ToSliceNum=None, Sigma=(1,1,1), ThreshLevel=0.75,
            AddText='', LogToConsole=False):
                                
    """
    Inputs:
    ******
    
    SrcSeg : Pydicom object
        Source DICOM SEG object.
                             
    FromSearchString : string
        All or part of the Source Segment Label of the segment containing the 
        segmentation to be copied.
                     
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                 
    TrgSeg : Pydicom object (optional; None by default)
        Target SEG object that the segmentation is to be copied to. TrgSeg is  
        only non-None if a Direct copy of the segmentation is to be made and 
        added to existing Direct-copied segmentation(s).
               
    ToSliceNum : integer (optional; None by default)
        Slice index within the Target DICOM stack where the segmentation is to
        be copied to (counting from 0).  This only applies for Direct copies,
        hence the default value None.
        
    Sigma : tuple (optional; (1,1,1) by default)
        The sigma values that will define the Gaussian if Gaussian image 
        blurring is used.  The tuple will define the sigma along x, y and z
        directions.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if used.  The threshold value will be ThreshLevel*M, where
        M is the maximum intensity value of the image.
    
    AddText : string (optional; '' by default)
        Sting of text to pass to ModifyRts/ModifySeg to be used when generating
        a filename for the new Target RTS/SEG.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    *******
        
    TrgSeg : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy of a segmentation to an existing Direct-copied 
        contour) Target SEG object.
    """
    
    import importlib
    
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    import ImageTools
    importlib.reload(ImageTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ConversionTools
    importlib.reload(ConversionTools)
    
    from copy import deepcopy
    import numpy as np
    from DicomTools import GetRoiNum
    from SegTools import GetFrameFromPixArr
    from GeneralTools import PrintTitle
    from GeneralTools import UniqueItems
    from ImageTools import ImportImage
    from SegTools import InitialiseSeg
    from SegTools import ModifySeg
    from ImageTools import GetImageInfo
    from GeneralTools import ProportionOfPixArrInExtent
    
    
    # Verify that a segment exists whose label matches FromSearchString:
    #SrcSearchStrings = GetRoiLabels(SrcSeg)
    FromSegNum = GetRoiNum(Roi=SrcSeg, SearchString=FromSearchString)
    
    if LogToConsole:
        print(f'\n\n\n***FromSegNum = {FromSegNum}')
    
    #ToSegNum = GetRoiNum(Roi=TrgSeg, SearchString=FromSearchString) # 04/12
        
        
    # Print the list of PFFGStoSliceInds:
    """This isn't required but useful."""
    from DicomTools import GetDicomSOPuids
    from SegTools import GetPFFGStoSliceIndsBySeg
    
    SegSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    
    SegPFFGStoSliceIndsBySeg = GetPFFGStoSliceIndsBySeg(Seg=SrcSeg, 
                                                        SOPuids=SegSOPuids)
    
    SegPFFGStoSliceIndsInSeg = SegPFFGStoSliceIndsBySeg[FromSegNum]
    
    if LogToConsole:
        print(f'\nThe Source SEG matching "{FromSearchString}" has'
              f'segmentations on slices {SegPFFGStoSliceIndsInSeg}.')
    """End of not-required code."""
    
    # The pixel array that only contains the frame to be copied:
    PixArrToCopy = GetFrameFromPixArr(Seg=SrcSeg, DicomDir=SrcDcmDir, 
                                      SearchString=FromSearchString, 
                                      SliceNum=FromSliceNum,
                                      LogToConsole=LogToConsole)
    
    """
    First must determine whether the pixels that make up the mask coincide with
    the Target image extent (14/01/2021). 
    """
    FracProp = ProportionOfPixArrInExtent(PixArr=PixArrToCopy, 
                                          FrameToSliceInds=[FromSliceNum],
                                          SrcDicomDir=SrcDcmDir,
                                          TrgDicomDir=TrgDcmDir)
    
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the voxels in the segmentation',
              'to be copied lie within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the voxels in the segmentation to be copied lie within'\
              + ' the physical extent of the Target image.'
              
        raise Exception(msg)
        
        return None
        
    
    # Determine which Use Case applies:
    UseCase = WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole)
    
    if LogToConsole:
        PrintTitle(f'Case {UseCase} applies.')
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    #if 'a' in UseCase:
    #    """ 
    #    For UseCases 1a/b, 2a/b or 3a, there will always be a single 
    #    PerFrameFunctionalGroupsSequence. 
    #    """
    #    
    #    TrgPFFGStoSliceInds = [FromSliceNum] # moved to end
        
    
    if 'a' in UseCase:
        """ 
        The segmentation will be copied to a different slice location, so the
        pixel shift between the indices at FromSliceNum and ToSliceNum needs to
        be determined.  This is common for all Direct copy sub-use cases.
        """
        
        from GeneralTools import GetPixelShiftBetweenSlices
        from GeneralTools import ShiftFrame
        from SegTools import GetPixArrInSeg
        from SegTools import AddCopiedPixArr
        
        if LogToConsole:
            print(f'\nCopying slice {FromSliceNum} from Source SEG to slice',
                  f'{ToSliceNum} in Target SEG..')
        
        # Import the Source 3D image:
        SrcIm = ImportImage(SrcDcmDir)
        TrgIm = ImportImage(TrgDcmDir)
        
        # Get pixel shift between FromSliceNum in SrcIm and ToSliceNum in
        # TrgIm in the Target image domain:
        PixShift = GetPixelShiftBetweenSlices(Image0=SrcIm, 
                                              SliceNum0=FromSliceNum, 
                                              Image1=TrgIm,
                                              SliceNum1=ToSliceNum,
                                              RefImage=TrgIm)
        
        if LogToConsole:
            print(f'\nPixShift = {PixShift}')
        
        
    if UseCase in ['1a', '2a']:
        # Shift the the in-plane elements in PixArrToCopy:
        PixArrToCopy = ShiftFrame(Frame=PixArrToCopy, 
                                  PixShift=PixShift)
        
        if TrgSeg:
            """
            Since this is a Direct copy, any existing segmentations in Target 
            with segment label matching FromSearchString are to be preserved.
            """
            
            # Get all frames in the Target segment containing the segmentation 
            # to be copied:
            TrgPixArrInSeg,\
            TrgPFFGStoSliceIndsInSeg = GetPixArrInSeg(Seg=TrgSeg, 
                                                      DicomDir=TrgDcmDir,
                                                      SearchString=FromSearchString)
            
            """
            PixArrToCopy needs to be added to previously Direct-copied 
            segmentations in TrgPixArrInSeg.
            """
            
            TrgPixArr,\
            TrgPFFGStoSliceInds = AddCopiedPixArr(OrigPixArr=TrgPixArrInSeg, 
                                                  OrigPFFGStoSliceInds=TrgPFFGStoSliceIndsInSeg, 
                                                  PixArrToAdd=PixArrToCopy, 
                                                  PFFGStoSliceIndsToAdd=[ToSliceNum])
            
        else:
            """
            A Target SEG was not provided.
            """
            TrgPixArr = deepcopy(PixArrToCopy)
            
            TrgPFFGStoSliceInds = deepcopy([ToSliceNum]) # 15/12 verify?
    
    
    if UseCase in ['1b', '2b']:
        
        TrgPixArr = deepcopy(PixArrToCopy)
        
        TrgPFFGStoSliceInds = deepcopy([FromSliceNum]) # 04/12 verify?
    
        
        
    if UseCase in ['3a', '3b', '4', '5']:
        from ConversionTools import PixArr2Image
        from ConversionTools import Image2PixArr
        from ImageTools import ResampleImage
        from ImageTools import GaussianBlurImage
        #from ImageTools import RecursiveGaussianBlurImage
        from ImageTools import BinaryThresholdImage
        from ConversionTools import ConvertImagePixelType
        
        # Import the 3D images:
        if not 'a' in UseCase:
            SrcIm = ImportImage(SrcDcmDir)
            TrgIm = ImportImage(TrgDcmDir)
        
        # Convert PixArrToCopy to a labelmap image:
        LabmapImToCopy = PixArr2Image(PixArr=PixArrToCopy,
                                      FrameToSliceInds=[FromSliceNum],
                                      RefIm=SrcIm)
        
        if LogToConsole:
            #from ImageTools import ImageMax
            #from ImageTools import ImageMin
            unique = UniqueItems(Items=PixArrToCopy, IgnoreZero=False)
            print('\nPrior to resampling:')
            #print(f'type(LabmapImToCopy) = {type(LabmapImToCopy)}')
            #print(f'LabmapImToCopy.GetPixelIDTypeAsString() = {LabmapImToCopy.GetPixelIDTypeAsString()}')
            #print(f'PixArrToCopy.shape = {PixArrToCopy.shape}')
            #print(f'There are {len(unique)} unique items in PixArrToCopy')
            #print(f'Max of LabmapImToCopy = {ImageMax(LabmapImToCopy)}')
            #print(f'Min of LabmapImToCopy = {ImageMin(LabmapImToCopy)}')
            #print(f'LabmapImToCopy.GetSize() = {LabmapImToCopy.GetSize()}')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(LabmapImToCopy, LogToConsole)
            
        
        
    if UseCase in ['3a', '3b', '4', '5']:
        
        #import SimpleITK as sitk
        
        
        """ First try resampling the image using a NearestNeighbor 
        interpolation. """
        interp = 'NearestNeighbor'
        
        if LogToConsole:
            print('\nUsing', interp, 'interpolation...')
            
        # Resample LabmapImToCopy to the Target image's grid:
        ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                          RefImage=TrgIm, 
                                          Interpolation=interp,
                                          LogToConsole=LogToConsole)
        
        if LogToConsole:
            print('\nAfter resampling:')
        
        
        # PixID required for if statement below.
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
            
        
        
        """ 
        If resampling to a smaller grid size aliasing can occur, leading to
        the appearent disappearance of segmentations (i.e. F2Sinds = []). 
        
        Try suggestions from Ziv Yaniv here:
            https://github.com/SimpleITK/SimpleITK/issues/1277
                
        1. Use the sitkLabelGaussian interpolator.
        
        2. Gaussian blur the segmentation, then resample using a linear 
        interpolator, and then threshold so that values in [a<1<b] are mapped
        to 1. You will need to select appropriate values for a, b.
        
        3. Compute the distance map of the segmentation, then resample using a
        linear interpolator, and then threshold the absolute value of the 
        distance map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you 
        need to select an appropriate dist_threshold.

        Tried:
        
        1. Resulting Gaussian interpolated resampled image had zeros only 
        (the image was converted to float prior to resampling).
        
        2. Works for UseCase 3B-i_SEG with FromSliceNum = 26 with binary
        threshold value of 0.04.
        
        
        Comment from Bradley Lowecamp:
            The implementation of the sitkLabelGaussian interpolator in the 
            ResampleImageFilter sets the sigma to the pixel spacing of the 
            input image to the resample filter. This is why it does not do a 
            better job in reducing aliasing.
        """
        #Workaround = 1
        Workaround = 2
        #Workaround = 0 # i.e. no workaround
        
        if not F2Sinds and Workaround:
            if LogToConsole:
                print('\nThe resampled image has no non-zero frames (e.g. due',
                      f'to aliasing).  \nTry workaround #{Workaround}.')
            
            # Convert LabmapImToCopy from 32-bit unsigned integer to a float:
            LabmapImToCopy = ConvertImagePixelType(Image=LabmapImToCopy, 
                                                   NewPixelType='sitkFloat32')
            
            if LogToConsole:
                print('\nAfter converting to float:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(LabmapImToCopy, LogToConsole)
        
            if Workaround == 1:
                if LogToConsole:
                    print(f'\nTrying workaround #{Workaround}...')
                
                interp = 'Gaussian' # <-- get zeros only (image was converted 
                # to float prior to resampling)
                
                if LogToConsole:
                    print('\nUsing', interp, 'interpolation...')
                    
                # Resample LabmapImToCopy to the Target image's grid:
                ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                                  RefImage=TrgIm, 
                                                  Interpolation=interp,
                                                  LogToConsole=LogToConsole)
                
                if LogToConsole:
                    print('\nAfter resampling:')
                
                
                # PixID required for if statement below.
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
                    
                    
                """ Will need to convert back to unsigned int and may need to
                do some binary thresholding, but only proceed if there are
                non-zero frames in the resampled image. """
                
                if not F2Sinds:
                    Workaround = 2 # try Workaround #2
                    
                    
            if Workaround == 2:
                if LogToConsole:
                    print(f'\nTrying workaround #{Workaround}...')
                      
                SrcSpacings = SrcIm.GetSpacing()
                TrgSpacings = TrgIm.GetSpacing()
                
                SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, TrgSpacings))
                
                """ 
                The FWHM of the Gaussian is:
                    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
                    
                Due to the Nyquist-Shannon sampling theorem, the FWHM should be
                at least 2x the voxel spacings of the Target image. For 
                symmetry considerations, it makes sense for one Source voxel 
                width to be blurred to three Target voxels widths.  :
                    
                    sigma = (3*TrgSpacings)/2.355
                    
                """
                
                if LogToConsole:
                    print(f'\nSpacingsRatio = {SpacingsRatio}')
                      
                # Gaussian blur the image:
                """
                For UseCase3B-i_SEG:
                    different image sizes:
                       Source: (208, 256, 50) 
                       Target: (378, 448, 21)
                
                   different voxel sizes:
                       Source: (0.8984375, 0.8984375, 2.9999990463256836) 
                       Target: (0.51339286565781, 0.51339286565781, 7.499998569488525)
                
                   different slice thickness:
                       Source: 3.0 
                       Target: 5.0
                
                Segmentation on slice 26.
                
                
                With Variance = (0, 0, 0.1):
                    --> blurred to 3 slices: [25, 26, 27] with 3 unique values: [0, 0.0055, 0.989]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.0051 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                    
                With Variance = (0, 0, 0.5):
                    --> blurred to 3 slices: [25, 26, 27] with 3 unique values: [0, 0.026, 0.947]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.0247 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                With Variance = (1, 1, 1) (default):
                    --> blurred to 3 slices: [25, 26, 27] with 727 unique values: [0, ..., 0.90]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.047 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                    
                With Variance = (0, 0, 1):
                    --> blurred to 3 slices: [25, 26, 27] with 3 unique values: [0, 0.050, 0.9001]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.047 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                    
                With Variance = (0, 0, 1.5):
                    --> blurred to 3 slices: [25, 26, 27] with 3 unique values: [0, 0.0712, 0.858]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.0663 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                    
                With Variance = (0, 0, 2):
                    --> blurred to 5 slices: [24, 25, 26, 27, 28] with 4 unique values: [0, 0.005, 0.09, 0.811]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.0837 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                    
                With Variance = (0, 0, 2.5):
                    --> blurred to 5 slices: [24, 25, 26, 27, 28] with 4 unique values: [0, 0.00736, 0.106, 0.773]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.10 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                With Variance = (0, 0, 3):
                    --> blurred to 5 slices: [24, 25, 26, 27, 28] with 4 unique values: [0, 0.010, 0.12, 0.737]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.114 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                With Variance = (0, 0, 4):
                    --> blurred to 5 slices: [24, 25, 26, 27, 28] with 4 unique values: [0, 0.016, 0.15, 0.675]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.137 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                With Variance = (0, 0, 5):
                    --> blurred to 5 slices: [24, 25, 26, 27, 28] with 4 unique values: [0, 0.023, 0.166, 0.622]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.156 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                With Variance = (0, 0, 6):
                    --> blurred to 5 slices: [24, 25, 26, 27, 28] with 4 unique values: [0, 0.030, 0.182, 0.576]
                    --> lin resampled to 2 slices: [10, 11] with Max = 0.171 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                With Variance = (0, 0, 7):
                    --> blurred to 7 slices: [23, 24, 25, 26, 27, 28, 29] with 5 unique values: [0, 0.0047, 0.037, 0.192, 0.532]
                    --> lin resampled to 4 slices: [9, 10, 11, 12] with Max = 0.183 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                
                
                With Variance = (0, 0, 10):
                    --> blurred to 7 slices: [23, 24, 25, 26, 27, 28, 29] with 5 unique values: [0, 0.010, 0.056, 0.213, 0.440]
                    --> lin resampled to 4 slices: [9, 10, 11, 12] with Max = 0.203 
                    --> binary thresholded to 2 slices: [10, 11] with thresh = max/2
                
                
                """
                #LabmapImToCopy = GaussianBlurImage(LabmapImToCopy, 
                #                                   Variance=(0, 0, 7))
                LabmapImToCopy = GaussianBlurImage(Im=LabmapImToCopy,
                                                   Variance=Sigma,
                                                   LogToConsole=LogToConsole)
                
                # Use the RecursiveGaussian image filter:
                """
                For UseCase4_SEG_S3_to_S7:
                    
                With SrcVoxelSize[2] = 6.25 and TrgVoxelSize[2] = 3.30, 
                i.e. SrcVoxelSize[2] ~= TrgVoxelSize[2]/2 and input slice = 3:
                    Sigma = 0.1 -->             1 blurred slice: [3]                                 --> linearly resampled to 0 slices: []
                    Sigma = 0.2 -->             3 blurred slices: [2, 3, 4]                          --> linearly resampled to 0 slices: []
                    Sigma = 0.3 -->             3 blurred slices: [2, 3, 4]                          --> linearly resampled to 0 slices: []
                    Sigma = 0.4 -->             5 blurred slices: [0, 1, 3, 5, 6]                    --> linearly resampled to 8 slices: []
                    Sigma = 0.5 --> Max = 1 and 6 blurred slices: [0, 2, 3, 4, 6, 7]                 --> linearly resampled to 0 slices: []
                    Sigma = 0.6 -->             7 blurred slices: [0, 1, 2, 3, 4, 5, 6]              --> linearly resampled to 0 slices: []
                    Sigma = 0.7 -->             3 blurred slices: [2, 3, 4]                          --> linearly resampled to 0 slices: []
                    Sigma = 0.8 -->             4 blurred slices: [3, 7, 8, 9]                       --> linearly resampled to 4 slices: [0, 1, 2, 3]
                    Sigma = 0.9 -->             6 blurred slices: [0, 3, 6, 7, 10, 11]               --> linearly resampled to 8 slices: [0, 1, 2, 3, 4, 5, 6, 7]
                    Sigma = 1.0 --> Max = 1 and 10 blurred slices: [0, 1, 3, 5, 6, 8, 9, 11, 12, 14] --> linearly resampled to 13 slices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
                    Sigma = 1.1 --> 9 blurred slices: [1, 3, 5, 7, 8, 10, 13, 15, 16]    --> linearly resampled to 16 slices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
                    
                    Sigma = 1.5 --> Max = 1    and  9 blurred slices: [0, 3, 6, 8, 10, 12, 15, 17, 19]                  --> lin resampled to 23 slices with Max = 3e-13
                    Sigma = 1.7 --> Max = 1    and 12 blurred slices: [0, 2, 3, 4, 6, 8, 9, 11, 14, 16, 19, 21]         --> lin resampled to 23 slices with Max = 3e-12
                    Sigma = 1.8 --> Max = 1    and 13 blurred slices: [0, 2, 3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22]    --> lin resampled to 23 slices with Max = 4e-13
                    Sigma = 1.9 --> Max = 0.99 and 14 blurred slices: [0, 2, 3, 4, 6, 7, 9, 12, 15, 17, 18, 20, 21, 23] --> lin resampled to 23 slices with Max = 3e-12
                    Sigma = 2.0 --> Max = 0.99 and 17 blurred slices: [0, 2, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24] --> lin resampled to 23 slices with Max = 6e-12
                    
                    Sigma = 3.0 --> Max = 0.82 and 11 blurred slices: [2, 3, 4, 8, 9, 12, 13, 17, 18, 21, 22] --> lin resampled to 23 slices with Max = 1e-6
                    
                    Sigma = 4.0 --> Max = 0.62 and 14 blurred slices: [1, 2, 3, 4, 5, 9, 10, 11, 15, 16, 17, 21, 22, 23] --> lin resampled to 23 slices with Max = 1e-6
                    
                """
                #LabmapImToCopy = RecursiveGaussianBlurImage(LabmapImToCopy, 
                #                                            Sigma=3,
                #                                            Direction=2)
                
                if LogToConsole:
                    print('\nAfter applying Gaussian blur:')
                
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(LabmapImToCopy, LogToConsole)
                
                #if not F2Sinds:
                #    Workaround = 3 # try Workaround #3
                    
                # Linearly resample the blurred image:
                interp = 'Linear'
                
                if LogToConsole:
                    print('\nUsing', interp, 'interpolation...')
                    print(f'   LabmapImToCopy.GetSize() = {LabmapImToCopy.GetSize()}')
                    print(f'   LabmapImToCopy.GetSpacing() = {LabmapImToCopy.GetSpacing()}')
                    print(f'   TrgIm.GetSize() = {TrgIm.GetSize()}')
                    print(f'   TrgIm.GetSpacing() = {TrgIm.GetSpacing()}')
                    
                # Resample LabmapImToCopy to the Target image's grid:
                ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                                  RefImage=TrgIm, 
                                                  Interpolation=interp,
                                                  LogToConsole=LogToConsole)
                
                if LogToConsole:
                    print('\nAfter resampling:')
                
                
                # UniqueVals is required below:
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
                    
                
                
                """ Prior to converting back to unsigned integer, must first
                binary threshold the image. """
                # Thresh values with results for Case3B-i_SEG and 
                # FromSliceNum = 26 (max value in Pix array prior to binary
                # thresholding = 0.0468482):
                #Thresh = 0.5 # no frames
                #Thresh = 0.005 # two frames (slices 10 & 11)
                #Thresh = 0.04 # one frame (slice 11)
                #Thresh = 0.5*max(UniqueVals) # 13/01/2021
                #Thresh = 0.75*max(UniqueVals) # 13/01/2021
                
                """
                Get 2 slices with 0.5*max
                Get 1 slice with 0.75*max
                """
                
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                if LogToConsole:
                    #print(f'\nBinary thresholding with Thresh = {Thresh} ',
                    #      f'(max = {max(UniqueVals)})...')
                    print('\nBinary thresholding at',
                          f'{ThreshLevel}*{max(UniqueVals)} = {ThreshValue}...')
                
                #ResLabmapImToCopy = BinaryThresholdImage(ResLabmapImToCopy, 
                #                                         Thresh)
                ResLabmapImToCopy = BinaryThresholdImage(Im=ResLabmapImToCopy, 
                                                         Thresh=ThreshValue)
                
                if LogToConsole:
                    print('\nAfter binary thresholding:')
                    
                
                # PixID will be needed for if statement below.
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
                    
                
                """ Following code was moved outside loop since useful to run
                under all conditions. """
                #""" ... then convert back to unsigned integer. """
                #if not PixID == 5:
                #    print(f'\nPixelID (= {PixID}) != 5 (sitkUInt32 = Unsigned',
                #          '32 bit integer).  Converting to sitkUInt32..')
                #    
                #    # Convert ResLabmapImToCopy from float to 32-bit unsigned 
                #    # integer:
                #    ResLabmapImToCopy = ConvertImagePixelType(Image=ResLabmapImToCopy, 
                #                                              NewPixelType='sitkUInt32')
                #    
                #    if True:#LogToConsole:
                #        print('\nAfter converting to unsigned int:')
                #        PixID, PixIDTypeAsStr, UniqueVals,\
                #        F2Sinds = GetImageInfo(ResLabmapImToCopy, 
                #                               LogToConsole=True)
                
                
            #if Workaround == 3:
            #    print(f'\nTrying workaround #{Workaround}...')
            #          
            #    print('\nStill need to write Workaround 3...')
                
        
        
        #""" 16/12:
        #    Instead of using NearestNeighbor, use a BSpline interpolation, then
        #    binary threshold the image using a low threshold (much lower than
        #    0.5).
        #"""
        #interp = 'BSpline'
        #interp = 'Linear'
        #interp = 'NearestNeighbor'
        #interp = 'Gaussian' # <-- get zeros only (image was converted to float 
        # prior to resampling)
        
        
            
        
        #""" For some reason when binary thresholding an empty image, the 
        #threshold value SetLowerThreshold is being rounded down to 0, resulting
        #in a non-binary image.  So only apply binary thresholding if there are
        #at least 2 unique values in the input image. """
        if False:
            PixArr, F2Sinds = Image2PixArr(LabmapIm=ResLabmapImToCopy)
            unique = UniqueItems(Items=PixArr, IgnoreZero=False)
            
            if len(unique) < 2:
                ResLabmapImToCopy = BinaryThresholdImage(ResLabmapImToCopy, Thresh=0.5)
                
                if LogToConsole:
                    print('\nAfter binary thresholding:')
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
            
        
        
        """ Ensure that ResLabmapImToCopy is a 32-bit unsigned integer. """
        if not PixID == 5:
            if LogToConsole:
                print(f'\nPixelID = {PixID} ({PixIDTypeAsStr})). Converting',
                      f'to unsigned 32-bit integer (sitkUInt32)..')
            
            # Convert ResLabmapImToCopy from float to 32-bit unsigned 
            # integer:
            ResLabmapImToCopy = ConvertImagePixelType(Image=ResLabmapImToCopy, 
                                                      NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print('\nAfter converting to 32-bit unsigned int:')
            
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
                        
                        
        # Convert ResLabmapImToCopy back to a pixel array:
        ResPixArrToCopy,\
        PFFGStoSliceIndsToCopy = Image2PixArr(LabmapIm=ResLabmapImToCopy)
        
        if LogToConsole:
            print('\nAfter converting to a pixel array:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
            
            for i in range(ResPixArrToCopy.shape[0]):
                maxval = np.amax(ResPixArrToCopy[i])
                minval = np.amin(ResPixArrToCopy[i])
                
                print(f'Frame {i} has max = {maxval}, min = {minval}')
        
        #from ImageTools import BinaryThresholdImage
        #
        #if not 'earest' in interp:
        #    ResLabmapImToCopy = BinaryThresholdImage(ResLabmapImToCopy, 
        #                                             Thresh=0.5)
        
    
        
        
    if UseCase == '3a':
        """
        For a Direct copy, one segmentation is always copied to a single
        segmentation, so if after resampling there is more than one frame, the
        frames will need to be averaged.
        """
        from GeneralTools import MeanPixArr
        
        if LogToConsole:
            unique = UniqueItems(Items=ResPixArrToCopy, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy')
        
        F = ResPixArrToCopy.shape[0]
        
        if LogToConsole:
            print(f'\nThe segmentation on FromSliceNum = {FromSliceNum} has',
                  f'been resampled to {F} frames: {PFFGStoSliceIndsToCopy}')
        
        if F > 1:
            if LogToConsole:
                print(f'\nResPixArrToCopy has {F} frames so the pixel arrays',
                      'will be averaged...')
            
            # Take pixel-by-pixel mean across all frames in ResPixArrToCopy:
            ResPixArrToCopy = MeanPixArr(ResPixArrToCopy)
            
            if LogToConsole:
                F = ResPixArrToCopy.shape[0]
            
                print(f'\nResPixArrToCopy now has {F} frames')
        else:
            if LogToConsole:
                print(f'\nResPixArrToCopy has F frames')
            
        
        # Shift the in-plane elements in ResPixArrToCopy:
        ResPixArrToCopy = ShiftFrame(Frame=ResPixArrToCopy, PixShift=PixShift)
        
        if LogToConsole:
            unique = UniqueItems(Items=ResPixArrToCopy, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy',
                  'after shifting the frame.')
        
        
        if TrgSeg:
            """
            Since this is a Direct copy, any existing segmentations in Target
            with segment label matching FromSearchString are to be preserved.
            """
            
            # Get all frames in the Target segment containing the segmentation 
            # to be copied:
            TrgPixArrInSeg,\
            TrgPFFGStoSliceIndsInSeg = GetPixArrInSeg(Seg=TrgSeg, 
                                                      DicomDir=TrgDcmDir,
                                                      SearchString=FromSearchString)
            
            """
            ResPixArrToCopy needs to be added to previously Direct-copied 
            segmentations in TrgPixArrInSeg.
            """
            
            TrgPixArr,\
            TrgPFFGStoSliceInds = AddCopiedPixArr(OrigPixArr=TrgPixArrInSeg, 
                                                  OrigPFFGStoSliceInds=TrgPFFGStoSliceIndsInSeg, 
                                                  PixArrToAdd=ResPixArrToCopy, 
                                                  PFFGStoSliceIndsToAdd=[ToSliceNum])
        else:
            """
            A Target SEG was not provided.
            """
            TrgPixArr = deepcopy(ResPixArrToCopy)
            
            TrgPFFGStoSliceInds = deepcopy([ToSliceNum]) # 04/12 verify?
            
    
    
    if UseCase in ['3b', '4']:
        TrgPixArr = deepcopy(ResPixArrToCopy)
        
        TrgPFFGStoSliceInds = deepcopy(PFFGStoSliceIndsToCopy)
        
        
    
    #if UseCase == '4':
    
        
    if UseCase == '5':
        from ImageTools import RegisterImages
        from ImageTools import TransformImage
        from ImageTools import BinaryThresholdImage
    
        # Register the images using an affine transformation:
        RegSrcIm, RegImFilt = RegisterImages(FixIm=TrgIm, MovIm=SrcIm, 
                                             Tx='affine',
                                             LogToConsole=LogToConsole)
        
        # Transform LabmapImToCopy using RegImFilt:
        TxLabmapImToCopy = TransformImage(Im=LabmapImToCopy, 
                                          RegImFilt=RegImFilt,
                                          Interpolation='Nearestneighbor')
        
        if LogToConsole:
            print('\nAfter transforming LabmapImToCopy by the registration',
                  'transformation:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabmapImToCopy, LogToConsole)
        
        """
        21/12/20:  
            Following the use of a nearestneighbor interpolator when 
            transforming LabmapImToCopy, i.e.:
            
            TxMap[0]["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]
            
            it's no longer necessary to perform the binary thresholding step.
        
        # Binary threshold TxLabmapImToCopy:
        Thresh = 0.1
        
        TxLabmapImToCopy = BinaryThresholdImage(Im=TxLabmapImToCopy, 
                                                Thresh=Thresh)
        
        if LogToConsole:
                print('\nAfter binary thresholding with Thresh = {Thresh}:')
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(TxLabmapImToCopy, LogToConsole)
        """
        
        
        # Convert TxLabmapImToCopy to a pixel array:
        TrgPixArr,\
        TrgPFFGStoSliceInds = Image2PixArr(LabmapIm=TxLabmapImToCopy)
        
        if LogToConsole:
            unique = UniqueItems(Items=TrgPixArr, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TrgPixArr')
            
            F = TrgPixArr.shape[0]
            
            print(f'\nThe segmentation on slice {FromSliceNum} has been',
                  f'registered to {F} frames: {TrgPFFGStoSliceInds}.')
    
    
    
    #if 'a' in UseCase:
    #    """ 
    #    For UseCases 1a/b, 2a/b or 3a, there will always be a single 
    #    PerFrameFunctionalGroupsSequence, and the single index will be slice
    #    ToSliceNum.  Over-write previous TrgPFFGStoSliceInds. 
    #    """
    #    
    #    TrgPFFGStoSliceInds = [ToSliceNum]
    #    
    #    print(f'\nTrgPFFGStoSliceInds over-written to {TrgPFFGStoSliceInds}')
        
        
    """ 
    If UseCase is 1a/2a/3a (i.e. a Direct copy is to be made) and a SEG has 
    been provided for Target this is not a clean copy, i.e. the copied 
    segmentation is to be added to an existing segment containing previously 
    copied segmentation(s). If a Relationship-preserving copy is to be made or 
    if a Target SEG has not been provided, the Source SEG will be used as a 
    template.
    """
    #if not 'a' in UseCase or not TrgSeg:
    #    TrgSeg = InitialiseSeg(SegTemplate=SrcSeg, SegNum=FromSegNum,
    #                           PFFGStoSliceInds=TrgPFFGStoSliceInds, 
    #                           DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    if LogToConsole:
        print('\n\nInputs to InitialiseSeg:')
        print(f'   SegNum = {FromSegNum}')
        print(f'   TrgPFFGStoSliceInds = {TrgPFFGStoSliceInds}')
        
    """ 
    09/12: When downsampling some resulting images have no segmentations,
    hence TrgPFFGStoSliceInds = [].  This results in 
    Seg.NumberOfFrames = f"{len(PFFGStoSliceInds)}"
    in InitialiseSeg.
    Until the cause of the above is determined, raise an exception if 
    TrgPFFGStoSliceInds = [].
    """
    if TrgPFFGStoSliceInds == []:
        msg = 'TrgPFFGStoSliceInds = []. There must be at least one slice '\
              + 'number in TrgPFFGStoSliceInds.'
        
        raise Exception(msg)
    
        
    if TrgSeg:
        """ Use TrgSeg to initialise the Target SEG. """
        #TrgSeg = InitialiseSeg(SegTemplate=TrgSeg, SegNum=ToSegNum,
        TrgSeg = InitialiseSeg(SegTemplate=TrgSeg, SearchString=FromSearchString,
                               PFFGStoSliceInds=TrgPFFGStoSliceInds, 
                               DicomDir=TrgDcmDir, NamePrefix=AddText,
                               LogToConsole=LogToConsole)
    else:
        """ Use SrcSeg to initialise the Target SEG. """
        TrgSeg = InitialiseSeg(SegTemplate=SrcSeg, SearchString=FromSearchString,
                               PFFGStoSliceInds=TrgPFFGStoSliceInds, 
                               DicomDir=TrgDcmDir, NamePrefix=AddText,
                               LogToConsole=LogToConsole)
        
    if LogToConsole:
        print(f'\n\n\n***FromSegNum = {FromSegNum}')
    
        print('\n\nInputs to ModifySeg:')
        print(f'   TrgPixArr.shape = {TrgPixArr.shape}')
        [print(f'   np.amax(TrgPixArr[{i}]) = {np.amax(TrgPixArr[i])}') for i in range(TrgPixArr.shape[0])]
        print(f'   TrgPFFGStoSliceInds = {TrgPFFGStoSliceInds}')
        
    # Modify the tags in TrgSeg:
    TrgSeg = ModifySeg(Seg=TrgSeg, PixArr=TrgPixArr,
                       PFFGStoSliceInds=TrgPFFGStoSliceInds, 
                       DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    
    """ Had 1040 lines of code. Now 207. """
    
    return TrgSeg






"""
******************************************************************************
******************************************************************************
COPY A CONTOUR / SEGMENTATION
******************************************************************************
******************************************************************************
"""

def CopyRoi(SrcRoi, FromSearchString, SrcDcmDir, FromSliceNum, TrgDcmDir,
            TrgRoi=None, ToSliceNum=None, Sigma=(1,1,1), ThreshLevel=0.75,
            AddText='', LogToConsole=False):
    """
    
    Inputs:
    ******
    
    SrcRoi : Pydicom object
        Source RTS/SEG object.
                             
    FromSearchString : string
        All or part of the Source ROIName/SegmentLabel of the ROI/segment 
        containing the contour/segmentation to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied (counting from 0).
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
        
    TrgRts : Pydicom object (optional; None by default)
        Target RTS object that the contour(s)/segmentation(s) is to be copied 
        to. TrgRts is only non-None if a Direct copy of the contour/
        segmentation is to be made and added to existing Direct-copied 
        contour(s)/segmentation(s).     
        
    ToSliceNum : integer (optional; None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to (counting from 0).  This only applies 
        for Direct copies, hence the default value None.
        
    Sigma : tuple (optional; (1,1,1) by default)
        The sigma values that will define the Gaussian if Gaussian image 
        blurring is used.  The tuple will define the sigma along x, y and z
        directions.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if used.  The threshold value will be ThreshLevel*M, where
        M is the maximum intensity value of the image.
    
    AddText : string (optional, '' by default)
        Sting of text to pass to ModifyRts/ModifySeg to be used when generating
        a filename for the new Target RTS/SEG.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    *******
        
    TrgRts : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing Direct-
        copied contour(s)/segmentation(s)) Target RTS/SEG object.
    """

    import time
    from DicomTools import IsSameModalities
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Compare the modalities of the RTSs/SEGs:
    SrcModality = SrcRoi.Modality
    
    if TrgRoi:
        TrgModality = TrgRoi.Modality
        
        if not IsSameModalities(SrcRoi, TrgRoi):
            msg = f"The Source ({SrcModality}) and Target ({TrgModality}) " \
                  + "modalities are different."
            
            raise Exception(msg)
            
    
    if SrcModality == 'RTSTRUCT':
        NewTrgRoi = CopyRts(SrcRoi, FromSearchString, SrcDcmDir, FromSliceNum, 
                            TrgDcmDir, TrgRoi, ToSliceNum, Sigma, ThreshLevel,
                            AddText, LogToConsole)
    
    elif SrcModality == 'SEG':
        NewTrgRoi = CopySeg(SrcRoi, FromSearchString, SrcDcmDir, FromSliceNum, 
                            TrgDcmDir, TrgRoi, ToSliceNum, Sigma, ThreshLevel,
                            AddText, LogToConsole)
        
    else:
        msg = f'The Source modality ({SrcModality}) must be either "RTS" or '\
              + '"SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n\nDone.  Took {Dtime} s to run.')
        
    return NewTrgRoi







"""
******************************************************************************
******************************************************************************
CHECK FOR ERRORS IN NEW TARGET RTS / SEG
******************************************************************************
******************************************************************************
"""

def ErrorCheckRoi(Roi, DicomDir, LogToConsole=False):
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
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe info or any errors that are found.
    
    Nerrors : integer
        The number of errors found.
    """
    
    import time
    import os
    
    if Roi.Modality == 'RTSTRUCT':
        import RtsTools
        import importlib
        importlib.reload(RtsTools)
        from RtsTools import ErrorCheckRts
        
        LogList, Nerrors = ErrorCheckRts(Roi, DicomDir, LogToConsole)
        
    elif Roi.Modality == 'SEG':
        import importlib
        import SegTools
        importlib.reload(SegTools)
        from SegTools import ErrorCheckSeg
        
        LogList, Nerrors = ErrorCheckSeg(Roi, DicomDir, LogToConsole)
        
    else:
        LogList = None
        Nerrors = None
    
    
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    Fname = CurrentDateTime + '_RoiCopy.log'
    
    CWD = os.getcwd()
    
    Fpath = os.path.join(CWD, 'logs', Fname)
    
    with open(Fpath, 'w') as f:
        for line in LogList:
            f.write(line)
        
    print('\nLog file saved to:\n\n', Fpath)
    
    return LogList, Nerrors






"""
******************************************************************************
******************************************************************************
EXPORT NEW TARGET RTS / SEG TO DISK
******************************************************************************
******************************************************************************
"""

def ExportTrgRoi(TrgRoi, SrcRoiFpath, ExportDir, NamePrefix=''):
    """
    Export RTS/SEG to disk.  
    
    Inputs:
    ******
    
    TrgRoi : Pydicom object
        Target RTS/SEG ROI object to be exported.
        
    SrcRoiFpath : string
        Full path of the Source RTS/SEG file (used to generate the filename of
        the new RTS/SEG file).
                              
    ExportDir : string
        Directory where the new RTS/SEG is to be exported.
    
    NamePrefix : string (optional; '' by default)
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
                           
    
    Outputs:
    *******
    
    TrgRoiFpath : string
        Full path of the exported Target RTS/SEG file.
    """
    
    import os
    import time
    
    # Get the filename of the original RTS/SEG file:
    #SrcRoiFname = os.path.split(SrcRoiFpath)[1]
    
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #FnamePrefix = CurrentDateTime + '_' + NamePrefix + '_from_' 
    FnamePrefix = CurrentDateTime + '_' + NamePrefix
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = FnamePrefix + SrcRoiFname
    #TrgRoiFname = FnamePrefix
    TrgRoiFname = FnamePrefix + '.dcm' # 'dcm' was missing (found 18/12)
    
    TrgRoiFpath = os.path.join(ExportDir, TrgRoiFname)
    
    TrgRoi.save_as(TrgRoiFpath)
        
    print('\nNew Target RTS/SEG exported to:\n\n', TrgRoiFpath)
    
    return TrgRoiFpath





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
        
    print('\nNew DICOM Registration Object exported to:\n\n', DroFpath)
    
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
        print(f'\nTook {Dtime} s to import the 3D images.')
        
    # The Transform parameter file:
    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    # Parse the Transform parameter file:
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\nTook {Dtime} s to parse the transform parameter file.')

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
        print(f'\n   Took {Dtime} s to get the transform matrix type.')
        
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
        print(f'\n   Took {Dtime} s to add sequences to',
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
        print(f'\n   Took {Dtime} s to update the ReferencedInstanceSequences',
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
        print(f'\n   Took {Dtime} s to add sequences to',
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
        print(f'\n   Took {Dtime} s to update the ReferencedInstanceSequences',
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
        print(f'\n   Took {Dtime} s to add sequences to',
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
        print(f'\n   Took {Dtime} s to modify the',
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
        print(f'\n   Took {Dtime} s to add sequences to',
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
        print(f'\n   Took {Dtime} s to modify the',
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
        print(f'\n   Took {Dtime} s to delete all GroupLength tags.')
        
    
    Dtime = round(times[-1] - times[0], 1)
    if True:#LogToConsole:
        print(f'\n   Took a total of {Dtime} s to create the DRO.')
        
        
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
        
    print('\nNew DICOM Registration Object exported to:\n\n', DroFpath)
    
    return DroFpath