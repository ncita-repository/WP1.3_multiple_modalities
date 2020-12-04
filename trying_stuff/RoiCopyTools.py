# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:12:23 2020

@author: ctorti
"""


"""
ROI copying tools
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
    ------
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    ToSliceNum : integer
        Slice index within the Target DICOM stack where the segmentation is to
        be copied to (counting from 0).  This only applies for Direct copies.
                        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.

              
    Output:
    ------
        
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
            print('\n   same image size')
        else:
            print('\n   different image sizes:\n',
                  f'       Source: {SrcSize} \n',
                  f'       Target: {TrgSize}')
        
        if SameSpacings:
            print('\n   same voxel sizes')
        else:
            print('\n   different voxel sizes:\n',
                  f'       Source: {SrcSpacing} \n',
                  f'       Target: {TrgSpacing}')
        
        if SameOrigins:
            print('\n   same origin')
        else:
            print('\n   different origins:\n',
                  f'       Source: {SrcIPPs[0]} \n',
                  f'       Target: {TrgIPPs[0]}')
            
        if SameDirs:
            print('\n   same patient orientation')
        else:
            print('\n   different patient orientations:\n',
                  f'      Source: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},',
                  f'               {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},',
                  f'               {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]\n',
                  f'      Target: [{TrgDirs[0]}, {TrgDirs[1]}, {TrgDirs[2]},',
                  f'               {TrgDirs[3]}, {TrgDirs[4]}, {TrgDirs[5]},',
                  f'               {TrgDirs[6]}, {TrgDirs[7]}, {TrgDirs[8]}]')
        
        
    return UseCase






"""
******************************************************************************
******************************************************************************
COPY A CONTOUR
******************************************************************************
******************************************************************************
"""

def CopyRts(SrcRts, FromRoiLabel, SrcDcmDir, FromSliceNum, TrgDcmDir,
            TrgRts=None, ToSliceNum=None, LogToConsole=False):
    """
    Make a Direct or Relationship-preserving copy of a single contour in a
    given ROI on a given slice to a new ROI in an existing RTS, or 
    in a newly created RTS (if an existing RTS is not provided). 
    
    For Direct copies the in-plane coordinates (currently the (x,y) coords) 
    will be preserved only, whilst for a Relationship-preserving copy, all
    coordinates are spatially preserved.
    
    There are five possible Main Use Cases (indicated by the numbers 1 to 5)
    and two Sub-cases (indicated by the letters a or b) that apply for a subset
    of the Main Use Cases that are accommodated:
    
    1a. Direct copy of the Source contour in FromRoiLabel on slice  
        FromSliceNum to a new ROI for slice ToSliceNum in Target. 
        The Source and Target images are the same. 
    
    1b. Relationship-preserving copy of the Source segmentation in FromRoiLabel 
        on slice FromSliceNum to a new ROI for slice ToSliceNum in Target. 
        The Source and Target images are the same. 
        
    2a. Direct copy of the Source contour in FromRoiLabel on slice 
        FromSliceNum to a new ROI for slice ToSliceNum in Target. 
        The Source and Target images have the same FOR, IOP, IPP, and VS.
    
    2b. Relationship-preserving copy of the Source contour in FromRoiLabel 
        on slice FromSliceNum to a new ROI for slice ToSliceNum in Target. 
        The Source and Target images have the same FOR, IOP, IPP, and VS. 
    
    3a. Direct copy of the Source contour in FromRoiLabel on slice 
        FromSliceNum to a new ROI for slice ToSliceNum in Target. 
        The Source and Target images have the same FOR and IOP but have 
        different VS (they will likely also have different IPP).
        
    3b. Relationship-preserving copy of the Source contour in FromRoiLabel 
        on slice FromSliceNum to a new ROI for any number of slices in 
        Target.
        The Source and Target images have the same FOR and IOP but have 
        different VS (they will likely also have different IPP).
    
    4.  Relationship-preserving copy of the Source contour in FromRoiLabel 
        on slice FromSliceNum to a new ROI for any number of slices in 
        Target.
        The Source and Target images have the same FOR but have different IOP 
        and IPP, and possibly different VS.
       
    5.  Relationship-preserving copy of the Source contour in FromRoiLabel 
        on slice FromSliceNum to a new ROI for any number of slices in 
        Target.
        The Source and Target images have different FOR, IPP and IOP, and 
        possibly different VS.
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, VS = Voxel spacings
    
    
    
    Inputs:
    ------
    
    SrcRts : Pydicom object
        Source RTS object.
                             
    FromRoiLabel : string
        All or part of the Source ROI label containing the contour to be
        copied.
    
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
                        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    -------
        
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
    
    import time
    from copy import deepcopy
    from GeneralTools import PrintTitle
    from GeneralTools import UniqueItems
    from DicomTools import GetRoiNum
    from RtsTools import GetPtsInContour
    from DicomTools import IsSameModalities
    from GeneralTools import GetPixelShiftBetweenSlices
    from GeneralTools import ShiftFrame
    from GeneralTools import MeanPixArr
    from RtsTools import InitialiseRts
    from RtsTools import ModifyRts
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Set the threshold level for binarising of the transformed labelmap (for
    # UseCase 5):
    Thresh = 0.1
    
    # Verify that FromRoiLabel exists:
    #SrcRoiLabels = GetRoiLabels(SrcRts)
    FromRoiNum = GetRoiNum(Roi=SrcRts, Label=FromRoiLabel)
    
        
    # Compare the modalities of the ROIs:
    if TrgRts:
        if not IsSameModalities(SrcRts, TrgRts):
            msg = f"The Source ({SrcRts.Modality}) and Target " \
                  + "({TrgRts.Modality} modalities are different."
            
            raise Exception(msg)
    
    # The contour to be copied:
    PtsToCopy, CStoSliceIndsToCopy = GetPtsInContour(Rts=SrcRts, 
                                                     DicomDir=SrcDcmDir, 
                                                     RoiLabel=FromRoiLabel, 
                                                     SliceNum=FromSliceNum)
    
    print(f'\nThe contour to be copied from slice {FromSliceNum} has',
          f'{len(PtsToCopy)} points and CStoSliceIndsToCopy =',
          f'{CStoSliceIndsToCopy}.')
    
    
    # Determine which Use Case applies:
    UseCase = WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole)
    
    if True:#LogToConsole:
        PrintTitle(f'Case {UseCase} applies.')
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    if 'a' in UseCase:
        from RtsTools import GetPtsInRoi
        from RtsTools import AddCopiedPtsByCnt
        from ImageTools import GetImageAttributes
        from ConversionTools import ContourToIndices
        from GeneralTools import ChangeZinds
        from ConversionTools import Inds2Pts
        
        
    if UseCase in ['1a', '2a']:
        """
        The contour will be copied to a different slice location, so the
        z-components of the coordinates will need to be modified.  For 
        UseCase 3a this operation will be done in index space (i.e. by
        shifting pixels).
        """
        
        SrcSize, SrcSpacings, SrcST,\
        SrcIPPs, SrcDirs = GetImageAttributes(DicomDir=SrcDcmDir)
        
        # Convert PtsToCopy to indices:
        IndsToCopy = ContourToIndices(Points=PtsToCopy, Origin=SrcIPPs[0], 
                                      Directions=SrcDirs, 
                                      Spacings=SrcSpacings)
        
        # Modify the z-component (k) of the indices to ToSliceNum:
        IndsToCopy = ChangeZinds(Indices=IndsToCopy, NewZind=ToSliceNum)
        
        # Convert back to physical points:
        PtsToCopy = Inds2Pts(Indices=IndsToCopy, Origin=SrcIPPs[0], 
                             Directions=SrcDirs, Spacings=SrcSpacings)
        
        """
        Since this is a Direct copy, any existing contours in Target with ROI
        label matching FromRoiLabel are to be preserved.
        """
        
        # Get all points (by contour) in the Target ROI containing the contour
        # to be copied:
        TrgPtsInRoiByCnt,\
        TrgCStoSliceIndsInRoi = GetPtsInRoi(Rts=TrgRts, 
                                            DicomDir=TrgDcmDir,
                                            RoiLabel=FromRoiLabel)
        
        print('\nThe Target ROI containing the contour to be copied has',
              f'{len(TrgPtsInRoiByCnt)} contours. \nTrgCStoSliceIndsInRoi =',
              f'{TrgCStoSliceIndsInRoi}.')
        
        
        """
        PtsToCopy need to be added to previously Direct-copied contour points
        in TrgPtsInRoiByCnt.
        
        Note:
            Since PtsToCopy is a list of points it needs to be enclosed in []
            so that it is a list (of length 1) of a list of points so that it
            is in the form expected of PtsToAddByCnt.
        """
        
        TrgCntDataByCnt,\
        TrgPtsByCnt,\
        TrgCStoSliceInds = AddCopiedPtsByCnt(OrigPtsByCnt=TrgPtsInRoiByCnt, 
                                             OrigCStoSliceInds=TrgCStoSliceIndsInRoi, 
                                             PtsToAddByCnt=[PtsToCopy], 
                                             CStoSliceIndsToAddByCnt=CStoSliceIndsToCopy)
        
    
    if UseCase in ['1b', '2b']:
        """
        TrgPtsByCnt is simply [PtsToCopy], so that it is a list (of length 1, 
        for the single contour).  Likewise for TrgCntDataByCnt.
        """
        
        from ConversionTools import Pts2ContourData
        
        TrgPtsByCnt = [deepcopy(PtsToCopy)]
        
        TrgCntDataByCnt = [Pts2ContourData(PtsToCopy)]
        
        
    
    if UseCase in ['3a', '3b', '4', '5']:
        from ImageTools import ImportImage
        from RtsTools import GetMaskFromRts
        from ConversionTools import PixArr2Image
        from ConversionTools import Image2PixArr
        #from ConversionTools import PixArr2Contours
        from ConversionTools import PixArr2PtsByContour
        
        # Import the 3D images:
        SrcIm = ImportImage(SrcDcmDir)
        TrgIm = ImportImage(TrgDcmDir)
        
        
        """
        ***********************************************************************
        Convert the contour to be copied to a 3D labelmap image.
        ***********************************************************************
        """    
        # Convert the contour data to be copied to a 2D mask:
        PixArrToCopy = GetMaskFromRts(Rts=SrcRts, 
                                      FromRoiLabel=FromRoiLabel, 
                                      FromSliceNum=FromSliceNum, 
                                      DicomDir=SrcDcmDir, 
                                      RefImage=SrcIm)
        
        print(f'\nPixArrToCopy.shape = {PixArrToCopy.shape}')
        
        # Convert PixArrToCopy to a labelmap image:
        LabmapImToCopy = PixArr2Image(PixArr=PixArrToCopy,
                                      FrameToSliceInds=[FromSliceNum],
                                      RefIm=SrcIm)
        
        print(f'\nLabmapImToCopy.GetSize() = {LabmapImToCopy.GetSize()}')
        
    
    if '3' in UseCase:
        from ImageTools import ResampleImage
        
        # Resample LabmapImToCopy to the Target image's grid:
        ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                          RefImage=TrgIm, 
                                          Interpolation='NearestNeighbor')
        
        print(f'\nResLabmapImToCopy.GetSize() = {ResLabmapImToCopy.GetSize()}')
        
        # Convert ResLabmapImToCopy to a pixel array:
        ResPixArrToCopy,\
        TrgCStoSliceInds = Image2PixArr(LabmapIm=ResLabmapImToCopy)
        
        print(f'\nResPixArrToCopy.shape = {ResPixArrToCopy.shape}')
        
        unique = UniqueItems(Items=ResPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy')
        
        F = ResPixArrToCopy.shape[0]
        
        print(f'\nThe segmentation (from the contour) on slice {FromSliceNum}',
              f'has been resampled to {F} frames (contours).\nTrgCStoSliceInds',
              f'= {TrgCStoSliceInds}.')
        
    
        
    if UseCase == '3a':
        """
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled labelmap-image-to-copy).  So if after resampling there is 
        more than one frame, the frames will be averaged.
        """
        
        if F > 1: 
            print(f'The {F} frames will be averaged...')
            
            # The pixel-by-pixel mean across all frames in ResPixArrToCopy:
            ResPixArrToCopy = MeanPixArr(ResPixArrToCopy)
            
            print(f'\nResPixArrToCopy.shape = {ResPixArrToCopy.shape}')
            
            F = ResPixArrToCopy.shape[0]
            
            print(f'\nResPixArrToCopy now has {F} frames')
        
        
        # Get pixel shift between FromSliceNum in SrcIm and ToSliceNum in 
        # TrgIm in the Target image domain:
        PixShift = GetPixelShiftBetweenSlices(Image0=SrcIm, 
                                              SliceNum0=FromSliceNum, 
                                              Image1=TrgIm,
                                              SliceNum1=ToSliceNum,
                                              RefImage=TrgIm)
        
        print(f'\nPixShift = {PixShift}')
        
        # Shift the the in-plane elements in ResPixArrToCopy:
        ResPixArrToCopy = ShiftFrame(Frame=ResPixArrToCopy, PixShift=PixShift)
        
        unique = UniqueItems(Items=ResPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy',
              'after shifting the frame.')
        
        # Modify TrgCStoSliceInds to reflect the slice to which ResPixArrToCopy
        # relates to:
        TrgCStoSliceInds = deepcopy([ToSliceNum])
        
        
            
    if '3' in UseCase:
        # Convert ResPixArrToCopy to a list of ContourData by contour and a
        # list of points by contour: 
        #TrgCntDataByCnt, TrgPtsByCnt = PixArr2Contours(PixArr=ResPixArrToCopy)
        TrgPtsByCnt,\
        TrgCntDataByCnt = PixArr2PtsByContour(PixArr=ResPixArrToCopy,
                                              FrameToSliceInds=TrgCStoSliceInds,
                                              DicomDir=TrgDcmDir)
        
        """
        03/12/20: The z-coordinates are suspiciously all the same.  Should look
        into this.
        """
        
        print(f'\nAfter converting ResPixArrToCopy to ContourData and points',
              f'TrgCStoSliceInds = {TrgCStoSliceInds}.')
        
        print(f'\nlen(TrgCntDataByCnt) = {len(TrgCntDataByCnt)}')
        print(f'\nlen(TrgPtsByCnt) = {len(TrgPtsByCnt)}')
        
    
    
    if UseCase == '3a':
        if TrgRts:
            """
            Since this is a Direct copy, any existing contours in Target with 
            ROI label matching FromRoiLabel are to be preserved.
            """
            
            # Get all points (by contour) in the Target ROI containing the 
            # contour to be copied:
            TrgPtsInRoiByCnt,\
            TrgCStoSliceIndsInRoi = GetPtsInRoi(Rts=TrgRts, 
                                                DicomDir=TrgDcmDir,
                                                RoiLabel=FromRoiLabel)
            
            print('\nThe Target ROI containing the contour to be copied has',
                  f'{len(TrgPtsInRoiByCnt)} contours. \nTrgCStoSliceIndsInRoi',
                  f'= {TrgCStoSliceIndsInRoi}.')
            
            
            """
            The contour points that arose from ResPixArrToCopy (TrgPtsByCnt)  
            needs to be added to previously Direct-copied contour points 
            TrgPtsInRoiByCnt.
            """
            
            TrgCntDataByCnt,\
            TrgPtsByCnt,\
            TrgCStoSliceInds = AddCopiedPtsByCnt(OrigPtsByCnt=TrgPtsInRoiByCnt, 
                                                 OrigCStoSliceInds=TrgCStoSliceIndsInRoi, 
                                                 PtsToAddByCnt=TrgPtsByCnt, 
                                                 CStoSliceIndsToAddByCnt=TrgCStoSliceInds)
        else:
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
                                          RegImFilt=RegImFilt)
        
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
            
            unique = UniqueItems(Items=TxPixArrPreBinary, NonZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArrToCopy',
                  'prior to binary thresholding TxLabmapImToCopy')
        """ end """
        
        # Binary threshold TxLabmapImToCopy:
        TxLabmapImToCopy = BinaryThresholdImage(Im=TxLabmapImToCopy,
                                                Thresh=Thresh)
            
        # Convert TxLabmapImToCopy to a pixel array:
        TxPixArrToCopy,\
        TrgCStoSliceInds = Image2PixArr(LabmapIm=TxLabmapImToCopy)
        
        unique = UniqueItems(Items=TxPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in TxPixArrToCopy',
              'after binary thresholding TxLabmapImToCopy')
        
        F = TxPixArrToCopy.shape[0]
        
        print(f'\nThe segmentation (from the contour) on slice {FromSliceNum}',
              f'has been registered to {F} frames (contours):',
              f'{TrgCStoSliceInds}.')
        
        
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
        TrgPtsByCnt,\
        TrgCntDataByCnt = PixArr2PtsByContour(PixArr=TxPixArrToCopy,
                                              FrameToSliceInds=TrgCStoSliceInds,
                                              DicomDir=TrgDcmDir)
    
    
    
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
    
    if TrgRts:
        """ Use TrgRts to initialise the Target RTS. """
        TrgRts = InitialiseRts(RtsTemplate=TrgRts, RoiNum=FromRoiNum,
                               CStoSliceInds=TrgCStoSliceInds, 
                               DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    else:
        """ Use SrcRts to initialise the Target RTS. """
        TrgRts = InitialiseRts(RtsTemplate=SrcRts, RoiNum=FromRoiNum,
                               CStoSliceInds=TrgCStoSliceInds, 
                               DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    # Modify the tags in TrgRts:
    TrgRts = ModifyRts(Rts=TrgRts, CntDataByCnt=TrgCntDataByCnt, 
                       PtsByCnt=TrgPtsByCnt, CStoSliceInds=TrgCStoSliceInds,
                       DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n\nDone.  Took {Dtime} s to run.')
    
    
    return TrgRts






"""
******************************************************************************
******************************************************************************
COPY A SEGMENTATION
******************************************************************************
******************************************************************************
"""



def CopySeg(SrcSeg, FromSegLabel, SrcDcmDir, FromSliceNum, TrgDcmDir,
            TrgSeg=None, ToSliceNum=None, LogToConsole=False):
                                
    """
    Make a direct or relationship-preserving copy of a single segmentation in a
    given segment on a given slice to a new segment in an existing SEG, or 
    in a newly created SEG (if an existing SEG is not provided). 
    
    For direct copies the in-plane coordinates (currently the (x,y) coords) 
    will be preserved only, whilst for a relationship-preserving copy, all
    coordinates are spatially preserved.
    
    There are five possible Main Use Cases (indicated by the numbers 1 to 5)
    and two Sub-cases (indicated by the letters a or b) that apply for a subset
    of the Main Use Cases that are accommodated:
    
    1a. Direct copy of the Source segmentation in FromSegLabel on slice  
        FromSliceNum to a new segment for slice ToSliceNum in Target. 
        The Source and Target images are the same. 
    
    1b. Relationship-preserving copy of the Source segmentation in FromSegLabel 
        on slice FromSliceNum to a new segment for slice ToSliceNum in Target. 
        The Source and Target images are the same. 
        
    2a. Direct copy of the Source segmentation in FromSegLabel on slice 
        FromSliceNum to a new segment for slice ToSliceNum in Target. 
        The Source and Target images have the same FOR, IOP, IPP, and VS.
    
    2b. Relationship-preserving copy of the Source segmentation in FromSegLabel 
        on slice FromSliceNum to a new segment for slice ToSliceNum in Target. 
        The Source and Target images have the same FOR, IOP, IPP, and VS. 
    
    3a. Direct copy of the Source segmentation in FromSegLabel on slice 
        FromSliceNum to a new segment for slice ToSliceNum in Target. 
        The Source and Target images have the same FOR and IOP but have 
        different VS (they will likely also have different IPP).
        
    3b. Relationship-preserving copy of the Source segmentation in FromSegLabel 
        on slice FromSliceNum to a new segment for any number of slices in 
        Target.
        The Source and Target images have the same FOR and IOP but have 
        different VS (they will likely also have different IPP).
    
    4.  Relationship-preserving copy of the Source segmentation in FromSegLabel 
        on slice FromSliceNum to a new segment for any number of slices in 
        Target.
        The Source and Target images have the same FOR but have different IOP 
        and IPP, and possibly different VS.
       
    5.  Relationship-preserving copy of the Source segmentation in FromSegLabel 
        on slice FromSliceNum to a new segment for any number of slices in 
        Target.
        The Source and Target images have different FOR, IPP and IOP, and 
        possibly different VS.
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, VS = Voxel spacings
    
    
    
    Inputs:
    ------
    
    SrcSeg : Pydicom object
        Source DICOM SEG object.
                             
    FromSegLabel : string
        All or part of the Source segment label containing the segmentation to 
        be copied.
                     
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
               
    ToSliceNum : integer (default = None)
        Slice index within the Target DICOM stack where the segmentation is to
        be copied to (counting from 0).  This only applies for Direct copies,
        hence the default value None.
                        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    -------
        
    TrgSeg : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy of a segmentation to an existing Direct-copied 
        contour) Target SEG object.
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    
    import time
    from copy import deepcopy
    from DicomTools import IsSameModalities
    from DicomTools import GetRoiNum
    from SegTools import GetFrameFromPixArr
    from GeneralTools import PrintTitle
    from GeneralTools import UniqueItems
    from ImageTools import ImportImage
    from SegTools import InitialiseSeg
    from SegTools import ModifySeg
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Set the threshold level for binarising of the transformed labelmap (for
    # UseCase 5):
    Thresh = 0.1
    
    # Verify that FromSegLabel exists:
    #SrcSegLabels = GetRoiLabels(SrcSeg)
    FromSegNum = GetRoiNum(Roi=SrcSeg, Label=FromSegLabel)
    
        
    # Compare the modalities of the segments:
    if TrgSeg:
        if not IsSameModalities(SrcSeg, TrgSeg):
            msg = f"The Source ({SrcSeg.Modality}) and Target " \
                  + "({TrgSeg.Modality} modalities are different."
            
            raise Exception(msg)
            
        ToSegNum = GetRoiNum(Roi=TrgSeg, Label=FromSegLabel) # 04/12
    
    # The pixel array that only contains the frame to be copied:
    PixArrToCopy = GetFrameFromPixArr(Seg=SrcSeg, DicomDir=SrcDcmDir, 
                                      SegLabel=FromSegLabel, 
                                      SliceNum=FromSliceNum,
                                      LogToConsole=LogToConsole)
    
    
    # Determine which Use Case applies:
    UseCase = WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole)
    
    if True:#LogToConsole:
        PrintTitle(f'Case {UseCase} applies.')
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    if 'a' in UseCase:
        """ 
        For UseCases 1a/b, 2a/b or 3a, there will always be a single 
        PerFrameFunctionalGroupsSequence. 
        """
        
        #TrgPFFGStoSliceInds = [FromSliceNum] # moved to end
        
    
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
        
        print(f'\nPixShift = {PixShift}')
        
        
    if UseCase in ['1a', '2a']:
        # Shift the the in-plane elements in PixArrToCopy:
        PixArrToCopy = ShiftFrame(Frame=PixArrToCopy, 
                                  PixShift=PixShift)
        
        if TrgSeg:
            """
            Since this is a Direct copy, any existing segmentations in Target 
            with segment label matching FromRoiLabel are to be preserved.
            """
            
            # Get all frames in the Target segment containing the segmentation 
            # to be copied:
            TrgPixArrInSeg,\
            TrgPFFGStoSliceIndsInSeg = GetPixArrInSeg(Seg=TrgSeg, 
                                                      DicomDir=TrgDcmDir,
                                                      SegLabel=FromSegLabel)
            
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
    
    
    if UseCase in ['1b', '2b']:
        
        TrgPixArr = deepcopy(PixArrToCopy)
        
        TrgPFFGStoSliceInds = deepcopy([FromSliceNum]) # 04/12 verify?
    
        
        
    if UseCase in ['3a', '3b', '4', '5']:
        from ConversionTools import PixArr2Image
        from ConversionTools import Image2PixArr
        
        # Import the 3D images:
        if not 'a' in UseCase:
            SrcIm = ImportImage(SrcDcmDir)
            TrgIm = ImportImage(TrgDcmDir)
        
        # Convert PixArrToCopy to a labelmap image:
        LabmapImToCopy = PixArr2Image(PixArr=PixArrToCopy,
                                      FrameToSliceInds=[FromSliceNum],
                                      RefIm=SrcIm)
        
        
    if '3' in UseCase:
        from ImageTools import ResampleImage
        
        # Resample LabmapImToCopy to the Target image's grid:
        ResLabmapImToCopy = ResampleImage(Image=LabmapImToCopy, 
                                          RefImage=TrgIm, 
                                          Interpolation='NearestNeighbor')
        
        # Convert ResLabmapImToCopy back to a pixel array:
        ResPixArrToCopy,\
        PFFGStoSliceIndsToCopy = Image2PixArr(LabmapIm=ResLabmapImToCopy)
        
        
    if UseCase == '3a':
        """
        For a Direct copy, one segmentation is always copied to a single
        segmentation, so if after resampling there is more than one frame, the
        frames will need to be averaged.
        """
        from GeneralTools import MeanPixArr
        
        unique = UniqueItems(Items=ResPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy')
        
        F = ResPixArrToCopy.shape[0]
        
        print(f'\nThe segmentation on FromSliceNum = {FromSliceNum} has been',
              f'resampled to {F} frames: {PFFGStoSliceIndsToCopy}')
        
        if F > 1:
            print(f'\nResPixArrToCopy has {F} frames so the pixel arrays will',
                  'be averaged...')
            
            # Take pixel-by-pixel mean across all frames in ResPixArrToCopy:
            ResPixArrToCopy = MeanPixArr(ResPixArrToCopy)
            
            F = ResPixArrToCopy.shape[0]
            
            print(f'\nResPixArrToCopy now has {F} frames')
        else:
            print(f'\nResPixArrToCopy has F frames')
            
        
        # Shift the in-plane elements in ResPixArrToCopy:
        ResPixArrToCopy = ShiftFrame(Frame=ResPixArrToCopy, PixShift=PixShift)
        
        unique = UniqueItems(Items=ResPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResPixArrToCopy',
              'after shifting the frame.')
        
        
        if TrgSeg:
            """
            Since this is a Direct copy, any existing segmentations in Target
            with segment label matching FromRoiLabel are to be preserved.
            """
            
            # Get all frames in the Target segment containing the segmentation 
            # to be copied:
            TrgPixArrInSeg,\
            TrgPFFGStoSliceIndsInSeg = GetPixArrInSeg(Seg=TrgSeg, 
                                                      DicomDir=TrgDcmDir,
                                                      SegLabel=FromSegLabel)
            
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
            
    
    
    if UseCase == '3b':
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
                                          RegImFilt=RegImFilt)
        
        # Binary threshold TxLabmapImToCopy:
        TxLabmapImToCopy = BinaryThresholdImage(Im=TxLabmapImToCopy, 
                                                Thresh=Thresh)
        
        # Convert TxLabmapImToCopy to a pixel array:
        TrgPixArr,\
        TrgPFFGStoSliceInds = Image2PixArr(LabmapIm=TxLabmapImToCopy)
        
        unique = UniqueItems(Items=TrgPixArr, NonZero=False)
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
    
    if TrgSeg:
        """ Use TrgSeg to initialise the Target SEG. """
        TrgSeg = InitialiseSeg(SegTemplate=TrgSeg, SegNum=ToSegNum,
                               PFFGStoSliceInds=TrgPFFGStoSliceInds, 
                               DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    else:
        """ Use SrcSeg to initialise the Target SEG. """
        TrgSeg = InitialiseSeg(SegTemplate=SrcSeg, SegNum=FromSegNum,
                               PFFGStoSliceInds=TrgPFFGStoSliceInds, 
                               DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    # Modify the tags in TrgSeg:
    TrgSeg = ModifySeg(Seg=TrgSeg, PixArr=TrgPixArr,
                       PFFGStoSliceInds=TrgPFFGStoSliceInds, 
                       DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n\nDone.  Took {Dtime} s to run.')
    
    """ Had 1040 lines of code. Now 207. """
    
    return TrgSeg

    