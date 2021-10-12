# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:01:08 2020

@author: ctorti
"""



"""
Old functions from RoiCopyTools.py and others
"""








def CopyContourWithinSeries_OLD(SrcDcmDir, SrcRtsFpath, FromSliceNum, ToSliceNum, 
                            LogToConsole=False):
    """
    Copy a single contour on a given slice to another slice within the same
    DICOM series.
    
    
    Inputs:
    ------
    
    SrcDcmDir : string
    
        Directory containing the Source (= Target) DICOMs
                            
    SrcRtsFpath : string
    
        Full path of the Source (= Target) DICOM-RTSTRUCT file 
                             
    FromSliceNum : integer
    
        Slice index of the Source (= Target) DICOM stack corresponding to the 
        contour to be copied (counting from 0)
                             
    ToSliceNum : integer
    
        Slice index of the Target (= Source) DICOM stack where the contour will
        be copied to (counting from 0)
                            
    LogToConsole : boolean
    
        Denotes whether some intermediate results will be logged to the console
                            
    
    Output:
    ------
        
    TrgRts : Pydicom object
    
        Modified Target (= Source) ROI object 
    
    """
    
    # Import the Source RTSTRUCT ROI:
    SrcRts = dcmread(SrcRtsFpath)
    
    # Import the Source (= Target) DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir)
    
    # Get the Source (= Target) DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices, 
    # and ContourSequence-to-DICOM slice indices for the Source DICOMs:
    SrcCIStoSliceInds = GetCIStoSliceInds(SrcRts, SrcSOPuids)
    
    SrcCStoSliceInds = GetCStoSliceInds(SrcRts, SrcSOPuids)
    
    
    # Verify that a contour exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcCIStoSliceInds:
        print(f'\nA contour does not exist on slice {FromSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a contour:
    """
    Note:
        If slice number ToSliceNum already has a contour, it will be 
        over-written with the contour data of slice number FromSliceNum.
        If it doesn't have a contour, ContourImageSequence and 
        ContourSequence need to be extended by one.
    """
    if ToSliceNum in SrcCIStoSliceInds:
        # Use SrcRts as a template for TrgRts: 
        TrgRts = deepcopy(SrcRts)
        
        # Create a new list of CIStoSliceInds:
        TrgCIStoSliceInds = deepcopy(SrcCIStoSliceInds)
        
    else:
        # Increase the length of the sequences:
        TrgRts = AddToRtsSequences(SrcRts)
        
        # Create a new list of CIStoSliceInds by adding the index ToSliceNum 
        # (representing the contour to be added to ContourImageSequence and to
        # ContourSequence) to SrcCIStoSliceInds:
        TrgCIStoSliceInds = AppendItemToListAndSort(SrcCIStoSliceInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nSrcCIStoSliceInds =', SrcCIStoSliceInds)
        print('\nTrgCIStoSliceInds =', TrgCIStoSliceInds)
    
    # Modify select tags in TrgRts:
    TrgRts = ModifyRtsTagVals(SrcRts, TrgRts, FromSliceNum, ToSliceNum, 
                              SrcDicoms, SrcCIStoSliceInds, TrgCIStoSliceInds,
                              LogToConsole)
    
    
    
    # Generate a new SOP Instance UID:
    TrgRts.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    TrgRts.SeriesInstanceUID = generate_uid()
    
    
    return TrgRts





def DirectCopyContourAcrossSeries_OLD(SrcDcmDir, SrcRtsFpath, FromSliceNum, 
                                  TrgDcmDir, TrgRtsFpath, ToSliceNum, 
                                  LogToConsole=False):
    """
    Make a direct copy of a single contour on a given slice to another slice in 
    a different DICOM series.  A direct copy implies that the z-coordinates of
    the Target points (copy of Source points) do not relate to the 
    z-coordinates of the DICOM slice, since the voxels corresponding to the 
    Target ROI volume do not relate to the voxels of the Source ROI volume.
    
    
    Inputs:
    ------
        
    SrcDcmDir : string
    
        Directory containing the Source DICOMs
                            
    SrcRtsFpath : string
    
        Full path of the Source DICOM-RTSTRUCT file 
                             
    FromSliceNum : integer
    
        Slice index of the Source DICOM stack corresponding to the contour to 
        be copied (counting from 0)
                       
    TrgDcmDir : string
    
        Directory containing the Target DICOMs
                            
    TrgRtsFpath : string
    
        Full path of the Target DICOM-RTSTRUCT file 
                             
    ToSliceNum : integer
    
        Slice index of the Target DICOM stack where the contour will be copied 
        to (counting from 0)
                            
    LogToConsole : boolean
    
        Denotes whether some intermediate results will be logged to the console
         
                   
    Output:
    ------
        
    NewTrgRts : Pydicom object
    
        Modified Target ROI object 
    
    """
    
    # Import the Source and Target RTSTRUCT ROI:
    SrcRts = dcmread(SrcRtsFpath)
    TrgRts = dcmread(TrgRtsFpath)
    
    # Import the Target DICOMs:
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices, 
    # and ContourSequence-to-DICOM slice indices for Source and Target:
    SrcCIStoSliceInds = GetCIStoSliceInds(SrcRts, SrcSOPuids)
    
    SrcCStoSliceInds = GetCStoSliceInds(SrcRts, SrcSOPuids)
    
    TrgCIStoSliceInds = GetCIStoSliceInds(TrgRts, TrgSOPuids)
    
    TrgCStoSliceInds = GetCStoSliceInds(TrgRts, TrgSOPuids)
    
    
        
    
    # Verify that a contour exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcCIStoSliceInds:
        print(f'\nA contour does not exist on slice {FromSliceNum} in the',
              'Source DICOM series.')
        
        return
    
    # Verify that the slice that a contour is to be copied to exists (this 
    # check won't be necessary in the OHIF-Viewer but a bad input could occur
    # here):
    if not ToSliceNum in list(range(len(TrgDicoms))):
        print(f'\nThe Target DICOM series does not have slice {ToSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a contour:
    """
    Note:
        If slice number ToSliceNum already has a contour, it will be 
        over-written with the contour data of slice number FromSliceNum.
        If it doesn't have a contour, ContourImageSequence and 
        ContourSequence need to be extended by one.
    """
    if ToSliceNum in TrgCIStoSliceInds:
        # Use TrgRts as a template for NewTrgRts: 
        NewTrgRts = deepcopy(TrgRts)
        
        # Create a new list of NewTrgCIStoSliceInds:
        NewTrgCIStoSliceInds = deepcopy(TrgCIStoSliceInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgRts = AddToRtsSequences(TrgRts)
        
        # Create a new list of TrgCIStoSliceInds by adding the index ToSliceNum 
        # (representing the contour to be added to ContourImageSequence and to
        # ContourSequence) to TrgCIStoSliceInds:
        NewTrgCIStoSliceInds = AppendItemToListAndSort(TrgCIStoSliceInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nSrcCIStoSliceInds    =', SrcCIStoSliceInds)
        print('\nTrgCIStoSliceInds    =', TrgCIStoSliceInds)
        print('\nNewTrgCIStoSliceInds =', NewTrgCIStoSliceInds)
    
    # Modify select tags in TrgRts:
    NewTrgRts = ModifyRtsTagVals(TrgRts, NewTrgRts, FromSliceNum, 
                                 ToSliceNum, TrgDicoms, 
                                 TrgCIStoSliceInds, NewTrgCIStoSliceInds,
                                 LogToConsole)
    
    
    
    # Generate a new SOP Instance UID:
    NewTrgRts.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRts.SeriesInstanceUID = generate_uid()
    
    
    return NewTrgRts





def MappedCopyContourAcrossSeries_OLD(SrcDcmDir, SrcRoiFpath, FromSliceNum, 
                                 TrgDcmDir, TrgRoiFpath, LogToConsole):
    """
    I need to work with labelmaps!!!  See MappedCopySegmentAcrossSeries...
    
    Make a relationship-preserving copy of a single contour on a given slice to  
    a different DICOM series.  Since the relationship is preserved, the voxels 
    corresponding to the points copied to the Target ROI volume spatially
    relate to the voxels corresponding to the points copied from the Source ROI
    volume.
    
    There are four possible cases that will need to be accommodated*:
        
        2. The Source and Target images have the same FOR, IOP, IPP, PS and ST 
        
        3. The Source and Target images have the same FOR and IOP but have
           different PS and/or ST (if different ST they will likely have 
           different IPP)
        
        4. The Source and Target images have the same FOR but have different
           IOP and IPP, and possibly different PS and/or ST
           
        5. The Source and Target images have different FOR (hence different IPP
           and IOP), and possibly different ST and/or PS
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, ST = SliceThickness, PS = PixelSpacing
    
    * The numbering deliberately starts at 2 for consistency with the 5 overall
    cases.  The first case is covered by another function.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcRoiFpath  - (String) Full path of the Source DICOM RTSTRUCT/SEG file
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the contour/segment to be copied 
                       (counting from 0)
                       
        TrgDcmDir    - (String) Directory containing the Target DICOMs
                            
        TrgRoiFpath  - (String) Full path of the Target DICOM RTSTRUCT/SEG file 
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewTrgRoi  - Modified Target ROI object 
    
    """
    
    # Import the Source and Target RTSTRUCT/SEG ROI:
    SrcRoi = pydicom.dcmread(SrcRoiFpath)
    TrgRoi = pydicom.dcmread(TrgRoiFpath)
    
    # Get the modality of the ROI:
    SrcModality = SrcRoi.Modality
    TrgModality = TrgRoi.Modality
    
    # Both modalities should be the same.  Return if not:
    if SrcModality == TrgModality:
        Modality = SrcModality
        
    else:
        print(f'\nThe Source ({SrcModality}) and Target ({TrgModality})',
              'modalities are different.')
        
        return
    
    # Import the Source and Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir, SortMethod='slices', Debug=False)
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir, SortMethod='slices', Debug=False)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices and
    # ContourSequence-to-DICOM slice indices, or the 
    # ReferencedImageSequence-to-DICOM slice indices and 
    # PerFrameFunctionalGroupsSequence-to-DICOM slice indices for Source and 
    # Target:
    if Modality == 'RTSTRUCT':
        SrcCIStoSliceInds = GetCIStoSliceInds(SrcRoi, SrcSOPuids)
    
        SrcCStoSliceInds = GetCStoSliceInds(SrcRoi, SrcSOPuids)
        
        TrgCIStoSliceInds = GetCIStoSliceInds(TrgRoi, TrgSOPuids)
        
        TrgCStoSliceInds = GetCStoSliceInds(TrgRoi, TrgSOPuids)
        
        # Both indices should be the same.  Return if not:
        if not VerifyCISandCS(SrcCIStoSliceInds, SrcCStoSliceInds, LogToConsole):
            return
        
        if not VerifyCISandCS(TrgCIStoSliceInds, TrgCStoSliceInds, LogToConsole):
            return
        
        # Verify that a contour exists for FromSliceNum (this check won't be 
        # necessary in the OHIF-Viewer but a bad input could occur here):
        if not FromSliceNum in SrcCIStoSliceInds:
            print(f'\nA contour does not exist on slice {FromSliceNum} in the',
                  'Source DICOM series.')
            
            return
        
    if Modality == 'SEG':
        SrcRIStoSliceInds = GetRIStoSliceInds(SrcRoi, SrcSOPuids)
    
        SrcPFFGStoSliceInds = GetPFFGStoSliceInds(SrcRoi, SrcSOPuids)
        
        TrgRIStoSliceInds = GetRIStoSliceInds(TrgRoi, TrgSOPuids)
        
        TrgPFFGStoSliceInds = GetPFFGStoSliceInds(TrgRoi, TrgSOPuids)
        
        # Verify that a segment exists for FromSliceNum (this check won't be 
        # necessary in the OHIF-Viewer but a bad input could occur here):
        if not FromSliceNum in SrcPFFGStoSliceInds:
            print(f'\nA segment does not exist on slice {FromSliceNum} in the',
                  'Source DICOM series.')
            
            return
    
    
    
    # Get the image attributes:
    SrcSize, SrcSpacings, SrcST, SrcIPPs,\
    SrcDirs = GetImageAttributes(SrcDcmDir)
    
    TrgSize, TrgSpacings, TrgST, TrgIPPs,\
    TrgDirs = GetImageAttributes(TrgDcmDir)
    
    
    # Check which case follows:
    if SrcFORuid == TrgFORuid:
        if SrcDirs == TrgDirs:
            if SrcSpacings == TrgSpacings:
                Case = 2
                
            else:
                Case = 3
            
        else:
            Case = 4
    
    else:
        Case = 5
    
    
    
    
    
    
    
    # Check if ToSliceNum already has a contour:
    """
    Note:
        If slice number ToSliceNum already has a contour, it will be 
        over-written with the contour data of slice number FromSliceNum.
        If it doesn't have a contour, ContourImageSequence and 
        ContourSequence need to be extended by one.
    """
    if ToSliceNum in TrgCIStoSliceInds:
        # Use TrgRtsRoi as a template for NewTrgRtsRoi: 
        NewTrgRtsRoi = deepcopy(TrgRtsRoi)
        
        # Create a new list of NewTrgCIStoSliceInds:
        NewTrgCIStoSliceInds = deepcopy(TrgCIStoSliceInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgRtsRoi = AddToRtsSequences(TrgRtsRoi)
        
        # Create a new list of TrgCIStoSliceInds by adding the index ToSliceNum 
        # (representing the contour to be added to ContourImageSequence and to
        # ContourSequence) to TrgCIStoSliceInds:
        NewTrgCIStoSliceInds = AppendItemToListAndSort(TrgCIStoSliceInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nSrcCIStoSliceInds    =', SrcCIStoSliceInds)
        print('\nTrgCIStoSliceInds    =', TrgCIStoSliceInds)
        print('\nNewTrgCIStoSliceInds =', NewTrgCIStoSliceInds)
    
    # Modify select tags in TrgRtsRoi:
    NewTrgRtsRoi = ModifyRtsTagVals(TrgRtsRoi, NewTrgRtsRoi, FromSliceNum, 
                                    ToSliceNum, TrgDicoms, 
                                    TrgCIStoSliceInds, NewTrgCIStoSliceInds,
                                    LogToConsole)
    
    
    
    # Generate a new SOP Instance UID:
    NewTrgRtsRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRtsRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    return NewTrgRtsRoi













def CopyRts_OLD(SrcDcmDir, SrcRtsFpath, FromSliceNum, FromRoiLabel,
            TrgDcmDir, TrgRtsFpath=None, ToSliceNum=None, 
            LogToConsole=False):
                                
    """
    Modified on 30/11/20.
    
    25/11: This was copied from CopySeg as a template so much of the 
    code has not yet been modified...
    
    
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
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    SrcRoiFpath : string
        Full path of the Source DICOM RTS file.
                             
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour to
        be copied (counting from 0).
                       
    FromRoiLabel : string
        All or part of the Source ROI label containing the contour to be
        copied.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    TrgRtsFpath : string (default = None)
        Full path of the Target DICOM RTS file. Default value is None, e.g. 
        if there is no DICOM RTS for the Target Series.
        
    ToSliceNum : integer (default = None)
        Slice index within the Target DICOM stack where the contour is to
        be copied to (counting from 0).  This only applies for Direct copies,
        hence the default value None.
                        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    -------
        
    SrcImage : SimpleITK object
        Source 3D image.
        
    _SrcImage : SimpleITK object
        Source 3D case-dependent image. If Main use case:
             = 1 or 2: SrcImage is returned unchanged.
             = 3 or 4: SrcImage resampled to TrgImage (ResSrcImage is returned).
             = 5: SrcImage registered to TrgImage (RegSrcImage is returned).
        
    TrgImage : SimpleITK object
        Target 3D image.
        
    SrcLabmapIms : List of SimpleITK objects
        List of Source 3D labelmaps - one labelmap per segment.
    
    _SrcLabmapIms : List of SimpleITK objects
        List of Source 3D case-dependent image - one labelmap per segment. If 
        Main use case:
            = 1 or 2: SrcLabmapIms is returned unchanged.
            = 3 or 4: Each SrcLabmapIm in ScrLabmapIms is resampled to TrgImage 
            (ResSrcLabmapIms is returned).
            = 5: Each SrcLabmapIm in SrcLabmapIms is registered to TrgImage 
            (RegSrcLabmapIms is returned).
    
    SrcLabmapToCopyIm : SimpleITK object
        Source 3D labelmap of the segmentation in FromSegLabel to be copied.
        
    _SrcLabmapToCopyIm : SimpleITK object
        Source 3D case-dependent image - one labelmap per segment. If Main use 
        case:
            = 1 or 2: SrcLabmapToCopyIm is returned unchanged.
            = 3 or 4: SrcLabmapToCopyIm is resampled to TrgImage 
            (ResSrcLabmapToCopyIm is returned).
            = 5: SrcLabmapToCopyIm is registered to TrgImage 
            (RegSrcLabmapToCopyIm is returned).
    
    TrgLabmapIms : List of SimpleITK objects
        List of Target 3D labelmaps - one labelmap per segment, if a Target SEG
        exists.  If not a list with a single empty 3D image is returned.
    
    NewTrgLabmapIms : List of SimpleITK objects
        List of new Target 3D labelmaps - one labelmap per segment.  If a 
        Target exists, the first n labelmaps are the same labelmaps in 
        TrgLabmapIms, and the final item is the labelmap that was copied.  If 
        not, a list with a single 3D labelmap is returned - the labelmap that
        was copied.  The new labelmap (last item in NewTrgLabmapIms) is case-
        dependent.  If Main use case:
            = 1 or 2: The new labelmap is simply SrcLabmapToCopyIm unchanged.
            = 3 or 4: The new labelmap is SrcLabmapToCopyIm resampled to 
            TrgImage (ResSrcLabmapToCopyIm is returned).
            = 5: The new labelmap is SrcLabmapToCopyIm transformed to the
            Target domain (TxSrcLabmapToCopyIm is returned).
    
    SrcRts : Pydicom object
        Source RTS object.
        
    TrgRts : Pydicom object
        Original Target RTS object or a copy of SrcRts if TrgRtsFpath = ''.
        
    NewTrgRts : Pydicom object
        New/modified Target RTS object.
    
    SrcDicoms : List of Pydicom objects
        List of Source DICOM objects.
    
    TrgDicoms : List of Pydicom objects
        List of Target DICOM objects.
    """
    
    
    # Determine which Use Case applies:
    UseCase = WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole)
    
    if True:#LogToConsole:
        #print(f'\nCase {UseCase} applies.')
        title = f'Case {UseCase} applies.'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
    
    """ 
    **************************************************************************
    Import the list of DICOM objects, 3D images, RTSs and Target ContourData. 
    **************************************************************************
    """
    SrcDicoms = ImportDicoms(SrcDcmDir) # not needed but useful
    TrgDicoms = ImportDicoms(TrgDcmDir) # needed for ModifyRtsTagVals
    
    SrcImage = ImportImage(SrcDcmDir)
    TrgImage = ImportImage(TrgDcmDir)
    
    SrcRts = dcmread(SrcRtsFpath)
    
    if TrgRtsFpath: # if TrgRtsFpath is not empty
        TrgRts = dcmread(TrgRtsFpath)
        
        TrgPtsByRoi = GetPtsByRoi(RtsFpath=TrgRtsFpath, DicomDir=TrgDcmDir)
    else:
        TrgPtsByRoi = []
    
    
    # Compare the modalities of the ROIs:
    if TrgRtsFpath:
        if not IsSameModalities(SrcRts, TrgRts):
            msg = f"The Source ({SrcRts.Modality}) and Target " \
                  + "({TrgRts.Modality} modalities are different."
            
            raise Exception(msg)
    
    
    """ 
    **************************************************************************
    Get the number of ROIs, their labels and the ROI number containing the ROI
    to be copied. 
    **************************************************************************
    """
    SrcNumOfRois = len(SrcRts.ROIContourSequence)
    
    SrcRoiLabels = GetRoiLabels(SrcRts)
    
    if not FromRoiLabel in SrcRoiLabels:
        raise Exception("There is no Source ROI with label " + FromRoiLabel\
                        + ".")
        
    FromRoiNum = SrcRoiLabels.index(FromRoiLabel)
    
    if TrgRtsFpath:
        TrgNumOfRois = len(TrgRts.ROIContourSequence)
    
        TrgRoiLabels = GetRoiLabels(TrgRts)
        
        #if 'a' in UseCase:
        #    """ 
        #    This is a direct copy. Is this the first time a contour is being
        #    copied from FromRoiLabel to TrgRts? If so ToRoiNum = 0 (the copied
        #    contour will be the first and only contour in the first and only
        #    ROI).  If not, ToRoiNum will be ... will always be 0. """ 
        #    # Has 
        #    if FromRoiLabel in TrgRoiLabels:
                
    
    if True:#LogToConsole:
        if TrgRtsFpath:
            msg = f'\nThe Source RTS has {SrcNumOfRois} ROIs with ' \
                  + f'label(s) {SrcRoiLabels} \nand the Target RTS has ' \
                  + f'{TrgNumOfRois} ROIs with label(s) {TrgRoiLabels}'
        else:
            msg = f'\nThe Source RTS has {SrcNumOfRois} ROIs with ' \
                  + f'label(s) {SrcRoiLabels}'
            
        print(msg)
    
    
        
    
    """
    **************************************************************************
    Get the DICOM SOP UIDs, the ContourImageSequence-to-slice indices
    (CIStoSliceInds), the ContourSequence-to-slice indices (CStoSliceInds), 
    and the CStoSliceInds grouped by ROI (CStoSliceIndsByRoi).
    **************************************************************************
    """
    
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    #SrcCIStoSliceInds = GetCIStoSliceInds(Rts=SrcRts, SOPuids=SrcSOPuids)
    
    SrcCStoSliceInds = GetCStoSliceInds(Rts=SrcRts, SOPuids=SrcSOPuids)
    
    SrcNumOfConsByRoi = GetNumOfContoursByRoi(Rts=SrcRts)
    
    SrcCStoSliceIndsByRoi = GroupListByRoi(ListToGroup=SrcCStoSliceInds, 
                                           NumOfConsByRoi=SrcNumOfConsByRoi)
    
    
    if TrgRtsFpath:
        TrgCIStoSliceInds = GetCIStoSliceInds(Rts=TrgRts, SOPuids=TrgSOPuids)
    
        TrgCStoSliceInds = GetCStoSliceInds(Rts=TrgRts, SOPuids=TrgSOPuids)
        
        TrgNumOfConsByRoi = GetNumOfContoursByRoi(Rts=TrgRts)
        
        TrgCStoSliceIndsByRoi = GroupListByRoi(ListToGroup=TrgCStoSliceInds, 
                                               NumOfConsByRoi=TrgNumOfConsByRoi)
        
        if True:#LogToConsole:
            print('\nThe Source series has contours on slices    ',
                  f'{SrcCStoSliceIndsByRoi} \nand the Target series has',
                  f'contours on slices {TrgCStoSliceIndsByRoi}',
                  '\n(grouped by ROI) \n(***Note: Counting from 0.***)')
    
    else:
        TrgCIStoSliceInds = []
    
        TrgCStoSliceInds = []
        
        TrgNumOfConsByRoi = []
        
        TrgCStoSliceIndsByRoi = []
        
        if True:#LogToConsole:
            print('The Source series has contours on slices',
                  f'{SrcCStoSliceIndsByRoi} \n(grouped by ROI)',
                  '\n(***Note: Counting from 0.***)')
    
    
    """ 
    Verify that a contour exists for FromSliceNum (this check won't be 
    necessary in the OHIF-Viewer but a bad input could occur here).
    """
    if not FromSliceNum in SrcCStoSliceInds:
        msg = "The Source series does not have a contour on slice " \
              + f"{FromSliceNum}"
        
        raise Exception(msg)
        
        return
    
    
    
    
    """
    **************************************************************************
    Create a 3D labelmap (as a SimpleITK image) for the segmentation to be 
    copied.
    **************************************************************************
    """    
    
    # Convert the contour data to be copied to a 2D mask:
    SrcPixArrToCopy = GetMaskFromRts(RtsFpath=SrcRtsFpath, 
                                     FromRoiLabel=FromRoiLabel, 
                                     FromSliceNum=FromSliceNum, 
                                     DicomDir=SrcDcmDir, 
                                     RefImage=SrcImage)
    
    # Convert the pixel array (frame) to be to be copied into an SimpleITK
    # image:
    """ This isn't required for Use Cases 1a/b and 2a/b. """
    SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                     FrameToSliceInds=[FromSliceNum],
                                     RefIm=SrcImage)
    
    
    """
    **************************************************************************
    Create a list of 3D labelmaps (as SimpleITK images) for all segmentations
    in all segments.
    
    This isn't needed but is created since the function returns SrcLabmapIms.
    **************************************************************************
    """
    
    SrcLabmapIms = GetLabmapImsByRoi(RtsFpath=SrcRtsFpath, 
                                     DicomDir=SrcDcmDir, 
                                     CStoSliceIndsByRoi=SrcCStoSliceIndsByRoi, 
                                     RefIm=SrcImage)
    
    if TrgRtsFpath:
        TrgLabmapIms = GetLabmapImsByRoi(RtsFpath=TrgRtsFpath, 
                                     DicomDir=TrgDcmDir, 
                                     CStoSliceIndsByRoi=TrgCStoSliceIndsByRoi, 
                                     RefIm=TrgImage)
    else:
        """ 
        This isn't needed but create a list of a single empty image since the 
        function returns TrgLabmapIms 
        """
        # Initialise the Target labelmap:
        TrgLabmapIms = []
        TrgLabmapIms.append(InitialiseImage(TrgImage))
    
    
    
    
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    if 'a' in UseCase:
        """ 
        This is common for all Direct copy sub-use cases.
        """
        
        # The contour points to be copied:
        SrcPtsToCopy = GetPtsInContour(RtsFpath=SrcRtsFpath, 
                                       DicomDir=SrcDcmDir, 
                                       RoiNum=FromRoiNum, 
                                       SliceNum=FromSliceNum)
        
        
        
    if UseCase == '1a' or UseCase == '2a':
        """
        Main use case = 1 or 2:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            voxel spacings, and either are the same series (Main use case = 1) 
            or different series (Main use case = 2).
            
        Sub-case = a:
            Direct copy of the contour on slice FromSliceNum in ROI 
            FromRoiLabel for Source to a new ROI on slice ToSliceNum for Target.
        
        """
        
        #print(f'\nApplying Case 2a..\n')
        if UseCase == '1a':
            title = 'Applying Case 1a..'
        else:
            title = 'Applying Case 2a..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        SrcSize, SrcSpacings, SrcST,\
        SrcIPPs, SrcDirs = GetImageAttributes(DicomDir=SrcDcmDir)
        
        # Convert SrcPtsToCopy to indices:
        SrcIndsToCopy = ContourToIndices(Points=SrcPtsToCopy, 
                                         Origin=SrcIPPs[0], 
                                         Directions=SrcDirs, 
                                         Spacings=SrcSpacings)
        
        # Change the z-component (k) of the indices to ToSliceNum:
        SrcIndsToCopy = ChangeZinds(Indices=SrcIndsToCopy, NewZind=ToSliceNum)
        
        # Convert back to physical points:
        SrcPtsToCopy = Inds2Pts(Indices=SrcIndsToCopy, 
                                Origin=SrcIPPs[0], 
                                Directions=SrcDirs, 
                                Spacings=SrcSpacings)
        
        # The new Target ContourImageSequence-to-slice indices are the original
        # list with ToSliceNum appended to it and sorted if this image was not
        # already in ContourImageSequence:
        if not ToSliceNum in TrgCIStoSliceInds:
            NewTrgCIStoSliceInds = AppendItemToListAndSort(OrigList=TrgCIStoSliceInds, 
                                                           ItemToAppend=ToSliceNum)
        
        # The new Target ContourSequence-to-slice indices are the original list
        # with ToSliceNum appended to it:
        NewTrgCStoSliceInds = deepcopy(TrgCStoSliceInds)
        NewTrgCStoSliceInds.append(ToSliceNum)
        
        # and for ContourSequence-to-slice indices, append [ToSliceNum]:
        NewTrgCStoSliceIndsByRoi = deepcopy(TrgCStoSliceIndsByRoi)
        NewTrgCStoSliceIndsByRoi.append([ToSliceNum])
        
        
        
    
    if UseCase == '1b' or UseCase == '2b':
        """
        Main use case = 1 or 2:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            voxel spacings, and either are the same series (Main use case = 1) 
            or different series (Main use case = 2).
            
        Sub-case = b:
            Relationship-preserving copy of the contour on slice FromSliceNum 
            (specified) in ROI FromRoiLabel for Source to a new ROI on 
            FromSliceNum (since ToSliceNum = FromSliceNum) for Target.
        """
        
        #print('\nApplying Case 2b..\n')
        if UseCase == '2a':
            title = 'Applying Case 2a..'
        else:
            title = 'Applying Case 2b..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        if LogToConsole:
            #ToSliceNum = deepcopy(FromSliceNum)
            
            print(f'\nCopying slice {FromSliceNum} from Source RTS to slice',
                  #f'{ToSliceNum} in Target SEG..')
                  f'{FromSliceNum} in Target RTS..')
        
        # The new Target ContourSequence-to-slice indices is the original one
        # with [FromSliceNum] appended to it:
        NewTrgCStoSliceInds = deepcopy(TrgCStoSliceInds)
        NewTrgCStoSliceInds.append(FromSliceNum)
        
        NewTrgCStoSliceIndsByRoi = deepcopy(TrgCStoSliceIndsByRoi)
        NewTrgCStoSliceIndsByRoi.append([FromSliceNum])
        
        
    if '1' in UseCase or '2' in UseCase:
        """
        What follows is common to UseCase 1a, 1b, 2a and 2b.
        """
        
        # Change SrcPtsToCopy to ContourData format (i.e. flat list of 
        # coordinates):
        SrcContourDataToCopy = Pts2ContourData(Points=SrcPtsToCopy)
        
        # Add SrcContourDataToCopy to the TrgContourData:
        NewTrgPtsByRoi = deepcopy(TrgPtsByRoi)
        NewTrgPtsByRoi.append(SrcPtsToCopy)
        
        # The new Target ROI numbers:
        NewTrgRoiNums = []
        
        print('')
        
        for RoiNum in range(len(NewTrgCStoSliceIndsByRoi)):
            # The number of contours for this ROI number:
            Ncontours = len(NewTrgCStoSliceIndsByRoi[RoiNum])
            
            print(f'RoiNum = {RoiNum}, Ncontours = {Ncontours}')
            
            NewTrgRoiNums.extend([RoiNum]*Ncontours)
        
        
        if LogToConsole:
            print('\nTrgCStoSliceInds    =', TrgCStoSliceInds)
            print('\nNewTrgCStoSliceInds =', NewTrgCStoSliceInds)
            
        print(f'\nThere are {len(NewTrgCStoSliceInds)} contours in',
              'NewTrgContourData')
        print(f'NewTrgCStoSliceInds: {NewTrgCStoSliceInds}')
        print(f'NewTrgRoiNums: {NewTrgRoiNums}')
        
        """
        26/11: Continue modifying this
        """
        
        if not TrgRtsFpath:
            # Target does not have an existing RTS, so initialise a new RTS 
            # object:
            TrgRts = InitialiseRts(Dicoms=TrgDicoms, Image=TrgImage,
                                   RtsTemplate=SrcRts, 
                                   LogToConsole=LogToConsole)
        
        # Modify select tags in NewTrgRts:
        NewTrgRts = ModifyRtsTagVals(TrgRts, TrgDicoms, NewTrgPtsByRoi, 
                                     NewTrgRoiNums, NewTrgCIStoSliceInds, 
                                     NewTrgCStoSliceIndsByRoi, FromRoiNum, 
                                     FromSliceNum, SrcRts, LogToConsole)
        
        
        """ 30/11: Need to continue modifying below... """
        
        
        """ Note:
            
        The following variables:
            
        ResSrcImage
        SrcLabmapIms
        ResSrcLabmapIms
        SrcLabmapToCopyIm
        ResSrcLabmapToCopyIm
        NewTrgLabmapIms
        
        were not used for this case but this main function outputs:
        
        return SrcImage, ResSrcImage, TrgImage,\
               SrcLabmapIms, ResSrcLabmapIms,\
               SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               NewTrgSeg
               
        So create:
            
        SrcLabmapToCopyIm
        NewTrgLabmapIms
        
        even though they're not needed for this case, and return 
        SrcLabmapToCopyIm for ResSrcLabmapToCopyIm, and SrcLabmapIms for
        ResSrcLabmapIms.
        """
        
        # This is already done prior to spitting to different use cases:
        #SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
        #                                 FrameToSliceInds=[ToSliceNum], 
        #                                 RefIm=TrgImage)
        
        ResSrcLabmapIms = []
        
        
        ResSrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                            #FrameToSliceInds=F2SInds,
                                            FrameToSliceInds=[ToSliceNum], 
                                            RefIm=TrgImage)
        
        # Create NewTrgLabmapIms even though it's not needed:
        NewTrgLabmapIms = []
        
        if TrgSegFpath:
            # Append the existing labelmaps for Target:
            for TrgLabmapIm in TrgLabmapIms:
                NewTrgLabmapIms.append(TrgLabmapIm)
            
        # then append the Source labelmap to be copied:    
        NewTrgLabmapIms.append(SrcLabmapToCopyIm)
        
        
        return SrcImage, SrcImage, TrgImage,\
               SrcLabmapIms, SrcLabmapIms,\
               SrcLabmapToCopyIm, SrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               SrcSeg, TrgSeg, NewTrgSeg,\
               SrcDicoms, TrgDicoms
               
        
        """
        23/10:  Still need to verify that Case 2 provides the expected results.
        
        26/10:  Resampling will need to be performed if the origins are not the
        same... 19/11: Now WhichUseCase() checks the origins and if they are
        not the same it's considered UseCase 3.
        """
        
    
    
    
    if '3' in UseCase:
        """
        What follows is common to both UseCase 3a and 3b.
        """
        
        # Resample the Source image:
        """ This is not required but useful info """
        ResSrcImage = ResampleImage(Image=SrcImage, 
                                    RefImage=TrgImage, 
                                    Interpolation='Linear')
        
        if LogToConsole:                            
            print('Results after resampling Source image:\n')
            
            print('The Source and Target image size:\n',
                          f'   Original Source:  {SrcImage.GetSize()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSize()} \n',
                          f'   Target:           {TrgImage.GetSize()}')
                
            print('The Source and Target voxel spacings:\n',
                          f'   Original Source:  {SrcImage.GetSpacing()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSpacing()} \n',
                          f'   Target:           {TrgImage.GetSpacing()}')
                
            print('The Source and Target Origin:\n',
                          f'   Original Source:  {SrcImage.GetOrigin()} \n',
                          f'   Resampled Source: {ResSrcImage.GetOrigin()} \n',
                          f'   Target:           {TrgImage.GetOrigin()}')
                
            print('The Source and Target Direction:\n',
                          f'   Original Source:  {SrcImage.GetDirection()} \n',
                          f'   Resampled Source: {ResSrcImage.GetDirection()} \n',
                          f'   Target:           {TrgImage.GetDirection()}')
                                   
                                    
        # Resample the Source labelmap images:
        """ This is not required but useful info """
        ResSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            ResSrcLabmapIm = ResampleImage(Image=SrcLabmapIm, 
                                           RefImage=TrgImage, 
                                           Interpolation='NearestNeighbor')
            
            ResSrcLabmapIms.append(ResSrcLabmapIm)
        
        
        
        # Resample the Source labelmap to be copied:
        ResSrcLabmapToCopyIm = ResampleImage(Image=SrcLabmapToCopyIm, 
                                             RefImage=TrgImage, 
                                             Interpolation='NearestNeighbor')
        
        #nda = sitk.GetArrayFromImage(ResSrcLabmapToCopyIm)
        #unique = UniqueItems(Nda=nda, NonZero=False)
        #print(f'\nThere are {len(unique)} unique items in ResSrcLabmapToCopyIm')
        
        
    if UseCase == '3a':
        """
        Main use case = 3:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            but different voxel spacings.  
            
        Sub-case = a:
            Direct copy of the segmentation on slice FromSliceNum (specified) 
            in segment FromSegLabel for Source to a new segment on slice  
            ToSliceNum (specified) for Target.
        """
        
        #print(f'\nApplying Case 3a..\n')
        title = 'Applying Case 3a..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        # Convert the resampled Source labelmap-to-be-copied image to a Numpy
        # array:
        ResSrcPixArrToCopy,\
        F2SInds = Image2PixArr(LabmapIm=ResSrcLabmapToCopyIm)
        
        unique = UniqueItems(Nda=ResSrcPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResSrcPixArrToCopy')
        
        print(f'\nThe segmentation on FromSliceNum = {FromSliceNum} has been',
              f'resampled to {ResSrcPixArrToCopy.shape[0]} frames:',
              f'{F2SInds}')
        
        """
        For Direct copy 1 segment is always copied to 1 frame, so if after
        resampling there is more than one frame, either all but one of frames 
        will need to be ignored, or all frames will need to be averaged or a
        logical OR operation performed on all frames, etc.
        """
        if ResSrcPixArrToCopy.shape[0] > 1:
            print(f'\nResSrcPixArrToCopy has {ResSrcPixArrToCopy.shape[0]}',
                  'frames so the pixel arrays will be averaged...')
            
            # Take pixel-by-pixel mean across all frames in ResSrcPixArrToCopy:
            ResSrcPixArrToCopy = MeanPixArr(ResSrcPixArrToCopy)
            
            print(f'\nResSrcPixArrToCopy now has {ResSrcPixArrToCopy.shape[0]}',
                  'frames')
        else:
            print(f'\nResSrcPixArrToCopy has {ResSrcPixArrToCopy.shape[0]}',
                  'frames')
            
        
        # Shift the the in-plane elements in the resampled Source frame that is
        # to be copied:
        ResSrcPixArrToCopy = ShiftFrame(Frame=ResSrcPixArrToCopy, 
                                        PixShift=PixShift)
        
        print(f'\nAfter shifting the frame, ResSrcPixArrToCopy has',
              f'{ResSrcPixArrToCopy.shape[0]} frames')
        
        unique = UniqueItems(Nda=ResSrcPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResSrcPixArrToCopy',
              'after shifting the frame.')
        
        # Convert the shifted resampled pixel array back to a SimpleITK image:
        """
        For Direct copy 1 segment is always copied to 1 frame, hence
        FrameToSliceInds=[ToSliceNum].
        """
        ResSrcLabmapToCopyIm = PixArr2Image(PixArr=ResSrcPixArrToCopy,
                                            #FrameToSliceInds=F2SInds,
                                            FrameToSliceInds=[ToSliceNum], 
                                            RefIm=TrgImage)
        
        
        
    if UseCase == '3b':
        """
        Main use case = 3:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            but different voxel spacings.  
            
        Sub-case = b:
            Relationship-preserving copy of segmentation on slice FromSliceNum 
            (specified) in segment FromSegLabel for Source to a new segment 
            with any number of segmentations on any number of slices for Target. 
        """
        
        #print(f'\nApplying Case 3b..')
        title = 'Applying Case 3a..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        
    if '3' in UseCase:
        """
        The following is common to UseCase 3a and 3b.
        """
        
        # Create a list of new Target Labelmap images including the existing 
        # labelmaps...:
        NewTrgLabmapIms = []
        
        if TrgSegFpath:
            # Append the existing labelmaps for Target:
            for TrgLabmapIm in TrgLabmapIms:
                NewTrgLabmapIms.append(TrgLabmapIm)
            
        # ... and append the resampled Source labelmap to be copied:    
        NewTrgLabmapIms.append(ResSrcLabmapToCopyIm)
        
        print(f'\nlen(NewTrgLabmapIms) = {len(NewTrgLabmapIms)}')
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrBySeg = []
        NewTrgPFFGStoSliceIndsBySeg = []
        NewTrgPFFGStoSliceInds = []
        NewTrgSegNumsBySeg = []
        NewTrgSegNums = []
        
        print('')
        
        for SegNum in range(len(NewTrgLabmapIms)):
            PixArr, F2SInds = Image2PixArr(LabmapIm=NewTrgLabmapIms[SegNum])
            
            unique = UniqueItems(Nda=PixArr, NonZero=False)

            print(f'   SegNum = {SegNum}, F2SInds = {F2SInds}')
            print(f'   There are {len(unique)} unique items in PixArr')
            
            NewTrgPixArrBySeg.append(PixArr)
            
            NewTrgPFFGStoSliceIndsBySeg.append(F2SInds)
            
            NewTrgPFFGStoSliceInds.extend(F2SInds)
            
            NewTrgSegNumsBySeg.append([SegNum]*len(F2SInds))
            
            NewTrgSegNums.extend([SegNum]*len(F2SInds))
            
        
        # Combine all frames:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoSliceInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for SegNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoSliceIndsBySeg[SegNum])):
                NewTrgPixArr[n] = NewTrgPixArrBySeg[SegNum][i]
                
                n += 1
        
        
        print(f'\nThe new Target pixel array has shape {NewTrgPixArr.shape}')
        print(f'\nNewTrgPFFGStoSliceInds: {NewTrgPFFGStoSliceInds}')
        print(f'NewTrgPFFGStoSliceIndsBySeg: {NewTrgPFFGStoSliceIndsBySeg}')
        print(f'\nNewTrgSegNums: {NewTrgSegNums}')
        print(f'NewTrgSegNumsBySeg: {NewTrgSegNumsBySeg}')
        
        
        unique = UniqueItems(Nda=NewTrgPixArr, NonZero=False)
        #unique = UniqueItems(Nda=NewTrgPixArr, NonZero=True)
        
        #print(f'\nUnique items in NewTrgPixArr = {unique}')
        print(f'\nThere are {len(unique)} unique items in NewTrgPixArr')
        
        if not TrgSegFpath:
            # Target does not have an existing SEG, so initialise a new SEG 
            # object:
            TrgSeg = InitialiseSeg(TrgDicoms, NewTrgPFFGStoSliceInds,
                                   TrgImage, SrcSeg, LogToConsole)
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoSliceInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole)
        
        
        if LogToConsole:
            print('\n\n')
            title = 'Results after resampling Source LabelMap image and ' \
                    + 'Source-to-copy LabelMap image:'
            ul = '*' * len(title)
            
            print('\n' + title)
            print(ul)
            
            print('\n*The Source and Target labelmap image size:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetSize()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetSize()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetSize()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSize()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSize()}')
                
            print('\n*The Source and Target labelmap image voxel spacings:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetSpacing()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetSpacing()} \n',
                          f'   Source-to-copy            {SrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSpacing()}')
                
            print('\n*The Source and Target labelmap image Origin:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetOrigin()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetOrigin()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Target:                   {TrgLabmapIm.GetOrigin()}')
                
            print('\n*The Source and Target labelmap image Direction:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetDirection()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetDirection()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Target:                   {TrgLabmapIm.GetDirection()}')
        
        
        
        if LogToConsole:
            #SrcLabmapImMax = ImageMax(SrcLabmapIm)
            #SrcLabmapImMin = ImageMin(SrcLabmapIm)
            
            #ResSrcLabmapImMax = ImageMax(ResSrcLabmapIm)
            #ResSrcLabmapImMin = ImageMin(ResSrcLabmapIm)
            
            SrcLabmapToCopyImMax = ImageMax(SrcLabmapToCopyIm)
            SrcLabmapToCopyImMin = ImageMin(SrcLabmapToCopyIm)
            
            ResSrcLabmapToCopyImMax = ImageMax(ResSrcLabmapToCopyIm)
            ResSrcLabmapToCopyImMin = ImageMin(ResSrcLabmapToCopyIm)
            
            TrgLabmapImMax = ImageMax(TrgLabmapIm)
            TrgLabmapImMin = ImageMin(TrgLabmapIm)
            
            print('\n*The Source and Target labelmap image minimum:*\n',
                          #f'   Source:                   {SrcLabmapImMin} \n',
                          #f'   Resampled Source:         {ResSrcLabmapImMin} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMin} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMin} \n',
                          f'   Target:                   {TrgLabmapImMin}')
            
            print('\n*The Source and Target labelmap image maximum:*\n',
                          #f'   Source:                   {SrcLabmapImMax} \n',
                          #f'   Resampled Source:         {ResSrcLabmapImMax} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMax} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMax} \n',
                          f'   Target:                   {TrgLabmapImMax}')
            
        
        return SrcImage, ResSrcImage, TrgImage,\
               SrcLabmapIms, ResSrcLabmapIms,\
               SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               SrcSeg, TrgSeg, NewTrgSeg,\
               SrcDicoms, TrgDicoms
    
         
        
        
    if UseCase == '4':
        """
        Main use case = 4:
            Source and Target have the same FrameOfReferenceUID and origins, 
            but different IOP and voxel spacings.  
            
        Does UseCase 4 require sub-cases?...
        
        Don't have test data set to develop this yet.
        On 17/11 Simon said he would ask Jess for "data with different 
        angulations".
        """
        
        """ need to continue """
        
        #print(f'\nApplying Case 4..')
        title = 'Applying Case 4..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        
    if UseCase == '5':
        """
        Main use case = 5:
            Source and Target have the different FrameOfReferenceUID.  
            
        Does UseCase 5 require sub-cases?...
        """
        
        
        #print(f'\nApplying Case 5..')
        title = 'Applying Case 5..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        # Register the images using an affine transformation:
        RegSrcImage, RegImFilt = RegisterImages(FixIm=TrgImage, MovIm=SrcImage, 
                                                Tx='affine',
                                                LogToConsole=LogToConsole)
        
        # Transform SrcLabmapToCopyIm:
        TxSrcLabmapToCopyIm = TransformImage(Im=SrcLabmapToCopyIm, 
                                             RegImFilt=RegImFilt)
        
        # Binary threshold TxSrcLabmapToCopyIm:
        TxSrcLabmapToCopyIm = BinaryThresholdImage(Im=TxSrcLabmapToCopyIm, 
                                                   Threshold=0.5)
        
        # Transform all Source labelmap images:
        """ This is not required but useful info """
        TxSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            TxSrcLabmapIm = TransformImage(Im=SrcLabmapIm, 
                                           RegImFilt=RegImFilt)
            
            TxSrcLabmapIm = BinaryThresholdImage(Im=TxSrcLabmapIm, 
                                                 Threshold=0.5)
            
            TxSrcLabmapIms.append(TxSrcLabmapIm)
        
        
        """
        24/11: The labelmaps in TxSrcLabmapIms look odd.
        """
        
        # Create a list of new Target Labelmap images including the existing 
        # labelmaps...:
        NewTrgLabmapIms = []
        
        if TrgSegFpath:
            # Append the existing labelmaps for Target:
            for TrgLabmapIm in TrgLabmapIms:
                NewTrgLabmapIms.append(TrgLabmapIm)
            
        # ... and append the transformed Source labelmap to be copied:    
        NewTrgLabmapIms.append(TxSrcLabmapToCopyIm)
        
        print(f'\nlen(NewTrgLabmapIms) = {len(NewTrgLabmapIms)}')
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrBySeg = []
        NewTrgPFFGStoSliceIndsBySeg = []
        NewTrgPFFGStoSliceInds = []
        NewTrgSegNumsBySeg = []
        NewTrgSegNums = []
        
        print('')
        
        for SegNum in range(len(NewTrgLabmapIms)):
            PixArr, F2SInds = Image2PixArr(LabmapIm=NewTrgLabmapIms[SegNum])
            
            unique = UniqueItems(Nda=PixArr, NonZero=False)

            print(f'   SegNum = {SegNum}, F2SInds = {F2SInds}')
            print(f'   There are {len(unique)} unique items in PixArr')
            
            NewTrgPixArrBySeg.append(PixArr)
            
            NewTrgPFFGStoSliceIndsBySeg.append(F2SInds)
            
            NewTrgPFFGStoSliceInds.extend(F2SInds)
            
            NewTrgSegNumsBySeg.append([SegNum]*len(F2SInds))
            
            NewTrgSegNums.extend([SegNum]*len(F2SInds))
            
        
        # Combine all frames:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoSliceInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for SegNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoSliceIndsBySeg[SegNum])):
                NewTrgPixArr[n] = NewTrgPixArrBySeg[SegNum][i]
                
                n += 1
        
        
        print(f'\nThe new Target pixel array has shape {NewTrgPixArr.shape}')
        print(f'\nNewTrgPFFGStoSliceInds: {NewTrgPFFGStoSliceInds}')
        print(f'NewTrgPFFGStoSliceIndsBySeg: {NewTrgPFFGStoSliceIndsBySeg}')
        print(f'\nNewTrgSegNums: {NewTrgSegNums}')
        print(f'NewTrgSegNumsBySeg: {NewTrgSegNumsBySeg}')
        
        
        unique = UniqueItems(Nda=NewTrgPixArr, NonZero=False)
        #unique = UniqueItems(Nda=NewTrgPixArr, NonZero=True)
        
        #print(f'\nUnique items in NewTrgPixArr = {unique}')
        print(f'\nThere are {len(unique)} unique items in NewTrgPixArr')
        
        if not TrgSegFpath:
            # Target does not have an existing SEG, so initialise a new SEG 
            # object:
            TrgSeg = InitialiseSeg(TrgDicoms, NewTrgPFFGStoSliceInds,
                                   TrgImage, SrcSeg, LogToConsole)
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoSliceInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole)
        
        #NewTrgSeg = deepcopy(TrgSeg)
        
        return SrcImage, RegSrcImage, TrgImage,\
               SrcLabmapIms, TxSrcLabmapIms,\
               SrcLabmapToCopyIm, TxSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               SrcSeg, TrgSeg, NewTrgSeg,\
               SrcDicoms, TrgDicoms









def CopySegWithinSeries(SrcDcmDir, SrcSegFpath, FromSliceNum, ToSliceNum, 
                        LogToConsole=False):
    """
    17/11:  This hasn't been updated to account for changes to various 
    functions - namely ModifySegTagVals.
    
    
    Copy a single segmentation on a given slice to another slice within the
    same DICOM series.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcSegFpath  - (String) Full path of the Source DICOM SEG file 
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the segment to be copied (counting from 
                       0)
                             
        ToSliceNum   - (Integer) Slice index in the Target DICOM stack where
                       the segment will be copied to (counting from 0)
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewSrcSeg - (Pydicom Object) Modified Source SEG ROI object
    
    """
    
    # Import the Source SEG ROI:
    Seg = dcmread(SrcSegFpath)
    
    # Import the Source DICOMs:
    Dicoms = ImportDicoms(DicomDir=SrcDcmDir)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, and the
    # Per-frameFunctionalGroupsSequence-to-DICOM slice indices:
    RIStoDcmInds = GetRIStoDcmInds(Seg, SOPuids)
    
    PFFGStoDcmInds = GetPFFGStoDcmInds(Seg, SOPuids)
    
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in PFFGStoDcmInds:
        print(f'\nA segment does not exist on slice {FromSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a segment:
    """
    Note:
        If slice number ToSliceNum already has a segment, it will be 
        over-written with the segment data of slice number FromSliceNum.
        If it doesn't have a segment, ReferencedImageSequence and 
        PerFrameFunctionalGroupsSequence need to be extended by one.
    """
    if ToSliceNum in RIStoDcmInds:
        # Use Seg as a template for NewSeg: 
        NewSeg = deepcopy(Seg)
        
        # Create a new list of RIStoDcmInds: <-- not needed?
        #NewRIStoDcmInds = copy.deepcopy(RIStoDcmInds)
        
        # Create a new list of PFFGStoDcmInds:
        NewPFFGStoDcmInds = deepcopy(PFFGStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewSeg = AddToSegSequences(Seg)
        
        # Create a new list of RIStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to ReferencedImageSequence)
        # to RIStoDcmInds: <-- not needed?
        #NewRIStoDcmInds = AddItemToListAndSort(RIStoDcmInds, ToSliceNum)
        
        # Create a new list of PFFGStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to 
        # PerFrameFunctionalGroupsSequence) to PFFGStoDcmInds:
        NewPFFGStoDcmInds = AddItemToListAndSort(PFFGStoDcmInds, ToSliceNum)
    
    
    if LogToConsole:
        print('\nPFFGStoDcmInds    =', PFFGStoDcmInds)
        print('\nNewPFFGStoDcmInds =', NewPFFGStoDcmInds)
        
    # Modify select tags in NewSegRoi:
    NewSeg = ModifySegTagVals(SrcSeg=Seg, SrcDicoms=Dicoms, 
                              SrcPFFGStoDcmInds=PFFGStoDcmInds, 
                              FromSliceNum=FromSliceNum, 
                              TrgSeg=Seg, TrgDicoms=Dicoms, 
                              TrgPFFGStoDcmInds=PFFGStoDcmInds, 
                              ToSliceNum=ToSliceNum, 
                              NewTrgSeg=NewSeg, 
                              NewTrgPFFGStoDcmInds=NewPFFGStoDcmInds, 
                              LogToConsole=LogToConsole)
                                 
    
    
    
    # Generate a new SOP Instance UID:
    NewSeg.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewSeg.SeriesInstanceUID = generate_uid()
    
    
    return NewSeg









def DirectCopySegAcrossSeries_OLD(SrcDcmDir, SrcSegFpath, FromSliceNum, 
                              TrgDcmDir, TrgSegFpath, ToSliceNum, 
                              LogToConsole=False):
    """
    17/11:  This hasn't been updated to account for changes to various 
    functions - namely ModifySegTagVals.
    
    
    Copy a single segmentation on a given slice to another slice in a different 
    DICOM series.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcSegFpath  - (String) Full path of the Source DICOM SEG file 
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the segment to be copied (counting from 
                       0)
                       
        TrgDcmDir    - (String) Directory containing the Target DICOMs
                            
        TrgSegFpath  - (String) Full path of the Target DICOM SEG file 
                             
        ToSliceNum   - (Integer) Slice index of the Target DICOM stack where
                       the segment will be copied to (counting from 0)
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewTrgSeg    - (Pydicom Object) Modified Target SEG object 
    
    """
    
    # Import the Source and Target SEG ROI:
    SrcSeg = dcmread(SrcSegFpath)
    TrgSeg = dcmread(TrgSegFpath)
    
    # Import the Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir) # <-- not needed but a 
    # required input to ModifySegTagVals 
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, and the
    # Per-frameFunctionalGroupsSequence-to-DICOM slice indices for Source and 
    # Target:
    #SrcRIStoDcmInds = GetRIStoDcmInds(SrcSeg, SrcSOPuids)
    
    SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcSeg, SrcSOPuids)
    
    TrgRIStoDcmInds = GetRIStoDcmInds(TrgSeg, TrgSOPuids)
    
    TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgSeg, TrgSOPuids)
        
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcPFFGStoDcmInds:
        print(f'\nA segment does not exist on slice {FromSliceNum} in the',
              'Source DICOM series.')
        
        return
    
    # Verify that the slice that a contour is to be copied to exists (this 
    # check won't be necessary in the OHIF-Viewer but a bad input could occur
    # here):
    if not ToSliceNum in list(range(len(TrgDicoms))):
        print(f'\nThe Target DICOM series does not have slice {ToSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a segment:
    """
    Note:
        If slice number ToSliceNum already has a segment, it will be 
        over-written with the segment data of slice number FromSliceNum.
        If it doesn't have a segment, ReferencedImageSequence and 
        PerFrameFunctionalGroupsSequence need to be extended by one.
    """
    if ToSliceNum in TrgRIStoDcmInds:
        # Use TrgSegRoi as a template for NewTrgSeg:
        NewTrgSeg = deepcopy(TrgSeg)
        
        # Create a new list of TrgRIStoDcmInds: <-- not needed?
        #NewTrgRIStoDcmInds = copy.deepcopy(TrgRIStoDcmInds)
        
        # Create a new list of PFFGStoDcmInds:
        NewTrgPFFGStoDcmInds = deepcopy(TrgPFFGStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgSeg = AddToSegSequences(TrgSeg)
        
        # Create a new list of TrgRIStoDcmInds by adding the index ToSliceNum 
        # (representing the contour to be added to ReferencedImageSequence) to 
        # TrgRIStoDcmInds: <-- not needed?
        #NewTrgRIStoDcmInds = AddItemToListAndSort(TrgRIStoDcmInds, ToSliceNum)
        
        # Create a new list of TrgPFFGStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to 
        # PerFrameFunctionalGroupsSequence) to TrgPFFGStoDcmInds:
        NewTrgPFFGStoDcmInds = AddItemToListAndSort(TrgPFFGStoDcmInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nTrgPFFGStoDcmInds    =', TrgPFFGStoDcmInds)
        print('\nNewTrgPFFGStoDcmInds =', NewTrgPFFGStoDcmInds)
    
    # Modify select tags in NewTrgSeg:
    NewTrgSeg = ModifySegTagVals(SrcSeg, SrcDicoms, SrcPFFGStoDcmInds, 
                                 FromSliceNum, TrgSeg, TrgDicoms, 
                                 TrgPFFGStoDcmInds, ToSliceNum, NewTrgSeg, 
                                 NewTrgPFFGStoDcmInds, LogToConsole)
    
    
    # Generate a new SOP Instance UID:
    NewTrgSeg.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    
    return NewTrgSeg







def DirectCopySegAcrossSeries(SrcDcmDir, SrcSegFpath, FromSliceNum, FromSegLabel,
                              TrgDcmDir, TrgSegFpath, ToSliceNum, 
                              LogToConsole=False):
    """
    17/11:  This hasn't been updated to account for changes to various 
    functions - namely ModifySegTagVals.
    
    
    Copy a single segmentation on a given slice to another slice in a different 
    DICOM series.
    
    
    Inputs:
    ------
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    SrcSegFpath : string
        Full path of the Source DICOM SEG file.
                             
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
                       
    FromSegLabel : string
        All or part of the Source segment label containing the segmentation to 
        be copied.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    TrgSegFpath : string
        Full path of the Target DICOM SEG file. May also be empty string ('') 
        if there is no DICOM SEG for the Target Series.
        
    ToSliceNum : integer
        Slice index of the Target DICOM stack corresponding to the segmentation
        that is to be copied (counting from 0).
                        
    LogToConsole : boolean
        Denotes whether some intermediate results will be logged to the console
          
                  
    Output:
    ------
        
    NewTrgSeg : Pydicom object
        Modified Target ROI object 
    """
    
    # Import the Source and Target SEG ROI:
    SrcSeg = dcmread(SrcSegFpath)
    TrgSeg = dcmread(TrgSegFpath)
    
    # Import the Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir) # <-- not needed but a 
    # required input to ModifySegTagVals 
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, and the
    # Per-frameFunctionalGroupsSequence-to-DICOM slice indices for Source and 
    # Target:
    #SrcRIStoDcmInds = GetRIStoDcmInds(SrcSeg, SrcSOPuids)
    
    SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcSeg, SrcSOPuids)
    
    TrgRIStoDcmInds = GetRIStoDcmInds(TrgSeg, TrgSOPuids)
    
    TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgSeg, TrgSOPuids)
        
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcPFFGStoDcmInds:
        print(f'\nA segment does not exist on slice {FromSliceNum} in the',
              'Source DICOM series.')
        
        return
    
    # Verify that the slice that a contour is to be copied to exists (this 
    # check won't be necessary in the OHIF-Viewer but a bad input could occur
    # here):
    if not ToSliceNum in list(range(len(TrgDicoms))):
        print(f'\nThe Target DICOM series does not have slice {ToSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a segment:
    """
    Note:
        If slice number ToSliceNum already has a segment, it will be 
        over-written with the segment data of slice number FromSliceNum.
        If it doesn't have a segment, ReferencedImageSequence and 
        PerFrameFunctionalGroupsSequence need to be extended by one.
    """
    if ToSliceNum in TrgRIStoDcmInds:
        # Use TrgSegRoi as a template for NewTrgSeg:
        NewTrgSeg = deepcopy(TrgSeg)
        
        # Create a new list of TrgRIStoDcmInds: <-- not needed?
        #NewTrgRIStoDcmInds = copy.deepcopy(TrgRIStoDcmInds)
        
        # Create a new list of PFFGStoDcmInds:
        NewTrgPFFGStoDcmInds = deepcopy(TrgPFFGStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgSeg = AddToSegSequences(TrgSeg)
        
        # Create a new list of TrgRIStoDcmInds by adding the index ToSliceNum 
        # (representing the contour to be added to ReferencedImageSequence) to 
        # TrgRIStoDcmInds: <-- not needed?
        #NewTrgRIStoDcmInds = AddItemToListAndSort(TrgRIStoDcmInds, ToSliceNum)
        
        # Create a new list of TrgPFFGStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to 
        # PerFrameFunctionalGroupsSequence) to TrgPFFGStoDcmInds:
        NewTrgPFFGStoDcmInds = AddItemToListAndSort(TrgPFFGStoDcmInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nTrgPFFGStoDcmInds    =', TrgPFFGStoDcmInds)
        print('\nNewTrgPFFGStoDcmInds =', NewTrgPFFGStoDcmInds)
    
    # Modify select tags in NewTrgSeg:
    NewTrgSeg = ModifySegTagVals(SrcSeg, SrcDicoms, SrcPFFGStoDcmInds, 
                                 FromSliceNum, TrgSeg, TrgDicoms, 
                                 TrgPFFGStoDcmInds, ToSliceNum, NewTrgSeg, 
                                 NewTrgPFFGStoDcmInds, LogToConsole)
    
    
    # Generate a new SOP Instance UID:
    NewTrgSeg.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    
    return NewTrgSeg







def MappedCopySegmentAcrossSeries_OLD(SrcDcmDir, SrcSegFpath, FromSliceNum, 
                                  FromSegLabel,
                                  TrgDcmDir, TrgSegFpath, ToSegLabel, 
                                  LogToConsole=False):
                                
    """
    16/11: Removed ToSegLabel from inputs. Instead of copying to an existing
    segment, create a new segment.
    
    Make a relationship-preserving copy of a single segment on a given slice to  
    a different DICOM series.  Since the relationship is preserved, the voxels 
    corresponding to the voxels copied to the Target ROI volume spatially
    relate to the voxels corresponding to the voxels copied from the Source ROI
    volume.
    
    There are four possible cases that will need to be accommodated*:
        
        2b. The Source and Target images have the same FOR, IOP, IPP, PS and ST 
        
        3b. The Source and Target images have the same FOR and IOP but have
           different PS and/or ST (if different ST they will likely have 
           different IPP)
        
        4. The Source and Target images have the same FOR but have different
           IOP and IPP, and possibly different PS and/or ST
           
        5. The Source and Target images have different FOR (hence different IPP
           and IOP), and possibly different ST and/or PS
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, ST = SliceThickness, PS = PixelSpacing
    
    * The numbering deliberately starts at 2b for consistency with the 5 
    overall cases.  Cases 1, 2a and 3a (=Direct copies) are covered by other 
    functions (CopySegmentWithinSeries for Case 1, and 
    DirectCopySegmentAcrossSeries for Cases 2a and 3a).
    
    
    Inputs:
    ------
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    SrcSegFpath : string
        Full path of the Source DICOM SEG file.
                             
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the 
        contour/segment to be copied (counting from 0).
                       
    FromSegLabel : string
        All or part of the segment label of the ROI in the Source SEG 
        containing the segment to be copied.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    TrgSegFpath : string
        Full path of the Target DICOM SEG file. May also be empty string ('') 
        if there is no DICOM SEG for the Target Series.
                   
    ToSegLabel : string
        All or part of the segment label of the ROI in the Target SEG that the 
        segment is to be copied to.
                        
    LogToConsole : boolean
        Denotes whether some intermediate results will be logged to the console
          
                  
    Output:
        
    NewTrgSeg : Pydicom object
        Modified Target ROI object 
    
    """
    
    # Import the Source and Target SEG ROI:
    SrcSeg = dcmread(SrcSegFpath)
    if TrgSegFpath: # if TrgSegFpath is not empty
        TrgSeg = dcmread(TrgSegFpath)
    
    # Get the modality of the ROIs:
    SrcModality = SrcSeg.Modality
    if TrgSegFpath:
        TrgModality = TrgSeg.Modality
    
    # Both modalities should be SEG.  Return if not:
    if TrgSegFpath:
        if SrcModality != 'SEG' or TrgModality != 'SEG':
            print(f'\nThe Source ({SrcModality}) and Target ({TrgModality})',
                  'modalities are different.')
            
            return
    
    # Import the Source and Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir)
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target FrameOfReferenceUID:
    SrcFOR = SrcDicoms[0].FrameOfReferenceUID
    TrgFOR = TrgDicoms[0].FrameOfReferenceUID
    
    # Get the Image Attributes for Source and Target using SimpleITK:
    SrcSitkSize, SrcSitkSpacing, SrcSitkST,\
    SrcSitkIPP, SrcSitkDir = GetImageAttributes(DicomDir=SrcDcmDir,
                                                Package='sitk')
    
    TrgSitkSize, TrgSitkSpacing, TrgSitkST,\
    TrgSitkIPP, TrgSitkDir = GetImageAttributes(DicomDir=TrgDcmDir,
                                                Package='sitk')
    
    if LogToConsole:
        CompareSourceTargetImageAttributes(SrcDcmDir, TrgDcmDir)
    
    
    """
    I'm not sure if the items in DimensionIndexSequence,
    e.g. DimensionDescriptionLabel = 'ReferencedSegmentNumber'
                                   = 'ImagePositionPatient'
                                   = 'StackID'
                                   = 'InStackPositionNumber'
                                   
    are important.  While the label 'ReferenedIndexSequence' are common to SEGs
    both generated by the OHIF-Viewer and RTS-to-SEG conversions, it seems that
    the label 'ImagePositionPatient' is only found in OHIF-generated SEGs, and
    the labels 'StackID' and 'InStackPositionNumber' from RTS-to-SEG 
    conversions.
    """
    
    # Get the number of ROIs in Source and Target:
    SrcNumOfROIs = len(SrcSeg.SegmentSequence)
    if TrgSegFpath:
        TrgNumOfROIs = len(TrgSeg.SegmentSequence)
    
    # Get the labels of all ROIs in the Source and Target SEGs:
    SrcSegLabels = GetSegLabels(SrcSeg)
    #SrcNumOfROIs = len(SrcSegLabels)
    
    if TrgSegFpath:
        TrgSegLabels = GetSegLabels(TrgSeg)
        #TrgNumOfROIs = len(TrgSegLabels)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            print(f'\nThe Source SEG has {SrcNumOfROIs} ROIs with label(s)',
                  f'{SrcSegLabels} \nand the Target SEG has {TrgNumOfROIs} ROIs',
                  f'with label(s) {TrgSegLabels}')
        else:
            print(f'\nThe Source SEG has {SrcNumOfROIs} ROIs with label(s)',
                  f'{SrcSegLabels}')
    
    
    
    # Get the DimensionIndexValues (which relate the segments to the item in 
    # the Referenced Instance Sequence):
    """
    Note:
        The DimensionalIndexValues are integers that start from 1 (unlike the
        RIS or PFFGS which begin at 0)!
    """
    SrcDIVs = GetDIVs(SrcSeg)
    if TrgSegFpath:
        TrgDIVs = GetDIVs(TrgSeg)
    
    # Group the DIVs by ROI:
    """
    Note:
        The DimensionalIndexValues are integers that start from 1 (unlike the
        RIS or PFFGS which begin at 0)!
    """
    SrcDIVsByRoi = GroupListBySegment(ListToGroup=SrcDIVs, DIVs=SrcDIVs)
    if TrgSegFpath:
        TrgDIVsByRoi = GroupListBySegment(ListToGroup=TrgDIVs, DIVs=TrgDIVs)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            print('\nThe Source SEG DimensionIndexValues grouped by ROI are',
                  f'\n{SrcDIVsByRoi} \n\n and the Target SEG',
                  f'DimensionIndexValues grouped by ROI are \n{TrgDIVsByRoi}')
        else:
            print('\nThe Source SEG DimensionIndexValues grouped by ROI are',
                  f'\n{SrcDIVsByRoi}')
        
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, 
    # the same grouped by ROI,
    # the PerFrameFunctionalGroupsSequence-to-DICOM slice indices, 
    # the same grouped by ROI, and 
    # the frame numbers grouped by ROI for Source and Target:
    """
    Note:
        The RIS and PFFGS are integers that start from 0 (unlike the DIVs which
        begin at 1)!
    """
    SrcRIStoDcmInds = GetRIStoDcmInds(SrcSeg, SrcSOPuids)

    SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcSeg, SrcSOPuids)
    
    SrcFrameNums = list(range(len(SrcPFFGStoDcmInds)))
    
    # Group the above by ROIs:
    #SrcRIStoDcmIndsByRoi = GroupListBySegment(ListToGroup=SrcRIStoDcmInds, 
    #                                      DIVs=SrcDIVs) # <-- RIS has more elements than DIV
    SrcPFFGStoDcmIndsByRoi = GroupListBySegment(ListToGroup=SrcPFFGStoDcmInds, 
                                                DIVs=SrcDIVs)
    
    SrcFrameNumsByRoi = GroupListBySegment(ListToGroup=SrcFrameNums, 
                                           DIVs=SrcDIVs)
    
    if TrgSegFpath:
        TrgRIStoDcmInds = GetRIStoDcmInds(TrgSeg, TrgSOPuids)
        
        TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgSeg, TrgSOPuids)
        
        TrgFrameNums = list(range(len(TrgPFFGStoDcmInds)))
        
        # Group the above by ROIs:
        #TrgRIStoDcmIndsByRoi = GroupListBySegment(ListToGroup=TrgRIStoDcmInds, 
        #                                      DIVs=TrgDIVs) # <-- RIS has more elements than DIV
        TrgPFFGStoDcmIndsByRoi = GroupListBySegment(ListToGroup=TrgPFFGStoDcmInds, 
                                                    DIVs=TrgDIVs)
        
        TrgFrameNumsByRoi = GroupListBySegment(ListToGroup=TrgFrameNums, 
                                               DIVs=TrgDIVs)
        
        if True:#LogToConsole:
            print('\nThe Source series has segments on slices',
                  f'{SrcPFFGStoDcmIndsByRoi}, \nand the Target series on',
                  f'slices {TrgPFFGStoDcmIndsByRoi} \n(grouped by ROI)')
            
            #print('\nThe RIS-to-DICOM slice num for the Source series is',
            #      f'\n{SrcRIStoDcmInds}, \nand for the Target series is',
            #      f'\n{TrgRIStoDcmInds}')
    
    else:
        if True:#LogToConsole:
            print('The Source series has segments on slices',
                  f'{SrcPFFGStoDcmIndsByRoi} \n(grouped by ROI)')
            
            #print('\nThe RIS-to-DICOM slice num for the Source series is',
            #      f'\n{SrcRIStoDcmInds}')
            
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcPFFGStoDcmInds:
        print('\nThe Source series does not have any segments on slice',
              f'{FromSliceNum}')
        
        return
    
    
    # The 3D pixel array:
    #SrcSegPixArr = SrcSeg.pixel_array
    #TrgSegPixArr = TrgSeg.pixel_array
    
    if False:
        # Create 3D labelmaps (= pixel arrays in same frame position as 
        # corresponding DICOM slices, with 2D zero masks in between):
        SrcLabmap = PixArr2Labmap(PixelArray=SrcSeg.pixel_array,
                                  NumOfSlices=SrcSitkSize[2],
                                  FrameToSliceInds=SrcPFFGStoDcmInds)
        
        if TrgSegFpath:
            TrgLabmap = PixArr2Labmap(PixelArray=TrgSeg.pixel_array,
                                      NumOfSlices=TrgSitkSize[2],
                                      FrameToSliceInds=TrgPFFGStoDcmInds)
        else:
            # Initialise the Target LabelMap:
            TrgLabmap = np.zeros((TrgSitkSize[2], TrgSitkSize[1], TrgSitkSize[0]), 
                                 dtype='uint')
        
        """ More useful to convert Pixel Array to sitk image.. """
    
    # Read in the Source and Target SimpleITK 3D images:
    SrcImage = ImportImage(SrcDcmDir)
    TrgImage = ImportImage(TrgDcmDir)
    
    if LogToConsole:
        title = 'Creating 3D labelmaps for *Source*...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    # Create SimpleITK images of the 3D labelmaps (= pixel arrays in same frame 
    # position as corresponding DICOM slices, with 2D zero masks in between)
    # for each ROI:
    SrcLabmapIms = []
    
    for RoiNum in range(SrcNumOfROIs):
        # The starting (= a) and ending (= b) frame numbers:
        a = SrcFrameNumsByRoi[RoiNum][0]
        b = SrcFrameNumsByRoi[RoiNum][-1]
        
        # The (partial) pixel array for this ROI:
        PixArr = SrcSeg.pixel_array[a:b+1]
        
        Inds = SrcPFFGStoDcmIndsByRoi[RoiNum]
        
        if LogToConsole:
            print(f'len(PixArr) = {len(PixArr)}')
            print(f'len(Inds) = {len(Inds)}')
        
        SrcLabmapIm = PixArr2Image(PixelArray=PixArr,
                                   Image=SrcImage,
                                   FrameToSliceInds=Inds)
        
        SrcLabmapIms.append(SrcLabmapIm)
        
    
    if LogToConsole:
        title = 'Creating 3D labelmaps for *Target*...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    if TrgSegFpath:
        TrgLabmapIms = []
        
        for RoiNum in range(TrgNumOfROIs):
            # The starting (= a) and ending (= b) frame numbers:
            a = TrgFrameNumsByRoi[RoiNum][0]
            b = TrgFrameNumsByRoi[RoiNum][-1]
            
            # The (partial) pixel array for this ROI:
            PixArr = TrgSeg.pixel_array[a:b+1]
            
            Inds = TrgPFFGStoDcmIndsByRoi[RoiNum]
            
            TrgLabmapIm = PixArr2Image(PixelArray=PixArr,
                                       Image=TrgImage,
                                       FrameToSliceInds=Inds)
            
            TrgLabmapIms.append(TrgLabmapIm)
            
        
    else:
        # Initialise the Target LabelMap:
        TrgLabmap = np.zeros((TrgSitkSize[2], TrgSitkSize[1], TrgSitkSize[0]), 
                             dtype='uint')
        
        TrgLabmapIm = sitk.GetImageFromArray(TrgLabmap)
        
        # Set the Origin, Spacing and Direction of TrgLabmapIm to that of 
        # TrgImage:
        TrgLabmapIm.SetOrigin(TrgImage.GetOrigin())
        TrgLabmapIm.SetSpacing(TrgImage.GetSpacing())
        TrgLabmapIm.SetDirection(TrgImage.GetDirection())
        
        # For consistency with other IF case:
        TrgLabmapIms = []
        TrgLabmapIms.append(TrgLabmapIm)
        
        
    """ 
    SrcLabmapIm is the 3D labelmap of all segments within SrcSeg.
    
    But when looking to copy a single segment from Source to Target, I'll need
    to know how that single segment is resampled.  So create a new labelmap
    containing the segment to be copied only.
    """
    
    if LogToConsole:
        title = 'Creating 3D labelmap for *Source* containing segment to copy only...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    # Create SimpleITK image of the 3D labelmap for the segment to be copied 
    # only:
    """
    There may be more than one frame number in PixelArray that corresponds to 
    FromSliceNum...
    # The frame number in PixelArray that corresponds to FromSliceNum:
    FromFrameNum = SrcPFFGStoDcmInds.index(FromSliceNum)
    """
    
    # Determine which Source ROI contains the segment to be copied:
    """
    What I call RoiNum is equivalent to Segment Number in the Segment Sequence.
    """
    FromRoiNum = [i for i, s in enumerate(SrcSegLabels) if FromSegLabel in s]
    FromRoiNum = FromRoiNum[0]
    
    if TrgSegFpath:
        # Determine which Target ROI the Source segment is to be copied to:
        ToRoiNum = [i for i, s in enumerate(TrgSegLabels) if ToSegLabel in s]
        ToRoiNum = ToRoiNum[0]
    
        print(f'\nFromRoiNum = {FromRoiNum} \nToRoiNum = {ToRoiNum}')
        
    else:
        print(f'\nFromRoiNum = {FromRoiNum}')
    
    # Determine which frame number in PixelArray corresponds to FromSliceNum:
    """
    FromFrameNum = SrcPFFGStoDcmIndsByRoi[FromRoiNum].index(FromSliceNum)  
    
    The problem with the above is that it's the frame number within the subset
    of frames pertaining to the ROI given by FromRoiNum - it doesn't relate to 
    the frame number in PixelArray (which contains all frames for all ROIs).
    So need something a bit more clever...
    """
    n = 0 # initialise frame number counter
    
    for RoiNum in range(SrcNumOfROIs):
        if RoiNum == FromRoiNum:
            # This ROI contains the frame to be copied, whose position is given
            # by the index of FromSliceNum is in SrcPFFGStoDcmIndsByRoi[RoiNum]
            # plus any frames that preceeded it (i.e. the frame counter n):
            FromFrameNum = n + SrcPFFGStoDcmIndsByRoi[RoiNum].index(FromSliceNum)
            
        else:
            # Add to the frame counter the number of frames in this ROI:
            n += len(SrcPFFGStoDcmIndsByRoi[RoiNum])
    
    
    print(f'\nFromSliceNum = {FromSliceNum} relates to frame number',
          f'{FromFrameNum} in PixelArray')
    
    # Determine which frame number(s) in PixelArray correspond to FromSliceNum: 
    #FromFrameNums = []
    #
    #for RoiNum in range(SrcNumOfROIs):
    #    if FromSliceNum in SrcPFFGStoDcmIndsByRoi[RoiNum]:
    #        ind = SrcPFFGStoDcmIndsByRoi[RoiNum].index(FromSliceNum)
    #        
    #        FromFrameNums.append(ind)
    #     
    #print(f'\nFromSliceNum = {FromSliceNum} relates to frame number(s)',
    #      f'{FromFrameNums} in PixelArray')
    
    # Knowing which ROI contains the frame to be copied and the frame number
    # defines the pixel array to be copied:
    FromPixelArr = SrcSeg.pixel_array[FromFrameNum]
    
    # Modify the pixel array so that all but the frame to be copied has zeros:
    """
    SrcPixArrToCopy = PixArrForOneFrame(PixelArray=FromPixelArr, 
                                        FrameNum=FromFrameNum)
    
    I either need to pass the entire PixelArray and specify the FrameNum to be
    copied (as was my original intention with the function PixArrForOneFrame),
    or pass the frame to be copied (FromPixelArr) and modify PixArrForOneFrame
    accordingly!
    """
    SrcPixArrToCopy = PixArrForOneFrame(PixelArray=SrcSeg.pixel_array, 
                                        FrameNum=FromFrameNum)
    
    if LogToConsole:
        print(f'\nSrcPixArrToCopy.max() = {SrcPixArrToCopy.max()}')
    
    SrcLabmapToCopyIm = PixArr2Image(PixelArray=SrcPixArrToCopy,
                                     Image=SrcImage,
                                     FrameToSliceInds=[FromSliceNum])
    
    if LogToConsole:
        print(f'\nMax value of SrcLabmapToCopyIm = {ImageMax(SrcLabmapToCopyIm)}')
    
    
    # Check which case follows:
    if SrcFOR == TrgFOR:
        if LogToConsole:
            print('\nThe Source and Target FrameOfReferenceUID are the same')
        
        # Compare the Source and Target directions. Get the maximum value of 
        # their absolute differences:
        AbsDiffs = abs(np.array(SrcSitkDir) - np.array(TrgSitkDir))
        MaxAbsDiff = max(AbsDiffs)
        
        # Consider the directions different if any of the vector differences 
        # are greater than epsilon:
        epsilon = 1e-06
        
        #if SrcDirs == TrgDirs: 
        """ The line above will result in insignificant differences (e.g. due
        to quantisation/pixelation) resulting in False, so instead consider
        MaxAbsDiff and epsilon. 
        """
        if MaxAbsDiff < epsilon:
            if LogToConsole:
                print('The Source and Target directions are the same')
            
            if SrcSitkSpacing == TrgSitkSpacing:
                if LogToConsole:
                    print('The Source and Target voxel sizes are the same')
                
                CaseNum = '2b'
                
            else:
                if LogToConsole:
                    print('The Source and Target voxel sizes are different:\n',
                          f'   Source: {SrcSitkSpacing} \n',
                          f'   Target: {TrgSitkSpacing}')
                    
                CaseNum = '3b'
            
        else:
            if LogToConsole:
                print('The Source and Target directions are different:\n',
                      f'   Source: [{SrcSitkDir[0]}, {SrcSitkDir[1]}, {SrcSitkDir[2]},',
                      f'            {SrcSitkDir[3]}, {SrcSitkDir[4]}, {SrcSitkDir[5]},',
                      f'            {SrcSitkDir[6]}, {SrcSitkDir[7]}, {SrcSitkDir[8]}]\n',
                      f'   Target: [{TrgSitkDir[0]}, {TrgSitkDir[1]}, {TrgSitkDir[2]},',
                      f'            {TrgSitkDir[3]}, {TrgSitkDir[4]}, {TrgSitkDir[5]},',
                      f'            {TrgSitkDir[6]}, {TrgSitkDir[7]}, {TrgSitkDir[8]}]')
                
                print(f'AbsDiffs = {AbsDiffs}')
            
            
            
            CaseNum = '4'
    
    else:
        if LogToConsole:
            print('The Source and Target FrameOfReferenceUID are different')
        
        CaseNum = '5'
        
    
    if LogToConsole:    
        print('The Source and Target origins:\n',
                          f'   Source: {SrcSitkIPP[0]} \n',
                          f'   Target: {TrgSitkIPP[0]}')
        
        print('The Source and Target image dims:\n',
                          f'   Source: {SrcSitkSize} \n',
                          f'   Target: {TrgSitkSize}')
    
    
    
    if LogToConsole:
        print(f'\nCase {CaseNum} applies.')
    
    
    if CaseNum == '2b':
        """
        The origins for Source and Target are the same since they have the same
        FrameOfReferenceUID.  And since the PixelSpacing and SliceThickness are
        also the same, ToSliceNum = FromSliceNum, and the relationships will be
        preserved.
        
        26/10: It's seems that having the same FOR does not guarantee the same
        origin.  Series 9 and 6 for Subject 011, MR4 have different origins!
        """
        
        print('\nApplying Case 2b..')
    
        ToSliceNum = deepcopy(FromSliceNum)
        
        if LogToConsole:
            print(f'\nCopying slice {FromSliceNum} from Source SEG to slice',
                  f'{ToSliceNum} in Target SEG..')
    
        # Check if ToSliceNum already has a segment:
        """
        Note:
            If slice number ToSliceNum already has a segment, it will be 
            over-written with the segment data of slice number FromSliceNum.
            If it doesn't have a segment, ReferencedImageSequence and 
            PerFrameFunctionalGroupsSequence need to be extended by one.
        """
        if ToSliceNum in TrgRIStoDcmInds:
            # Use TrgSeg as a template for NewTrgSegRoi:
            NewTrgSeg = deepcopy(TrgSeg)
            
            # Create a new list of TrgRIStoDcmInds: <-- not needed?
            #NewTrgRIStoDcmInds = deepcopy(TrgRIStoDcmInds)
            
            # Create a new list of PFFGStoDcmInds:
            NewTrgPFFGStoDcmInds = deepcopy(TrgPFFGStoDcmInds)
            
        else:
            # Increase the length of the sequences:
            NewTrgSeg = AddToSegSequences(TrgSeg)
            
            # Create a new list of TrgRIStoDcmInds by adding the index  
            # ToSliceNum (representing the contour to be added to 
            # ReferencedImageSequence) to TrgRIStoDcmInds: <-- not needed?
            #NewTrgRIStoDcmInds = AddItemToListAndSort(TrgRIStoDcmInds, 
            #                                          ToSliceNum)
            
            # Create a new list of TrgPFFGStoDcmInds by adding the index 
            # ToSliceNum (representing the segment to be added to 
            # PerFrameFunctionalGroupsSequence) to TrgPFFGStoDcmInds:
            NewTrgPFFGStoDcmInds = AddItemToListAndSort(TrgPFFGStoDcmInds, 
                                                        ToSliceNum)
            
        
        if LogToConsole:
            print('\nTrgPFFGStoDcmInds    =', TrgPFFGStoDcmInds)
            print('\nNewTrgPFFGStoDcmInds =', NewTrgPFFGStoDcmInds)
        
        # Modify select tags in NewTrgSeg:
        #NewTrgSeg = ModifySegTagVals(TrgSeg, NewTrgSeg, FromSliceNum, 
        #                             ToSliceNum, TrgDicoms, 
        #                             #TrgPFFGStoDcmInds, NewTrgPFFGStoDcmInds,
        #                             SrcPFFGStoDcmInds, NewTrgPFFGStoDcmInds,
        #                             LogToConsole)
        
        NewTrgSeg = ModifySegTagVals(SrcSeg, SrcDicoms, FromSliceNum, 
                                     TrgSeg, NewTrgSeg, TrgDicoms, ToSliceNum, 
                                     SrcPFFGStoDcmInds, 
                                     TrgPFFGStoDcmInds, NewTrgPFFGStoDcmInds, 
                                     LogToConsole)
        
        
        
        """
        23/10:  Still need to verify that Case 2 provides the expected results.
        
        26/10:  Resampling will need to be performed if the origins are not the
        same...
        """
        
        
    if CaseNum == '3b':
        print(f'\nApplying Case 3b..')
        
        # Read in the Source and Target SimpleITK 3D images:
        #SrcImage = ImportImage(SrcDcmDir)
        #TrgImage = ImportImage(TrgDcmDir)
        
        #SrcSitkOrigin = SrcSitkIm.GetOrigin() 
        #TrgSitkOrigin = SrcSitkIm.GetOrigin() 
        #
        #SrcSitkDirs = SrcSitkIm.GetDirection() 
        #
        #SrcSitkSpacings = SrcSitkIm.GetSpacing() 
        #
        #SrcSitkDims = SrcSitkIm.GetSize() 
        
        
        # Get the Image Attributes for Source and Target using SimpleITK:
        #SrcSitkSize, SrcSitkSpacing,\
        #SrcSitkIPPs, SrcSitkDirs = GetImageAttributes(DicomDir=SrcDcmDir,
        #                                              Package='sitk')
        #
        #TrgSitkSize, TrgSitkSpacing,\
        #TrgSitkIPPs, TrgSitkDirs = GetImageAttributes(DicomDir=TrgDcmDir,
        #                                              Package='sitk')
        
        
        """
        Although it's the Source Labelmap that needs to be resampled, first
        prove that I can resample the image.
        """
        
        
        # Resample the Source image:
        ResSrcImage = ResampleImage(Image=SrcImage, RefImage=TrgImage, 
                                    Interpolation='Linear')
        
        if LogToConsole:                            
            print('Results after resampling Source image:\n')
            
            print('The Source and Target image size:\n',
                          f'   Original Source:  {SrcImage.GetSize()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSize()} \n',
                          f'   Target:           {TrgImage.GetSize()}')
                
            print('The Source and Target voxel spacings:\n',
                          f'   Original Source:  {SrcImage.GetSpacing()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSpacing()} \n',
                          f'   Target:           {TrgImage.GetSpacing()}')
                
            print('The Source and Target Origin:\n',
                          f'   Original Source:  {SrcImage.GetOrigin()} \n',
                          f'   Resampled Source: {ResSrcImage.GetOrigin()} \n',
                          f'   Target:           {TrgImage.GetOrigin()}')
                
            print('The Source and Target Direction:\n',
                          f'   Original Source:  {SrcImage.GetDirection()} \n',
                          f'   Resampled Source: {ResSrcImage.GetDirection()} \n',
                          f'   Target:           {TrgImage.GetDirection()}')
                                   
                                    
        # Resample the Source labelmap images:
        ResSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            ResSrcLabmapIm = ResampleImage(Image=SrcLabmapIm, RefImage=TrgImage, 
                                           Interpolation='NearestNeighbor')
            
            ResSrcLabmapIms.append(ResSrcLabmapIm)
            
        
        """
        Now need to work with the labelmap only containing the segment to be 
        copied.
        """
        # Resample the Source labelmap:
        ResSrcLabmapToCopyIm = ResampleImage(Image=SrcLabmapToCopyIm, 
                                             RefImage=TrgImage, 
                                             Interpolation='NearestNeighbor')
        
        # Perform a pixel-wise OR of ResSrcLabmapToCopyIm and the labelmap of
        # Target that the segment is to be copied to:
        ORLabmapIm = OrImages(ResSrcLabmapToCopyIm, TrgLabmapIms[ToRoiNum])
        
        """
        I can't tell if OR has worked since my ROIs (tumour and brain) overlap!
        So rather than performing OR operation, add the labelmaps.
        """
        # Perform a pixel-wise addition of ResSrcLabmapToCopyIm and the 
        # labelmap of Target that the segment is to be copied to: 
        #ADDLabmapIm = AddImages(ResSrcLabmapToCopyIm, TrgLabmapIms[ToRoiNum])
        
        # Create a new list of Target Labelmaps to include the new/modified 
        # ROI/labelmap as well as any others (unmodified):
        NewTrgLabmapIms = []
        
        for RoiNum in range(len(TrgLabmapIms)):
            if RoiNum == ToRoiNum:
                NewTrgLabmapIms.append(ORLabmapIm)
                #NewTrgLabmapIms.append(ADDLabmapIm) # <-- useful for visualisation but not correct for PixelData
            
            else:
                NewTrgLabmapIms.append(TrgLabmapIms[RoiNum])
                
                
        """
        Now need to convert NewTrgLabmapIms into a pixel array.
        """
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrByRoi = []
        NewTrgPFFGStoDcmIndsByRoi = []
        NewTrgPFFGStoDcmInds = []
        NewTrgRoiNumsByRoi = []
        NewTrgRoiNums = []
        
        for RoiNum in range(len(NewTrgLabmapIms)):
            PixArr, FTSInds = Image2PixArr(LabmapImage=NewTrgLabmapIms[RoiNum])
            
            NewTrgPixArrByRoi.append(PixArr)
            
            NewTrgPFFGStoDcmIndsByRoi.append(FTSInds)
            
            NewTrgPFFGStoDcmInds.extend(FTSInds)
            
            NewTrgRoiNumsByRoi.append([RoiNum]*len(FTSInds))
            
            NewTrgRoiNums.extend([RoiNum]*len(FTSInds))
            
        
        # Combine all frames PixArrByRoi:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoDcmInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for RoiNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoDcmIndsByRoi[RoiNum])):
                NewTrgPixArr[n] = NewTrgPixArrByRoi[RoiNum][i]
                
                n += 1
        
        
        
        print(f'\nThere are {NewTrgPixArr.shape[0]} frames in NewTrgPixArr')
        print(f'NewTrgPFFGStoDcmIndsByRoi: {NewTrgPFFGStoDcmIndsByRoi}')
        print(f'NewTrgRoiNumsByRoi: {NewTrgRoiNumsByRoi}')
        
        
        
        
        # Create a list of frame numbers for NewTrg that indicate the frame
        # numbers that the masks originated from, and a similar list indicating
        # where they came from:
        #NewTrgFrameNums = []
        #NewTrgFrameSrcs = []
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgRoiNums, NewTrgPFFGStoDcmInds, 
                                     FromRoiNum, FromSliceNum, SrcSeg, 
                                     ToRoiNum, LogToConsole=False)
        
        
        
        
        if LogToConsole:
            print('\n\n')
            title = 'Results after resampling Source LabelMap image and ' \
                    + 'Source-to-copy LabelMap image:'
            ul = '*' * len(title)
            
            print('\n' + title)
            print(ul)
            
            print('\n*The Source and Target labelmap image size:*\n',
                          f'   Source:                   {SrcLabmapIm.GetSize()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetSize()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetSize()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSize()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSize()}')
                
            print('\n*The Source and Target labelmap image voxel spacings:*\n',
                          f'   Source:                   {SrcLabmapIm.GetSpacing()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetSpacing()} \n',
                          f'   Source-to-copy            {SrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSpacing()}')
                
            print('\n*The Source and Target labelmap image Origin:*\n',
                          f'   Source:                   {SrcLabmapIm.GetOrigin()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetOrigin()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Target:                   {TrgLabmapIm.GetOrigin()}')
                
            print('\n*The Source and Target labelmap image Direction:*\n',
                          f'   Source:                   {SrcLabmapIm.GetDirection()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetDirection()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Target:                   {TrgLabmapIm.GetDirection()}')
        
        
        
        if LogToConsole:
            SrcLabmapImMax = ImageMax(SrcLabmapIm)
            SrcLabmapImMin = ImageMin(SrcLabmapIm)
            
            ResSrcLabmapImMax = ImageMax(ResSrcLabmapIm)
            ResSrcLabmapImMin = ImageMin(ResSrcLabmapIm)
            
            SrcLabmapToCopyImMax = ImageMax(SrcLabmapToCopyIm)
            SrcLabmapToCopyImMin = ImageMin(SrcLabmapToCopyIm)
            
            ResSrcLabmapToCopyImMax = ImageMax(ResSrcLabmapToCopyIm)
            ResSrcLabmapToCopyImMin = ImageMin(ResSrcLabmapToCopyIm)
            
            TrgLabmapImMax = ImageMax(TrgLabmapIm)
            TrgLabmapImMin = ImageMin(TrgLabmapIm)
            
            print('\n*The Source and Target labelmap image minimum:*\n',
                          f'   Source:                   {SrcLabmapImMin} \n',
                          f'   Resampled Source:         {ResSrcLabmapImMin} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMin} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMin} \n',
                          f'   Target:                   {TrgLabmapImMin}')
            
            print('\n*The Source and Target labelmap image maximum:*\n',
                          f'   Source:                   {SrcLabmapImMax} \n',
                          f'   Resampled Source:         {ResSrcLabmapImMax} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMax} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMax} \n',
                          f'   Target:                   {TrgLabmapImMax}')
        
        
        
        
        
        # Plot results:
        #PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcImage, 
        #                              ResSrcIm=ResSrcImage, 
        #                              TrgIm=TrgImage, 
        #                              LogToConsole=False)
        
        
        #return SrcImage, ResSrcImage, TrgImage
        #return SrcLabmapIm, ResSrcLabmapIm, TrgLabmapIm
        #return SrcLabmapIm, ResSrcLabmapIm, TrgLabmapIm,\
        #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm
        #return SrcImage, ResSrcImage, TrgImage,\
        #       SrcLabmapIm, ResSrcLabmapIm,\
        #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm, TrgLabmapIm
        #return SrcImage, ResSrcImage, TrgImage,\
        #       SrcLabmapIms, ResSrcLabmapIms,\
        #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
        #       TrgLabmapIms, NewTrgLabmapIms,\
        #       NewTrgSeg
         
         
        
        
    if CaseNum == '4':
        print(f'\nApplying Case 4..')
        
        
        
    if CaseNum == '5':
        print(f'\nApplying Case 5..')
        
        
    
        
        
    """
    The following was moved to ModifySegTagVals:
    # Check if NewTrgSeg exists (temporary check):
    if 'NewTrgSeg' in locals():
    
        # Generate a new SOP Instance UID:
        NewTrgSeg.SOPInstanceUID = pydicom.uid.generate_uid()
        
        # Generate a new Series Instance UID:
        NewTrgSeg.SeriesInstanceUID = pydicom.uid.generate_uid()
        
        
        return NewTrgSeg
    
    #else:
        #return []
        #return ResSrcImage
    """
        
        
        
    return SrcImage, ResSrcImage, TrgImage,\
           SrcLabmapIms, ResSrcLabmapIms,\
           SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
           TrgLabmapIms, NewTrgLabmapIms,\
           NewTrgSeg
           
           
           
           
           
           



def MappedCopySegAcrossSeries(SrcDcmDir, SrcSegFpath, FromSliceNum, FromSegLabel,
                              TrgDcmDir, TrgSegFpath, LogToConsole=False):
                                
    """
    16/11:  Removed ToSegLabel from list of inputs. Rather than copying to an
    existing segment copy to a new one.
    
    18/11:  See CopySegAcrossSeries = DirectCopySegAcrossSeries 
                                      + MappedCopySegAcrossSeries
    and the presence of ToSliceNum will indicate that a direct copy is to be
    made. The default value of ToSliceNum will be None, which will indicate
    that a relationship-preserving copy is to be made.
    
    
    Make a relationship-preserving copy of a single segmentation on a given   
    slice to a different DICOM series.  Since the relationship is preserved,  
    the voxels corresponding to the voxels copied to the Target ROI volume 
    spatially relate to the voxels corresponding to the voxels copied from the 
    Source ROI volume.
    
    There are four possible cases that will need to be accommodated*:
        
        2b. The Source and Target images have the same FOR, IOP, IPP, PS and ST 
        
        3b. The Source and Target images have the same FOR and IOP but have
            different PS and/or ST (if different ST they will likely have 
            different IPP)
        
        4. The Source and Target images have the same FOR but have different
           IOP and IPP, and possibly different PS and/or ST
           
        5. The Source and Target images have different FOR (hence different IPP
           and IOP), and possibly different ST and/or PS
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, ST = SliceThickness, PS = PixelSpacing
    
    * The numbering deliberately starts at 2b for consistency with the 5 
    overall cases.  Cases 1, 2a and 3a (=Direct copies) are covered by other 
    functions (CopySegWithinSeries for Case 1, and 
    DirectCopySegAcrossSeries for Cases 2a and 3a).
    
    
    Inputs:
    ------
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    SrcSegFpath : string
        Full path of the Source DICOM SEG file.
                             
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
                       
    FromSegLabel : string
        All or part of the Source segment label containing the segmentation to 
        be copied.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    TrgSegFpath : string
        Full path of the Target DICOM SEG file. May also be empty string ('') 
        if there is no DICOM SEG for the Target Series.
                        
    LogToConsole : boolean
        Denotes whether some intermediate results will be logged to the console
          
                  
    Output:
        
    NewTrgSeg : Pydicom object
        Modified Target ROI object
    """
    
    
    # Import the Source and Target SEG ROI:
    SrcSeg = dcmread(SrcSegFpath)
    if TrgSegFpath: # if TrgSegFpath is not empty
        TrgSeg = dcmread(TrgSegFpath)
    
    # Compare the modalities of the ROIs:
    if TrgSegFpath:
        if not IsSameModalities(SrcSeg, TrgSeg):
            msg = f"The Source ({SrcSeg.Modality}) and Target " \
                  + "({TrgSeg.Modality} modalities are different."
            
            raise Exception(msg)
            
    
    # Import the Source and Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir)
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target FrameOfReferenceUID:
    SrcFORuid = SrcDicoms[0].FrameOfReferenceUID
    TrgFORuid = TrgDicoms[0].FrameOfReferenceUID
    
    # Get the Image Attributes for Source and Target using SimpleITK:
    SrcFORuid, SrcSitkSize, SrcSitkSpacing, SrcSitkST,\
    SrcSitkIPP, SrcSitkDir = GetImageAttributes(DicomDir=SrcDcmDir,
                                                Package='sitk')
    
    TrgFORuid, TrgSitkSize, TrgSitkSpacing, TrgSitkST,\
    TrgSitkIPP, TrgSitkDir = GetImageAttributes(DicomDir=TrgDcmDir,
                                                Package='sitk')
    
    if LogToConsole:
        CompareSourceTargetImageAttributes(SrcDcmDir, TrgDcmDir)
    
    
    """
    I'm not sure if the items in DimensionIndexSequence,
    e.g. DimensionDescriptionLabel = 'ReferencedSegmentNumber'
                                   = 'ImagePositionPatient'
                                   = 'StackID'
                                   = 'InStackPositionNumber'
                                   
    are important.  While the label 'ReferencedIndexSequence' are common to 
    SEGs both generated by the OHIF-Viewer and RTS-to-SEG conversions, it seems 
    that the label 'ImagePositionPatient' is only found in OHIF-generated SEGs, 
    and the labels 'StackID' and 'InStackPositionNumber' from RTS-to-SEG 
    conversions.
    """
    
    # Get the number of segments in Source and Target:
    SrcNumOfSegs = len(SrcSeg.SegmentSequence)
    if TrgSegFpath:
        TrgNumOfSegs = len(TrgSeg.SegmentSequence)
    
    # Get the labels of all segments in the Source and Target SEGs:
    SrcSegLabels = GetRoiLabels(SrcSeg)
    #SrcNumOfSegs = len(SrcSegLabels)
    
    if TrgSegFpath:
        TrgSegLabels = GetRoiLabels(TrgSeg)
        #TrgNumOfSegs = len(TrgSegLabels)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            print(f'\nThe Source SEG has {SrcNumOfSegs} segments with',
                  f'label(s) {SrcSegLabels} \nand the Target SEG has',
                  f'{TrgNumOfSegs} segments with label(s) {TrgSegLabels}')
        else:
            print(f'\nThe Source SEG has {SrcNumOfSegs} segments with',
                  f'label(s) {SrcSegLabels}')
    
    
    
    # Get the DimensionIndexValues (which relate the segments to the item in 
    # the Referenced Instance Sequence):
    """
    Note:
        The DimensionalIndexValues are integers that start from 1 (unlike the
        RIS or PFFGS which begin at 0)!
    """
    SrcDIVs = GetDIVs(SrcSeg)
    if TrgSegFpath:
        TrgDIVs = GetDIVs(TrgSeg)
    
    # Group the DIVs by ROI:
    """
    Note:
        The DimensionalIndexValues are integers that start from 1 (unlike the
        RIS or PFFGS which begin at 0)!
    """
    SrcDIVsBySeg = GroupListBySegment(ListToGroup=SrcDIVs, DIVs=SrcDIVs)
    if TrgSegFpath:
        TrgDIVsBySeg = GroupListBySegment(ListToGroup=TrgDIVs, DIVs=TrgDIVs)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            print('\nThe Source SEG DimensionIndexValues grouped by segment',
                  f'are \n{SrcDIVsBySeg} \n(***Note: Counting from 1.***)\n\n',
                  'and the Target SEG DimensionIndexValues grouped by segment',
                  f'are \n{TrgDIVsBySeg} \n(***Note: Counting from 1.***)')
        else:
            print('\nThe Source SEG DimensionIndexValues grouped by segment',
                  f'are \n{SrcDIVsBySeg} \n(***Note: Counting from 1.***)')
        
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, 
    # the same grouped by ROI,
    # the PerFrameFunctionalGroupsSequence-to-DICOM slice indices, 
    # the same grouped by ROI, and 
    # the frame numbers grouped by ROI for Source and Target:
    """
    Note:
        The RIStoDcmInds and PFFGStoDcmInds are integers that start from 0 
        (unlike the DIVs which begin at 1)!
    """
    #SrcRIStoDcmInds = GetRIStoDcmInds(SrcSeg, SrcSOPuids)

    SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcSeg, SrcSOPuids)
    
    SrcFrameNums = list(range(len(SrcPFFGStoDcmInds)))
    
    # Group the above by segment:
    #SrcRIStoDcmIndsBySeg = GroupListBySegment(ListToGroup=SrcRIStoDcmInds, 
    #                                      DIVs=SrcDIVs) # <-- RIS has more elements than DIV
    SrcPFFGStoDcmIndsBySeg = GroupListBySegment(ListToGroup=SrcPFFGStoDcmInds, 
                                                DIVs=SrcDIVs)
    
    SrcFrameNumsBySeg = GroupListBySegment(ListToGroup=SrcFrameNums, 
                                           DIVs=SrcDIVs)
    
    if TrgSegFpath:
        #TrgRIStoDcmInds = GetRIStoDcmInds(TrgSeg, TrgSOPuids)
        
        TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgSeg, TrgSOPuids)
        
        TrgFrameNums = list(range(len(TrgPFFGStoDcmInds)))
        
        # Group the above by segment:
        #TrgRIStoDcmIndsBySeg = GroupListBySegment(ListToGroup=TrgRIStoDcmInds, 
        #                                      DIVs=TrgDIVs) # <-- RIS has more elements than DIV
        TrgPFFGStoDcmIndsBySeg = GroupListBySegment(ListToGroup=TrgPFFGStoDcmInds, 
                                                    DIVs=TrgDIVs)
        
        TrgFrameNumsBySeg = GroupListBySegment(ListToGroup=TrgFrameNums, 
                                               DIVs=TrgDIVs)
        
        if True:#LogToConsole:
            print('\nThe Source series has segmentations on slices    ',
                  f'{SrcPFFGStoDcmIndsBySeg} \nand the Target series has',
                  f'segmentations on slices {TrgPFFGStoDcmIndsBySeg}',
                  '\n(grouped by segment) \n(***Note: Counting from 0.***)')
            
            #print('\nThe RIS-to-DICOM slice num for the Source series is',
            #      f'\n{SrcRIStoDcmInds}, \nand for the Target series is',
            #      f'\n{TrgRIStoDcmInds}')
    
    else:
        if True:#LogToConsole:
            print('The Source series has segmentations on slices',
                  f'{SrcPFFGStoDcmIndsBySeg} \n(grouped by segment)',
                  '\n(***Note: Counting from 0.***)')
            
            #print('\nThe RIS-to-DICOM slice num for the Source series is',
            #      f'\n{SrcRIStoDcmInds}')
    
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcPFFGStoDcmInds:
        msg = "The Source series does not have a segmentation on slice " \
              + f"{FromSliceNum}"
        
        raise Exception(msg)
        
        return
    
    
    
    # Read in the Source and Target SimpleITK 3D images:
    SrcImage = ImportImage(SrcDcmDir)
    TrgImage = ImportImage(TrgDcmDir)
    
    if LogToConsole:
        title = 'Creating 3D labelmaps for *Source*...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    
    SrcPixArr = SrcSeg.pixel_array
    
    print(f'\nThe Source pixel array has shape {SrcPixArr.shape}')
    
    # Create SimpleITK images of the 3D labelmaps (= pixel arrays in same frame 
    # position as corresponding DICOM slices, with 2D zero masks in between)
    # for each segment:
    SrcLabmapIms = PixArr2ImagesBySeg(PixArr=SrcPixArr, 
                                      Image=SrcImage, 
                                      FrameToSliceIndsBySeg=SrcFrameNumsBySeg)
        
    
    if LogToConsole:
        title = 'Creating 3D labelmaps for *Target*...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    if TrgSegFpath:
        TrgPixArr = TrgSeg.pixel_array
        
        print(f'The Target pixel array has shape {TrgPixArr.shape}')
        
        TrgLabmapIms = PixArr2ImagesBySeg(PixArr=TrgPixArr, 
                                          Image=TrgImage, 
                                          FrameToSliceIndsBySeg=TrgFrameNumsBySeg)
            
        
    else:
        # Initialise the Target labelmap:
        #TrgLabmapIm = InitialiseLabmapIm(TrgImage)
        TrgLabmapIm = InitialiseImage(TrgImage)
        
        # For consistency with other IF case:
        TrgLabmapIms = []
        TrgLabmapIms.append(TrgLabmapIm)
        
        
    
    
    
    
    """ 
    SrcLabmapIm is the 3D labelmap of all segments within SrcSeg.
    
    But when looking to copy a single segment from Source to Target, I'll need
    to know how that single segment is resampled.  So create a new labelmap
    containing the segment to be copied only.
    """
    
    if LogToConsole:
        title = 'Creating 3D labelmap for *Source* containing segment to copy only...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    # Create SimpleITK image of the 3D labelmap for the segment to be copied 
    # only:
    """
    There may be more than one frame number in PixelArray that corresponds to 
    FromSliceNum...
    # The frame number in PixelArray that corresponds to FromSliceNum:
    FromFrameNum = SrcPFFGStoDcmInds.index(FromSliceNum)
    """
    
    # Determine which Source segment contains the segmentation to be copied:
    """
    SegNum is equivalent to SegmentNumber in SegmentSequence.
    """
    FromSegNum = [i for i, s in enumerate(SrcSegLabels) if FromSegLabel in s]
    FromSegNum = FromSegNum[0]
    
    if TrgSegFpath:
        # Determine which Target segment the Source segmentation is to be 
        # copied to:
        #ToSegNum = [i for i, s in enumerate(TrgSegLabels) if ToSegLabel in s]
        #ToSegNum = ToSegNum[0]
        
        # The segment will be added to a new segment:
        ToSegNum = len(TrgSegLabels) + 1
    
        print(f'\nFromSegNum = {FromSegNum} \nToSegNum = {ToSegNum}')
        
    else:
        print(f'\nFromSegNum = {FromSegNum}')
    
    # Determine which frame number in PixelArray corresponds to FromSliceNum:
    """
    FromFrameNum = SrcPFFGStoDcmIndsByRoi[FromRoiNum].index(FromSliceNum)  
    
    The problem with the above is that it's the frame number within the subset
    of frames pertaining to the segment given by FromSegNum - it doesn't relate 
    to the frame number in PixelArray (which contains all frames for all 
    segments).
    """
    n = 0 # initialise frame number counter
    
    for SegNum in range(SrcNumOfSegs):
        if SegNum == FromSegNum:
            # This segment contains the frame to be copied, whose position is
            # given by the index of FromSliceNum is in 
            # SrcPFFGStoDcmIndsBySeg[SegNum]
            # plus any frames that preceeded it (i.e. the frame counter n):
            FromFrameNum = n + SrcPFFGStoDcmIndsBySeg[SegNum].index(FromSliceNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(SrcPFFGStoDcmIndsBySeg[SegNum])
    
    
    print(f'\nFromSliceNum = {FromSliceNum} relates to frame number',
          f'{FromFrameNum} in PixelArray (***Note: Counting from 0.***)')
    
    # Determine which frame number(s) in PixelArray correspond to FromSliceNum: 
    #FromFrameNums = []
    #
    #for RoiNum in range(SrcNumOfROIs):
    #    if FromSliceNum in SrcPFFGStoDcmIndsByRoi[RoiNum]:
    #        ind = SrcPFFGStoDcmIndsByRoi[RoiNum].index(FromSliceNum)
    #        
    #        FromFrameNums.append(ind)
    #     
    #print(f'\nFromSliceNum = {FromSliceNum} relates to frame number(s)',
    #      f'{FromFrameNums} in PixelArray')
    
    # Knowing which ROI contains the frame to be copied and the frame number
    # defines the pixel array to be copied:
    #FromPixelArr = SrcSeg.pixel_array[FromFrameNum]
    
    # Modify the pixel array so that all but the frame to be copied has zeros:
    """
    SrcPixArrToCopy = GetFrameFromPixArr(PixArr=FromPixelArr, 
                                         FrameNum=FromFrameNum)
    
    I either need to pass the entire PixelArray and specify the FrameNum to be
    copied (as was my original intention with the function GetFrameFromPixArr),
    or pass the frame to be copied (FromPixelArr) and modify GetFrameFromPixArr
    accordingly!
    """
    SrcPixArrToCopy = GetFrameFromPixArr(PixArr=SrcSeg.pixel_array, 
                                         FrameNum=FromFrameNum)
    
    if LogToConsole:
        print(f'\nSrcPixArrToCopy.max() = {SrcPixArrToCopy.max()}')
    
    
    """ Moved to Case 3 since not required for case 2: 
    SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                     Image=SrcImage,
                                     FrameToSliceInds=[FromSliceNum])
    
    if LogToConsole:
        print(f'\nMax value of SrcLabmapToCopyIm = {ImageMax(SrcLabmapToCopyIm)}')
    """
    
    
    # Check which case follows:
    if SrcFORuid == TrgFORuid:
        if LogToConsole:
            print('\nThe Source and Target FrameOfReferenceUID are the same')
        
        # Compare the Source and Target directions. Get the maximum value of 
        # their absolute differences:
        AbsDiffs = abs(np.array(SrcSitkDir) - np.array(TrgSitkDir))
        MaxAbsDiff = max(AbsDiffs)
        
        # Consider the directions different if any of the vector differences 
        # are greater than epsilon:
        epsilon = 1e-06
        
        #if SrcDirs == TrgDirs: 
        """ The line above will result in insignificant differences (e.g. due
        to quantisation/pixelation) resulting in False, so instead consider
        MaxAbsDiff and epsilon. 
        """
        if MaxAbsDiff < epsilon:
            if LogToConsole:
                print('The Source and Target directions are the same')
            
            if SrcSitkSpacing == TrgSitkSpacing:
                if LogToConsole:
                    print('The Source and Target voxel sizes are the same')
                
                CaseNum = '2b'
                
            else:
                if LogToConsole:
                    print('The Source and Target voxel sizes are different:\n',
                          f'   Source: {SrcSitkSpacing} \n',
                          f'   Target: {TrgSitkSpacing}')
                    
                CaseNum = '3b'
            
        else:
            if LogToConsole:
                print('The Source and Target directions are different:\n',
                      f'   Source: [{SrcSitkDir[0]}, {SrcSitkDir[1]}, {SrcSitkDir[2]},',
                      f'            {SrcSitkDir[3]}, {SrcSitkDir[4]}, {SrcSitkDir[5]},',
                      f'            {SrcSitkDir[6]}, {SrcSitkDir[7]}, {SrcSitkDir[8]}]\n',
                      f'   Target: [{TrgSitkDir[0]}, {TrgSitkDir[1]}, {TrgSitkDir[2]},',
                      f'            {TrgSitkDir[3]}, {TrgSitkDir[4]}, {TrgSitkDir[5]},',
                      f'            {TrgSitkDir[6]}, {TrgSitkDir[7]}, {TrgSitkDir[8]}]')
                
                print(f'AbsDiffs = {AbsDiffs}')
            
            
            
            CaseNum = '4'
    
    else:
        if LogToConsole:
            print('The Source and Target FrameOfReferenceUID are different')
        
        CaseNum = '5'
        
    
    if LogToConsole:    
        print('The Source and Target origins:\n',
                          f'   Source: {SrcSitkIPP[0]} \n',
                          f'   Target: {TrgSitkIPP[0]}')
        
        print('The Source and Target image dims:\n',
                          f'   Source: {SrcSitkSize} \n',
                          f'   Target: {TrgSitkSize}')
    
    
    
    if LogToConsole:
        print(f'\nCase {CaseNum} applies.')
    
    
    if CaseNum == '2b':
        """
        The origins for Source and Target are the same since they have the same
        FrameOfReferenceUID.  And since the PixelSpacing and SliceThickness are
        also the same, ToSliceNum = FromSliceNum, and the relationships will be
        preserved.
        
        26/10: It's seems that having the same FOR does not guarantee the same
        origin.  Series 9 and 6 for Subject 011, MR4 have different origins!
        So need to check whether origins are the same?...
        
        17/11:  This case is up-to-date now but maintaining the same outputs
        for different cases is very messy and inefficient (e.g. see creation of
        unnecessary variables for Case 2b).
        """
        
        print('\nApplying Case 2b..')
    
        ToSliceNum = deepcopy(FromSliceNum)
        
        if LogToConsole:
            print(f'\nCopying slice {FromSliceNum} from Source SEG to slice',
                  f'{ToSliceNum} in Target SEG..')
        
        
        # The new list of Target PerFrameFunctionalGroupsSequence-to-DICOM
        # slice index = the original list with FromSliceNum appended to it:
        NewTrgPFFGStoDcmInds = deepcopy(TrgPFFGStoDcmInds)
        NewTrgPFFGStoDcmInds.append(FromSliceNum)
        
        # Similarly for the new Target PerFrameFunctionalGroupsSequence-to-
        # DICOM slice index grouped by segment:
        NewTrgPFFGStoDcmIndsBySeg = deepcopy(TrgPFFGStoDcmIndsBySeg)
        NewTrgPFFGStoDcmIndsBySeg.append([FromSliceNum])
        
        # The new Target segment numbers are:
        NewTrgSegNums = []
        
        for SegNum in range(len(NewTrgPFFGStoDcmIndsBySeg)):
            # The number of frames for this segment number:
            Nframes = len(NewTrgPFFGStoDcmIndsBySeg[SegNum])
            
            print(f'\nSegNum = {SegNum}, Nframes = {Nframes}')
            
            NewTrgSegNums.extend([SegNum]*Nframes)
            
            
        # The new Target pixel array = original Target pixel array with
        # SrcPixArrToCopy added to it:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoDcmInds), 
                                 TrgPixArr.shape[1], 
                                 TrgPixArr.shape[2]))
        
        NewTrgPixArr[:TrgPixArr.shape[0]] = TrgPixArr
        NewTrgPixArr[-1] = SrcPixArrToCopy
        
        
        if LogToConsole:
            print('\nTrgPFFGStoDcmInds    =', TrgPFFGStoDcmInds)
            print('\nNewTrgPFFGStoDcmInds =', NewTrgPFFGStoDcmInds)
            
        print(f'\nThere are {NewTrgPixArr.shape[0]} frames in NewTrgPixArr')
        print(f'NewTrgPFFGStoDcmInds: {NewTrgPFFGStoDcmInds}')
        print(f'NewTrgSegNums: {NewTrgSegNums}')
        
        # Modify select tags in NewTrgSeg:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoDcmInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole=False)
        
        
        """ Note:
            
        The following variables:
            
        ResSrcImage
        SrcLabmapIms
        ResSrcLabmapIms
        SrcLabmapToCopyIm
        ResSrcLabmapToCopyIm
        NewTrgLabmapIms
        
        were not used for this case but this main function outputs them:
        
        return SrcImage, ResSrcImage, TrgImage,\
               SrcLabmapIms, ResSrcLabmapIms,\
               SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               NewTrgSeg
               
        So return:
            
        SrcImage for ResSrcImage
        SrcLabmapIms for ResSrcLabmapIms
        
        and create:
            
        SrcLabmapToCopyIm
        NewTrgLabmapIms
        
        even though they're not needed for this case, and return:
            
        SrcLabmapToCopyIm for ResSrcLabmapToCopyIm
        """
        
        # Create SrcLabmapToCopyIm even though it's not needed:
        SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                         Image=SrcImage,
                                         FrameToSliceInds=[FromSliceNum])
        
        # Create NewTrgLabmapIms even though it's not needed:
        NewTrgLabmapIms = []
        
        for TrgLabmapIm in TrgLabmapIms:
            NewTrgLabmapIms.append(TrgLabmapIm)
            
        # ... and append the Source labelmap to be copied:    
        NewTrgLabmapIms.append(SrcLabmapToCopyIm)
        
        return SrcImage, SrcImage, TrgImage,\
               SrcLabmapIms, SrcLabmapIms,\
               SrcLabmapToCopyIm, SrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               NewTrgSeg
               
        
        """
        23/10:  Still need to verify that Case 2 provides the expected results.
        
        26/10:  Resampling will need to be performed if the origins are not the
        same...
        """
        
        
    if CaseNum == '3b':
        print(f'\nApplying Case 3b..')
        
        # Resample the Source image:
        """ This is not required but useful info """
        ResSrcImage = ResampleImage(Image=SrcImage, RefImage=TrgImage, 
                                    Interpolation='Linear')
        
        if LogToConsole:                            
            print('Results after resampling Source image:\n')
            
            print('The Source and Target image size:\n',
                          f'   Original Source:  {SrcImage.GetSize()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSize()} \n',
                          f'   Target:           {TrgImage.GetSize()}')
                
            print('The Source and Target voxel spacings:\n',
                          f'   Original Source:  {SrcImage.GetSpacing()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSpacing()} \n',
                          f'   Target:           {TrgImage.GetSpacing()}')
                
            print('The Source and Target Origin:\n',
                          f'   Original Source:  {SrcImage.GetOrigin()} \n',
                          f'   Resampled Source: {ResSrcImage.GetOrigin()} \n',
                          f'   Target:           {TrgImage.GetOrigin()}')
                
            print('The Source and Target Direction:\n',
                          f'   Original Source:  {SrcImage.GetDirection()} \n',
                          f'   Resampled Source: {ResSrcImage.GetDirection()} \n',
                          f'   Target:           {TrgImage.GetDirection()}')
                                   
                                    
        # Resample the Source labelmap images:
        """ This is not required but useful info """
        ResSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            ResSrcLabmapIm = ResampleImage(Image=SrcLabmapIm, RefImage=TrgImage, 
                                           Interpolation='NearestNeighbor')
            
            ResSrcLabmapIms.append(ResSrcLabmapIm)
        
        
        # Convert pixel array to be to be copied to image:
        SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                         Image=SrcImage,
                                         FrameToSliceInds=[FromSliceNum])
    
        if LogToConsole:
            print(f'\nMax value of SrcLabmapToCopyIm = {ImageMax(SrcLabmapToCopyIm)}')
        
        
        # Resample the Source labelmap to be copied:
        ResSrcLabmapToCopyIm = ResampleImage(Image=SrcLabmapToCopyIm, 
                                             RefImage=TrgImage, 
                                             Interpolation='NearestNeighbor')
        
        # Create a new list of Target Labelmaps including the existing 
        # labelmaps...:
        NewTrgLabmapIms = []
        
        for TrgLabmapIm in TrgLabmapIms:
            NewTrgLabmapIms.append(TrgLabmapIm)
            
        # ... and append the resampled Source labelmap to be copied:    
        NewTrgLabmapIms.append(ResSrcLabmapToCopyIm)
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrBySeg = []
        NewTrgPFFGStoDcmIndsBySeg = []
        NewTrgPFFGStoDcmInds = []
        NewTrgSegNumsBySeg = []
        NewTrgSegNums = []
        
        for SegNum in range(len(NewTrgLabmapIms)):
            PixArr, F2SInds = Image2PixArr(LabmapIm=NewTrgLabmapIms[SegNum])
            
            NewTrgPixArrBySeg.append(PixArr)
            
            NewTrgPFFGStoDcmIndsBySeg.append(F2SInds)
            
            NewTrgPFFGStoDcmInds.extend(F2SInds)
            
            NewTrgSegNumsBySeg.append([SegNum]*len(F2SInds))
            
            NewTrgSegNums.extend([SegNum]*len(F2SInds))
            
        
        # Combine all frames:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoDcmInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for SegNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoDcmIndsBySeg[SegNum])):
                NewTrgPixArr[n] = NewTrgPixArrBySeg[SegNum][i]
                
                n += 1
        
        
        print(f'\nThe new Target pixel array has shape {NewTrgPixArr.shape}')
        print(f'\nNewTrgPFFGStoDcmIndsBySeg: {NewTrgPFFGStoDcmIndsBySeg}')
        print(f'NewTrgSegNumsBySeg: {NewTrgSegNumsBySeg}')
        print(f'NewTrgPFFGStoDcmInds: {NewTrgPFFGStoDcmInds}')
        print(f'NewTrgSegNums: {NewTrgSegNums}')
        
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoDcmInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole=False)
        
        
        if LogToConsole:
            print('\n\n')
            title = 'Results after resampling Source LabelMap image and ' \
                    + 'Source-to-copy LabelMap image:'
            ul = '*' * len(title)
            
            print('\n' + title)
            print(ul)
            
            print('\n*The Source and Target labelmap image size:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetSize()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetSize()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetSize()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSize()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSize()}')
                
            print('\n*The Source and Target labelmap image voxel spacings:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetSpacing()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetSpacing()} \n',
                          f'   Source-to-copy            {SrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSpacing()}')
                
            print('\n*The Source and Target labelmap image Origin:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetOrigin()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetOrigin()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Target:                   {TrgLabmapIm.GetOrigin()}')
                
            print('\n*The Source and Target labelmap image Direction:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetDirection()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetDirection()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Target:                   {TrgLabmapIm.GetDirection()}')
        
        
        
        if LogToConsole:
            #SrcLabmapImMax = ImageMax(SrcLabmapIm)
            #SrcLabmapImMin = ImageMin(SrcLabmapIm)
            
            #ResSrcLabmapImMax = ImageMax(ResSrcLabmapIm)
            #ResSrcLabmapImMin = ImageMin(ResSrcLabmapIm)
            
            SrcLabmapToCopyImMax = ImageMax(SrcLabmapToCopyIm)
            SrcLabmapToCopyImMin = ImageMin(SrcLabmapToCopyIm)
            
            ResSrcLabmapToCopyImMax = ImageMax(ResSrcLabmapToCopyIm)
            ResSrcLabmapToCopyImMin = ImageMin(ResSrcLabmapToCopyIm)
            
            TrgLabmapImMax = ImageMax(TrgLabmapIm)
            TrgLabmapImMin = ImageMin(TrgLabmapIm)
            
            print('\n*The Source and Target labelmap image minimum:*\n',
                          #f'   Source:                   {SrcLabmapImMin} \n',
                          #f'   Resampled Source:         {ResSrcLabmapImMin} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMin} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMin} \n',
                          f'   Target:                   {TrgLabmapImMin}')
            
            print('\n*The Source and Target labelmap image maximum:*\n',
                          #f'   Source:                   {SrcLabmapImMax} \n',
                          #f'   Resampled Source:         {ResSrcLabmapImMax} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMax} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMax} \n',
                          f'   Target:                   {TrgLabmapImMax}')
            
        
        return SrcImage, ResSrcImage, TrgImage,\
               SrcLabmapIms, ResSrcLabmapIms,\
               SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               NewTrgSeg
    
         
        
        
    if CaseNum == '4':
        print(f'\nApplying Case 4..')
        
        
        
    if CaseNum == '5':
        print(f'\nApplying Case 5..')
    
        
        
    """
    The following was moved to each CaseNum since some of the outputs don't
    exist for some cases (e.g. resampled Source image for Case 2b). """
    #return SrcImage, ResSrcImage, TrgImage,\
    #       SrcLabmapIms, ResSrcLabmapIms,\
    #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
    #       TrgLabmapIms, NewTrgLabmapIms,\
    #       NewTrgSeg
    
    
    
    

    



def CopySeg_OLD(SrcDcmDir, SrcSegFpath, FromSliceNum, FromSegLabel,
            TrgDcmDir, TrgSegFpath=None, ToSliceNum=None, 
            LogToConsole=False):
                                
    """
    New version created on 30/11/20.
    
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
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                            
    SrcSegFpath : string
        Full path of the Source DICOM SEG file.
                             
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
                       
    FromSegLabel : string
        All or part of the Source segment label containing the segmentation to 
        be copied.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
                            
    TrgSegFpath : string (default = None)
        Full path of the Target DICOM SEG file. Default value is None, e.g. 
        if there is no DICOM SEG for the Target Series.
        
    ToSliceNum : integer (default = None)
        Slice index within the Target DICOM stack where the segmentation is to
        be copied to (counting from 0).  This only applies for Direct copies,
        hence the default value None.
                        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
          
                  
    Outputs:
    -------
        
    SrcImage : SimpleITK object
        Source 3D image.
        
    _SrcImage : SimpleITK object
        Source 3D case-dependent image. If Main use case:
             = 1 or 2: SrcImage is returned unchanged.
             = 3 or 4: SrcImage resampled to TrgImage (ResSrcImage is returned).
             = 5: SrcImage registered to TrgImage (RegSrcImage is returned).
        
    TrgImage : SimpleITK object
        Target 3D image.
        
    SrcLabmapIms : List of SimpleITK objects
        List of Source 3D labelmaps - one labelmap per segment.
    
    _SrcLabmapIms : List of SimpleITK objects
        List of Source 3D case-dependent image - one labelmap per segment. If 
        Main use case:
            = 1 or 2: SrcLabmapIms is returned unchanged.
            = 3 or 4: Each SrcLabmapIm in ScrLabmapIms is resampled to TrgImage 
            (ResSrcLabmapIms is returned).
            = 5: Each SrcLabmapIm in SrcLabmapIms is registered to TrgImage 
            (RegSrcLabmapIms is returned).
    
    SrcLabmapToCopyIm : SimpleITK object
        Source 3D labelmap of the segmentation in FromSegLabel to be copied.
        
    _SrcLabmapToCopyIm : SimpleITK object
        Source 3D case-dependent image - one labelmap per segment. If Main use 
        case:
            = 1 or 2: SrcLabmapToCopyIm is returned unchanged.
            = 3 or 4: SrcLabmapToCopyIm is resampled to TrgImage 
            (ResSrcLabmapToCopyIm is returned).
            = 5: SrcLabmapToCopyIm is registered to TrgImage 
            (RegSrcLabmapToCopyIm is returned).
    
    TrgLabmapIms : List of SimpleITK objects
        List of Target 3D labelmaps - one labelmap per segment, if a Target SEG
        exists.  If not a list with a single empty 3D image is returned.
    
    NewTrgLabmapIms : List of SimpleITK objects
        List of new Target 3D labelmaps - one labelmap per segment.  If a 
        Target exists, the first n labelmaps are the same labelmaps in 
        TrgLabmapIms, and the final item is the labelmap that was copied.  If 
        not, a list with a single 3D labelmap is returned - the labelmap that
        was copied.  The new labelmap (last item in NewTrgLabmapIms) is case-
        dependent.  If Main use case:
            = 1 or 2: The new labelmap is simply SrcLabmapToCopyIm unchanged.
            = 3 or 4: The new labelmap is SrcLabmapToCopyIm resampled to 
            TrgImage (ResSrcLabmapToCopyIm is returned).
            = 5: The new labelmap is SrcLabmapToCopyIm transformed to the
            Target domain (TxSrcLabmapToCopyIm is returned).
    
    SrcSeg : Pydicom object
        Source SEG object.
        
    TrgSeg : Pydicom object
        Original Target SEG object or a copy of SrcSeg if TrgSegFpath = ''.
        
    NewTrgSeg : Pydicom object
        New/modified Target SEG object.
    
    SrcDicoms : List of Pydicom objects
        List of Source DICOM objects.
    
    TrgDicoms : List of Pydicom objects
        List of Target DICOM objects.
    
    
        
    Notes:
    -----
    
    16/11:  Removed ToSegLabel from list of inputs. Rather than copying to an
    existing segment copy to a new one.
    
    18/11:  The presence of ToSliceNum will indicate that a direct copy is to 
    be made. The default value of ToSliceNum will be None, which will indicate
    that a relationship-preserving copy is to be made.
    """
    
    
    # Determine which Use Case applies:
    UseCase = WhichUseCase(SrcDcmDir, TrgDcmDir, ToSliceNum, LogToConsole)
    
    if True:#LogToConsole:
        #print(f'\nCase {UseCase} applies.')
        title = f'Case {UseCase} applies.'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
    
    """ 
    **************************************************************************
    Import the list of DICOM objects, 3D images and SEGs. 
    **************************************************************************
    """
    SrcDicoms = ImportDicoms(SrcDcmDir) # not needed but useful
    TrgDicoms = ImportDicoms(TrgDcmDir) # needed for ModifySegTagVals
    
    SrcImage = ImportImage(SrcDcmDir)
    TrgImage = ImportImage(TrgDcmDir)
    
    SrcSeg = dcmread(SrcSegFpath)
    
    SrcPixArr = SrcSeg.pixel_array
    
    if TrgSegFpath: # if TrgSegFpath is not empty
        TrgSeg = dcmread(TrgSegFpath)
    
    
    # Compare the modalities of the ROIs:
    if TrgSegFpath:
        if not IsSameModalities(SrcSeg, TrgSeg):
            msg = f"The Source ({SrcSeg.Modality}) and Target " \
                  + "({TrgSeg.Modality} modalities are different."
            
            raise Exception(msg)
    
    
    """ 
    **************************************************************************
    Get the number of segments and the segment labels. 
    **************************************************************************
    """
    SrcNumOfSegs = len(SrcSeg.SegmentSequence)
    if TrgSegFpath:
        TrgNumOfSegs = len(TrgSeg.SegmentSequence)
    
    SrcSegLabels = GetRoiLabels(SrcSeg)
    #SrcNumOfSegs = len(SrcSegLabels)
    
    if TrgSegFpath:
        TrgSegLabels = GetRoiLabels(TrgSeg)
        #TrgNumOfSegs = len(TrgSegLabels)
    else:
        TrgSegLabels = [] # 26/11
    
    if True:#LogToConsole:
        if TrgSegFpath:
            msg = f'\nThe Source SEG has {SrcNumOfSegs} segments with ' \
                  + f'label(s) {SrcSegLabels} \nand the Target SEG has ' \
                  + f'{TrgNumOfSegs} segments with label(s) {TrgSegLabels}'
        else:
            msg = f'\nThe Source SEG has {SrcNumOfSegs} segments with ' \
                  + f'label(s) {SrcSegLabels}'
            
        print(msg)
    
    
    
    """ 
    **************************************************************************
    Get the DimensionIndexValues (which relate the segments to the item in the
    Referenced Instance Sequence) and group by segment (required for 
    SrcFrameNumsBySeg).
    **************************************************************************
    
    Note:
        The DIVs are integers that start from 1 (unlike the RIS or PFFGS which 
        begin at 0)!
    """
    SrcDIVs = GetDIVs(SrcSeg)
    if TrgSegFpath:
        TrgDIVs = GetDIVs(TrgSeg)
    
    SrcDIVsBySeg = GroupListBySegment(ListToGroup=SrcDIVs, DIVs=SrcDIVs)
    if TrgSegFpath:
        TrgDIVsBySeg = GroupListBySegment(ListToGroup=TrgDIVs, DIVs=TrgDIVs)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            msg = '\nThe Source SEG DimensionIndexValues grouped by segment '\
                  + f'are \n{SrcDIVsBySeg} \n(***Note: Counting from 1.***)'\
                  + '\n\nand the Target SEG DimensionIndexValues grouped by '\
                  + f'segment are \n{TrgDIVsBySeg} \n(***Note: Counting from '\
                  + '1.***)'
        else:
            msg = '\nThe Source SEG DimensionIndexValues grouped by segment '\
                  + f'are \n{SrcDIVsBySeg} \n(***Note: Counting from 1.***)'
            
        print(msg)
        
    
    """
    **************************************************************************
    Get the DICOM SOP UIDs, PerFrameFunctionalGroupsSequence-to-slice 
    indices (PFFGStoSliceInds) and PFFGStoSliceInds grouped by segment
    (PFFGStoSliceIndsBySeg), and the frame numbers grouped by segment 
    (required for PixArr2ImagesBySeg()).
    **************************************************************************
    
    Note:
        The PFFGStoSliceInds and FrameNums are integers that start from 0 (unlike 
        the DIVs which begin at 1)!
    """
    
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    SrcPFFGStoSliceInds = GetPFFGStoSliceInds(SrcSeg, SrcSOPuids)
    
    SrcPFFGStoSliceIndsBySeg = GroupListBySegment(ListToGroup=SrcPFFGStoSliceInds, 
                                                  DIVs=SrcDIVs)
    
    SrcFrameNums = list(range(len(SrcPFFGStoSliceInds)))
    
    SrcFrameNumsBySeg = GroupListBySegment(ListToGroup=SrcFrameNums, 
                                           DIVs=SrcDIVs)
    
    if TrgSegFpath:
        TrgPFFGStoSliceInds = GetPFFGStoSliceInds(TrgSeg, TrgSOPuids)
        
        TrgPFFGStoSliceIndsBySeg = GroupListBySegment(ListToGroup=TrgPFFGStoSliceInds, 
                                                      DIVs=TrgDIVs)
        
        TrgFrameNums = list(range(len(TrgPFFGStoSliceInds)))
        
        TrgFrameNumsBySeg = GroupListBySegment(ListToGroup=TrgFrameNums, 
                                               DIVs=TrgDIVs)
        
        if True:#LogToConsole:
            print('\nThe Source series has segmentations on slices    ',
                  f'{SrcPFFGStoSliceIndsBySeg} \nand the Target series has',
                  f'segmentations on slices {TrgPFFGStoSliceIndsBySeg}',
                  '\n(grouped by segment) \n(***Note: Counting from 0.***)')
    
    else:
        TrgPFFGStoSliceInds = []
        
        TrgPFFGStoSliceIndsBySeg = []
        
        TrgFrameNums = []
        
        TrgFrameNumsBySeg = []
        
        if True:#LogToConsole:
            print('The Source series has segmentations on slices',
                  f'{SrcPFFGStoSliceIndsBySeg} \n(grouped by segment)',
                  '\n(***Note: Counting from 0.***)')
    
    
    """ 
    Verify that a segment exists for FromSliceNum (this check won't be 
    necessary in the OHIF-Viewer but a bad input could occur here).
    """
    if not FromSliceNum in SrcPFFGStoSliceInds:
        msg = "The Source series does not have a segmentation on slice " \
              + f"{FromSliceNum}"
        
        raise Exception(msg)
        
        return
    
    
    
    
    """
    **************************************************************************
    Create a 3D labelmap (as a SimpleITK image) for the segmentation to be 
    copied.
    **************************************************************************
    """    
    
    if LogToConsole:
        msg = 'Creating labelmap for *Source* segmentation to be copied...'
        PrintTitle(msg)
    
    # Determine which Source segment contains the segmentation to be copied:
    """
    Note:
        SegNum is equivalent to SegmentNumber in SegmentSequence.
    """
    FromSegNum = [i for i, s in enumerate(SrcSegLabels) if FromSegLabel in s]
    FromSegNum = FromSegNum[0]
    
    if TrgSegFpath:
        # Determine which Target segment the Source segmentation is to be 
        # copied to:
        #ToSegNum = [i for i, s in enumerate(TrgSegLabels) if ToSegLabel in s]
        #ToSegNum = ToSegNum[0]
        
        # The segment will be added to a new segment:
        ToSegNum = len(TrgSegLabels) + 1
    
        #print(f'\nFromSegNum = {FromSegNum} \nToSegNum = {ToSegNum}') # 26/11
        
    else:
        #print(f'\nFromSegNum = {FromSegNum}') # 26/11
        ToSegNum = 0 # 26/11
        
    print(f'\nFromSegNum = {FromSegNum} \nToSegNum = {ToSegNum}') # 26/11
    
    # Determine which frame number in PixelArray corresponds to FromSliceNum:
    #FromFrameNum = SrcPFFGStoSliceIndsByRoi[FromRoiNum].index(FromSliceNum)  
    """
    The problem with the above is that it's the frame number within the subset
    of frames pertaining to the segmentation given by FromSegNum - it doesn't  
    relate to the frame number in PixelArray (which contains all frames for all 
    segmentations).
    """
    n = 0 # initialise frame number counter
    
    for SegNum in range(SrcNumOfSegs):
        if SegNum == FromSegNum:
            """ 
            This segment contains the frame to be copied, whose position is
            given by the index of FromSliceNum is in SrcPFFGStoSliceIndsBySeg[SegNum]
            plus any frames that preceeded it (i.e. the frame counter n).
            """
            FromFrameNum = n + SrcPFFGStoSliceIndsBySeg[SegNum].index(FromSliceNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(SrcPFFGStoSliceIndsBySeg[SegNum]) 
            #n -= 1 # 19/11: Subtract 1 so that FromFrameNum begins from 0
    
    
    print(f'\nFromSliceNum = {FromSliceNum} relates to FromFrameNum =',
          f'{FromFrameNum} in PixelArray (***Note: Counting from 0.***)')

    
    # Create a new pixel array that only contains the frame to be copied:
    SrcPixArrToCopy = GetFrameFromPixArr(PixArr=SrcPixArr, 
                                         FrameNum=FromFrameNum)
    
    if LogToConsole:
        print(f'\nSrcPixArrToCopy.max() = {SrcPixArrToCopy.max()}')
    
    
    """ 
    Not needed for Case 2 but kept here for now since the function outputs
    SrcLabmapToCopyIm.
    """
    # Convert the pixel array (frame) to be to be copied into an SimpleITK
    # image:
    SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                     FrameToSliceInds=[FromSliceNum],
                                     RefIm=SrcImage)
    
    if LogToConsole:
        print(f'\nMax value of SrcLabmapToCopyIm = {ImageMax(SrcLabmapToCopyIm)}')
        
    #print(f'\n---> SrcLabmapToCopyIm.GetSize() = {SrcLabmapToCopyIm.GetSize()}')
    
    
    
    
    
    
    """
    **************************************************************************
    Create a list of 3D labelmaps (as SimpleITK images) for all segmentations
    in all segments.
    
    This isn't needed but is created since the function returns SrcLabmapIms. 
    **************************************************************************
    """
    
    if LogToConsole:
        PrintTitle('Creating 3D labelmaps for *Source*...')
    print(f'\nThe Source pixel array has shape {SrcPixArr.shape}')
    
    SrcLabmapIms = PixArr2ImagesBySeg(PixArr=SrcPixArr, 
                                      FrameToSliceIndsBySeg=SrcPFFGStoSliceIndsBySeg,
                                      FrameNumsBySeg=SrcFrameNumsBySeg,
                                      RefIm=SrcImage)
    
    if TrgSegFpath:
        TrgPixArr = TrgSeg.pixel_array
        
        if LogToConsole:
            PrintTitle('Creating 3D labelmaps for *Target*...')
        print(f'The Target pixel array has shape {TrgPixArr.shape}')
        
        TrgLabmapIms = PixArr2ImagesBySeg(PixArr=TrgPixArr, 
                                          FrameToSliceIndsBySeg=TrgPFFGStoSliceIndsBySeg,
                                          FrameNumsBySeg=TrgFrameNumsBySeg, 
                                          RefIm=TrgImage)
    else:
        """ 
        This isn't needed but create a list of a single empty image since the 
        function returns TrgLabmapIms 
        """
        # Initialise the Target labelmap:
        TrgLabmapIms = []
        TrgLabmapIms.append(InitialiseImage(TrgImage))
        
        
    
    
    
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    if 'a' in UseCase:
        """ 
        This is common for all Direct copy sub-use cases.
        """
        
        print(f'\nCopying slice {FromSliceNum} from Source SEG to slice',
                  f'{ToSliceNum} in Target SEG..')
        
        # Get pixel shift between FromSliceNum and ToSliceNum:
        PixShift = GetPixelShiftBetweenSlices(SrcImage, FromSliceNum, 
                                              ToSliceNum)
        
        print(f'\nPixShift = {PixShift}')
        
        
        
    if UseCase == '1a' or UseCase == '2a':
        """
        Main use case = 1 or 2:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            voxel spacings, and either are the same series (Main use case = 1) 
            or different series (Main use case = 2).
            
        Sub-case = a:
            Direct copy of the segmentation on slice FromSliceNum in segment 
            FromSegLabel for Source to a new segment on slice ToSliceNum for 
            Target.
        
        """
        
        #print(f'\nApplying Case 2a..\n')
        if UseCase == '1a':
            title = 'Applying Case 1a..'
        else:
            title = 'Applying Case 2a..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        # The new list of Target PerFrameFunctionalGroupsSequence-to-DICOM
        # slice index = the original list with ToSliceNum appended to it:
        NewTrgPFFGStoSliceInds = deepcopy(TrgPFFGStoSliceInds)
        NewTrgPFFGStoSliceInds.append(ToSliceNum)
        
        # Similarly for the new Target PerFrameFunctionalGroupsSequence-to-
        # DICOM slice index grouped by segment:
        NewTrgPFFGStoSliceIndsBySeg = deepcopy(TrgPFFGStoSliceIndsBySeg)
        NewTrgPFFGStoSliceIndsBySeg.append([ToSliceNum])
        
        
        # Shift the the in-plane elements in the Source frame that is to be 
        # copied:
        SrcPixArrToCopy = ShiftFrame(Frame=SrcPixArrToCopy, 
                                     PixShift=PixShift)
        
        
    
    if UseCase == '1b' or UseCase == '2b':
        """
        Main use case = 1 or 2:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            voxel spacings, and either are the same series (Main use case = 1) 
            or different series (Main use case = 2).
            
        Sub-case = b:
            Relationship-preserving copy of segmentation on slice FromSliceNum 
            in segment FromSegLabel for Source to a new segment on FromSliceNum 
            (since ToSliceNum = FromSliceNum) for Target.
        
        
        26/10: It's seems that having the same FOR does not guarantee the same
        origin.  Series 9 and 6 for Subject 011, MR4 have different origins!
        So need to check whether origins are the same?...
        
        17/11:  This case is up-to-date now but maintaining the same outputs
        for different cases is very messy and inefficient (e.g. see creation of
        unnecessary variables for Case 2b).
        
        18/11:  Added check if origins are the same in WhichUseCase().
        """
        
        #print('\nApplying Case 2b..\n')
        if UseCase == '2a':
            title = 'Applying Case 2a..'
        else:
            title = 'Applying Case 2b..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        if LogToConsole:
            #ToSliceNum = deepcopy(FromSliceNum)
            
            print(f'\nCopying slice {FromSliceNum} from Source SEG to slice',
                  #f'{ToSliceNum} in Target SEG..')
                  f'{FromSliceNum} in Target SEG..')
        
        # The new list of Target PerFrameFunctionalGroupsSequence-to-slice
        # index = the original list with FromSliceNum appended to it:
        NewTrgPFFGStoSliceInds = deepcopy(TrgPFFGStoSliceInds)
        NewTrgPFFGStoSliceInds.append(FromSliceNum)
        
        # Similarly for the new Target PerFrameFunctionalGroupsSequence-to-
        # slice index grouped by segment:
        NewTrgPFFGStoSliceIndsBySeg = deepcopy(TrgPFFGStoSliceIndsBySeg)
        NewTrgPFFGStoSliceIndsBySeg.append([FromSliceNum])
        
        
        
    if '1' in UseCase or '2' in UseCase:
        """
        What follows is common to UseCase 1a, 1b, 2a and 2b.
        """
        
        # The new Target segment numbers are:
        NewTrgSegNums = []
        
        print('')
        
        for SegNum in range(len(NewTrgPFFGStoSliceIndsBySeg)):
            # The number of frames for this segment number:
            Nframes = len(NewTrgPFFGStoSliceIndsBySeg[SegNum])
            
            print(f'SegNum = {SegNum}, Nframes = {Nframes}')
            
            NewTrgSegNums.extend([SegNum]*Nframes)
            
        
        # The new Target pixel array = original Target pixel array with
        # SrcPixArrToCopy added to it:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoSliceInds), 
                                 TrgPixArr.shape[1], 
                                 TrgPixArr.shape[2]))
        
        #print(f'\nTrgPixArr.shape = {TrgPixArr.shape}')
        #print(f'NewTrgPixArr.shape = {NewTrgPixArr.shape}')
        #print(f'SrcPixArrToCopy.shape = {SrcPixArrToCopy.shape}')
        
        NewTrgPixArr[:TrgPixArr.shape[0]] = TrgPixArr
        NewTrgPixArr[-1] = SrcPixArrToCopy
        
        
        if LogToConsole:
            print('\nTrgPFFGStoSliceInds    =', TrgPFFGStoSliceInds)
            print('\nNewTrgPFFGStoSliceInds =', NewTrgPFFGStoSliceInds)
            
        print(f'\nThere are {NewTrgPixArr.shape[0]} frames in NewTrgPixArr')
        print(f'NewTrgPFFGStoSliceInds: {NewTrgPFFGStoSliceInds}')
        print(f'NewTrgSegNums: {NewTrgSegNums}')
        
        if not TrgSegFpath:
            # Target does not have an existing SEG, so initialise a new SEG 
            # object:
            TrgSeg = InitialiseSeg(TrgDicoms, NewTrgPFFGStoSliceInds,
                                   TrgImage, SrcSeg, LogToConsole)
        
        # Modify select tags in NewTrgSeg:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoSliceInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole=False)
        
        
        
        """ Note:
            
        The following variables:
            
        ResSrcImage
        SrcLabmapIms
        ResSrcLabmapIms
        SrcLabmapToCopyIm
        ResSrcLabmapToCopyIm
        NewTrgLabmapIms
        
        were not used for this case but this main function outputs:
        
        return SrcImage, ResSrcImage, TrgImage,\
               SrcLabmapIms, ResSrcLabmapIms,\
               SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               NewTrgSeg
               
        So create:
            
        SrcLabmapToCopyIm
        NewTrgLabmapIms
        
        even though they're not needed for this case, and return 
        SrcLabmapToCopyIm for ResSrcLabmapToCopyIm, and SrcLabmapIms for
        ResSrcLabmapIms.
        """
        
        # This is already done prior to spitting to different use cases:
        #SrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
        #                                 FrameToSliceInds=[ToSliceNum], 
        #                                 RefIm=TrgImage)
        
        ResSrcLabmapIms = []
        
        
        ResSrcLabmapToCopyIm = PixArr2Image(PixArr=SrcPixArrToCopy,
                                            #FrameToSliceInds=F2SInds,
                                            FrameToSliceInds=[ToSliceNum], 
                                            RefIm=TrgImage)
        
        # Create NewTrgLabmapIms even though it's not needed:
        NewTrgLabmapIms = []
        
        if TrgSegFpath:
            # Append the existing labelmaps for Target:
            for TrgLabmapIm in TrgLabmapIms:
                NewTrgLabmapIms.append(TrgLabmapIm)
            
        # then append the Source labelmap to be copied:    
        NewTrgLabmapIms.append(SrcLabmapToCopyIm)
        
        
        return SrcImage, SrcImage, TrgImage,\
               SrcLabmapIms, SrcLabmapIms,\
               SrcLabmapToCopyIm, SrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               SrcSeg, TrgSeg, NewTrgSeg,\
               SrcDicoms, TrgDicoms
               
        
        """
        23/10:  Still need to verify that Case 2 provides the expected results.
        
        26/10:  Resampling will need to be performed if the origins are not the
        same... 19/11: Now WhichUseCase() checks the origins and if they are
        not the same it's considered UseCase 3.
        """
        
    
    
    
    if '3' in UseCase:
        """
        What follows is common to both UseCase 3a and 3b.
        """
        
        # Resample the Source image:
        """ This is not required but useful info """
        ResSrcImage = ResampleImage(Image=SrcImage, 
                                    RefImage=TrgImage, 
                                    Interpolation='Linear')
        
        if LogToConsole:                            
            print('Results after resampling Source image:\n')
            
            print('The Source and Target image size:\n',
                          f'   Original Source:  {SrcImage.GetSize()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSize()} \n',
                          f'   Target:           {TrgImage.GetSize()}')
                
            print('The Source and Target voxel spacings:\n',
                          f'   Original Source:  {SrcImage.GetSpacing()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSpacing()} \n',
                          f'   Target:           {TrgImage.GetSpacing()}')
                
            print('The Source and Target Origin:\n',
                          f'   Original Source:  {SrcImage.GetOrigin()} \n',
                          f'   Resampled Source: {ResSrcImage.GetOrigin()} \n',
                          f'   Target:           {TrgImage.GetOrigin()}')
                
            print('The Source and Target Direction:\n',
                          f'   Original Source:  {SrcImage.GetDirection()} \n',
                          f'   Resampled Source: {ResSrcImage.GetDirection()} \n',
                          f'   Target:           {TrgImage.GetDirection()}')
                                   
                                    
        # Resample the Source labelmap images:
        """ This is not required but useful info """
        ResSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            ResSrcLabmapIm = ResampleImage(Image=SrcLabmapIm, 
                                           RefImage=TrgImage, 
                                           Interpolation='NearestNeighbor')
            
            ResSrcLabmapIms.append(ResSrcLabmapIm)
        
        
        
        # Resample the Source labelmap to be copied:
        ResSrcLabmapToCopyIm = ResampleImage(Image=SrcLabmapToCopyIm, 
                                             RefImage=TrgImage, 
                                             Interpolation='NearestNeighbor')
        
        #nda = sitk.GetArrayFromImage(ResSrcLabmapToCopyIm)
        #unique = UniqueItems(Nda=nda, NonZero=False)
        #print(f'\nThere are {len(unique)} unique items in ResSrcLabmapToCopyIm')
        
        
    if UseCase == '3a':
        """
        Main use case = 3:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            but different voxel spacings.  
            
        Sub-case = a:
            Direct copy of the segmentation on slice FromSliceNum in segment 
            FromSegLabel for Source to a new segment on slice ToSliceNum for 
            Target.
        """
        
        #print(f'\nApplying Case 3a..\n')
        title = 'Applying Case 3a..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        # Convert the resampled Source labelmap-to-be-copied image to a Numpy
        # array:
        ResSrcPixArrToCopy,\
        F2SInds = Image2PixArr(LabmapIm=ResSrcLabmapToCopyIm)
        
        unique = UniqueItems(Nda=ResSrcPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResSrcPixArrToCopy')
        
        print(f'\nThe segmentation on FromSliceNum = {FromSliceNum} has been',
              f'resampled to {ResSrcPixArrToCopy.shape[0]} frames:',
              f'{F2SInds}')
        
        """
        For Direct copy 1 segment is always copied to 1 frame, so if after
        resampling there is more than one frame, either all but one of frames 
        will need to be ignored, or all frames will need to be averaged or a
        logical OR operation performed on all frames, etc.
        """
        if ResSrcPixArrToCopy.shape[0] > 1:
            print(f'\nResSrcPixArrToCopy has {ResSrcPixArrToCopy.shape[0]}',
                  'frames so the pixel arrays will be averaged...')
            
            # Take pixel-by-pixel mean across all frames in ResSrcPixArrToCopy:
            ResSrcPixArrToCopy = MeanPixArr(ResSrcPixArrToCopy)
            
            print(f'\nResSrcPixArrToCopy now has {ResSrcPixArrToCopy.shape[0]}',
                  'frames')
        else:
            print(f'\nResSrcPixArrToCopy has {ResSrcPixArrToCopy.shape[0]}',
                  'frames')
            
        
        # Shift the the in-plane elements in the resampled Source frame that is
        # to be copied:
        ResSrcPixArrToCopy = ShiftFrame(Frame=ResSrcPixArrToCopy, 
                                        PixShift=PixShift)
        
        print(f'\nAfter shifting the frame, ResSrcPixArrToCopy has',
              f'{ResSrcPixArrToCopy.shape[0]} frames')
        
        unique = UniqueItems(Nda=ResSrcPixArrToCopy, NonZero=False)
        print(f'\nThere are {len(unique)} unique items in ResSrcPixArrToCopy',
              'after shifting the frame.')
        
        # Convert the shifted resampled pixel array back to a SimpleITK image:
        """
        For Direct copy 1 segment is always copied to 1 frame, hence
        FrameToSliceInds=[ToSliceNum].
        """
        ResSrcLabmapToCopyIm = PixArr2Image(PixArr=ResSrcPixArrToCopy,
                                            #FrameToSliceInds=F2SInds,
                                            FrameToSliceInds=[ToSliceNum], 
                                            RefIm=TrgImage)
        
        
        
    if UseCase == '3b':
        """
        Main use case = 3:
            Source and Target have the same FrameOfReferenceUID, IOP, origins, 
            but different voxel spacings.  
            
        Sub-case = b:
            Relationship-preserving copy of segmentation on slice FromSliceNum 
            in segment FromSegLabel for Source to a new segment with any number 
            of segmentations on any number of slices for Target. 
        """
        
        #print(f'\nApplying Case 3b..')
        title = 'Applying Case 3a..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        
    if '3' in UseCase:
        """
        The following is common to UseCase 3a and 3b.
        """
        
        # Create a list of new Target Labelmap images including the existing 
        # labelmaps...:
        NewTrgLabmapIms = []
        
        if TrgSegFpath:
            # Append the existing labelmaps for Target:
            for TrgLabmapIm in TrgLabmapIms:
                NewTrgLabmapIms.append(TrgLabmapIm)
            
        # ... and append the resampled Source labelmap to be copied:    
        NewTrgLabmapIms.append(ResSrcLabmapToCopyIm)
        
        print(f'\nlen(NewTrgLabmapIms) = {len(NewTrgLabmapIms)}')
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrBySeg = []
        NewTrgPFFGStoSliceIndsBySeg = []
        NewTrgPFFGStoSliceInds = []
        NewTrgSegNumsBySeg = []
        NewTrgSegNums = []
        
        print('')
        
        for SegNum in range(len(NewTrgLabmapIms)):
            PixArr, F2SInds = Image2PixArr(LabmapIm=NewTrgLabmapIms[SegNum])
            
            unique = UniqueItems(Nda=PixArr, NonZero=False)

            print(f'   SegNum = {SegNum}, F2SInds = {F2SInds}')
            print(f'   There are {len(unique)} unique items in PixArr')
            
            NewTrgPixArrBySeg.append(PixArr)
            
            NewTrgPFFGStoSliceIndsBySeg.append(F2SInds)
            
            NewTrgPFFGStoSliceInds.extend(F2SInds)
            
            NewTrgSegNumsBySeg.append([SegNum]*len(F2SInds))
            
            NewTrgSegNums.extend([SegNum]*len(F2SInds))
            
        
        # Combine all frames:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoSliceInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for SegNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoSliceIndsBySeg[SegNum])):
                NewTrgPixArr[n] = NewTrgPixArrBySeg[SegNum][i]
                
                n += 1
        
        
        print(f'\nThe new Target pixel array has shape {NewTrgPixArr.shape}')
        print(f'\nNewTrgPFFGStoSliceInds: {NewTrgPFFGStoSliceInds}')
        print(f'NewTrgPFFGStoSliceIndsBySeg: {NewTrgPFFGStoSliceIndsBySeg}')
        print(f'\nNewTrgSegNums: {NewTrgSegNums}')
        print(f'NewTrgSegNumsBySeg: {NewTrgSegNumsBySeg}')
        
        
        unique = UniqueItems(Nda=NewTrgPixArr, NonZero=False)
        #unique = UniqueItems(Nda=NewTrgPixArr, NonZero=True)
        
        #print(f'\nUnique items in NewTrgPixArr = {unique}')
        print(f'\nThere are {len(unique)} unique items in NewTrgPixArr')
        
        if not TrgSegFpath:
            # Target does not have an existing SEG, so initialise a new SEG 
            # object:
            TrgSeg = InitialiseSeg(TrgDicoms, NewTrgPFFGStoSliceInds,
                                   TrgImage, SrcSeg, LogToConsole)
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoSliceInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole)
        
        
        if LogToConsole:
            print('\n\n')
            title = 'Results after resampling Source LabelMap image and ' \
                    + 'Source-to-copy LabelMap image:'
            ul = '*' * len(title)
            
            print('\n' + title)
            print(ul)
            
            print('\n*The Source and Target labelmap image size:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetSize()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetSize()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetSize()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSize()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSize()}')
                
            print('\n*The Source and Target labelmap image voxel spacings:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetSpacing()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetSpacing()} \n',
                          f'   Source-to-copy            {SrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSpacing()}')
                
            print('\n*The Source and Target labelmap image Origin:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetOrigin()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetOrigin()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Target:                   {TrgLabmapIm.GetOrigin()}')
                
            print('\n*The Source and Target labelmap image Direction:*\n',
                          #f'   Source:                   {SrcLabmapIm.GetDirection()} \n',
                          #f'   Resampled Source:         {ResSrcLabmapIm.GetDirection()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Target:                   {TrgLabmapIm.GetDirection()}')
        
        
        
        if LogToConsole:
            #SrcLabmapImMax = ImageMax(SrcLabmapIm)
            #SrcLabmapImMin = ImageMin(SrcLabmapIm)
            
            #ResSrcLabmapImMax = ImageMax(ResSrcLabmapIm)
            #ResSrcLabmapImMin = ImageMin(ResSrcLabmapIm)
            
            SrcLabmapToCopyImMax = ImageMax(SrcLabmapToCopyIm)
            SrcLabmapToCopyImMin = ImageMin(SrcLabmapToCopyIm)
            
            ResSrcLabmapToCopyImMax = ImageMax(ResSrcLabmapToCopyIm)
            ResSrcLabmapToCopyImMin = ImageMin(ResSrcLabmapToCopyIm)
            
            TrgLabmapImMax = ImageMax(TrgLabmapIm)
            TrgLabmapImMin = ImageMin(TrgLabmapIm)
            
            print('\n*The Source and Target labelmap image minimum:*\n',
                          #f'   Source:                   {SrcLabmapImMin} \n',
                          #f'   Resampled Source:         {ResSrcLabmapImMin} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMin} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMin} \n',
                          f'   Target:                   {TrgLabmapImMin}')
            
            print('\n*The Source and Target labelmap image maximum:*\n',
                          #f'   Source:                   {SrcLabmapImMax} \n',
                          #f'   Resampled Source:         {ResSrcLabmapImMax} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMax} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMax} \n',
                          f'   Target:                   {TrgLabmapImMax}')
            
        
        return SrcImage, ResSrcImage, TrgImage,\
               SrcLabmapIms, ResSrcLabmapIms,\
               SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               SrcSeg, TrgSeg, NewTrgSeg,\
               SrcDicoms, TrgDicoms
    
         
        
        
    if UseCase == '4':
        """
        Main use case = 4:
            Source and Target have the same FrameOfReferenceUID and origins, 
            but different IOP.  The voxel spacings may be the same or different.  
            
        Does UseCase 4 require sub-cases?...
        
        Don't have test data set to develop this yet.
        On 17/11 Simon said he would ask Jess for "data with different 
        angulations".
        """
        
        """ need to continue """
        
        #print(f'\nApplying Case 4..')
        title = 'Applying Case 4..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        
        
    if UseCase == '5':
        """
        Main use case = 5:
            Source and Target have the different FrameOfReferenceUID (hence
            different IOP and origins) and may have the same or different voxel 
            spacings.  
            
            Relationship-preserving copy of segmentation on slice FromSliceNum 
            in segment FromSegLabel for Source to a new segment with any number 
            of segmentations on any number of slices for Target. 
    
            
        Does UseCase 5 require sub-cases?...
        """
        
        
        #print(f'\nApplying Case 5..')
        title = 'Applying Case 5..'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
        
        # Register the images using an affine transformation:
        RegSrcImage, RegImFilt = RegisterImages(FixIm=TrgImage, MovIm=SrcImage, 
                                                Tx='affine',
                                                LogToConsole=LogToConsole)
        
        # Transform SrcLabmapToCopyIm:
        TxSrcLabmapToCopyIm = TransformImage(Im=SrcLabmapToCopyIm, 
                                             RegImFilt=RegImFilt)
        
        # Binary threshold TxSrcLabmapToCopyIm:
        TxSrcLabmapToCopyIm = BinaryThresholdImage(Im=TxSrcLabmapToCopyIm, 
                                                   Threshold=0.5)
        
        # Transform all Source labelmap images:
        """ This is not required but useful info """
        TxSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            TxSrcLabmapIm = TransformImage(Im=SrcLabmapIm, 
                                           RegImFilt=RegImFilt)
            
            TxSrcLabmapIm = BinaryThresholdImage(Im=TxSrcLabmapIm, 
                                                 Threshold=0.5)
            
            TxSrcLabmapIms.append(TxSrcLabmapIm)
        
        
        """
        24/11: The labelmaps in TxSrcLabmapIms look odd.
        """
        
        # Create a list of new Target Labelmap images including the existing 
        # labelmaps...:
        NewTrgLabmapIms = []
        
        if TrgSegFpath:
            # Append the existing labelmaps for Target:
            for TrgLabmapIm in TrgLabmapIms:
                NewTrgLabmapIms.append(TrgLabmapIm)
            
        # ... and append the transformed Source labelmap to be copied:    
        NewTrgLabmapIms.append(TxSrcLabmapToCopyIm)
        
        print(f'\nlen(NewTrgLabmapIms) = {len(NewTrgLabmapIms)}')
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrBySeg = []
        NewTrgPFFGStoSliceIndsBySeg = []
        NewTrgPFFGStoSliceInds = []
        NewTrgSegNumsBySeg = []
        NewTrgSegNums = []
        
        print('')
        
        for SegNum in range(len(NewTrgLabmapIms)):
            PixArr, F2SInds = Image2PixArr(LabmapIm=NewTrgLabmapIms[SegNum])
            
            unique = UniqueItems(Nda=PixArr, NonZero=False)

            print(f'   SegNum = {SegNum}, F2SInds = {F2SInds}')
            print(f'   There are {len(unique)} unique items in PixArr')
            
            NewTrgPixArrBySeg.append(PixArr)
            
            NewTrgPFFGStoSliceIndsBySeg.append(F2SInds)
            
            NewTrgPFFGStoSliceInds.extend(F2SInds)
            
            NewTrgSegNumsBySeg.append([SegNum]*len(F2SInds))
            
            NewTrgSegNums.extend([SegNum]*len(F2SInds))
            
        
        # Combine all frames:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoSliceInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for SegNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoSliceIndsBySeg[SegNum])):
                NewTrgPixArr[n] = NewTrgPixArrBySeg[SegNum][i]
                
                n += 1
        
        
        print(f'\nThe new Target pixel array has shape {NewTrgPixArr.shape}')
        print(f'\nNewTrgPFFGStoSliceInds: {NewTrgPFFGStoSliceInds}')
        print(f'NewTrgPFFGStoSliceIndsBySeg: {NewTrgPFFGStoSliceIndsBySeg}')
        print(f'\nNewTrgSegNums: {NewTrgSegNums}')
        print(f'NewTrgSegNumsBySeg: {NewTrgSegNumsBySeg}')
        
        
        unique = UniqueItems(Nda=NewTrgPixArr, NonZero=False)
        #unique = UniqueItems(Nda=NewTrgPixArr, NonZero=True)
        
        #print(f'\nUnique items in NewTrgPixArr = {unique}')
        print(f'\nThere are {len(unique)} unique items in NewTrgPixArr')
        
        if not TrgSegFpath:
            # Target does not have an existing SEG, so initialise a new SEG 
            # object:
            TrgSeg = InitialiseSeg(TrgDicoms, NewTrgPFFGStoSliceInds,
                                   TrgImage, SrcSeg, LogToConsole)
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgSegNums, NewTrgPFFGStoSliceInds, 
                                     FromSegNum, FromSliceNum, SrcSeg, 
                                     LogToConsole)
        
        #NewTrgSeg = deepcopy(TrgSeg)
        
        return SrcImage, RegSrcImage, TrgImage,\
               SrcLabmapIms, TxSrcLabmapIms,\
               SrcLabmapToCopyIm, TxSrcLabmapToCopyIm,\
               TrgLabmapIms, NewTrgLabmapIms,\
               SrcSeg, TrgSeg, NewTrgSeg,\
               SrcDicoms, TrgDicoms
    
        
        
    """
    The following was moved to each UseCase since some of the outputs don't
    exist for some cases (e.g. resampled Source image for Case 2b). """
    #return SrcImage, ResSrcImage, TrgImage,\
    #       SrcLabmapIms, ResSrcLabmapIms,\
    #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
    #       TrgLabmapIms, NewTrgLabmapIms,\
    #       NewTrgSeg






def CopyContour(SrcRts, FromSearchString, SrcDcmDir, FromSliceNum, TrgDcmDir,
                TrgRts=None, ToSliceNum=None, Sigma=(1,1,1), ThreshLevel=0.75,
                AddText='', LogToConsole=False):
    """
    Note 01/02/2021:
        This function (which can only copy a single contour) was formerly 
        called CopyRts. Now CopyRts will, depending on the inputs, copy either:
            1. A specific contour
            2. All contours in a specific ROI
            3. All contours in all ROIs
    
        
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
    from RtsTools import GetPtsByCntForSliceNum, ProportionOfContourInExtent
    from GeneralTools import GetPixelShiftBetweenSlices
    from GeneralTools import ShiftFrame
    from RtsTools import InitialiseRts
    from RtsTools import ModifyRts
    
    
    # Verify that an ROI exists whose name matches FromSearchString:
    #SrcRoiLabels = GetRoiLabels(SrcRts)
    FromRoiNum = GetRoiNum(Roi=SrcRts, SearchString=FromSearchString)
    
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
    
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                if LogToConsole:
                    print('\nBinary thresholding at',
                      f'{ThreshLevel}*{max(UniqueVals)} = {ThreshValue}...')
                
                ResLabmapImToCopy = BinaryThresholdImage(Im=ResLabmapImToCopy, 
                                                         Thresh=ThreshValue)
                
                if LogToConsole:
                    print(f'\nAfter binary thresholding:')
                
                
                # PixID required for if statement below.
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapImToCopy, LogToConsole)
        
        
        
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
        TrgRts = InitialiseRts(RtsTemplate=TrgRts, RoiNum=FromRoiNum,
                               CStoSliceInds=TrgCStoSliceInds,
                               CntDataByObjByFrame=TrgCntDataByObjByFrame,
                               DicomDir=TrgDcmDir, NamePrefix=AddText,
                               LogToConsole=LogToConsole)
    else:
        """ Use SrcRts to initialise the Target RTS. """
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
    TrgRts = ModifyRts(Rts=TrgRts, CntDataByObjByFrame=TrgCntDataByObjByFrame, 
                       PtsByObjByFrame=TrgPtsByObjByFrame, 
                       CStoSliceInds=TrgCStoSliceInds,
                       DicomDir=TrgDcmDir, LogToConsole=LogToConsole)
    
    
    return TrgRts








def CopySegmentation(SrcSeg, FromSearchString, SrcDcmDir, FromSliceNum, TrgDcmDir,
            TrgSeg=None, ToSliceNum=None, Sigma=(1,1,1), ThreshLevel=0.75,
            AddText='', LogToConsole=False):
                                
    """
    Note 08/02/2021:
        This function (which can only copy a single segmentation) was formerly 
        called CopySeg. Now CopySeg will, depending on the inputs, copy either:
            1. A specific segmentation
            2. All segmentations in a specific segment
            3. All segmentations in all segments
        
        
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






def GetInputsToCopyRoi(TestNum, UseOrigTrgRoi=False):
    import os
    from copy import deepcopy
    
    """ Subject TCGA-BB-A5HY in XNAT_TEST """
    
    Sub1_Dir = r"C:\Data\XNAT_TEST\TCGA-BB-A5HY"
    
    """ Series 2 """
    # IPP =
    # [1.0, 2.1e-10, 0.0,
    # -2.1e-10, 1.0, 0.0,
    # 0.0, 0.0, 1.0]


    Ses2_Ser2_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\2_Head Routine  5.0  H30s\DICOM")
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0

    Ses2_Ser4_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\4_Head Routine  5.0  H60s\DICOM")
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    Ses2_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\5_Head Routine  0.75  H20s\DICOM")
    # 512, 512, 330
    # 0.4, 0.4, 0.5
    # 0.75

    Ses2_Ser2_LeftEyeCT_RtsFpath = os.path.join(Sub1_Dir,
                                                r"Session2\assessors",
                                                r"RoiCollection_hkhGnxQ_2kX4w5b39s0\RTSTRUCT",
                                                r"AIM_20210201_000731.dcm") 
    # Left eye CT
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    Ses2_Ser2_RightEyeCT_RtsFpath = os.path.join(Sub1_Dir,
                                                 r"Session2\assessors",
                                                 r"RoiCollection_hkhGnxQ_2NZvk8B9pwA\RTSTRUCT",
                                                 r"AIM_20210201_000842.dcm") 
    # Right eye CT (single slice) <-- on slice 5 (0-indexed)
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    
    """ Series 4 """
    
    Ses4_Ser10_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans",
                                     r"10_T1  AXIAL 2MM POST repeat_S7_DIS3D\DICOM")
    # 192, 192, 80
    # 1.2, 1.2, 2.0 
    # 2.0
    
    Ses4_Ser11_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\11_T1 3D AX POST_S5_DIS3D\DICOM")
    # 512, 512, 192
    # 0.5, 0.5, 1.0 
    # 1.0

    Ses4_Ser12_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\12_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    # 204, 256, 57
    # 0.9, 0.9, 3.0 
    # 3.0

    Ses4_Ser13_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\13_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    # 192, 192, 80
    # 1.2, 1.2, 2.0 
    # 2.0

    Ses4_Ser14_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\14_T2_FLAIR 3MM_S2_DIS3D\DICOM")
    # 512, 512, 46
    # 0.45, 0.45, 3.9 
    # 3.0

    Ses4_Ser11_Ventricles_RtsFpath = os.path.join(Sub1_Dir,
                                                  r"Session4\assessors",
                                                  r"RoiCollection_hkhGnxQ_BNths5fNmA\RTSTRUCT",
                                                  r"AIM_20210131_230845.dcm") 
    # Ventricles
    # 512, 512, 192
    # 0.5, 0.5, 1.0 
    # 1.0
    
    Ses4_Ser13_NasalCavity_RtsFpath = os.path.join(Sub1_Dir,
                                                   r"Session4\assessors",
                                                   r"RoiCollection_hkhGnxQ_6daB7S8Eckp\RTSTRUCT",
                                                   r"AIM_20210131_214821.dcm") 
    # Nasal cavity
    # 192, 192, 80
    # 1.2, 1.2, 2.0 
    # 2.0
    
    """ Series 6"""
    # IPP =
    # [1.0, 2.1e-10, 0.0,
    # -2.1e-10, 1.0, 0.0,
    # 0.0, 0.0, 1.0]


    Ses6_Ser3_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session6\scans\3_Head Routine  5.0  H30f\DICOM")
    # 512, 512, 38
    # 0.43, 0.43, 5.0 
    # 5.0

    Ses6_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session6\scans\5_Head Routine  0.75  H20f\DICOM")
    # 512, 512, 375
    # 0.43, 0.43, 0.5 
    # 0.75
    
    
    
    
    """ Subject ACRIN-FMISO-Brain-011 in ACRIN-FMISO-Brain """
    
    Sub2_Dir = r"C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011"
    
    Sub2_RoiDir = r"C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011"
    Sub2_RtsDir = os.path.join(Sub2_RoiDir, "Fixed_RTS_Collections")
    Sub2_SegDir = os.path.join(Sub2_RoiDir, "Fixed_SEG_Collections")

    MR4_Dir = r"04-10-1960-MRI Brain wwo Contrast-69626"
    # IPP =
    # [1.00, 0.00, 0.07,
    # -0.02, 0.95, 0.30,
    # -0.07, -0.30, 0.95]
    
    MR12_Dir = r"06-11-1961-MRI Brain wwo Contrast-79433"
    # IPP =
    # [0.99, 0.11, 0.07,
    # -0.10, 0.99, -0.03,
    # -0.07, 0.02, 1.00]
    
    
    """ DICOM directories for MR4 """
    
    MR4_Ser3_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"3-T2 TSE AXIAL-07507") 
    # 378, 448, 21
    # 0.51, 0.51, 7.5 mm
    # 5.0 mm
    
    MR4_Ser4_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"4-T2 AXIAL FLAIR DARK FL-48830") 
    # 416, 512, 21 pix
    # 0.45, 0.45, 5.0 mm
    # 5.0 mm
    
    MR4_Ser5_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"5-T1 SE AXIAL 3MM-81246")        
    # 208, 256, 50 
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    MR4_Ser7_DcmDir = os.path.join(Sub2_Dir, MR4_Dir, 
                                   r"7-ep2ddiff3scantracep2ADC-22197") 
    # 192, 192, 21 pix
    # 1.3, 1.3, 7.5 mm
    # 5.0 mm
    
    MR4_Ser8_DcmDir = os.path.join(Sub2_Dir, MR4_Dir, 
                                   r"8-T1 SE AXIAL POST FS FC-59362") 
    # 212, 256, 30 pix
    # 0.9, 0.9, 5.0 mm
    # 5.0 mm
    
    MR4_Ser9_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"9-T1 SE AXIAL POST 3MM 5 MIN DELAY-07268") 
    # 208, 256, 50 pix
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    
    """ RTS/SEG filepaths for MR4 """
    
    MR4_Ser8_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    # tumour
    MR4_Ser8_SegFpath = os.path.join(Sub2_SegDir, 
                                     "Seg_from_MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    # tumour
    
    MR4_Ser9_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") 
    # tumour drawn on S9 stack
    #MR4_Ser9_RtsFpath = os.path.join(Sub2_RtsDir, 
    #                                 "MR4_S9_brain_RTS_AIM_20201112_142559.dcm")
    MR4_Ser9_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S9_tumour_and_brain_RTS_AIM_20201116_122858.dcm")
    
    #MR4_Ser9_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "Seg_from_MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") 
    # tumour drawn on S9 stack
    #MR4_Ser9_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S9_tumour_SEG_20201105_115408.dcm")
    MR4_Ser9_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S9_tumour_and_brain_SEG_20201105_115939.dcm")
    
    MR4_Ser5_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    # tumour drawn on S5 stack
    #MR4_Ser5_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "Seg_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    # tumour drawn on S5 stack
    MR4_Ser5_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S5_tumour_SEG_20201216_084441.dcm") 
    # single frame on slice 24 or 26
    
    #MR4_Ser3_RtsFpath = os.path.join(Sub2_RtsDir, 
    #                                 "MR4_S3_tumour_and_brain_RTS_AIM_20201112_162427.dcm") 
    # <-- WARNING: tumour & brain are not separate ROIs!
    MR4_Ser3_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm") 
    # separate ROIs
    
    MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S3_brain_SEG_20201030_133110.dcm")
    #MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") 
    # one ROI!
    #MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm") 
    # two ROIs
    MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") 
    # two ROIs
    #MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S3_tumour_SEG_20201216_085421.dcm") 
    # single frame on slice 10
    
    
    
    """ DICOM directories for MR12 """
    
    MR12_Ser3_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"3-T2 TSE AX IPAT2-08659") 
    # 336, 448, 24
    # 0.54, 0.54, 7.5 mm
    # 5.0 mm
    
    MR12_Ser4_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"4-T2 AXIAL FLAIR DARK FL-94212") 
    # 384, 512, 24
    # 0.47, 0.47, 7.5 mm
    # 5.0 mm
    
    MR12_Ser5_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"5-T1 SE AXIAL-43742") 
    # 192, 256, 35
    # 0.94, 0.94, 5.0 mm
    # 5.0 mm
    
    
    MR12_Ser8_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"8-T1 SE AXIAL POST FS FC-81428") 
    # 256, 192, 35 
    # 0.94, 0.94, 5.0 mm
    # 5.0 mm
    
    
    """ ******************************************************************* """
    """ TestNums for Testing                                                """
    """ ******************************************************************* """
    
    if TestNum == 'RR1':
        """ Rts Relationship-preserving copy test 1 """
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_LeftEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses2_Ser4_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Alternative inputs to test other possible combinations:
        #FromRoiLabel = None # copy all ROIs
        #FromSliceNum = 10 # there are no contours on slice 10
        #FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1 # direct copy to next slice
        #ToSliceNum = 0
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'RR2':
        """ Rts Relationship-preserving copy test 2 """
        
        SrcDcmDir = deepcopy(Ses4_Ser13_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser13_NasalCavity_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser10_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR3':
        """ Rts Relationship-preserving copy test 3 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_LeftEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses2_Ser5_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RR4':
        """ Rts Relationship-preserving copy test 4 """
        
        SrcDcmDir = deepcopy(Ses4_Ser13_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser13_NasalCavity_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'RR5':
        """ Rts Relationship-preserving copy test 5 """
        
        SrcDcmDir = deepcopy(Ses4_Ser13_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser13_NasalCavity_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser14_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'RR6':
        """ Rts Relationship-preserving copy test 6 """
        
        SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser10_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR7':
        """ Rts Relationship-preserving copy test 7 """
        
        SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser12_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'RR8':
        """ Rts Relationship-preserving copy test 8 """
        
        SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser13_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'RR9':
        """ Rts Relationship-preserving copy test 9 """
        
        SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses4_Ser14_DcmDir)
        TrgRoiFpath = None
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'RR10':
        """ Rts Relationship-preserving copy test 10 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_LeftEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser9_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' MR4 S9'
    
    elif TestNum == 'RR11':
        """ Rts Relationship-preserving copy test 11 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_LeftEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser7_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'MR4_Ser2'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' MR4 S7'
    
    elif TestNum == 'RR12':
        """ Rts Relationship-preserving copy test 12 """
        
        SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser8_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'MR4_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' MR4 S8'
        
    elif TestNum == 'RR13':
        """ Rts Relationship-preserving copy test 13 """
        
        SrcDcmDir = deepcopy(Ses4_Ser13_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser13_NasalCavity_RtsFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser4_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'MR4_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' MR4 S4'
    
    elif TestNum == 'RR14':
        """ Rts Relationship-preserving copy test 14 """
        
        SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
        
        TrgDcmDir = deepcopy(MR12_Ser8_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'MR12_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' MR12 S8'
        
    elif TestNum == 'RR15':
        """ Rts Relationship-preserving copy test 15 """
        
        SrcDcmDir = deepcopy(Ses4_Ser13_DcmDir)
        SrcRoiFpath = deepcopy(Ses4_Ser13_NasalCavity_RtsFpath)
        
        TrgDcmDir = deepcopy(MR12_Ser8_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'MR12_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' MR12 S4'
    
    #elif TestNum == 'RR16':
    # """ Rts Relationship-preserving copy test 16 """
    #
    #    SrcDcmDir = deepcopy(Ses4_Ser13_DcmDir)
    #    SrcRoiFpath = deepcopy(Ses4_Ser13_NasalCavity_RtsFpath)
    #    
    #    TrgDcmDir = deepcopy(MR12_Ser3_DcmDir)
    #    TrgDcmDir = deepcopy(MR12_Ser4_DcmDir)
    #    TrgDcmDir = deepcopy(MR12_Ser8_DcmDir)
    #    TrgRoiFpath = None
    #    
    #    # Inputs required for this test:
    #    FromRoiLabel = 'Nasal cavity'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' MR12 S3'
    #    TxtToAddToRoiLabel = ' MR12 S4'
    #    TxtToAddToRoiLabel = ' MR12 S8'
    
    #elif TestNum == 'RR16':
    # """ Rts Relationship-preserving copy test 16 """
    #
    #    SrcDcmDir = deepcopy(Ses4_Ser11_DcmDir)
    #    SrcRoiFpath = deepcopy(Ses4_Ser11_Ventricles_RtsFpath)
    #    
    #    TrgDcmDir = deepcopy(MR12_Ser3_DcmDir)
    #    TrgDcmDir = deepcopy(MR12_Ser4_DcmDir)
    #    TrgRoiFpath = None
    #    
    #    # Inputs required for this test:
    #    FromRoiLabel = 'Ventricles'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' MR12 S3'
    #    TxtToAddToRoiLabel = ' MR12 S4'
     
    if TestNum == 'RD1':
        """ Rts Direct copy test 1 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_RightEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 6
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'RD2':
        """ Rts Direct copy test 2 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_RightEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses2_Ser4_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'RD3':
        """ Rts Direct copy test 3 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_RightEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses2_Ser5_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RD4':
        """ Rts Direct copy test 4 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_RightEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses6_Ser3_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'RD5':
        """ Rts Direct copy test 5 """
        
        SrcDcmDir = deepcopy(Ses2_Ser2_DcmDir)
        SrcRoiFpath = deepcopy(Ses2_Ser2_RightEyeCT_RtsFpath)
        
        TrgDcmDir = deepcopy(Ses6_Ser5_DcmDir)
        TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_5'
        
    
    
    
        """ ************************************************************** """
        """ ************************************************************** """
        """ TestNums for Development                                       """
        """ ************************************************************** """
        """ ************************************************************** """
    
    elif TestNum in ['1R', '1S']:
        """ UseCase 1 (direct copy) Rts or Seg """
        
        SrcDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = 29
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
        
    elif TestNum in ['2aR', '2aS']:
        """ UseCase 2a (direct copy) Rts or Seg """
        
        SrcDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser5_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser5_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser5_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = 22
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
        
    elif TestNum in ['2bR', '2bS']:
        """ UseCase 2b (relationship-preserving copy) Rts or Seg """
        
        SrcDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser5_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser5_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser5_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 27
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    elif TestNum in ['3aiR', '3aiS']:
        """ UseCase 3ai (direct copy) Rts or Seg """
        
        SrcDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser3_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser3_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser3_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 25
        ToSliceNum = 20
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
    
    elif TestNum in ['3aiiR', '3aiiS']:
        """ UseCase 3aii (direct copy) Rts or Seg """
        
        SrcDcmDir = deepcopy(MR4_Ser3_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser3_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser3_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 11
        ToSliceNum = 26
        
        SrcLabel = 'MR4_Ser3'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
    
    elif TestNum in ['3biR', '3biS']:
        """ UseCase 3bi (relationship-preserving copy) Rts or Seg"""
        
        SrcDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser3_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser3_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser3_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 26
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    elif TestNum in ['3biiR', '3biiS']:
        """ UseCase 3bii (relationship-preserving copy) Rts or Seg """
        
        SrcDcmDir = deepcopy(MR4_Ser3_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser3_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser3_SegFpath)
        
        TrgDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
            else:
                TrgRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 11
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser3'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    elif TestNum in ['5bR', '5bS']:
        """ UseCase 5b (relationship-preserving copy) Rts or Seg"""
        
        SrcDcmDir = deepcopy(MR4_Ser9_DcmDir)
        if 'R' in TestNum:
            SrcRoiFpath = deepcopy(MR4_Ser9_RtsFpath)
        else:
            SrcRoiFpath = deepcopy(MR4_Ser9_SegFpath)
        
        TrgDcmDir = deepcopy(MR12_Ser8_DcmDir)
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                raise Exception('A RTS does not yet exist for MR12 Ser8.')
            else:
                raise Exception('A SEG does not yet exist for MR12 Ser8.')
        else:
            TrgRoiFpath = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR12_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    
    
    
    
    
    return SrcRoiFpath, FromSliceNum, FromRoiLabel, ToSliceNum, SrcDcmDir,\
           TrgDcmDir, TrgRoiFpath, TxtToAddToRoiLabel, SrcLabel, TrgLabel
        





def CopyRoi_OLD(SrcRoi, FromSearchString, SrcDcmDir, FromSliceNum, TrgDcmDir,
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
        """ If a single contour is to be copied FromSliceNum will be an integer.
        If FromSliceNum is None, all contours on all ROIs are to be copied
        (01/02/2021). """
        
        if isinstance(FromSliceNum, int):
            NewTrgRoi = CopyContour(SrcRoi, FromSearchString, SrcDcmDir, 
                                    FromSliceNum, TrgDcmDir, TrgRoi, ToSliceNum, 
                                    Sigma, ThreshLevel, AddText, LogToConsole)
        else:
            NewTrgRoi = CopyRts_OLD(SrcRoi, FromSearchString, SrcDcmDir, 
                                FromSliceNum, TrgDcmDir, TrgRoi, ToSliceNum, 
                                Sigma, ThreshLevel, AddText, LogToConsole)
    
    elif SrcModality == 'SEG':
        if isinstance(FromSliceNum, int):
            NewTrgRoi = CopySegmentation(SrcRoi, FromSearchString, SrcDcmDir, 
                                         FromSliceNum, TrgDcmDir, TrgRoi, 
                                         ToSliceNum, Sigma, ThreshLevel,
                                         AddText, LogToConsole)
        else:
            NewTrgRoi = CopySeg_OLD()
        
    else:
        msg = f'The Source modality ({SrcModality}) must be either "RTS" or '\
              + '"SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n\nDone.  Took {Dtime} s to run.')
        
    return NewTrgRoi






def CreateDictOfPaths_OLD():
    import os
    
    """ WARNING:
        Don't forget to add new directories/filepaths to the dictionary Dict
        below!
    """
    
    """ Subject TCGA-BB-A5HY in XNAT_TEST """
    
    Sub1_Dir = r"C:\Data\XNAT_TEST\TCGA-BB-A5HY"
    
    """ Session 1 """
    
    Ses1_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session1\scans\5_AX T2 BRAIN  PROPELLER\DICOM")
    
    Ses1_Ser6_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session1\scans\6_AX GRE BRAIN\DICOM")
    
    Ses1_Ser7_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session1\scans\7_AX T1 BRAIN\DICOM")
    
    Ses1_Ser8_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session1\scans\8_Ax T2 FLAIR BRAIN POST\DICOM")
    
    Ses1_Ser9_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session1\scans\9_AX T1 POST BRAIN\DICOM")
    
    
    """ Session 2 """
    # IPP =
    # [1.0, 2.1e-10, 0.0,
    # -2.1e-10, 1.0, 0.0,
    # 0.0, 0.0, 1.0]


    Ses2_Ser2_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\2_Head Routine  5.0  H30s\DICOM")
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    Ses2_Ser3_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\3_Head Routine  5.0  J30s  3\DICOM")
    
    Ses2_Ser4_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\4_Head Routine  5.0  H60s\DICOM")
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    Ses2_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session2\scans\5_Head Routine  0.75  H20s\DICOM")
    # 512, 512, 330
    # 0.4, 0.4, 0.5
    # 0.75

    Ses2_Ser2_LeftEyeCT_RtsFpath = os.path.join(Sub1_Dir,
                                                r"Session2\assessors",
                                                r"RoiCollection_hkhGnxQ_2kX4w5b39s0\RTSTRUCT",
                                                r"AIM_20210201_000731.dcm")
    
    Ses2_Ser2_LeftEyeCT_SegFpath = os.path.join(Sub1_Dir,
                                                r"Session2\assessors",
                                                r"RoiCollection_hlDxpyQ_7cJDZotLuP_icr_roiCollectionData\SEG",
                                                r"SEG_20210213_140215.dcm") 
    # Left eye CT
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    Ses2_Ser2_RightEyeCT_RtsFpath = os.path.join(Sub1_Dir,
                                                 r"Session2\assessors",
                                                 r"RoiCollection_hkhGnxQ_2NZvk8B9pwA\RTSTRUCT",
                                                 r"AIM_20210201_000842.dcm") 
    
    Ses2_Ser2_RightEyeCT_SegFpath = os.path.join(Sub1_Dir,
                                                 r"Session2\assessors",
                                                 r"RoiCollection_hlDxpyQ_AqEfVVmmysA_icr_roiCollectionData\SEG",
                                                 r"SEG_20210213_135857.dcm") 
    # Right eye CT (single slice) <-- on slice 5 (0-indexed)
    # 512, 512, 33
    # 0.4, 0.4, 5.0 
    # 5.0
    
    
    """ Session 3 """
    
    Ses3_Ser6_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session3\scans",
                                     r"6_T2 FLAIR Ax\DICOM")
    
    Ses3_Ser7_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session3\scans",
                                     r"7_Diff Ax\DICOM")
    
    Ses3_Ser13_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session3\scans",
                                     r"13_T1 Ax Post fs IAC\DICOM")
    
    Ses3_Ser14_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session3\scans",
                                     r"14_T1 Ax Post Brain\DICOM")
    
    
    """ Session 4 """
    
    Ses4_Ser10_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans",
                                     r"10_T1  AXIAL 2MM POST repeat_S7_DIS3D\DICOM")
    # 192, 192, 80
    # 1.2, 1.2, 2.0 
    # 2.0
    
    Ses4_Ser11_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\11_T1 3D AX POST_S5_DIS3D\DICOM")
    # 512, 512, 192
    # 0.5, 0.5, 1.0 
    # 1.0

    Ses4_Ser12_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\12_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    # 204, 256, 57
    # 0.9, 0.9, 3.0 
    # 3.0

    Ses4_Ser13_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\13_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    # 192, 192, 80
    # 1.2, 1.2, 2.0 
    # 2.0

    Ses4_Ser14_DcmDir = os.path.join(Sub1_Dir,
                                     r"Session4\scans\14_T2_FLAIR 3MM_S2_DIS3D\DICOM")
    # 512, 512, 46
    # 0.45, 0.45, 3.9 
    # 3.0

    Ses4_Ser11_Ventricles_RtsFpath = os.path.join(Sub1_Dir,
                                                  r"Session4\assessors",
                                                  r"RoiCollection_hkhGnxQ_BNths5fNmA\RTSTRUCT",
                                                  r"AIM_20210131_230845.dcm") 
    
    Ses4_Ser11_Ventricles_SegFpath = os.path.join(Sub1_Dir,
                                                  r"Session4\assessors",
                                                  r"RoiCollection_hlDxpyQ_5Dz9OulRH6n_icr_roiCollectionData\SEG",
                                                  r"SEG_20210213_144855.dcm") 
    # Ventricles
    # 512, 512, 192
    # 0.5, 0.5, 1.0 
    # 1.0
    
    Ses4_Ser13_NasalCavity_RtsFpath = os.path.join(Sub1_Dir,
                                                   r"Session4\assessors",
                                                   r"RoiCollection_hkhGnxQ_6daB7S8Eckp\RTSTRUCT",
                                                   r"AIM_20210131_214821.dcm") 
    
    Ses4_Ser13_NasalCavity_SegFpath = os.path.join(Sub1_Dir,
                                                   r"Session4\assessors",
                                                   r"RoiCollection_hlDxpyQ_AJ47Vd6ikeE_icr_roiCollectionData\SEG",
                                                   r"SEG_20210213_152014.dcm") 
    # Nasal cavity
    # 192, 192, 80
    # 1.2, 1.2, 2.0 
    # 2.0
    
    
    """ Session 5"""

    Ses5_Ser6_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session5\scans\6_T2 FLAIR Ax\DICOM")
    
    Ses5_Ser7_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session5\scans\7_CISS Ax\DICOM")
    
    Ses5_Ser12_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session5\scans\12_T2 Ax fs Post Brain\DICOM")
    
    Ses5_Ser14_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session5\scans\14_T1 Ax Post fs IAC\DICOM")
    
    Ses5_Ser16_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session5\scans\16_T1 Ax Post fs IAC repeat\DICOM")
    
    Ses5_Ser12_LeftEye_RtsFpath = os.path.join(Sub1_Dir,
                                               r"Session5\assessors",
                                               r"RoiCollection_hkhGnxQ_4rFYFKE1r6r\RTSTRUCT",
                                               r"AIM_20210131_212509.dcm")
    
    Ses5_Ser12_LeftEye_SegFpath = os.path.join(Sub1_Dir,
                                               r"Session5\assessors",
                                               r"RoiCollection_hlDxpyQ_8CVLpZmvRD8_icr_roiCollectionData\SEG",
                                               r"SEG_20210213_153745.dcm")
    
    Ses5_Ser16_Tumour_RtsFpath = os.path.join(Sub1_Dir,
                                               r"Session5\assessors",
                                               r"RoiCollection_hkhGnxQ_79YdSO5eYPo\RTSTRUCT",
                                               r"AIM_20210131_211622.dcm")
    
    Ses5_Ser16_Tumour_SegFpath = os.path.join(Sub1_Dir,
                                               r"Session5\assessors",
                                               r"RoiCollection_hlDxpyQ_18XInHEwVjJ_icr_roiCollectionData\SEG",
                                               r"SEG_20210213_155524.dcm")
    
    
    """ Session 6"""
    # IPP =
    # [1.0, 2.1e-10, 0.0,
    # -2.1e-10, 1.0, 0.0,
    # 0.0, 0.0, 1.0]


    Ses6_Ser3_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session6\scans\3_Head Routine  5.0  H30f\DICOM")
    # 512, 512, 38
    # 0.43, 0.43, 5.0 
    # 5.0
    
    Ses6_Ser4_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session6\scans\4_Head Routine  5.0  J30f  3\DICOM")
    
    Ses6_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session6\scans\5_Head Routine  0.75  H20f\DICOM")
    # 512, 512, 375
    # 0.43, 0.43, 0.5 
    # 0.75
    
    
    """ Session 7"""
    
    Ses7_Ser21_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session7\scans\21_T1  AXIAL 2MM POST_S6_DIS3D\DICOM")
    
    Ses7_Ser22_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session7\scans\22_T1 3D AX POST_S5_DIS3D\DICOM")
    
    Ses7_Ser23_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session7\scans\23_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    
    Ses7_Ser24_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session7\scans\24_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    
    Ses7_Ser25_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session7\scans\25_T2_FLAIR 3MM_S2_DIS3D\DICOM")
    
    
    """ Session 8"""
    
    Ses8_Ser2_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session8\scans\2_Head Routine  5.0  H30f\DICOM")
    
    Ses8_Ser4_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session8\scans\4_Head Routine  5.0  J30f  3\DICOM")
    
    Ses8_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session8\scans\5_Head Routine  0.75  H20f\DICOM")
    
    
    
    """ Session 9"""
    
    Ses9_Ser3_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session9\scans\3_T2 FLAIR Ax\DICOM")
    
    Ses9_Ser4_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session9\scans\4_Diff Ax\DICOM")
    
    Ses9_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session9\scans\5_Diff Ax_ADC\DICOM")
    
    Ses9_Ser14_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session9\scans\14_T2spc Ax\DICOM")
    
    Ses9_Ser15_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session9\scans\15_T1 3D Ax Post\DICOM")
    
    
    
    """ Session 10"""
    
    Ses10_Ser3_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session10\3\DICOM")
    
    Ses10_Ser6_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session10\6\DICOM")
    
    Ses10_Ser19_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session10\19\DICOM")
    
    Ses10_Ser22_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session10\22\DICOM")
    
    Ses10_Ser3_RightEye_RtsFpath = os.path.join(Sub1_Dir,
                                                r"Session10",
                                                r"AIM_20210131_183753\RTSTRUCT",
                                                r"AIM_20210131_183753.dcm")
    
    Ses10_Ser3_LowerBrain_RtsFpath = os.path.join(Sub1_Dir,
                                                  r"Session10",
                                                  r"AIM_20210131_232127\RTSTRUCT",
                                                  r"AIM_20210131_232127.dcm")
    
    
    
    """ Session 11"""
    
    Ses11_Ser2_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session11\scans\2_Head Routine  0.75  H20s\DICOM")
    
    Ses11_Ser3_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session11\scans\3_Head Routine  5.0  H30s\DICOM")
    
    Ses11_Ser5_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session11\scans\5_Head Routine  5.0  J30s  3\DICOM")
    
    Ses11_Ser9_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session11\scans\9_TEMPORAL BONES RIGHT  5.0  H30s\DICOM")
    
    Ses11_Ser10_DcmDir = os.path.join(Sub1_Dir,
                                    r"Session11\scans\10_TEMPORAL BONES RIGHT  0.6  H30s\DICOM")
    
    
    
    
    """ Subject ACRIN-FMISO-Brain-011 in ACRIN-FMISO-Brain """
    
    Sub2_Dir = r"C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011"
    
    Sub2_RoiDir = r"C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011"
    Sub2_RtsDir = os.path.join(Sub2_RoiDir, "Fixed_RTS_Collections")
    Sub2_SegDir = os.path.join(Sub2_RoiDir, "Fixed_SEG_Collections")

    MR4_Dir = r"04-10-1960-MRI Brain wwo Contrast-69626"
    # IPP =
    # [1.00, 0.00, 0.07,
    # -0.02, 0.95, 0.30,
    # -0.07, -0.30, 0.95]
    
    MR12_Dir = r"06-11-1961-MRI Brain wwo Contrast-79433"
    # IPP =
    # [0.99, 0.11, 0.07,
    # -0.10, 0.99, -0.03,
    # -0.07, 0.02, 1.00]
    
    
    """ DICOM directories for MR4 """
    
    MR4_Ser3_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"3-T2 TSE AXIAL-07507") 
    # 378, 448, 21
    # 0.51, 0.51, 7.5 mm
    # 5.0 mm
    
    MR4_Ser4_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"4-T2 AXIAL FLAIR DARK FL-48830") 
    # 416, 512, 21 pix
    # 0.45, 0.45, 5.0 mm
    # 5.0 mm
    
    MR4_Ser5_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"5-T1 SE AXIAL 3MM-81246")        
    # 208, 256, 50 
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    MR4_Ser7_DcmDir = os.path.join(Sub2_Dir, MR4_Dir, 
                                   r"7-ep2ddiff3scantracep2ADC-22197") 
    # 192, 192, 21 pix
    # 1.3, 1.3, 7.5 mm
    # 5.0 mm
    
    MR4_Ser8_DcmDir = os.path.join(Sub2_Dir, MR4_Dir, 
                                   r"8-T1 SE AXIAL POST FS FC-59362") 
    # 212, 256, 30 pix
    # 0.9, 0.9, 5.0 mm
    # 5.0 mm
    
    MR4_Ser9_DcmDir = os.path.join(Sub2_Dir, MR4_Dir,
                                   r"9-T1 SE AXIAL POST 3MM 5 MIN DELAY-07268") 
    # 208, 256, 50 pix
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    
    """ RTS/SEG filepaths for MR4 """
    
    MR4_Ser8_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    # tumour
    MR4_Ser8_SegFpath = os.path.join(Sub2_SegDir, 
                                     "Seg_from_MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    # tumour
    
    MR4_Ser9_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") 
    # tumour drawn on S9 stack
    #MR4_Ser9_RtsFpath = os.path.join(Sub2_RtsDir, 
    #                                 "MR4_S9_brain_RTS_AIM_20201112_142559.dcm")
    MR4_Ser9_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S9_tumour_and_brain_RTS_AIM_20201116_122858.dcm")
    
    #MR4_Ser9_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "Seg_from_MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") 
    # tumour drawn on S9 stack
    #MR4_Ser9_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S9_tumour_SEG_20201105_115408.dcm")
    MR4_Ser9_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S9_tumour_and_brain_SEG_20201105_115939.dcm")
    
    MR4_Ser5_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    # tumour drawn on S5 stack
    #MR4_Ser5_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "Seg_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    # tumour drawn on S5 stack
    MR4_Ser5_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S5_tumour_SEG_20201216_084441.dcm") 
    # single frame on slice 24 or 26
    
    #MR4_Ser3_RtsFpath = os.path.join(Sub2_RtsDir, 
    #                                 "MR4_S3_tumour_and_brain_RTS_AIM_20201112_162427.dcm") 
    # <-- WARNING: tumour & brain are not separate ROIs!
    MR4_Ser3_RtsFpath = os.path.join(Sub2_RtsDir, 
                                     "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm") 
    # separate ROIs
    
    MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S3_brain_SEG_20201030_133110.dcm")
    #MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") 
    # one ROI!
    #MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm") 
    # two ROIs
    MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
                                     "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") 
    # two ROIs
    #MR4_Ser3_SegFpath = os.path.join(Sub2_SegDir, 
    #                                 "MR4_S3_tumour_SEG_20201216_085421.dcm") 
    # single frame on slice 10
    
    
    
    """ DICOM directories for MR12 """
    
    MR12_Ser3_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"3-T2 TSE AX IPAT2-08659") 
    # 336, 448, 24
    # 0.54, 0.54, 7.5 mm
    # 5.0 mm
    
    MR12_Ser4_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"4-T2 AXIAL FLAIR DARK FL-94212") 
    # 384, 512, 24
    # 0.47, 0.47, 7.5 mm
    # 5.0 mm
    
    MR12_Ser5_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"5-T1 SE AXIAL-43742") 
    # 192, 256, 35
    # 0.94, 0.94, 5.0 mm
    # 5.0 mm
    
    
    MR12_Ser8_DcmDir = os.path.join(Sub2_Dir, MR12_Dir,
                                    r"8-T1 SE AXIAL POST FS FC-81428") 
    # 256, 192, 35 
    # 0.94, 0.94, 5.0 mm
    # 5.0 mm
    
    
    Dict = {}
    
    Dict['Sub1_Dir'] = Sub1_Dir
    
    
    Dict['Ses1_Ser5_DcmDir'] = Ses1_Ser5_DcmDir
    Dict['Ses1_Ser6_DcmDir'] = Ses1_Ser6_DcmDir
    Dict['Ses1_Ser7_DcmDir'] = Ses1_Ser7_DcmDir
    Dict['Ses1_Ser8_DcmDir'] = Ses1_Ser8_DcmDir
    Dict['Ses1_Ser9_DcmDir'] = Ses1_Ser9_DcmDir
    
    
    Dict['Ses2_Ser2_DcmDir'] = Ses2_Ser2_DcmDir
    Dict['Ses2_Ser3_DcmDir'] = Ses2_Ser3_DcmDir
    Dict['Ses2_Ser4_DcmDir'] = Ses2_Ser4_DcmDir
    Dict['Ses2_Ser5_DcmDir'] = Ses2_Ser5_DcmDir
    
    Dict['Ses2_Ser2_LeftEyeCT_RtsFpath'] = Ses2_Ser2_LeftEyeCT_RtsFpath
    Dict['Ses2_Ser2_RightEyeCT_RtsFpath'] = Ses2_Ser2_RightEyeCT_RtsFpath
    Dict['Ses2_Ser2_LeftEyeCT_SegFpath'] = Ses2_Ser2_LeftEyeCT_SegFpath
    Dict['Ses2_Ser2_RightEyeCT_SegFpath'] = Ses2_Ser2_RightEyeCT_SegFpath
    
    
    Dict['Ses3_Ser6_DcmDir'] = Ses3_Ser6_DcmDir
    Dict['Ses3_Ser7_DcmDir'] = Ses3_Ser7_DcmDir
    Dict['Ses3_Ser13_DcmDir'] = Ses3_Ser13_DcmDir
    Dict['Ses3_Ser14_DcmDir'] = Ses3_Ser14_DcmDir
    
    
    Dict['Ses4_Ser10_DcmDir'] = Ses4_Ser10_DcmDir
    Dict['Ses4_Ser11_DcmDir'] = Ses4_Ser11_DcmDir
    Dict['Ses4_Ser12_DcmDir'] = Ses4_Ser12_DcmDir
    Dict['Ses4_Ser13_DcmDir'] = Ses4_Ser13_DcmDir
    Dict['Ses4_Ser14_DcmDir'] = Ses4_Ser14_DcmDir
    
    Dict['Ses4_Ser11_Ventricles_RtsFpath'] = Ses4_Ser11_Ventricles_RtsFpath
    Dict['Ses4_Ser13_NasalCavity_RtsFpath'] = Ses4_Ser13_NasalCavity_RtsFpath
    Dict['Ses4_Ser11_Ventricles_SegFpath'] = Ses4_Ser11_Ventricles_SegFpath
    Dict['Ses4_Ser13_NasalCavity_SegFpath'] = Ses4_Ser13_NasalCavity_SegFpath
    
    
    Dict['Ses5_Ser6_DcmDir'] = Ses5_Ser6_DcmDir
    Dict['Ses5_Ser7_DcmDir'] = Ses5_Ser7_DcmDir
    Dict['Ses5_Ser12_DcmDir'] = Ses5_Ser12_DcmDir
    Dict['Ses5_Ser14_DcmDir'] = Ses5_Ser14_DcmDir
    Dict['Ses5_Ser16_DcmDir'] = Ses5_Ser16_DcmDir
    
    Dict['Ses5_Ser12_LeftEye_RtsFpath'] = Ses5_Ser12_LeftEye_RtsFpath
    Dict['Ses5_Ser16_Tumour_RtsFpath'] = Ses5_Ser16_Tumour_RtsFpath
    Dict['Ses5_Ser12_LeftEye_SegFpath'] = Ses5_Ser12_LeftEye_SegFpath
    Dict['Ses5_Ser16_Tumour_SegFpath'] = Ses5_Ser16_Tumour_SegFpath
    
    
    Dict['Ses6_Ser3_DcmDir'] = Ses6_Ser3_DcmDir
    Dict['Ses6_Ser4_DcmDir'] = Ses6_Ser4_DcmDir
    Dict['Ses6_Ser5_DcmDir'] = Ses6_Ser5_DcmDir
    
    
    Dict['Ses7_Ser21_DcmDir'] = Ses7_Ser21_DcmDir
    Dict['Ses7_Ser22_DcmDir'] = Ses7_Ser22_DcmDir
    Dict['Ses7_Ser23_DcmDir'] = Ses7_Ser23_DcmDir
    Dict['Ses7_Ser24_DcmDir'] = Ses7_Ser24_DcmDir
    Dict['Ses7_Ser25_DcmDir'] = Ses7_Ser25_DcmDir
    
    
    Dict['Ses8_Ser2_DcmDir'] = Ses8_Ser2_DcmDir
    Dict['Ses8_Ser4_DcmDir'] = Ses8_Ser4_DcmDir
    Dict['Ses8_Ser5_DcmDir'] = Ses8_Ser5_DcmDir
    
    
    Dict['Ses9_Ser3_DcmDir'] = Ses9_Ser3_DcmDir
    Dict['Ses9_Ser4_DcmDir'] = Ses9_Ser4_DcmDir
    Dict['Ses9_Ser5_DcmDir'] = Ses9_Ser5_DcmDir
    Dict['Ses9_Ser14_DcmDir'] = Ses9_Ser14_DcmDir
    Dict['Ses9_Ser15_DcmDir'] = Ses9_Ser15_DcmDir
    
    
    Dict['Ses10_Ser3_DcmDir'] = Ses10_Ser3_DcmDir
    Dict['Ses10_Ser6_DcmDir'] = Ses10_Ser6_DcmDir
    Dict['Ses10_Ser19_DcmDir'] = Ses10_Ser19_DcmDir
    Dict['Ses10_Ser22_DcmDir'] = Ses10_Ser22_DcmDir
    
    Dict['Ses10_Ser3_RightEye_RtsFpath'] = Ses10_Ser3_RightEye_RtsFpath
    Dict['Ses10_Ser3_LowerBrain_RtsFpath'] = Ses10_Ser3_LowerBrain_RtsFpath
    
    
    Dict['Ses11_Ser2_DcmDir'] = Ses11_Ser2_DcmDir
    Dict['Ses11_Ser3_DcmDir'] = Ses11_Ser3_DcmDir
    Dict['Ses11_Ser5_DcmDir'] = Ses11_Ser5_DcmDir
    Dict['Ses11_Ser9_DcmDir'] = Ses11_Ser9_DcmDir
    Dict['Ses11_Ser10_DcmDir'] = Ses11_Ser10_DcmDir
    
    
    
    Dict['Sub2_Dir'] = Sub2_Dir
    Dict['Sub2_RoiDir'] = Sub2_RoiDir
    Dict['Sub2_RtsDir'] = Sub2_RtsDir
    Dict['Sub2_SegDir'] = Sub2_SegDir
    
    Dict['MR4_Dir'] = MR4_Dir
    
    Dict['MR4_Ser3_DcmDir'] = MR4_Ser3_DcmDir
    Dict['MR4_Ser4_DcmDir'] = MR4_Ser4_DcmDir
    Dict['MR4_Ser5_DcmDir'] = MR4_Ser5_DcmDir
    Dict['MR4_Ser7_DcmDir'] = MR4_Ser7_DcmDir
    Dict['MR4_Ser8_DcmDir'] = MR4_Ser8_DcmDir
    Dict['MR4_Ser9_DcmDir'] = MR4_Ser9_DcmDir
    
    Dict['MR4_Ser3_RtsFpath'] = MR4_Ser3_RtsFpath
    Dict['MR4_Ser3_SegFpath'] = MR4_Ser3_SegFpath
    Dict['MR4_Ser5_RtsFpath'] = MR4_Ser5_RtsFpath
    Dict['MR4_Ser5_SegFpath'] = MR4_Ser5_SegFpath
    Dict['MR4_Ser8_RtsFpath'] = MR4_Ser8_RtsFpath
    Dict['MR4_Ser8_SegFpath'] = MR4_Ser8_SegFpath
    Dict['MR4_Ser9_RtsFpath'] = MR4_Ser9_RtsFpath
    Dict['MR4_Ser9_SegFpath'] = MR4_Ser9_SegFpath
    
    Dict['MR12_Dir'] = MR12_Dir
    
    Dict['MR12_Ser3_DcmDir'] = MR12_Ser3_DcmDir
    Dict['MR12_Ser4_DcmDir'] = MR12_Ser4_DcmDir
    Dict['MR12_Ser5_DcmDir'] = MR12_Ser5_DcmDir
    Dict['MR12_Ser8_DcmDir'] = MR12_Ser8_DcmDir
    
    """
    Dict[''] = 
    Dict[''] = 
    Dict[''] = 
    Dict[''] = 
    Dict[''] = 
    """
    
    return Dict







def CheckValidityOfInputs_OLD(SrcRoiName, SrcSliceNum, TrgRoiName, TrgSliceNum, 
                          LogToConsole=False):
    """
    Establish whether the inputs SrcRoiName, SrcSliceNum, TrgRoiName and
    TrgSliceNum are valid combinations. If not raise exception.
    
    Inputs:
    ******
        
    SrcRoiName : string
        The Source ROIName or SegmentLabel.
    
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
    
    TrgRoiName : string
        The Target ROIName or SegmentLabel.
        
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
            
            
    if SrcRoiName == None:
        msg = f'WARNING: SrcRoiName = {SrcRoiName}, so all contours/'\
              + 'segmentations in all ROIs/segments in the RTS/SEG are to'\
              + f' be copied. \nBut '
                  
        if SrcSliceNum != None and TrgSliceNum != None:
            msg += f'SrcSliceNum (= {SrcSliceNum} != None) and TrgSliceNum '\
                   + f'(= {TrgSliceNum} != None) are defined. '
            msg += 'This is not a valid combination of inputs.'
            
            raise Exception(msg)
            
        elif SrcSliceNum != None:
            msg += f'SrcSliceNum (= {SrcSliceNum} != None) is defined.'
            msg += 'This is not a valid combination of inputs.'
            
            raise Exception(msg)
            
        elif TrgSliceNum != None:
            msg += f'TrgSliceNum (= {TrgSliceNum} != None) is defined.'
            msg += 'This is not a valid combination of inputs.'
        
            raise Exception(msg)
    
    if SrcSliceNum == None and TrgSliceNum != None:
        msg = f'WARNING: SrcSliceNum = {SrcSliceNum}, so all contours/'\
              + 'segmentations in the ROI/segment whose label matches '\
              + f'{SrcRoiName} are to be copied. \nBut TrgSliceNum (= '\
              + f'{TrgSliceNum} != None) is defined. This is not a valid '\
              + 'combination of inputs.'
              
        raise Exception(msg)
    
    return True






def RunCopyXnatRoi_INCOMPLETE(TestNums, LogToConsole=False, XnatUrl=None, ProjId=None, 
                   SubjLabel=None, SrcStudyLabel=None, SrcSeriesLabel=None, 
                   SrcAsrMod=None, SrcAsrName=None, SrcSliceNum=None, 
                   TrgStudyLabel=None, TrgSeriesLabel=None, TrgAsrMod=None, 
                   TrgAsrName=None, TrgSliceNum=None, TxtToAddToTrgAsrName='', 
                   XnatSession=None, PathsDict=None, XnatDownloadDir='default', 
                   ExportLogFiles=False):
    """
    
    Inputs:
    ******
    
    TestNums : list of strings or None
        Either a list of strings denoting the tests to run (as defined in 
        CopyRoiTestConfig.py), or None. Several tests can be run in series by 
        providing a list of comma-separated strings, e.g. ['SR1', 'SR2', 'SR3'].
        To run a single test a single-itemed list must be provided, e.g. ['RD4'].
        To run the algorithm on a new set or different combination of inputs,
        either update CopyRoiTestConfig.py or set TestNums = None and provide
        the necessary optional inputs as described below.
    
    LogToConsole : boolean (optional if TestNums != None; False by default)
        If True, intermediate results will be logged to the console during the
        running of CopyRoi().
    
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
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
        
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
        
    
    Notes:
    *****
    
    If running one or more pre-configured tests, CopyRoiTestConfig.py must be
    copied to the current working directory.
    """

    import time
    import os
    from pypref import Preferences
    from pydicom import dcmread
    #from DicomTools import IsSameModalities
    from GeneralTools import PrintTitle, ExportListToTxt, ExportDictionaryToJson
    
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
    
    
    #if TestNums != None:
    
    for TestNum in TestNums:
        TestRunTimes = [] # this this Test run only
        TestRunTimes.append(time.time())
        
        PrintTitle(f'Running Test {TestNum}:')
        
        """ Get the pre-configured inputs for this test: """
        TestInputs = Preferences(directory=os.getcwd(), 
                                 filename='CopyRoiTestConfig.py')
        
        XnatUrl = TestCfg.get('XnatUrl')
        ProjId = TestCfg.get('ProjId')
        SubjLabel = TestCfg.get('SubjLabel')
        SrcStudyLabel = TestCfg.get(TestNum)['SrcStudyLabel']
        SrcSeriesLabel = TestCfg.get(TestNum)['SrcSeriesLabel']
        SrcAsrMod = TestCfg.get(TestNum)['SrcAsrMod']
        SrcAsrName = TestCfg.get(TestNum)['SrcAsrName']
        SrcSliceNum = TestCfg.get(TestNum)['SrcSliceNum']
        TrgSeriesLabel = TrgStudyLabel = TestCfg.get(TestNum)['TrgStudyLabel']
        TrgSeriesLabel = TestCfg.get(TestNum)['TrgSeriesLabel']
        TrgAsrMod = TestCfg.get(TestNum)['TrgAsrMod']
        TrgAsrName = TestCfg.get(TestNum)['TrgAsrName']
        TrgSliceNum = TestCfg.get(TestNum)['TrgSliceNum']
        TxtToAddToTrgAsrName = TestCfg.get(TestNum)['TxtToAddToTrgAsrName']
    
        """ Get the pre-defined defaults for remaining required inputs: """
        Defaults = Preferences(directory=os.getcwd(), 
                               filename='CopyRoiDefaults.py')
    
        ResInterp = Defaults.get('ResInterp')
        PreResVar = Defaults.get('PreResVar')
        ForceReg = Defaults.get('ForceReg')
        Tx = Defaults.get('Tx')
        TxMaxIters = Defaults.get('TxMaxIters')
        TxInterp = Defaults.get('TxInterp'),
        ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
        PostTxVar = Defaults.get('PostTxVar')
        ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
        
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
        
        """ New assessor name (StructureSetLabel or SeriesDescription) and 
        RTS/SEG filename: """
        if 'RR' in TestNum or 'RD' in TestNum:
            NewAsrName = SrcRoi.StructureSetLabel + TxtToAddTrgRoiName
        else:
            """ 'SR' in TestNum or 'SD' in TestNum: """
            NewAsrName = SrcRoi.SeriesDescription + TxtToAddTrgRoiName
            
        
        """ Text to add to file names: """
        TxtToAddToFname = f'Test_{TestNum}_{NewAsrName}'
        
        
        """ Copy Contours/Segmentations """
        
        NewTrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
        = CopyRoi(SrcRoiFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir, 
                  TrgRoiFpath, TrgSliceNum, ResInterp, PreResVar, #PostResThresh, 
                  ForceReg, Tx, TxMaxIters, TxInterp, 
                  ApplyPostTxBlur, PostTxVar,
                  ApplyPostTxBin, #ThreshPostTx, 
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



