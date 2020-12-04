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