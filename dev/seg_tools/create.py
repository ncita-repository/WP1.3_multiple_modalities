# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 13:47:10 2021

@author: ctorti
"""



def create_seg(SrcSegFpath, TrgSegFpath, TrgSegLabel, TrgPixArrBySeg, 
               TrgF2SindsBySeg, TrgDicomDir, SrcSegLabel, AddTxtToSegLabel='', 
               LogToConsole=False):
    """
    Create an SEG object for the target dataset.
    
    Parameters
    ----------
    SrcSegFpath : str
        Filepath of the Source SEG file. If TrgSegFpath = None, the Source SEG
        will be used as a template for the new Target SEG (NewTrgSeg).
    TrgSegFpath : str or None
        If TrgSegFpath != None, the Target SEG will be modified to create 
        NewTrgSeg.  If TrgSegFpath = None, the Source SEG will be used as a 
        template.
    TrgSegLabel : str or None
        All or part of the SegmentLabel of the destination segment.
    TrgPixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment for Target.  The pixel array is a sub-array of 
        the (entire) SEG's pixel array. 
    TrgF2SindsBySeg : list of a list of ints
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in TrgPixArrBySeg.
    TrgDicomDir : str
        Directory containing the Target DICOMs.
    SrcSegLabel : str
        All or part of the Source SegmentLabel of the segment containing the 
        segmentations(s) that were copied.
    AddTxtToSegLabel : str, optional ('' by default)
        String to be added to SeriesDescription for the new SEG. 
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
        
    Returns
    -------   
    NewTrgSeg : Pydicom Object
        New SEG object.
    
    Note
    ----
    Performance of adding/removing sequences using append()/pop():
    
    - Took 168.4 ms to add 10000 sequences to ReferencedInstanceSequence 
    (16.8 us/sequence)
    - Took 81.7 ms to remove 10000 sequences to ReferencedInstanceSequence 
    (8.2 us/sequence)
    
    - Took 92.4 ms to add 10000 sequences to PerFrameFunctionalGroupsSequence 
    (9.2 us/sequence)
    -Took 45.8 ms to remove 10000 sequences to PerFrameFunctionalGroupsSequence 
    (4.6 us/sequence)
    
    19/02:
        After avoiding the use of deepcopy to increase sequences in RtsSeg()
        I've had to re-instate them since the tag values were being duplicated 
        despite calls to modify them.  For some reason the removal of deepcopy
        from CreateSeg() hasn't resulted in the unexpected behaviours found 
        with CreateRts().  Nevertheless deepcopy was re-instated in CreateSeg()
        as a precaution.  Oddly the re-introduction of deepcopy doesn't seem to
        have lengthened the time to run CreateSeg() or to perform error 
        checking.
    """
    
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    #import time
    import datetime
    from copy import deepcopy
    import numpy as np
    from general_tools.general import flatten_list, get_unique_items
    from general_tools.general import reduce_list_of_str_floats_to_16
    from dicom_tools.imports import import_dicoms
    from dicom_tools.metadata import get_roicol_labels
    from image_tools.attrs_info import get_im_attrs
    from seg_tools.metadata import get_RSOP_uids_in_RIS
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of CreateSeg():')
        print('\n\n', '-'*120)
        print(f'   TrgF2SindsBySeg = {TrgF2SindsBySeg}')
        print(f'   len(TrgPixArrBySeg) = {len(TrgPixArrBySeg)}')
    
    SrcSeg = dcmread(SrcSegFpath)
    
    """ Use TrgSeg or SrcSeg as a template for NewTrgSeg. """
    if TrgSegFpath:
        NewTrgSeg = dcmread(TrgSegFpath)
    else:
        NewTrgSeg = deepcopy(SrcSeg)
    
    #print(f'\nTrgDicomDir = {TrgDicomDir}')
    
    TrgDicoms = import_dicoms(TrgDicomDir)
    
    Size, Spacings, ST, IPPs, Dirs, ListOfWarnings\
        = get_im_attrs(DicomDir=TrgDicomDir, Package='sitk')
    
    NumOfFramesBySeg = []
    for s in range(len(TrgPixArrBySeg)):
        NumOfFramesBySeg.append(TrgPixArrBySeg[s].shape[0])
    
    #print(f'\n\n\nNumOfFramesBySeg = {NumOfFramesBySeg}')
    
    NumOfFrames = sum(NumOfFramesBySeg)
    
    #NumOfSegs = len(TrgF2SindsBySeg)
    
    """ The list of segment labels in SrcSeg: """
    AllSrcSegLabels = get_roicol_labels(SrcSeg)
    
    if SrcSegLabel:
        """ The list of segment numbers for the segments of interest: """
        SrcSegNums = get_roicol_labels(SrcSeg, SrcSegLabel)
    else:
        SrcSegNums = list(range(len(AllSrcSegLabels)))
    
    UniqueTrgF2Sinds = get_unique_items(Items=TrgF2SindsBySeg, IgnoreZero=False, 
                                        MaintainOrder=True)
    
    FlattenedTrgF2Sinds = flatten_list(TrgF2SindsBySeg)
    
    if LogToConsole:
        print(f'   SrcSegNums = {SrcSegNums}')
        print(f'   NumOfFramesBySeg = {NumOfFramesBySeg}')
        print(f'   NumOfFrames = {NumOfFrames}')
    
    """ Generate a new SOPInstanceUID. """
    NewSOPuid = generate_uid()
    NewTrgSeg.SOPInstanceUID = deepcopy(NewSOPuid)
    NewTrgSeg.file_meta.MediaStorageSOPInstanceUID = deepcopy(NewSOPuid)
    
    """ Generate a new SeriesInstanceUID. """
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    #NewDate = time.strftime("%Y%m%d", time.gmtime())
    #NewTime = time.strftime("%H%M%S", time.gmtime())
    TimeNow = datetime.datetime.now()
    NewDate = TimeNow.strftime('%Y%m%d')
    NewTime = TimeNow.strftime('%H%M%S.%f')
    
    """ If TrgSegFpath != None, some tags will not need to be replaced.
    If TrgSegFpath = None, use the corresponding values in the first Target 
    DICOM. """
    
    if TrgSegFpath == None:
        NewTrgSeg.StudyDate = TrgDicoms[0].StudyDate
        try:
            NewTrgSeg.SeriesDate = TrgDicoms[0].SeriesDate
        except AttributeError:
            pass
        #NewTrgSeg.ContentDate = TrgDicoms[0].ContentDate
        NewTrgSeg.ContentDate = NewDate
        NewTrgSeg.StudyTime = TrgDicoms[0].StudyTime
        try:
            NewTrgSeg.SeriesTime = TrgDicoms[0].SeriesTime
        except AttributeError:
            pass
        #NewTrgSeg.ContentTime = TrgDicoms[0].ContentTime
        NewTrgSeg.ContentTime = NewTime
        NewTrgSeg.Manufacturer = TrgDicoms[0].Manufacturer
        try:
            NewTrgSeg.SeriesDescription += AddTxtToSegLabel
        except AttributeError:
            NewTrgSeg.SeriesDescription = AddTxtToSegLabel
        NewTrgSeg.PatientName = TrgDicoms[0].PatientName
        NewTrgSeg.PatientID = TrgDicoms[0].PatientID
        NewTrgSeg.PatientBirthDate = TrgDicoms[0].PatientBirthDate
        NewTrgSeg.PatientSex = TrgDicoms[0].PatientSex
        try:
            NewTrgSeg.PatientAge = TrgDicoms[0].PatientAge
        except AttributeError:
            pass
        NewTrgSeg.StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
        NewTrgSeg.StudyID = TrgDicoms[0].StudyID
        NewTrgSeg.SeriesNumber = TrgDicoms[0].SeriesNumber
        NewTrgSeg.InstanceNumber = TrgDicoms[0].InstanceNumber
        NewTrgSeg.FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
        NewTrgSeg.PositionReferenceIndicator = TrgDicoms[0].PositionReferenceIndicator
        NewTrgSeg.Rows = TrgDicoms[0].Rows
        NewTrgSeg.Columns = TrgDicoms[0].Columns
        
        IOP = TrgDicoms[0].ImageOrientationPatient
        
        # Ensure that the character length is limited to 16:
        IOP = reduce_list_of_str_floats_to_16(IOP)
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PlaneOrientationSequence[0]\
                 .ImageOrientationPatient = IOP
           
        """ 
        Notes:
        1. SliceThickness in PixelMeasuresSequence appears to be the z-Spacing 
        rather than SliceThickness from the DICOM metadata. 
        2. The character limit of VR DS is 16 characters.
        """
        Zspacing = f"{Spacings[2]}"
        
        #if len(Zspacing) > 16:
        #    Zspacing = Zspacing[:16]
        
        Zspacing = reduce_list_of_str_floats_to_16([Zspacing])[0]
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .SliceThickness = deepcopy(Zspacing)
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .SpacingBetweenSlices = deepcopy(Zspacing)
    
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .PixelSpacing = TrgDicoms[0].PixelSpacing
    
    """ Modify NumberOfFrames. """
    
    NewTrgSeg.NumberOfFrames = f"{len(FlattenedTrgF2Sinds)}"
    
    """ Modify ReferencedSeriesSequence. """
    
    NewTrgSeg.ReferencedSeriesSequence[0]\
             .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    """ Modify ReferencedInstanceSequence in ReferencedSeriesSequence. """
        
    for i in range(len(UniqueTrgF2Sinds)):
        """ The slice index s determines the SOPInstanceUID to reference. """
        s = UniqueTrgF2Sinds[i]
        
        """ The number of sequences in ReferencedSeriesSequence: """
        N = len(NewTrgSeg.ReferencedSeriesSequence[0]\
                         .ReferencedInstanceSequence)
        
        if i > N - 1:
            """ Increase the sequence by one. """
            #NewTrgSeg.ReferencedSeriesSequence[0]\
            #         .ReferencedInstanceSequence.append(NewTrgSeg.ReferencedSeriesSequence[0]\
            #                                                     .ReferencedInstanceSequence[-1])
                     
            LastItem = deepcopy(NewTrgSeg.ReferencedSeriesSequence[0]\
                                         .ReferencedInstanceSequence[-1])
            
            NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence.append(LastItem)
        
        """ Update the sequence. """
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPClassUID = TrgDicoms[s].SOPClassUID
        
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
    
    """ Check if there are more sequences in ReferencedInstanceSequence than 
    required. """
    N = len(NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    if N > len(UniqueTrgF2Sinds):
        for i in range(N - len(UniqueTrgF2Sinds)):
            """ Remove the last sequence. """
            NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence.pop()
    
    """ Modify SegmentSequence. """
    
    """ Loop through all SrcSegNums. """
    for s in range(len(SrcSegNums)):
        SegNum = SrcSegNums[s]
        
        """ Replace the s^th SegmentSequence in NewTrgSeg with the SegNum^th 
        SegmentSequence in SrcSeg."""
        NewTrgSeg.SegmentSequence[s] = deepcopy(SrcSeg.SegmentSequence[SegNum])
        
        NewTrgSeg.SegmentSequence[s].SegmentNumber = int(s+1)
        
        NewTrgSeg.SegmentSequence[s]\
                 .SegmentLabel = SrcSeg.SegmentSequence[SegNum].SegmentLabel
        
    """ Check if there are more sequences in SegmentSequence than required. """
    N = len(NewTrgSeg.SegmentSequence)
    
    if N > len(SrcSegNums):
        for i in range(N - len(SrcSegNums)):
            """ Remove the last sequence: """
            NewTrgSeg.SegmentSequence.pop()
    
    """ Modify PerFrameFunctionalGroupsSequence. """
    
    """ Get the ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence."""
    RefSOPsInRIS = get_RSOP_uids_in_RIS(NewTrgSeg)
    
    i = 0 # total frame counter (for all segments)
    
    for s in range(len(TrgF2SindsBySeg)):
        for f in range(len(TrgF2SindsBySeg[s])):  
            N = len(NewTrgSeg.PerFrameFunctionalGroupsSequence)
            
            if i > N - 1:
                """ Increase the sequence by one. """
                #NewTrgSeg.PerFrameFunctionalGroupsSequence\
                #         .append(NewTrgSeg.PerFrameFunctionalGroupsSequence[-1])
                         
                LastItem = deepcopy(NewTrgSeg.PerFrameFunctionalGroupsSequence[-1])
                
                NewTrgSeg.PerFrameFunctionalGroupsSequence\
                         .append(LastItem)
            
            """ The DICOM slice number: """
            d = TrgF2SindsBySeg[s][f]
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .DerivationImageSequence[0]\
                     .SourceImageSequence[0]\
                     .ReferencedSOPClassUID = TrgDicoms[d].SOPClassUID
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .DerivationImageSequence[0]\
                     .SourceImageSequence[0]\
                     .ReferencedSOPInstanceUID = TrgDicoms[d].SOPInstanceUID
            
            SOPuid = TrgDicoms[d].SOPInstanceUID
            
            ind = RefSOPsInRIS.index(SOPuid)
            
            """ Note: While s (segment number) and ind (the index of SOPuid
            within the ReferencedSOPInstanceUIDs in ReferencedInstanceSequence)
            are are both 0-indexed, DimensionIndexValues are 1-indexed. """
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .FrameContentSequence[0]\
                     .DimensionIndexValues = [s + 1, ind + 1]
            
            IPP = TrgDicoms[d].ImagePositionPatient
            
            # Ensure that the characters are limited to 16:
            IPP = reduce_list_of_str_floats_to_16(IPP)
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .PlanePositionSequence[0]\
                     .ImagePositionPatient = IPP
               
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .SegmentIdentificationSequence[0]\
                     .ReferencedSegmentNumber = s + 1
            
            i += 1 # increment the total frame count
    
    """ Check if there are more sequences in PerFrameFunctionalGroupsSequence 
    than required. """
    N = len(NewTrgSeg.PerFrameFunctionalGroupsSequence)
    
    if N > len(FlattenedTrgF2Sinds):
        for i in range(N - len(FlattenedTrgF2Sinds)):
            """ Remove the last sequence: """
            NewTrgSeg.PerFrameFunctionalGroupsSequence.pop()
    
    """ Modify PixelData. """
    
    """ Start by copying the first (or only) frame: """
    PixArr = deepcopy(TrgPixArrBySeg[0])
    
    if len(TrgPixArrBySeg) > 1:
        """ If there are more than one frame, vertically stack all pixel arrays 
        into a single pixel array. """
        for s in range(1, len(TrgPixArrBySeg)):
            PixArr = np.vstack((PixArr, TrgPixArrBySeg[s]))
    
    """ Note: The following doesn't work for binary arrays:
    NewTrgSeg.PixelData = PixArr.tobytes()
    """
    
    """ Ravel (to 1D array) PixArr and pack bits (see
    https://github.com/pydicom/pydicom/issues/1230). """
    packed = pack_bits(PixArr.ravel())
    
    """ Convert PixArr to bytes. """
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    if LogToConsole:
        print('-'*120)
        
    return NewTrgSeg