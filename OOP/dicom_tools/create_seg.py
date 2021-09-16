# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 13:47:10 2021

@author: ctorti
"""



def create_seg(srcDataset, trgDataset, newDataset, params):
    # TODO update docstrings
    """
    Create an SEG object for the target dataset.
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    newDataset : Propagator Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    newSeg : Pydicom Object
        The new ROI Collection for the target dataset.
    
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
    
    from importlib import reload
    import general_tools.general
    reload(general_tools)
    
    #from pydicom import dcmread
    from pydicom.uid import generate_uid
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    #import time
    import datetime
    from copy import deepcopy
    import numpy as np
    from general_tools.general import (
        flatten_list, get_unique_items, reduce_list_of_str_floats_to_16
        )
    #from dicom_tools.imports import import_dicoms
    #from dicom_tools.metadata import get_roicol_labels
    #from image_tools.attrs_info import get_im_attrs
    from seg_tools.metadata import get_RSOPuids_in_RIS
    from general_tools.console_printing import (
        print_indsByRoi, print_pixarrBySeg
    )
    
    srcSeg = srcDataset.roicol
    srcSegLabel = srcDataset.roiName
    allSrcSegLabels = srcDataset.roiNames
    
    """ 
    Since the source SEG has been copied/propagated to the target domain, it's
    the tags from the target DICOMs that are relevant. The image attributes of
    newDataset are the same as that of trgDataset.
    """
    trgDicoms = trgDataset.dicoms
    trgSeg = trgDataset.roicol
    #trgImSize = trgDataset.imSize
    trgImSpacings = trgDataset.imSpacings
    #trgSlcThick = trgDataset.imSlcThick
    #trgIPPs = trgDataset.imPositions
    #trgDirs = trgDataset.imDirections 
    
    newF2SindsBySeg = newDataset.f2sIndsBySeg
    newPixarrBySeg = newDataset.pixarrBySeg
    
    addToRoicolLab = params.cfgDict['addToRoicolLab']
    p2c = params.cfgDict['p2c']
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of create_seg():')
        print('\n\n', '-'*120)
        print_indsByRoi(newF2SindsBySeg)
        print_pixarrBySeg(newPixarrBySeg)
    
    """ Use trgSeg or srcSeg as a template for newSeg. """
    if trgSeg:
        newSeg = deepcopy(trgSeg)
    else:
        newSeg = deepcopy(srcSeg)
    
    numOfFramesBySeg = []
    for s in range(len(newPixarrBySeg)):
        numOfFramesBySeg.append(newPixarrBySeg[s].shape[0])
    
    #print(f'\n\n\nnumOfFramesBySeg = {numOfFramesBySeg}')
    
    numOfFrames = sum(numOfFramesBySeg)
    
    #NumOfSegs = len(newF2SindsBySeg)
    
    # The list of segment labels in srcSeg:
    #allSrcSegLabels = get_roicol_labels(srcSeg)
    
    if srcSegLabel:
        # The list of segment numbers for the segments of interest:
        #srcSegNums = get_roicol_labels(srcSeg, srcSegLabel)
        srcSegNums = srcDataset.roiNums
    else:
        srcSegNums = list(range(len(allSrcSegLabels)))
    
    uniqueTrgF2Sinds = get_unique_items(
        items=newF2SindsBySeg, ignoreZero=False, maintainOrder=True
        )
    
    #print(newF2SindsBySeg)
    flattenedTrgF2Sinds = flatten_list(newF2SindsBySeg)
    
    if p2c:
        print(f'   srcSegNums = {srcSegNums}')
        print(f'   numOfFramesBySeg = {numOfFramesBySeg}')
        print(f'   numOfFrames = {numOfFrames}')
    
    # Generate a new SOPInstanceUID.
    newSOPuid = generate_uid()
    newSeg.SOPInstanceUID = deepcopy(newSOPuid)
    newSeg.file_meta.MediaStorageSOPInstanceUID = deepcopy(newSOPuid)
    
    # Generate a new SeriesInstanceUID.
    newSeg.SeriesInstanceUID = generate_uid()
    
    #currentDate = time.strftime("%Y%m%d", time.gmtime())
    #currentTime = time.strftime("%H%M%S", time.gmtime())
    timeNow = datetime.datetime.now()
    currentDate = timeNow.strftime('%Y%m%d')
    currentTime = timeNow.strftime('%H%M%S.%f')
    
    if addToRoicolLab == '':
        addToRoicolLab = timeNow.strftime(' %Y%m%d %H%M%S')
    
    """ 
    If trgSeg != None, some tags will not need to be replaced. 
    If trgSeg = None, use the corresponding values in the first Target DICOM.
    """
    
    if trgSeg == None:
        newSeg.StudyDate = trgDicoms[0].StudyDate
        try:
            newSeg.SeriesDate = trgDicoms[0].SeriesDate
        except AttributeError:
            pass
        #newSeg.ContentDate = trgDicoms[0].ContentDate
        newSeg.ContentDate = currentDate
        newSeg.StudyTime = trgDicoms[0].StudyTime
        try:
            newSeg.SeriesTime = trgDicoms[0].SeriesTime
        except AttributeError:
            pass
        #newSeg.ContentTime = trgDicoms[0].ContentTime
        newSeg.ContentTime = currentTime
        newSeg.Manufacturer = trgDicoms[0].Manufacturer
        try:
            newSeg.SeriesDescription += addToRoicolLab
        except AttributeError:
            newSeg.SeriesDescription = addToRoicolLab
        newSeg.PatientName = trgDicoms[0].PatientName
        newSeg.PatientID = trgDicoms[0].PatientID
        newSeg.PatientBirthDate = trgDicoms[0].PatientBirthDate
        newSeg.PatientSex = trgDicoms[0].PatientSex
        try:
            newSeg.PatientAge = trgDicoms[0].PatientAge
        except AttributeError:
            pass
        newSeg.StudyInstanceUID = trgDicoms[0].StudyInstanceUID
        newSeg.StudyID = trgDicoms[0].StudyID
        newSeg.SeriesNumber = trgDicoms[0].SeriesNumber
        newSeg.InstanceNumber = trgDicoms[0].InstanceNumber
        newSeg.FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
        newSeg.PositionReferenceIndicator = trgDicoms[0].PositionReferenceIndicator
        newSeg.Rows = trgDicoms[0].Rows
        newSeg.Columns = trgDicoms[0].Columns
        
        IOP = trgDicoms[0].ImageOrientationPatient
        
        #print(f'type(IOP[0]) = {type(IOP[0])}\n')
        
        """
        ImageOrientationPatient should have VR Decimal String (DS) but it
        appears to be of type pydicom.valuerep.DSfloat.
        """
        
        # Convert IOP to str:
        IOP = [str(item) for item in IOP]
        
        #print(f'IOP = {IOP}\n')
        #print(f'type(IOP[0]) = {type(IOP[0])}\n')
        
        # Ensure that the character length is limited to 16:
        IOP = reduce_list_of_str_floats_to_16(IOP)
        
        # Convert to string with maximum:
        #IOP = [float(item) for item in IOP]
        
        
        newSeg.SharedFunctionalGroupsSequence[0]\
              .PlaneOrientationSequence[0]\
              .ImageOrientationPatient = IOP
           
        """ 
        Notes:
        1. SliceThickness in PixelMeasuresSequence appears to be the z-Spacing 
        rather than SliceThickness from the DICOM metadata. 
        2. The character limit of VR DS is 16 characters.
        """
        zSpacing = f"{trgImSpacings[2]}"
        
        #if len(zSpacing) > 16:
        #    zSpacing = zSpacing[:16]
        
        zSpacing = reduce_list_of_str_floats_to_16([zSpacing])[0]
        
        newSeg.SharedFunctionalGroupsSequence[0]\
              .PixelMeasuresSequence[0]\
              .SliceThickness = deepcopy(zSpacing)
        
        newSeg.SharedFunctionalGroupsSequence[0]\
              .PixelMeasuresSequence[0]\
              .SpacingBetweenSlices = deepcopy(zSpacing)
    
        newSeg.SharedFunctionalGroupsSequence[0]\
              .PixelMeasuresSequence[0]\
              .PixelSpacing = trgDicoms[0].PixelSpacing
    
    # Modify NumberOfFrames.
    newSeg.NumberOfFrames = f"{len(flattenedTrgF2Sinds)}"
    
    # Modify ReferencedSeriesSequence.
    newSeg.ReferencedSeriesSequence[0]\
          .SeriesInstanceUID = trgDicoms[0].SeriesInstanceUID
    
    # Modify ReferencedInstanceSequence in ReferencedSeriesSequence.
    for i in range(len(uniqueTrgF2Sinds)):
        """ The slice index s determines the SOPInstanceUID to reference. """
        s = uniqueTrgF2Sinds[i]
        
        # The number of sequences in ReferencedSeriesSequence:
        N = len(newSeg.ReferencedSeriesSequence[0]\
                      .ReferencedInstanceSequence)
        
        if i > N - 1:
            # Increase the sequence by one.
            #newSeg.ReferencedSeriesSequence[0]\
            #      .ReferencedInstanceSequence.append(
            #          newSeg.ReferencedSeriesSequence[0]\
            #                 .ReferencedInstanceSequence[-1]
            #                 )
                     
            lastItem = deepcopy(newSeg.ReferencedSeriesSequence[0]\
                                      .ReferencedInstanceSequence[-1])
            
            newSeg.ReferencedSeriesSequence[0]\
                  .ReferencedInstanceSequence.append(lastItem)
        
        # Update the sequence.
        newSeg.ReferencedSeriesSequence[0]\
              .ReferencedInstanceSequence[i]\
              .ReferencedSOPClassUID = trgDicoms[s].SOPClassUID
        
        newSeg.ReferencedSeriesSequence[0]\
              .ReferencedInstanceSequence[i]\
              .ReferencedSOPInstanceUID = trgDicoms[s].SOPInstanceUID
    
    # Check if there are more sequences in ReferencedInstanceSequence than 
    # required.
    N = len(newSeg.ReferencedSeriesSequence[0]\
                  .ReferencedInstanceSequence)
    
    if N > len(uniqueTrgF2Sinds):
        for i in range(N - len(uniqueTrgF2Sinds)):
            # Remove the last sequence.
            newSeg.ReferencedSeriesSequence[0]\
                  .ReferencedInstanceSequence.pop()
    
    """ Modify SegmentSequence. """
    
    # Loop through all srcSegNums.
    for s in range(len(srcSegNums)):
        SegNum = srcSegNums[s]
        
        """Replace the s^th SegmentSequence in newSeg with the SegNum^th 
        SegmentSequence in srcSeg."""
        newSeg.SegmentSequence[s] = deepcopy(srcSeg.SegmentSequence[SegNum])
        
        newSeg.SegmentSequence[s].SegmentNumber = int(s+1)
        
        newSeg.SegmentSequence[s]\
              .SegmentLabel = srcSeg.SegmentSequence[SegNum].SegmentLabel
        
    # Check if there are more sequences in SegmentSequence than required.
    N = len(newSeg.SegmentSequence)
    
    if N > len(srcSegNums):
        for i in range(N - len(srcSegNums)):
            # Remove the last sequence: 
            newSeg.SegmentSequence.pop()
    
    """ Modify PerFrameFunctionalGroupsSequence. """
    
    # Get the ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence.
    refSOPsInRIS = get_RSOPuids_in_RIS(newSeg)
    
    i = 0 # total frame counter (for all segments)
    
    for s in range(len(newF2SindsBySeg)):
        for f in range(len(newF2SindsBySeg[s])):  
            N = len(newSeg.PerFrameFunctionalGroupsSequence)
            
            if i > N - 1:
                # Increase the sequence by one.
                #newSeg.PerFrameFunctionalGroupsSequence\
                #      .append(newSeg.PerFrameFunctionalGroupsSequence[-1])
                         
                lastItem = deepcopy(newSeg.PerFrameFunctionalGroupsSequence[-1])
                
                newSeg.PerFrameFunctionalGroupsSequence\
                      .append(lastItem)
            
            # The DICOM slice number:
            d = newF2SindsBySeg[s][f]
            
            newSeg.PerFrameFunctionalGroupsSequence[i]\
                  .DerivationImageSequence[0]\
                  .SourceImageSequence[0]\
                  .ReferencedSOPClassUID = trgDicoms[d].SOPClassUID
            
            newSeg.PerFrameFunctionalGroupsSequence[i]\
                  .DerivationImageSequence[0]\
                  .SourceImageSequence[0]\
                  .ReferencedSOPInstanceUID = trgDicoms[d].SOPInstanceUID
            
            SOPuid = trgDicoms[d].SOPInstanceUID
            
            ind = refSOPsInRIS.index(SOPuid)
            
            """ Note: While s (segment number) and ind (the index of SOPuid
            within the ReferencedSOPInstanceUIDs in ReferencedInstanceSequence)
            are are both 0-indexed, DimensionIndexValues are 1-indexed. """
            newSeg.PerFrameFunctionalGroupsSequence[i]\
                  .FrameContentSequence[0]\
                  .DimensionIndexValues = [s + 1, ind + 1]
            
            IPP = trgDicoms[d].ImagePositionPatient
            
            # Convert IPP to str:
            IPP = [str(item) for item in IPP]
            
            # Ensure that the characters are limited to 16:
            IPP = reduce_list_of_str_floats_to_16(IPP)
            
            newSeg.PerFrameFunctionalGroupsSequence[i]\
                  .PlanePositionSequence[0]\
                  .ImagePositionPatient = IPP
               
            newSeg.PerFrameFunctionalGroupsSequence[i]\
                  .SegmentIdentificationSequence[0]\
                  .ReferencedSegmentNumber = s + 1
            
            i += 1 # increment the total frame count
    
    # Check if there are more sequences in PerFrameFunctionalGroupsSequence 
    # than required.
    N = len(newSeg.PerFrameFunctionalGroupsSequence)
    
    if N > len(flattenedTrgF2Sinds):
        for i in range(N - len(flattenedTrgF2Sinds)):
            # Remove the last sequence:
            newSeg.PerFrameFunctionalGroupsSequence.pop()
    
    """ Modify PixelData. """
    
    # Start by copying the first (or only) frame:
    pixarr = deepcopy(newPixarrBySeg[0])
    
    if len(newPixarrBySeg) > 1:
        """ If there are more than one frame, vertically stack all pixel arrays 
        into a single pixel array. """
        for s in range(1, len(newPixarrBySeg)):
            pixarr = np.vstack((pixarr, newPixarrBySeg[s]))
    
    """ Note: The following doesn't work for binary arrays:
    newSeg.PixelData = pixarr.tobytes()
    """
    
    """ Ravel (to 1D array) pixarr and pack bits (see
    https://github.com/pydicom/pydicom/issues/1230). """
    packed = pack_bits(pixarr.ravel())
    
    """ Convert pixarr to bytes. """
    newSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    if p2c:
        print('-'*120)
    
    return newSeg