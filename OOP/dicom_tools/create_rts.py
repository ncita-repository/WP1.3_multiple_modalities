# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 16:59:11 2021

@author: ctorti
"""



def create_rts(srcDataset, trgDataset, newDataset, params):
    """
    srcRtsFpath, TrgRtsFpath, TrgCntDataByCntByRoi, TrgPtsByCntByRoi, 
    TrgC2SindsByRoi, TrgDicomDir, #SrcRoiName, 
    AddTxtToSSLabel='', p2c=False
    
    Create an RTS object for the target dataset. 
    
    Parameters
    ----------
    srcDataset : DataImporter Object
        DataImporter Object for the source DICOM series.
    trgDataset : DataImporter Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    newRts : Pydicom object
        New RTS object for the Target dataset.
    
    Notes
    -----
    Performance of adding/removing sequences using append()/pop():
        
    * Took 329.9 ms to add 10000 sequences to ContourImageSequence 
    (33.0 us/sequence)
    * Took 165.4 ms to remove 10000 sequences to ContourImageSequence 
    (16.5 us/sequence)
    
    * Took 177.7 ms to add 10000 sequences to ContourSequence 
    (17.8 us/sequence)
    * Took 81.7 ms to remove 10000 sequences to ContourSequence 
    (8.2 us/sequence)
    
    19/02:
        After avoiding the use of deepcopy to increase sequences I've had to
        re-instate them since the tag values were being duplicated despite
        calls to modify them.  For Test RR3, creating the new RTS with 45 
        contours on 45 slices previously took 2 s to execute without using 
        deepcopy.  Using deepcopy it takes 11 s.
        
        This has also lengthened the time it takes to error check the RTS. The
        same RTS object formerly took 9 s to error check.  Now it takes 31 s.
        It's not clear why error checking should take longer.
        
        Test RR4 previously took 1 s to create the RTS and 5 s to error check
        it.  Now it takes 115 s and 343 s.
    """
    
    from importlib import reload
    import general_tools.general
    reload(general_tools)
    
    from pydicom.uid import generate_uid
    #import time
    import datetime
    from copy import deepcopy
    from general_tools.general import (
        get_unique_items#, flatten_list, reduce_list_of_str_floats_to_16
        )
    from general_tools.console_printing import (
        print_indsByRoi, print_ptsByCntByRoi
    )
    
    srcRts = srcDataset.roicol
    #srcSegLabel = srcDataset.roiName
    #allSrcSegLabels = srcDataset.roiNames
    
    """ 
    Since the source RTS has been copied/propagated to the target domain, it's
    the tags from the target DICOMs that are relevant. The image attributes of
    newDataset are the same as that of trgDataset.
    """
    trgDicoms = trgDataset.dicoms
    trgRts = trgDataset.roicol
    #trgImSize = trgDataset.imSize
    #trgImSpacings = trgDataset.imSpacings
    #trgSlcThick = trgDataset.imSlcThick
    #trgIPPs = trgDataset.imPositions
    #trgDirs = trgDataset.imDirections 
    
    newC2SindsByRoi = newDataset.c2sIndsByRoi
    newPtsByCntByRoi = newDataset.ptsByCntByRoi
    newCntdataByCntByRoi = newDataset.cntdataByCntByRoi
    
    addToRoicolLab = params.cfgDict['addToRoicolLab']
    p2c = params.cfgDict['p2c']
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of create_rts():')
        print('\n\n', '-'*120)
        print_indsByRoi(newC2SindsByRoi)
        print_ptsByCntByRoi(newPtsByCntByRoi)
        print_ptsByCntByRoi(newCntdataByCntByRoi)
    
    """ Use trgRts or srcRts as a template for newRts. """
    if trgRts:
        newRts = deepcopy(trgRts)
    else:
        newRts = deepcopy(srcRts)
    
    uniqueTrgC2Sinds = get_unique_items(
        items=newC2SindsByRoi, ignoreZero=False, maintainOrder=True
        )
    
    # Generate a new SOPInstanceUID.
    newSOPuid = generate_uid()
    newRts.SOPInstanceUID = deepcopy(newSOPuid)
    newRts.file_meta.MediaStorageSOPInstanceUID = deepcopy(newSOPuid)
    
    # Generate a new SeriesInstanceUID.
    newRts.SeriesInstanceUID = generate_uid()
    
    #currentDate = time.strftime("%Y%m%d", time.gmtime())
    #currentTime = time.strftime("%H%M%S", time.gmtime())
    
    timeNow = datetime.datetime.now()
    currentDate = timeNow.strftime('%Y%m%d')
    currentTime = timeNow.strftime('%H%M%S.%f')
    
    if addToRoicolLab == '':
        addToRoicolLab = timeNow.strftime(' %Y%m%d %H%M%S')
    
    """ 
    If trgRts != None, some tags will not need to be replaced.
    If trgRts = None, use the corresponding values in the first Target DICOM. 
    """
    
    if trgRts == None:
        newRts.StudyDate = trgDicoms[0].StudyDate
        newRts.StudyTime = trgDicoms[0].StudyTime
        newRts.Manufacturer = trgDicoms[0].Manufacturer
        newRts.ManufacturerModelName = trgDicoms[0].ManufacturerModelName
        #newRts.SeriesDate = Dicoms[0].SeriesDate
        newRts.PatientName = trgDicoms[0].PatientName
        newRts.PatientID = trgDicoms[0].PatientID
        newRts.PatientBirthDate = trgDicoms[0].PatientBirthDate
        newRts.PatientSex = trgDicoms[0].PatientSex
        try:
            newRts.PatientAge = trgDicoms[0].PatientAge
        except AttributeError:
            pass
        newRts.StudyInstanceUID = trgDicoms[0].StudyInstanceUID
        newRts.StudyID = trgDicoms[0].StudyID
        newRts.SeriesNumber = trgDicoms[0].SeriesNumber
        newRts.FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
        newRts.PositionReferenceIndicator = trgDicoms[0].PositionReferenceIndicator
        
        newRts.ReferencedFrameOfReferenceSequence[0]\
              .FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
          
        newRts.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .ReferencedSOPInstanceUID = trgDicoms[0].StudyInstanceUID
        
        newRts.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .SeriesInstanceUID = trgDicoms[0].SeriesInstanceUID
       
    try:
        newRts.StructureSetLabel += addToRoicolLab
    except AttributeError:
        newRts.StructureSetLabel = addToRoicolLab
    # Modify StructureSetDate and StructureSetTime to the present:
    newRts.StructureSetDate = currentDate
    newRts.StructureSetTime = currentTime
    
    # Start timing:
    #times = []
    #times.append(time.time())
    
    """ Modify ContourImageSequence. """
    for i in range(len(uniqueTrgC2Sinds)):
        """ The slice index s determines the SOPInstanceUID to reference. """
        s = uniqueTrgC2Sinds[i]
        
        N = len(newRts.ReferencedFrameOfReferenceSequence[0]\
                      .RTReferencedStudySequence[0]\
                      .RTReferencedSeriesSequence[0]\
                      .ContourImageSequence)
        
        if i > N - 1:
            # Increase the sequence by one.
            #newRts.ReferencedFrameOfReferenceSequence[0]\
            #      .RTReferencedStudySequence[0]\
            #      .RTReferencedSeriesSequence[0]\
            #      .ContourImageSequence.append(
            #          newRts.ReferencedFrameOfReferenceSequence[0]\
            #                .RTReferencedStudySequence[0]\
            #                .RTReferencedSeriesSequence[0]\
            #                .ContourImageSequence[-1]
            #                )
                     
            lastItem = deepcopy(newRts.ReferencedFrameOfReferenceSequence[0]\
                                      .RTReferencedStudySequence[0]\
                                      .RTReferencedSeriesSequence[0]\
                                      .ContourImageSequence[-1])
            
            newRts.ReferencedFrameOfReferenceSequence[0]\
                  .RTReferencedStudySequence[0]\
                  .RTReferencedSeriesSequence[0]\
                  .ContourImageSequence.append(lastItem)
        
        # Update the sequence:
        newRts.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .ContourImageSequence[i]\
              .ReferencedSOPClassUID = trgDicoms[s].SOPClassUID
           
        newRts.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .ContourImageSequence[i]\
              .ReferencedSOPInstanceUID = trgDicoms[s].SOPInstanceUID
    
    # Check if there are more sequences in ContourImageSequence than required:
    N = len(newRts.ReferencedFrameOfReferenceSequence[0]\
                  .RTReferencedStudySequence[0]\
                  .RTReferencedSeriesSequence[0]\
                  .ContourImageSequence)
    
    if N > len(uniqueTrgC2Sinds):
        for i in range(N - len(uniqueTrgC2Sinds)):
            # Remove the last sequence:
            newRts.ReferencedFrameOfReferenceSequence[0]\
                  .RTReferencedStudySequence[0]\
                  .RTReferencedSeriesSequence[0]\
                  .ContourImageSequence.pop()
    
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#p2c:
    #    print(f'\nTook {Dtime} s to modify 'ContourImageSequence.')
    
    
    """ Modify StructureSetROISequence. """
    for r in range(len(newC2SindsByRoi)):
        # The number of contours in this ROI:
        #C = len(newC2SindsByRoi[r])
        N = len(newRts.StructureSetROISequence)
        
        if r > N - 1:
            # Increase the sequence by one.
            lastItem = deepcopy(newRts.StructureSetROISequence[-1])
            
            newRts.StructureSetROISequence.append(lastItem)
        
        # Update the sequence:
        newRts.StructureSetROISequence[r].ROINumber = f"{r+1}"
        
        newRts.StructureSetROISequence[r]\
              .ReferencedFrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
           
        newRts.StructureSetROISequence[r]\
              .ROIName = newRts.StructureSetROISequence[r]\
                               .ROIName + addToRoicolLab
    
    
    # Check if there are more sequences in StructureSetROISequence than 
    # required:
    N = len(newRts.StructureSetROISequence)
    
    if N > len(newC2SindsByRoi):
        for i in range(N - len(newC2SindsByRoi)):
            # Remove the last sequence:
            newRts.StructureSetROISequence.pop()
               
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#p2c:
    #    print(f'\nTook {Dtime} s to modify 'StructureSetROISequence.')
        
        
    """ 
    Modify ROIContourSequence. 
    """
    #print(f'\nnewC2SindsByRoi = {newC2SindsByRoi}')
    for r in range(len(newC2SindsByRoi)):
        
        n = 0 # total contour counter for this ROI
        
        # The number of ROIContourSequences:
        N = len(newRts.ROIContourSequence)
        
        if r > N - 1:
            # Increase the sequence by one.
            lastItem = deepcopy(newRts.ROIContourSequence[-1])
            
            newRts.ROIContourSequence.append(lastItem)
        
        for c in range(len(newC2SindsByRoi[r])):
            #print(f'\n\n\nn = {n}')
            
            # The DICOM slice number:
            s = newC2SindsByRoi[r][c]
            
            # The number of ContourSequences:
            N = len(newRts.ROIContourSequence[r].ContourSequence)
        
            if n > N - 1:
                # Increase the sequence by one.
                #newRts.ROIContourSequence[r]\
                #      .ContourSequence.append(newRts.ROIContourSequence[r]\
                #                                    .ContourSequence[-1])
                         
                lastItem = deepcopy(newRts.ROIContourSequence[r]\
                                          .ContourSequence[-1])
                
                newRts.ROIContourSequence[r]\
                      .ContourSequence.append(lastItem)
                   
            #print(f'\nnewRts.ROIContourSequence[{r}].ContourSequence[{i}].ContourImageSequence[0].ReferencedSOPInstanceUID =',
            #      f'{newRts.ROIContourSequence[r].ContourSequence[i].ContourImageSequence[0].ReferencedSOPInstanceUID}')
            #print(f'\nDicoms[{s}].SOPInstanceUID = {Dicoms[s].SOPInstanceUID}')
            
            # Update the sequence:
            newRts.ROIContourSequence[r]\
                  .ContourSequence[n]\
                  .ContourImageSequence[0]\
                  .ReferencedSOPClassUID = trgDicoms[s].SOPClassUID
               
            newRts.ROIContourSequence[r]\
                  .ContourSequence[n]\
                  .ContourImageSequence[0]\
                  .ReferencedSOPInstanceUID = trgDicoms[s].SOPInstanceUID
            
            newRts.ROIContourSequence[r]\
                  .ContourSequence[n]\
                  .NumberOfContourPoints = f"{len(newPtsByCntByRoi[r][c])}"
             
            #print(f'\nn = {n}')
            
            newRts.ROIContourSequence[r]\
                  .ContourSequence[n]\
                  .ContourNumber = f"{n+1}"
                     
            #print(f'newRts.ROIContourSequence[{r}].ContourSequence[{n}].ContourNumber =',
            #      f'{newRts.ROIContourSequence[r].ContourSequence[n].ContourNumber}')
               
            newRts.ROIContourSequence[r]\
                  .ContourSequence[n]\
                  .ContourData = newCntdataByCntByRoi[r][c]
            
            newRts.ROIContourSequence[r]\
                  .ReferencedROINumber = f"{r+1}"
            
            n += 1 # increment the total contour counter
        
        # Check if there are more sequences in ContourSequence than required:
        N = len(newRts.ROIContourSequence[r].ContourSequence)
        
        if N > len(newC2SindsByRoi[r]):
            #print(f'There are {N} sequences in ROIContourSequence[{r}].ContourSequence',
            #      f'but only {len(newC2SindsByRoi[r])} contours in this ROI.')
            for i in range(N - len(newC2SindsByRoi[r])):
                # Remove the last sequence:
                newRts.ROIContourSequence[r].ContourSequence.pop()
    
    #print('\nCheck ContourNumbers:')
    #for r in range(len(newRts.ROIContourSequence)):
    #    for n in range(len(newRts.ROIContourSequence[r].ContourSequence)):
    #        print(f'\nn = {n}')
    #                 
    #        print(f'newRts.ROIContourSequence[{r}].ContourSequence[{n}].ContourNumber =',
    #              f'{newRts.ROIContourSequence[r].ContourSequence[n].ContourNumber}')
    
    # Check if there are more sequences in ROIContourSequence than required:
    N = len(newRts.ROIContourSequence)
    
    if N > len(newC2SindsByRoi):
        for i in range(N - len(newC2SindsByRoi)):
            # Remove the last sequence:
            newRts.ROIContourSequence.pop()
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#p2c:
    #    print(f'\nTook {Dtime} s to modify 'ROIContourSequence.')
    
    
    """ 
    Modify RTROIObservationsSequence. 
    """
    for r in range(len(newC2SindsByRoi)):
        #""" The number of contours in this ROI: """
        #C = len(newC2SindsByRoi[r])
        N = len(newRts.RTROIObservationsSequence)
        
        if r > N - 1:
            # Increase the sequence by one.
            lastItem = deepcopy(newRts.RTROIObservationsSequence[-1])
            
            newRts.RTROIObservationsSequence.append(lastItem)
    
        # Update the sequence:
        newRts.RTROIObservationsSequence[r].ObservationNumber = f"{r+1}"
        
        newRts.RTROIObservationsSequence[r].ReferencedROINumber = f"{r+1}"
    
    # Check if there are more sequences in RTROIObservationsSequence than 
    # required:
    N = len(newRts.RTROIObservationsSequence)
    
    if N > len(newC2SindsByRoi):
        for i in range(N - len(newC2SindsByRoi)):
            # Remove the last sequence:
            newRts.RTROIObservationsSequence.pop()
            
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#p2c:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'RTROIObservationsSequence.')
    
    if p2c:
        print('-'*120)
    
    return newRts