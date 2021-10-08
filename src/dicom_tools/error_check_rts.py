# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:03:21 2021

@author: ctorti
"""


def error_check_rts(rts, trgDataset, params):
    """
    Check a RTS for errors in dependencies based on provided directory of
    the DICOMs that relate to the RTS.  
    
    Parameters
    ----------
    rts : Pydicom object
        RTS object to be error checked.
    trgDataset : DataImporter Object
        DataImporter Object for the target DICOM series.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    logList : list of strs
        A list of strings that describe any errors that are found.  If no 
        errors are found an empty list ([]) will be returned.
    Nerrors : int
        The number of errors found.
    
    Notes
    -----
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
    
    from copy import deepcopy
    
    """ 
    Since the source RTS has been copied/propagated to the target domain, it's
    the tags from the target DICOMs that are relevant. The image attributes of
    newDataset are the same as that of trgDataset.
    """
    dicom = trgDataset.dicoms[0]
    SOPuids = trgDataset.sopuids
    #spacings = trgDataset.imSpacings
    #ST = trgDataset.imSlcThick
    #IPPs = trgDataset.imPositions
    
    p2c = params.cfgDict['p2c']
    
    logList = []
    Nerrors = 0
    
    if p2c:
        print('\n\n', '-'*120)
        #PrintTitle('Error checking RTS..')
        print('Running of running error_check_rts():')
        print('\n\n', '-'*120)
    
    # Determine whether MediaStorageSOPInstanceUID matches SOPInstanceUID
    MSSOPuid = deepcopy(rts.file_meta.MediaStorageSOPInstanceUID)
    SOPuid = deepcopy(rts.SOPInstanceUID)
    
    if MSSOPuid == SOPuid:
        msg = f'INFO:  MediaStorageSOPInstanceUID {MSSOPuid} does not ' +\
            f'match SOPInstanceUID {SOPuid}.\n'
    else:
        msg = f'ERROR:  MediaStorageSOPInstanceUID {MSSOPuid} does not ' +\
            f'match SOPInstanceUID {SOPuid}.\n'   
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if rts.PatientName == dicom.PatientName:
        msg = f'INFO:  RTS PatientName {rts.PatientName} matches DICOM '\
              + f'PatientName {dicom.PatientName}.\n'
    else:
        msg = f'ERROR:  RTS PatientName {rts.PatientName} does not match '\
              + f'DICOM PatientName {dicom.PatientName}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
        
    if rts.PatientID == dicom.PatientID:
        msg = f'INFO:  RTS PatientID {rts.PatientID} matches ' +\
            f'DICOM PatientID {dicom.PatientID}.\n'
    else:
        msg = f'ERROR:  RTS PatientID {rts.PatientID} does not match ' +\
            f'DICOM PatientID {dicom.PatientID}.\n'    
        Nerrors +=  1 
    logList.append(msg)
    if p2c:
        print(msg)

    if rts.StudyInstanceUID == dicom.StudyInstanceUID:
        msg = f'INFO:  RTS StudyInstanceUID {rts.StudyInstanceUID} ' +\
            f'matches DICOM StudyInstanceUID {dicom.StudyInstanceUID}.\n'
    else:
        msg = f'ERROR:  RTS StudyInstanceUID {rts.StudyInstanceUID} does ' +\
            f'not match DICOM StudyInstanceUID {dicom.StudyInstanceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
            
    #if rts.StudyID == dicom.StudyID:
    #    msg = f'INFO:  The RTS StudyID {rts.StudyID} matches '\
    #          + f'the DICOM StudyID {dicom.StudyID}.\n'
    #else:
    #    msg = f'ERROR:  The RTS StudyID {rts.StudyID} does not match '\
    #          + f'the DICOM StudyID {dicom.StudyID}.\n'
    #
    #    Nerrors +=  1
    #    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
            
    if rts.FrameOfReferenceUID == dicom.FrameOfReferenceUID:
        msg = f'INFO:  RTS FrameOfReferenceUID {rts.FrameOfReferenceUID}' +\
            '  matches the DICOM FrameOfReferenceUID ' +\
            f'{dicom.FrameOfReferenceUID}.\n'
    else:
        msg = f'ERROR:  RTS FrameOfReferenceUID {rts.FrameOfReferenceUID}' +\
            '  does not match the DICOM FrameOfReferenceUID ' +\
            f'{dicom.FrameOfReferenceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    #RFORS = deepcopy(rts.ReferencedFrameOfReferenceSequence[0])
    RFORS = deepcopy(rts.ReferencedFrameOfReferenceSequence)
    
    #FORuidInRFORS = deepcopy(RFORS.FrameOfReferenceUID)
    FORuidInRFORS = deepcopy(RFORS[0].FrameOfReferenceUID)
    
    if FORuidInRFORS == dicom.FrameOfReferenceUID:
        msg = 'INFO:  FrameOfReferenceUID in ' +\
            f'ReferencedFrameOfReferenceSequence {FORuidInRFORS} ' +\
            'matches DICOM FrameOfReferenceUID ' +\
            f'{dicom.FrameOfReferenceUID}.\n'
    else:
        msg = 'ERROR:  FrameOfReferenceUID in ' +\
            f'ReferencedFrameOfReferenceSequence {FORuidInRFORS} does ' +\
            'not match DICOM FrameOfReferenceUID ' +\
            f'{dicom.FrameOfReferenceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Verify that the ReferencedSOPInstanceUID in 
    # ReferencedFrameOfReferenceSequence matches the RTS StudyInstanceUID.
    #refSOPuid = RFORS.RTReferencedStudySequence[0].ReferencedSOPInstanceUID
    refSOPuid = RFORS[0].RTReferencedStudySequence[0].ReferencedSOPInstanceUID
    
    if refSOPuid == rts.StudyInstanceUID:
        msg = 'INFO:  ReferencedSOPInstanceUID in ' +\
            f'ReferencedFrameOfReferenceSequence {refSOPuid} matches ' +\
            f'the RTS StudyInstanceUID {rts.StudyInstanceUID}.\n'
    else:
        msg = 'ERROR:  ReferencedSOPInstanceUID in ' +\
            f'ReferencedFrameOfReferenceSequence {refSOPuid} does not ' +\
            f'match the RTS StudyInstanceUID {rts.StudyInstanceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Verify that the SeriesInstanceUID in ReferencedFrameOfReferenceSequence
    # matches the DICOM SeriesInstanceUID.
    #rtsSeriesuid = deepcopy(RFORS.RTReferencedStudySequence[0]\
    #                             .RTReferencedSeriesSequence[0]\
    #                             .SeriesInstanceUID)
    rtsSeriesuid = deepcopy(RFORS[0].RTReferencedStudySequence[0]\
                                    .RTReferencedSeriesSequence[0]\
                                    .SeriesInstanceUID)
    
    if rtsSeriesuid == dicom.SeriesInstanceUID:
        msg = 'INFO:  SeriesInstanceUID in ' +\
            f'ReferencedFrameOfReferenceSequence {rtsSeriesuid} matches ' +\
            f' DICOM SeriesInstanceUID {dicom.SeriesInstanceUID}.\n'
    else:
        msg = 'ERROR:  SeriesInstanceUID in ' +\
            f'ReferencedFrameOfReferenceSequence {rtsSeriesuid} does not ' +\
            f'match DICOM SeriesInstanceUID {dicom.SeriesInstanceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Verify that the ReferencedSOPInstanceUIDs in 
    # ReferencedFrameOfReferenceSequence match the DICOM SOPInstanceUIDs.
    #CIS = deepcopy(RFORS.RTReferencedStudySequence[0]\
    #                    .RTReferencedSeriesSequence[0]\
    #                    .ContourImageSequence)
    CIS = deepcopy(RFORS[0].RTReferencedStudySequence[0]\
                           .RTReferencedSeriesSequence[0]\
                           .ContourImageSequence)
    
    refSOPuidsInRFORS = [
        CIS[i].ReferencedSOPInstanceUID for i in range(len(CIS))
        ]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    isMatch = [refSOPuid in SOPuids for refSOPuid in refSOPuidsInRFORS]
    
    #numOfMatches = isMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    inds = [i for i, x in enumerate(isMatch) if x==False]
    
    if inds:
        nonMatchingUids = [refSOPuidsInRFORS[i] for i in inds]
        
        for i in range(len(nonMatchingUids)):
            uid = nonMatchingUids[i]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in ' +\
                'ContourImageSequence {i+1} does not match any DICOM ' +\
                'SOPInstanceUID.\n'
            Nerrors +=  1
    else:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs in ' +\
            'ContourImageSequence match the DICOM SOPInstanceUIDs.\n'
    logList.append(msg)
    if p2c:
        print(msg)
    
    ## Verify that the StructureSetROISequence has only 1 sequence.
    #N = len(rts.StructureSetROISequence)
    #
    #if N == 1:
    #    msg = f'INFO:  There are {N} sequences in StructureSetROISequence.'\
    #          + ' There should be 1.\n'
    #else:
    #    msg = f'ERROR:  There are {N} sequences in StructureSetROISequence.'\
    #          + ' There should be 1.\n'
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
    
    ## Verify that the ROINumber in StructureSetROISequence is 1.
    #N = int(deepcopy(rts.StructureSetROISequence[0].ROINumber))
    #
    #if N == 1:
    #    msg = f'INFO:  ROINumber in StructureSetROISequence is {N}.  It '\
    #          + 'should be 1.\n'
    #else:
    #    msg = f'ERROR:  ROINumber in StructureSetROISequence is {N}.  It '\
    #          + 'should be 1.\n'
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
    
    # Verify that the ReferencedFrameOfReferenceUID in StructureSetROISequence
    # matches the DICOM FrameOfReferenceUID.
    refFORuidInSSRS = deepcopy(rts.StructureSetROISequence[0]\
                                  .ReferencedFrameOfReferenceUID)
    
    if refFORuidInSSRS == dicom.FrameOfReferenceUID:
        msg = 'INFO:  ReferencedFrameOfReferenceUID in ' +\
            f'StructureSetROISequence {refFORuidInSSRS} matches DICOM ' +\
            f'FrameOfReferenceUID {dicom.FrameOfReferenceUID}.\n'
    else:
        msg = 'ERROR:  ReferencedFrameOfReferenceUID in ' +\
            f'StructureSetROISequence {refFORuidInSSRS} does not match ' +\
            f'DICOM FrameOfReferenceUID {dicom.FrameOfReferenceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    ## Verify that the ROIContourSequence has only 1 sequence.
    #N = len(rts.ROIContourSequence)
    #
    #if N == 1:
    #    msg = f'INFO:  There are {N} sequences in ROIContourSequence. There '\
    #          + 'should be 1.\n'
    #else:
    #    msg = f'ERROR:  There are {N} sequences in ROIContourSequence. There '\
    #          + 'should be 1.\n'
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
    
    # Verify that the ReferencedSOPInstanceUIDs in ROIContourSequence
    # match the ReferencedSOPInstanceUIDs in ContourImageSequence.
    CS = deepcopy(rts.ROIContourSequence[0].ContourSequence)
    
    refSOPuidsInRCS = [
        CS[i].ContourImageSequence[0].ReferencedSOPInstanceUID for i in range(len(CS))
        ]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    isMatch = [uid in refSOPuidsInRFORS for uid in refSOPuidsInRCS]
    
    #numOfMatches = isMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    inds = [i for i, x in enumerate(isMatch) if x==False]
    
    if inds:
        nonMatchingUids = [refSOPuidsInRCS[i] for i in inds]
        
        for i in range(len(nonMatchingUids)):
            uid = nonMatchingUids[i]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in ' +\
                f'ContourSequence {i+1} does not match any ' +\
                'ReferencedSOPInstanceUID in ContourImageSequence.\n'
            Nerrors +=  1
    else:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs in ContourSequence ' +\
            'match the ReferencedSOPInstanceUIDs in ContourImageSequence.\n'
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Verify that the ContourNumbers in each ContourSequence range from 1 to
    # len(CS).
    expectedNums = list(range(1, len(CS)+1))
    
    contourNums = [int(CS[i].ContourNumber) for i in range(len(CS))]
    
    #print(f'\n\n\ncontourNums = {contourNums}')
    
    # Determine whether the contour numbers match the expected contour numbers:
    isMatch = [contourNums[i] == expectedNums[i] for i in range(len(CS))]
    #isMatch = [contourNums[i] == i+1 for i in range(len(CS))]
    
    # Find the indices of any non-matching contour numbers:
    inds = [i for i, x in enumerate(isMatch) if x==False]
    
    if inds:        
        for i in range(len(inds)):
            cn = contourNums[inds[i]]
            en = expectedNums[inds[i]]
            
            msg = f'ERROR:  ContourNumber {i+1} is {cn} but is expected ' +\
                f'to be {en}.\n'
            logList.append(msg)      
            Nerrors +=  1
    else:
        msg = 'INFO:  The ContourNumbers are as expected.\n'
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Verify that the number of elements in ContourData is 3x 
    # NumberOfContourPoints for each ContourSequence.
    NCP = [int(CS[i].NumberOfContourPoints) for i in range(len(CS))]
    
    NCD = [len(CS[i].ContourData) for i in range(len(CS))]
    
    isMatch = [NCD[i] == 3*NCP[i] for i in range(len(CS))]
    
    # Find the indices of any non-zero modulos:
    inds = [i for i, x in enumerate(isMatch) if x==False]
    
    if inds:        
        for i in range(len(inds)):
            ncp = NCP[inds[i]]
            
            ncd = NCD[inds[i]]
            
            msg = f'ERROR:  The number of elements in ContourData {ncd} is ' +\
                f'not 3x NumberOfContourPoints {ncp} for ContourSequence ' +\
                f'{i+1}.\n'
            Nerrors +=  1
    else:
        msg = 'INFO:  The number of elements in ContourData are 3x ' +\
            'NumberOfContourPoints for ContourSequence.\n'
    logList.append(msg)
    if p2c:
        print(msg)
    
    ## Verify that the ReferencedROINumber in ROIContourSequence is 1.
    #N = int(deepcopy(rts.ROIContourSequence[0].ReferencedROINumber))
    #
    #if N == 1:
    #    msg = f'INFO:  ReferencedROINumber in ROIContourSequence is {N}. '\
    #          + 'It should be 1.\n'
    #else:
    #    msg = f'ERROR:  ReferencedROINumber in ROIContourSequence is {N}. '\
    #          + 'It should be 1.\n'
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
    
    ## Verify that the ROIObservationsSequence has only 1 sequence.
    #N = len(rts.RTROIObservationsSequence)
    #
    #if N == 1:
    #    msg = f'INFO:  There are {N} sequences in RTROIObservationsSequence.'\
    #          + ' There should be 1.\n'
    #else:
    #    msg = f'ERROR:  There are {N} sequences in RTROIObservationsSequence.'\
    #          + ' There should be 1.\n'  
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
    
    ## Verify that the ObservationNumber is 1.
    #N = int(deepcopy(rts.RTROIObservationsSequence[0].ObservationNumber))
    #
    #if N == 1:
    #    msg = 'INFO:  ObservationNumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #else:
    #    msg = 'ERROR:  ObservationNumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
            
    ## Verify that the ReferencedROINumber is 1.
    #N = int(deepcopy(rts.RTROIObservationsSequence[0].ReferencedROINumber))
    #
    #if N == 1:
    #    msg = 'INFO:  ReferencedROINumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #else:
    #    msg = 'ERROR:  ReferencedROINumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #    Nerrors +=  1
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #logList.append(msg)
    #if p2c:
    #    print(msg)
    
    #print(f'There were {Nerrors} errors found in the RTS.\n')
    
    if p2c:
        print('-'*120)
    
    return logList, Nerrors