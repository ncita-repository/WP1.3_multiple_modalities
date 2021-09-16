# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:02:33 2021

@author: ctorti
"""


def error_check_seg(seg, trgDataset, params):
    # TODO update docstrings
    """
    Check a SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the SEG.  
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object to be error checked.
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
    """
    
    from copy import deepcopy
    from general_tools.general import are_items_equal_to_within_eps
    
    """ 
    Since the source SEG has been copied/propagated to the target domain, it's
    the tags from the target DICOMs that are relevant. The image attributes of
    newDataset are the same as that of trgDataset.
    """
    dicom = trgDataset.dicoms[0]
    SOPuids = trgDataset.sopuids
    spacings = trgDataset.imSpacings
    #ST = trgDataset.imSlcThick
    IPPs = trgDataset.imPositions
    
    p2c = params.cfgDict['p2c']
    
    # The maximum allowed error between two lists/items:
    epsilon = 1e-05
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of error_check_seg():')
        print('\n\n', '-'*120)
    
    pixarr = seg.pixel_array
    
    # The number of frames, rows and columns in the SEG's pixel array:
    if len(pixarr.shape) > 2:
        numOfFrames, numOfRows, numOfCols = pixarr.shape
    else:
        numOfRows, numOfCols = pixarr.shape
        numOfFrames = 1
    
    #NumOfDicoms = len(SOPuids)
    
    if p2c:
        print(f'numOfFrames = {numOfFrames}')
    
    logList = []
    Nerrors = 0
    
    # Determine whether MediaStorageSOPInstanceUID matches SOPInstanceUID
    MSSOPuid = seg.file_meta.MediaStorageSOPInstanceUID
    SOPuid = seg.SOPInstanceUID
    
    if MSSOPuid == SOPuid:
        msg = f'INFO:  MediaStorageSOPInstanceUID {MSSOPuid} matches ' +\
            f'SOPInstanceUID {SOPuid}.\n'
    else:
        msg = f'ERROR:  MediaStorageSOPInstanceUID {MSSOPuid} does not ' +\
            f'match SOPInstanceUID {SOPuid}.\n'
        Nerrors +=  1 
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Determine whether any of the ReferencedSOPInstanceUIDs do not match the
    # SOPInstanceUIDs        
    RIS = seg.ReferencedSeriesSequence[0]\
             .ReferencedInstanceSequence
    
    refSOPuidsInRIS = [RIS[i].ReferencedSOPInstanceUID for i in range(len(RIS))]

    isMatch = [RefSOPuid in SOPuids for RefSOPuid in refSOPuidsInRIS]
    
    #numOfMatches = isMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    inds = [i for i, x in enumerate(isMatch) if x==False]
    
    if inds:
        for i in range(len(inds)):
            uid = refSOPuidsInRIS[inds[i]]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in ' +\
                f'ReferencedInstanceSequence[{inds[i]}] does not match any' +\
                'SOPInstanceUID.\n'
            logList.append(msg)
            Nerrors +=  1
            if p2c:
                print(msg)
    else:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs in ' +\
            'ReferencedInstanceSequence matches the SOPInstanceUIDs.\n'
        logList.append(msg)
        if p2c:
            print(msg)
                           
    if seg.ReferencedSeriesSequence[0].SeriesInstanceUID == dicom.SeriesInstanceUID:
        msg = 'INFO:  SEG SeriesInstanceUID in ReferencedSeriesSequence, ' +\
            f'{seg.ReferencedSeriesSequence[0].SeriesInstanceUID}, matches ' +\
            f'DICOM SeriesInstanceUID, {dicom.SeriesInstanceUID}.\n'
    else:
        msg = 'ERROR:  SEG SeriesInstanceUID in ReferencedSeriesSequence, ' +\
            f'{seg.ReferencedSeriesSequence[0].SeriesInstanceUID}, does not' +\
            f' match DICOM SeriesInstanceUID {dicom.SeriesInstanceUID}.\n'
        Nerrors +=  1 
    logList.append(msg)
    if p2c:
        print(msg)
    
    if seg.PatientName == dicom.PatientName:
        msg = f'INFO:  SEG PatientName, {seg.PatientName}, matches DICOM ' +\
            f'PatientName, {dicom.PatientName}.\n'
    else:
        msg = f'ERROR:  SEG PatientName, {seg.PatientName}, does not match '+\
            f'DICOM PatientName, {dicom.PatientName}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if seg.PatientID == dicom.PatientID:
        msg = f'INFO:  SEG PatientID, {seg.PatientID}, matches the DICOM ' +\
            f'PatientID, {dicom.PatientID}.\n'
    else:
        msg = f'ERROR:  SEG PatientID, {seg.PatientID}, does not match the' +\
            f' DICOM PatientID, {dicom.PatientID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if seg.StudyInstanceUID == dicom.StudyInstanceUID:
        msg = f'INFO:  SEG StudyInstanceUID, {seg.StudyInstanceUID}, ' +\
            f'matches the DICOM StudyInstanceUID, {dicom.StudyInstanceUID}.\n'
    else:
        msg = f'ERROR:  SEG StudyInstanceUID, {seg.StudyInstanceUID}, does '+\
            'not match the DICOM StudyInstanceUID, {dicom.StudyInstanceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
            
    #if seg.StudyID == dicom.StudyID:
    #    msg = f'INFO:  The SEG StudyID {seg.StudyID} matches the'\
    #          + f' DICOM StudyID {dicom.StudyID}.'
    #else:
    #    msg = f'ERROR:  The SEG StudyID {seg.StudyID} does not match the'\
    #          + f' DICOM StudyID {dicom.StudyID}.'
    #
    #    Nerrors +=  1
    #    
    #logList.append(msg)
    #if p2c:
    #    print(msg)
            
    if seg.FrameOfReferenceUID == dicom.FrameOfReferenceUID:
        msg = f'INFO:  SEG FrameOfReferenceUID, {seg.FrameOfReferenceUID}, ' +\
            f'matches DICOM FrameOfReferenceUID, {dicom.FrameOfReferenceUID}.\n'
    else:
        msg = f'ERROR:  SEG FrameOfReferenceUID, {seg.FrameOfReferenceUID}, '+\
            'does not match DICOM FrameOfReferenceUID, '\
            f'{dicom.FrameOfReferenceUID}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if int(seg.NumberOfFrames) == numOfFrames:
        msg = f'INFO:  NumberOfFrames, {seg.NumberOfFrames}, matches the ' +\
            f'number of frames in PixelData, {numOfFrames}.\n'
    else:
        msg = f'ERROR:  NumberOfFrames, {seg.NumberOfFrames}, does not ' +\
            f'match the number of frames in PixelData, {numOfFrames}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if seg.Rows == dicom.Rows:
        msg = f'INFO:  SEG Rows, {seg.Rows}, matches DICOM Rows, ' +\
            f'{dicom.Rows}.\n'
    else:
        msg = f'ERROR:  SEG Rows, {seg.Rows}, does not match DICOM Rows, ' +\
            f'{dicom.Rows}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if seg.Columns == dicom.Columns:
        msg = f'INFO:  SEG Columns, {seg.Columns}, matches DICOM ' +\
            f'Columns, {dicom.Columns}.\n'
    else:
        msg = f'ERROR:  SEG Columns, {seg.Columns}, does not match DICOM ' +\
            f'Columns, {dicom.Columns}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    """ This is already checked below
    if seg.Rows == numOfRows:
        msg = f'INFO:  SEG Rows, {seg.Rows}, matches the number of rows in '\
              + f'PixelData, {numOfRows}.\n'
    else:
        msg = f'ERROR:  SEG Rows, {seg.Rows}, does not match the number of '\
              + f'rows in PixelData, {numOfRows}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    if seg.Columns == numOfCols:
        msg = f'INFO:  SEG Columns, {seg.Columns}, matches the number of '\
              + f'columns in PixelData, {numOfCols}.\n'
    else:
        msg = f'ERROR:  SEG Columns, {seg.Columns}, does not match the '\
              + f'number of columns in PixelData, {numOfCols}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    """
    
    # Keep track of the number of sequences in SegmentSequence.
    numSS = len(seg.SegmentSequence)
    
    # Verify that the SegmentNumber in each SegmentSequence increments from 1.
    SS = deepcopy(seg.SegmentSequence)
    
    for i in range(numSS):
        N = int(SS[i].SegmentNumber)
        
        #print(f'\n\nN = {N}, type(N) = {type(N)}')
        
        if N == i + 1:
            msg = f'INFO:  SegmentNumber in SegmentSequence[{i}] is {N}.\n'
            logList.append(msg)
        else:
            msg = f'ERROR:  SegmentNumber in SegmentSequence[{i}] is {N}.' +\
                f' It should be {i + 1}.\n'
            Nerrors +=  1
            logList.append(msg)
    if p2c:
        print(msg)
    
    # Check various tags in SharedFunctionalGroupsSequence.
    segIOP = seg.SharedFunctionalGroupsSequence[0]\
                .PlaneOrientationSequence[0]\
                .ImageOrientationPatient
    
    segIOP = [float(item) for item in segIOP]
    
    dcmIOP = [float(item) for item in dicom.ImageOrientationPatient]
    
    if segIOP == dcmIOP:
        msg = 'INFO:  SEG ImageOrientationPatient in ' +\
            f'SharedFunctionalGroupsSequence, {segIOP} matches ' +\
                f'DICOM ImageOrientationPatient, {dcmIOP}.\n'
    else:
        if are_items_equal_to_within_eps(segIOP, dcmIOP, epsilon):
            msg = 'INFO:  SEG ImageOrientationPatient in ' +\
                f'SharedFunctionalGroupsSequence, {segIOP}, is within ' +\
                f'{epsilon} of DICOM ImageOrientationPatient, {dcmIOP}.\n'
        else:
            msg = 'ERROR:  SEG ImageOrientationPatient in ' +\
                f'SharedFunctionalGroupsSequence, {segIOP}, is not within ' +\
                f'{epsilon} of DICOM ImageOrientationPatient, {dcmIOP}.\n'
            Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
       
    segST = seg.SharedFunctionalGroupsSequence[0]\
               .PixelMeasuresSequence[0]\
               .SliceThickness
    
    """ The SliceThickness appears to be the z-Spacing rather than the
    SliceThickness from the DICOM metadata. """
    
    if float(segST) == spacings[2]:
        msg = 'INFO:  SEG SliceThickness in SharedFunctionalGroupsSequence,' +\
            f' {segST}, matches the slice thickness calculated from DICOM ' +\
            'ImagePositionPatient and ' +\
            f'{spacings[2]}.\n'
    else:
        if are_items_equal_to_within_eps(float(segST), spacings[2], epsilon):
            msg = 'INFO:  SEG SliceThickness in ' +\
                f'SharedFunctionalGroupsSequence, {segST}, is {epsilon} of ' +\
                'the slice thickness calculated from DICOM ' +\
                'ImagePositionPatient and ImageOrientationPatient, ' +\
                f'{spacings[2]}.\n'
        else:
            print(f'\n{float(segST) - spacings[2]}')
            
            msg = 'ERROR:  SEG SliceThickness in ' +\
                f'SharedFunctionalGroupsSequence, {segST}, is not within ' +\
                f'{epsilon} of the slice thickness calculated from DICOM ' +\
                'ImagePositionPatient and ImageOrientationPatient, ' +\
                f'{spacings[2]}.\n'
            Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    segSBS = seg.SharedFunctionalGroupsSequence[0]\
                .PixelMeasuresSequence[0]\
                .SpacingBetweenSlices
    
    if float(segSBS) == spacings[2]:
        msg = 'INFO:  SEG SpacingBetweenSlices in ' +\
            f'SharedFunctionalGroupsSequence, {segSBS}, matches the slice' +\
            ' thickness calculated from DICOM ImagePositionPatient and ' +\
            f'ImageOrientationPatient, {spacings[2]}.\n'
    else:
        if are_items_equal_to_within_eps(float(segSBS), spacings[2], epsilon):
            msg = 'INFO:  SEG SpacingBetweenSlices in ' +\
                f'SharedFunctionalGroupsSequence, {segSBS}, is within ' +\
                f'{epsilon} of the slice thickness calculated from DICOM ' +\
                'ImagePositionPatient and ImageOrientationPatient, ' +\
                f'{spacings[2]}.\n'
        else:
            #print(f'\n{float(segSBS) - spacings[2]}')
            msg = 'ERROR:  SEG SpacingBetweenSlices in ' +\
                f'SharedFunctionalGroupsSequence, {segSBS}, is not within ' +\
                f'{epsilon} of the slice thickness calculated from DICOM ' +\
                'ImagePositionPatient and ImageOrientationPatient, ' +\
                f'{spacings[2]}.\n'
            Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)

    segPS = seg.SharedFunctionalGroupsSequence[0]\
               .PixelMeasuresSequence[0]\
               .PixelSpacing
    
    segPS = [float(item) for item in segPS]
    
    dcmPS = [float(item) for item in dicom.PixelSpacing]
    
    if segPS == dcmPS:
        msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence, ' +\
            f'{segPS}, matches the DICOM PixelSpacing, {dcmPS}.\n'
    else:
        if are_items_equal_to_within_eps(segPS, dcmPS, epsilon):
            msg = 'INFO:  SEG PixelSpacing in ' +\
                f'SharedFunctionalGroupsSequence, {segPS} is within ' +\
                f'{epsilon} of the DICOM PixelSpacing, {dcmPS}.\n'
        else:
            msg = 'ERROR:  SEG PixelSpacing in ' +\
                f'SharedFunctionalGroupsSequence, {segPS}, is not within ' +\
                f'{epsilon} of the DICOM PixelSpacing, {dcmPS}.\n'
            Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Verify that the number of sequences in PerFrameFunctionalGroupsSequence
    # is equal to NumberOfFrames.
    PFFGS = seg.PerFrameFunctionalGroupsSequence
    
    if len(PFFGS) == int(seg.NumberOfFrames):
        msg = 'INFO:  The number of sequences in ' +\
            f'PerFrameFunctionGroupsSequence, {len(PFFGS)}, ' +\
            f'matches NumberOfFrames {seg.NumberOfFrames}.\n'
    else:
        msg = 'ERROR:  The number of sequences in ' +\
            f'PerFrameFunctionGroupsSequence, {len(PFFGS)}, does not ' +\
            f'match NumberOfFrames {seg.NumberOfFrames}.\n'
        Nerrors +=  1
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Determine whether any of the ReferencedSOPInstanceUIDs in
    # PerFrameFunctionGroupsSequence do not match the SOPInstanceUIDs.       
    RefSOPuidsInPFFGS = [
        PFFGS[i].DerivationImageSequence[0]\
                .SourceImageSequence[0]\
                .ReferencedSOPInstanceUID for i in range(len(PFFGS))
                ]
    
    isMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInPFFGS]
    
    #numOfMatches = isMatch.count(True)
    
    # Find the indices of any non-matching ReferencedSOPInstanceUID:
    inds = [i for i, x in enumerate(isMatch) if x==False]
    
    if inds:        
        for i in range(len(inds)):
            uid = RefSOPuidsInPFFGS[inds[i]]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in ' +\
                f'PerFrameFunctionalGroupsSequence[{inds[i]}] does not ' +\
                'match any DICOM SOPInstanceUID.\n'
            logList.append(msg)
            Nerrors +=  1
    else:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs in ' +\
            'PerFrameFunctionalGroupsSequence match the DICOM ' +\
            'SOPInstanceUIDs.\n'
        logList.append(msg)
    if p2c:
        print(msg)
    
    # Determine whether the DimensionIndexValues are sensible integers,
    # and if the indexed ReferencedSOPInstanceUID agrees with the 
    # ReferencedSOPInstanceUID within the SourceImageSequence.
    DIVs = [PFFGS[i].FrameContentSequence[0]\
                    .DimensionIndexValues for i in range(len(PFFGS))]
    
    firstElements = [int(DIVs[i][0]) for i in range(len(DIVs))]
    
    # Find the indices of any elements that exceed the number of sequences in
    # SegmentSequence, numSS:
    inds = [i for i, x in enumerate(firstElements) if x > numSS]
    
    if inds:        
        for i in range(len(inds)):
            DIV = DIVs[inds[i]]
            
            msg = 'ERROR:  The first element in ' +\
                f'DimensionalIndexValue[{inds[i]}], {DIV}, exceeds the ' +\
                f'number of sequences in  SegmentSequence, {numSS}.\n'
            logList.append(msg)
            Nerrors +=  1
    else:
        msg = 'INFO:  The first element in each DimensionalIndexValue are ' +\
            'consistent with the number of sequences in SegmentSequence, ' +\
            f'{numSS}.\n'
        logList.append(msg)
    if p2c:
        print(msg)
                
    # Determine if any of the second elements in the DIVs exceed the number
    #of sequences in ReferencedInstanceSequence:
    secondElements = [int(DIVs[i][1]) for i in range(len(DIVs))]
    
    # Find the indices of any elements that exceed len(RIS):
    inds = [i for i, x in enumerate(secondElements) if x > len(RIS)]

    if inds:       
        for i in range(len(inds)):
            DIV = DIVs[inds[i]]
            
            msg = 'ERROR:  The second element in ' +\
                f'DimensionalIndexValue[{inds[i]}], {DIV}, exceeds the ' +\
                'number of sequences in ReferencedInstanceSequence, ' +\
                '{len(RIS)}.\n'
            logList.append(msg)
            Nerrors +=  1
    else:
        msg = 'INFO:  The second elements in DimensionalIndexValue matches '+\
            'the number of sequences in ReferencedInstanceSequence.\n'
        logList.append(msg)
    if p2c:
        print(msg)
    
    #print(f'\n\n\nlen(RefSOPuidsInPFFGS) = {RefSOPuidsInPFFGS}')
    #print(f'len(secondElements) = {len(secondElements)}')
    #print(f'DIVs = {DIVs}')
    
    # Determine if the ReferencedSOPInstanceUID in the 
    # ReferencedInstanceSequence indexed by the second elements in the DIVs 
    # match the ReferencedSOPInstanceUID in the SourceImageSequence.
    n = 0
    for i in range(len(secondElements)):
        refSOPinPFFGS = RefSOPuidsInPFFGS[i]
        
        """ -1 since i is zero-indexed: """
        ind = secondElements[i] - 1
        
        #print(f'ind = {ind}')
        
        refSOPinRIS = refSOPuidsInRIS[ind]
        
        if refSOPinPFFGS == refSOPinRIS:
            msg = 'INFO:  ReferencedSOPInstanceUID referenced in ' +\
                f'SourceImageSequence[{i}] matches ReferencedSOPInstanceUID' +\
                f'referenced in ReferencedInstanceSequence[{DIVs[i]}], as ' +\
                f'indexed by DimensionalIndexValue[{i}], {DIVs[i]}.\n'
        else:
            msg = 'ERROR:  ReferencedSOPInstanceUID referenced in ' +\
                f'SourceImageSequence[{i}] does not match ' +\
                'ReferencedSOPInstanceUID referenced in ' +\
                f'ReferencedInstanceSequence[{DIVs[i]}], as indexed by ' +\
                f'DimensionalIndexValue[{i}], {DIVs[i]}.\n'
            logList.append(msg)
            Nerrors +=  1
            n += 1
                  
    if n == 0:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs referenced in ' +\
            'SourceImageSequence match the ReferencedSOPInstanceUIDs ' +\
            'referenced in ReferencedInstanceSequence, as indexed by ' +\
            'DimensionalIndexValue.\n'
        logList.append(msg)
    if p2c:
        print(msg)
    
    # Determine if ImagePositionPatient in PerFrameFunctionalGroupsSequence 
    # matches ImagePositionPatient of the DICOM as indexed by the DIVs.
    IPPsInPFFGS = [PFFGS[i].PlanePositionSequence[0]\
                           .ImagePositionPatient for i in range(len(PFFGS))]
    
    IPPsInPFFGS = [[float(item) for item in IPP] for IPP in IPPsInPFFGS]
    
    c = 0
    for i in range(len(IPPsInPFFGS)):
        if RefSOPuidsInPFFGS[i] in SOPuids:
            ind = SOPuids.index(RefSOPuidsInPFFGS[i])
            
            IPP = IPPs[ind]
            
            #if IPPsInPFFGS[i] != IPP:
            if are_items_equal_to_within_eps(IPPsInPFFGS[i], IPP, epsilon):
                msg = f'INFO: ImagePositionPatient, {IPPsInPFFGS[i]}, ' +\
                    f'referenced in PlanePositionSequence[{i}] is within ' +\
                    f'{epsilon} of the DICOM ImagePositionPatient, {IPP}, ' +\
                    f'indexed by the second element, {secondElements[i]}, ' +\
                    f'in DimensionalIndexValue[{i}], {DIVs[i]}.\n'
            else:
                msg = f'ERROR: ImagePositionPatient, {IPPsInPFFGS[i]}, ' +\
                    f'referenced in PlanePositionSequence[{i+1}] is not ' +\
                    f'within {epsilon} of the DICOM ImagePositionPatient, ' +\
                    f'{IPP}, indexed by the second element, ' +\
                    f'{secondElements[i]}, in DimensionalIndexValue[{i}],' +\
                    f' {DIVs[i]}.\n'
                logList.append(msg)
                Nerrors +=  1
                c += 1
        
        if c == 0:
            msg = 'INFO: All ImagePositionPatient referenced in ' +\
                f'PlanePositionSequence are within {epsilon} of the DICOM ' +\
                'ImagePositionPatient indexed by the second element in ' +\
                'DimensionalIndexValue.\n'
            logList.append(msg)
            
        #else:
        #    msg = f'.\n'
        #    
        #    Nerrors +=  1
                  
    logList.append(msg)
    if p2c:
        print(msg)
    
    # Check that the ReferencedSegmentNumber match the first elements in the
    # DimensionIndexValues.
    RSNs = [PFFGS[i].SegmentIdentificationSequence[0]\
                    .ReferencedSegmentNumber for i in range(len(PFFGS))]
    
    diffs = [firstElements[i] - RSNs[i] for i in range(len(RSNs))]
    
    inds = [i for i, x in enumerate(diffs) if x > 0]
    
    if inds:       
        for i in range(len(inds)):
            msg = 'ERROR:  ReferencedSegmentNumber in ' +\
                'SegmentIdentificationSequence of ' +\
                f'PerFrameFunctionalGroupsSequence[{inds[i]}], ' +\
                f'{RSNs[inds[i]]}, does not match the first element in ' +\
                f'DimensionIndexValues[{inds[i]}], {DIVs[inds[i]]}.\n'
            logList.append(msg)
            Nerrors +=  1
    else:
        msg = 'INFO:  Each ReferencedSegmentNumber in ' +\
            'SegmentIdentificationSequence of ' +\
            'PerFrameFunctionalGroupsSequence matches the first element ' +\
            'in DimensionIndexValues.\n'
        logList.append(msg)
    if p2c:
        print(msg)
    
    # Determine if the dimensions of the pixel array in PixelData do not
    # match the expected dimensions (the number of frames (pixarr) has 
    # already been compared to NumberOfFrames).
    if numOfRows == seg.Rows:
        msg = f'INFO:  The number of rows in PixelData, {numOfRows}, ' +\
            f'matches the number of SEG Rows, {seg.Rows}.\n'
        logList.append(msg)
    else:
        msg = f'ERROR:  The number of rows in PixelData, {numOfRows}, ' +\
            f'does not match the number of SEG Rows, {seg.Rows}.\n'
        logList.append(msg)
        Nerrors +=  1
    if p2c:
        print(msg)
            
    if numOfCols == seg.Columns:
        msg = f'INFO:  The number of columns in PixelData, {numOfCols}, ' +\
            f'matches the number of SEG Columns, {seg.Columns}.\n'
        logList.append(msg)
    else:
        msg = f'ERROR:  The number of columns in PixelData, {numOfCols}, ' +\
            f'does not match the number of SEG Columns, {seg.Columns}.\n'
        logList.append(msg)
        Nerrors +=  1
    if p2c:
        print(msg)
    
    #print(f'There were {Nerrors} errors found in the SEG.\n')
    
    if p2c:
        print('-'*120)
            
    return logList, Nerrors