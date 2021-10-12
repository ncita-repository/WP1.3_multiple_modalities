# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 18:21:35 2020

@author: ctorti
"""

def ErrorCheckSeg(Seg, DicomDir, LogToConsole=False):
    """
    Check a SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the SEG.  
    
    Inputs:
    ------
    
    Seg : Pydicom object
        SEG object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Seg.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
                           
    
    Outputs:
    -------
    
    ErrorList : list of strings
        A list of strings that describe any errors that are found.  If no 
        errors are found an empty list ([]) will be returned.
    """
    
    import importlib
    import GeneralTools
    importlib.reload(GeneralTools)
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from DicomTools import ImportDicom
    from ImageTools import GetImageAttributes
    from GeneralTools import AreItemsEqualToWithinEpsilon
    from GeneralTools import AreListsEqualToWithinEpsilon
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    Dicom = ImportDicom(DicomDir)
    
    #Size, Spacings, ST,\
    #IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    Size, Spacings, ST,\
    IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='pydicom')
    
    # The number of frames, rows and columns in the SEG's pixel array:
    if len(Seg.pixel_array.shape) > 2:
        PixArrF, PixArrR, PixArrC = Seg.pixel_array.shape
    else:
        PixArrR, PixArrC = Seg.pixel_array.shape
        PixArrF = 1
    
    #NumOfDicoms = len(SOPuids)
    
    ErrorList = []
    
    """Determine whether the number of sequences in 
    ReferencedInstanceSequence matches the number of DICOMs."""
    
    RIS = deepcopy(Seg.ReferencedSeriesSequence[0]\
                      .ReferencedInstanceSequence)
    
    if len(SOPuids) != len(RIS):
        msg = f'The number of "SOPInstanceUID"s ({len(SOPuids)}) does not'\
              + ' match the number of sequences in '\
              + f'"ReferencedInstanceSequence" ({len(RIS)}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs do not
    match the SOPInstanceUIDs."""        
    RefSOPuidsInRIS = [RIS[i].ReferencedSOPInstanceUID for i in range(len(RIS))]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInRIS]
    
    NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if NumOfMatches == 0:
        NonMatchingUids = [RefSOPuidsInRIS[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'"ReferencedSOPInstanceUID" {uid} in '\
                  + f'"ReferencedInstanceSequence" {i+1} does not match any '\
                  + '"SOPInstanceUID".'
        
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
                
                
    
    """Determine whether any of the SOPInstanceUIDs do not match the
    ReferencedSOPInstanceUIDs."""        
    
    # Determine whether the SOP UIDs match the Referenced SOP UIDs:
    IsMatch = [SOPuid in RefSOPuidsInRIS for SOPuid in SOPuids]
    
    NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if NumOfMatches == 0:
        NonMatchingUids = [SOPuids[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'"SOPInstanceUID" {uid} is not referenced in any '\
                  + '"ReferencedSOPInstanceUID" in '\
                  + '"ReferencedInstanceSequence".'
        
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
    
    
                           
    if Seg.ReferencedSeriesSequence[0].SeriesInstanceUID != Dicom.SeriesInstanceUID:
        msg = f'SEG "SeriesInstanceUID" in "ReferencedSeriesSequence" '\
              + f'({Seg.PatientName}) does not match DICOM "SeriesInstanceUID"'\
              + f'({Dicom.SeriesInstanceUID}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
            
       
    if Seg.PatientName != Dicom.PatientName:
        msg = f'SEG "PatientName" ({Seg.PatientName}) does not match DICOM '\
              + f'"PatientName" ({Dicom.PatientName}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
        
    if Seg.PatientID != Dicom.PatientID:
        msg = f'SEG "PatientID" ({Seg.PatientID}) does not match the DICOM '\
              + f'"PatientID" ({Dicom.PatientID}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)

    if Seg.StudyInstanceUID != Dicom.StudyInstanceUID:
        msg = f'SEG "StudyInstanceUID" ({Seg.StudyInstanceUID}) does not '\
              + 'match the DICOM tag "StudyInstanceUID" '\
              + f'({Dicom.StudyInstanceUID}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
    #if Seg.StudyID != Dicom.StudyID:
    #    msg = f'The SEG "StudyID" ({Seg.StudyID}) does not match the '\
    #          + f'DICOM "StudyID" ({Dicom.StudyID}).'
    #    
    #    ErrorList.append(msg)
    #    
    #    if LogToConsole:
    #        print('\n' + msg)
            
    if Seg.FrameOfReferenceUID != Dicom.FrameOfReferenceUID:
        msg = 'SEG "FrameOfReferenceUID" ({Seg.FrameOfReferenceUID}) '\
              + 'does not match DICOM "FrameOfReferenceUID" '\
              + f'({Dicom.FrameOfReferenceUID}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
    if int(Seg.NumberOfFrames) != PixArrF:
        msg = f'"NumberOfFrames" ({Seg.NumberOfFrames}) does not match the '\
              + f'number of frames in "PixelData" ({PixArrF}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
    if Seg.Rows != Dicom.Rows:
        msg = f'SEG "Rows" ({Seg.Rows}) does not match DICOM "Rows" '\
              + f'({Dicom.Rows}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
    if Seg.Columns != Dicom.Columns:
        msg = f'SEG "Columns" ({Seg.Columns}) does not match DICOM "Columns" '\
              + f'({Dicom.Columns}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
    """There should only be one sequence in SegmentSequence."""
    N = len(Seg.SegmentSequence)
    
    if N > 1:
        msg = 'There are {N} sequences in "SegmentSequence".  There should '\
              + f'be only 1.'
              
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
    
    """Verify that the SegmentNumber in SegmentSequence is 1."""
    N = deepcopy(Seg.SegmentSequence[0].SegmentNumber)
    
    if N != 1:
        msg = '"SegmentNumber" in "SegmentSequence" is {N}.  It should be 1.'
              
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
            
    """Check various tags in SharedFunctionalGroupsSequence."""
    SegIOP = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                         .PlaneOrientationSequence[0]\
                         .ImageOrientationPatient)
    
    SegIOP = [float(item) for item in SegIOP]
    
    DcmIOP = [float(item) for item in Dicom.ImageOrientationPatient]
    
    #if SegIOP != Dicom.ImageOrientationPatient:
    if not AreListsEqualToWithinEpsilon(SegIOP, DcmIOP):
        msg = 'SEG "ImageOrientationPatient" in '\
              + f'"SharedFunctionalGroupsSequence" ({SegIOP}) does not match '\
              + f'DICOM "ImageOrientationPatient" ({DcmIOP}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
       
    
    SegST = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                        .PixelMeasuresSequence[0]\
                        .SliceThickness)
    
    """ The SliceThickness appears to be the z-Spacing rather than the
    SliceThickness from the DICOM metadata. """
    
    #if float(SegST) != Spacings[2]:
    if not AreItemsEqualToWithinEpsilon(float(SegST), Spacings[2], epsilon=1e-05):
        print(f'\n{float(SegST) - Spacings[2]}')
        
        msg = 'SEG "SliceThickness" in "SharedFunctionalGroupsSequence" '\
              + f'({SegST}) does not match the slice thickness calculated '\
              + 'from DICOM "ImagePositionPatient" and '\
              + f'"ImageOrientationPatient" ({Spacings[2]}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
    SegSBS = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                         .PixelMeasuresSequence[0]\
                         .SpacingBetweenSlices)
    
    #if float(SegSBS) != Spacings[2]:
    if not AreItemsEqualToWithinEpsilon(float(SegSBS), Spacings[2], epsilon=1e-05):
        print(f'\n{float(SegSBS) - Spacings[2]}')
        
        msg = 'SEG "SpacingBetweenSlices" in "SharedFunctionalGroupsSequence"'\
              + f' ({SegSBS}) does not match the slice thickness calculated '\
              + 'from the DICOM "ImagePositionPatient" and '\
              + f'"ImageOrientationPatient" ({Spacings[2]}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)

    
    SegPS = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                        .PixelMeasuresSequence[0]\
                        .PixelSpacing)
    
    SegPS = [float(item) for item in SegPS]
    
    DcmPS = [float(item) for item in Dicom.PixelSpacing]
    
    #if SegPS != Dicom.PixelSpacing:
    if not AreListsEqualToWithinEpsilon(SegPS, DcmPS):
        msg = 'SEG "PixelSpacing" in "SharedFunctionalGroupsSequence '\
              + f'({SegPS}) does not match DICOM "PixelSpacing" ({DcmPS}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
    
    
    """Verify that the number of sequences in PerFrameFunctionalGroupsSequence
    is equal to NumberOfFrames."""
    PFFGS = deepcopy(Seg.PerFrameFunctionalGroupsSequence)
    
    if len(PFFGS) != int(Seg.NumberOfFrames):
        msg = f'The number of sequences in "PerFrameFunctionGroupsSequence" '\
              + f'({len(PFFGS)}) does not match "NumberOfFrames" '\
              + f'({Seg.NumberOfFrames}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs in the
    PerFrameFunctionGroupsSequence do not match the SOPInstanceUIDs."""        
    RefSOPuidsInPFFGS = [PFFGS[i].DerivationImageSequence[0]\
                                 .SourceImageSequence[0]\
                                 .ReferencedSOPInstanceUID for i in range(len(PFFGS))]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInPFFGS]
    
    NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if NumOfMatches == 0:
        NonMatchingUids = [RefSOPuidsInPFFGS[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'"ReferencedSOPInstanceUID" {uid} in '\
                  + '"PerFrameFunctionalGroupsSequence" {i+1} does not match '\
                  + 'any DICOM "SOPInstanceUID".'
        
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
                
    
    
    """Determine whether the DimensionIndexValues are sensible integers,
    and if the indexed ReferencedSOPInstanceUID agrees with the 
    ReferencedSOPInstance UID within the SourceImageSequence."""
    DIVs = [PFFGS[i].FrameContentSequence[0]\
                    .DimensionIndexValues for i in range(len(PFFGS))]
    
    # Determine whether any of the first elements in the DIVs are not 1  
    # (there should only be one segment, so the first element should always
    # be 1):
    FirstElements = [int(DIVs[i][0]) for i in range(len(DIVs))]
    
    #NumOfMatches = FirstElements.count('1')
    
    # Find the indices of any non-matching elements:
    Inds = [i for i, x in enumerate(FirstElements) if x != 1]
    
    if Inds:
        divs = [DIVs[i] for i in Inds]
        
        for i in range(len(Inds)):
            msg = f'The first element in "DimensionalIndexValue" {i+1} '\
                  + f'({divs[i]}) is not "1".'
        
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
                
    # Determine if any of the second elements in the DIVs exceed the number
    # of sequences in ReferencedInstanceSequence:
    SecondElements = [int(DIVs[i][1]) for i in range(len(DIVs))]
    
    # Find the indices of any elements that exceed len(RIS):
    Inds = [i for i, x in enumerate(SecondElements) if x > len(RIS)]

    if Inds:
        #divs = [DIVs[i] for i in Inds]
        
        for i in range(len(Inds)):
            msg = f'The second element in "DimensionalIndexValue" {i+1} '\
                  + f'({DIVs[i]}) exceeds the number of sequences in '\
                  + f'"ReferencedInstanceSequence" (N = {len(RIS)}).'
        
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
                
    """Determine if the ReferencedSOPInstanceUID in the 
    ReferencedInstanceSequence indexed by the second elements in the DIVs 
    match the ReferencedSOPInstanceUID in the SourceImageSequence."""
    for i in range(len(SecondElements)):
        RefSOPinPFFGS = RefSOPuidsInPFFGS[i]
        
        # Need to -1 since i is zero-indexed:
        ind = SecondElements[i] - 1
        
        RefSOPinRIS = RefSOPuidsInRIS[ind]
        
        if RefSOPinPFFGS != RefSOPinRIS:
            msg = '"ReferencedSOPInstanceUID" referenced in '\
                  + f'"SourceImageSequence" {i+1} does not match '\
                  + '"ReferencedSOPInstanceUID" referenced in '\
                  + '"ReferencedInstanceSequence" as indexed by the second '\
                  + f'element in "DimensionalIndexValue" {i+1} ({DIVs[i]}).'
                  
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
    
    
    """Determine if ImagePositionPatient in PerFrameFunctionalGroupsSequence 
    matches ImagePositionPatient of the DICOM as indexed by the DIVs."""
    IPPsInPFFGS = [PFFGS[i].PlanePositionSequence[0]\
                           .ImagePositionPatient for i in range(len(PFFGS))]
    
    # Convert from strings to floats:
    IPPsInPFFGS = [[float(item) for item in IPP] for IPP in IPPsInPFFGS]
    
    for i in range(len(IPPsInPFFGS)):
        ind = SOPuids.index(RefSOPuidsInPFFGS[i])
        
        IPP = IPPs[ind]
        
        #if IPPsInPFFGS[i] != IPP:
        if not AreListsEqualToWithinEpsilon(IPPsInPFFGS[i], IPP):
            msg = f'"ImagePositionPatient" ({IPPsInPFFGS[i]}) referenced in '\
                  + f'"PlanePositionSequence" {i+1} does not match DICOM '\
                  + f'"ImagePositionPatient" ({IPP}) as indexed by the second'\
                  + f' element ({SecondElements[i]}) in '\
                  + f'"DimensionalIndexValue" {i+1} ({DIVs[i]}).'
                  
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
    
    
    """Determine whether any of the ReferencedSegmentNumber are not 1 (there
    should only be one segment)."""
    RSNs = [PFFGS[i].SegmentIdentificationSequence[0]\
                    .ReferencedSegmentNumber for i in range(len(PFFGS))]
    
    # Find the indices of any non-matching elements:
    Inds = [i for i, x in enumerate(RSNs) if x != 1]
    
    if Inds:        
        for i in range(len(Inds)):
            msg = f'"ReferencedSegmentNumber" {i+1} is not 1.'\
        
            ErrorList.append(msg)
            
            if LogToConsole:
                print('\n' + msg)
    
    
    """Determine if the dimensions of the pixel array in PixelData do not
    match the expected dimensions (the number of frames (PixArrF) has 
    already been compared to NumberOfFrames)."""
    if PixArrR != Seg.Rows:
        msg = f'The number of rows in "PixelData" ({PixArrR}) does not '\
              + f'match the number of SEG "Rows" ({Seg.Rows}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
            
    if PixArrC != Seg.Columns:
        msg = f'The number of columns in "PixelData" ({PixArrC}) does not'\
              + f' match the number of SEG "Columns" ({Seg.Columns}).'
        
        ErrorList.append(msg)
        
        if LogToConsole:
            print('\n' + msg)
    
    
    print(f'\nThere were {len(ErrorList)} errors found in the SEG.')
            
    return ErrorList