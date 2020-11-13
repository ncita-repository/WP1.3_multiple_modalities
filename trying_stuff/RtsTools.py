# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:55:17 2020

@author: ctorti
"""




# Import packages and functions:
import os
import time
from copy import deepcopy


"""
******************************************************************************
******************************************************************************
DICOM RTSTRUCT FUNCTIONS
******************************************************************************
******************************************************************************
"""


def GetCIStoDcmInds(RtsRoi, SOPuids):
    """
    Get the DICOM slice numbers that correspond to each Contour Image Sequence
    in an RTSTRUCT ROI.
    
    Inputs:
        RtsRoi       - ROI Object from an RTSTRUCT file
        
        SOPuids      - List of strings of the SOP UIDs of the DICOMs
        
    Returns:
        CIStoDcmInds - List of integers of the DICOM slice numbers that
                       correspond to each Contour Image Sequence in RtsRoi
    """
    
    CIStoDcmInds = [] 
    
    # Loop through each Contour Image Sequence:
    sequences = RtsRoi.ReferencedFrameOfReferenceSequence[0]\
                      .RTReferencedStudySequence[0]\
                      .RTReferencedSeriesSequence[0]\
                      .ContourImageSequence
    
    for sequence in sequences:
        # Get the Reference SOP Instance UID:
        uid = sequence.ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in SOPuids:
        CIStoDcmInds.append(SOPuids.index(uid))
        
    return CIStoDcmInds




def GetCStoDcmInds(RtsRoi, SOPuids):
    """
    Get the DICOM slice numbers that correspond to each Contour Sequence
    in an RTSTRUCT ROI.
    
    Inputs:
        RtsRoi      - ROI Object from an RTSTRUCT file
        
        SOPuids     - List of strings of the SOP UIDs of the DICOMs
        
    Returns:
        CStoDcmInds - List of integers of the DICOM slice numbers that 
                      correspond to each Contour Sequence in RtsRoi 
    """
    
    CStoDcmInds = [] 
    
    # Loop through each Contour Image Sequence:
    sequences = RtsRoi.ROIContourSequence[0].ContourSequence
    
    for sequence in sequences:
        # Get the Reference SOP Instance UID:
        uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in SOPuids:
        CStoDcmInds.append(SOPuids.index(uid))
        
    return CStoDcmInds




def VerifyCISandCS(CIStoDcmInds, CStoDcmInds, LogToConsole=False):
    """
    Verify that the ContourImageSequence-to-DICOM slice indeces and the
    ContourSequence-to-DICOM slice indeces are the same.
    
    Inputs:
        CIStoDcmInds - List of integers of the DICOM slice numbers that 
                      correspond to each Contour Image Sequence
        
        CStoDcmInds  - List of integers of the DICOM slice numbers that 
                      correspond to each Contour Sequence
        
    Returns:
        (Boolean) True/False if they match/do not match  
    """
    
    if CIStoDcmInds != CStoDcmInds:
        if LogToConsole:
            print('\nThe ContourImageSequence-to-DICOM slice indices are not',
                  'the same as the ContourSequence-to-DICOM slice indices.')
            
        return False
    
    else:
        return True
    
    






def AddToRtsSequences(RtsRoi):
    """
    Append the last item in the following sequences in the RTS Object:
        Contour Image Sequence
        Contour Sequence
    
    Inputs:
        RtsRoi    - ROI Object from an RTSTRUCT file
        
    Returns:
        NewRtsRoi - Modified ROI Object with sequences lengthened by one
    """
    
    # Use RtsRoi as a template for NewRtsRoi: 
    NewRtsRoi = deepcopy(RtsRoi)
        
    # The last item in Contour Image Sequence:
    last = deepcopy(RtsRoi.ReferencedFrameOfReferenceSequence[0]\
                          .RTReferencedStudySequence[0]\
                          .RTReferencedSeriesSequence[0]\
                          .ContourImageSequence[-1])
    
    # Append to the Contour Image Sequence:
    NewRtsRoi.ReferencedFrameOfReferenceSequence[0]\
             .RTReferencedStudySequence[0]\
             .RTReferencedSeriesSequence[0]\
             .ContourImageSequence\
             .append(last)

    # The last item in Contour Sequence:
    last = deepcopy(RtsRoi.ROIContourSequence[0]\
                          .ContourSequence[-1])
    
    # Append to the Contour Sequence:
    NewRtsRoi.ROIContourSequence[0]\
             .ContourSequence\
             .append(last)
             
    return NewRtsRoi






def ModifyRtsTagVals(OrigRtsRoi, NewRtsRoi, FromSliceNum, ToSliceNum, Dicoms,
                     OrigCIStoDcmInds, NewCIStoDcmInds, LogToConsole=False):
    """
    Modify various tag values of the new RTS ROI so that they are consistent
    with the corresponding tag values of the DICOMs (e.g. SOPClassUID, 
    SOPInstanceUID), and so that the ROIContourSequence values are consistent
    with those of the original RTS ROI for existing contours, and with the
    contour that has been copied from the original RTS ROI.
         
    
    Inputs:
        OrigRtsRoi       - Original ROI Object 
        
        NewRtsRoi        - New/modified ROI Object
        
        FromSliceNum     - (Integer) Slice index in the DICOM stack
                           corresponding to the contour to be copied 
                           (counting from 0)
                          
        ToSliceNum       - (Integer) Slice index in the DICOM stack
                           where the contour will be copied to (counting from 
                           0)
                          
        Dicoms           - A list of DICOM Objects that include the DICOMs that
                           NewRtsRoi relate to
                          
        OrigCIStoDcmInds - List of integers of the DICOM slice numbers that 
                           correspond to each Contour Image Sequence in 
                           OrigRtsRoi
                          
        NewCIStoDcmInds  - List of integers of the DICOM slice numbers that 
                           correspond to each Contour Image Sequence in 
                           NewRtsRoi
                           
        LogToConsole     - (Boolean) Denotes whether or not some intermediate
                            results will be logged to the console
        
    Returns:
        NewRtsRoi        - Modified new ROI Object for the Target series
    """
    
    # Get the sequence number that corresponds to the slice number to be 
    # copied, and the slice number to be copied to:
    FromOrigSeqNum = OrigCIStoDcmInds.index(FromSliceNum)
    ToNewSeqNum = NewCIStoDcmInds.index(ToSliceNum)
    
    if LogToConsole:
            print(f'\nFromOrigSeqNum = {FromOrigSeqNum} -->',
                  f'OrigCIStoDcmInds[{FromOrigSeqNum}] =',
                  f'{OrigCIStoDcmInds[FromOrigSeqNum]}')
            
            print(f'ToNewSeqNum   = {ToNewSeqNum} -->',
                  f'NewCIStoDcmInds[{ToNewSeqNum}] =',
                  f'{NewCIStoDcmInds[ToNewSeqNum]}')
    
    # Loop through each index in NewCIStoDcmInds and modify the values:
    for i in range(len(NewCIStoDcmInds)):
        # The DICOM slice index for the i^th sequence:
        s = NewCIStoDcmInds[i]
        
        if LogToConsole:
            print(f'\nNewCIStoDcmInds[{i}] = {NewCIStoDcmInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Class UID and Referenced SOP Instance UID 
        # so that they are consistent with the SOP Class UID and SOP Instance  
        # UID of their corresponding DICOMs:
        """ SOPClassUIDs will need to be modified for cases where the contour
        is being copied across different modality images. Not required for
        copying within the same series:
        NewRtsRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewRtsRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        """ Ditto:
        NewRtsRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourImageSequence[0]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewRtsRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourImageSequence[0]\
                 .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
                 
        
        
        # Modify the Contour Number to match the value of the Source contour
        #  being copied:
        NewRtsRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourNumber = f"{i + 1}"
                     
                     
        # Modify select values in ROIContourSequence:         
        if i == ToNewSeqNum: 
            if LogToConsole:
                print('\nThis is a new sequence.')
                print(f'Copying the {FromOrigSeqNum}th sequence from Source',
                      f'to the {i}th sequence in Target.')
                
            # This is the new sequence. Copy from FromSeqNum.
            # Modify the Number of Contour Points to match the value of the 
            # Source contour being copied:
            val = deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                     .ContourSequence[FromOrigSeqNum]\
                                     .NumberOfContourPoints)
            
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .NumberOfContourPoints = deepcopy(val)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            val = deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                     .ContourSequence[FromOrigSeqNum]\
                                     .ContourData)
            
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .ContourData = deepcopy(val)
        
        else:
            # This is an existing sequence. Find the index of OrigCIStoDcmInds 
            # for this sequence number (of NewCIStoDcmInds):
            ind = OrigCIStoDcmInds.index(NewCIStoDcmInds[i])
            
            if LogToConsole:
                print('\nThis is an existing sequence.')
                print(f'Copying the {ind}th sequence from Source to the {i}th',
                      f'sequence in Target.')
                                          
            # Modify the Number of Contour Points to match the value of the 
            # Source contour being copied:
            val = deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                     .ContourSequence[ind]\
                                     .NumberOfContourPoints)
            
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .NumberOfContourPoints = deepcopy(val)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            val = deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                     .ContourSequence[ind]\
                                     .ContourData)
            
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .ContourData = deepcopy(val)
                                          
    return NewRtsRoi









def ExportRtsRoi(TrgRtsRoi, SrcRtsFpath, NamePrefix, ExportDir):
    """
    Export contour ROI Object to disc.
    
    
    Inputs:
        TrgRtsRoi   - Target RT-STRUCT ROI object to be exported
        
        SrcRtsFpath - (String) Full path of the Source DICOM-RTSTRUCT file
                      (used to generate the filename of the new RTS file)
        
        NamePrefix  - (String) Prefix to be added to the assigned filename and
                      ROI Name (e.g. 'MR4_S9_s23_to_MR4_S9_s22') 
                            
        ExportDir   - (String) Directory where the RTSTRUCT is to be exported
                           
    
                            
    Returns:
        TrgRtsFpath - (String) Full path of the exported Target DICOM-RTSTRUCT 
                      file
    
    """
    
    """
    The following was moved into the main function CopyContourWithinSeries:
    # Generate a new SOP Instance UID:
    RtsRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    RtsRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    """
    
    # Get the filename of the original RT-Struct file:
    SrcRtsFname = os.path.split(SrcRtsFpath)[1]
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    TrgRtsRoi.StructureSetDate = NewDate
    TrgRtsRoi.StructureSetTime = NewTime
    
    FnamePrefix = NewDate + '_' + NewTime + '_' + NamePrefix + '_from_' 
    RoiNamePrefix = NamePrefix + '_from_'
    
    # Modify the Structure Set Label (this appears under Name in XNAT):
    #TrgRtsRoi.StructureSetLabel = 'Copy_of_' + TrgRtsRoi.StructureSetLabel
    TrgRtsRoi.StructureSetLabel = RoiNamePrefix + TrgRtsRoi.StructureSetLabel
    
    # Modify the ROI Name for the Structure Set ROI Sequence:             
    TrgRtsRoi.StructureSetROISequence[0].ROIName = RoiNamePrefix \
                                                  + TrgRtsRoi.StructureSetLabel
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = 'New_' + SrcRtsFname + '_' + NewDate + '_' + NewTime
    TrgRtsFname = FnamePrefix + SrcRtsFname
    
    TrgRtsFpath = os.path.join(ExportDir, TrgRtsFname)
    
    TrgRtsRoi.save_as(TrgRtsFpath)
        
    print('\nRTS ROI exported to:\n\n', TrgRtsFpath)
    
    return TrgRtsFpath