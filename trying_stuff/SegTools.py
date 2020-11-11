# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:54:12 2020

@author: ctorti
"""




# Import packages and functions:
import os
import time
import copy
import pydicom
from pydicom.pixel_data_handlers.numpy_handler import pack_bits
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
from DicomTools import GetDicomSOPuids




"""
******************************************************************************
******************************************************************************
SEGMENT FUNCTIONS
******************************************************************************
******************************************************************************
"""      


def GetSegLabels(SegRoi):
    """
    Get the list of segment labels in a SEG ROI.
    
    Inputs:
        SegRoi    - ROI Object from a SEG file
    
        
    Returns:
        SegLabels - (List of strings) Segment labels
    """
    
    SegLabels = [] 
    
    sequences = SegRoi.SegmentSequence
    
    for sequence in sequences:
        label = sequence.SegmentLabel
        
        SegLabels.append(label)
        
    return SegLabels






def GetRSOPuidsInRIS(SegRoi):
    """
    Get the list of referenced SOP Instance UIDs in the Referenced Instance
    Sequence of a SEG ROI.
    
    Inputs:
        SegRoi   - ROI Object from a SEG file
    
        
    Returns:
        RSOPuids - (List of strings) Referenced SOP UIDs
    """
    
    RSOPuids = [] 
    
    sequences = SegRoi.ReferencedSeriesSequence[0]\
                      .ReferencedInstanceSequence
    
    for sequence in sequences:
        uid = sequence.ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids





def GetRIStoDcmInds(Roi, SOPuids):
    """
    Get the DICOM slice numbers that correspond to each Referenced Instance
    Sequence in a RTSTRUCT or SEG ROI.
    
    Inputs:
        Roi          - ROI Object from a RTS/SEG file
        
        SOPuids      - (List of strings) SOP UIDs of the DICOMs
        
    Returns:
        RIStoDcmInds - (List of integers) DICOM slice numbers that correspond
                       to each Referenced Instance Sequence in Roi
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    RSOPuids = GetRSOPuidsInRIS(Roi)
    
    RIStoDcmInds = [] 
    
    for Ruid in RSOPuids:
        # Find the matching index of Ruid in SOPuids:
        RIStoDcmInds.append(SOPuids.index(Ruid))
        
    return RIStoDcmInds





def GetRSOPuidsInPFFGS(SegRoi):
    """
    Get the list of referenced SOP Instance UIDs in the Per-Frame Functional
    Groups Sequence of a SEG ROI.
    
    Inputs:
        SegRoi   - ROI Object from a SEG file
    
        
    Returns:
        RSOPuids - (List of strings) Referenced SOP UIDs
    """
    
    RSOPuids = [] 
    
    sequences = SegRoi.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        uid = sequence.DerivationImageSequence[0]\
                      .SourceImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids





def GetPFFGStoDcmInds(SegRoi, SOPuids):
    """
    Get the DICOM slice numbers that correspond to each Per-frame Functional
    Groups Sequence in a SEG ROI.
    
    Inputs:
        SegRoi         - ROI Object from a SEG file
        
        SOPuids        - (List of strings) SOP UIDs of the DICOMs
        
    Returns:
        PFFGStoDcmInds - (List of integers) DICOM slice numbers that correspond
                         to each Per-frame Functional Groups Sequence in SegRoi
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    RSOPuids = GetRSOPuidsInPFFGS(SegRoi)
    
    PFFGStoDcmInds = [] 
    
    for Ruid in RSOPuids:
        # Find the matching index of Ruid in SOPuids:
        PFFGStoDcmInds.append(SOPuids.index(Ruid))
        
    return PFFGStoDcmInds





def GetDIVs(SegRoi):
    """
    Get the DimensionIndexValues in a SEG ROI.
    
    Inputs:
        SegRoi - ROI Object from a SEG file
        
        
    Returns:
        DIVs   - (List of lists of integers) Dimension Index Values
        
    """
    
    # Initialise the list of DIVs:
    DIVs = [] 
    
    sequences = SegRoi.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        DIV = sequence.FrameContentSequence[0]\
                      .DimensionIndexValues
        
        DIVs.append(DIV)
        
    return DIVs






def GroupListBySegment(ListToGroup, DIVs):
    """
    Group a list of items by common segment/ROI using the Dimension Index 
    Values.
    
    Inputs:
        ListToGroup - List to be grouped
        
        DIVs        - (List of lists of integers) List of Dimension Index 
                      Values
        
        
    Returns:
        GroupedList - (List of lists) ListToGroup grouped by ROI
        
        
    Note:
        DIVs is a list of lists of the form:
            
        [ [1, 10], [1, 11], [1, 12], ... [2, 6], [2, 7], [2, 8], ... ]
        
        where each list of index values that belong to a common ROI are 
        separated into separate lists.
        
        The goal is to group a list (of anything) by ROI.  In the case of the
        DIVs, that is:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ],
          [ [2, 6], [2, 7], [2, 8], ... ],
          ... 
         ]
          
         
        If there's only one ROI, ListToGroup will still be nested in a list to 
        maintain continuity.  In the case of the DIVs with only one ROI:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ]
         ]
    """
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in DIVs:
        RoiNums.append(DIV[0])
        
        
    UniqueRoiNums = list(set(RoiNums))
    
    
    if len(UniqueRoiNums) == 1:
        GroupedList = [ListToGroup]
        
    else:
        GroupedList = []
        
        for RoiNum in UniqueRoiNums:
            ListThisRoi = [] # initialise
            
            for i in range(len(ListToGroup)):
                if DIVs[i][0] == RoiNum:
                    ListThisRoi.append(ListToGroup[i])
        
            GroupedList.append(ListThisRoi)
             
    return GroupedList




def AddToSegSequences_OLD(SegRoi):
    """
    Append the last item in the following sequences in the SEG Object: 
        Referenced Instance Sequence
        Per-frame Functional Groups Sequence
    
    Inputs:
        SegRoi    - ROI Object from a SEG file
        
    Returns:
        NewSegRoi - Modified ROI Object with sequences lengthened by one
    """
    
    # Use SegRoi as a template for NewSegRoi: 
    NewSegRoi = copy.deepcopy(SegRoi)
        
    # The last item in Referenced Instance Sequence:
    last = copy.deepcopy(SegRoi.ReferencedSeriesSequence[0]\
                               .ReferencedInstanceSequence[-1])
    
    # Append to the Referenced Instance Sequence:
    NewSegRoi.ReferencedSeriesSequence[0]\
             .ReferencedInstanceSequence\
             .append(last)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = copy.deepcopy(SegRoi.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    NewSegRoi.PerFrameFunctionalGroupsSequence\
             .append(last)
             
    return NewSegRoi





def AddToPFFGS(SegRoi):
    """
    Append the last item in the Per-Frame Functional Groups Sequence of the
    SEG Object to itself (i.e. increase the length of the sequence by one). 
    
    Inputs:
        SegRoi    - ROI Object from a SEG file
        
    Returns:
        NewSegRoi - Modified ROI Object with sequences lengthened by one
    """
    
    # Use SegRoi as a template for NewSegRoi: 
    NewSegRoi = copy.deepcopy(SegRoi)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = copy.deepcopy(SegRoi.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    NewSegRoi.PerFrameFunctionalGroupsSequence\
             .append(last)
             
    return NewSegRoi






def PixArr2Labmap(PixelArray, NumOfSlices, FrameToSliceInds):
    """
    Convert a 3D SEG pixel array to a 3D SEG labelmap.  The labelmap will 
    contain the 2D pixel arrays at the same frame position as the slice
    position of the corresponding DICOM series (with 2D zero masks in-between).
    
    Inputs:
        PixelArray       - (Numpy array) The SEG's PixelData loaded as a 
                           Numpy array with dtype uint8
        
        NumOfSlices      - (Integer) The number of slices in the DICOM series
        
        FrameToSliceInds - (List of integers) The DICOM slice numbers that
                           correspond to each frame in PixelArray,
                           e.g. PFFGStoDcmInds 
                           = PerFrameFunctionalGroupsSequence-to-DICOM slice
                           number
        
    Returns:
        LabelMap         - (Numpy array) A zero-padded version of PixelArray
                           with dtype uint

    """
    
    #print('\nRunning PixArr2Labmap()...')
    
    NumOfFrames, NumOfRows, NumOfCols = PixelArray.shape
    
    #print(f'\n\nThe maximum value in PixelArray is {PixelArray.max()}')
    
    # Initialise LabelMap:
    #LabelMap = np.zeros((NumOfSlices, NumOfRows, NumOfCols), dtype='bool')
    LabelMap = np.zeros((NumOfSlices, NumOfRows, NumOfCols), dtype='uint')
    
    # Over-write all frames in LabelMap with non-zero 2D masks from PixelArray:
    for i in range(len(FrameToSliceInds)):
        # The DICOM slice number for the i^th frame:
        s = FrameToSliceInds[i]
        
        #print(f'\nThe DICOM slice number for the {i}^th frame is {s}')
        
        LabelMap[s] = PixelArray[i]
        
    
    #print(f'\n\nThe maximum value in LabelMap is {LabelMap.max()}')
    
    #print('\nCompleted running PixArr2Labmap.')
    
    return LabelMap





def PixArr2Image(PixelArray, Image, FrameToSliceInds):
    """
    Convert a 3D SEG pixel array to a 3D SimpleITK image using an intermediate
    conversion from 3D SEG pixel array to a Numpy array (using PixArr2Labmap).  
    
    Inputs:
        PixelArray       - (Numpy array) The SEG's PixelData loaded as a 
                           Numpy array with dtype uint8
                           
        Image            - (SimpleITK image) The 3D image that PixelArray
                           relates to
        
        FrameToSliceInds - (List of integers) The DICOM slice numbers that
                           correspond to each frame in PixelArray,
                           e.g. PFFGStoDcmInds 
                           = PerFrameFunctionalGroupsSequence-to-DICOM slice
                           number
        
    Returns:
        LabelMapImage    - (SimpleITK image) A zero-padded version of 
                           PixelArray

    """
    
    #print('\nRunning PixArr2Image()...')
    
    ImSize = Image.GetSize()
    
    nda = PixArr2Labmap(PixelArray, ImSize[2], FrameToSliceInds)
    
    #print('\nContinuing with running PixArr2Image()..')
    
    LabelMapImage = sitk.GetImageFromArray(nda)
    
    # Set the Origin, Spacing and Direction of LabelMapImage to that of Image:
    LabelMapImage.SetOrigin(Image.GetOrigin())
    LabelMapImage.SetSpacing(Image.GetSpacing())
    LabelMapImage.SetDirection(Image.GetDirection())
    
    #print(f'\n\nThe maximum value in nda is {nda.max()}')
    #print(f'\n\nThe maximum value in LabelMapImage is {ImageMax(LabelMapImage)}')
    
    #print('\nCompleted running PixArr2Image().')
    
        
    return LabelMapImage




def Image2PixArr(LabmapImage):
    """
    Convert a 3D SimpleITK image to a 3D SEG pixel array.  
    
    Inputs:
        LabelMapImage    - (SimpleITK image) A zero-padded version of 
                           PixelArray
        
        
    Returns:
        PixelArray       - (Numpy array) The SEG's PixelData as a Numpy array
                           with dtype uint8
                           
        FrameToSliceInds - (List of integers) The DICOM slice numbers that
                           correspond to each frame in PixelArray,
                           e.g. PFFGStoDcmInds 
                           = PerFrameFunctionalGroupsSequence-to-DICOM slice
                           number

    """
    
    # Convert from SimpleITK image to a Numpy array:
    LabmapNda = sitk.GetArrayFromImage(LabmapImage)
    
    FrameToSliceInds = []
    
    # Loop through all frames in LabmapNda:
    for i in range(LabmapNda.shape[0]):
        Sum = np.amax(LabmapNda[i])
        
        if Sum:
            # This is a non-zero frame:
            FrameToSliceInds.append(i)
            
            
    
    # Initialise PixelArray:
    PixelArray = np.zeros((len(FrameToSliceInds), 
                           LabmapNda.shape[1], 
                           LabmapNda.shape[2]))
        
    # Use FrameToSliceInds to index the non-zero frames in LabmapNda:
    for i in range(len(FrameToSliceInds)):
        PixelArray[i] = LabmapNda[FrameToSliceInds[i]]
        
    
    return PixelArray, FrameToSliceInds
        




def PixArrForOneFrame(PixelArray, FrameNum):
    """
    Modify a 3D SEG pixel array so that all but one frame has non-zero 2D masks
    (i.e. set all other frames to zero).  
    
    Inputs:
        PixelArray    - (Numpy array) The SEG's PixelData loaded as a Numpy
                        array with dtype uint8
        
        FrameNum      - (Integer) The number of slices in the DICOM series
        
    Returns:
        NewPixelArray - (Numpy array) New PixelArray

    """
    
    NumOfFrames, NumOfRows, NumOfCols = PixelArray.shape
    
    # Initialise NewPixelArray:
    NewPixelArray = np.zeros((1, NumOfRows, NumOfCols), dtype='uint')
    
    NewPixelArray[0] = PixelArray[FrameNum]
            
    return NewPixelArray





def CropNonZerosIn2dMask(Orig2dMask):
    
    #print('\nShape of Orig2dMask =', Orig2dMask.shape)
    
    # The non-zero elements in the Orig2dMask:
    non0 = np.nonzero(Orig2dMask)
    
    # If there are non-zero elements:
    if non0[0].size:
    
        a = non0[0][0]
        b = non0[0][-1]
        
        c = non0[1][0]
        d = non0[1][-1]
    
        # The non-zero (cropped) array:
        Cropped2dMask = Orig2dMask[a:b, c:d]
        
        #print('\nShape of Cropped2dMask =', Cropped2dMask.shape)
        
        return Cropped2dMask
    
    else:
        #print('There were no non-zero elements in the 2D mask.')
        
        return non0[0]
    




def Compare2dMasksFrom3dMasks(OrigSegRoi, NewSegRoi, OrigDicomDir, NewDicomDir):
    """ Compare cropped masks of non-zero elements only. """ 
    
    # Get the DICOM SOP UIDs:
    OrigSOPuids = GetDicomSOPuids(DicomDir=OrigDicomDir)
    NewSOPuids = GetDicomSOPuids(DicomDir=NewDicomDir)
    
    # Get the Per-frameFunctionalGroupsSequence-to-DICOM slice indices:
    OrigPFFGStoDcmInds = GetPFFGStoDcmInds(OrigSegRoi, OrigSOPuids)
    NewPFFGStoDcmInds = GetPFFGStoDcmInds(NewSegRoi, NewSOPuids)
    
    # Combined indices from OrigPFFGStoDcmInds and NewPFFGStoDcmInds:
    AllInds = OrigPFFGStoDcmInds + NewPFFGStoDcmInds
    
    # Remove duplicates:
    AllInds = list(set(AllInds))
    
    
    # Get the 3D SEG masks:
    Orig3dMask = OrigSegRoi.pixel_array
    New3dMask = NewSegRoi.pixel_array
    
    OrigShape = Orig3dMask.shape
    NewShape = New3dMask.shape
    
    print(f'Segments exist in OrigSegRoi on slices {OrigPFFGStoDcmInds}')
    print(f'Shape of Orig3dMask = {OrigShape}')
    print(f'\nSegments exist in NewSegRoi on slices {NewPFFGStoDcmInds}')
    print(f'Shape of New3dMask = {NewShape}\n')
    
    # Initialise the 3D cropped SEG masks:
    Orig3dMaskCropped = []
    New3dMaskCropped = []
    
    for i in range(OrigShape[0]):
        cropped = CropNonZerosIn2dMask(Orig3dMask[i])
        
        Orig3dMaskCropped.append(cropped)
        
    for i in range(NewShape[0]):
        cropped = CropNonZerosIn2dMask(New3dMask[i])
        
        New3dMaskCropped.append(cropped)
    
    
    Nrows = len(AllInds)
    
    Ncols = 2
    
    n = 1 # initialised sub-plot number
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows))
    
    for i in range(len(AllInds)):
        SliceNum = AllInds[i]
        
        # Does slice SliceNum have a segment in OrigSegRoi or NewSegRoi?
        if SliceNum in OrigPFFGStoDcmInds:
            OrigFrameNum = OrigPFFGStoDcmInds.index(SliceNum)
        
            ax = plt.subplot(Nrows, Ncols, n, aspect='equal')
            ax.imshow(Orig3dMaskCropped[OrigFrameNum])
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Orig slice {SliceNum}')
            
        n += 1 # increment sub-plot number
        
        if SliceNum in NewPFFGStoDcmInds:
            NewFrameNum = NewPFFGStoDcmInds.index(SliceNum)
    
            ax = plt.subplot(Nrows, Ncols, n, aspect='equal')
            ax.imshow(New3dMaskCropped[NewFrameNum])
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'New slice {SliceNum}')
            
        n += 1 # increment sub-plot number
    
    return






def ModifySegTagVals_OLD(SegRoi, NewSegRoi, FromSliceNum, ToSliceNum, Dicoms,
                     OrigPFFGStoDcmInds, NewPFFGStoDcmInds, LogToConsole):
    """
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
        SegRoi             - Original ROI Object from a SEG file
        
        NewSegRoi          - New SEG ROI Object
        
        FromSliceNum       - (Integer) Slice index in the Source DICOM stack
                             corresponding to the segmentation to be copied 
                             (counting from 0)
                          
        ToSliceNum         - (Integer) Slice index in the Target DICOM stack 
                             where the segmentation will be copied to (counting 
                             from 0)
                          
        Dicoms             - A list of DICOM Objects that include the DICOMs 
                             that SegRoi relate to
                            
        OrigPFFGStoDcmInds - (List of integers) DICOM slice numbers that 
                             correspond to each Per-frame Functional Groups  
                             Sequence in SegRoi
                          
        NewPFFGStoDcmInds  - (List of integers) DICOM slice numbers that 
                             correspond to each Per-frame Functional Groups 
                             Sequence in NewSegRoi
                             
        LogToConsole       - (Boolean) Denotes whether or not some intermediate
                             results will be logged to the console
        

        *** I don't need the following but kept here for now:    
              
        OrigRIStoDcmInds   - (List of integers) DICOM slice numbers that 
                             correspond to each Referenced Instance Sequence in 
                             SegRoi
                            
        NewRIStoDcmInds    - (List of integers) DICOM slice numbers that 
                             correspond to each Referenced Instance Sequence in 
                             NewSegRoi
        ***
        
    Returns:
        NewSegRoi         - Modified new SEG ROI Object
    """
    
    # Get the sequence number that corresponds to the slice number to be 
    # copied, and the slice number to be copied to:
    FromOrigSeqNum = OrigPFFGStoDcmInds.index(FromSliceNum)
        
    ToNewSeqNum = NewPFFGStoDcmInds.index(ToSliceNum)
    
    if LogToConsole:
        print(f'\nFromOrigSeqNum = {FromOrigSeqNum} -->',
              f'OrigPFFGStoDcmInds[{FromOrigSeqNum}] =',
              f'{OrigPFFGStoDcmInds[FromOrigSeqNum]}')
        
        print(f'ToNewSeqNum    = {ToNewSeqNum} -->',
              f'NewPFFGStoDcmInds[{ToNewSeqNum}] =',
              f'{NewPFFGStoDcmInds[ToNewSeqNum]}')
    
    
    # Modify the Number of Frames in NewSegRoi:
    #NewSegRoi.NumberOfFrames = str(int(SegRoi.NumberOfFrames) + 1)
    
    # The 2D array to be copied is:
    #MaskToCopy = SegRoi.pixel_array[FromSliceNum]
    
    # The original 3D mask is:
    #Orig3dMask = NewSegRoi.pixel_array
    Orig3dMask = SegRoi.pixel_array
    
    # Let the new 3D mask be the original 3D mask with the last 2D mask 
    # appended to it:
    #New3dMask = np.append(Orig3dMask, Orig3dMask[-1], axis=0)
    
    # The shape of Orig3dMask:
    OrigShape = Orig3dMask.shape
    
    print('\nThe shape of the original 3D mask is', OrigShape)
    
    """
    # Re-shape New3DMask from (NumOfFrames, Rows, Columns) to 
    # (Rows, Columns, NumOfFrames):
    New3dMask = np.reshape(Orig3dMask, (OrigShape[1], OrigShape[2], OrigShape[0]))
    
    # Stack the last 2D array to the 3D array:
    New3dMask = np.dstack((New3dMask, New3dMask[:,:,-1]))
    
    # The shape of New3dMask:
    NewShape = New3dMask.shape
    
    # Reshape New3dMask to the original form:
    New3dMask = np.reshape(New3dMask, (NewShape[2], NewShape[0], NewShape[1]))
    """
    
    OrigNframes = OrigShape[0]
    OrigNrows = OrigShape[1]
    OrigNcols = OrigShape[2]
    
    NewNframes = OrigNframes + 1
    
    # Modify the Number of Frames in NewSegRoi:
    NewSegRoi.NumberOfFrames = str(NewNframes)
    
    # Create an empty 3D array with required shape:
    #New3dMask = np.zeros((NewNframes, OrigNrows, OrigNcols), dtype=int)
    New3dMask = np.zeros((NewNframes, OrigNrows, OrigNcols), dtype=bool)
    
    
    # The shape of New3dMask:
    #NewShape = New3dMask.shape
    #
    #print('\nThe shape of the new 3D mask is', NewShape)
    #
    # Convert New3dMask to bytes:
    #NewSegRoi.PixelArray = ConvertNumpyArrayToPixelArray(New3dMask) 
    #""" Note: The function ConvertNumpyArrayToPixelArray doesn't yet exist! """
    #
    #NewSegRoi.PixelData = New3dMask.tobytes()
    
    
    # Loop through each frame/sequence in NewSegRoi, modifying the relevant tag 
    # values and New3dMask:
    for i in range(NewNframes):
        # The DICOM slice index for the i^th sequence:
        s = NewPFFGStoDcmInds[i]
        
        if LogToConsole:
            print(f'\nNewPFFGStoDcmInds[{i}] = {NewPFFGStoDcmInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Class UID and Referenced SOP Instance UID 
        # so that they are consistent with the SOP Class UID and SOP Instance  
        # UID of their corresponding DICOMs:
        """ SOPClassUIDs will need to be modified for cases where the contour
        is being copied across different modality images. Not required for
        copying within the same series:
        NewSegRoi.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewSegRoi.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        """ Ditto:
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
                 
        
        # Modify the tags in FrameContentSequence:
        #NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
        #         .FrameContentSequence[0]\
        #         .StackID = f'{i + 1}'
                 
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .InStackPositionNumber = i + 1
                 
        
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .DimensionIndexValues = [1, 1, i + 1]
                     
        # Modify the ImagePositionPatient:  
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .PlanePositionSequence[0]\
                 .ImagePositionPatient = copy.deepcopy(Dicoms[s].ImagePositionPatient)
        
        
        # Modify New3dMask:         
        if i == ToNewSeqNum:
            if LogToConsole:
                print('\nThis is a new sequence.')
                print(f'Copying the {FromOrigSeqNum}th sequence from Source',
                      f'to the {i}th sequence in Target.')
                
            # This is the new sequence. Set the i^th 2D mask in New3DMask 
            # equal to the FromSeqNum^th 2D mask in Orig3dMask:
            New3dMask[i] = Orig3dMask[FromOrigSeqNum]
            #New3dMask[i,:,:] = Orig3dMask[FromOrigSeqNum,:,:]
        
        else:
            # This is an existing sequence. Find the index of 
            # OrigPFFGStoDcmInds for this sequence number (of 
            # NewPFFGStoDcmInds):
            ind = OrigPFFGStoDcmInds.index(NewPFFGStoDcmInds[i])
            
            if LogToConsole:
                print('\nThis is an existing sequence.')
                print(f'Copying the {ind}th sequence from Source to the {i}th',
                      f'sequence in Target.')
                                          
            # Set the i^th 2D mask in New3DMask equal to the ind^th 2D mask in 
            # Orig3dMask:
            New3dMask[i] = Orig3dMask[ind]
            #New3dMask[i,:,:] = Orig3dMask[ind,:,:]
        
    
    # The shape of New3dMask:
    NewShape = New3dMask.shape
    
    print('\nThe shape of the new 3D mask is', NewShape)
    
    # Ravel (to 1D array) New3dMask and pack bits 
    # (see https://github.com/pydicom/pydicom/issues/1230):
    packed = pack_bits(New3dMask.ravel())
    
    # Convert New3dMask to bytes:
    #NewSegRoi.PixelData = New3dMask.tobytes() <-- doesn't work
    NewSegRoi.PixelData = packed + b'\x00' if len(packed) % 2 else packed
                                          
    return NewSegRoi





def ModifySegTagVals_OLD2(SrcSeg, SrcDicoms, SrcPFFGStoDcmInds, FromSliceNum, 
                     TrgSeg, TrgDicoms, TrgPFFGStoDcmInds, ToSliceNum, 
                     NewTrgSeg, NewTrgPFFGStoDcmInds, LogToConsole=False):
    """
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
        SrcSeg               - (Pydicom Object) Source SEG ROI Object
        
        SrcDicoms            - (List of Pydicom Objects) List of DICOM Objects
                               that include the DICOMs that SrcSeg relate to
        
        SrcPFFGStoDcmInds    - (List of integers) The DICOM slice numbers that 
                               correspond to each Per-frame Functional Groups  
                               Sequence in SrcSeg
        
        FromSliceNum         - (Integer) Slice index in the Source DICOM stack
                               corresponding to the segmentation to be copied 
                               (counting from 0)
                               
        TrgSeg               - (Pydicom Object) Original Target SEG ROI Object
        
        TrgDicoms            - (List of Pydicom Objects) List of DICOM Objects
                               that include the DICOMs that TrgSeg relate to
        
        TrgPFFGStoDcmInds    - (List of integers) DICOM slice numbers that 
                               correspond to each Per-frame Functional Groups  
                               Sequence in TrgSeg
                               
        ToSliceNum           - (Integer) Slice index in the Target DICOM stack 
                               where the segmentation will be copied to 
                               (counting from 0)
                               
        NewTrgSeg            - (Pydicom Object) New Target SEG ROI Object
                          
        NewTrgPFFGStoDcmInds - (List of integers) DICOM slice numbers that 
                               correspond to each Per-frame Functional Groups 
                               Sequence in NewTrgSeg
                             
        LogToConsole         - (Boolean) Denotes whether or not intermediate
                               results will be logged to the console
        

        *** I don't need the following but kept here for now:    
              
        OrigRIStoDcmInds     - (List of integers) DICOM slice numbers that 
                               correspond to each Referenced Instance Sequence 
                               in SegRoi
                            
        NewRIStoDcmInds      - (List of integers) DICOM slice numbers that 
                               correspond to each Referenced Instance Sequence 
                               in NewSegRoi
        ***
        
    Returns:
        NewTrgSeg            - Modified new SEG ROI Object
    """
    
    # Get the sequence number that corresponds to the slice number to be 
    # copied, and the slice number to be copied to:
    FromSrcSeqNum = SrcPFFGStoDcmInds.index(FromSliceNum)
        
    ToNewTrgSeqNum = NewTrgPFFGStoDcmInds.index(ToSliceNum)
    
    if LogToConsole:
        print(f'\nFromSrcSeqNum     = {FromSrcSeqNum} -->',
              f'ScrPFFGStoDcmInds[{FromSrcSeqNum}] =',
              f'{SrcPFFGStoDcmInds[FromSrcSeqNum]}')
        
        print(f'ToNewTrgSeqNum    = {ToNewTrgSeqNum} -->',
              f'NewTrgPFFGStoDcmInds[{ToNewTrgSeqNum}] =',
              f'{NewTrgPFFGStoDcmInds[ToNewTrgSeqNum]}')
    
    
    # Modify the Number of Frames in NewSegRoi:
    #NewSegRoi.NumberOfFrames = str(int(SegRoi.NumberOfFrames) + 1)
    
    # The 2D array to be copied is:
    #MaskToCopy = SegRoi.pixel_array[FromSliceNum]
    
    # The original 3D mask is:
    Src3dMask = SrcSeg.pixel_array
    Trg3dMask = TrgSeg.pixel_array
    
    # Let the new 3D mask be the original 3D mask with the last 2D mask 
    # appended to it:
    #New3dMask = np.append(Orig3dMask, Orig3dMask[-1], axis=0)
    
    # The shape of Trg3dMask:
    OrigShape = Trg3dMask.shape
    
    print('\nThe shape of the original Target 3D mask is', OrigShape)
    
    OrigNframes = OrigShape[0]
    OrigNrows = OrigShape[1]
    OrigNcols = OrigShape[2]
    
    NewNframes = OrigNframes + 1
    
    # Modify the Number of Frames in NewTrgSeg:
    NewTrgSeg.NumberOfFrames = str(NewNframes)
    
    # Create an empty 3D array with required shape:
    NewTrg3dMask = np.zeros((NewNframes, OrigNrows, OrigNcols), dtype=int) # 05/11/20
    #NewTrg3dMask = np.zeros((NewNframes, OrigNrows, OrigNcols), dtype=bool) # 05/11/20
    
    
    # Loop through each frame/sequence in NewTrgSeg, modifying the relevant tag 
    # values and NewTrg3dMask:
    for i in range(NewNframes):
        # The DICOM slice index for the i^th sequence:
        s = NewTrgPFFGStoDcmInds[i]
        
        if LogToConsole:
            print(f'\nNewTrgPFFGStoDcmInds[{i}] = {NewTrgPFFGStoDcmInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Class UID and Referenced SOP Instance UID 
        # so that they are consistent with the SOP Class UID and SOP Instance  
        # UID of their corresponding DICOMs:
        """ SOPClassUIDs will need to be modified for cases where the contour
        is being copied across different modality images. Not required for
        copying within the same series:
        NewSegRoi.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
        
        """ Ditto:
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
                 
        
        # Modify the tags in FrameContentSequence:
        #NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
        #         .FrameContentSequence[0]\
        #         .StackID = f'{i + 1}'
            
        """ 
        I can't find InStackPositionNumber in a SEG file created in the
        OHIF-Viewer.  It seems that this tag is generated from RTS-to-SEG 
        conversions using James' Java tool. 
        Commenting this out for now.
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .InStackPositionNumber = i + 1
        """
                 
        
        """ This is wrong - the DIVs need to be indexed against the items in 
        the Referenced Instance Sequence (and the first dimension needs to 
        relate to the ROI number as there can be more than one)!
        """
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .DimensionIndexValues = [1, 1, i + 1]
                     
        # Modify the ImagePositionPatient:  
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .PlanePositionSequence[0]\
                 .ImagePositionPatient = copy.deepcopy(TrgDicoms[s].ImagePositionPatient)
        
        
        # Modify NewTrg3dMask:         
        if i == ToNewTrgSeqNum:
            if LogToConsole:
                print('\nThis is a new sequence.')
                print(f'Copying the {FromSrcSeqNum}th sequence from Source',
                      f'to the {i}th sequence in Target.')
                
            # This is the new sequence. Set the i^th 2D mask in NewTrg3DMask 
            # equal to the FromSrcSeqNum^th 2D mask in Src3dMask:
            NewTrg3dMask[i] = Src3dMask[FromSrcSeqNum]
            #New3dMask[i,:,:] = Orig3dMask[FromOrigSeqNum,:,:]
        
        else:
            # This is an existing sequence. Find the index of 
            # TrgPFFGStoDcmInds for this sequence number (of 
            # NewTrgPFFGStoDcmInds):
            ind = TrgPFFGStoDcmInds.index(NewTrgPFFGStoDcmInds[i])
            
            if LogToConsole:
                print('\nThis is an existing sequence.')
                print(f'Copying the {ind}th sequence from Target to the {i}th',
                      f'sequence in Target.')
                                          
            # Set the i^th 2D mask in NewTrg3DMask equal to the ind^th 2D mask
            # in Trg3dMask:
            NewTrg3dMask[i] = Trg3dMask[ind]
            #New3dMask[i,:,:] = Orig3dMask[ind,:,:]
        
    
    # The shape of NewTrg3dMask:
    NewShape = NewTrg3dMask.shape
    
    print('\nThe shape of the new 3D mask is', NewShape)
    
    # Ravel (to 1D array) NewTrg3dMask and pack bits 
    # (see https://github.com/pydicom/pydicom/issues/1230):
    packed = pack_bits(NewTrg3dMask.ravel())
    
    # Convert NewTrg3dMask to bytes:
    #NewSegRoi.PixelData = New3dMask.tobytes() <-- doesn't work
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
                                          
    return NewTrgSeg






def ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, NewTrgRoiNums,  
                     NewTrgPFFGStoDcmInds, FromRoiNum, FromSliceNum, SrcSeg,
                     ToRoiNum, LogToConsole=False):
    """
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:                              
        TrgSeg               - (Pydicom Object) Original Target SEG ROI Object
        
        TrgDicoms            - (List of Pydicom Objects) List of DICOM Objects
                               that include the DICOMs that TrgSeg relate to
        
        NewTrgPixArr         - (Numpy array) New Target SEG pixel array
        
        NewTrgRoiNums        - (List of integers) ROI numbers that correspond 
                               to each frame in NewTrgSeg
        
        NewTrgPFFGStoDcmInds - (List of integers) DICOM slice numbers that 
                               correspond to each frame in NewTrgSeg (which 
                               will also become the Per-frame Functional Groups 
                               Sequence to DICOM slice numbers)
                               
        FromRoiNum           - (Integer) Source ROI number of the segmentation 
                               that was copied (counting from 0) (this is only 
                               used for generating a new Series Description) 
        
        FromSliceNum         - (Integer) Source slice index in the DICOM stack
                               corresponding to the segmentation that was  
                               copied (counting from 0) (this is only used for  
                               generating a new Series Description)
                               
        SrcSeg               - (Pydicom Object) Source SEG ROI Object (this is 
                               only used for generating a new Series 
                               Description)
        
        ToRoiNum             - (Integer) Target ROI number that the 
                               segmentation was copied to (counting from 0) 
                               (this is only used for generating a new Series
                               Description) 
                             
        LogToConsole         - (Boolean) Denotes whether or not intermediate
                               results will be logged to the console
        
        ***
        
    Returns:
        NewTrgSeg            - (Pydicom Object) New Target SEG ROI Object
    """
    
    # Use TrgSeg as a template for NewTrgSeg:
    NewTrgSeg = copy.deepcopy(TrgSeg)
    
    # The original and new number of frames in pixel array:
    OrigN = TrgSeg.pixel_array.shape[0]
    NewN = NewTrgPixArr.shape[0]
    
    # Increase the length of the Per-Frame Functional Groups Sequence in 
    # NewTrgSeg:
    for i in range(NewN - OrigN):
        NewTrgSeg = AddToPFFGS(NewTrgSeg)
            
    
    # Loop through each frame/sequence in NewTrgSeg, modifying the relevant tag 
    # values and NewTrg3dMask:
    for i in range(NewN):
        # The DICOM slice index for the i^th sequence:
        s = NewTrgPFFGStoDcmInds[i]
        
        # The ROI Number (Segment number):
        r = NewTrgRoiNums[i]
        
        if LogToConsole:
            print(f'\nNewTrgPFFGStoDcmInds[{i}] = {NewTrgPFFGStoDcmInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Class UID and Referenced SOP Instance UID 
        # so that they are consistent with the SOP Class UID and SOP Instance  
        # UID of their corresponding DICOMs:
        """ SOPClassUIDs will need to be modified for cases where the contour
        is being copied across different modality images. Not required for
        copying within the same series:
        NewSegRoi.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
        
        """ Ditto:
        NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
                 
        
        # Modify the tags in FrameContentSequence:
        #NewSegRoi.PerFrameFunctionalGroupsSequence[i]\
        #         .FrameContentSequence[0]\
        #         .StackID = f'{i + 1}'
            
        """ 
        I can't find InStackPositionNumber in a SEG file created in the
        OHIF-Viewer.  It seems that this tag is generated from RTS-to-SEG 
        conversions using James' Java tool. 
        Commenting this out for now.
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .InStackPositionNumber = i + 1
        """
        
        # Modify the Dimension Index Values using the RoiNum and SliceNum:
        """
        Note: r and s are integers starting from 0.  Need to add 1 to them so 
        that the integers in DIVs start from 1.
        """        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .DimensionIndexValues = [r + 1, s + 1]
                     
        # Modify the ImagePositionPatient: 
        IPP = copy.deepcopy(TrgDicoms[s].ImagePositionPatient)
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .PlanePositionSequence[0]\
                 .ImagePositionPatient = IPP
                 
        # Modify the Referenced Segment Number (= RoiNum + 1):
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .SegmentIdentificationSequence[0]\
                 .ReferencedSegmentNumber = r + 1
                 
        
        print(f'\nFrame number = {i}')
        print(f'DICOM slice number = {s}')
        print(f'Roi/Segment number = {r}')
        
        
    
    # The shape of NewTrgPixArr:
    NewShape = NewTrgPixArr.shape
    
    print('\nThe shape of the new 3D pixel array is', NewShape)
    
    # Ravel (to 1D array) NewTrgPixArr and pack bits 
    # (see https://github.com/pydicom/pydicom/issues/1230):
    packed = pack_bits(NewTrgPixArr.ravel())
    
    # Convert NewTrgPixArr to bytes:
    #NewSegRoi.PixelData = NewTrgPixArr.tobytes() <-- doesn't work
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    # Modify the Number of Frames:
    NewTrgSeg.NumberOfFrames = f"{len(NewTrgPFFGStoDcmInds)}"
    
    # Generate a new SOP Instance UID:
    NewTrgSeg.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgSeg.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    # Modify the Series Description:
    NewTrgSegSD = 'RPCopy_of_R' + str(FromRoiNum) + '_s' + str(FromSliceNum) \
                  + '_in_' + SrcSeg.SeriesDescription + '_to_R' + str(ToRoiNum)\
                  + '_in_' + TrgSeg.SeriesDescription
                  
    #NewTrgSegSD = 'Case3d_RPCopy'
    
    NewTrgSeg.SeriesDescription = NewTrgSegSD
                                          
    return NewTrgSeg






def ExportSegRoi(NewTrgSeg, SrcSegFpath, TrgSegFpath, NamePrefix, FromSliceNum,
                 ExportDir):
    """
    Export contour ROI Object to disc.
    
    
    Inputs:
        NewTrgSeg      - New Target SEG ROI object to be exported
        
        SrcSegFpath    - (String) Full path of the Source DICOM SEG file
        
        TrgSegFpath    - (String) Full path of the Target DICOM SEG file
        
        NamePrefix     - (String) Prefix to be added to the assigned filename
                         after the DateTime stamp (e.g. 'Case3d_RPCopy') 
                            
        ExportDir      - (String) Directory where the SEG is to be exported
                           
    
                            
    Returns:
        NewTrgSegFpath - (String) Full path of the exported Target DICOM SEG 
                         file
    
    """
    
    # Modify the Series Description to a much shorter one (the one assigned in
    # ModifySegTagVals is very long):    
    NewTrgSeg.SeriesDescription = NamePrefix
    
    # The current DateTime:
    CDT = time.strftime("%Y%m%d_%H%M%S_", time.gmtime())
    
    # The filename of the Source and Target SEGs:
    SrcSegFname = os.path.split(SrcSegFpath)[1]
    TrgSegFname = os.path.split(TrgSegFpath)[1]
    
    NewTrgSegFname = CDT + NamePrefix + '_s' + str(FromSliceNum) + '_from_' \
                     + os.path.splitext(SrcSegFname)[0] + '_to_' + TrgSegFname
    
    NewTrgSegFpath = os.path.join(ExportDir, NewTrgSegFname)
    
    
    NewTrgSeg.save_as(NewTrgSegFpath)
        
    print('\nNew Target SEG ROI exported to:\n\n', NewTrgSegFpath)
    
    return NewTrgSegFpath

