# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:31:16 2020

@author: ctorti
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







def GetFrameFromPixArr_OLD(PixArr, FrameNum):
    """
    Extract a single (2D) frame from a 3D pixel array.  
    
    Inputs:
    ------
    
    PixArr : Numpy array
        PixelData from a SEG file loaded as a Numpy array
        
    FrameNum : integer
        The number of slices in the DICOM series.
        
    
    Outputs:
    -------
    
    Frame : Numpy array
        Frame from PixArr.

    """
    
    import numpy as np
    
    F, R, C = PixArr.shape
    
    # Initialise Frame:
    Frame = np.zeros((1, R, C), dtype='uint')
    
    Frame[0] = PixArr[FrameNum]
            
    return Frame





def GetDIVsGroupedByRoi_OLD1(SegRoi):
    """
    Get the DimensionIndexValues in a SEG ROI grouped by common ROI.
    
    Inputs:
        SegRoi - ROI Object from a SEG file
        
        
    Returns:
        DIVs   - (List of lists of integers) Dimension Index Values grouped by
                 ROI
        
        
    Note:
        DIVs is a nested list of lists of the form:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ],
          [ [2, 6], [2, 7], [2, 8], ... ],
          ... 
         ]
          
         
        where each list of index values that belong to a common ROI are 
        separated into separate lists.  If there's only one ROI, the list of
        index values will still be nested in a list to maintain continuity:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ]
         ]
    """
    
    # Get the flat list of DIVs:
    FlatList = GetDIVs(SegRoi) 
        
    # Group the DIVs by common ROI.
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in FlatList:
        RoiNums.append(DIV[0])
        
        
    UniqueRoiNums = list(set(RoiNums))
    
    
    if len(UniqueRoiNums) == 1:
        DimIndVals = [FlatList]
        
    else:
        DimIndVals = []
        
        RoiNum = 1 # initial value
        
        ListThisRoi = [] # initialise
        
        for i in range(len(FlatList)):
            
            if FlatList[i][0] == RoiNum:
                ListThisRoi.append(FlatList[i])
                
            else:
                DimIndVals.append(ListThisRoi)
                
                RoiNum += 1 # increment
                
                ListThisRoi = [] # re-initialise
                
                ListThisRoi.append(FlatList[i])
        
            
        DimIndVals.append(ListThisRoi)
             
    return DimIndVals







def GetDIVsGroupedByRoi_OLD2(SegRoi):
    """
    Get the DimensionIndexValues in a SEG ROI grouped by common ROI.
    
    Inputs:
        SegRoi - ROI Object from a SEG file
        
        
    Returns:
        DIVs   - (List of lists of integers) Dimension Index Values grouped by
                 ROI
        
        
    Note:
        DIVs is a nested list of lists of the form:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ],
          [ [2, 6], [2, 7], [2, 8], ... ],
          ... 
         ]
          
         
        where each list of index values that belong to a common ROI are 
        separated into separate lists.  If there's only one ROI, the list of
        index values will still be nested in a list to maintain continuity:
            
        [ [ [1, 10], [1, 11], [1, 12], ... ]
         ]
    """
    
    # Get the flat list of DIVs:
    FlatList = GetDIVs(SegRoi) 
        
    # Group the DIVs by common ROI.
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in FlatList:
        RoiNums.append(DIV[0])
        
        
    UniqueRoiNums = list(set(RoiNums))
    
    
    if len(UniqueRoiNums) == 1:
        DimIndVals = [FlatList]
        
    else:
        DimIndVals = []
        
        for n in range(len(UniqueRoiNums)):
            ListThisRoi = [] # initialise
            
            for i in range(len(FlatList)):
                if FlatList[i][0] == n:
                    ListThisRoi.append(FlatList[i])
        
            DimIndVals.append(ListThisRoi)
             
    return DimIndVals






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
    NewSegRoi = deepcopy(SegRoi)
        
    # The last item in Referenced Instance Sequence:
    last = deepcopy(SegRoi.ReferencedSeriesSequence[0]\
                          .ReferencedInstanceSequence[-1])
    
    # Append to the Referenced Instance Sequence:
    NewSegRoi.ReferencedSeriesSequence[0]\
             .ReferencedInstanceSequence\
             .append(last)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = deepcopy(SegRoi.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    NewSegRoi.PerFrameFunctionalGroupsSequence\
             .append(last)
             
    return NewSegRoi






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
                 .ImagePositionPatient = deepcopy(Dicoms[s].ImagePositionPatient)
        
        
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
                 .ImagePositionPatient = deepcopy(TrgDicoms[s].ImagePositionPatient)
        
        
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







def ModifySegTagVals_OLD3(TrgSeg, TrgDicoms, NewTrgPixArr, NewTrgRoiNums,  
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
    NewTrgSeg = deepcopy(TrgSeg)
    
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
        IPP = deepcopy(TrgDicoms[s].ImagePositionPatient)
        
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








def ExportSegRoi(TrgSegRoi, SrcSegFpath, NamePrefix, ExportDir):
    """
    Export contour ROI Object to disc.
    
    
    Inputs:
        TrgSegRoi   - Target SEG ROI object to be exported
        
        SrcSegFpath - (String) Full path of the Source DICOM SEG file (used to
                      generate the filename of the new SEG file)
        
        NamePrefix  - (String) Prefix to be added to the assigned filename and
                      Content Label (e.g. 'MR4_S9_s23_to_MR4_S9_s22') 
                            
        ExportDir   - (String) Directory where the SEG is to be exported
                           
    
                            
    Returns:
        TrgSegFpath - (String) Full path of the exported Target DICOM SEG file
    
    """
    
    """
    The following was moved into the main function ModifySegTagVals:
    # Generate a new SOP Instance UID:
    SegRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    SegRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    """
    
    # Get the filename of the original SEG file:
    SrcSegFname = os.path.split(SrcSegFpath)[1]
    
    # Modify the Content Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    #TrgSegRoi.ContentDate = NewDate
    #TrgSegRoi.ContentTime = NewTime
    
    FnamePrefix = NewDate + '_' + NewTime + '_' + NamePrefix + '_from_' 
    
    ContentLabelPrefix = NamePrefix + '_from_'
    
    # Modify the Content Label (this appears under Name in XNAT):
    #TrgSegRoi.ContentLabel = 'Copy_of_' + TrgSegRoi.ContentLabel
    TrgSegRoi.ContentLabel = ContentLabelPrefix + TrgSegRoi.ContentLabel
    
    # Modify Content Description:
    TrgSegRoi.ContentDescription = ContentLabelPrefix \
                                   + TrgSegRoi.ContentDescription
                                   
    # Modify the Series Description as with the Content Label:
    TrgSegRoi.SeriesDescription = ContentLabelPrefix \
                                  + TrgSegRoi.SeriesDescription
    
    
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = 'New_' + SrcSegFname + '_' + NewDate + '_' + NewTime
    #TrgRoiFname = FnamePrefix + '_SEG_' + NewDate + '_' + NewTime + '_' \
    #              + SrcSegFname
    TrgSegFname = FnamePrefix + SrcSegFname
    
    TrgSegFpath = os.path.join(ExportDir, TrgSegFname)
    
    TrgSegRoi.save_as(TrgSegFpath)
        
    print('\nSEG ROI exported to:\n\n', TrgSegFpath)
    
    return TrgSegFpath







def Compare2dMasksFrom3dMasksOLD(OrigSegRoi, NewSegRoi, OrigSliceNum, NewSliceNum,
                              OrigDicomDir, NewDicomDir):
    """ Compare cropped masks of non-zero elements only. """ 
    
    # Get the DICOM SOP UIDs:
    OrigSOPuids = GetDicomSOPuids(DicomDir=OrigDicomDir)
    NewSOPuids = GetDicomSOPuids(DicomDir=NewDicomDir)
    
    # Get the Per-frameFunctionalGroupsSequence-to-DICOM slice indices:
    OrigPFFGStoDcmInds = GetPFFGStoDcmInds(OrigSegRoi, OrigSOPuids)
    NewPFFGStoDcmInds = GetPFFGStoDcmInds(NewSegRoi, NewSOPuids)
    
    # Get the 3D SEG masks:
    Orig3dMask = OrigSegRoi.pixel_array
    New3dMask = NewSegRoi.pixel_array
    
    OrigShape = Orig3dMask.shape
    NewShape = New3dMask.shape
    
    print(f'Segments exist in OrigSegRoi on slices {OrigPFFGStoDcmInds}')
    print(f'Shape of Orig3dMask = {OrigShape}')
    print(f'\nSegments exist in NewSegRoi on slices {NewPFFGStoDcmInds}')
    print(f'Shape of New3dMask = {NewShape}\n')
    
    
    # Initialise Orig2dMaskCropped and New2dMaskCropped:
    Orig2dMaskCropped = np.zeros((1, 1))
    New2dMaskCropped = np.zeros((1, 1))
    
    
    # Does slice OrigSliceNum in OrigSegRoi have a segment?
    if OrigSliceNum in OrigPFFGStoDcmInds:
        # A segment frame exists for OrigSliceNum:
        OrigFrame = True
        
        # Get the frame number for this slice number:
        OrigFrameNum = OrigPFFGStoDcmInds.index(OrigSliceNum)
        
        Orig2dMask = OrigSegRoi.pixel_array[OrigFrameNum]
        
        print(f'Shape of Orig2dMask = Orig3dMask[{OrigFrameNum}] = {Orig2dMask.shape}')
    
        #print('Cropping Orig2dMask to array of non-zero elements..')
    
        Orig2dMaskCropped = CropNonZerosIn2dMask(Orig2dMask)
        #Orig2dMaskCropped = CropNonZerosIn2dMask(Orig3dMask[OrigFrameNum])
        
    else:
        OrigFrame = False
        
        print(f'A segment does not exist on slice {OrigSliceNum} in OrigSegRoi.')
        
    
    # Does slice NewSliceNum in NewSegRoi have a segment?
    if NewSliceNum in NewPFFGStoDcmInds:
        # A segment frame exists for NewSliceNum:
        NewFrame = True
        
        # Get the frame number for this slice number:
        NewFrameNum = NewPFFGStoDcmInds.index(NewSliceNum)
        
        New2dMask = NewSegRoi.pixel_array[NewFrameNum]
    
        print(f'Shape of New2dMask = New3dMask[{NewFrameNum}]  = {New2dMask.shape}')
    
        #print('Cropping New2dMask to array of non-zero elements..')
        
        New2dMaskCropped = CropNonZerosIn2dMask(New2dMask)
        #New2dMaskCropped = CropNonZerosIn2dMask(New3dMask[FrameNum])
        
        
    else:
        NewFrame = False
        
        print(f'A segment does not exist on slice {NewSliceNum} in NewSegRoi.')
    
    
    # Proceed only if either cropped arrays are not empty (i.e. if at least one
    # of the 2D masks had non-zero elements):
    if Orig2dMaskCropped.any() or New2dMaskCropped.any():
    
        print('')
        
        fig, ax = plt.subplots(1, 2)
        
        if Orig2dMaskCropped.any():
            print('Shape of Orig2dMaskCropped =', Orig2dMaskCropped.shape)
            
            ax = plt.subplot(1, 2, 1, aspect='equal')
            ax.imshow(Orig2dMaskCropped)
        
        else:
            if OrigFrame:
                print(f'Frame {OrigFrameNum} in OrigSegRoi does not have any',
                      'non-zero elements.')
        
        if New2dMaskCropped.any():
            print('Shape of New2dMaskCropped  =', New2dMaskCropped.shape)
        
            ax = plt.subplot(1, 2, 2, aspect='equal')
            ax.imshow(New2dMaskCropped)
        
        else:
            if NewFrame:
                print(f'Frame {NewFrameNum} in NewSegRoi does not have any',
                      'non-zero elements.')
            
    else:
        if OrigFrame and NewFrame:
            print(f'Frame {OrigFrameNum} and {NewFrameNum} does not have any',
                  'non-zero elements in OrigSegRoi and NewSegRoi.')
    
    return







def DELETE_THIS_AddToRIStoDcmInds(RIStoDcmInds, ToSliceNum):
    """
    Append the DICOM slice number (ToSliceNum) representing the segmentation to 
    be added to the Referenced Instance Sequence to RIStoDcmInds (the 
    ReferencedInstanceSequence-to-DICOM slice indeces), and sort the list of 
    indeces.
    
    Inputs:
        RIStoDcmInds    - List of integers of the DICOM slice numbers that 
                          correspond to each Referenced Instance Sequence
        
        ToSliceNum      - (Integer) Slice index in the Target DICOM stack where
                          the segmentation will be copied to (counting from 0)
        
    Returns:
        NewRIStoDcmInds - Modified and sorted list of RIStoDcmInds with the  
                          slice index to be copied
    """
    
    NewRIStoDcmInds = copy.deepcopy(RIStoDcmInds)

    NewRIStoDcmInds.append(ToSliceNum)
    
    NewRIStoDcmInds.sort()
    
    return NewRIStoDcmInds




def DELETE_THIS_AddToPFFGStoDcmInds(PFFGStoDcmInds, ToSliceNum):
    """
    Append the DICOM slice number (ToSliceNum) representing the segmentation to 
    be added to the Per-frame Functional Groups Sequence to PFFGStoDcmInds (the 
    PerFrameFunctionalGroupsSequence-to-DICOM slice indeces), and sort the list 
    of indeces.
    
    Inputs:
        PFFGStoDcmInds    - List of integers of the DICOM slice numbers that 
                            correspond to each Per-frame Functional Groups 
                            Sequence
        
        ToSliceNum        - (Integer) Slice index in the Target DICOM stack 
                            where the segmentation will be copied to (counting 
                            from 0)
        
    Returns:
        NewPFFGStoDcmInds - Modified and sorted list of PFFGStoDcmInds with the 
                            slice index to be copied
    """
    
    NewPFFGStoDcmInds = copy.deepcopy(PFFGStoDcmInds)

    NewPFFGStoDcmInds.append(ToSliceNum)
    
    NewPFFGStoDcmInds.sort()
    
    return NewPFFGStoDcmInds







def ModifySegTagVals_OLD(TrgSeg, TrgDicoms, NewTrgPixArr, NewTrgSegNums,  
                     NewTrgPFFGStoSliceInds, FromSegNum, FromSliceNum, SrcSeg,
                     LogToConsole=False):
    """
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
    ------
                              
    TrgSeg : Pydicom object
        Original Target SEG ROI Object.
    
    TrgDicoms : List of Pydicom objects
        List of DICOM Objects that include the DICOMs that TrgSeg relate to.
    
    NewTrgPixArr : Numpy array
        New Target SEG pixel array.
    
    NewTrgSegNums : List of integers
        Segment numbers that correspond to each frame in NewTrgSeg.
    
    NewTrgPFFGStoSliceInds : List of integers
        Slice numbers that correspond to each frame in NewTrgSeg.
                           
    FromSegNum : integer
        Source segment number of the segmentation that was copied (counting  
        from 0) (this is only used for generating a new Series Description).
    
    FromSliceNum : integer
        Source slice index in the DICOM stack corresponding to the segmentation 
        that was copied (counting from 0) (this is only used for generating a 
        new Series Description).
                           
    SrcSeg : Pydicom object
        Source SEG ROI object (this is only used for generating a new Series 
        Description).
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    NewTrgSeg : Pydicom object
        Modified Target SEG ROI object.
    """
    
    from copy import deepcopy
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    from pydicom.uid import generate_uid
    
    # Use TrgSeg as a template for NewTrgSeg:
    NewTrgSeg = deepcopy(TrgSeg)
    
    # The original and new number of frames in pixel array:
    OrigF = TrgSeg.pixel_array.shape[0]
    NewF = NewTrgPixArr.shape[0]
    
    # Increase the length of the Per-Frame Functional Groups Sequence in 
    # NewTrgSeg:
    for i in range(NewF - OrigF):
        NewTrgSeg = AddToPFFGS(NewTrgSeg)
            
    
    # Loop through each frame/sequence in NewTrgSeg, modifying the relevant tag 
    # values:
    for i in range(NewF):
        # The slice index for the i^th sequence:
        s = NewTrgPFFGStoSliceInds[i]
        
        # The new Target segment number:
        r = NewTrgSegNums[i]
        
        if LogToConsole:
            print(f'\nNewTrgPFFGStoSliceInds[{i}] = {NewTrgPFFGStoSliceInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Instance UID:
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
        
        # Modify the Dimension Index Values using the RoiNum and SliceNum:
        """
        Note: r and s are integers starting from 0.  Need to add 1 to them so 
        that the integers in DIVs start from 1.
        """        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .DimensionIndexValues = [r + 1, s + 1]
                     
        # Modify the ImagePositionPatient: 
        IPP = deepcopy(TrgDicoms[s].ImagePositionPatient)
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .PlanePositionSequence[0]\
                 .ImagePositionPatient = IPP
                 
        # Modify the Referenced Segment Number (= RoiNum + 1):
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .SegmentIdentificationSequence[0]\
                 .ReferencedSegmentNumber = r + 1
                 
        
        print(f'\nFrame number       = {i}')
        print(f'DICOM slice number = {s}')
        print(f'Segment number     = {r}')
        
        
    
    
    # The original and new number of segments:
    """
    Since the copied segmentation will be copied to a new segment NewS is 
    always OrigS + 1.
    """
    OrigS = len(TrgSeg.SegmentSequence)
    NewS = OrigS + 1
    
    # Increase the length of the Segment Sequence in NewTrgSeg:
    for i in range(NewS - OrigS):
        NewTrgSeg = AddToSS(NewTrgSeg)
        
        
        
    #print(f'The new Target pixel array has shape {NewTrgPixArr.shape}')
    
    # Ravel (to 1D array) NewTrgPixArr and pack bits 
    # (see https://github.com/pydicom/pydicom/issues/1230):
    packed = pack_bits(NewTrgPixArr.ravel())
    
    # Convert NewTrgPixArr to bytes:
    #NewSegRoi.PixelData = NewTrgPixArr.tobytes() <-- doesn't work
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    # Modify the Number of Frames:
    NewTrgSeg.NumberOfFrames = f"{len(NewTrgPFFGStoSliceInds)}"
    
    # Generate a new SOP Instance UID:
    #NewTrgSeg.SOPInstanceUID = pydicom.uid.generate_uid()
    NewTrgSeg.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    #NewTrgSeg.SeriesInstanceUID = pydicom.uid.generate_uid()
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    # The segmentation was copied to a new segment with number:
    ToSegNum = NewTrgSegNums[-1] 
    """ verify the above """
    
    # Modify the Series Description:
    NewTrgSegSD = 'RPCopy_of_Seg' + str(FromSegNum) + '_s' + str(FromSliceNum) \
                  + '_in_' + SrcSeg.SeriesDescription + '_to_Seg' \
                  + str(ToSegNum) + '_in_' + TrgSeg.SeriesDescription
                  
    #NewTrgSegSD = 'Case3d_RPCopy'
    
    NewTrgSeg.SeriesDescription = NewTrgSegSD
    
    # Modify the Segment Label of the last (new) segment:
    FromSegLabel = deepcopy(NewTrgSeg.SegmentSequence[FromSegNum].SegmentLabel)
    
    NewSegLabel = 'NEW_' + FromSegLabel
    
    #print(f'\n--> FromSegNum = {FromSegNum}, ToSegNum = {ToSegNum}')
    #print(f'\n--> len(NewTrgSeg.SegmentSequence) = {len(NewTrgSeg.SegmentSequence)}')
    
    NewTrgSeg.SegmentSequence[ToSegNum].SegmentLabel = NewSegLabel
                                          
    return NewTrgSeg






def ExportSeg(TrgSeg, SrcSegFpath, NamePrefix, ExportDir):
    """
    Export SEG to disc.
    
    Inputs:
    ------
    
    TrgSeg : Pydicom object
        New Target SEG ROI object to be exported.
        
    SrcSegFpath : string
        Full path of the Source DICOM SEG file (used to generate the 
        filename of the new RTS file).
        
    NamePrefix : string
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
                            
    ExportDir : string
        Directory where the SEG is to be exported.
                           
                            
    Outputs:
    -------
    
    TrgSegFpath : string
        Full path of the exported Target DICOM SEG file.
    
    """
    
    import os
    import time
    
    # The filename of the Source SEG:
    SrcSegFname = os.path.split(SrcSegFpath)[1]
    
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    FnamePrefix = CurrentDateTime + '_' + NamePrefix + '_from_' 

    # Modify the Series Description to a much shorter one (the one assigned in
    # ModifySegTagVals is very long):    
    #NewTrgSeg.SeriesDescription = NamePrefix
    
    TrgSegFname = FnamePrefix + SrcSegFname
    
    TrgSegFpath = os.path.join(ExportDir, TrgSegFname)
    
    TrgSeg.save_as(TrgSegFpath)
        
    print('\nSEG ROI exported to:\n\n', TrgSegFpath)
    
    return TrgSegFpath






def InitialiseSeg(SegTemplate, SearchString, PFFGStoSliceInds, DicomDir, 
                  NamePrefix='', LogToConsole=False):
    """
    Initialise a SEG object based on an existing SEG object (SegTemplate).
         
    
    Inputs:
    ******
    
    SegTemplate : Pydicom object
        SEG object that will be used as a template to create the new object.

    SearchString : string
        All or part of the SegmentLabel of the segment of interest. 
        
    PFFGStoSliceInds : List of integers
        List (for each segmentation) of the slice numbers that correspond to 
        each frame in the SEG's pixel array.
                    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    NamePrefix : string (optional; '' by default)
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    Seg : Pydicom object
        New SEG object.
    """
    
    from pydicom.uid import generate_uid
    import time
    from copy import deepcopy
    from DicomTools import GetRoiNum
    from DicomTools import GetDicomFpaths
    from DicomTools import ImportDicom
    from ImageTools import GetImageAttributes
    
    Dicom = ImportDicom(DicomDir)
    
    NumOfDicoms = len(GetDicomFpaths(DicomDir))
    
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    
    SegNum = GetRoiNum(SegTemplate, SearchString)
    
    #print(f'\n\n***SegNum = {SegNum}')
    
    # The Segment Label of the segment containing the segmentation to be 
    # copied:
    SegLabel = SegTemplate.SegmentSequence[SegNum].SegmentLabel
    
    #print(f'\n\n***SegLabel = {SegLabel}')
    
    # Use SegTemplate as a template for Seg:
    Seg = deepcopy(SegTemplate)
    
    # Generate a new SOPInstance UID:
    Seg.SOPInstanceUID = generate_uid()
    
    # Generate a new SeriesInstance UID:
    Seg.SeriesInstanceUID = generate_uid()
    
    Seg.file_meta.MediaStorageSOPInstanceUID = deepcopy(Seg.SOPInstanceUID) # 18/12
    
    
    """Added on 16/12 to see if it solves the issues for UseCases 2b and 3b-i.
    It didn't resolve the issue.
    
    # Generate a new DimensionOrganizationUID:
    Seg.DimensionOrganizationSequence[0].DimensionOrganizationUID = generate_uid()
    
    # Set the DimensionOrganizationUIDs in DimensionIndexSequence to the new
    # uid:
    DOuid = deepcopy(Seg.DimensionOrganizationSequence[0].DimensionOrganizationUID)
    
    DIS = deepcopy(Seg.DimensionIndexSequence)
    
    for i in range(len(DIS)):
        Seg.DimensionIndexSequence[i].DimensionOrganizationUID = DOuid
    
    end of code added on 16/12."""
    
    # Modify the Series Description:
    #if 'SEG' in NamePrefix:
    #    SD = NamePrefix + '_from_' + Seg.SeriesDescription
    #else:
    #    SD = NamePrefix + '_SEG_from_' + Seg.SeriesDescription
    
    #print(f'\n\n***NamePrefix = {NamePrefix}')
        
    #Seg.SeriesDescription = SD
    Seg.SeriesDescription = NamePrefix
    
    #print(f'\n\n***Seg.SeriesDescription = {Seg.SeriesDescription}')
    
    # Modify the Content Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    Seg.ContentDate = NewDate
    Seg.ContentTime = NewTime
    
    
    
    # Modify various tags to match the values in Dicoms:
    Seg.StudyDate = Dicom.StudyDate
    Seg.SeriesDate = Dicom.SeriesDate
    #Seg.ContentDate = Dicom.ContentDate
    Seg.StudyTime = Dicom.StudyTime
    Seg.SeriesTime = Dicom.SeriesTime
    #Seg.ContentTime = Dicom.ContentTime
    
    Seg.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = Dicom.SeriesInstanceUID # 15/12/20
    
    Seg.PatientName = Dicom.PatientName
    Seg.PatientID = Dicom.PatientID
    Seg.PatientBirthDate = Dicom.PatientBirthDate
    Seg.PatientSex = Dicom.PatientSex
    Seg.PatientAge = Dicom.PatientAge
    Seg.StudyInstanceUID = Dicom.StudyInstanceUID
    
    #Seg.StudyID = Dicom.StudyID # Leave StudyID unchanged (22/12/20)
    
    Seg.FrameOfReferenceUID = Dicom.FrameOfReferenceUID
    Seg.NumberOfFrames = f"{len(PFFGStoSliceInds)}"
    Seg.Rows = Dicom.Rows
    Seg.Columns = Dicom.Columns
    
    Seg.SharedFunctionalGroupsSequence[0]\
       .PlaneOrientationSequence[0]\
       .ImageOrientationPatient = Dicom.ImageOrientationPatient
       
    """ 
    Notes:
        
    1. The SliceThickness appears to be the z-Spacing rather than the
    SliceThickness from the DICOM metadata. 
    
    2. The character limit of VR DS is 16 characters.
    """
    Zspacing = f"{Spacings[2]}"
    
    if len(Zspacing) > 16:
        Zspacing = Zspacing[:16]
    
    Seg.SharedFunctionalGroupsSequence[0]\
       .PixelMeasuresSequence[0]\
       .SliceThickness = Zspacing
    
    Seg.SharedFunctionalGroupsSequence[0]\
       .PixelMeasuresSequence[0]\
       .SpacingBetweenSlices = Zspacing

    Seg.SharedFunctionalGroupsSequence[0]\
       .PixelMeasuresSequence[0]\
       .PixelSpacing = Dicom.PixelSpacing
       
       
    """ 
    Modify the number of sequences in ReferencedInstanceSequence (RIS). 
    The number of sequences RIS should be equal to the number of DICOMs.
    """
    # The original and required number of sequences in 
    # ReferencedInstanceSequence:
    OrigRIS = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)

    ReqRIS = NumOfDicoms
    
    # Add/remove sequences by appending/popping the last sequence N times, 
    # where N is the difference between ReqRIS and OrigRIS:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            Seg.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence\
               .append(deepcopy(Seg.ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1]))          
    else:
        for i in range(OrigRIS - ReqRIS):
            Seg.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    
    """ 
    Modify the number of sequences in SegmentSequence. 
    There must be only one. Preserve the sequence that relates to the segment
    of interest and remove any others.
    """
    
    # Over-write the first SegmentSequence with the sequence to be preserved:
    Seg.SegmentSequence[0] = deepcopy(Seg.SegmentSequence[SegNum])
    
    OrigSS = len(Seg.SegmentSequence)
    
    ReqSS = 1
    
    if OrigSS > ReqSS:
        for i in range(OrigSS - ReqSS):
            Seg.SegmentSequence.pop()
    
    
    """ 
    Modify the number of sequences in Per-FrameFunctionGroupsSequence. 
    The number of sequences must be equal to the length of PFFGStoSliceInds.
    """

    OrigPFFGS = len(Seg.PerFrameFunctionalGroupsSequence)
    
    ReqPFFGS = len(PFFGStoSliceInds)
    
    if ReqPFFGS > OrigPFFGS:
        for i in range(ReqPFFGS - OrigPFFGS):
            Seg.PerFrameFunctionalGroupsSequence\
               .append(deepcopy(Seg.PerFrameFunctionalGroupsSequence[-1]))          
    else:
        for i in range(OrigPFFGS - ReqPFFGS):
            Seg.PerFrameFunctionalGroupsSequence.pop()

    
    
    #print(f'\n\n***Seg.SegmentSequence[0].SegmentLabel before change = {Seg.SegmentSequence[0].SegmentLabel}')
    
    """ 
    Modify the Segment Label to SegLabel.
    """
    Seg.SegmentSequence[0].SegmentLabel = SegLabel # 09/12/2020
    
    #print(f'\n\n***Seg.SegmentSequence[0].SegmentLabel after change = {Seg.SegmentSequence[0].SegmentLabel}')
    
    
    
    """
    Output some results to the console.
    """
    if LogToConsole:       
        NewRIS = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)
        
        NewSS = len(Seg.SegmentSequence)
        
        NewPFFGS = len(Seg.PerFrameFunctionalGroupsSequence)
        
        print(f'\nThere are {ReqRIS} DICOMs and there were {OrigRIS}',
              f'sequences in ReferencedSeriesSequence. Now there are {NewRIS}',
              'sequences.')
        
        print(f'\nThere were {OrigSS} sequences in SegmentSequence.',
              f'Now there are {NewSS} sequences.')
        
        print(f'\nThere are {ReqPFFGS} segmentations and there were',
              f'{OrigPFFGS} sequences in PerFrameFunctionalGroupsSequence.',
              f'Now there are {NewPFFGS} sequences.')
        
    return Seg

    
    
    







def ModifySeg(Seg, PixArr, PFFGStoSliceInds, DicomDir, LogToConsole=False):
    """
    Modify various tag values of a SEG so that they are consistent with the
    corresponding tag values of the DICOMs.
         
    
    Inputs:
    ******
                              
    Seg : Pydicom object
        The SEG ROI object to modify.
     
    PixArr : Numpy array
        The SEG's pixel array.
    
    PFFGStoSliceInds : List of integers
        List of slice numbers that correspond to each frame in PixArr.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    Seg : Pydicom object
        Modified SEG ROI object.
    """
    
    from DicomTools import ImportDicoms
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    import numpy as np
    
    F, R, C = PixArr.shape
    
    #if F == 0:
    #    PixArr = np.reshape(PixArr, (1, R, C))
    #    
    #    F, R, C = PixArr.shape
    
    if LogToConsole:
        print(f'\n*****PixArr.shape = {PixArr.shape}')
    
    P = len(PFFGStoSliceInds)
    
    if not F == P:
        raise Exception(f"There are {F} frames in PixArr and {P} indices in "\
                        + "PFFGStoSliceInds. They must be equal.")
    
    
    Dicoms = ImportDicoms(DicomDir)
    
    RIS = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)
             
    if not RIS == len(Dicoms):
        raise Exception(f"There are {len(Dicoms)} DICOMs and {RIS} sequences "\
                        + "in ReferencedInstanceSequence. They must be equal.")
    
    # Modify ReferencedInstanceSequence:
    for i in range(len(Dicoms)):
        Seg.ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
    
    
    Seg.NumberOfFrames = f"{F}"
    
    # Modify SegmentSequence:
    Seg.SegmentSequence[0].SegmentNumber = 1
    
    # Get the ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence:
    RefSOPs = GetRSOPuidsInRIS(Seg)
    
    # Modify PerFrameFunctionalGroupsSequence:
    for i in range(P):
        # The slice index for the i^th sequence:
        s = PFFGStoSliceInds[i]
        
        if LogToConsole:
            print('\n\nResults of ModifySeg:')
            print(f'   PFFGStoSliceInds[{i}] = {PFFGStoSliceInds[i]}')
            print(f'   DICOM slice number = {s}')
        
        Seg.PerFrameFunctionalGroupsSequence[i]\
           .DerivationImageSequence[0]\
           .SourceImageSequence[0]\
           .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        # Modify the DimensionIndexValues (DIVs), ImagePositionPatient and
        # ReferencedSegmentNumber:
        """
        Since there is always only one segment, the first integer in the DIV 
        will always be 1, as is the Referenced Segment Number. The second 
        integer will be the 1 + the index of the SOP Instance UID for the s^th 
        DICOM within RefSOPs (+1 since ind is an integer counting from 0, 
        whereas the DIVs are integers counting from 1).
        """
        ind = RefSOPs.index(Dicoms[s].SOPInstanceUID)
        
        Seg.PerFrameFunctionalGroupsSequence[i]\
           .FrameContentSequence[0]\
           .DimensionIndexValues = [1, ind + 1]
        
        Seg.PerFrameFunctionalGroupsSequence[i]\
           .PlanePositionSequence[0]\
           .ImagePositionPatient = Dicoms[s].ImagePositionPatient
           
        Seg.PerFrameFunctionalGroupsSequence[i]\
           .SegmentIdentificationSequence[0]\
           .ReferencedSegmentNumber = 1
           
        if LogToConsole:
              print(f'   DimensionIndexValues = {[1, ind + 1]}')
              print(f'   ReferencedSegmentNumber = {1}')
              
        
        
    
    
    if LogToConsole:
        print('\n   Prior to ravelling and packing bits:')
        print(f'   PixArr.shape = {PixArr.shape}')
        [print(f'   np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
    
    

    # Ravel (to 1D array) PixArr and pack bits (see
    # https://github.com/pydicom/pydicom/issues/1230):
    packed = pack_bits(PixArr.ravel())
    
    # Convert PixArr to bytes:
    #Seg.PixelData = PixArr.tobytes() <-- doesn't work
    Seg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    
    if False:
        F, R, C = PixArr.shape
        
        # If there are more than one frames in PixArr:
        if F > 1:
            # Ravel (to 1D array) PixArr and pack bits (see
            # https://github.com/pydicom/pydicom/issues/1230):
            packed = pack_bits(PixArr.ravel())
            
            # Convert PixArr to bytes:
            Seg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
            
            if LogToConsole:
                print(f'\n   len(packed) = {len(packed)}') 
                """ 
                = 39936 for UseCase 1a
                = 6656 for UseCase 1b
                """
        else:
            # Convert PixArr to bytes:
            #Seg.PixelData = PixArr.tobytes()
            
            #if LogToConsole:
                #print(f'\n   len(PixArr.tobytes()) = {len(PixArr.tobytes())}') 
                #""" 
                #=  for UseCase 1a
                #= 212992 for UseCase 1b
                #"""
            
            #PixArr = np.reshape(PixArr, (R, C))
            
            if LogToConsole:
                print('\n   There is only one frame in PixArr, so it was reshaped',
                      f'to shape {PixArr.shape}')
                print(f'   np.amax(PixArr) = {np.amax(PixArr)}')
            
            # Ravel (to 1D array) PixArr and pack bits (see
            # https://github.com/pydicom/pydicom/issues/1230):
            packed = pack_bits(PixArr.ravel())
            
            # Convert PixArr to bytes:
            Seg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
            
            if LogToConsole:
                print(f'\n   len(packed) = {len(packed)}') 
                """ 
                = 39936 for UseCase 1a
                = 6656 for UseCase 1b
                """
        
        
    
    
    
    if LogToConsole:
        print(f'\n   len(packed) = {len(packed)}') 
        """ 
        = 39936 for UseCase 1a
        = 6656 for UseCase 1b
        """
                
        PixArr2 = Seg.pixel_array
        
        print('\n   After packing bits:')
        print(f'   PixArr2.shape = {PixArr2.shape}')
        if len(PixArr2.shape) > 2:
            [print(f'   np.amax(PixArr2[{i}]) = {np.amax(PixArr2[i])}') for i in range(PixArr2.shape[0])]
        else:
            print(f'   np.amax(PixArr2) = {np.amax(PixArr2)}')
    
                                          
    return Seg







def ErrorCheckSeg_OLD(Seg, DicomDir, LogToConsole=False):
    """
    Check a SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the SEG.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Seg.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe any errors that are found.  If no 
        errors are found an empty list ([]) will be returned.
        
    Nerrors : integer
        The number of errors found.
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
    
    # The maximum allowed error between two lists/items:
    epsilon = 1e-05
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    Dicom = ImportDicom(DicomDir)
    
    #Size, Spacings, ST,\
    #IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=DicomDir, Package='pydicom')
    
    # The number of frames, rows and columns in the SEG's pixel array:
    if len(Seg.pixel_array.shape) > 2:
        PixArrF, PixArrR, PixArrC = Seg.pixel_array.shape
    else:
        PixArrR, PixArrC = Seg.pixel_array.shape
        PixArrF = 1
    
    #NumOfDicoms = len(SOPuids)
    
    LogList = []
    Nerrors = 0
    
    
    """Determine whether MediaStorageSOPInstanceUID matches SOPInstanceUID."""
    MSSOPuid = deepcopy(Seg.file_meta.MediaStorageSOPInstanceUID)
    SOPuid = deepcopy(Seg.SOPInstanceUID)
    
    if MSSOPuid == SOPuid:
        msg = f'INFO:  MediaStorageSOPInstanceUID {MSSOPuid} matches '\
              + f'SOPInstanceUID {SOPuid}.\n'
    else:
        msg = f'ERROR:  MediaStorageSOPInstanceUID {MSSOPuid} does not '\
              + f'match SOPInstanceUID {SOPuid}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether the number of sequences in 
    ReferencedInstanceSequence matches the number of DICOMs."""
    
    RIS = deepcopy(Seg.ReferencedSeriesSequence[0]\
                      .ReferencedInstanceSequence)
    
    if len(SOPuids) == len(RIS):
        msg = f'INFO:  The number of SOPInstanceUIDs {len(SOPuids)}'\
              + ' matches the number of sequences in '\
              + f'ReferencedInstanceSequence {len(RIS)}.\n'
    else:
        msg = f'ERROR:  The number of SOPInstanceUIDs {len(SOPuids)} does'\
              + ' not match the number of sequences in '\
              + f'ReferencedInstanceSequence {len(RIS)}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs do not
    match the SOPInstanceUIDs."""        
    RefSOPuidsInRIS = [RIS[i].ReferencedSOPInstanceUID for i in range(len(RIS))]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInRIS]
    
    NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if Inds:
        NonMatchingUids = [RefSOPuidsInRIS[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + f'ReferencedInstanceSequence {i+1} does not match any '\
                  + 'SOPInstanceUID.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The ReferencedSOPInstanceUIDs in '\
                  + f'ReferencedInstanceSequence matches the '\
                  + 'SOPInstanceUIDs.\n'
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
                
                
    
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
            
            msg = f'ERROR:  SOPInstanceUID {uid} is not referenced in any '\
                  + 'ReferencedSOPInstanceUID in '\
                  + 'ReferencedInstanceSequence.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The SOPInstanceUIDs match the '\
                  + 'ReferencedSOPInstanceUIDs in '\
                  + 'ReferencedInstanceSequence.\n'
                  
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
                           
    if Seg.ReferencedSeriesSequence[0].SeriesInstanceUID == Dicom.SeriesInstanceUID:
        msg = f'INFO:  SEG SeriesInstanceUID in ReferencedSeriesSequence '\
              + f'{Seg.ReferencedSeriesSequence[0].SeriesInstanceUID} '\
              + 'matches DICOM SeriesInstanceUID '\
              + f'{Dicom.SeriesInstanceUID}.\n'
    else:
        msg = f'ERROR:  SEG SeriesInstanceUID in ReferencedSeriesSequence '\
              + f'{Seg.ReferencedSeriesSequence[0].SeriesInstanceUID} does '\
              + 'not match DICOM SeriesInstanceUID '\
              + f'{Dicom.SeriesInstanceUID}.\n'
        
        Nerrors = Nerrors + 1 
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
            
       
    if Seg.PatientName == Dicom.PatientName:
        msg = f'INFO:  SEG PatientName {Seg.PatientName} matches DICOM '\
              + f'PatientName {Dicom.PatientName}.\n'
    else:
        msg = f'ERROR:  SEG PatientName {Seg.PatientName} does not match '\
              + f'DICOM PatientName {Dicom.PatientName}.\n'
            
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
        
    if Seg.PatientID == Dicom.PatientID:
        msg = f'INFO:  SEG PatientID {Seg.PatientID} matches the '\
              + f'DICOM PatientID {Dicom.PatientID}.\n'
    else:
        msg = f'ERROR:  SEG PatientID ({Seg.PatientID}) does not match the '\
              + f'DICOM PatientID ({Dicom.PatientID}).\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)

    if Seg.StudyInstanceUID == Dicom.StudyInstanceUID:
        msg = f'INFO:  SEG StudyInstanceUID {Seg.StudyInstanceUID} '\
              + 'matches the DICOM StudyInstanceUID '\
              + f'{Dicom.StudyInstanceUID}.\n'
    else:
        msg = f'ERROR:  SEG StudyInstanceUID {Seg.StudyInstanceUID} does '\
              + 'not match the DICOM StudyInstanceUID '\
              + f'{Dicom.StudyInstanceUID}.\n'
                      
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    #if Seg.StudyID == Dicom.StudyID:
    #    msg = f'INFO:  The SEG StudyID {Seg.StudyID} matches the'\
    #          + f' DICOM StudyID {Dicom.StudyID}.'
    #else:
    #    msg = f'ERROR:  The SEG StudyID {Seg.StudyID} does not match the'\
    #          + f' DICOM StudyID {Dicom.StudyID}.'
    #
    #    Nerrors = Nerrors + 1
    #    
    #LogList.append(msg)
    #    
    #if LogToConsole:
    #    print(msg)
            
    if Seg.FrameOfReferenceUID == Dicom.FrameOfReferenceUID:
        msg = f'INFO:  SEG FrameOfReferenceUID {Seg.FrameOfReferenceUID} '\
              + 'matches DICOM FrameOfReferenceUID '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
    else:
        msg = f'ERROR:  SEG FrameOfReferenceUID {Seg.FrameOfReferenceUID} '\
              + 'does not match DICOM FrameOfReferenceUID '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    if int(Seg.NumberOfFrames) == PixArrF:
        msg = f'INFO:  NumberOfFrames {Seg.NumberOfFrames} matches the '\
              + f'number of frames in PixelData {PixArrF}.\n'
    else:
        msg = f'ERROR:  NumberOfFrames {Seg.NumberOfFrames} does not '\
              + f'match the number of frames in PixelData {PixArrF}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    if Seg.Rows == Dicom.Rows:
        msg = f'INFO:  SEG Rows {Seg.Rows} matches DICOM Rows '\
              + f'{Dicom.Rows}.\n'
    else:
        msg = f'ERROR:  SEG Rows {Seg.Rows} does not match DICOM Rows '\
              + f'{Dicom.Rows}.\n'
              
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
        
    if Seg.Columns == Dicom.Columns:
        msg = f'INFO:  SEG Columns {Seg.Columns} matches DICOM '\
              + f'Columns {Dicom.Columns}.\n'
    else:
        msg = f'ERROR:  SEG Columns {Seg.Columns} does not match DICOM '\
              + f'Columns {Dicom.Columns}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """There should only be one sequence in SegmentSequence."""
    N = len(Seg.SegmentSequence)
    
    if N == 1:
        msg = f'INFO:  There are {N} sequences in SegmentSequence.  There '\
              + f'should be only 1.\n'
    else:
        msg = f'ERROR:  There are {N} sequences in SegmentSequence.  There '\
              + f'should be only 1.\n'
        
        Nerrors = Nerrors + 1
              
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Verify that the SegmentNumber in SegmentSequence is 1."""
    N = deepcopy(Seg.SegmentSequence[0].SegmentNumber)
    
    if N == 1:
        msg = f'INFO:  SegmentNumber in SegmentSequence is {N}.  It '\
              + 'should be 1.\n'
    else:
        msg = f'ERROR:  SegmentNumber in SegmentSequence is {N}.  It '\
              + 'should be 1.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
            
    """Check various tags in SharedFunctionalGroupsSequence."""
    SegIOP = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                         .PlaneOrientationSequence[0]\
                         .ImageOrientationPatient)
    
    SegIOP = [float(item) for item in SegIOP]
    
    DcmIOP = [float(item) for item in Dicom.ImageOrientationPatient]
    
    if SegIOP == DcmIOP:
        msg = 'INFO:  SEG ImageOrientationPatient in '\
              + f'SharedFunctionalGroupsSequence {SegIOP} matches '\
              + f'DICOM ImageOrientationPatient {DcmIOP}.\n'
    else:
        if AreListsEqualToWithinEpsilon(SegIOP, DcmIOP, epsilon):
            msg = 'INFO:  SEG ImageOrientationPatient in '\
                  + f'SharedFunctionalGroupsSequence {SegIOP} is within '\
                  + f'{epsilon} of DICOM ImageOrientationPatient {DcmIOP}.\n'
        else:
            msg = 'ERROR:  SEG ImageOrientationPatient in '\
                  + f'SharedFunctionalGroupsSequence {SegIOP} is not within '\
                  + f'{epsilon} of DICOM ImageOrientationPatient {DcmIOP}.\n'
            
            Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
       
    
    SegST = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                        .PixelMeasuresSequence[0]\
                        .SliceThickness)
    
    """ The SliceThickness appears to be the z-Spacing rather than the
    SliceThickness from the DICOM metadata. """
    
    if float(SegST) == Spacings[2]:
        msg = 'INFO:  SEG SliceThickness in SharedFunctionalGroupsSequence'\
              + f' {SegST} matches the slice thickness calculated '\
              + 'from DICOM ImagePositionPatient and '\
              + f'ImageOrientationPatient {Spacings[2]}.\n'
    else:
        if AreItemsEqualToWithinEpsilon(float(SegST), Spacings[2], epsilon):
            msg = 'INFO:  SEG SliceThickness in '\
                  + f'SharedFunctionalGroupsSequence {SegST} is within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM'\
                  + ' ImagePositionPatient and ImageOrientationPatient '\
                  + f'{Spacings[2]}.\n'
        else:
            print(f'\n{float(SegST) - Spacings[2]}')
            
            msg = 'ERROR:  SEG SliceThickness in '\
                  + f'SharedFunctionalGroupsSequence {SegST} is not within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient '\
                  + f'{Spacings[2]}.\n'
            
            Nerrors = Nerrors + 1
              
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    SegSBS = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                         .PixelMeasuresSequence[0]\
                         .SpacingBetweenSlices)
    
    if float(SegSBS) == Spacings[2]:
        msg = 'INFO:  SEG SpacingBetweenSlices in '\
              + f'SharedFunctionalGroupsSequence {SegSBS} matches the slice '\
              + 'thickness calculated from DICOM ImagePositionPatient and '\
              + f'ImageOrientationPatient {Spacings[2]}.\n'
    else:
        if AreItemsEqualToWithinEpsilon(float(SegSBS), Spacings[2], epsilon):
            msg = 'INFO:  SEG SpacingBetweenSlices in '\
                  + f'SharedFunctionalGroupsSequence {SegSBS} is within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient '\
                  + f'{Spacings[2]}.\n'
        else:
            print(f'\n{float(SegSBS) - Spacings[2]}')
            
            msg = 'ERROR:  SEG SpacingBetweenSlices in '\
                  + f'SharedFunctionalGroupsSequence {SegSBS} is not within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient '\
                  + f'{Spacings[2]}.\n'
        
            Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)

    
    SegPS = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                        .PixelMeasuresSequence[0]\
                        .PixelSpacing)
    
    SegPS = [float(item) for item in SegPS]
    
    DcmPS = [float(item) for item in Dicom.PixelSpacing]
    
    if SegPS == DcmPS:
        msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence '\
              + f'{SegPS} matches the DICOM PixelSpacing {DcmPS}.\n'
    else:
        if AreListsEqualToWithinEpsilon(SegPS, DcmPS, epsilon):
            msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence '\
                  + f'{SegPS} is within {epsilon} of the DICOM PixelSpacing '\
                  + f'{DcmPS}.\n'
        else:
            msg = 'ERROR:  SEG PixelSpacing in SharedFunctionalGroupsSequence'\
                  + f'{SegPS} is not within {epsilon} of the DICOM '\
                  + f'PixelSpacing {DcmPS}.\n'
            
            Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    
    
    """Verify that the number of sequences in PerFrameFunctionalGroupsSequence
    is equal to NumberOfFrames."""
    PFFGS = deepcopy(Seg.PerFrameFunctionalGroupsSequence)
    
    if len(PFFGS) == int(Seg.NumberOfFrames):
        msg = f'INFO:  The number of sequences in '\
              + f'PerFrameFunctionGroupsSequence {len(PFFGS)} '\
              + f'matches NumberOfFrames {Seg.NumberOfFrames}.\n'
    else:
        msg = f'ERROR:  The number of sequences in '\
              + f'PerFrameFunctionGroupsSequence {len(PFFGS)} does not '\
              + f'match NumberOfFrames {Seg.NumberOfFrames}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
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
    
    if Inds:
        NonMatchingUids = [RefSOPuidsInPFFGS[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + 'PerFrameFunctionalGroupsSequence {i+1} does not match '\
                  + 'any DICOM SOPInstanceUID.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The ReferencedSOPInstanceUIDs in '\
              + 'PerFrameFunctionalGroupsSequence match the DICOM '\
              + 'SOPInstanceUIDs.\n'
                  
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
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
            msg = f'ERROR:  The first element in DimensionalIndexValue {i+1},'\
                  + f' {divs[i]}, is not "1".\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The first elements in DimensionalIndexValue are 1.\n'
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
                
    # Determine if any of the second elements in the DIVs exceed the number
    # of sequences in ReferencedInstanceSequence:
    SecondElements = [int(DIVs[i][1]) for i in range(len(DIVs))]
    
    # Find the indices of any elements that exceed len(RIS):
    Inds = [i for i, x in enumerate(SecondElements) if x > len(RIS)]

    if Inds:
        #divs = [DIVs[i] for i in Inds]
        
        for i in range(len(Inds)):
            msg = 'ERROR:  The second element in DimensionalIndexValue '\
                  + f'{i+1}, {DIVs[i]} exceeds the number of sequences in '\
                  + f'ReferencedInstanceSequence {len(RIS)}.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The second elements in DimensionalIndexValue '\
                  + f'matches the number of sequences in '\
                  + f'ReferencedInstanceSequence.\n'
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
                
    """Determine if the ReferencedSOPInstanceUID in the 
    ReferencedInstanceSequence indexed by the second elements in the DIVs 
    match the ReferencedSOPInstanceUID in the SourceImageSequence."""
    for i in range(len(SecondElements)):
        RefSOPinPFFGS = RefSOPuidsInPFFGS[i]
        
        # Need to -1 since i is zero-indexed:
        ind = SecondElements[i] - 1
        
        RefSOPinRIS = RefSOPuidsInRIS[ind]
        
        if RefSOPinPFFGS == RefSOPinRIS:
            msg = 'INFO:  ReferencedSOPInstanceUID referenced in '\
                  + f'SourceImageSequence {i+1} matches '\
                  + 'ReferencedSOPInstanceUID referenced in '\
                  + 'ReferencedInstanceSequence as indexed by the second '\
                  + f'element in DimensionalIndexValue {i+1}, {DIVs[i]}.\n'
        else:
            msg = 'ERROR:  ReferencedSOPInstanceUID referenced in '\
                  + f'SourceImageSequence {i+1} does not match '\
                  + 'ReferencedSOPInstanceUID referenced in '\
                  + 'ReferencedInstanceSequence as indexed by the second '\
                  + f'element in DimensionalIndexValue {i+1}, {DIVs[i]}.\n'
            
            Nerrors = Nerrors + 1
                  
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
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
        if AreListsEqualToWithinEpsilon(IPPsInPFFGS[i], IPP, epsilon):
            msg = f'INFO: ImagePositionPatient {IPPsInPFFGS[i]} '\
                  + f'referenced in PlanePositionSequence {i+1} is within '\
                  + f'{epsilon} of the DICOM ImagePositionPatient {IPP} as '\
                  + f'indexed by the second element {SecondElements[i]} in '\
                  + f'DimensionalIndexValue {i+1}, {DIVs[i]}.\n'
        else:
            msg = f'ERROR: ImagePositionPatient {IPPsInPFFGS[i]} '\
                  + f'referenced in PlanePositionSequence {i+1} is not within'\
                  + f' {epsilon} of the DICOM ImagePositionPatient {IPP} as '\
                  + f'indexed by the second element {SecondElements[i]} in '\
                  + f'DimensionalIndexValue {i+1}, {DIVs[i]}.\n'
            
            Nerrors = Nerrors + 1
                  
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether any of the ReferencedSegmentNumber are not 1 (there
    should only be one segment)."""
    RSNs = [PFFGS[i].SegmentIdentificationSequence[0]\
                    .ReferencedSegmentNumber for i in range(len(PFFGS))]
    
    # Find the indices of any non-matching elements:
    Inds = [i for i, x in enumerate(RSNs) if x != 1]
    
    if Inds:        
        for i in range(len(Inds)):
            msg = f'ERROR:  ReferencedSegmentNumber is {i+1}. '\
                  + 'It should be 1.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  ReferencedSegmentNumber is {i+1}. It should be 1.'
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine if the dimensions of the pixel array in PixelData do not
    match the expected dimensions (the number of frames (PixArrF) has 
    already been compared to NumberOfFrames)."""
    if PixArrR == Seg.Rows:
        msg = f'INFO:  The number of rows in PixelData {PixArrR} '\
              + f'matches the number of SEG Rows {Seg.Rows}.\n'
    else:
        msg = f'ERROR:  The number of rows in PixelData {PixArrR} does '\
              + f'not match the number of SEG Rows {Seg.Rows}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    if PixArrC == Seg.Columns:
        msg = f'INFO:  The number of columns in PixelData {PixArrC} '\
              + f'matches the number of SEG Columns {Seg.Columns}.\n'
    else:
        msg = f'ERROR:  The number of columns in PixelData {PixArrC} does'\
              + f' not match the number of SEG Columns {Seg.Columns}.\n'
        
        Nerrors = Nerrors + 1
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    print(f'\nThere were {Nerrors} errors found in the SEG.')
            
    return LogList, Nerrors







def GetSegDataFromListOfSegs_OLD(ListOfSegs, ListOfDicomDirs, SearchString, 
                             LogToConsole=False):
    """ 
    Note:
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    """
    
    import importlib
    import SegTools
    importlib.reload(SegTools)
    import DicomTools
    importlib.reload(DicomTools)
    
    import numpy as np
    from SegTools import GetPFFGStoSliceIndsBySeg
    from DicomTools import GetDicomFpaths
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    from SegTools import GetFramesFromPixArr
    
    ListOfPixArrs = []
    ListOfF2Sinds = []
    ListOfSegNums = []
    ListOfDcmFpaths = []
    
    for i in range(len(ListOfSegs)):
        if LogToConsole:
            print(f'\n\nListOfSegs[{i}]:')
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        ListOfDcmFpaths.append(DcmFpaths)
        
        """
        11/12/20:  
            There should only be one segment so perhaps this should be checked
            and raise exception if the number of segments is greater than one.
            And if not, SegNum will always be = 0 (no need to GetRoiNum).
        """
        
        if ListOfSegs[i]:
            SegNum = GetRoiNums(ListOfSegs[i], SearchString)
            
            SOPuids = GetDicomSOPuids(ListOfDicomDirs[i])
        
            PFFGS2SindsBySeg = GetPFFGStoSliceIndsBySeg(ListOfSegs[i], 
                                                        SOPuids)
            
            ListOfF2Sinds.append(PFFGS2SindsBySeg[SegNum])
            
            if LogToConsole:
                print(f'\nListOfF2Sinds[{i}] = {ListOfF2Sinds[i]}')
            
            
            
            #PixArr = ListOfSegs[i].pixel_array # this contains all segments!
            PixArr = GetFramesFromPixArr(ListOfSegs[i], ListOfDicomDirs[i], 
                                         SearchString, LogToConsole)
            
            if LogToConsole:      
                print(f'\nPixArr.shape = {PixArr.shape}')
                [print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
                
            """
            If PixArr has a single frame it will have shape (NumOfRows, NumOfCols).
            For compatibility with multi-framed pixel arrays, convert single-framed
            pixel arrays to shape (1, NumOfRows, NumOfCols).
            """
            if len(PixArr.shape) < 3:
                R, C = PixArr.shape
                
                PixArr = np.reshape(PixArr, (1, R, C))
                
                if LogToConsole:      
                    print(f'\nPixArr was a single frame so reshaped to {PixArr.shape}')
                    
        else:
            SegNum = None
            ListOfF2Sinds.append(None)
            PixArr = None
        
        ListOfSegNums.append(SegNum)
        
        ListOfPixArrs.append(PixArr)
        
        #if LogToConsole: 
        #    print('')
        #    print(f'\nListOfPixArrs[{i}].shape = {ListOfPixArrs[i].shape}')
        
    
    return ListOfPixArrs, ListOfF2Sinds, ListOfSegNums, ListOfDcmFpaths





