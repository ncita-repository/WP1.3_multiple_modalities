# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:01:08 2020

@author: ctorti
"""



"""
Old functions from RoiCopyTools.py and others
"""


def ModifyRtsTagVals(TrgRts, TrgDicoms, NewTrgContourData, NewTrgRoiNums,
                     NewTrgCStoSliceInds, FromRoiNum, FromSliceNum, SrcRts,
                     LogToConsole=False):
    """
    Modify various tag values of the new RTS ROI so that they are consistent
    with the corresponding tag values of the DICOMs (e.g. SOPClassUID, 
    SOPInstanceUID), and so that the ROIContourSequence values are consistent
    with those of the original RTS ROI for existing contours, and with the
    contour that has been copied from the original RTS ROI.
         
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
    ------
                              
    TrgRts : Pydicom object
        Original Target RTS ROI Object.
    
    TrgDicoms : List of Pydicom objects
        List of DICOM Objects that include the DICOMs that TrgRts relate to.
    
    NewTrgContourData : list of floats
        Flat list of coordinates for Target containing the new ContourData.
    
    NewTrgRoiNums : List of integers
        ROI numbers that correspond to each ROI in NewTrgRts.
    
    NewTrgCStoSliceInds : List of integers
        Slice numbers that correspond to each contour in NewTrgRts.
                           
    FromRoiNum : integer
        Source ROI number of the contour that was copied (counting from 0) 
        (this is only used for generating a new Series Description).
    
    FromSliceNum : integer
        Source slice index in the DICOM stack corresponding to the contour 
        that was copied (counting from 0) (this is only used for generating a 
        new Series Description).
                           
    SrcRts : Pydicom object
        Source RTS ROI object (this is only used for generating a new Series 
        Description).
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    NewTrgRts : Pydicom object
        New Target RTS ROI object.
    """
    
    # Get the sequence number that corresponds to the slice number to be 
    # copied, and the slice number to be copied to:
    FromSeqNum = OrigCIStoDcmInds.index(FromSliceNum)
    ToSeqNum = NewCIStoDcmInds.index(ToSliceNum)
    
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
        
        NewRts.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .ContourImageSequence[i]\
              .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        """ Ditto:
        NewRts.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourImageSequence[0]\
              .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewRts.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourImageSequence[0]\
              .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
                 
        
        
        # Modify the Contour Number to match the value of the Source contour
        #  being copied:
        NewRts.ROIContourSequence[0]\
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
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[FromOrigSeqNum]\
                                  .NumberOfContourPoints)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .NumberOfContourPoints = deepcopy(val)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[FromOrigSeqNum]\
                                  .ContourData)
            
            NewRts.ROIContourSequence[0]\
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
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[ind]\
                                  .NumberOfContourPoints)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .NumberOfContourPoints = deepcopy(val)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[ind]\
                                  .ContourData)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourData = deepcopy(val)
                  
                  
    # Generate a new SOP Instance UID:
    NewTrgRts.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRts.SeriesInstanceUID = generate_uid()
                                          
    return NewRts










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








def VerifyCISandCS(CIStoDcmInds, CStoDcmInds, LogToConsole=False):
    """
    CIS and CS will only be the same if there is only one ROI so this function
    is of limited use.
    
    
    
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





def ResampleImage(Image, RefImage, Method, Interpolation):#, RefImageSpacings):
    
    #NewSrcImage = sitk.Resample(image=SrcImage, 
    #                            size=TrgSitkSize,
    #                            transform=sitk.Transform(),
    #                            interpolator=sitk.sitkLinear,
    #                            outputOrigin=TrgSitkOrigin,
    #                            outputSpacing=TrgSitkSpacing,
    #                            outputDirection=TrgSitkDirection,
    #                            defaultPixelValue=0,
    #                            outputPixelType=TrgImage.GetPixelID())
    
    
    # Define which interpolator to use:
    """ 
    Use the linear (or BSpline) interpolator for intensity images, and
    NearestNeighbour for binary images (e.g. segmentation) so that no new
    labels are introduced.
    """
    
    if 'inear' in Interpolation:
        Interpolator = sitk.sitkLinear
    if 'pline' in Interpolation:
        Interpolator = sitk.sitkBSpline
    if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:    
        Interpolator = sitk.sitkNearestNeighbor
            
        
    #print('\nUsing', Interpolation, 'interpolation\n')
    
    Dimension = Image.GetDimension()
    

    
    if 'ResampleImageFilter' in Method:
        Resampler = sitk.ResampleImageFilter()
        
        Resampler.SetReferenceImage(RefImage)
        
        Resampler.SetInterpolator(Interpolator)
        
        if 'Identity' in Method:
            Resampler.SetTransform(sitk.Transform())
        if 'Affine' in Method:
            Resampler.SetTransform(sitk.AffineTransform(Dimension))
            
            
        #print(f'\nImage.GetPixelIDValue() = {Image.GetPixelIDValue()}')
        
        Resampler.SetOutputSpacing(RefImage.GetSpacing())
        Resampler.SetSize(RefImage.GetSize())
        Resampler.SetOutputDirection(RefImage.GetDirection())
        Resampler.SetOutputOrigin(RefImage.GetOrigin())
        #Resampler.SetDefaultPixelValue(Image.GetPixelIDValue())
        Resampler.SetDefaultPixelValue(0)
        
        #print('\n')
        #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: print("\rProgress: {0:03.1f}%...".format(100*Resampler.GetProgress()),end=''))
        #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: sys.stdout.flush())
        
        ResImage = Resampler.Execute(Image)
    
    
    
    if 'Resample2' in Method:
        """
        Use the 2nd method listed in the SimpleITK Transforms and Resampling 
        tutorial (the new Size, Spacing, Origin and Direction will be 
        determined intrinsically from RefImage):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        
        Also follow the example here:
            
        https://gist.github.com/zivy/79d7ee0490faee1156c1277a78e4a4c4
        """
        
        # Define the transform which maps RefImage to Image with the 
        # translation mapping the image origins to each other:
        Transform = sitk.AffineTransform(Dimension)
        #Transform.SetMatrix(Image.GetDirection())                                      # <-- 9:55 on 30/10
        Transform.SetMatrix(RefImage.GetDirection())                                    # <-- 9:55 on 30/10
        #Transform.SetTranslation(np.array(Image.GetOrigin()) - RefImage.GetOrigin())   # <-- 9:55 on 30/10
        Transform.SetTranslation(np.array(RefImage.GetOrigin()) - Image.GetOrigin())    # <-- 9:55 on 30/10
        
        if 'Centres' in Method:
            # Use the TransformContinuousIndexToPhysicalPoint to compute an 
            # indexed point's physical coordinates as this takes into account  
            # size, spacing and direction cosines. 
            """
            For the vast majority of images the direction cosines are the 
            identity matrix, but when this isn't the case simply multiplying 
            the central index by the spacing will not yield the correct 
            # coordinates resulting in a long debugging session. 
            """
            RefCentre = np.array(RefImage.TransformContinuousIndexToPhysicalPoint(np.array(RefImage.GetSize())/2.0))
            
            # Modify the transformation to align the centres of the original 
            # and reference image instead of their origins:
            CenteringTransform = sitk.TranslationTransform(Dimension)
            
            ImageCentre = np.array(Image.TransformContinuousIndexToPhysicalPoint(np.array(Image.GetSize())/2.0))
            
            CenteringTransform.SetOffset(np.array(Transform.GetInverse().TransformPoint(ImageCentre) - RefCentre))
            
            CenteredTransform = sitk.Transform(Transform)
            
            CenteredTransform.AddTransform(CenteringTransform)

            ResImage = sitk.Resample(Image, RefImage, CenteredTransform, 
                                     Interpolator, Image.GetPixelIDValue())
            
            
        if 'Origins' in Method:
            ResImage = sitk.Resample(Image, RefImage, Transform, Interpolator, 
                                     Image.GetPixelIDValue())
        
        
    
    
    if 'Resample3' in Method:
        """
        Use the 3rd method listed in the SimpleITK Transforms and Resampling 
        tutorial (explicitly define the new Size, Spacing, Origin and 
        Direction):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        """
        
        NewSize = RefImage.GetSize()
        NewOrigin = RefImage.GetOrigin()
        NewSpacing = RefImage.GetSpacing()
        NewDirection = RefImage.GetDirection()
        
        
        # Define the transform which maps RefImage to Image with the 
        # translation mapping the image origins to each other:
        Transform = sitk.AffineTransform(Dimension)
        #Transform.SetMatrix(Image.GetDirection())                                     # <-- 10:03 on 30/10
        Transform.SetMatrix(RefImage.GetDirection())                                   # <-- 10:03 on 30/10
        #Transform.SetTranslation(np.array(Image.GetOrigin()) - RefImage.GetOrigin())  # <-- 10:03 on 30/10
        Transform.SetTranslation(np.array(RefImage.GetOrigin()) - Image.GetOrigin())   # <-- 10:03 on 30/10
        
        if 'Centres' in Method:
            # Use the TransformContinuousIndexToPhysicalPoint to compute an 
            # indexed point's physical coordinates as this takes into account  
            # size, spacing and direction cosines. 
            """
            For the vast majority of images the direction cosines are the 
            identity matrix, but when this isn't the case simply multiplying 
            the central index by the spacing will not yield the correct 
            # coordinates resulting in a long debugging session. 
            """
            RefCentre = np.array(RefImage.TransformContinuousIndexToPhysicalPoint(np.array(RefImage.GetSize())/2.0))
            
            # Modify the transformation to align the centres of the original 
            # and reference image instead of their origins:
            CenteringTransform = sitk.TranslationTransform(Dimension)
            
            ImageCentre = np.array(Image.TransformContinuousIndexToPhysicalPoint(np.array(Image.GetSize())/2.0))
            
            CenteringTransform.SetOffset(np.array(Transform.GetInverse().TransformPoint(ImageCentre) - RefCentre))
            
            CenteredTransform = sitk.Transform(Transform)
            
            CenteredTransform.AddTransform(CenteringTransform)
        
            ResImage = sitk.Resample(Image, NewSize, CenteredTransform, 
                                     Interpolator, NewOrigin, NewSpacing, 
                                     NewDirection, Image.GetPixelIDValue())  
            
            
        if 'Origins' in Method:
            ResImage = sitk.Resample(Image, NewSize, Transform, 
                                     Interpolator, NewOrigin, NewSpacing, 
                                     NewDirection, Image.GetPixelIDValue())  
            
            
            
            
    if Method == 'NotCorrect':
        """ 
        Use the 3rd method listed in the SimpleITK Transforms and Resampling 
        tutorial (explicitly define the new Size, Spacing, Origin and 
        Direction):
            
        https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynb
        """
        
        ExtremePts = GetExtremePoints(RefImage)
        
        Affine = sitk.Euler3DTransform(RefImage.TransformContinuousIndexToPhysicalPoint([size/2 for size in RefImage.GetSize()]),
                                       RefImage.GetAngleX(), 
                                       RefImage.GetAngleY(), 
                                       RefImage.GetAngleZ(), 
                                       RefImage.GetTranslation())
        
        InvAffine = Affine.GetInverse()
        
        # Transformed ExtremePts:
        TxExtremePts = [InvAffine.TransformPoint(pt) for pt in ExtremePts]
        
        min_x = min(TxExtremePts)[0]
        min_y = min(TxExtremePts, key=lambda p: p[1])[1]
        min_z = min(TxExtremePts, key=lambda p: p[2])[2]
        max_x = max(TxExtremePts)[0]
        max_y = max(TxExtremePts, key=lambda p: p[1])[1]
        max_z = max(TxExtremePts, key=lambda p: p[2])[2]
        
        NewSpacing = RefImage.GetSpacing()
        #NewDirection = np.eye(3).flatten()
        NewDirection = RefImage.GetDirection()
        #NewOrigin = [min_x, min_y, min_z]
        NewOrigin = RefImage.GetOrigin()
        
        #NewSize = [math.ceil((max_x - min_x)/NewSpacing[0]), 
        #           math.ceil((max_y - min_y)/NewSpacing[1]),
        #           math.ceil((max_z - min_z)/NewSpacing[2])]
        NewSize = RefImage.GetSize()
        
        ResImage = sitk.Resample(Image, NewSize, Affine, Interpolator, 
                                 NewOrigin, NewSpacing, NewDirection)
    
    
    #sitk.Show(ResImage)

    return ResImage









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
    
    
    
    

    

    
def Plot_Src_ResSrc_Trg_Images_OLD(SrcIm, ResSrcIm, TrgIm, SrcImLabel, TrgImLabel,
                               ImageType, ExportPlot, ExportDir, 
                               LogToConsole=False, dpi=80):
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm.GetSize()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetSize()} \n',
                      f'\n   Target:           {TrgIm.GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm.GetSpacing()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetSpacing()} \n',
                      f'\n   Target:           {TrgIm.GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm.GetOrigin()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetOrigin()} \n',
                      f'\n   Target:           {TrgIm.GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm.GetDirection()} \n',
                      f'\n   Resampled Source: {ResSrcIm.GetDirection()} \n',
                      f'\n   Target:           {TrgIm.GetDirection()}')

    
    
    # Convert from SimpleITK image to Numpy arrays:
    SrcNpa = sitk.GetArrayFromImage(SrcIm)
    ResSrcNpa = sitk.GetArrayFromImage(ResSrcIm)
    TrgNpa = sitk.GetArrayFromImage(TrgIm)
    
    # Prepare the figure:
    Nslices = max(SrcIm.GetSize()[2], TrgIm.GetSize()[2])
    
    # Set the number of subplot rows and columns:
    Ncols = 3
    
    Nrows = Nslices
    
    
    # Create a figure with two subplots and the specified size:
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Loop through each slice: 
    for i in range(Nslices):
        #if LogToConsole:
        #        print('\ni =', i)
        
        #plt.subplot(Nrows, Ncols, n)
        #plt.axis('off')
        #ax = plt.subplot(Nrows, Ncols, n)

        #ax = plt.subplot(Nrows, Ncols, n, aspect=AR)
         
        # Plot the Original Source slice:
        if i < SrcIm.GetSize()[2]:
            #maximum = SrcNpa[i].max()
            maximum = np.amax(SrcNpa[i])
            if maximum <= 1:
                maximum = 1
            
            ax = plt.subplot(Nrows, Ncols, n)
            im = ax.imshow(SrcNpa[i], cmap=plt.cm.Greys_r)
            cbar = fig.colorbar(im, ax=ax)
            cbar.mappable.set_clim(0, maximum)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Original Source slice {i}')
            
        n += 1 # increment sub-plot number
            
        # Plot the Resampled Source slice:
        if i < ResSrcIm.GetSize()[2]:
            #maximum = ResSrcNpa[i].max()
            maximum = np.amax(ResSrcNpa[i])
            if maximum <= 1:
                maximum = 1
                
            ax = plt.subplot(Nrows, Ncols, n)
            im = ax.imshow(ResSrcNpa[i], cmap=plt.cm.Greys_r)
            cbar = fig.colorbar(im, ax=ax)
            cbar.mappable.set_clim(0, maximum)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Resampled Source slice {i}')
            
        n += 1 # increment sub-plot number
        
        # Plot the Target slice:
        if i < TrgIm.GetSize()[2]:
            #maximum = TrgNpa[i].max()
            maximum = np.amax(TrgNpa[i])
            if maximum <= 1:
                maximum = 1
                
            ax = plt.subplot(Nrows, Ncols, n)
            im = ax.imshow(TrgNpa[i], cmap=plt.cm.Greys_r)
            cbar = fig.colorbar(im, ax=ax)
            cbar.mappable.set_clim(0, maximum)
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Target slice {i}')
            
        n += 1 # increment sub-plot number
        
      
    
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
    
    if ExportPlot:
        # Short version of Interpolation for exported filename:
        #if 'inear' in Interpolation:
        #    Interp = 'Lin'
        #if 'pline' in Interpolation:
        #    Interp = 'BSp'
        #if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:    
        #    Interp = 'NN'
            
        # Short version of ResampleImageFilter:
        #if Method == 'ResampleImageFilter':
        #    Method = 'ResampleImFilt'
            
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
    
        #ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
        #              + SrcImLabel + '_' + TrgImLabel + '_images_' \
        #              + Method + '_' + Interp + '.jpg'
                      
        ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
                      + SrcImLabel + '_' + TrgImLabel + '_' + ImageType \
                      + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return







def Plot_Src_ResSrc_Trg_NewTrg_Images(SrcIm, ResSrcIm, TrgIm, NewTrgIm,
                                      SrcImLabel, TrgImLabel,
                                      ImageType, ExportPlot, ExportDir, 
                                      LogToConsole=False, dpi=80):
    
    """ 
    SrcIm/ResSrcIm/TrgIm could either be a SimpleITK image,
    i.e. type(SrcIm).__name__ == 'Image', or a list of SimpleITK images,
    i.e. type(SrcIm).__name__ == 'list'.
    
    e.g. SrcIm can be a list of sitk images of the labelmaps for each of 
    multiple ROIs.
    """
    
    # If ___Im is an Image (i.e. not a list), redefine it as a list with one 
    # item:  
    if type(SrcIm).__name__ == 'Image':
        SrcIm = [SrcIm]
        
    if type(ResSrcIm).__name__ == 'Image':
        ResSrcIm = [ResSrcIm]
        
    if type(TrgIm).__name__ == 'Image':
        TrgIm = [TrgIm]
        
    if type(NewTrgIm).__name__ == 'Image':
        NewTrgIm = [NewTrgIm]
    
    
    if LogToConsole:
        print(f'\nScrIm has {len(SrcIm)} segments')
        print(f'ResScrIm has {len(ResSrcIm)} segments')
        print(f'TrgIm has {len(TrgIm)} segments')
        print(f'NewTrgIm has {len(NewTrgIm)} segments')
    
    
    if LogToConsole:
        print('\nThe Source and Target image size:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSize()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetSize()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetSize()}',
                      f'\n   New Target:       {NewTrgIm[0].GetSize()}')
            
        print('\nThe Source and Target voxel spacings:\n',
                      f'\n   Original Source:  {SrcIm[0].GetSpacing()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetSpacing()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetSpacing()}',
                      f'\n   New Target:       {NewTrgIm[0].GetSpacing()}')
            
        print('\nThe Source and Target Origin:\n',
                      f'\n   Original Source:  {SrcIm[0].GetOrigin()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetOrigin()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetOrigin()}',
                      f'\n   New Target:       {NewTrgIm[0].GetOrigin()}')
            
        print('\nThe Source and Target Direction:\n',
                      f'\n   Original Source:  {SrcIm[0].GetDirection()} \n',
                      f'\n   Resampled Source: {ResSrcIm[0].GetDirection()} \n',
                      f'\n   Original Target:  {TrgIm[0].GetDirection()}',
                      f'\n   New Target:       {NewTrgIm[0].GetDirection()}')

    
    """
    The code below works but it's less efficient than what follows, i.e.
    instead of converting from sitk images to numpy and adding the segments,
    it's more efficient to add the segments as sitk images the convert to 
    numpy.
    
    # Convert from SimpleITK image to Numpy arrays, assigning different values
    # for each ROI (i.e. 1, 2, 3, ...):
    SrcNpas = []
    ResSrcNpas = []
    TrgNpas = []
    NewTrgNpas = []
    
    # Loop through all ROIs:
    for r in range(len(SrcIm)):
        SrcNpas.append(sitk.GetArrayFromImage(SrcIm[r]))
        #SrcNpas.append((r + 1)*sitk.GetArrayFromImage(SrcIm[r]))
        
    for r in range(len(ResSrcIm)):
        ResSrcNpas.append(sitk.GetArrayFromImage(ResSrcIm[r]))
        #ResSrcNpas.append((r + 1)*sitk.GetArrayFromImage(ResSrcIm[r]))
        
    for r in range(len(TrgIm)):
        TrgNpas.append(sitk.GetArrayFromImage(TrgIm[r]))
        #TrgNpas.append((r + 1)*sitk.GetArrayFromImage(TrgIm[r]))
        
    for r in range(len(NewTrgIm)):
        NewTrgNpas.append(sitk.GetArrayFromImage(NewTrgIm[r]))
            
    
    if LogToConsole:
        print(f'\nlen(SrcNpas)    = {len(SrcNpas)}')
        print(f'len(ResSrcNpas) = {len(ResSrcNpas)}')
        print(f'len(TrgNpas)    = {len(TrgNpas)}')
        print(f'len(NewTrgNpas) = {len(NewTrgNpas)}')
        
        print(f'\nSrcNpas[0].shape    = {SrcNpas[0].shape}')
        print(f'ResSrcNpas[0].shape = {ResSrcNpas[0].shape}')
        print(f'TrgNpas[0].shape    = {TrgNpas[0].shape}')
        print(f'NewTrgNpas[0].shape = {NewTrgNpas[0].shape}')
        
    # Take the sum of the labelmaps for all ROIs:
    SrcNpaSum = np.zeros((SrcNpas[0].shape))
    ResSrcNpaSum = np.zeros((ResSrcNpas[0].shape))
    TrgNpaSum = np.zeros((TrgNpas[0].shape))
    NewTrgNpaSum = np.zeros((NewTrgNpas[0].shape))
    
    for r in range(len(SrcNpas)):
        SrcNpaSum += SrcNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(SrcNpaSum) = {np.amax(SrcNpaSum)}')
        
    for r in range(len(ResSrcNpas)):
        ResSrcNpaSum += ResSrcNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(ResSrcNpaSum) = {np.amax(ResSrcNpaSum)}')
        
    for r in range(len(TrgNpas)):
        TrgNpaSum += TrgNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(TrgNpaSum) = {np.amax(TrgNpaSum)}')
        
    for r in range(len(NewTrgNpas)):
        NewTrgNpaSum += NewTrgNpas[r]
        
        if r == 0:
            print('')
        print(f'r = {r}, np.amax(NewTrgNpaSum) = {np.amax(NewTrgNpaSum)}')
        
        
    print(f'\nSrcNpaSum.shape    = {SrcNpaSum.shape}')
    print(f'ResSrcNpaSum.shape = {ResSrcNpaSum.shape}')
    print(f'TrgNpaSum.shape    = {TrgNpaSum.shape}')
    print(f'NewTrgNpaSum.shape = {NewTrgNpaSum.shape}')
    """
    
    
    # Add the segments contained in the list of sitk images into a single sitk
    # image:
    SrcImSum = InitialiseImage(SrcIm[0])
    for r in range(len(SrcIm)):
        SrcImSum = AddImages(SrcImSum, SrcIm[r])
        
    ResSrcImSum = InitialiseImage(ResSrcIm[0])
    for r in range(len(ResSrcIm)):
        ResSrcImSum = AddImages(ResSrcImSum, ResSrcIm[r])
        
    TrgImSum = InitialiseImage(TrgIm[0])
    for r in range(len(TrgIm)):
        TrgImSum = AddImages(TrgImSum, TrgIm[r])
       
    NewTrgImSum = InitialiseImage(NewTrgIm[0])
    for r in range(len(NewTrgIm)):
        NewTrgImSum = AddImages(NewTrgImSum, NewTrgIm[r])
    
        
    # Convert to numpy arrays:
    SrcNpaSum = sitk.GetArrayFromImage(SrcImSum)
    ResSrcNpaSum = sitk.GetArrayFromImage(ResSrcImSum)
    TrgNpaSum = sitk.GetArrayFromImage(TrgImSum)
    NewTrgNpaSum = sitk.GetArrayFromImage(NewTrgImSum)
    
    
    #Nslices = max(SrcIm.GetSize()[2], TrgIm.GetSize()[2])

    #SrcNslices = SrcNpas[0].shape[0]
    #ResSrcNslices = ResSrcNpas[0].shape[0]
    #TrgNslices = TrgNpas[0].shape[0]
    #NewTrgNslices = NewTrgNpas[0].shape[0]
    
    SrcNslices = SrcNpaSum.shape[0]
    ResSrcNslices = ResSrcNpaSum.shape[0]
    TrgNslices = TrgNpaSum.shape[0]
    NewTrgNslices = NewTrgNpaSum.shape[0]
    
    Nslices = max(SrcNslices, TrgNslices)

    
    # Get the maxima for all slices for all labelmaps:
    if 'Labelmaps' in ImageType:
        SrcMaximaBySlice = np.amax(SrcNpaSum, axis=(1,2))
        
        ResSrcMaximaBySlice = np.amax(ResSrcNpaSum, axis=(1,2))
            
        TrgMaximaBySlice = np.amax(TrgNpaSum, axis=(1,2))
        
        NewTrgMaximaBySlice = np.amax(NewTrgNpaSum, axis=(1,2))
           
        if LogToConsole:
            print(f'\nSrcMaximaBySlice   = {SrcMaximaBySlice}')
            print(f'\nResSrcMaximaBySlice = {ResSrcMaximaBySlice}')
            print(f'\nTrgMaximaBySlice    = {TrgMaximaBySlice}')
            print(f'\nNewTrgMaximaBySlice = {NewTrgMaximaBySlice}')
        
        
        #AllMaxima = [max(SrcMaximaBySlice), max(ResSrcMaximaBySlice),
        #             max(TrgMaximaBySlice), max(NewTrgMaximaBySlice)]
            
        
        # Determine which slices have at least one non-zero labelmap across all
        # labelmaps:
        SumOfMaximaBySlice = []
        
        for s in range(Nslices):
            count = 0
            
            if s < SrcNslices:
                count += SrcMaximaBySlice[s]
                
            if s < ResSrcNslices:
                count += ResSrcMaximaBySlice[s]
                
            if s < TrgNslices:
                count += TrgMaximaBySlice[s]
                
            if s < NewTrgNslices:
                count += NewTrgMaximaBySlice[s]
                
            SumOfMaximaBySlice.append(count)
        
        
        #print(f'\nSumOfMaximaBySlice = {SumOfMaximaBySlice}')
        
    
    # Prepare the figure.
    
    # Set the number of subplot rows and columns:
    Ncols = 4
    
    if 'Labelmaps' in ImageType:
        Nrows = np.count_nonzero(np.array(SumOfMaximaBySlice))
    else:
        Nrows = Nslices
    
    
    # Create a figure with two subplots and the specified size:
    #fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows))
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows), dpi=dpi)
        
    
    n = 1 # initialised sub-plot number
    
    # Define the colourmaps to use for each ROI:
    #cmaps = [plt.cm.spring, plt.cm.cool, plt.cm.autumn, plt.cm.winter, 
    #         plt.cm.hsv, plt.cm.jet]
    
    # Loop through each slice: 
    for s in range(Nslices):
        # Only proceed if ImageType is not 'Labelmaps' or (if it is) at least 
        # one of the labelmaps is non-zero for this slice number:
        if not 'Labelmaps' in ImageType or SumOfMaximaBySlice[s]:
            # Plot the Original Source slice:
            if s < SrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                # Loop through all ROIs:
                #maxima = []
                #for r in range(len(SrcNpas)):
                #    im = ax.imshow(SrcNpas[r][s], cmap=plt.cm.Greys_r)
                #    #im = ax.imshow(SrcNpas[r][s], cmap=cmaps[r])
                #    #maxima.append(np.amax(SrcNpas[r][s]))
                
                im = ax.imshow(SrcNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if SrcMaximaBySlice[s] < 1 else SrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Original Source slice {s}')
                
            n += 1 # increment sub-plot number
                
            # Plot the Resampled Source slice:
            if s < ResSrcNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                #maxima = []
                #npa = np.zeros()
                #for r in range(len(ResSrcNpas)):
                #    im = ax.imshow(ResSrcNpas[r][s], cmap=plt.cm.Greys_r)
                #    #im = ax.imshow(ResSrcNpas[r][s], cmap=cmaps[r])
                #    #maxima.append(np.amax(ResSrcNpas[r][s]))
                
                im = ax.imshow(ResSrcNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if ResSrcMaximaBySlice[s] < 1 else ResSrcMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Resampled Source slice {s}')
                
            n += 1 # increment sub-plot number
            
            # Plot the Target slice:
            if s < TrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                #maxima = []
                #for r in range(len(TrgNpas)):
                #    im = ax.imshow(TrgNpas[r][s], cmap=plt.cm.Greys_r)
                #    #im = ax.imshow(TrgNpas[r][s], cmap=cmaps[r])
                #    #maxima.append(np.amax(TrgNpas[r][s]))
                
                im = ax.imshow(TrgNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if TrgMaximaBySlice[s] < 1 else TrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'Target slice {s}')
                
            n += 1 # increment sub-plot number
            
            # Plot the New Target slice:
            if s < NewTrgNslices:
                ax = plt.subplot(Nrows, Ncols, n)
                
                #for r in range(len(NewTrgNpas)):
                #    im = ax.imshow(NewTrgNpas[r][s], cmap=plt.cm.Greys_r)
                
                im = ax.imshow(NewTrgNpaSum[s], cmap=plt.cm.Greys_r)
                
                if 'Labelmaps' in ImageType:
                    m = 1 if NewTrgMaximaBySlice[s] < 1 else NewTrgMaximaBySlice[s]
                    cbar = fig.colorbar(im, ax=ax)
                    cbar.mappable.set_clim(0, m)
                    
                ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                ax.set_title(f'New Target slice {s}')
                
            n += 1 # increment sub-plot number
        
    
    if ExportPlot:
        CurrentDate = time.strftime("%Y%m%d", time.gmtime())
        CurrentTime = time.strftime("%H%M%S", time.gmtime())
        CurrentDateTime = CurrentDate + '_' + CurrentTime
                      
        ExportFname = CurrentDateTime + '_Orig' + SrcImLabel + '_Resamp' \
                      + SrcImLabel + '_' + TrgImLabel + '_New' + TrgImLabel \
                      + '_' + ImageType + '.jpg'
        
        ExportFpath = os.path.join(ExportDir, ExportFname)
        
        plt.savefig(ExportFpath, bbox_inches='tight')
        
        print('\nPlot exported to:\n\n', ExportFpath)
        
    return







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