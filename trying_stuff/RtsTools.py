# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:55:17 2020

@author: ctorti
"""




# Import packages and functions:
import os
import time
import copy
from ImageTools import GetImageAttributes




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
    NewRtsRoi = copy.deepcopy(RtsRoi)
        
    # The last item in Contour Image Sequence:
    last = copy.deepcopy(RtsRoi.ReferencedFrameOfReferenceSequence[0]\
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
    last = copy.deepcopy(RtsRoi.ROIContourSequence[0]\
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
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .NumberOfContourPoints = copy.deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                                                      .ContourSequence[FromOrigSeqNum]\
                                                                      .NumberOfContourPoints)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .ContourData = copy.deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                                            .ContourSequence[FromOrigSeqNum]\
                                                            .ContourData)
        
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
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .NumberOfContourPoints = copy.deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                                                      .ContourSequence[ind]\
                                                                      .NumberOfContourPoints)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            NewRtsRoi.ROIContourSequence[0]\
                     .ContourSequence[i]\
                     .ContourData = copy.deepcopy(OrigRtsRoi.ROIContourSequence[0]\
                                                            .ContourSequence[ind]\
                                                            .ContourData)
                                          
    return NewRtsRoi






def GetExtremePoints(Image):

    pts = [Image.TransformIndexToPhysicalPoint((0, 0, 0)), 
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), 0, 0)),
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), Image.GetHeight(), 0)),
           Image.TransformIndexToPhysicalPoint((0, Image.GetHeight(), 0)),
           Image.TransformIndexToPhysicalPoint((0, 0, Image.GetDepth())), 
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), 0, Image.GetDepth())),
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), Image.GetHeight(), Image.GetDepth())),
           Image.TransformIndexToPhysicalPoint((0, Image.GetHeight(), Image.GetDepth()))]
    
    return pts






def ConvertPcsPointToIcs(PointPCS, Origin, Directions, Spacings):
    """
    Convert a point from the Patient Coordinate System (PCS) to the Image 
    Coordinate System (ICS).
    
    Input:
        PointPCS   - (List of floats) Points in the Patient Coordinate System, 
                     e.g. [x, y, z]
        
        Origin     - (List of floats) The 3D image origin 
                     (= ImagePositionPatient of the first slice in the DICOM 
                     series), e.g. [x0, y0, z0]
        
        Directions - (List of floats) The direction cosine along x (rows), 
                     y (columns) and z (slices) (= the cross product of the 
                     x and y direction cosines from ImageOrientationPatient 
                     appended to ImageOrientationPatient),
                     e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz]
        
        Spacings   - (List of floats) The pixel spacings along x, y and z 
                     (= SliceThickness appended to PixelSpacing), 
                     e.g. [di, dj, dk]
    
        
    Returns:
        PointICS   - (List of floats) PointPCS converted to the ICS,
                     e.g. [x, y, z]
        
    
    Equations:
    
    P = S + (di*i).X + (dj*j).Y + (dk*k).Z
    
    where P = (Px, Py, Pz) is the point in the Patient Coordinate System
    
          S = (Sx, Sy, Sz) is the origin in the Patient Coordinate System (i.e.
          the ImagePositionPatient)
          
          X = (Xx, Xy, Xz) is the direction cosine vector along rows (x) 
          
          Y = (Yx, Yy, Yz) is the direction cosine vector along columns (y)
          
          Z = (Zx, Zy, Zz) is the direction cosine vector along slices (z)
          
          di, dj and dk are the pixel spacings along x, y and z
          
          i, j, and k are the row (x), column (y) and slice (z) indices
          
          
    Solve for i, j and k, i.e.:
        
        Vx = Xx*di*i + Yx*dj*j + Zx*dk*k
        
        Vy = Xy*di*i + Yy*dj*j + Zy*dk*k
        
        Vz = Xz*di*i + Yz*dk*k + Zz*dk*k
        
    where Vx = Px - Sx
    
          Vy = Py - Sy
          
          Vz = Pz - Sz
          
    Solving for i, j and k results in the expressions:
        
        i = (1/di)*d/e
            
        j = (1/dj)*( (a/b)*di*i + c/b )
        
        k = (1/dk)*(1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
        or
        k = (1/dk)*(1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
        
        (either expression for k should be fine)
        
    In the case where X = [1, 0, 0]
                      Y = [0, 1, 0]
                      Z = [0, 0, 1]
                      
    The expressions simplify to:
        
        a = 0
        
        b = 1
        
        c = Vy
        
        d = Vx
        
        e = 1
        
        i = Vx/di
        
        j = Vy/dj
        
        k = Vz/dk
        
    """
    
    # Define S, X, Y and Z:
    S = Origin # the origin
    
    X = Directions[0:3] # the row (x) direction cosine vector
    Y = Directions[3:6] # the column (y) direction cosine vector
    Z = Directions[6:] # the slice (z) direction cosine vector
    
    # The indices of the largest direction cosines along rows, columns and 
    # slices:
    ind_i = X.index(max(X)) # typically = 0
    ind_j = Y.index(max(Y)) # typically = 1
    ind_k = Z.index(max(Z)) # typically = 2

    # The pixel spacings:
    di = Spacings[0]
    dj = Spacings[1]
    dk = Spacings[2] 
    
    # Simplifying expressions:
    Xx = X[ind_i]
    Xy = X[ind_j]
    Xz = X[ind_k]
    
    Yx = Y[ind_i]
    Yy = Y[ind_j]
    Yz = Y[ind_k]
    
    Zx = Z[ind_i]
    Zy = Z[ind_j]
    Zz = Z[ind_k]
    
    # Define simplifying expressions:
    Vx = PointPCS[0] - S[0]
    Vy = PointPCS[1] - S[1]
    Vz = PointPCS[2] - S[2]
    
    a = Xz*Zy/Zz - Xy
    
    b = Yy - Yz*Zy/Zz
    
    c = Vy - Vz*Zy/Zz
    
    d = Vx - (c/b)*Yx - (Zx/Zz)*(Vz - (c/b)*Yz)
    
    e = Xx + (a/b)*Yx - (Zx/Zz)*(Xz + (a/b)*Yz)
    
    # Solve for i, j and k:
    i = (1/di)*d/e
    
    j = (1/dj)*( (a/b)*di*i + c/b )
    
    k = (1/dk)*(1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
    #k = (1/dk)*(1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
    
    PointICS = [i, j, k]
    
    return PointICS






def ConvertContourDataToIcsPoints(ContourData, DicomDir):
    """
    Convert contour data from a flat list of points in the Patient Coordinate
    System (PCS) to a list of [x, y, z] points in the Image Coordinate System. 
    
    Inputs:
        ContourData - (List of floats) Flat list of points in the PCS,
                      e.g. [x0, y0, z0, x1, y1, z1, ...]
        
        DicomDir    - (String) Directory containing DICOMs
        
        
    Returns:
        PointsICS   - (List of floats) List of a list of [x, y, z] coordinates
                      of each point in ContourData converted to the ICS,
                      e.g. [[x0, y0, z0], [x1, y1, z1], ...]
    """
    
    
    # Get the image attributes:
    #Origin, Directions, Spacings, Dimensions = GetImageAttributes(DicomDir)
    IPPs, Directions, Spacings, Dimensions = GetImageAttributes(DicomDir)
    
    # Initialise the array of contour points:
    PointsICS = []
    
    # Iterate for all contour points in threes:
    for p in range(0, len(ContourData), 3):
        # The point in Patient Coordinate System:
        pointPCS = [ContourData[p], ContourData[p+1], ContourData[p+2]]
        
        # Convert ptPCS to Image Coordinate System:
        #pointICS = ConvertPcsPointToIcs(pointPCS, Origin, Directions, Spacings)
        pointICS = ConvertPcsPointToIcs(pointPCS, IPPs[0], Directions, Spacings)
        
        PointsICS.append(pointICS)
        
    return PointsICS





def UnpackPoints(Points):
    """
    Unpack list of a list of [x, y, z] coordinates into lists of x, y and z
    coordinates.
    
    Input:
        Points - (List of floats) List of a list of [x, y, z] coordinates,
                 e.g. [[x0, y0, z0], [x1, y1, z1], ...]
        
    Returns:
        X      - (List of floats) List of all x coordinates in Points,
                 e.g. [x0, x1, x2, ...]
        
        Y      - (List of floats) List of all y coordinates in Points,
                 e.g. [y0, y1, y2, ...]
        
        Z      - (List of floats) List of all z coordinates in Points,
                 e.g. [z0, z1, z2, ...]
    """
    
    # Initialise unpacked lists:
    X = []
    Y = []
    Z = []
    
    # Unpack tuple and store in X, Y, Z:
    for x, y, z in Points:
        X.append(x)
        Y.append(y)
        Z.append(z)
        
    return X, Y, Z





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