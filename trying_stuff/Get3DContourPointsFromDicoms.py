# -*- coding: utf-8 -*-
"""
Created on Fri May 22 12:08:21 2020

@author: ctorti
Parts of code below was adapted from James Petts' Java code to convert from 
3D coordinates in RTSTRUCTs to 2D in-image coordinates (copied from and 
provided by James Darcy via Slack on 04/02/2020).
"""


"""
Function:
    Get3DContourPointsFromDicoms()
    
Purpose:
    Get array of 3D contour points from an ROI object that correspond to a
    series of DICOM objects contained in a directory, e.g.:
        [
         [[x1_1, y1_1, z1_1], [x2_1, y2_1, z2_1], [x3_1, y3_1, z3_1], ... ],
         [[x1_2, y1_2, z1_2], [x2_2, y2_2, z2_2], [x3_2, y3_2, z3_2], ... ],
         ...
         [[x1_N, y1_N, z1_N], [x2_N, y2_N, z2_N], [x3_N, y3_N, z3_N], ... ]
        ]
    
        where each sub-array contains the [x, y, z] points contained in a 
        contour, of which there are N number of contours.
        
        The contour points are converted from the patient coordinate system to
        the image coordinate system.

Input:
    DicomDir      - Directory containing DICOM objects
    
    RoiFpath      - Filepath of an ROI Collection
    
Returns:
    ContourPoints - array of contour points in Roi that belong to all DICOMs
    
"""


def Get3DContourPointsFromDicoms(DicomDir, RoiFpath):
    # Import packages:
    import pydicom
    from GetDicoms import GetDicoms
    
    # Get the DICOMs:
    DicomFpaths, Dicoms = GetDicoms(DicomDir, 'slices', Debug=False)

    # Read in the ROI Collection:
    Roi = pydicom.dcmread(RoiFpath)
    
    # Create an array of contour points for all DICOMs:
    AllContourPoints = []
    
    # Loop through each DICOM and get the contour points that relate to it:
    for Dicom in Dicoms:
        
        # Get the SOP Instance UID, Image Position Patient, Orientation and 
        # Pixel Spacing from Dicom:
        SopUid = Dicom.SOPInstanceUID
        Position = Dicom.ImagePositionPatient 
        # e.g. ['-104.92378045297', '-152.49065317178', '-9.2214860452719']
        Orientation = Dicom.ImageOrientationPatient 
        # e.g. ['0.99726487384986', '-0.0222629979573', '-0.070477871046', 
        #       '1.863426e-010', '0.95355613378866', '-0.301215370979']
        PixSpacing = Dicom.PixelSpacing # e.g. ['0.8984375', '0.8984375']
        
        # The patient coordinates:
        x0 = Position[0]
        y0 = Position[1]
        z0 = Position[2]
        
        # The direction cosines along rows and columns:
        Theta_r = Orientation[0:3]
        Theta_c = Orientation[3:6]
        
        # The indeces of the largest direction cosines along rows and columns:
        IndMaxTheta_r = Theta_r.index(max(Theta_r)) # typically = 0
        IndMaxTheta_c = Theta_c.index(max(Theta_c)) # typically = 1
    
        # The pixel spacings:
        dx = PixSpacing[0]
        dy = PixSpacing[1]
            
        
        # Get a list of the Contour Sequences in the Roi:
        ContourSequences = Roi.ROIContourSequence[0].ContourSequence
        
        # Get a list of all Referenced SOP Instance UIDs from the ContourSequences:
        RefSopUids = [sequence.ContourImageSequence[0]\
                      .ReferencedSOPInstanceUID for sequence in ContourSequences]
        
        # Create an array of (x,y) points from the flattened contour data:
        ContourPoints = []
            
        # Get the index of the Contour Sequence whose Referenced SOP Instance UID
        # matches the DICOM's SOP Instance UID:
        """ Note that a contour sequence may not exist for this Dicom, so check if
        SopUid is in RefSopUids.  Proceed only if true. """
        if SopUid in RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first match,
            # and there may be more than one!
            
            # Get the indeces of all matching RefSopUids:
            inds = [i for i, e in enumerate(RefSopUids) if e==SopUid]
            
            # Iterate for each index in inds:
            for ind in inds:
                # Get the Contour Sequence that corresponds to this ind:
                ContourSequence = ContourSequences[ind]
                
                # Get the Contour Data and number of contour points for this 
                # ContourSequence:
                ContourData = [float(item) for item in ContourSequence.ContourData]
                #NoContourPts = int(ContourSequence.NumberOfContourPoints)
            
                # Store the contour points for each contour ("_c"):
                """ Note: There will be multiple contours if len(inds) > 1. """
                ContourPoints_c = []
                
                # Iterate for all contour points in threes:
                for p in range(0, len(ContourData), 3):
                    # Define vector v as the the patient position subtracted from 
                    # the contour coordinates:
                    v = ContourData[p]   - x0, \
                        ContourData[p+1] - y0, \
                        ContourData[p+2] - z0
                    
                    # Convert from patient coord system to image coord system:   
                    i = (v[IndMaxTheta_r] - v[IndMaxTheta_c]*Theta_c[IndMaxTheta_r]/Theta_c[IndMaxTheta_c]) / \
                    (
                    Theta_r[IndMaxTheta_r]*dx *
                    (1 - (Theta_c[IndMaxTheta_r]*Theta_r[IndMaxTheta_c]) / (Theta_r[IndMaxTheta_r]*Theta_c[IndMaxTheta_c]))
                    )
                    
                    # Equation above copied below but without line breaks:
                    # i = (v[IndMaxTheta_r] - Theta_c[IndMaxTheta_r]*v[IndMaxTheta_c]/Theta_c[IndMaxTheta_c]) / ( Theta_r[IndMaxTheta_r]*dx * (1 - (Theta_c[IndMaxTheta_r]*Theta_r[IndMaxTheta_c]) / (Theta_r[IndMaxTheta_r]*Theta_c[IndMaxTheta_c])) )
                    
                    j = (v[IndMaxTheta_c] - Theta_r[IndMaxTheta_c]*i*dx) / (Theta_c[IndMaxTheta_c]*dy)
                    
                    """ 11/05/2020:  Add k for 3D: """
                    k = v[2]
        
                    """ For 2D: """
                    #ContourPoints_c.append([i,j])
                    """ 11/05/2020:  For 3D: """
                    ContourPoints_c.append([i,j,k])
    
                    
                # Append contourPts with ContourPts_c:
                #ContourPoints.append(ContourPoints_c)
                ContourPoints.extend(ContourPoints_c)
                
                
        # If ContourPoints is [] print warning:
        #if ContourPoints==[]:
        #    print('No contour points were found for this DICOM (none of',
        #          'the Referenced SOP Instance UIDs in the Contour Sequence',
        #          'match the \nDICOM SOP Instance UID).')
        
        
        # Append ContourPoints to AllContourPts:
        AllContourPoints.append(ContourPoints)
    
    return AllContourPoints
