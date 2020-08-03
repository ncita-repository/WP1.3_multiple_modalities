# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 15:33:41 2020

@author: ctorti
"""


"""
Function:
    GetContourPoints()
    
Purpose:
    Get array of contour points from an ROI object that correspond to a given
    DICOM object, e.g.:
        [
         [[x1_1, y1_1], [x2_1, y2_1], [x3_1, y3_1], ... ],
         [[x1_2, y1_2], [x2_2, y2_2], [x3_2, y3_2], ... ],
         ...
         [[x1_N, y1_N], [x2_N, y2_N], [x3_N, y3_N], ... ]
        ]
    
        where each sub-array contains the [x, y] points contained in a contour,
        of which there are N number of contours.
        
        The contour points are converted from the patient coordinate system to
        the image coordinate system.

Input:
    Dicom         - a DICOM object
    
    Roi           - a ROI object
    
Returns:
    ContourPoints - array of contour points in Roi that belong to Dicom
    
"""


def GetContourPoints(Dicom, Roi):
    
    #print('\nRoi:\n', Roi)
    
    # Get the SOP Instance UID, Image Position Patient, Orientation and 
    # Pixel Spacing from Dicom:
    SopUid = Dicom.SOPInstanceUID
    Position = Dicom.ImagePositionPatient
    Orientation = Dicom.ImageOrientationPatient
    PixSpacing = Dicom.PixelSpacing
    
    # The patient coordinates:
    x0 = Position[0]
    y0 = Position[1]
    z0 = Position[2]
    
    # The direction cosines along rows and columns:
    Theta_r = Orientation[0:3]
    Theta_c = Orientation[3:6]
    
    # The indeces of the largest direction cosines along rows and columns:
    MaxTheta_r = Theta_r.index(max(Theta_r))
    MaxTheta_c = Theta_c.index(max(Theta_c))

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
        
            # Store the contour points for each contour (more than one if 
            # inds > 1):
            ContourPoints_c = []
            
            # Iterate for all contour points in threes:
            for p in range(0, len(ContourData), 3):
                # Define vector v as the the patient position subtracted from 
                # the contour coordinates:
                v = ContourData[p] - x0, \
                    ContourData[p+1] - y0, \
                    ContourData[p+2] - z0
                
                # Convert from patient coord system to image coord system:          
                i = (v[MaxTheta_r] - Theta_c[MaxTheta_r]*v[MaxTheta_c]/Theta_c[MaxTheta_c]) / \
                (
                Theta_r[MaxTheta_r]*dx * \
                (1 - (Theta_c[MaxTheta_r]*Theta_r[MaxTheta_c]) / (Theta_r[MaxTheta_r]*Theta_c[MaxTheta_c]))
                )
                
                j = (v[MaxTheta_c] - Theta_r[MaxTheta_c]*i*dx) / (Theta_c[MaxTheta_c]*dy)
    
                ContourPoints_c.append([i,j])
                
            # Append contourPts with ContourPts_c:
            ContourPoints.append(ContourPoints_c)
            
            
    # If ContourPoints is [] print warning:
    if ContourPoints==[]:
        print('No contour points were found for this DICOM (none of', \
              'the Referenced SOP Instance UIDs in the Contour Sequence', \
              'match the \nDICOM SOP Instance UID).')
            
        
    return ContourPoints