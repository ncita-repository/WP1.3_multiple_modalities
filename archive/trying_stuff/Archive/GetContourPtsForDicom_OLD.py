# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:10:58 2020

@author: ctorti

Parts of code below was adapted from James Petts' Java code to convert from 
3D coordinates in RTSTRUCTs to 2D in-image coordinates (copied from and 
provided by James Darcy via Slack on 04/02/2020).
"""


"""
Function:
    GetContourPtsForDicom()
    
Purpose:
    Get contour points belonging to a DICOM image from an ROI object.
    
    , e.g.:
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

Inputs:
    DicomFpath     - filepath to a DICOM file
    
    RoiFpath       - filepath to a ROI Collection file
    
    CoordSys       - the coordinate system of the contour points returned;
                     acceptable values are:
                        - 'PCS' for Patient Coordinate System (i.e. ContourData
                        in Roi is in PCS)
                        - 'ICS' for Image Coordinate System (accounts for the
                        Image Plane attributes)
    
Returns:
    ContourPts     - array of all contour points in Roi that belong to DICOM
    
                    e.g. [[x0, y0, z0],
                          [x1, y1, z1],
                          ...
                          ]
    
    ContourPtsArr  - array of arrays of contour points for each contour in Roi
                     that belong to DICOM
                     
                     e.g. [[[x0_c0, y0_c0, z0_c0],
                            [x1_c0, y1_c0, z1_c0],
                            ...
                            ],
                           [[x0_c1, y0_c1, z0_c1],
                            [x1_c1, y1_c1, z1_c1],
                            ...
                           ],
                           ...
                          ]
                           
                           where _cN denotes the N^th contour
    
    ContourPtsDict - dictionary of all contour points in Roi that belong to
                     DICOM
    
                     e.g. {'x':[x0, x1, ...],
                           'y':[y0, y1, ...],
                           'z':[z0, z1, ...]
                           }
    
"""


#def GetContourPtsForDicom(DicomFpath, RoiFpath, CoordSys):
def GetContourPtsForDicom(DicomFpath, RoiFpath):
    # Import packages:
    import pydicom
    
    # Read in the DICOM and Roi:
    Dicom = pydicom.read_file(DicomFpath)
    Roi = pydicom.read_file(RoiFpath)

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
    
    # Initialise the array of contour points (ContourPts) and array of array of
    # contour points (ContourPtsArr):
    Pts_PCS = []
    PtsArr_PCS = []
    Pts_ICS = []
    PtsArr_ICS = []
        
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
            PtsArr_c_PCS = [] # array of array
            Pts_c_PCS = [] # flat array
            PtsArr_c_ICS = [] # array of array
            Pts_c_ICS = [] # flat array
            
            # Iterate for all contour points in threes:
            for p in range(0, len(ContourData), 3):
                #if CoordSys=='PCS':
                # Keep contour points in Patient Coordinate System:
                Pts_c_PCS.append([ContourData[p],
                                  ContourData[p+1],
                                  ContourData[p+2]])
    
                PtsArr_c_PCS.append([ContourData[p],
                                     ContourData[p+1],
                                     ContourData[p+2]])
        
                #if CoordSys=='ICS':
                # Convert contour points from Patient to Image Coordinate
                # System:
                
                # Define vector v as the the patient position subtracted 
                # from the contour coordinates:
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
                
                """ 11/05/2020:  Add k for 3D - note this doesn't take into
                account the direction cosine along the z-axis!!!: """
                k = v[2]
    
                """ For 2D: """
                #ContourPoints_c.append([i,j])
                """ 11/05/2020:  For 3D: """
                Pts_c_ICS.append([i, j, k])
                
                PtsArr_c_ICS.append([i, j, k])


            Pts_PCS.extend(Pts_c_PCS)
            PtsArr_PCS.extend(PtsArr_c_PCS)
            Pts_ICS.extend(Pts_c_ICS)
            PtsArr_ICS.extend(PtsArr_c_ICS)
            """ Note:  PtsArr will be the same as Pts if there is
            only one contour. """
            

            
            
    # If ContourPoints is [] print warning:
    #if ContourPoints==[]:
    #    print('No contour points were found for this DICOM (none of',
    #          'the Referenced SOP Instance UIDs in the Contour Sequence',
    #          'match the \nDICOM SOP Instance UID).')
     
    
    # Create a dictionary:
    
    # Initialise an array to store the x, y and z coordinates:
    X = []
    Y = []
    Z = []  
    
    # Loop through each array of points in ContourPts and append the x, y and z
    # coordinates to X, Y and Z:
    #for ArrayOfPts in ContourPoints:
    #    if ArrayOfPts: # if there are points for this contour
    #        for Point in ArrayOfPts:
    #            X.append(Point[0])
    #            Y.append(Point[1])
    #            Z.append(Point[2])
    
    #print('\nContourPts:\n', ContourPts)
    
    # Loop through each array of points in ContourPts and append the x, y and z
    # coordinates to X, Y and Z:
    for Point in Pts_PCS:
        X.append(Point[0])
        Y.append(Point[1])
        Z.append(Point[2])
    
    PtsDict_PCS = {'x':X, 'y':Y, 'z':Z}  
    
    # Initialise an array to store the x, y and z coordinates:
    X = []
    Y = []
    Z = []
    
    for Point in Pts_ICS:
        X.append(Point[0])
        Y.append(Point[1])
        Z.append(Point[2])
    
    PtsDict_ICS = {'x':X, 'y':Y, 'z':Z}      
        
    #return ContourPts, ContourPtsArr , ContourPtsDict
    return Pts_PCS, PtsArr_PCS, PtsDict_PCS, Pts_ICS, PtsArr_ICS, PtsDict_ICS