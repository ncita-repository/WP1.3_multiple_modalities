# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:56:06 2020

@author: ctorti

Note: 
    GetContourPtsForDicom() was adapted from James Petts' Java code to convert 
    from 3D coordinates in RTSTRUCTs to 2D in-image coordinates (copied from 
    and provided by James Darcy via Slack on 04/02/2020).
    
    For this function the derivations were made from scratch.
"""


"""
Function:
    GetContourPtsForDicom_v2()
    
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



"""
r = s + i.X + j.Y + k.Z

"""

def GetContourPtsForDicom_v2(DicomFpath, RoiFpath):
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
    
    
    # Define S as the origin along x, y and z, and X, Y and Z as the direction 
    # cosines along x, y, and z:
    S = [float(item) for item in Position]
    
    X = [float(item) for item in Orientation[0:3]]
    Y = [float(item) for item in Orientation[3:6]] 
    Z = [float(item) for item in Orientation[6:9]]
    
    # The indeces of the largest direction cosines along x, y and z:
    ind_i = X.index(max(X)) # typically = 0
    ind_j = Y.index(max(Y)) # typically = 1
    ind_k = Z.index(max(Z)) # typically = 2

    # The pixel spacings:
    di = float(PixSpacing[0])
    dj = float(PixSpacing[1])
    dk = float(PixSpacing[2])
    
    
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
            R = [float(item) for item in ContourSequence.ContourData]
            #NoContourPts = int(ContourSequence.NumberOfContourPoints)
        
            # Store the contour points for each contour ("_c"):
            """ Note: There will be multiple contours if len(inds) > 1. """
            PtsArr_c_PCS = [] # array of array
            Pts_c_PCS = [] # flat array
            PtsArr_c_ICS = [] # array of array
            Pts_c_ICS = [] # flat array
            
            # Iterate for all contour points in threes:
            for p in range(0, len(R), 3):
                #if CoordSys=='PCS':
                # Keep contour points in Patient Coordinate System:
                Pts_c_PCS.append([R[p],
                                  R[p+1],
                                  R[p+2]])
    
                PtsArr_c_PCS.append([R[p],
                                     R[p+1],
                                     R[p+2]])
        
                #if CoordSys=='ICS':
                # Convert contour points from Patient to Image Coordinate
                # System:
                
                # Define vector v as the the patient position subtracted 
                # from the contour coordinates:
                #V = R[p]   - S[0], \
                #    R[p+1] - S[1], \
                #    R[p+2] - S[2]
                    
                Vx = R[p]   - S[0]
                Vy = R[p+1] - S[1]
                Vz = R[p+2] - S[2]
                    
                # Define some helpful simplifying expressions:
                a = 1 - X[ind_k]*Z[ind_j] / ( X[ind_j]*Z[ind_k] )
                
                b = Y[ind_k]*Z[ind_j] / ( Y[ind_j]*Z[ind_k] ) - 1
                
                c = Vz*Z[ind_j]/Z[ind_k] - Vy
                
                d = 1 - Y[ind_k]*Z[ind_i]*Z[ind_k]/Y[ind_i]
                
                e = 1 - X[ind_k]*Z[ind_i] / ( X[ind_i]*Z[ind_k] )
                
                f = 1 - Y[ind_k]*Z[ind_i]/Y[ind_i]
                
                # Convert from patient coord system to image coord system:   
                i = ( Vx - (c*d/b)*Y[ind_i]/ ( Y[ind_j]*Z[ind_k] ) - Vz*Z[ind_i]/Z[ind_k] ) / ( di*( e*X[ind_i] + (a*f/b)*Y[ind_i]/Y[ind_j] ) )
                
                j = (1 / ( dj*Y[ind_j] ) ) * ( (c/b) + (a/b)*di*i )
                
                k = (1 / ( dk*Z[ind_k] ) ) * ( Vz - di*X[ind_k]*i - dj*Y[ind_k]*j )
    
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
    x = []
    y = []
    z = []  
    
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
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    PtsDict_PCS = {'x':x, 'y':y, 'z':z}  
    
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []
    
    for Point in Pts_ICS:
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    PtsDict_ICS = {'x':x, 'y':y, 'z':z}      
        
    #return ContourPts, ContourPtsArr , ContourPtsDict
    return Pts_PCS, PtsArr_PCS, PtsDict_PCS, Pts_ICS, PtsArr_ICS, PtsDict_ICS