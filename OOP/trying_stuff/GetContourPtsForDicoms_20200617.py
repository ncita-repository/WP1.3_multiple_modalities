# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:04:37 2020

@author: ctorti
"""

""" 
Since use of pydicom to get the ImageOrientationPatient only returns the
direction cosines along x and y, I've had to use SimpleElastix (SimpleITK)
to read in the 3D volume to get all direction cosines.
"""


def GetContourPtsForDicoms(DicomDir, RoiFpath):
    # Import packages:
    import pydicom
    import SimpleITK as sitk
    from GetDicomFpaths import GetDicomFpaths

    
    # Get the DICOM filepaths:
    DicomFpaths = GetDicomFpaths(DirPath=DicomDir, 
                                 SortMethod='slices', 
                                 Debug=False)

    # Read in the 3D stack of DICOMs:
    sitkReader = sitk.ImageSeriesReader()
    sitkNames = sitkReader.GetGDCMSeriesFileNames(DicomDir)
    sitkReader.SetFileNames(sitkNames)
    sitkIm = sitkReader.Execute()
    sitkIm = sitk.Cast(sitkIm, sitk.sitkFloat32)
    
    # Get the Origin (Image Position Patient), Direction (Image Orientation
    # Patient) and Spacing from sitkIm:
    Origin = sitkIm.GetOrigin() 
    # e.g. (-104.92378234863281, -152.4906463623047, -9.22148609161377)
    Direction = sitkIm.GetDirection() 
    # e.g. (0.9972648738498597, 1.8634260454120088e-10, 0.07391056342109259,
    #       -0.022262997957353512, 0.9535561337886633, 0.3003915089278791,
    #       -0.07047787104598309, -0.301215370979001, 0.9509480374756598)
    Spacing = sitkIm.GetSpacing() # e.g. (0.8984375, 0.8984375, 5.0)


    # Define S, X, Y and Z:
    S = Origin # the origin
    
    X = Direction[0:3] # the direction cosine along rows (X)
    Y = Direction[3:6] # the direction cosine along columns (Y)
    Z = Direction[6:9] # the direction cosine along slices (Z)
    
    # The indeces of the largest direction cosines along rows, columns and 
    # slices:
    ind_i = X.index(max(X)) # typically = 0
    ind_j = Y.index(max(Y)) # typically = 1
    ind_k = Z.index(max(Z)) # typically = 2

    # The pixel spacings:
    di = Spacing[0]
    dj = Spacing[1]
    dk = Spacing[2] 
    
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
    
    
    # Read in the Roi:
    Roi = pydicom.read_file(RoiFpath)
    
    # Get a list of the Contour Sequences in the Roi:
    ContourSequences = Roi.ROIContourSequence[0].ContourSequence
    
    # Get a list of all Referenced SOP Instance UIDs from the ContourSequences:
    RefSopUids = [sequence.ContourImageSequence[0]\
                  .ReferencedSOPInstanceUID for sequence in ContourSequences]
    
    
    
    # Initialise the array of contour points (AllContourPts) and array of array 
    # of contour points (AllContourPtsArr):
    AllPts_PCS = []
    AllPtsArr_PCS = []
    AllPts_ICS = []
    AllPtsArr_ICS = []

    
    # Loop through each DICOM filepath and get the array of contour points 
    # (ContourPts), array of array of contour points (ContourPtsArr) and 
    # dictionary of contour points (ContourPtsDict) that relate to it:
    for DicomFpath in DicomFpaths:
        # Read in the DICOM:
        Dicom = pydicom.read_file(DicomFpath)
    
        # Get the SOP Instance UID:
        SopUid = Dicom.SOPInstanceUID
    
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
                    #a = Yz*Zy/(Yy*Zz) - 1
                    
                    #b = Vz*Zy/(Vy*Zz) - 1
                    
                    #c = 1 - Xz*Zy/(Xy*Zz)
                    
                    #d = 1 - b*Vy*Yz/(a*Vz*Yy)
                    
                    #e = 1 + c*Xy*Yz/(a*Xz*Yy)
                    
                    #f = 1 + c*Xy*Yx/(a*Xx*Yy) - e*Xz*Zx/(Xx*Zz)
                    
                    #g = 1 - b*Vy*Yx/(a*Vx*Yy) - d*Vz*Zx/(Vx*Zz)
                    
                    a = (Xy - Xz/Zz) / (Yz/Zz - Yy)
                    
                    b = (Vz/Zz - Vy) / (Yz/Zz - Yy)
                    
                    c = Vx - b*Yx - (Zx/Zz)*(Vz - b*Yz)
                    
                    d = Xx + a*Yx - (Zx/Zz)*(Xz + a*Yz)
                    
                    
                    # Convert from patient coord system to image coord system:   
                    #i = g*Vx / (f*Xx*di)
                    
                    #j = 1/(a*dj*Yy) * (b*Vy + c*Xy*di*i)
                    
                    #k = ( 1/(dk*Zz) ) * (Vz - Xz*di*i - Yz*dj*j)
                    #k = ( 1/(dk*Zz) ) * (d*Vz - e*Xz*di*i) # similar results to above
                    
                    i = (1/di)*c/d
                    
                    j = (1/dj)*(a*di*i + b)
                    
                    k = (1/dk)*(1/Zz)*(Vz - b*Yz - (Xz + a*Yz)*di*i)
                    
        
                    Pts_c_ICS.append([i, j, k])
                    
                    PtsArr_c_ICS.append([i, j, k])
    
    
                Pts_PCS.extend(Pts_c_PCS)
                PtsArr_PCS.extend(PtsArr_c_PCS)
                Pts_ICS.extend(Pts_c_ICS)
                PtsArr_ICS.extend(PtsArr_c_ICS)
                """ Note:  PtsArr will be the same as Pts if there is
                only one contour. """
            
            
        
        
        AllPts_PCS.extend(Pts_PCS)
        AllPtsArr_PCS.append(PtsArr_PCS)
        
        AllPts_ICS.extend(Pts_ICS)
        AllPtsArr_ICS.append(PtsArr_ICS)
        
    
        
    # Create a dictionary:
    
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []  
    
    # Loop through each array of points in AllContourPts and append the x, y 
    # and z coordinates to X, Y and Z:
    for Point in AllPts_PCS:
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    AllPtsDict_PCS = {'x':x, 'y':y, 'z':z}      
    
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []  
    
    # Loop through each array of points in AllContourPts and append the x, y 
    # and z coordinates to X, Y and Z:
    for Point in AllPts_ICS:
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    AllPtsDict_ICS = {'x':x, 'y':y, 'z':z} 
        
        
    return AllPts_PCS, AllPtsArr_PCS, AllPtsDict_PCS, AllPts_ICS, AllPtsArr_ICS, AllPtsDict_ICS