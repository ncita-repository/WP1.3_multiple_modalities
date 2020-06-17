# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:04:37 2020

@author: ctorti
"""



def GetContourPtsForDicoms(DicomDir, RoiFpath):
    # Import packages:
    import pydicom
    from GetDicomFpaths import GetDicomFpaths

    
    # Get the DICOM filepaths:
    DicomFpaths = GetDicomFpaths(DirPath=DicomDir, 
                                 SortMethod='slices', 
                                 Debug=False)
    
    # Read in the Roi:
    Roi = pydicom.read_file(RoiFpath)
    
    # Get a list of the Contour Sequences in the Roi:
    ContourSequences = Roi.ROIContourSequence[0].ContourSequence
    
    # Get a list of all Referenced SOP Instance UIDs from the ContourSequences:
    RefSopUids = [sequence.ContourImageSequence[0]\
                  .ReferencedSOPInstanceUID for sequence in ContourSequences]
    
    
    
    # Initialise the array of contour points (AllPts) and array of array 
    # of contour points (AllPtsArr):
    AllPts = []
    AllPtsArr = []

    
    # Loop through each DICOM filepath and get the array of contour points 
    # (ContourPts), array of array of contour points (ContourPtsArr) and 
    # dictionary of contour points (ContourPtsDict) that relate to it:
    for DicomFpath in DicomFpaths:
        # Read in the DICOM:
        Dicom = pydicom.read_file(DicomFpath)
    
        # Get the SOP Instance UID:
        SopUid = Dicom.SOPInstanceUID
    
        # Initialise the array of contour points (Pts_PCS) and array of array
        # of contour points (Pts_PCS_Arr):
        Pts = []
        PtsArr = []

            
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
                PtsArr_c = [] # array of array
                Pts_c = [] # flat array

                
                # Iterate for all contour points in threes:
                for p in range(0, len(R), 3):
                    #if CoordSys=='PCS':
                    # Keep contour points in Patient Coordinate System:
                    Pts_c.append([R[p],
                                  R[p+1],
                                  R[p+2]])
        
                    PtsArr_c.append([R[p],
                                     R[p+1],
                                     R[p+2]])
    
    
                Pts.extend(Pts_c)
                PtsArr.extend(PtsArr_c)
                """ Note:  PtsArr will be the same as Pts if there is
                only one contour. """
            
            
        
        
        AllPts.extend(Pts)
        AllPtsArr.append(PtsArr)
        
    
        
    # Create a dictionary:
    
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []  
    
    # Loop through each array of points in AllPts and append the x, y and z
    # coordinates to X, Y and Z:
    for Point in AllPts:
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    AllPtsDict = {'x':x, 'y':y, 'z':z}      
        
        
    return AllPts, AllPtsArr, AllPtsDict