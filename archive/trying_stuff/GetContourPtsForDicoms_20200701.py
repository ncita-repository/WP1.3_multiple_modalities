# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:04:37 2020

@author: ctorti
"""


"""
Details of LUT:

    LUT is a P x 12 array (P rows for all P points, 11 columns) containing:
        
    - the point number: integer starting from 0, 
    - the (fixed image) contour number it belongs to: integer starting from 0,
    - the (fixed image) slice number it came from: integer starting from 0,
    - the (fixed image) index of the point: array of intergers (N=3),
    - the (fixed image) point in PCS: array of floats (N=3),
    - the (fixed image) point in ICS: array of floats (N=3),
    - the (moving image) contour number it belongs to: integer starting from 0,
    - the (moving image) index of the transformed point: array of intergers (N=3),
    - the (moving image) transformed point in PCS: array of floats (N=3),
    - the (moving image) transformed point in ICS: array of floats (N=3),
    - the deformation in PCS: array of floats (N=3),
    - the deformation in ICS: array of floats (N=3)
    
    
    e.g. [125, 
          2, 
          14, 
          [-19.904, -10.938, 2.551],
          [89.940, 146.288, 12.000],
          ...
          ]
    
    => The 126th point is part of the 3rd contour (not the 3rd contour in the 
    15th slice but the total running contour number) that came from the 15th 
    slice.
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
    
    
    
    # Initialise the array of contour points (Pts) and array of array 
    # of contour points organised on a slice-by-slice basis (PtsBySlice):
    Pts = []
    PtsBySlice = []

    # Initialise the LUT:
    LUT = []
    
    PtNo = 0 # iterate the point number
    CntNo = 0 # iterate the contour number
    
    # Loop through each DICOM filepath and get the array of contour points 
    # (ContourPts), array of array of contour points (ContourPtsArr) and 
    # dictionary of contour points (ContourPtsDict) that relate to it:
    for SliceNo in range(len(DicomFpaths)):
        # Read in the DICOM:
        Dicom = pydicom.read_file(DicomFpaths[SliceNo])
    
        # Get the SOP Instance UID:
        SopUid = Dicom.SOPInstanceUID
    
        # Initialise the array of contour points (Pts_s) and array of array
        # of contour points (PtsArr_s) for this slice:
        #Pts_s = [] # flat array of points in this slice
        #PtsArr_s = [] # array of array of points in this slice for each contour in this slice

            
        # Get the index of the Contour Sequence whose Referenced SOP Instance UID
        # matches the DICOM's SOP Instance UID:
        """ Note that a contour sequence may not exist for this Dicom, so check if
        SopUid is in RefSopUids.  Proceed only if true. """
        if SopUid in RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first match,
            # and there may be more than one!
            
            # Initialise the array of contour points (Pts_s) and array of 
            # array of contour points (PtsBySlice_s) for this slice:
            PtsThisSlice = []
            #Pts = [] # flat array
            Pts_s = [] # flat array
            #PtsArr = [] # array of array
            PtsBySlice_s = [] # array of array
                
            
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
                # Subtract 1 from ContourNumber so that the counting starts 
                # from 0 as with all other variables:
                #ContourNo = int(ContourSequence.ContourNumber) - 1
            
                # Initialise the array of contour points (Pts_s) and array of 
                # array of contour points (PtsArr_s) for this slice:
                """ Note: There will be multiple contours if len(inds) > 1. """
                #PtsThisSlice = []
                ##Pts = [] # flat array
                #Pts_s = [] # flat array
                ##PtsArr = [] # array of array
                #PtsArr_s = [] # array of array

                
                # Iterate for all contour points in threes:
                for p in range(0, len(ContourData), 3):
                    Pt = [ContourData[p], ContourData[p+1], ContourData[p+2]]
                    
                    PtsThisSlice.append(Pt)
        
                    #LUT.append([PtNo, SliceNo, ContourNo, int(p/3)])
                    LUT.append([PtNo, CntNo, SliceNo, [], Pt])
                    """ Note that the empty array in the 3rd position will be
                    replaced with the point indeces later, as will the 
                    subsequent items that will finalise the LUT. """
                    
                    PtNo = PtNo + 1
                    
                    
                CntNo = CntNo + 1

    
                #Pts.extend(PtsThisSlice)
                Pts_s.extend(PtsThisSlice)
                #PtsArr.extend(PtsThisSlice)
                #PtsArr.append(PtsThisSlice) # 29/06
                PtsBySlice_s.append(PtsThisSlice) 
                """ Note:  PtsBySlice_s will be the same as Pts_s if there is
                only one contour. """
        
        
            #AllPts.extend(Pts)
            Pts.extend(Pts_s)
            #AllPtsArr.append(PtsArr)
            PtsBySlice.append(PtsBySlice_s)
                
        else: 
            Pts.extend([])
            PtsBySlice.append([])
        
    
      
    ## Create a dictionary:
    #
    ## Initialise an array to store the x, y and z coordinates:
    #x = []
    #y = []
    #z = []  
    #
    ## Loop through each array of points in AllPts and append the x, y and z
    ## coordinates to X, Y and Z:
    #for Point in Pts:
    #    x.append(Point[0])
    #    y.append(Point[1])
    #    z.append(Point[2])
    #
    #PtsDict = {'x':x, 'y':y, 'z':z}      
    #    
    #return Pts, PtsBySlice, PtsDict, LUT
    
    return Pts, PtsBySlice, LUT