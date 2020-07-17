# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:04:37 2020

@author: ctorti
"""





"""
Function:
    GetContourPtsForDicoms()
    
Purpose:
    Create two arrays of contour points (points organised in one of two ways)
    and a LUT.


Input:
    DicomDir        - Directory containing DICOMs 
    
    RoiFpath        - Filepath of the DICOM-RTSTRUCT ROI Collection
    
    
Returns:
    Pts - Array of points, e.g.:
                    
            [
             [x0, y0, z0],
             [x1, y1, z1],
             ...
             ]
    
    PtsBySliceAndContourAndContour - Array of array of array of points organised by 
                           contour and slice, e.g.
                
                [
                 [ <-- start of array of array of points in slice 0
                  [x0_c0_s0, y0_c0_s0, z0_c0_s0], [x1_c0_s0, y1_c0_s0, z1_c0_s0], ...], <-- array of points in contour 0 in slice 0
                  [x0_c1_s0, y0_c1_s0, z0_c1_s0], [x1_c1_s0, y1_c1_s0, z1_c1_s0], ...], <-- array of points in contour 1 in slice 0
                  ...
                  ], <-- end of array of array of contour in slice 0
                 
                 [ <-- start of array of array of points in slice 1
                  [x0_c0_s1, y0_c0_s1, z0_c0_s1], [x1_c0_s1, y1_c0_s1, z1_c0_s1], ...], <-- array of points in contour 0 in slice 1
                  [x0_c1_s1, y0_c1_s1, z0_c1_s1], [x1_c1_s1, y1_c1_s1, z1_c1_s1], ...], <-- array of points in contour 1 in slice 1
                  ...
                  ], <-- end of array of array of contour in slice 1
                 ...
                 ]
                 
    LUT - A P x 5 array (P rows for all P points, 5 columns) containing:
        
            - (col 0): the point number: integer starting from 0, 
            - (col 1): the (fixed image) slice number it came from: integer starting from 0,    (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE NEXT ROW)
            - (col 2): the (fixed image) contour number it belongs to: integer starting from 0, (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE PREVIOUS ROW)
            - (col 3): an empty list ([]) as a placeholder for the (fixed image) index of the point: array of integers (N=3),
            - (col 4): the (fixed image) point in PCS: array of floats (N=3)
            
            e.g. [0, 
                  0 
                  12, 
                  [],
                  [-19.904, -10.938, 2.551]]
            
            => The 1st point is part of the 1st contour (the 1st contour in
            absolute terms, not the 1st of possibly multiple contours in any
            given slice) that came from the 13th slice, whose coordinates in 
            the PCS are [-19,904, -10.938, 2.551].
            
       - Note that LUT will have items added to it so that it will ultimately
         be a P x 12 array (P rows for all P points, 12 columns) containing:
        
            - (col 0): the point number: integer starting from 0, 
            - (col 1): the (fixed image) slice number it came from: integer starting from 0,    (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE NEXT ROW)
            - (col 2): the (fixed image) contour number it belongs to: integer starting from 0, (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE PREVIOUS ROW)
            - (col 3): the input (fixed image) index of the point: array of integers (N=3),
            - (col 4): the input (fixed image) point in PCS: array of floats (N=3),
            - (col 5): the input (fixed image) point in ICS: array of floats (N=3),
            - (col 6): the (output / moving image) contour number the transformed point: integer starting from 0,
            - (col 7): the (output / moving image) index of the transformed point: array of intergers (N=3),
            - (col 8): the (output / moving image) transformed point in PCS: array of floats (N=3),
            - (col 9): the (output / moving image) transformed point in ICS: array of floats (N=3),
            - (col 10): the deformation in PCS: array of floats (N=3),
            - (col 11): the deformation in ICS: array of floats (N=3)

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
    
    
    
    # Initialise the array of contour points (Pts) and array of array of 
    # contour points organised slice-by-slice along rows and contour-by-contour
    # (PtsBySliceAndContourAndContour):
    Pts = []
    PtsBySliceAndContour = []

    # Initialise the LUT:
    LUT = []
    
    PtNo = 0 # iterate the point number
    CntNo = 0 # iterate the contour number
    
    # Loop through each DICOM:
    for SliceNo in range(len(DicomFpaths)):
        # Read in the DICOM:
        Dicom = pydicom.read_file(DicomFpaths[SliceNo])
    
        # Get the SOP Instance UID:
        SopUid = Dicom.SOPInstanceUID
    
        # Initialise the array of contour points for this slice:
        PtsThisSlice = []

            
        # Get the index of the Contour Sequence whose Referenced SOP Instance 
        # UID matches the DICOM's SOP Instance UID:
        """ Note that a contour sequence may not exist for this Dicom, so check
        if SopUid is in RefSopUids.  Proceed only if true. """
        if SopUid in RefSopUids:
            # Get the indeces of all matching RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first 
            # match, and there may be more than one!
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

                # Initialise the array of contour points for this contour:
                PtsThisContour = []
                
                # Iterate for all contour points in threes:
                for p in range(0, len(ContourData), 3):
                    Pt = [ContourData[p], ContourData[p+1], ContourData[p+2]]
                    
                    Pts.append(Pt)
                    
                    PtsThisContour.append(Pt)
        
                    #LUT.append([PtNo, CntNo, SliceNo, [], Pt])
                    #LUT.append([PtNo, SliceNo, CntNo, [], Pt]) # Changed on 08/07/2020
                    LUT.append([PtNo, SliceNo, CntNo, 1, [], Pt]) # Changed on 17/07/2020
                    """ 
                    Notes:
                        
                    1. The 1 in the 3rd position signals that the contour 
                    point is original (not interpolated, for example).
                    
                    2. The empty array in the 4th position will be replaced 
                    with the point indeces later, as will the subsequent items 
                    that will finalise the LUT. 
                    """
                    
                    PtNo = PtNo + 1
                    
                    
                CntNo = CntNo + 1

                PtsThisSlice.append(PtsThisContour)
        
        
            PtsBySliceAndContour.append(PtsThisSlice)
                
        else: 
            PtsBySliceAndContour.append([])
            
        
    
      
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
    #return Pts, PtsBySliceAndContour, PtsDict, LUT
    
    return Pts, PtsBySliceAndContour, LUT