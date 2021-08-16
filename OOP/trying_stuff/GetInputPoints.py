# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:09:22 2020

@author: ctorti
"""


"""
Function:
    GetInputPoints()
    

MODIFY THE STUFF BELOW
    
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
        
            - the point number: integer starting from 0, 
            - the input (fixed image) slice number it came from: integer starting from 0,    (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE NEXT ROW)
            - the input (fixed image) contour number it belongs to: integer starting from 0, (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE PREVIOUS ROW)
            - an empty list ([]) as a placeholder for the input (fixed image) index of the point: array of integers (N=3),
            - the input (fixed image) point in PCS: array of floats (N=3)
            
            e.g. [0, 
                  12, 
                  0, 
                  [],
                  [-19.904, -10.938, 2.551]]
            
            => The 1st point is part of the 1st contour (the 1st contour in
            absolute terms, not the 1st of possibly multiple contours in any
            given slice) that came from the 13th slice, whose coordinates in 
            the PCS are [-19,904, -10.938, 2.551].
            
       - Note that LUT will have items added to it so that it will ultimately
         be a P x 12 array (P rows for all P points, 12 columns) containing:
        
            - (col 0): the point number: integer starting from 0, 
            - (col 1): the input (fixed image) slice number it came from: integer starting from 0,    (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE NEXT ROW)
            - (col 2): the input (fixed image) contour number it belongs to: integer starting from 0, (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE PREVIOUS ROW)
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



#def GetInputPoints(DicomDir, RoiFpath, Origin, Directions, Spacings): # July 24
def GetInputPoints(FixedDicomDir, FixedRoiFpath):
    # Import packages:
    #import os
    #import SimpleITK as sitk
    #import pydicom
    #from GetDicoms import GetDicoms
    #import importlib
    #import PCStoICS
    #importlib.reload(PCStoICS)
    from PCStoICS import PCStoICS
    #import RegUtilityFuncs as ruf
    #import GetContourPtsForDicoms
    #importlib.reload(GetContourPtsForDicoms)
    from GetContourPtsForDicoms import GetContourPtsForDicoms
    from GetImageAttributes import GetImageAttributes
    
    
    # Get the contour points:    
    PtsPCS, PtsBySliceAndContourPCS, LUT = GetContourPtsForDicoms(FixedDicomDir, 
                                                                  FixedRoiFpath)
    
    P = len(PtsPCS)
    S = len(PtsBySliceAndContourPCS)
    
    # List of slice numbers:
    SliceNos = list(range(S))
    
    
    # Chose which package to use to get the Image Plane Attributes:
    #package = 'sitk'
    package = 'pydicom'

    # Get the Image Attributes for the Fixed image:
    FixOrigin, FixDirs,\
    FixSpacings, FixDims = GetImageAttributes(DicomDir=FixedDicomDir, 
                                              Package=package)


    # Convert points in Pts_PCS from PCS to ICS:
    PtsICS = PCStoICS(Pts_PCS=PtsPCS, Origin=FixOrigin, 
                      Directions=FixDirs, Spacings=FixSpacings)

    
    # Initialise list of converted points in PtsBySliceAndContourPCS from PCS 
    # to ICS:
    PtsBySliceAndContourICS = []
    
    # Initialise ContourTypeBySlice:
    ContourTypeBySlice = []
    
    # Loop through each slice of points:
    for PtsThisSlicePCS in PtsBySliceAndContourPCS:
        # If PtsThisSlicePCS is not empty:
        if PtsThisSlicePCS:
            # Initialise the converted points for this slice:
            PtsThisSliceICS = []
            
            # Loop through each contour of points:
            for PtsThisContourPCS in PtsThisSlicePCS:
                PtsThisContourICS = PCStoICS(Pts_PCS=PtsThisContourPCS, 
                                             Origin=FixOrigin, 
                                             Directions=FixDirs, 
                                             Spacings=FixSpacings)
        
                PtsThisSliceICS.append(PtsThisContourICS)
               
            PtsBySliceAndContourICS.append(PtsThisSliceICS)
            
            ContourTypeBySlice.append(1)
            
        else:
            PtsBySliceAndContourICS.append([])
            
            ContourTypeBySlice.append(0)
            
    
    
    # Convert each PCS point (last item) in LUT and append to the sub-array:
    for i in range(P):
        PtPCS = LUT[i][-1]
        
        #print('PtPCS =', PtPCS)
        
        #PtICS = PCStoICS(Pts_PCS=[PtPCS], Origin=Origin,
        #                 Directions=Directions, Spacings=Spacings)
        PtICS = PCStoICS(Pts_PCS=PtPCS, Origin=FixOrigin,
                         Directions=FixDirs, Spacings=FixSpacings)
        
        #print('PtICS =', PtICS[0])
        
        #LUT[i].append(PtICS[0])
        LUT[i].append(PtICS)
    
 
    
    # Convert the LUT to a dictionary:
    PointData = {
                 'PointNo'       : [LUT[i][0] for i in range(P)], 
                 'InSliceNo'     : [LUT[i][1] for i in range(P)], 
                 'InContourNo'   : [LUT[i][2] for i in range(P)],
                 'InContourType' : [LUT[i][3] for i in range(P)],
                 'InPointIndex'  : [LUT[i][4] for i in range(P)],
                 'InPointPCS'    : [LUT[i][5] for i in range(P)],
                 'InPointICS'    : [LUT[i][6] for i in range(P)]
                 } 
    
    # Create a dictionary to store the contour data for the Fixed image
    # (input points) arranged by slices and contours:
    FixContourData = {
                      'SliceNo'     : SliceNos, 
                      'ContourType' : ContourTypeBySlice,
                      'PointPCS'    : PtsBySliceAndContourPCS,
                      'PointICS'    : PtsBySliceAndContourICS
                      } 
    
    
    # Create a dictionary to store the contour data for both the Fixed (input) 
    # and Moving (output) image domains, arranged by slices and contours:
    ContourData = {
                   'SliceNo'        : SliceNos, 
                   'FixContourType' : ContourTypeBySlice,
                   'FixPtsPCS'      : PtsBySliceAndContourPCS#,
                   #'FixPointsICS'   : PtsBySliceAndContourICS
                   } 
    
    
          
    #print(len(SliceNos))
    #print(len(PtsBySliceAndContourPCS))     
    
    return PtsPCS, PtsBySliceAndContourPCS,\
           PtsICS, PtsBySliceAndContourICS,\
           LUT, PointData, FixContourData, ContourData