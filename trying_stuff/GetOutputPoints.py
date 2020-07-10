# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:57:32 2020

@author: ctorti
"""


"""
Function:
    GetOutputPoints()
    
Purpose:
    Parse the transformed points in outputpoints.txt, re-structure them 
    slice-by-slice (along rows) and contour-by-contour (along columns), 
    convert to Image Coordinate System, and add to existing LUT.


Input:
    FixPts_ICS       - Input points of Transformix in Image Coordinate System
    
    Origin          - Moving image origin
    
    Directions      - Moving image direction 
    
    Spacings        - Moving image pixel spacings
    
    Dimensions      - Moving image dimensions
    
    LUT   - A P x 6 array (P rows for all P points, 6 columns) created by 
            calling function GetInputPoints(), containing:
        
            - (col 0): the point number: integer starting from 0, 
            - (col 1): the input (fixed image) slice number it came from: integer starting from 0,    (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE NEXT ROW)
            - (col 2): the input (fixed image) contour number it belongs to: integer starting from 0, (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE PREVIOUS ROW)
            - (col 3): an empty list ([]) as a placeholder for the input (fixed image) index of the point: array of integers (N=3),
            - (col 4): the input (fixed image) point in PCS: array of floats (N=3),
            - (col 5): the input (fixed image) point in ICS: array of floats (N=3)
            
            e.g. [0, 
                  12, 
                  0, 
                  [],
                  [-19.904, -10.938, 2.551],
                  [89.940, 146.288, 12.000]]
            
            => The 1st point is part of the 1st contour (the 1st contour in
            absolute terms, not the 1st of possibly multiple contours in any
            given slice) that came from the 13th slice, whose coordinates in 
            the PCS and ICS are [-19,904, -10.938, 2.551] and 
            [89.940, 146.288, 12.000].
    
    
Returns:
    
    PtNos           - List of point numbers (starting from 0) parsed from 
                      outputpoints.txt ("Point")
    
    FixInds          - List of [x, y, z] indeces of the input (Fixed) points 
                      parsed from outputpoints.txt ("InputIndex")
    
    FixOutInds      - List of [x, y, z] ??? indeces parsed from 
                      outputpoints.txt ("OutputIndexFixed")
                      
    MovPts_PCS      - List of [x, y, z] transformed (Moving) points in Patient 
                      Coordinate System parsed from outputpoints.txt 
                      ("OutputPoint")
   
    MovPts_ICS      - List of [x, y, z] transformed (Moving) points converted
                      to Image Coordinate System
   
    Def_PCS         - List of [x, y, z] deformations from input-to-output
                      points in Patient Coordinate System parsed from 
                      outputpoints.txt ("Deformation")
   
    Def_ICS         - List of [x, y, z] deformations from input-to-output
                      points converted to Image Coordinate System 
   
    MovInds         - List of [x, y, z] indeces of the transformed (Moving) 
                      points parsed from outputpoints.txt ("OutputIndexMoving")
   
    MovPtsBySliceAndContour_PCS - List of [x, y, z] transformed (Moving) pts 
                                  in Patient Coordinate Syste arranged 
                                  slice-by-slice along rows and 
                                  contour-by-contour along rows
                                 
   
    MovPtsBySliceAndContour_ICS - List of [x, y, z] transformed (Moving) pts 
                                  in Image Coordinate Syste arranged 
                                  slice-by-slice along rows and 
                                  contour-by-contour along rows
   
    LUT   - A P x 12 array (P rows for all P points, 12 columns) containing:
        
            - (col 0): the point number: integer starting from 0, 
            - (col 1): the input (fixed image) slice number it came from: integer starting from 0,    (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE NEXT ROW)
            - (col 2): the input (fixed image) contour number it belongs to: integer starting from 0, (NOTE 08/07/20: THIS WAS PREVIOUSLY IN THE PREVIOUS ROW)
            - (col 3): the (fixed image) index of the point: array of integers (N=3),
            - (col 4): the (fixed image) point in PCS: array of floats (N=3),
            - (col 5): the (fixed image) point in ICS: array of floats (N=3),
            - (col 6): the (output / moving image) contour number of the transformed point: integer starting from 0,
            - (col 7): the (output / moving image) index of the transformed point: array of intergers (N=3),
            - (col 8): the (output / moving image) transformed point in PCS: array of floats (N=3),
            - (col 9): the (output / moving image) transformed point in ICS: array of floats (N=3),
            - (col 10): the deformation in PCS: array of floats (N=3),
            - (col 11): the deformation in ICS: array of floats (N=3)
            
            
            e.g. [125, 
                  14, 
                  2, 
                  [-19.904, -10.938, 2.551],
                  [89.940, 146.288, 12.000],
                  ...
                  ]
            
            => The 126th point is part of the 3rd contour (not the 3rd contour 
            in the 15th slice but the total running contour number) that came 
            from the 15th slice.
"""


def GetOutputPoints(FixPts_ICS, MovOrigin, MovDirs, MovSpacings, MovDims, LUT):
    """ Convert OutputPts generated by Transformix into an array of an array
    of points for each slice in MovingIm, using the z index from 
    MovingOutputIndex and the (x,y) points from OutputPoint. """
    from ParseTransformixOutput import ParseTransformixOutput
    import importlib
    import PCStoICS
    importlib.reload(PCStoICS)
    import GetMovContourNosByPoint
    importlib.reload(GetMovContourNosByPoint)
    
    from PCStoICS import PCStoICS
    from GetMovContourNosByPoint import GetMovContourNosByPoint
    
    PtNos, FixInds, FixPts_PCS,\
    FixOutInds, MovPts_PCS,\
    Def_PCS, MovInds = ParseTransformixOutput()
    
    
    # Convert MovPts_PCS to ICS:
    MovPts_ICS = PCStoICS(Pts_PCS=MovPts_PCS, Origin=MovOrigin, 
                          Directions=MovDirs, Spacings=MovSpacings)
    
    # Number of points:
    P = len(PtNos)
    
    # Loop through each point:
    for p in range(P):
        # The slice number (LUT[p][1]) should be equal to the z-index in 
        # FixInds (FixInds[2]):
        #if LUT[p][2] != FixInds[p][2]:
        if LUT[p][1] != FixInds[p][2]: # changed on 08/07/2020
            #print(f'\nLUT[{p}][2] (= {LUT[p][2]}) != FixInds[{p}][2] (= {FixInds[p][2]})')
            print(f'\nLUT[{p}][1] (= {LUT[p][1]}) != FixInds[{p}][2] (= {FixInds[p][2]})') # changed on 08/07/20
            
        # Update the LUT. 
        
        # Start by replacing the 4th item (empty array) by the array of 
        # [x, y, z] indeces:
        LUT[p][3] = FixInds[p]
        
        # Append the contour number of the transformed point. This will be 
        # inferred from the z-index of the transformed point (MovInds[p][2])
        # and the contour no that it came from (LUT[p][1]). 
        # For now append an empty array as a placeholder since this will be
        # done outside of this loop:
        LUT[p].append([]) # 7th item in list
        
        # Append the transformed point's indeces to LUT:
        LUT[p].append(MovInds[p]) # 8th item in list
        
        # Append the transformed point in PCS to LUT:
        LUT[p].append(MovPts_PCS[p]) # 9th item in list
        
        # Append the transformed point in ICS to LUT:
        LUT[p].append(MovPts_ICS[p]) # 10th item in list
        
        # Append the deformation in PCS to LUT:
        LUT[p].append(Def_PCS[p]) # 11th item in list
        
    
    Def_ICS = [[MovPts_ICS[p][i] - FixPts_ICS[p][i] for i in range(3)] for p in range(P)]
    
    
    # Get the list of contour numbers point-by-point (MovContourNosByPoint),
    # the indeces of the points slice-by-slice (along rows) and 
    # contour-by-contour (along columns) (PtIndsByMovSliceAndContour), and the
    # set of moving slice numbers that contain transformed contour points 
    # (MovSliceNosSet):
    MovContourNosByPoint,\
    PtIndsByMovSliceAndContour,\
    MovSliceNosSet = GetMovContourNosByPoint(MovInds=MovInds, LUT=LUT)
    
    
    # Append the deformation in ICS to LUT and modify the contour no of the 
    # transformed point using MovContourNosByPoint:
    for p in range(P):
        LUT[p].append(Def_ICS[p])  # 12th item in list
        
        LUT[p][6] = MovContourNosByPoint[p]
        
    
    
    # Use PtIndsByMovSliceAndContour to create a list of output (moving) 
    # contour points slice-by-slice along rows and contour-by-contour along 
    # columns.
      
    # The total number of slices in moving image (whether they contain contour
    # point or not):
    S = MovDims[2]
    
    # Initialise MovPtsBySliceAndContours_PCS and _ICS with a list of empty 
    # lists for each slice in moving image:
    MovPtsBySliceAndContour_PCS = []
    [MovPtsBySliceAndContour_PCS.append([]) for i in range(S)]
    
    MovPtsBySliceAndContour_ICS = []
    [MovPtsBySliceAndContour_ICS.append([]) for i in range(S)]
    
    # The number of slices in moving image that contain contours:
    #S = len(PtIndsByMovSliceAndContour)
    
    # Loop through each slice number in the set of slice numbers that contain
    # contour points (i.e. in MovSliceNosSet) and modify :
    for i in range(len(MovSliceNosSet)):
        #print(f'\ns = {s}')
        
        # The slice number for moving image:
        s = MovSliceNosSet[i]
        
        # The point indices grouped by contour for this slice:
        IndsByContours = PtIndsByMovSliceAndContour[i]
        
        # The number of groups of contour point indices in this slice:
        C = len(IndsByContours)
        
        # Initialise the contour-by-contour list of points for this slice:
        MovPtsThisSlice_PCS = []
        MovPtsThisSlice_ICS = []
        
        # Loop through each contour-by-contour list of points:
        for c in range(C):
            # The point indeces for this contour:
            IndsThisContour = IndsByContours[c]
            
            # This "contour" might consist of a single point.
            if isinstance(IndsThisContour, int):
                # IndsThisContour is an integer: 
                MovPtsThisContour_PCS = LUT[IndsThisContour][8]
                MovPtsThisContour_ICS = LUT[IndsThisContour][9]
                
            else:
                # IndsThisContour is a list of integers:
                MovPtsThisContour_PCS = [LUT[i][8] for i in IndsThisContour]
                MovPtsThisContour_ICS = [LUT[i][9] for i in IndsThisContour]
            
            MovPtsThisSlice_PCS.append(MovPtsThisContour_PCS)
            MovPtsThisSlice_ICS.append(MovPtsThisContour_ICS)
        
               
        # Modify the s^th items:
        MovPtsBySliceAndContour_PCS[s] = MovPtsThisSlice_PCS
        MovPtsBySliceAndContour_ICS[s] = MovPtsThisSlice_ICS
        
           
    return PtNos, FixInds, FixOutInds,\
           MovPts_PCS, MovPts_ICS,\
           Def_PCS, Def_ICS, MovInds,\
           MovPtsBySliceAndContour_PCS, MovPtsBySliceAndContour_ICS,\
           LUT
