# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 14:27:54 2020

@author: ctorti
"""



"""
Function:
    GetMovContourNosByPoint()
    
Purpose:
    To group the transformed (moving) points by contours and slices.  The 
    slice numbers are given by the z-index (parsed from outputpoints.txt) and
    the contour numbers are assigned based on the contour they were transformed
    from and which slice they were transformed to.


Input:
    MovInds - List of [x, y, z] indeces of the transformed (Moving) points
              parsed from outputpoints.txt ("OutputIndexMoving")
    
    LUT     - A P x 6 array (P rows for all P points, 6 columns) created by 
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
    
    MovContourNosByPoint       - List of contour numbers for each contour point
    
    PtIndsByMovSliceAndContour - List of point indeces grouped by slice (along
                                 rows) and contour (along columns) for every
                                 slice in moving image that has contour points
                                 (i.e. not for every slice in moving image, as
                                 some slices will not necessarily have contour
                                 points)
    
    MovSliceNosSet             - List of slice numbers for moving image that
                                 has contour points
   
  
"""



def GetMovContourNosByPoint(MovInds, LUT):
    
    # Get list of all moving slice numbers (z-indeces in MovInds):
    MovSliceNos = [MovInds[i][2] for i in range(len(MovInds))]
    
    # Get list of unique slice numbers:
    MovSliceNosSet = list(set(MovSliceNos))
    
    #print('MovSliceNosSet =', MovSliceNosSet)
    
    # Initialise indeces of the points that correspond to the slice numbers in
    # MovSliceNosSet, and the corresponding fixed contour numbers organised by 
    # moving slices:
    PtIndsByMovSlice = []
    FixContourNosByMovSlice = []
    
    # Get indeces of all slice nos in MovSliceNos that are equal to each 
    # index in MovSliceNosSet:
    for zind in MovSliceNosSet:
        inds = [i for i, e in enumerate(MovSliceNos) if e==zind]
        
        PtIndsByMovSlice.append(inds)
        
        #FixContourNosByMovSlice.append([LUT[ind][1] for ind in inds])
        FixContourNosByMovSlice.append([LUT[ind][2] for ind in inds]) # Changed on 08/07/2020
        
    #print('\nPtIndsByMovSlice =', PtIndsByMovSlice)
    
    #print('\nFixContourNosByMovSlice =', FixContourNosByMovSlice)
    
    # Initialise array of moving contour numbers organised by moving slices, 
    # and a moving contour number counter:
    MovContourNosByMovSlice = []
    MovContourNo = 0
    
    # Loop through each array in FixContourNoByMovSlice:
    for FixContourNosThisSlice in FixContourNosByMovSlice:
        
        # Initialise array of moving contour numbers for this slice:
        MovContourNosThisSlice = []
        
        for i in range(len(FixContourNosThisSlice)):
            ThisFixContourNo = FixContourNosThisSlice[i]
            
            # If i=0 this is the start of a new contour.
            if i == 0:
                MovContourNosThisSlice.append(MovContourNo)
            
            else:
                PreviousFixContourNo = FixContourNosThisSlice[i-1]
                
                # If this FixContourNo is the same as the previous:
                if ThisFixContourNo == PreviousFixContourNo:
                    # This point is part of the previous contour:
                    MovContourNosThisSlice.append(MovContourNo)
                    
                else:
                    # This is the start of a different contour:
                    MovContourNo = MovContourNo + 1
                    
                    MovContourNosThisSlice.append(MovContourNo)
                    
        
        MovContourNosByMovSlice.append(MovContourNosThisSlice)
        
        # Next slice will be start of next contour no:
        MovContourNo = MovContourNo + 1
            
    
    #print('\nMovContourNosByMovSlice =', MovContourNosByMovSlice)        
    
    
    # Initialise a list of moving contour nos in the order of the points 
    # themselves, by using the indeces in PtIndsByMovSlice to index 
    # MovContourNosBySlice:
    MovContourNosByPoint = []
    [MovContourNosByPoint.append([]) for i in range(len(MovInds))]
    
    # Number of slices:
    S = len(MovContourNosByMovSlice)
    
    # Loop through each array in MovContourNosByMovSlice:
    for s in range(S):
        
        MovContourNosThisSlice = MovContourNosByMovSlice[s]
        
        inds = PtIndsByMovSlice[s]
        
        #print('\nMovContourNosThisSlice =', MovContourNosThisSlice)
        
        #print('\ninds =', inds)
        
        # Number of contours numbers in this slice:
        C = len(MovContourNosThisSlice)
        
        # Loop through each array in MovContourNosThisSlice:
        for c in range(C):
            ThisMovContourNo = MovContourNosThisSlice[c]
            
            # The index in the list of points for this contour:
            ind = inds[c]
            
            MovContourNosByPoint[ind] = ThisMovContourNo
            
    
    #print('\nMovContourNosByPoint =', MovContourNosByPoint) 
    #print('\nMovContourNosByPoint[0:2] =', [MovContourNosByPoint[i] for i in range(2)]) 
    #N = sum([len(MovContourNosByMovSlice[i]) for i in range(2)])
    #print(f'\nMovContourNosByPoint[0:{N}] =', [MovContourNosByPoint[i] for i in range(N)]) 
    
    
    """
    Create a list of point indeces organised slice-by-slice along rows and
    contour-by-contour along columns using MovContourNosByPoint.
    """
    
    # Number of slices:
    S = len(PtIndsByMovSlice)
    
    # Initialise array of point indeces and moving contour numbers 
    # slice-by-slice along rows and contour by contour along columns:
    PtIndsByMovSliceAndContour = []
    MovContourNosByMovSliceAndContour = []
    
    # Loop through each slice in PtIndsByMovSlice:
    for s in range(S):
        #if s < 2:
        #    print('\ns =', s)
    
        # The point indices and moving contour numbers for this slice:
        PtIndsThisSlice = PtIndsByMovSlice[s]
        MovContourNosThisSlice = MovContourNosByMovSlice[s]
        
        #if s < 2:
        #    print('\nMovContourNosThisSlice =', MovContourNosThisSlice)
        #    print('\nPtIndsThisSlice        =', PtIndsThisSlice)
        
        # A list of unique contour numbers for this slice:
        UniqueContourNos = list(set(MovContourNosThisSlice))
        
        if len(UniqueContourNos) == 1:
            # All items in this slice belong to the same contour:
            PtIndsByMovSliceAndContour.append(PtIndsThisSlice)
            MovContourNosByMovSliceAndContour.append(MovContourNosThisSlice)
            
        else:
            # The point indices and contour numbers need to be split and grouped
            # by contour number.     
            # Get the indices of each distinct contour number in MovContourNosThisSlice:
            IndsByContour = [[i for i, e in enumerate(MovContourNosThisSlice) if e==c] for c in UniqueContourNos]
            
            # Use each set of indeces in IndsByContour to index the point indices and moving
            # contour numbers in PtIndsThisSlice and MovContourNosThisSlice:
            PtIndsThisSliceByContour = [[PtIndsThisSlice[i] for i in IndsByContour[j]] for j in range(len(IndsByContour))]
            MovContourNosThisSliceByContour = [[MovContourNosThisSlice[i] for i in IndsByContour[j]] for j in range(len(IndsByContour))]
            
            #if s < 2:
            #    print('\nMovContourNosThisSliceByContour =', MovContourNosThisSliceByContour)
            #    print('\nPtIndsThisSliceByContour =', PtIndsThisSliceByContour)
                
            # Append to PtIndsByMovSliceAndContour and MovContourNosByMovSliceAndContour:
            PtIndsByMovSliceAndContour.append(PtIndsThisSliceByContour)
            MovContourNosByMovSliceAndContour.append(MovContourNosThisSliceByContour)
            
        #if s < 2:
        #    print(f'\nMovContourNosByMovSliceAndContour[{s}] =', MovContourNosByMovSliceAndContour[s])
        #    print(f'\nPtIndsByMovSliceAndContour[{s}] =', PtIndsByMovSliceAndContour[s])     
        
    #print('\nMovContourNosByMovSliceAndContour =', MovContourNosByMovSliceAndContour) 
    #print('\nMovContourNosByMovSliceAndContour[0:2] =', [MovContourNosByMovSliceAndContour[i] for i in range(2)]) 
    
    #print('\nPtIndsByMovSliceAndContour =', PtIndsByMovSliceAndContour) 
    #print('\nPtIndsByMovSliceAndContour[0:2] =', [PtIndsByMovSliceAndContour[i] for i in range(2)]) 

    return MovContourNosByPoint, PtIndsByMovSliceAndContour, MovSliceNosSet       
