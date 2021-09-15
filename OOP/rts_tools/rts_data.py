# -*- coding: utf-8 -*-


#from copy import deepcopy

from importlib import reload
import general_tools.console_printing
reload(general_tools.console_printing)

from general_tools.console_printing import (
    print_indsByRoi, print_ptsByCntByRoi
    )
from conversion_tools.inds_pts_cntdata import (
    cntdata_to_pts, ptsByCntByRoi_to_cntdataByCntByRoi
    )

def get_rts_data_of_interest(
        rts, allPtsByCntByRoi, allC2SindsByRoi, roiNums, 
        allRoiNames, roiName, slcNum, p2c=False
        ):
    """
    Get data of interest from an RTS.
    
    Parameters
    ----------
    rts : Pydicom object
        Pydicom representation of an RTSTRUCT file.
    allPtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    allC2SindsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour. 
    roiNums : list of ints
        A list of indices corresponding to the ROI(s) of interest.
        The list may be of length 1.
    allRoiNames : list of strs
        List of all ROI names. Only used if printing results to console.
    roiName : str or None
        The ROIName of interest, or None if all ROIs are of interest.
    slcNum : int or None
        The slice indeces within the DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If slcNum = None, a relationship-preserving copy will be made.
    p2c : bool, optional
        If True some results will be printed to the console. The default value
        is False.
                       
    Returns
    -------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates to be copied from 
        the RTS.
    cntdataByCntByRoi : list of a list of a list of strs
        List (for each ROI) of a list (for all contours) of a flat list of 
        coordinates in ptsByCntByRoi converted from floats to strings.
    c2sIndsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour to be copied from the RTS.
    
    Notes
    -----
    There are 3 possible main use cases:
        1. Copy a single contour
        2. Copy an entire ROI
        3. Copy all ROIs
        
    Which main case applies depends on the inputs:
        slcNum
        roiName
    
    If slcNum != None (i.e. if a slice number is defined) a single contour
    will be copied from the ROI that matches roiName (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple ROIs and 
    roiName = None.  Possible outcomes:

        - Raise exception with message "There are multiple ROIs so the ROI
        label must be provided"
        - Copy the contour to all ROIs that have a contour on slice 
        slcNum (* preferred option?)
        
    If slcNum = None but roiName != None, all contours within the ROI given
    by roiName will be copied (--> Main Use Case 2).  
    
    If slcNum = None and roiName = None, all contours within all ROIs will be
    copied (--> Main Use Case 3).
    
    If slcNum = None and roiName = None, all ROIs will be copied.  
    If not, get the ROI number (zero-indexed) to be copied, and reduce 
    ptsByCntByRoi and c2sIndsByRoi to the corresponding item.
    """
    
    #import DicomTools
    #import importlib
    #importlib.reload(DicomTools)
    #from DicomTools import GetRoiNums, GetRoiLabels
    #from copy import deepcopy
    #from pydicom import dcmread
    
    
    #allPtsByCntByRoi, allC2SindsByRoi = get_ptsByCntByRoi(Rts, DicomDir, p2c)
    
    Nrois = len(allC2SindsByRoi) # number of ROIs
    
    if p2c:
        #from GeneralTools import PrintIndsByRoi, PrintPtsByCntByRoi
        
        print('\n\n', '-'*120)
        print('Running of get_rts_data_of_interest():')
        print('\n\n', '-'*120)
        print('   allC2SindsByRoi =')
        print_indsByRoi(allC2SindsByRoi)
        #print(f'   allC2SindsByRoi = {allC2SindsByRoi}')
    
    reducedBySlice = False
            
    if slcNum:
        # Limit data to those which belong to the chosen slice number.
        ptsByCntByRoi = []
        c2sIndsByRoi = []
    
        for r in range(Nrois):
            # All points-by-contour and contour-to-slice indices for this ROI:
            #allPtsByCnt = deepcopy(allPtsByCntByRoi[r])
            #allC2Sinds = deepcopy(allC2SindsByRoi[r])
            allPtsByCnt = list(allPtsByCntByRoi[r])
            allC2Sinds = list(allC2SindsByRoi[r])
                
            # The contour number(s) that relate to slcNum:
            cntNums = [i for i, e in enumerate(allC2Sinds) if e==slcNum]
            
            #print(f'\ncntNums = {cntNums}')
            #print(f'allC2Sinds = {allC2Sinds}')
            
            if len(cntNums) != len(allC2Sinds):
                # Some contours will be rejected.
                reducedBySlice = True
            
                if cntNums:
                    # Keep only the indeces and points that relate to cntNums:
                    ptsByCnt = [allPtsByCnt[c] for c in cntNums]
                    c2sInds = [allC2Sinds[c] for c in cntNums]
                else:
                    ptsByCnt = []
                    c2sInds = []
                
                ptsByCntByRoi.append(ptsByCnt)
                c2sIndsByRoi.append(c2sInds)
            #else:
                # All contours to remain.
                
        
        if reducedBySlice:
            """ 
            Replace allPtsByCntByRoi and allC2SindsByRoi with ptsByCntByRoi and 
            c2sIndsByRoi in case further restricting of data is required below:
            """
            #allPtsByCntByRoi = deepcopy(ptsByCntByRoi)
            #allC2SindsByRoi = deepcopy(c2sIndsByRoi)
            allPtsByCntByRoi = list(ptsByCntByRoi)
            allC2SindsByRoi = list(c2sIndsByRoi)
        #else:
            # allPtsByCntByRoi and allC2SindsByRoi remain unchanged.
        
        if p2c and reducedBySlice:
            print('\n   After limiting data to those that relate to slice',
                  f'number {slcNum}, the c2sIndsByRoi =')
            print_indsByRoi(c2sIndsByRoi)
            #print('\n   After limiting data to those that relate to slice',
            #      f'number {slcNum}, the c2sIndsByRoi = {c2sIndsByRoi}')
    
    
    #ReducedByRoi = False
    
    if roiName:
        # Limit data to those which belong to the chosen ROI(s).
        ptsByCntByRoi = []
        c2sIndsByRoi = []
        
        # Get the ROI number(s) whose name matches roiName:
        #roiNums = GetRoiNums(Roi=Rts, SearchString=roiName)
        
        #print(f'\nroiNums = {roiNums}')
        #print(f'allC2SindsByRoi = {allC2SindsByRoi}')
        
        if len(roiNums) != len(allC2SindsByRoi):
            #ReducedByRoi = True
            
            ptsByCntByRoi = [allPtsByCntByRoi[r] for r in roiNums]
            c2sIndsByRoi = [allC2SindsByRoi[r] for r in roiNums]
            
            if p2c:
                # Get the names of all ROIs:
                #allRoiNames = GetRoiLabels(Roi=Rts)
            
                # Limit the list of ROI names to those that belong to the 
                # chosen ROIs(s):
                roiNames = [allRoiNames[i] for i in roiNums]
            
                print('\n   After limiting data to those whose ROI name',
                      f'matches {roiNames}, the c2sIndsByRoi =')
                print_indsByRoi(c2sIndsByRoi)
                #print('\n   After limiting data to those whose ROI name',
                #      f'matches {roiNames}, the c2sIndsByRoi = {c2sIndsByRoi}')
                
        else:
            #ptsByCntByRoi = deepcopy(allPtsByCntByRoi)
            #c2sIndsByRoi = deepcopy(allC2SindsByRoi)
            ptsByCntByRoi = list(allPtsByCntByRoi)
            c2sIndsByRoi = list(allC2SindsByRoi)
    else:
            #ptsByCntByRoi = deepcopy(allPtsByCntByRoi)
            #c2sIndsByRoi = deepcopy(allC2SindsByRoi)
            ptsByCntByRoi = list(allPtsByCntByRoi)
            c2sIndsByRoi = list(allC2SindsByRoi)
    
    cntdataByCntByRoi = ptsByCntByRoi_to_cntdataByCntByRoi(
        ptsByCntByRoi
        )
    
    if p2c:
        print('\n   Final outputs of get_rts_data_of_interest():')
        print_ptsByCntByRoi(ptsByCntByRoi)
        print_indsByRoi(c2sIndsByRoi)
        #print(f'   c2sIndsByRoi = {c2sIndsByRoi}')
        #print(f'   ptsByCntByRoi = {ptsByCntByRoi}')
        print('-'*120)
        
    return ptsByCntByRoi, cntdataByCntByRoi, c2sIndsByRoi


def raise_error_if_no_rts_data_of_interest(c2sIndsByRoi, allRoiNames, slcNum):
    """
    Raise error if there is no matching RTSTRUCT data of interest.
    
    Parameters
    ----------
    c2sIndsByRoi : list of a list of int
        List (for each ROI) of the slice numbers that correspond to each
        contour in ContourData, i.e. the list is a sub-list of allC2SindsByRoi.
    allRoiNames : list of strs
        List (for each ROI) of ROI name(s).
    slcNum : int or None
        The slice number of interest, or None if all slices are of interest.
    
    Returns
    -------
    None.
    """
    
    if c2sIndsByRoi == []:
        msg = "There are no contours "
        
        if slcNum:
            msg += f"on slice {slcNum} "
              
        msg += f"in any ROI. There are {len(allRoiNames)} ROIs:"
        
        for i in range(len(allRoiNames)):
            msg += f"\nROI {i+1} with name '{allRoiNames[i]}' has "\
                   + f"{len(c2sIndsByRoi[i])} contours that correspond "\
                   + f"to slices {c2sIndsByRoi[i]}" 
        raise Exception(msg)