# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:50:51 2021

@author: ctorti
"""

def get_seg_data_of_interest(seg, allF2SindsBySeg, segNums, segLabs, 
                             segLab, slcNum, p2c=False):
    # 16/07 need to tidy docstrings
    """   
    Get data of interest from a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        The SEG object.
    allF2SindsBySeg : list of a list of ints
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in seg.
    segNums : list of ints
        List of the segment number(s) that match the segment of interest.  
        The list may have a single int.
    segLabs : list of strs
        List of all segment labels. Only used if printing results to console.
    segLab : str or None
        The SEG label of interest, or None (if all segments are of interest).
    slcNum : int or None
        The index within the DICOM stack of the segmentation of interest, or
        None if all segmentations within a segment(s) is of interest.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
    f2sIndsBySeg : list of a list of int
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg, i.e. the list is a sub-list 
        of allF2SindsBySeg.
        
    Notes
    -----
    There are 3 possible main use cases:
        1. Copy a single segmentation
        2. Copy an entire segment
        3. Copy all segments
        
    Which main case applies depends on the inputs:
        slcNum
        segLab
    
    If slcNum != None (i.e. if a slice number is defined) a single 
    segmentation will be copied from the segment that matches segLab 
    (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple segments and 
    segLab = None.  Possible outcomes:

        - Raise exception with message "There are multiple segments so the 
        description of the desired segment must be provided"
        - Copy the segmentation to all segments that have a segmentation on 
        slice slcNum (* preferred option?)
        
    If slcNum = None but segLab != None, all segmentations within the
    segment given by srcSegLab will be copied (--> Main Use Case 2).  
    
    If slcNum = None and segLab = None, all segmentations within all
    segments will be copied (--> Main Use Case 3).
    
    If slcNum = None and segLab = None, all segments will be copied.  
    If not, get the segment number (zero-indexed) to be copied, and reduce 
    pixarrBySeg and f2sIndsBySeg to the corresponding item.
    """
    
    import numpy as np
    from copy import deepcopy
    #from pydicom import dcmread
    from general_tools.console_printing import print_inds_by_roi
    from general_tools.console_printing import print_pixarr_shape_by_seg
    #from general_tools.metadata import get_roicol_nums, get_roicol_labels
    from seg_tools.pixarrs import get_pixarrBySeg
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_seg_data_of_interest():')
        print('\n\n', '-'*120)
        print(f'   segLab = {segLab}')
        print(f'   slcNum = {slcNum}')
    
    
    # Get all pixel arrays by seg:
    """ 16/07: making changes to get_pixarrBySeg... """
    allPixarrBySeg = get_pixarrBySeg(seg, allF2SindsBySeg, p2c)
    
    #print('\n\nallPixarrBySeg =', allPixarrBySeg)
    #print(f'\n\nallF2SindsBySeg = {allF2SindsBySeg}')
    
    Nsegs = len(allF2SindsBySeg)
    
    # The shape of the first pixel array:
    F, R, C = allPixarrBySeg[0].shape
    
    if p2c:
        print_inds_by_roi(allF2SindsBySeg)
    
    # Initialise variable that indicates if the SEG data was reduced by
    # frame/slice number (FromSliceNum):
    reducedBySlc = False
            
    if slcNum:
        # Limit data to those which belong to the chosen slice number:
        pixarrBySeg = []
        f2sIndsBySeg = []
    
        for s in range(Nsegs):
            # The pixel array and frame-to-slice indices for this segment:
            allPixarr = deepcopy(allPixarrBySeg[s])
            allF2Sinds = deepcopy(allF2SindsBySeg[s])
                
            # Get the frame number(s) that relate to slcNum:
            frameNums = [i for i, x in enumerate(allF2Sinds) if x==slcNum]
            
            if p2c:
                print(f'\nallF2Sinds = {allF2Sinds}')
                print(f'frameNums = {frameNums}')
            
            if len(frameNums) != len(allF2Sinds):
                # Some frames will be rejected:
                reducedBySlc = True
                
                # Initialise pixarr with the reduced number of frames matching
                # len(frameNums):
                pixarr = np.zeros((len(frameNums), R, C), dtype='uint')
                f2sInds = []
                
                if frameNums:
                    # Keep only the indeces and frames that relate to frameNums:
                    for f in range(len(frameNums)):
                        pixarr[f] = allPixarr[frameNums[f]]
                        
                        f2sInds.append(allF2Sinds[frameNums[f]])
                    
                    if p2c:
                        print(f'f2sInds = {f2sInds}')
                    
                    # Append non-empty pixel arrays and non-empty f2sInds:
                    pixarrBySeg.append(pixarr)
                    f2sIndsBySeg.append(f2sInds)
            """ 
            else 
            all frames to remain
            """
        
        #print('\n\npixarrBySeg =', pixarrBySeg)
        
        if reducedBySlc:
            # Replace allPixarrBySeg and allF2SindsBySeg with pixarrBySeg and 
            # f2sIndsBySeg in case further restricting of data is required 
            # below:
            allPixarrBySeg = deepcopy(pixarrBySeg)
            allF2SindsBySeg = deepcopy(f2sIndsBySeg)
        """ 
        else 
            allPixarrBySeg and allF2SindsBySeg remain unchanged
        """
        
        #print('\n\nallPixarrBySeg =', allPixarrBySeg)
        
        if p2c and reducedBySlc:
            print('\n   After limiting data to those that relate to slice',
                  f'number {slcNum}:')#, the f2sIndsBySeg =')
            print_inds_by_roi(f2sIndsBySeg)
            print('   pixarrBySeg = ', pixarrBySeg)
            print_pixarr_shape_by_seg(pixarrBySeg)
    
    # Initialise variable that indicates if the SEG data was reduced by
    # chosen segment label (FromSegLabel):
    #ReducedBySeg = False
    
    if segLab:
        # Limit data to those which belong to the chosen segment(s):
        pixarrBySeg = []
        f2sIndsBySeg = []
        
        #print(f'\nsegNums = {segNums}')
        #print(f'allF2SindsBySeg = {allF2SindsBySeg}')
        
        if len(segNums) != len(allF2SindsBySeg):
            #ReducedBySeg = True
            
            pixarrBySeg = [allPixarrBySeg[s] for s in segNums]
            f2sIndsBySeg = [allF2SindsBySeg[s] for s in segNums]
            
            if p2c:
                # Limit the list of segment names to those that belong to the
                # chosen segment(s):
                segLabs = [segLabs[i] for i in segNums]
            
                print('\n   After limiting data to those whose segment name',
                      f'matches {segLabs}:')
                print_inds_by_roi(f2sIndsBySeg)
                print_pixarr_shape_by_seg(pixarrBySeg)
                
        else:
            pixarrBySeg = deepcopy(allPixarrBySeg)
            f2sIndsBySeg = deepcopy(allF2SindsBySeg)
        
    else:
        pixarrBySeg = deepcopy(allPixarrBySeg)
        f2sIndsBySeg = deepcopy(allF2SindsBySeg)
    
    if p2c:
        print('\n   Final outputs of GetSegDataOfInterest():')
        print_pixarr_shape_by_seg(pixarrBySeg)
        print_inds_by_roi(f2sIndsBySeg)
        print('-'*120)
        
    return pixarrBySeg, f2sIndsBySeg

def raise_error_if_no_seg_data_of_interest(f2sIndsBySeg, segLabs, slcNum):
    """
    Raise error if there is no matching SEG data of interest.

    Parameters
    ----------
    f2sIndsBySeg : list of a list of int
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg, i.e. the list is a sub-list 
        of allF2SindsBySeg.
    segLabs : list of strs
        List (for each segment) of the segment label.
    slcNum : int or None

    Returns
    -------
    None.

    """
    
    if f2sIndsBySeg == []:
        msg = "There are no segmentations "
        
        if slcNum:
            msg += f"on slice {slcNum} "
              
        msg += f"in any segment. There are {len(segLabs)} segments:"
        
        for i in range(len(segLabs)):
            msg += f"\nSegment {i+1} with label '{segLabs[i]}' has "\
                   + f"{len(f2sIndsBySeg[i])} segmentations that correspond "\
                   + f"to slices {f2sIndsBySeg[i]}" 
        raise Exception(msg)


""" Old below """







def get_seg_data_of_interest_OLD(segFpath, slcNum, segLab, dicomDir, 
                             p2c=False):
    """
    Get data of interest from a SEG.
    
    Parameters
    ----------
    segFpath : str
        Filepath of the Source SEG file.
    slcNum : int or None
        The slice indeces within the DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a 
        single segmentation). The index is zero-indexed.
        If slcNum = None, all segmentations attributed to all slices will be
        parsed.
    segLab : str or None
        All or part of the SegmentLabel of the segment containing the 
        segmentation(s) to be copied.  If segLab = None, all segments will be
        parsed.
    dicomDir : str 
        Directory containing the corresponding DICOMs.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
    f2sIndsBySeg : List of a list of int
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg.
        
    Notes
    -----
    There are 3 possible main use cases:
        1. Copy a single segmentation
        2. Copy an entire segment
        3. Copy all segments
        
    Which main case applies depends on the inputs:
        slcNum
        segLab
    
    If slcNum != None (i.e. if a slice number is defined) a single 
    segmentation will be copied from the segment that matches segLab 
    (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple segments and 
    segLab = None.  Possible outcomes:

        - Raise exception with message "There are multiple segments so the 
        description of the desired segment must be provided"
        - Copy the segmentation to all segments that have a segmentation on 
        slice slcNum (* preferred option?)
        
    If slcNum = None but segLab != None, all segmentations within the
    segment given by srcSegLab will be copied (--> Main Use Case 2).  
    
    If slcNum = None and segLab = None, all segmentations within all
    segments will be copied (--> Main Use Case 3).
    
    If slcNum = None and segLab = None, all segments will be copied.  
    If not, get the segment number (zero-indexed) to be copied, and reduce 
    pixarrBySeg and f2sIndsBySeg to the corresponding item.
    """
    
    import numpy as np
    from copy import deepcopy
    from pydicom import dcmread
    from general_tools.console_printing import print_inds_by_roi
    from general_tools.console_printing import print_pixarr_shape_by_seg
    from general_tools.metadata import get_roicol_nums, get_roicol_labels
    from seg_tools.pixarrs import get_pixarrBySeg
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_seg_data_of_interest():')
        print('\n\n', '-'*120)
        print(f'   segLab = {segLab}')
        print(f'   slcNum = {slcNum}')
    
    seg = dcmread(segFpath)
    
    pixarrBySeg_all, f2sIndsBySeg_all = get_pixarrBySeg(seg, dicomDir, 
                                                      p2c)
    
    #print('\n\npixarrBySeg_all =', pixarrBySeg_all)
    #print(f'\n\nf2sIndsBySeg_all = {f2sIndsBySeg_all}')
    
    Nsegs = len(f2sIndsBySeg_all)
    
    """ Get the segment labels: """
    segLabs = get_roicol_labels(Roi=seg)
    
    """ The shape of the first pixel array: """
    F, R, C = pixarrBySeg_all[0].shape
    
    if p2c:
        print(f'   segLabs = {segLabs}')
        #print('   f2sIndsBySeg_all =')
        print_inds_by_roi(f2sIndsBySeg_all)
    
    """ Initialise variable that indicates if the SEG data was reduced by
    frame/slice number (FromSliceNum): """
    reducedBySlc = False
            
    if slcNum:
        """ Limit data to those which belong to the chosen slice number. """
        
        pixarrBySeg = []
        f2sIndsBySeg = []
    
        for s in range(Nsegs):
            """ The pixel array and frame-to-slice indices for this segment:"""
            pixarr_all = deepcopy(pixarrBySeg_all[s])
            f2sInds_all = deepcopy(f2sIndsBySeg_all[s])
                
            """ Get the frame number(s) that relate to slcNum: """
            FrmNums = [i for i, x in enumerate(f2sInds_all) if x==slcNum]
            
            if p2c:
                print(f'\nf2sInds_all = {f2sInds_all}')
                print(f'FrmNums = {FrmNums}')
            
            if len(FrmNums) != len(f2sInds_all):
                """ Some frames will be rejected: """
                reducedBySlc = True
                
                """ Initialise pixarr with the reduced number of frames 
                matching len(FrmNums): """
                pixarr = np.zeros((len(FrmNums), R, C), dtype='uint')
                f2sInds = []
                
                if FrmNums:
                    """ Keep only the indeces and frames that relate to 
                    FrmNums. """
                    #pixarr = [pixarr_all[f] for f in FrmNums]
                    
                    #""" Set pixarr to the first frame in FrmNums: """
                    #pixarr = pixarr_all[0]
                    #
                    #if len(FrmNums) > 1:
                    #    """ Add additional frames: """
                    #    for f in range(1, len(FrmNums)):
                    #        pixarr = np.vstack((pixarr, pixarr_all[f]))
                    
                    for f in range(len(FrmNums)):
                        pixarr[f] = pixarr_all[FrmNums[f]]
                        
                        f2sInds.append(f2sInds_all[FrmNums[f]])
                    
                    if p2c:
                        print(f'f2sInds = {f2sInds}')
                    
                    """ Append non-empty pixel arrays and non-empty f2sInds:"""
                    pixarrBySeg.append(pixarr)
                    f2sIndsBySeg.append(f2sInds)
                
                #else:
                #    #pixarr = []
                #    #pixarr = np.zeros((0, R, C), dtype='uint')
                #    f2sInds = []
                
                #""" Append empty pixel arrays and empty f2sInds: """
                #pixarrBySeg.append(pixarr)
                #f2sIndsBySeg.append(f2sInds)
            """ else all frames to remain. """
        
        #print('\n\npixarrBySeg =', pixarrBySeg)
        
        if reducedBySlc:
            """ Replace pixarr_allByRoi and f2sIndsBySeg_all with pixarrBySeg  
            and f2sIndsBySeg in case further restricting of data is required 
            below: """
            pixarrBySeg_all = deepcopy(pixarrBySeg)
            f2sIndsBySeg_all = deepcopy(f2sIndsBySeg)
        """ else pixarrBySeg_all and f2sIndsBySeg_all remain unchanged. """
        
        #print('\n\npixarrBySeg_all =', pixarrBySeg_all)
        
        if p2c and reducedBySlc:
            print('\n   After limiting data to those that relate to slice',
                  f'number {slcNum}:')#, the f2sIndsBySeg =')
            print_inds_by_roi(f2sIndsBySeg)
            print('   pixarrBySeg = ', pixarrBySeg)
            print_pixarr_shape_by_seg(pixarrBySeg)
    
    
    """ Initialise variable that indicates if the SEG data was reduced by
    chosen segment label (FromSegLabel): """
    #ReducedBySeg = False
    
    if segLab:
        """ Limit data to those which belong to the chosen segment(s). """
        
        pixarrBySeg = []
        f2sIndsBySeg = []
        
        """ Get the segment number(s) whose name matches segLab: """
        segNums = get_roicol_nums(Roi=seg, searchStr=segLab)
        
        #print(f'\nsegNums = {segNums}')
        #print(f'f2sIndsBySeg_all = {f2sIndsBySeg_all}')
        
        if len(segNums) != len(f2sIndsBySeg_all):
            #ReducedBySeg = True
            
            pixarrBySeg = [pixarrBySeg_all[s] for s in segNums]
            f2sIndsBySeg = [f2sIndsBySeg_all[s] for s in segNums]
            
            if p2c:
                """ Get the names of all ROIs: """
                segLabs = get_roicol_labels(Roi=seg)
            
                """ Limit the list of segment descriptions to those that belong to 
                the chosen segment(s): """
                segLabs = [segLabs[i] for i in segNums]
            
                print('\n   After limiting data to those whose segment name',
                      f'matches {segLabs}:')#, the f2sIndsBySeg =')
                print_inds_by_roi(f2sIndsBySeg)
                print_pixarr_shape_by_seg(pixarrBySeg)
                
        else:
            pixarrBySeg = deepcopy(pixarrBySeg_all)
            f2sIndsBySeg = deepcopy(f2sIndsBySeg_all)
        
    else:
        pixarrBySeg = deepcopy(pixarrBySeg_all)
        f2sIndsBySeg = deepcopy(f2sIndsBySeg_all)
    
    if p2c:
        print('\n   Final outputs of GetSegDataOfInterest():')
        print_pixarr_shape_by_seg(pixarrBySeg)
        print_inds_by_roi(f2sIndsBySeg)
        print('-'*120)
        
    return pixarrBySeg, f2sIndsBySeg