# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:54:12 2020

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
SEGMENT FUNCTIONS
******************************************************************************
******************************************************************************
"""      

def GetSegDataOfInterest(SegFpath, SliceNum, SegLabel, DicomDir, 
                         LogToConsole=False):
    """
    06/07/21:
        See get_seg_data_of_interest in seg_tools.get_seg_data_of_interest.
        
    Get data of interest from a SEG.
    
    Inputs:
    ******
    
    SegFpath : string
        Filepath of the Source SEG file.
        
    SliceNum : integer or None
        The slice indeces within the DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a 
        single segmentation). The index is zero-indexed.
        If SliceNum = None, all segmentations attributed to all slices will be
        parsed.
        
    SegLabel : string or None
        All or part of the SegmentLabel of the segment containing the 
        segmentation(s) to be copied.  If SegLabel = None, all segments will be
        parsed.
    
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
        
    F2SindsBySeg : List of a list of integers
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in PixArrBySeg.
        
        
    Notes:
    *****
    
    There are 3 possible main use cases:
        1. Copy a single segmentation
        2. Copy an entire segment
        3. Copy all segments
        
    Which main case applies depends on the inputs:
        SliceNum
        SegLabel
    
    If SliceNum != None (i.e. if a slice number is defined) a single 
    segmentation will be copied from the segment that matches SegLabel 
    (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple segments and 
    SegLabel = None.  Possible outcomes:

        - Raise exception with message "There are multiple segments so the 
        description of the desired segment must be provided"
        - Copy the segmentation to all segments that have a segmentation on 
        slice SliceNum (* preferred option?)
        
    If SliceNum = None but SegLabel != None, all segmentations within the
    segment given by FromSegLabel will be copied (--> Main Use Case 2).  
    
    If SliceNum = None and SegLabel = None, all segmentations within all
    segments will be copied (--> Main Use Case 3).
    
    If SliceNum = None and SegLabel = None, all segments will be copied.  
    If not, get the segment number (zero-indexed) to be copied, and reduce 
    PixArrBySeg and F2SindsBySeg to the corresponding item.
    """
    
    import numpy as np
    from DicomTools import GetRoiNums, GetRoiLabels
    from copy import deepcopy
    from pydicom import dcmread
    
    if LogToConsole:
        from GeneralTools import PrintIndsByRoi, PrintPixArrShapeBySeg
        
        print('\n\n', '-'*120)
        print('Running of GetSegDataOfInterest():')
        print('\n\n', '-'*120)
        print(f'   SegLabel = {SegLabel}')
        print(f'   SliceNum = {SliceNum}')
    
    
    Seg = dcmread(SegFpath)
    
    AllPixArrBySeg, AllF2SindsBySeg = GetPixArrBySeg(Seg, DicomDir, 
                                                     LogToConsole)
    
    #print('\n\nAllPixArrBySeg =', AllPixArrBySeg)
    #print(f'\n\nAllF2SindsBySeg = {AllF2SindsBySeg}')
    
    Nsegs = len(AllF2SindsBySeg)
    
    """ Get the segment labels: """
    SegLabels = GetRoiLabels(Roi=Seg)
    
    """ The shape of the first pixel array: """
    F, R, C = AllPixArrBySeg[0].shape
    
    if LogToConsole:
        print(f'   SegLabels = {SegLabels}')
        #print('   AllF2SindsBySeg =')
        PrintIndsByRoi(AllF2SindsBySeg)
    
    """ Initialise variable that indicates if the SEG data was reduced by
    frame/slice number (FromSliceNum): """
    ReducedBySlice = False
            
    if SliceNum:
        """ Limit data to those which belong to the chosen slice number. """
        
        PixArrBySeg = []
        F2SindsBySeg = []
    
        for s in range(Nsegs):
            """ The pixel array and frame-to-slice indices for this segment:"""
            AllPixArr = deepcopy(AllPixArrBySeg[s])
            AllF2Sinds = deepcopy(AllF2SindsBySeg[s])
                
            """ Get the frame number(s) that relate to SliceNum: """
            FrmNums = [i for i, x in enumerate(AllF2Sinds) if x==SliceNum]
            
            if LogToConsole:
                print(f'\nAllF2Sinds = {AllF2Sinds}')
                print(f'FrmNums = {FrmNums}')
            
            if len(FrmNums) != len(AllF2Sinds):
                """ Some frames will be rejected: """
                ReducedBySlice = True
                
                """ Initialise PixArr with the reduced number of frames 
                matching len(FrmNums): """
                PixArr = np.zeros((len(FrmNums), R, C), dtype='uint')
                F2Sinds = []
                
                if FrmNums:
                    """ Keep only the indeces and frames that relate to 
                    FrmNums. """
                    #PixArr = [AllPixArr[f] for f in FrmNums]
                    
                    #""" Set PixArr to the first frame in FrmNums: """
                    #PixArr = AllPixArr[0]
                    #
                    #if len(FrmNums) > 1:
                    #    """ Add additional frames: """
                    #    for f in range(1, len(FrmNums)):
                    #        PixArr = np.vstack((PixArr, AllPixArr[f]))
                    
                    for f in range(len(FrmNums)):
                        PixArr[f] = AllPixArr[FrmNums[f]]
                        
                        F2Sinds.append(AllF2Sinds[FrmNums[f]])
                    
                    if LogToConsole:
                        print(f'F2Sinds = {F2Sinds}')
                    
                    """ Append non-empty pixel arrays and non-empty F2Sinds:"""
                    PixArrBySeg.append(PixArr)
                    F2SindsBySeg.append(F2Sinds)
                
                #else:
                #    #PixArr = []
                #    #PixArr = np.zeros((0, R, C), dtype='uint')
                #    F2Sinds = []
                
                #""" Append empty pixel arrays and empty F2Sinds: """
                #PixArrBySeg.append(PixArr)
                #F2SindsBySeg.append(F2Sinds)
            """ else all frames to remain. """
        
        #print('\n\nPixArrBySeg =', PixArrBySeg)
        
        if ReducedBySlice:
            """ Replace AllPixArrByRoi and AllF2SindsBySeg with PixArrBySeg  
            and F2SindsBySeg in case further restricting of data is required 
            below: """
            AllPixArrBySeg = deepcopy(PixArrBySeg)
            AllF2SindsBySeg = deepcopy(F2SindsBySeg)
        """ else AllPixArrBySeg and AllF2SindsBySeg remain unchanged. """
        
        #print('\n\nAllPixArrBySeg =', AllPixArrBySeg)
        
        if LogToConsole and ReducedBySlice:
            print('\n   After limiting data to those that relate to slice',
                  f'number {SliceNum}:')#, the F2SindsBySeg =')
            PrintIndsByRoi(F2SindsBySeg)
            print('   PixArrBySeg = ', PixArrBySeg)
            PrintPixArrShapeBySeg(PixArrBySeg)
    
    
    """ Initialise variable that indicates if the SEG data was reduced by
    chosen segment label (FromSegLabel): """
    #ReducedBySeg = False
    
    if SegLabel:
        """ Limit data to those which belong to the chosen segment(s). """
        
        PixArrBySeg = []
        F2SindsBySeg = []
        
        """ Get the segment number(s) whose name matches SegLabel: """
        SegNums = GetRoiNums(Roi=Seg, SearchString=SegLabel)
        
        #print(f'\nSegNums = {SegNums}')
        #print(f'AllF2SindsBySeg = {AllF2SindsBySeg}')
        
        if len(SegNums) != len(AllF2SindsBySeg):
            #ReducedBySeg = True
            
            PixArrBySeg = [AllPixArrBySeg[s] for s in SegNums]
            F2SindsBySeg = [AllF2SindsBySeg[s] for s in SegNums]
            
            if LogToConsole:
                """ Get the names of all ROIs: """
                SegNames = GetRoiLabels(Roi=Seg)
            
                """ Limit the list of segment descriptions to those that belong to 
                the chosen segment(s): """
                SegNames = [SegNames[i] for i in SegNums]
            
                print('\n   After limiting data to those whose segment name',
                      f'matches {SegNames}:')#, the F2SindsBySeg =')
                PrintIndsByRoi(F2SindsBySeg)
                PrintPixArrShapeBySeg(PixArrBySeg)
                
        else:
            PixArrBySeg = deepcopy(AllPixArrBySeg)
            F2SindsBySeg = deepcopy(AllF2SindsBySeg)
        
    else:
        PixArrBySeg = deepcopy(AllPixArrBySeg)
        F2SindsBySeg = deepcopy(AllF2SindsBySeg)
    
    if LogToConsole:
        print('\n   Final outputs of GetSegDataOfInterest():')
        PrintPixArrShapeBySeg(PixArrBySeg)
        PrintIndsByRoi(F2SindsBySeg)
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg


def ProportionOfPixArrInExtent(PixArr, FrameToSliceInds, SrcDicomDir, 
                               TrgDicomDir, LogToConsole=False):
    """
    07/07/21:
        See prop_of_pixarr_in_extent in general_tools.geometry.py
        
    COMMENT 29/01/2021:
        I was using the same image (i.e. same DICOM directory) to generate the
        pixel array (i.e. PixArr2PtsByContour()) and to generate the vertices
        (i.e. GetImageVertices()).  The former had to be the source DICOMs and
        the latter the target DICOMs.
        
    
    Determine what proportion of voxels that make up a 3D pixel array intersect
    with the 3D volume of a 3D image.
    
    Inputs:
    ******
    
    PixArr : Numpy array with dimensions Z x Y x X where Z may be 1
        The pixel array that may represent a mask or labelmap.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    SrcDicomDir : string
        Directory containing the DICOMs that relate to PixArr.
    
    TrgDicomDir : string
        Directory containing the DICOMs that relate to the image grid/extent
        within which the test is for.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    FracProp : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArr2PtsByContour
    from ImageTools import ImportImage, GetImageVertices
    from GeneralTools import IsPointInPolygon
    
    F, R, C = PixArr.shape
    
    # Convert the pixel array to a list (for each frame in PixArr) of a list 
    # (for each point) of a list (for each dimension) of physical coordinates:
    PtsByObjByFrame,\
    CntDataByObjByFrame = PixArr2PtsByContour(PixArr, 
                                              FrameToSliceInds, 
                                              SrcDicomDir)
    
    Image = ImportImage(TrgDicomDir)
    
    Vertices = GetImageVertices(Image)
    
    IsInside = []
    
    for f in range(len(PtsByObjByFrame)):
        for o in range(len(PtsByObjByFrame[f])):
            for pt in PtsByObjByFrame[f][o]:
                IsInside.append(IsPointInPolygon(pt, Vertices))
    
    FracProp = sum(IsInside)/len(IsInside)
    
    if LogToConsole:
        print('\n\n\nFunction ProportionOfPixArrInExtent():')
        print(f'   PixArr has {F}x{R}x{C} (FxRxC)')
        print(f'   FrameToSliceInds = {FrameToSliceInds}')
        print(f'   len(PtsByObjByFrame) = {len(PtsByObjByFrame)}')
        print(f'   len(CntDataByObjByFrame) = {len(CntDataByObjByFrame)}')
        for f in range(len(PtsByObjByFrame)):
            print(f'   Frame {f} has {len(PtsByObjByFrame[f])} objects')
            for o in range(len(PtsByObjByFrame[f])):
                print(f'      Object {o} has {len(PtsByObjByFrame[f][o])} points')
                ##print(f'   len(CntDataByFrame[{i}]) = {len(CntDataByFrame[i])}')
        #print(f'   \nPtsByFrame = {PtsByFrame}')
        #print(f'   \nCntDataByFrame = {CntDataByFrame}')
        
        if False:
            print('\n   The vertices of the image grid:')
            for i in range(len(Vertices)):
                print(f'   {Vertices[i]}')
            
            print('\n   The points from the input indices:')
            for i in range(len(PtsByObjByFrame[f][o])):
                print(f'   {PtsByObjByFrame[f][o][i]}')
                
            print(f'\n   len(IsInside) = {len(IsInside)}')
            print(f'   sum(IsInside) = {sum(IsInside)}')
            
        if FracProp < 1:
            print('\n   The vertices of the image grid:')
            for i in range(len(Vertices)):
                print(f'   {Vertices[i]}')
            N = len(IsInside)
            limit = 10
            if N > limit:
                print(f'\n   The first {limit} points that lie outside of extent:')
                #print([f'\n   {PtsByObjByFrame[0][0][i]}' for i in range(5)])
                for i in range(limit):
                    print(f'   {PtsByObjByFrame[0][0][i]}')
            else:
                print(f'\n   The first {N} points that lie outside of extent:')
                #print([f'\n   {PtsByObjByFrame[0][0][i]}' for i in range(N)])
                for i in range(N):
                    print(f'   {PtsByObjByFrame[0][0][i]}')
    
    return FracProp


def ProportionOfSegsInExtent(PixArrBySeg, F2SindsBySeg, SrcDicomDir, 
                             TrgDicomDir, LogToConsole=False):
    """
    07/07/21:
        See prop_of_segs_in_extent in general_tools.geometry.py
        
    COMMENT 29/01/2021:
        I was using the same image (i.e. same DICOM directory) to generate the
        pixel array (i.e. PixArr2PtsByContour()) and to generate the vertices
        (i.e. GetImageVertices()).  The former had to be the source DICOMs and
        the latter the target DICOMs.
        
    
    Determine what proportion of voxels that make up a list of 3D pixel arrays 
    intersect with the 3D volume of a 3D image.
    
    Inputs:
    ******
    
    PixArrBySeg : list of Numpy arrays
        A list (for each segment) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    
    F2SindsBySeg : List of list of integers
        A list (for each segment) of a list (for each frame) of the slice 
        numbers that correspond to each frame in each pixel array in 
        PixArrBySeg.   
        
    SrcDicomDir : string
        Directory containing the DICOMs that relate to PixArrBySeg.
    
    TrgDicomDir : string
        Directory containing the DICOMs that relate to the image grid/extent
        within which the test is for.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    FracProp : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArrByRoi2PtsByCntByRoi
    from RtsTools import ProportionOfRoisInExtent
    
    # Convert the list of pixel arrays-by-segment to a list of 
    # points-by-contour-by-ROI:
    PtsByCntByRoi, CntDataByCntByRoi,\
    C2SindsByRoi = PixArrByRoi2PtsByCntByRoi(PixArrBySeg, F2SindsBySeg, 
                                             SrcDicomDir, Thresh=0.5)
    
    FracProp = ProportionOfRoisInExtent(PtsByCntByRoi, TrgDicomDir, 
                                        LogToConsole)
    
    return FracProp


def GetFrameNums(Seg, SearchString, DicomDir):
    """
    06/07/21:
        See get_frame_nums in seg_tools.metadata.py
    
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
        
    SearchString : string
        All or part of the Segment Label containing the segmentation of
        interest.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
                            
    Outputs:
    *******
    
    FrameNums : list of integer
        List of the frame numbers that correspond to all frames in the SEG's 
        pixel array that belongs to the segment of interest.  If only one frame 
        corresponds to the segment, the returned list will be of length one 
        (e.g. [3]).
        
    PFFGStoSliceInds : list of integers
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    # Get the DimensionIndexValues (DIVs) and group the DIVs by segment.
    """ Note:  The DIVs are integers that start from 1. """
    DIVs = GetDIVs(Seg)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir)
    
    # Get the PerFrameFunctionalGroupsSequence-to-slice indices
    # (PFFGStoSliceInds) and PFFGStoSliceInds grouped by segment
    # (PFFGStoSliceIndsBySeg):
    """ Note:  The PFFGStoSliceInds are integers that start from 0. """
    PFFGStoSliceInds = GetPFFGStoSliceInds(Seg, SOPuids)
    
    PFFGStoSliceIndsBySeg = GroupListBySegment(ListToGroup=PFFGStoSliceInds, 
                                               DIVs=DIVs)
    
    SegNum = GetRoiNum(Seg, SearchString)
    
    # The PFFGStoSliceInds for the segment of interest:
    PFFGStoSliceInds = PFFGStoSliceIndsBySeg[SegNum]
    
    # The number of segments:
    S = len(PFFGStoSliceIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSliceInds)
    
    n = 0 # initialise frame number counter
    
    FrameNums = []
    
    for s in range(S):
        if s == SegNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of PFFGStoSliceInds[f] in PFFGStoSliceIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                FrameNum = n + PFFGStoSliceIndsBySeg[s].index(PFFGStoSliceInds[f])
                
                FrameNums.append(FrameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSliceIndsBySeg[s])
    
    return FrameNums, PFFGStoSliceInds


def GetRSOPuidsInRIS(Seg):
    """
    06/07/21:
        See get_RSOP_uids_in_RIS in seg_tools.metadata.py
        
    Get the list of ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence
    of a SEG.
    
    Inputs:
    ******
        
    Seg : Pydicom object
        ROI object from a SEG file.
    
        
    Outputs:
    *******
        
    RSOPuids : list of strings
        List of ReferencedSOPInstanceUIDs.
    """
    
    RSOPuids = [] 
    
    sequences = Seg.ReferencedSeriesSequence[0]\
                   .ReferencedInstanceSequence
    
    for sequence in sequences:
        uid = sequence.ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids


def GetRIStoSliceInds(Seg, SOPuids):
    """
    06/07/21:
        See get_RIS_to_slice_inds in seg_tools.metadata.py
        
    Get the slice numbers that correspond to each ReferencedInstanceSequence
    in a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        List of SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
    
    RIStoSliceInds : list of integers
        Slice numbers that correspond to each ReferencedInstanceSequence in 
        Seg.
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    RSOPuids = GetRSOPuidsInRIS(Seg)
    
    RIStoSliceInds = [] 
    
    for Ruid in RSOPuids:
        # Find the matching index of Ruid in SOPuids:
        RIStoSliceInds.append(SOPuids.index(Ruid))
        
    return RIStoSliceInds


def GetRSOPuidsInPFFGS(Seg):
    """
    06/07/21:
        See get_RSOP_uids_in_PFFGS in seg_tools.metadata.py
        
    Get the list of ReferencedSOPInstanceUIDs in the 
    Per-FrameFunctionalGroupsSequence of a SEG.
    
    Inputs:
    ******
        
    Seg : Pydicom object
        ROI object from a SEG file.
    
        
    Outputs:
    *******
    
    RSOPuids : list of strings
        List of ReferencedSOPUIDs.
    """
    
    RSOPuids = [] 
    
    sequences = Seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        uid = sequence.DerivationImageSequence[0]\
                      .SourceImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids


def GetPFFGStoSliceInds(Seg, SOPuids):
    """
    06/07/21:
        See get_PFFGS_to_slice_inds in seg_tools.metadata.py
        
    Get a list of the slice numbers that correspond to each  
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Inputs:
    ******
        
    Seg : Pydcom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
        
    PFFGStoSliceInds : list of integers
        List (for each frame) of the slice numbers that correspond to each 
        Per-FrameFunctionalGroupsSequence in Seg.
    """
    
    # Get the list of ReferencedSOPInstanceUIDs:
    RSOPuids = GetRSOPuidsInPFFGS(Seg)
    
    PFFGStoSliceInds = [] 
    
    for i in range(len(RSOPuids)):
        if RSOPuids[i] in SOPuids:
            # Find the matching index of RSOPuids[i] in SOPuids:
            PFFGStoSliceInds.append(SOPuids.index(RSOPuids[i]))
        else:
            msg = f'ReferencedSOPInstanceUID[{i}], {RSOPuids[i]}, is not in '\
                  + 'the list of SOPInstanceUIDs.'
            
            raise Exception(msg)
        
    return PFFGStoSliceInds


def GetPFFGStoSliceIndsBySeg(Seg, SOPuids):
    """
    06/07/21:
        See get_PFFGS_to_slice_inds_by_seg in seg_tools.metadata.py
        
    Get a list (per segment) of the slice numbers that correspond to each 
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Inputs:
    ******
        
    Seg : Pydcom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
        
    PFFGStoSliceIndsBySeg : list of a list of integers
        List (for each segment) of a list (for each frame) of the slice numbers 
        that correspond to each Per-FrameFunctionalGroupsSequence in Seg.
    """
    
    PFFGStoSliceInds = GetPFFGStoSliceInds(Seg, SOPuids)
    
    DIVs = GetDIVs(Seg)
    
    PFFGStoSliceIndsBySeg = GroupListBySegment(ListToGroup=PFFGStoSliceInds, 
                                               DIVs=DIVs)
    
    return PFFGStoSliceIndsBySeg


def GetDIVs(Seg):
    """
    06/07/21:
        See get_DIVs in seg_tools.metadata.py
        
    Get the DimensionIndexValues in a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
        
    Outputs:
    *******
    
    DIVs : list of lists of integers
        List of DimensionIndexValues in Seg.
        
    """
    
    # Initialise the list of DIVs:
    DIVs = [] 
    
    sequences = Seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        DIV = sequence.FrameContentSequence[0]\
                      .DimensionIndexValues
        
        DIVs.append(DIV)
        
    return DIVs


def GroupListBySegment(ListToGroup, DIVs):
    """
    06/07/21:
        See group_list_by_seg in seg_tools.metadata.py
        
    Group a list of items by common segment/ROI using the DimensionIndexValues.
    
    Inputs:
    ******
    
    ListToGroup : list
        The list to be grouped.
        
    DIVs : list of lists of integers
        List of DimensionIndexValues.
        
        
    Outputs:
    *******
    
    GroupedList : list of lists
        ListToGroup grouped by ROI.
        
        
    Note:
    ----
    
    DIVs is a list of lists of the form:
            
    [ [1, 10], [1, 11], [1, 12], ... [2, 6], [2, 7], [2, 8], ... ]
    
    where each list of index values that belong to a common ROI are separated
    into separate lists.
    
    The goal is to group a list (of anything) by ROI.  In the case of the DIVs,
    that is:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ],
      [ [2, 6], [2, 7], [2, 8], ... ],
      ... 
     ]
      
     
    If there's only one ROI, ListToGroup will still be nested in a list to 
    maintain continuity.  In the case of the DIVs with only one ROI:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ]
     ]
    """
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in DIVs:
        RoiNums.append(DIV[0])
        
    UniqueRoiNums = list(set(RoiNums))
    
    if len(UniqueRoiNums) == 1:
        GroupedList = [ListToGroup]
        
    else:
        GroupedList = []
        
        for RoiNum in UniqueRoiNums:
            ListThisRoi = [] # initialise
            
            for i in range(len(ListToGroup)):
                if DIVs[i][0] == RoiNum:
                    ListThisRoi.append(ListToGroup[i])
        
            GroupedList.append(ListThisRoi)
             
    return GroupedList


def AddToPFFGS(Seg):
    """
    Append the last item in the Per-FrameFunctionalGroupsSequence of the
    SEG to itself (i.e. increase the length of the sequence by one). 
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    
    Outputs:
    *******
    
    NewSeg : Pydicom object
        Modified ROI object with Per-FrameFunctionalGroupsSequence lengthened
        by one.
    """
    
    from copy import deepcopy
    
    # Use Seg as a template for NewSeg: 
    NewSeg = deepcopy(Seg)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = deepcopy(Seg.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    NewSeg.PerFrameFunctionalGroupsSequence.append(last)
             
    return NewSeg


def AddToSS(Seg):
    """
    Append the last item in the SegmentSequence of the SEG Object to itself 
    (i.e. increase the length of the sequence by one). 
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    
    Outputs:
    *******
    
    NewSeg : Pydicom object
        Modified Seg with SegmentSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use Seg as a template for NewSeg: 
    NewSeg = deepcopy(Seg)
             
    # The last item in Segment Sequence:
    last = deepcopy(Seg.SegmentSequence[-1])
    
    # Append to the Segment Sequence:
    NewSeg.SegmentSequence.append(last)
             
    return NewSeg


def ChangeNumOfSegSequences(Seg, Dicoms, NumOfFrames, LogToConsole=False):
    """
    Adjust the number of sequences in a SEG object based so that they are 
    consistent with the number of DICOMs in the series and the number of
    frames in the pixel array.
         
    
    Inputs:
    ******
                              
    Seg : Pydicom object
        SEG object.
        
    Dicoms : List of Pydicom objects
        List of DICOM objects that relate to the SEG object to be created.
    
    NumOfFrames : integer
        The number of frames in the SEG object's pixel array.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
        
    Seg : Pydicom object
        SEG object with modified number of sequences.
    """
    
    
    """ 
    Modify the number of sequences in ReferencedSeriesSequence. 
    The number of sequences must be equal to the number of DICOMs.
    """
    
    from copy import deepcopy
    
    # The number of DICOMs:
    D = len(Dicoms)
    
    # The number of sequences in the ReferencedSeriesSequence:
    R = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)
    
    
    # Add/remove sequences by appending/popping the last sequence N times, 
    # where N is the difference between D and R:
    if D > R:
        for i in range(D - R):
            Seg.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence\
               .append(deepcopy(Seg.ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1]))          
    else:
        for i in range(R - D):
            Seg.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()

    
    """ 
    Modify the number of sequences in SegmentSequence. 
    There must be only one sequence - for the (one) new segment.
    """
    
    S = len(Seg.SegmentSequence)
    
    if S > 1:
        for i in range(S - 1):
            Seg.SegmentSequence.pop()
    
    
    """ 
    Modify the number of sequences in Per-FrameFunctionGroupsSequence. 
    The number of sequences must be equal to NumOfFrames.
    """
    
    # The number of sequences in the Per-FrameFunctionGroupsSequence:
    P = len(Seg.PerFrameFunctionalGroupsSequence)
    
    if NumOfFrames > P:
        for i in range(NumOfFrames - P):
            Seg.PerFrameFunctionalGroupsSequence\
               .append(deepcopy(Seg.PerFrameFunctionalGroupsSequence[-1]))          
    else:
        for i in range(P - NumOfFrames):
            Seg.PerFrameFunctionalGroupsSequence.pop()

    
    if True:#LogToConsole:       
        R2 = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)
        
        S2 = len(Seg.SegmentSequence)
        
        P2 = len(Seg.PerFrameFunctionalGroupsSequence)
        
        print(f'\nThere are {D} DICOMs and there were {R} sequences in',
              f'ReferencedSeriesSequence. Now there are {R2} sequences.')
        
        print(f'\nThere were {S} sequences in SegmentSequence. Now there are',
              f'{S2} sequences.')
        
        print(f'\nThere are {NumOfFrames} frames in the pixel array and there',
              f'were {P} sequences in PerFrameFunctionalGroupsSequence. Now',
              f'there are {P2} sequences.')
        
    
    return Seg


def GetFrameFromPixArr(Seg, DicomDir, SearchString, SliceNum, LogToConsole):
    """
    09/07/21:
        See get_frame_from_pixarr in general_tools.pixarr_ops.py
    
    Extract a single (2D) frame from a 3D pixel array that matches a given
    segment label.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    SearchString : string
        All or part of the Segment label containing the segmentation to be
        copied.
    
    SliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    Frame : Numpy array
        Frame from PixArr.
    """
    
    import numpy as np
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiLabels
    from DicomTools import GetRoiNum
    
    SegLabels = GetRoiLabels(Seg)
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  SegNum is equivalent to SegmentNumber in SegmentSequence. """
    SegNum = GetRoiNum(Seg, SearchString)
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    PFFGStoSliceIndsBySeg = GetPFFGStoSliceIndsBySeg(Seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[SegNum]
    
    if not SliceNum in PFFGStoSliceIndsInSeg:
        msg = f'There is no segmentation in the segment matching the term '\
              + f'"{SearchString}" for slice {SliceNum} in the SEG.'
              
        raise Exception(msg)
    
    # Get the corresponding FrameNum:
    
    n = 0 # initialise frame number counter
    
    for i in range(len(SegLabels)):
        if i == SegNum:
            """ 
            This segment contains the frame to be copied, whose position is
            given by the index of FromSliceNum in 
            SrcPFFGStoSliceIndsBySeg[i] plus any frames that preceeded it 
            (i.e. the frame counter n).
            """
            FrameNum = n + PFFGStoSliceIndsBySeg[i].index(SliceNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSliceIndsBySeg[i]) 
    
    
    if LogToConsole:
        print('\n\nResults of GetFrameFromPixArr:')
        print(f'   SliceNum = {SliceNum} relates to FrameNum = {FrameNum} in',
              'PixelArray')
    
    PixArr = Seg.pixel_array
    
    if len(PixArr.shape) < 3:
        # This is a single-frame PixArr, i.e. shape (R, C). Reshape it to 
        # shape (1, R, C):
        R, C = PixArr.shape
        
        PixArr = np.reshape(PixArr, (1, R, C))
    else:
        # This is a multi-frame PixArr, i.e. shape (AllF, R, C).
        AllF, R, C = PixArr.shape
    
    # Initialise Frame:
    Frame = np.zeros((1, R, C), dtype='uint')
    
    Frame[0] = PixArr[FrameNum]
            
    return Frame


def GetFramesFromPixArr(Seg, DicomDir, SearchString, LogToConsole):
    """
    09/07/21:
        See get_frames_from_pixarr in general_tools.pixarr_ops.py
    
    Extract all frames from a 3D pixel array that match a given segment label.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    SearchString : string
        All or part of the Segment Label of the segment containing the 
        segmentation to be copied.
    
    SliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    Frame : Numpy array
        Frame from PixArr.
    """
    
    import numpy as np
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  SegNum is equivalent to SegmentNumber in SegmentSequence. """
    SegNum = GetRoiNum(Seg, SearchString)
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    PFFGStoSliceIndsBySeg = GetPFFGStoSliceIndsBySeg(Seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[SegNum]
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSliceIndsInSeg)
    
    # Get the number of rows and columns required for Frames:
    PixArr = Seg.pixel_array
    
    if len(PixArr.shape) < 3:
        """ This is a single-frame PixArr, i.e. shape (1, R, C). """
        R, C = PixArr.shape
    else:
        """ This is a multi-frame PixArr, i.e. shape (AllF, R, C). """
        AllF, R, C = PixArr.shape
    
    Frames = np.zeros((F, R, C), dtype='uint')
    
    for f in range(F):
        SliceNum = PFFGStoSliceIndsInSeg[f]
        
        Frame = GetFrameFromPixArr(Seg, DicomDir, SearchString, SliceNum, 
                                   LogToConsole)
        
        Frames[f] = Frame
    
    if LogToConsole:
        print('\n\n***Result from GetFramesFromPixArr:')
        print(f'   There are {F} frames in the segment matching',
              f'"{SearchString}".')
        print(f'   PFFGStoSliceIndsInSeg = {PFFGStoSliceIndsInSeg}')
        print(f'   PixArr.shape = {PixArr.shape}')
        if len(PixArr.shape) > 2:
            [print(f'   np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        else:
            print(f'   np.amax(PixArr) = {np.amax(PixArr)}')
        print(f'   Frames.shape = {Frames.shape}')
        if len(Frames.shape) > 2:
            [print(f'   np.amax(Frames[{i}]) = {np.amax(Frames[i])}') for i in range(Frames.shape[0])]
        else:
            print(f'   np.amax(Frames) = {np.amax(Frames)}')
            
    return Frames



def GetPixArrInSeg(Seg, DicomDir, SearchString):
    """
    06/07/21:
        See get_pixarr_in_seg in seg_tools.pixarrs
        
    Get the frames in a SEG's pixel array that belong to a segment with
    specified label.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    SearchString : string
        All or part of the Segment Label containing the segment of interest.
                       
                            
    Outputs:
    *******
    
    PixArrInSeg : Numpy array
        Sub-array of the SEG's pixel array that contains the frames that belong
        to the segment of interest. 
        
    PFFGStoSliceIndsInSeg : List of integers
        List (for each segmentation) of the slice numbers that correspond to 
        each frame in PixArr.
    """
    
    import numpy as np
    
    AllPixArr = Seg.pixel_array
    
    AllF, R, C = AllPixArr.shape
    
    # Get the frame numbers of the SEG's pixel array that correspond to the 
    # segment of interest, and the corresponding Per-frame Functional Groups
    # Sequence-to-slice indices:
    FrameNumsInSeg, PFFGStoSliceIndsInSeg = GetFrameNums(Seg, SearchString, 
                                                         DicomDir)
    
    F = len(FrameNumsInSeg)
    
    PixArrInSeg = np.zeros((F, R, C), dtype='uint')
    
    for i in range(F):  
        PixArrInSeg[i] = AllPixArr[FrameNumsInSeg[i]]
    
        
    return PixArrInSeg, PFFGStoSliceIndsInSeg



def GetPixArrBySeg(Seg, DicomDir, LogToConsole=False):
    """
    06/07/21:
        See get_pixarr_by_seg in seg_tools.pixarrs
        
    Get a list of pixel arrays grouped by segment in a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
        
    F2SindsBySeg : List of a list of integers
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in PixArrBySeg.
    """
    
    import numpy as np
    from DicomTools import GetDicomSOPuids
    
    if LogToConsole:
        from GeneralTools import PrintIndsByRoi, PrintPixArrShapeBySeg
        
        print('\n\n', '-'*120)
        print('Running of GetPixArrBySeg():')
        print('\n\n', '-'*120)
        
    AllPixArr = Seg.pixel_array
    
    """ If the Seg has a single frame it will have shape (NumOfRows, NumOfCols).
    For compatibility with multi-framed pixel arrays, convert single-framed
    pixel arrays to shape (1, NumOfRows, NumOfCols). """
    
    if len(AllPixArr.shape) < 3:
        R, C = AllPixArr.shape
        
        AllPixArr = np.reshape(AllPixArr, (1, R, C))
    
    AllF, R, C = AllPixArr.shape
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    F2SindsBySeg = GetPFFGStoSliceIndsBySeg(Seg, SOPuids)
    
    Nsegs = len(F2SindsBySeg)
    
    if LogToConsole:
        print(f'   AllPixArr.shape = {AllPixArr.shape}')
        #print('   F2SindsBySeg =')
        PrintIndsByRoi(F2SindsBySeg)
        
    PixArrBySeg = []
    
    j = 0 # total frame counter for all frames in AllPixArr
    
    for s in range(Nsegs):
        # The number of frames in AllPixArr that belong to this segment:
        F = len(F2SindsBySeg[s])
        
        # Initialise the pixel array for this segment:
        PixArr = np.zeros((F, R, C), dtype='uint')
        
        if F > 0:
            for i in range(F):
                PixArr[i] = AllPixArr[j]
                
                j += 1 # increment the total frame counter
        
        PixArrBySeg.append(PixArr)
    
    if LogToConsole:
        PrintPixArrShapeBySeg(PixArrBySeg)
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg


def AddCopiedPixArr(OrigPixArr, OrigPFFGStoSliceInds, 
                    PixArrToAdd, PFFGStoSliceIndsToAdd):
    """
    Add frame(s) from a pixel array to an existing pixel array.
       
    
    Inputs:
    ******
    
    OrigPixArr : Numpy array
        A FxRxC (Frames x Rows x Columns) array of all frames that belong to
        the segment of interest. 
        
    OrigPFFGStoSliceInds : list of integers
        List (for each segmentation) of slice numbers that correspond to each 
        frame in OrigPixArr.
        
    PixArrToAdd : Numpy array
        A FxRxC (Frames x Rows x Columns) array of all frames that originated
        from a single segmentation.  F will be one if a single segmentation is
        to be copied, but may be greater than one, if, for example, the 
        segmentation to be copied was resampled to a number of segmentations.
        
    PFFGStoSliceIndsToAdd : list of integers
        List (for each segmentation) of slice numbers that correspond to each 
        frame in PixArrToAdd.
        
        
    Outputs:
    *******
    
    NewPixArr : Numpy array
        A FxRxC (Frames x Rows x Columns) array of all frames in 
        OrigPixArrInSeg and PixArrToAddToSeg. Note that F will not necessarily
        be the sum of the number of frames in OrigPixArrInSeg and 
        PixArrToAddToSeg.  Some frames in both pixel arrays may correspond
        to the same slice number, in which case, any duplicate frame(s) in
        OrigPixArrInSeg will be over-written with the frame(s) in
        PixArrToAddToSeg that correspond to the same slice number.
        
    NewPFFGStoSliceInds : list of integers
        Modified list (for each segmentation) of slice numbers that correspond   
        to each frame in NewPixArr.
    """
    
    import numpy as np
    from GeneralTools import UniqueItems
    
    # Combine OrigPFFGStoSliceInds and PFFGStoSliceIndsToAdd:
    NewPFFGStoSliceInds = OrigPFFGStoSliceInds + PFFGStoSliceIndsToAdd
    
    # Remove any duplicate slice numbers:
    NewPFFGStoSliceInds = UniqueItems(NewPFFGStoSliceInds)
    
    F = len(NewPFFGStoSliceInds)
    
    OrigF, R, C = OrigPixArr.shape
    
    NewPixArr = np.zeros((F, R, C), dtype='uint')
    
    for FrameNum in range(F):
        # The corresponding slice number for this frame:
        SliceNum = NewPFFGStoSliceInds[FrameNum]
        
        if SliceNum in PFFGStoSliceIndsToAdd:
            """
            This frame is to be over-written by the corresponding frame in
            PixArrToAdd.
            """
            
            ind = PFFGStoSliceIndsToAdd.index(SliceNum)
            
            NewPixArr[FrameNum] = PixArrToAdd[ind]
            
        else:
            """
            This existing frame (in OrigPixArr) is to be preserved.
            """
            
            ind = OrigPFFGStoSliceInds.index(SliceNum)
            
            NewPixArr[FrameNum] = OrigPixArr[ind]
    
    
    return NewPixArr, NewPFFGStoSliceInds


def CreateSeg(SrcSegFpath, TrgSegFpath, TrgSegLabel, TrgPixArrBySeg, 
              TrgF2SindsBySeg, TrgDicomDir, SrcSegLabel, AddTxtToSegLabel='', 
              LogToConsole=False):
    """
    09/07/21:
        See create_seg in seg_tools.create.py
    
    Create an SEG object for the target dataset.
         
    
    Inputs:
    ******
    
    SrcSegFpath : string
        Filepath of the Source SEG file. If TrgSegFpath = None, the Source SEG
        will be used as a template for the new Target SEG (NewTrgSeg).
    
    TrgSegFpath : string or None
        If TrgSegFpath != None, the Target SEG will be modified to create 
        NewTrgSeg.  If TrgSegFpath = None, the Source SEG will be used as a 
        template.
        
    TrgSegLabel : string or None
        All or part of the SegmentLabel of the destination segment.

    TrgPixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment for Target.  The pixel array is a sub-array of 
        the (entire) SEG's pixel array. 
        
    TrgF2SindsBySeg : List of a list of integers
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in TrgPixArrBySeg.
                    
    TrgDicomDir : string
        Directory containing the Target DICOMs.
    
    SrcSegLabel : string
        All or part of the Source SegmentLabel of the segment containing the 
        segmentations(s) that were copied.
    
    AddTxtToSegLabel : string (optional; '' by default)
        String to be added to SeriesDescription for the new SEG. 
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    NewTrgSeg : Pydicom object
        New SEG object.
        
    
    Notes:
    *****
    
    Performance of adding/removing sequences using append()/pop():
        
    - Took 168.4 ms to add 10000 sequences to ReferencedInstanceSequence 
    (16.8 us/sequence)
    - Took 81.7 ms to remove 10000 sequences to ReferencedInstanceSequence 
    (8.2 us/sequence)
    
    - Took 92.4 ms to add 10000 sequences to PerFrameFunctionalGroupsSequence 
    (9.2 us/sequence)
    -Took 45.8 ms to remove 10000 sequences to PerFrameFunctionalGroupsSequence 
    (4.6 us/sequence)
    
    19/02:
        After avoiding the use of deepcopy to increase sequences in RtsSeg()
        I've had to re-instate them since the tag values were being duplicated 
        despite calls to modify them.  For some reason the removal of deepcopy
        from CreateSeg() hasn't resulted in the unexpected behaviours found 
        with CreateRts().  Nevertheless deepcopy was re-instated in CreateSeg()
        as a precaution.  Oddly the re-introduction of deepcopy doesn't seem to
        have lengthened the time to run CreateSeg() or to perform error 
        checking.
    """
    
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    #import time
    import datetime
    from copy import deepcopy
    import numpy as np
    from GeneralTools import FlattenList, UniqueItems, ReduceListOfStringFloatsTo16
    from DicomTools import ImportDicoms, GetRoiNums, GetRoiLabels
    from ImageTools import GetImageAttributes
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of CreateSeg():')
        print('\n\n', '-'*120)
        print(f'   TrgF2SindsBySeg = {TrgF2SindsBySeg}')
        print(f'   len(TrgPixArrBySeg) = {len(TrgPixArrBySeg)}')
    
    SrcSeg = dcmread(SrcSegFpath)
    
    """ Use TrgSeg or SrcSeg as a template for NewTrgSeg. """
    if TrgSegFpath:
        NewTrgSeg = dcmread(TrgSegFpath)
    else:
        NewTrgSeg = deepcopy(SrcSeg)
    
    #print(f'\nTrgDicomDir = {TrgDicomDir}')
    
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=TrgDicomDir, Package='sitk')
    
    NumOfFramesBySeg = []
    for s in range(len(TrgPixArrBySeg)):
        NumOfFramesBySeg.append(TrgPixArrBySeg[s].shape[0])
    
    #print(f'\n\n\nNumOfFramesBySeg = {NumOfFramesBySeg}')
    
    NumOfFrames = sum(NumOfFramesBySeg)
    
    #NumOfSegs = len(TrgF2SindsBySeg)
    
    """ The list of segment labels in SrcSeg: """
    AllSrcSegLabels = GetRoiLabels(SrcSeg)
    
    if SrcSegLabel:
        """ The list of segment numbers for the segments of interest: """
        SrcSegNums = GetRoiNums(SrcSeg, SrcSegLabel)
    else:
        SrcSegNums = list(range(len(AllSrcSegLabels)))
    
    UniqueTrgF2Sinds = UniqueItems(Items=TrgF2SindsBySeg, IgnoreZero=False, 
                                   MaintainOrder=True)
    
    FlattenedTrgF2Sinds = FlattenList(TrgF2SindsBySeg)
    
    if LogToConsole:
        print(f'   SrcSegNums = {SrcSegNums}')
        print(f'   NumOfFramesBySeg = {NumOfFramesBySeg}')
        print(f'   NumOfFrames = {NumOfFrames}')
    
    
    """ Generate a new SOPInstanceUID. """
    NewSOPuid = generate_uid()
    NewTrgSeg.SOPInstanceUID = deepcopy(NewSOPuid)
    NewTrgSeg.file_meta.MediaStorageSOPInstanceUID = deepcopy(NewSOPuid)
    
    """ Generate a new SeriesInstanceUID. """
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    #NewDate = time.strftime("%Y%m%d", time.gmtime())
    #NewTime = time.strftime("%H%M%S", time.gmtime())
    
    TimeNow = datetime.datetime.now()
    NewDate = TimeNow.strftime('%Y%m%d')
    NewTime = TimeNow.strftime('%H%M%S.%f')
    
    """ If TrgSegFpath != None, some tags will not need to be replaced.
    If TrgSegFpath = None, use the corresponding values in the first Target 
    DICOM. """
    
    if TrgSegFpath == None:
        NewTrgSeg.StudyDate = TrgDicoms[0].StudyDate
        try:
            NewTrgSeg.SeriesDate = TrgDicoms[0].SeriesDate
        except AttributeError:
            pass
        #NewTrgSeg.ContentDate = TrgDicoms[0].ContentDate
        NewTrgSeg.ContentDate = NewDate
        NewTrgSeg.StudyTime = TrgDicoms[0].StudyTime
        try:
            NewTrgSeg.SeriesTime = TrgDicoms[0].SeriesTime
        except AttributeError:
            pass
        #NewTrgSeg.ContentTime = TrgDicoms[0].ContentTime
        NewTrgSeg.ContentTime = NewTime
        NewTrgSeg.Manufacturer = TrgDicoms[0].Manufacturer
        try:
            NewTrgSeg.SeriesDescription += AddTxtToSegLabel
        except AttributeError:
            NewTrgSeg.SeriesDescription = AddTxtToSegLabel
        NewTrgSeg.PatientName = TrgDicoms[0].PatientName
        NewTrgSeg.PatientID = TrgDicoms[0].PatientID
        NewTrgSeg.PatientBirthDate = TrgDicoms[0].PatientBirthDate
        NewTrgSeg.PatientSex = TrgDicoms[0].PatientSex
        try:
            NewTrgSeg.PatientAge = TrgDicoms[0].PatientAge
        except AttributeError:
            pass
        NewTrgSeg.StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
        NewTrgSeg.StudyID = TrgDicoms[0].StudyID
        NewTrgSeg.SeriesNumber = TrgDicoms[0].SeriesNumber
        NewTrgSeg.InstanceNumber = TrgDicoms[0].InstanceNumber
        NewTrgSeg.FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
        NewTrgSeg.PositionReferenceIndicator = TrgDicoms[0].PositionReferenceIndicator
        NewTrgSeg.Rows = TrgDicoms[0].Rows
        NewTrgSeg.Columns = TrgDicoms[0].Columns
        
        #IOP = TrgDicoms[0].ImageOrientationPatient # 01/09/21
        IOP = [str(item) for item in TrgDicoms[0].ImageOrientationPatient] # 01/09/21
        
        """
        print(f'\n\n\nIOP = {IOP}')
        print(f'type(IOP) = {type(IOP)}')
        print(f'type(IOP[0]) = {type(IOP[0])}')
        print(f'IOP[0] = {IOP[0]}')
        print(f'len(IOP[0]) = {len(IOP[0])}')
        print('\n\n\n')
        """
        
        # Ensure that the character length is limited to 16:
        IOP = ReduceListOfStringFloatsTo16(IOP)
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PlaneOrientationSequence[0]\
                 .ImageOrientationPatient = IOP
           
        """ 
        Notes:
            
        1. SliceThickness in PixelMeasuresSequence appears to be the z-Spacing 
        rather than SliceThickness from the DICOM metadata. 
        
        2. The character limit of VR DS is 16 characters.
        """
        Zspacing = f"{Spacings[2]}"
        
        #if len(Zspacing) > 16:
        #    Zspacing = Zspacing[:16]
        
        Zspacing = ReduceListOfStringFloatsTo16([Zspacing])[0]
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .SliceThickness = deepcopy(Zspacing)
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .SpacingBetweenSlices = deepcopy(Zspacing)
    
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .PixelSpacing = TrgDicoms[0].PixelSpacing
    
    
    #""" Modify SeriesDescription. <-- this was done above """
    #
    #if len(SegNums) == 1:
    #    """ There is only one segment. The new Target SeriesDescription will be
    #    the Source SeriesDescription with the Source SegmentLabel that was 
    #    copied + AddTxtToSegLabel. """
    #    NewSD = SrcSeg.SeriesDescription + '_' \
    #            + SrcSeg.SegmentSequence[SegNums[0]].SegmentLabel \
    #            + AddTxtToSegLabel
    #            
    #    NewTrgSeg.SeriesDescription = NewSD
    #else:
        
    
    
    """ Modify NumberOfFrames. """
    
    NewTrgSeg.NumberOfFrames = f"{len(FlattenedTrgF2Sinds)}"
    
    
    """ Modify ReferencedSeriesSequence. """
    
    NewTrgSeg.ReferencedSeriesSequence[0]\
             .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    """ Modify ReferencedInstanceSequence in ReferencedSeriesSequence. """
        
    for i in range(len(UniqueTrgF2Sinds)):
        """ The slice index s determines the SOPInstanceUID to reference. """
        s = UniqueTrgF2Sinds[i]
        
        """ The number of sequences in ReferencedSeriesSequence: """
        N = len(NewTrgSeg.ReferencedSeriesSequence[0]\
                         .ReferencedInstanceSequence)
        
        if i > N - 1:
            """ Increase the sequence by one. """
            #NewTrgSeg.ReferencedSeriesSequence[0]\
            #         .ReferencedInstanceSequence.append(NewTrgSeg.ReferencedSeriesSequence[0]\
            #                                                     .ReferencedInstanceSequence[-1])
                     
            LastItem = deepcopy(NewTrgSeg.ReferencedSeriesSequence[0]\
                                         .ReferencedInstanceSequence[-1])
            
            NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence.append(LastItem)
        
        """ Update the sequence. """
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPClassUID = TrgDicoms[s].SOPClassUID
        
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
    
    """ Check if there are more sequences in ReferencedInstanceSequence than 
    required. """
    N = len(NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    if N > len(UniqueTrgF2Sinds):
        for i in range(N - len(UniqueTrgF2Sinds)):
            """ Remove the last sequence. """
            NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence.pop()
    
    
    
    """ Modify SegmentSequence. """
    
    """ Loop through all SrcSegNums. """
    for s in range(len(SrcSegNums)):
        SegNum = SrcSegNums[s]
        
        """ Replace the s^th SegmentSequence in NewTrgSeg with the SegNum^th 
        SegmentSequence in SrcSeg."""
        NewTrgSeg.SegmentSequence[s] = deepcopy(SrcSeg.SegmentSequence[SegNum])
        
        NewTrgSeg.SegmentSequence[s].SegmentNumber = int(s+1)
        
        NewTrgSeg.SegmentSequence[s]\
                 .SegmentLabel = SrcSeg.SegmentSequence[SegNum].SegmentLabel
        
    """ Check if there are more sequences in SegmentSequence than required. """
    N = len(NewTrgSeg.SegmentSequence)
    
    if N > len(SrcSegNums):
        for i in range(N - len(SrcSegNums)):
            """ Remove the last sequence: """
            NewTrgSeg.SegmentSequence.pop()
    
    
    """ Modify PerFrameFunctionalGroupsSequence. """
    
    """ Get the ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence."""
    RefSOPsInRIS = GetRSOPuidsInRIS(NewTrgSeg)
    
    i = 0 # total frame counter (for all segments)
    
    for s in range(len(TrgF2SindsBySeg)):
        for f in range(len(TrgF2SindsBySeg[s])):  
            N = len(NewTrgSeg.PerFrameFunctionalGroupsSequence)
            
            if i > N - 1:
                """ Increase the sequence by one. """
                #NewTrgSeg.PerFrameFunctionalGroupsSequence\
                #         .append(NewTrgSeg.PerFrameFunctionalGroupsSequence[-1])
                         
                LastItem = deepcopy(NewTrgSeg.PerFrameFunctionalGroupsSequence[-1])
                
                NewTrgSeg.PerFrameFunctionalGroupsSequence\
                         .append(LastItem)
            
            """ The DICOM slice number: """
            d = TrgF2SindsBySeg[s][f]
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .DerivationImageSequence[0]\
                     .SourceImageSequence[0]\
                     .ReferencedSOPClassUID = TrgDicoms[d].SOPClassUID
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .DerivationImageSequence[0]\
                     .SourceImageSequence[0]\
                     .ReferencedSOPInstanceUID = TrgDicoms[d].SOPInstanceUID
            
            SOPuid = TrgDicoms[d].SOPInstanceUID
            
            ind = RefSOPsInRIS.index(SOPuid)
            
            """ Note: While s (segment number) and ind (the index of SOPuid
            within the ReferencedSOPInstanceUIDs in ReferencedInstanceSequence)
            are are both 0-indexed, DimensionIndexValues are 1-indexed. """
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .FrameContentSequence[0]\
                     .DimensionIndexValues = [s + 1, ind + 1]
            
            #IPP = TrgDicoms[d].ImagePositionPatient # 01/09/21
            IPP = [str(item) for item in TrgDicoms[d].ImagePositionPatient] # 01/09/21
            
            # Ensure that the characters are limited to 16:
            IPP = ReduceListOfStringFloatsTo16(IPP)
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .PlanePositionSequence[0]\
                     .ImagePositionPatient = IPP
               
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .SegmentIdentificationSequence[0]\
                     .ReferencedSegmentNumber = s + 1
            
            i += 1 # increment the total frame count
    
    """ Check if there are more sequences in PerFrameFunctionalGroupsSequence 
    than required. """
    N = len(NewTrgSeg.PerFrameFunctionalGroupsSequence)
    
    if N > len(FlattenedTrgF2Sinds):
        for i in range(N - len(FlattenedTrgF2Sinds)):
            """ Remove the last sequence: """
            NewTrgSeg.PerFrameFunctionalGroupsSequence.pop()        
    
    
    """ Modify PixelData. """
    
    """ Start by copying the first (or only) frame: """
    PixArr = deepcopy(TrgPixArrBySeg[0])
    
    if len(TrgPixArrBySeg) > 1:
        """ If there are more than one frame, vertically stack all pixel arrays 
        into a single pixel array. """
        for s in range(1, len(TrgPixArrBySeg)):
            PixArr = np.vstack((PixArr, TrgPixArrBySeg[s]))
    
    """ Note: The following doesn't work for binary arrays:
    NewTrgSeg.PixelData = PixArr.tobytes()
    """
    
    """ Ravel (to 1D array) PixArr and pack bits (see
    https://github.com/pydicom/pydicom/issues/1230). """
    packed = pack_bits(PixArr.ravel())
    
    """ Convert PixArr to bytes. """
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    
    if LogToConsole:
        print('-'*120)
        
    return NewTrgSeg



def ErrorCheckSeg(Seg, DicomDir, LogToConsole=False):
    """
    Check a SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the SEG.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Seg.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe any errors that are found.  If no 
        errors are found an empty list ([]) will be returned.
        
    Nerrors : integer
        The number of errors found.
    """
    
    #import importlib
    #import GeneralTools
    #importlib.reload(GeneralTools)
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from DicomTools import ImportDicom
    from ImageTools import GetImageAttributes
    from GeneralTools import AreItemsEqualToWithinEpsilon
    from GeneralTools import AreListsEqualToWithinEpsilon
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of ErrorCheckSeg():')
        print('\n\n', '-'*120)
        
    # The maximum allowed error between two lists/items:
    epsilon = 1e-05
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    Dicom = ImportDicom(DicomDir)
    
    #Size, Spacings, ST,\
    #IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=DicomDir, Package='pydicom')
    
    PixArr = Seg.pixel_array
    
    """ The number of frames, rows and columns in the SEG's pixel array: """
    if len(PixArr.shape) > 2:
        NumOfFrames, NumOfRows, NumOfCols = PixArr.shape
    else:
        NumOfRows, NumOfCols = PixArr.shape
        NumOfFrames = 1
    
    #NumOfDicoms = len(SOPuids)
    
    if LogToConsole:
        print(f'NumOfFrames = {NumOfFrames}')
    
    LogList = []
    Nerrors = 0
    
    
    """Determine whether MediaStorageSOPInstanceUID matches SOPInstanceUID."""
    MSSOPuid = Seg.file_meta.MediaStorageSOPInstanceUID
    SOPuid = Seg.SOPInstanceUID
    
    if MSSOPuid == SOPuid:
        msg = f'INFO:  MediaStorageSOPInstanceUID {MSSOPuid} matches '\
              + f'SOPInstanceUID {SOPuid}.\n'
    else:
        msg = f'ERROR:  MediaStorageSOPInstanceUID {MSSOPuid} does not '\
              + f'match SOPInstanceUID {SOPuid}.\n'
        Nerrors = Nerrors + 1 
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs do not
    match the SOPInstanceUIDs."""        
    
    RIS = Seg.ReferencedSeriesSequence[0]\
             .ReferencedInstanceSequence
    
    RefSOPuidsInRIS = [RIS[i].ReferencedSOPInstanceUID for i in range(len(RIS))]

    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInRIS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
    if Inds:
        for i in range(len(Inds)):
            uid = RefSOPuidsInRIS[Inds[i]]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + f'ReferencedInstanceSequence[{Inds[i]}] does not match '\
                  + 'any SOPInstanceUID.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
            if LogToConsole:
                print(msg)
    else:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs in '\
                  + 'ReferencedInstanceSequence matches the '\
                  + 'SOPInstanceUIDs.\n'
        LogList.append(msg)
        if LogToConsole:
            print(msg)
    
                           
    if Seg.ReferencedSeriesSequence[0].SeriesInstanceUID == Dicom.SeriesInstanceUID:
        msg = 'INFO:  SEG SeriesInstanceUID in ReferencedSeriesSequence, '\
              + f'{Seg.ReferencedSeriesSequence[0].SeriesInstanceUID}, '\
              + 'matches DICOM SeriesInstanceUID, '\
              + f'{Dicom.SeriesInstanceUID}.\n'
    else:
        msg = 'ERROR:  SEG SeriesInstanceUID in ReferencedSeriesSequence, '\
              + f'{Seg.ReferencedSeriesSequence[0].SeriesInstanceUID}, does '\
              + 'not match DICOM SeriesInstanceUID '\
              + f'{Dicom.SeriesInstanceUID}.\n'
        Nerrors = Nerrors + 1 
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
       
    if Seg.PatientName == Dicom.PatientName:
        msg = f'INFO:  SEG PatientName, {Seg.PatientName}, matches DICOM '\
              + f'PatientName, {Dicom.PatientName}.\n'
    else:
        msg = f'ERROR:  SEG PatientName, {Seg.PatientName}, does not match '\
              + f'DICOM PatientName, {Dicom.PatientName}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.PatientID == Dicom.PatientID:
        msg = f'INFO:  SEG PatientID, {Seg.PatientID}, matches the '\
              + f'DICOM PatientID, {Dicom.PatientID}.\n'
    else:
        msg = f'ERROR:  SEG PatientID, {Seg.PatientID}, does not match the '\
              + f'DICOM PatientID, {Dicom.PatientID}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.StudyInstanceUID == Dicom.StudyInstanceUID:
        msg = f'INFO:  SEG StudyInstanceUID, {Seg.StudyInstanceUID}, '\
              + 'matches the DICOM StudyInstanceUID, '\
              + f'{Dicom.StudyInstanceUID}.\n'
    else:
        msg = f'ERROR:  SEG StudyInstanceUID, {Seg.StudyInstanceUID}, does '\
              + 'not match the DICOM StudyInstanceUID, '\
              + f'{Dicom.StudyInstanceUID}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    #if Seg.StudyID == Dicom.StudyID:
    #    msg = f'INFO:  The SEG StudyID {Seg.StudyID} matches the'\
    #          + f' DICOM StudyID {Dicom.StudyID}.'
    #else:
    #    msg = f'ERROR:  The SEG StudyID {Seg.StudyID} does not match the'\
    #          + f' DICOM StudyID {Dicom.StudyID}.'
    #
    #    Nerrors = Nerrors + 1
    #    
    #LogList.append(msg)
    #    
    #if LogToConsole:
    #    print(msg)
            
    if Seg.FrameOfReferenceUID == Dicom.FrameOfReferenceUID:
        msg = f'INFO:  SEG FrameOfReferenceUID, {Seg.FrameOfReferenceUID}, '\
              + 'matches DICOM FrameOfReferenceUID, '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
    else:
        msg = f'ERROR:  SEG FrameOfReferenceUID, {Seg.FrameOfReferenceUID}, '\
              + 'does not match DICOM FrameOfReferenceUID, '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if int(Seg.NumberOfFrames) == NumOfFrames:
        msg = f'INFO:  NumberOfFrames, {Seg.NumberOfFrames}, matches the '\
              + f'number of frames in PixelData, {NumOfFrames}.\n'
    else:
        msg = f'ERROR:  NumberOfFrames, {Seg.NumberOfFrames}, does not '\
              + f'match the number of frames in PixelData, {NumOfFrames}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.Rows == Dicom.Rows:
        msg = f'INFO:  SEG Rows, {Seg.Rows}, matches DICOM Rows, '\
              + f'{Dicom.Rows}.\n'
    else:
        msg = f'ERROR:  SEG Rows, {Seg.Rows}, does not match DICOM Rows, '\
              + f'{Dicom.Rows}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.Columns == Dicom.Columns:
        msg = f'INFO:  SEG Columns, {Seg.Columns}, matches DICOM '\
              + f'Columns, {Dicom.Columns}.\n'
    else:
        msg = f'ERROR:  SEG Columns, {Seg.Columns}, does not match DICOM '\
              + f'Columns, {Dicom.Columns}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """ This is already checked below
    if Seg.Rows == NumOfRows:
        msg = f'INFO:  SEG Rows, {Seg.Rows}, matches the number of rows in '\
              + f'PixelData, {NumOfRows}.\n'
    else:
        msg = f'ERROR:  SEG Rows, {Seg.Rows}, does not match the number of '\
              + f'rows in PixelData, {NumOfRows}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
        
    if Seg.Columns == NumOfCols:
        msg = f'INFO:  SEG Columns, {Seg.Columns}, matches the number of '\
              + f'columns in PixelData, {NumOfCols}.\n'
    else:
        msg = f'ERROR:  SEG Columns, {Seg.Columns}, does not match the '\
              + f'number of columns in PixelData, {NumOfCols}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    """
    
    
    """Keep track of the number of sequences in SegmentSequence."""
    NumSS = len(Seg.SegmentSequence)
    
    """Verify that the SegmentNumber in each SegmentSequence increments from 
    1."""
    SS = deepcopy(Seg.SegmentSequence)
    
    for i in range(NumSS):
        N = int(SS[i].SegmentNumber)
        
        #print(f'\n\nN = {N}, type(N) = {type(N)}')
        
        if N == i + 1:
            msg = f'INFO:  SegmentNumber in SegmentSequence[{i}] is {N}.\n'
            LogList.append(msg)
        else:
            msg = f'ERROR:  SegmentNumber in SegmentSequence[{i}] is {N}.'\
                  + f' It should be {i + 1}.\n'
            Nerrors = Nerrors + 1
            LogList.append(msg)

    if LogToConsole:
        print(msg)
    
    
    """Check various tags in SharedFunctionalGroupsSequence."""
    SegIOP = Seg.SharedFunctionalGroupsSequence[0]\
                .PlaneOrientationSequence[0]\
                .ImageOrientationPatient
    
    SegIOP = [float(item) for item in SegIOP]
    
    DcmIOP = [float(item) for item in Dicom.ImageOrientationPatient]
    
    if SegIOP == DcmIOP:
        msg = 'INFO:  SEG ImageOrientationPatient in '\
              + f'SharedFunctionalGroupsSequence, {SegIOP} matches '\
              + f'DICOM ImageOrientationPatient, {DcmIOP}.\n'
    else:
        if AreListsEqualToWithinEpsilon(SegIOP, DcmIOP, epsilon):
            msg = 'INFO:  SEG ImageOrientationPatient in '\
                  + f'SharedFunctionalGroupsSequence, {SegIOP}, is within '\
                  + f'{epsilon} of DICOM ImageOrientationPatient, {DcmIOP}.\n'
        else:
            msg = 'ERROR:  SEG ImageOrientationPatient in '\
                  + f'SharedFunctionalGroupsSequence, {SegIOP}, is not within'\
                  + f' {epsilon} of DICOM ImageOrientationPatient, {DcmIOP}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
       
    
    SegST = Seg.SharedFunctionalGroupsSequence[0]\
               .PixelMeasuresSequence[0]\
               .SliceThickness
    
    """ The SliceThickness appears to be the z-Spacing rather than the
    SliceThickness from the DICOM metadata. """
    
    if float(SegST) == Spacings[2]:
        msg = 'INFO:  SEG SliceThickness in SharedFunctionalGroupsSequence,'\
              + f' {SegST}, matches the slice thickness calculated '\
              + 'from DICOM ImagePositionPatient and '\
              + f'ImageOrientationPatient, {Spacings[2]}.\n'
    else:
        if AreItemsEqualToWithinEpsilon(float(SegST), Spacings[2], epsilon):
            msg = 'INFO:  SEG SliceThickness in '\
                  + f'SharedFunctionalGroupsSequence, {SegST}, is within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM'\
                  + ' ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
        else:
            print(f'\n{float(SegST) - Spacings[2]}')
            
            msg = 'ERROR:  SEG SliceThickness in '\
                  + f'SharedFunctionalGroupsSequence, {SegST}, is not within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    SegSBS = Seg.SharedFunctionalGroupsSequence[0]\
                .PixelMeasuresSequence[0]\
                .SpacingBetweenSlices
    
    if float(SegSBS) == Spacings[2]:
        msg = 'INFO:  SEG SpacingBetweenSlices in '\
              + f'SharedFunctionalGroupsSequence, {SegSBS}, matches the slice'\
              + ' thickness calculated from DICOM ImagePositionPatient and '\
              + f'ImageOrientationPatient, {Spacings[2]}.\n'
    else:
        if AreItemsEqualToWithinEpsilon(float(SegSBS), Spacings[2], epsilon):
            msg = 'INFO:  SEG SpacingBetweenSlices in '\
                  + f'SharedFunctionalGroupsSequence, {SegSBS}, is within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
        else:
            #print(f'\n{float(SegSBS) - Spacings[2]}')
            msg = 'ERROR:  SEG SpacingBetweenSlices in '\
                  + f'SharedFunctionalGroupsSequence, {SegSBS}, is not within'\
                  + f' {epsilon} of the slice thickness calculated from DICOM'\
                  + ' ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)

    
    SegPS = Seg.SharedFunctionalGroupsSequence[0]\
               .PixelMeasuresSequence[0]\
               .PixelSpacing
    
    SegPS = [float(item) for item in SegPS]
    
    DcmPS = [float(item) for item in Dicom.PixelSpacing]
    
    if SegPS == DcmPS:
        msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence, '\
              + f'{SegPS}, matches the DICOM PixelSpacing, {DcmPS}.\n'
    else:
        if AreListsEqualToWithinEpsilon(SegPS, DcmPS, epsilon):
            msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence,'\
                  + f' {SegPS} is within {epsilon} of the DICOM PixelSpacing,'\
                  + f' {DcmPS}.\n'
        else:
            msg = 'ERROR:  SEG PixelSpacing in SharedFunctionalGroupsSequence'\
                  + f', {SegPS}, is not within {epsilon} of the DICOM '\
                  + f'PixelSpacing, {DcmPS}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Verify that the number of sequences in PerFrameFunctionalGroupsSequence
    is equal to NumberOfFrames."""
    PFFGS = Seg.PerFrameFunctionalGroupsSequence
    
    if len(PFFGS) == int(Seg.NumberOfFrames):
        msg = 'INFO:  The number of sequences in '\
              + f'PerFrameFunctionGroupsSequence, {len(PFFGS)}, '\
              + f'matches NumberOfFrames {Seg.NumberOfFrames}.\n'
    else:
        msg = 'ERROR:  The number of sequences in '\
              + f'PerFrameFunctionGroupsSequence, {len(PFFGS)}, does not '\
              + f'match NumberOfFrames {Seg.NumberOfFrames}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs in
    PerFrameFunctionGroupsSequence do not match the SOPInstanceUIDs."""        
    RefSOPuidsInPFFGS = [PFFGS[i].DerivationImageSequence[0]\
                                 .SourceImageSequence[0]\
                                 .ReferencedSOPInstanceUID for i in range(len(PFFGS))]
    
    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInPFFGS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching ReferencedSOPInstanceUID:
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
    if Inds:        
        for i in range(len(Inds)):
            uid = RefSOPuidsInPFFGS[Inds[i]]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + f'PerFrameFunctionalGroupsSequence[{Inds[i]}] does not '\
                  + 'match any DICOM SOPInstanceUID.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs in '\
              + 'PerFrameFunctionalGroupsSequence match the DICOM '\
              + 'SOPInstanceUIDs.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether the DimensionIndexValues are sensible integers,
    and if the indexed ReferencedSOPInstanceUID agrees with the 
    ReferencedSOPInstanceUID within the SourceImageSequence."""
    DIVs = [PFFGS[i].FrameContentSequence[0]\
                    .DimensionIndexValues for i in range(len(PFFGS))]
    
    FirstElements = [int(DIVs[i][0]) for i in range(len(DIVs))]
    
    """ Find the indices of any elements that exceed the number of sequences in
    SegmentSequence, NumSS: """
    Inds = [i for i, x in enumerate(FirstElements) if x > NumSS]
    
    if Inds:        
        for i in range(len(Inds)):
            DIV = DIVs[Inds[i]]
            
            msg = f'ERROR:  The first element in DimensionalIndexValue[{Inds[i]}],'\
                  + f' {DIV}, exceeds the number of sequences in '\
                  + f'SegmentSequene, {NumSS}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The first element in each DimensionalIndexValue are '\
              + 'consistent with the number of sequences in SegmentSequence, '\
              + f'{NumSS}.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
                
    """ Determine if any of the second elements in the DIVs exceed the number
    of sequences in ReferencedInstanceSequence: """
    SecondElements = [int(DIVs[i][1]) for i in range(len(DIVs))]
    
    """ Find the indices of any elements that exceed len(RIS): """
    Inds = [i for i, x in enumerate(SecondElements) if x > len(RIS)]

    if Inds:       
        for i in range(len(Inds)):
            DIV = DIVs[Inds[i]]
            
            msg = f'ERROR:  The second element in DimensionalIndexValue[{Inds[i]}],'\
                  + f' {DIV}, exceeds the number of sequences in '\
                  + f'ReferencedInstanceSequence, {len(RIS)}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The second elements in DimensionalIndexValue '\
                  + 'matches the number of sequences in '\
                  + 'ReferencedInstanceSequence.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    #print(f'\n\n\nlen(RefSOPuidsInPFFGS) = {RefSOPuidsInPFFGS}')
    #print(f'len(SecondElements) = {len(SecondElements)}')
    #print(f'DIVs = {DIVs}')
    
    
    """Determine if the ReferencedSOPInstanceUID in the 
    ReferencedInstanceSequence indexed by the second elements in the DIVs 
    match the ReferencedSOPInstanceUID in the SourceImageSequence."""
    n = 0
    for i in range(len(SecondElements)):
        RefSOPinPFFGS = RefSOPuidsInPFFGS[i]
        
        """ -1 since i is zero-indexed: """
        ind = SecondElements[i] - 1
        
        #print(f'ind = {ind}')
        
        RefSOPinRIS = RefSOPuidsInRIS[ind]
        
        if RefSOPinPFFGS == RefSOPinRIS:
            msg = 'INFO:  ReferencedSOPInstanceUID referenced in '\
                  + f'SourceImageSequence[{i}] matches '\
                  + 'ReferencedSOPInstanceUID referenced in '\
                  + f'ReferencedInstanceSequence[{DIVs[i]}], as indexed by '\
                  + f'DimensionalIndexValue[{i}], {DIVs[i]}.\n'
        else:
            msg = 'ERROR:  ReferencedSOPInstanceUID referenced in '\
                  + f'SourceImageSequence[{i}] does not match '\
                  + 'ReferencedSOPInstanceUID referenced in '\
                  + f'ReferencedInstanceSequence[{DIVs[i]}], as indexed by '\
                  + f'DimensionalIndexValue[{i}], {DIVs[i]}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
            n += 1
                  
    if n == 0:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs referenced in '\
                  + 'SourceImageSequence match the '\
                  + 'ReferencedSOPInstanceUIDs referenced in '\
                  + 'ReferencedInstanceSequence, as indexed by '\
                  + 'DimensionalIndexValue.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine if ImagePositionPatient in PerFrameFunctionalGroupsSequence 
    matches ImagePositionPatient of the DICOM as indexed by the DIVs."""
    IPPsInPFFGS = [PFFGS[i].PlanePositionSequence[0]\
                           .ImagePositionPatient for i in range(len(PFFGS))]
    
    IPPsInPFFGS = [[float(item) for item in IPP] for IPP in IPPsInPFFGS]
    
    c = 0
    for i in range(len(IPPsInPFFGS)):
        if RefSOPuidsInPFFGS[i] in SOPuids:
            ind = SOPuids.index(RefSOPuidsInPFFGS[i])
            
            IPP = IPPs[ind]
            
            #if IPPsInPFFGS[i] != IPP:
            if AreListsEqualToWithinEpsilon(IPPsInPFFGS[i], IPP, epsilon):
                msg = f'INFO: ImagePositionPatient, {IPPsInPFFGS[i]}, '\
                      + f'referenced in PlanePositionSequence[{i}] is within '\
                      + f'{epsilon} of the DICOM ImagePositionPatient, {IPP},'\
                      + f' indexed by the second element, {SecondElements[i]}'\
                      + f', in DimensionalIndexValue[{i}], {DIVs[i]}.\n'
            else:
                msg = f'ERROR: ImagePositionPatient, {IPPsInPFFGS[i]}, '\
                      + f'referenced in PlanePositionSequence[{i+1}] is not '\
                      + f'within {epsilon} of the DICOM ImagePositionPatient,'\
                      + f' {IPP}, indexed by the second element, '\
                      + f'{SecondElements[i]}, in DimensionalIndexValue[{i}],'\
                      + f' {DIVs[i]}.\n'
                LogList.append(msg)
                Nerrors = Nerrors + 1
                c += 1
        
        if c == 0:
            msg = 'INFO: All ImagePositionPatient referenced in '\
                      + f'PlanePositionSequence are within {epsilon} of the '\
                      + 'DICOM ImagePositionPatient indexed by the second '\
                      + 'element in DimensionalIndexValue.\n'
            LogList.append(msg)
            
        #else:
        #    msg = f'.\n'
        #    
        #    Nerrors = Nerrors + 1
                  
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Check that the ReferencedSegmentNumber match the first elements in the
    DimensionIndexValues.
    """
    RSNs = [PFFGS[i].SegmentIdentificationSequence[0]\
                    .ReferencedSegmentNumber for i in range(len(PFFGS))]
    
    Diffs = [FirstElements[i] - RSNs[i] for i in range(len(RSNs))]
    
    Inds = [i for i, x in enumerate(Diffs) if x > 0]
    
    if Inds:       
        for i in range(len(Inds)):
            msg = 'ERROR:  ReferencedSegmentNumber in '\
                  + 'SegmentIdentificationSequence of '\
                  + f'PerFrameFunctionalGroupsSequence[{Inds[i]}], '\
                  + f'{RSNs[Inds[i]]}, does not match the first element in '\
                  + f'DimensionIndexValues[{Inds[i]}], {DIVs[Inds[i]]}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  Each ReferencedSegmentNumber in '\
              + 'SegmentIdentificationSequence of '\
              + 'PerFrameFunctionalGroupsSequence matches the first element '\
              + 'in DimensionIndexValues.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine if the dimensions of the pixel array in PixelData do not
    match the expected dimensions (the number of frames (PixArrF) has 
    already been compared to NumberOfFrames)."""
    if NumOfRows == Seg.Rows:
        msg = f'INFO:  The number of rows in PixelData, {NumOfRows}, '\
              + f'matches the number of SEG Rows, {Seg.Rows}.\n'
        LogList.append(msg)
    else:
        msg = f'ERROR:  The number of rows in PixelData, {NumOfRows}, does '\
              + f'not match the number of SEG Rows, {Seg.Rows}.\n'
        LogList.append(msg)
        Nerrors = Nerrors + 1
    
    
    if LogToConsole:
        print(msg)
            
    if NumOfCols == Seg.Columns:
        msg = f'INFO:  The number of columns in PixelData, {NumOfCols}, '\
              + f'matches the number of SEG Columns, {Seg.Columns}.\n'
        LogList.append(msg)
    else:
        msg = f'ERROR:  The number of columns in PixelData, {NumOfCols}, does'\
              + f' not match the number of SEG Columns, {Seg.Columns}.\n'
        LogList.append(msg)
        Nerrors = Nerrors + 1
    
    if LogToConsole:
        print(msg)
    
    
    #print(f'There were {Nerrors} errors found in the SEG.\n')
    
    if LogToConsole:
        print('-'*120)
            
    return LogList, Nerrors


def GetSegDataFromListOfSegs_OLD(ListOfSegs, ListOfDicomDirs, SearchString, 
                             LogToConsole=False):
    """ 
    Note:
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    #import numpy as np
    from DicomTools import GetDicomFpaths
    #from DicomTools import GetDicomSOPuids
    #from DicomTools import GetRoiNums
    from ImageTools import GetImageAttributes
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of GetSegDataFromListOfSegs():')
        print('\n\n', '-'*120)
    
    ListOfPixArrBySeg = []
    ListOfF2SindsBySeg = []
    #ListOfSegNums = []
    ListOfDcmFpaths = []
    ListOfIPPs = []
    ListOfDirections = []
    ListOfSpacings = []
    
    for i in range(len(ListOfSegs)):
        if LogToConsole:
            print(f'\n\nListOfSegs[{i}]:')
        
        if ListOfSegs[i]:
            PixArrBySeg, F2SindsBySeg = GetPixArrBySeg(ListOfSegs[i], 
                                                       ListOfDicomDirs[i], 
                                                       LogToConsole)
            
            if LogToConsole:      
                print(f'F2SindsBySeg = {F2SindsBySeg}')
            
        else:
            PixArrBySeg = None
            F2SindsBySeg = None
            
            if LogToConsole:      
                print('No segmentations for this dataset')
        
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs,\
        ListOfWarnings = GetImageAttributes(ListOfDicomDirs[i])
        
        ListOfPixArrBySeg.append(PixArrBySeg)
        ListOfF2SindsBySeg.append(F2SindsBySeg)
        ListOfDcmFpaths.append(DcmFpaths)
        #ListOfOrigins.append(IPPs[0])
        ListOfIPPs.append(IPPs)
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
            
        if LogToConsole:      
            print(f'\nlen(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
    
    if LogToConsole:
        print('-'*120)
        
    return ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
           ListOfIPPs, ListOfDirections, ListOfSpacings



def GetSegDataFromListOfSegs(ListOfSegs, ListOfDicomDirs, SearchString, 
                             LogToConsole=False):
    """ 
    Note:
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    #import numpy as np
    #from DicomTools import GetDicomFpaths
    #from DicomTools import GetDicomSOPuids
    #from DicomTools import GetRoiNums
    #from ImageTools import GetImageAttributes
    from ImageTools import GetImAttributesFromListOfDcmDirs
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of GetSegDataFromListOfSegs():')
        print('\n\n', '-'*120)
    
    ListOfPixArrBySeg = []
    ListOfF2SindsBySeg = []
    #ListOfSegNums = []
    
    for i in range(len(ListOfSegs)):
        if LogToConsole:
            print(f'\n\nListOfSegs[{i}]:')
        
        if ListOfSegs[i]:
            PixArrBySeg, F2SindsBySeg = GetPixArrBySeg(ListOfSegs[i], 
                                                       ListOfDicomDirs[i], 
                                                       LogToConsole)
            
            if LogToConsole:      
                print(f'F2SindsBySeg = {F2SindsBySeg}')
            
        else:
            PixArrBySeg = None
            F2SindsBySeg = None
            
            if LogToConsole:      
                print('No segmentations for this dataset')
        
        ListOfPixArrBySeg.append(PixArrBySeg)
        ListOfF2SindsBySeg.append(F2SindsBySeg)
            
        if LogToConsole:      
            print(f'\nlen(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
    ListOfDcmFpaths, ListOfIPPs, ListOfDirections, ListOfSpacings\
        = GetImAttributesFromListOfDcmDirs(ListOfDicomDirs)
    
    if LogToConsole:
        print('-'*120)
        
    return ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
           ListOfIPPs, ListOfDirections, ListOfSpacings



def GetSegDataFromListOfLabImBySeg_OLD(ListOfLabImBySeg, ListOfDicomDirs, #SearchString, 
                                   LogToConsole=False):
    """ 
    Note:  Modelled from GetSegDataFromListOfSegs.
    
    ListOfLabImBySeg is a list (for each dataset, which could be of length 1) 
    of a list (for each segment) of label images (SimpleITK Images).
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    #import numpy as np
    from ConversionTools import ImageBySeg2PixArrBySeg
    from DicomTools import GetDicomFpaths
    #from DicomTools import GetDicomSOPuids
    #from DicomTools import GetRoiNums
    from ImageTools import GetImageAttributes
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of GetSegDataFromListOfSegs():')
        print('\n\n', '-'*120)
    
    ListOfPixArrBySeg = []
    ListOfF2SindsBySeg = []
    #ListOfSegNums = []
    ListOfDcmFpaths = []
    ListOfIPPs = []
    ListOfDirections = []
    ListOfSpacings = []
    
    for i in range(len(ListOfLabImBySeg)):
        if LogToConsole:
            print(f'\n\nListOfLabImBySeg[{i}]:')
        
        if ListOfLabImBySeg[i]:
            PixArrBySeg, F2SindsBySeg = ImageBySeg2PixArrBySeg(ListOfLabImBySeg[i])
            
            if LogToConsole:      
                print(f'F2SindsBySeg = {F2SindsBySeg}')
            
        else:
            PixArrBySeg = None
            F2SindsBySeg = None
            
            if LogToConsole:      
                print('No segmentations for this dataset')
        
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs,\
        ListOfWarnings = GetImageAttributes(ListOfDicomDirs[i])
        
        ListOfPixArrBySeg.append(PixArrBySeg)
        ListOfF2SindsBySeg.append(F2SindsBySeg)
        ListOfDcmFpaths.append(DcmFpaths)
        #ListOfOrigins.append(IPPs[0])
        ListOfIPPs.append(IPPs)
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
            
        if LogToConsole:      
            print(f'\nlen(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
    
    if LogToConsole:
        print('-'*120)
        
    return ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
           ListOfIPPs, ListOfDirections, ListOfSpacings




def GetSegDataFromListOfLabImBySeg(ListOfLabImBySeg, ListOfDicomDirs, #SearchString, 
                                   LogToConsole=False):
    """ 
    Note:  Modelled from GetSegDataFromListOfSegs.
    
    ListOfLabImBySeg is a list (for each dataset, which could be of length 1) 
    of a list (for each segment) of label images (SimpleITK Images).
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    #import numpy as np
    from ConversionTools import ImageBySeg2PixArrBySeg
    #from DicomTools import GetDicomFpaths
    #from DicomTools import GetDicomSOPuids
    #from DicomTools import GetRoiNums
    #from ImageTools import GetImageAttributes
    from ImageTools import GetImAttributesFromListOfDcmDirs
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of GetSegDataFromListOfSegs():')
        print('\n\n', '-'*120)
    
    ListOfPixArrBySeg = []
    ListOfF2SindsBySeg = []
    
    for i in range(len(ListOfLabImBySeg)):
        if LogToConsole:
            print(f'\n\nListOfLabImBySeg[{i}]:')
        
        if ListOfLabImBySeg[i]:
            PixArrBySeg, F2SindsBySeg = ImageBySeg2PixArrBySeg(ListOfLabImBySeg[i])
            
            if LogToConsole:      
                print(f'F2SindsBySeg = {F2SindsBySeg}')
            
        else:
            PixArrBySeg = None
            F2SindsBySeg = None
            
            if LogToConsole:      
                print('No segmentations for this dataset')
        
        
        ListOfPixArrBySeg.append(PixArrBySeg)
        ListOfF2SindsBySeg.append(F2SindsBySeg)
            
        if LogToConsole:      
            print(f'\nlen(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
    
    ListOfDcmFpaths, ListOfIPPs, ListOfDirections, ListOfSpacings\
        = GetImAttributesFromListOfDcmDirs(ListOfDicomDirs)
    
    if LogToConsole:
        print('-'*120)
        
    return ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
           ListOfIPPs, ListOfDirections, ListOfSpacings







"""
******************************************************************************
******************************************************************************
COPY A SEGMENTATION / SEGMENT / SEG
******************************************************************************
******************************************************************************
"""

def CopySeg(SrcSegFpath, SrcSliceNum, SrcSegLabel, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgSegFpath=None, TrgSegLabel=None, TrgSliceNum=None, 
            ResInterp='BlurThenLinear', PreResVariance=(1,1,1),
            ApplyPostResBlur=False, PostResVariance=(1,1,1),
            ForceReg=False, UseDroForTx=True, SelxOrSitk='Sitk', 
            Transform='affine', MaxIters='512', InitMethod='centerofgravity', 
            SrcFidsFpath='moving_fiducials.txt', 
            TrgFidsFpath='fixed_fiducials.txt', TxInterp='NearestNeighbor', 
            ApplyPostTxBin=True, ApplyPostTxBlur=True, PostTxVariance=(1,1,1), 
            #TxMatrix=None, GridDims=None, GridRes=None, VectGridData=None,
            Dro=None, TxtToAddToSegLabel='', #OverwriteDro=False, 
            LogToConsole=False, 
            DictOfInputs={}, ListOfInputs=[], ListOfTimings=[], 
            SampleDroDir='default'):
    """
    Note 09/02/2021:
        This function will, depending on the inputs, copy either:
            1. A specific segmentation
            2. All segmentations in a specific segment
            3. All segmentations in all segments
        
     
    Inputs:
    ******
    
    SrcSegFpath : string
        Filepath of the Source SEG file.
    
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a  
        single segmentation). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
        
    SrcSegLabel : string
        All or part of the Source SegmentLabel for the segment containing the
        segmentation(s) to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
        
    TrgRtsFpath : string or None (optional; None by default)
        Filepath to the Target SEG file that the segmentation(s) is/are to be 
        copied to. TrgSegFpath != None only if an existing Target SEG exists 
        and is to added to. An existing SEG will be added to only for Direct 
        copy operations.     
    
    TrgSegLabel : string or None (optional; None by default) 
        All or part of the SegmentLabel of the destination segment.
        
    TrgSliceNum : integer (optional; None by default)
        The slice index within the Target DICOM stack where the segmentation 
        will be copied to (applies only for the case of direct copies of single 
        segmentations). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the segmentation(s) will be copied to 
        will not depend on user input.
        
    ResInterp : string
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    PreResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to perform binary thresholding following 
    #    resampling if a non-label (e.g. linear) interpolator is used. Example 
    #    use is on a labelmap image if a Gaussian blurring + linearly resampling 
    #    approach is taken to combat aliasing effects.
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
        
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
    
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, even 
        if the Source and Target images have the same FrameOfReferenceUID.  
    
    UseDroForTx : boolean (optional; True by default)
        If True the transformation parameters from any suitable DRO (for the
        required Transform) will be used instead of performing image
        registration.
    
    SelxOrSitk : string (optional; 'Selx' by default)
        Denotes which package to use for image registration and transformation.
        Acceptable values include:
            - 'Selx' for SimpleElastix
            - 'Sitk' for SimpleITK
            
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxIters : string (optional; '512' by default)
        The maximum number of iterations used for the optimiser during image 
        registration (if applicable).
    
    InitMethod : string (optional, 'centerofgravity' by default)
        The initialisation method used for initial alignment prior to image
        registration (if applicable).
    
    SrcFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Source/Moving image.
    
    TrgFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Target/Fixed image.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
        
    #ThreshPostTx : float (optional; 0.05 by default)
    #    The threshold level used to perform binary thresholding after 
    #    registration transformation.
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
        
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
    
    TxtToAddToSegLabel : string (optional, '' by default)
        String of text to add to the segment label of the new Target SEG.
    
    #TxMatrix : list of float strings or None (optional; None by default)
    #    List of float strings representing the non-deformable transformation 
    #    that transforms the moving (Source) image to the fixed (Target) image;
    #    or None.  If not None, image registration will be skipped to save 
    #    computational time.
    #
    #GridDims : list of integers or None (optional; None by default)
    #    List of integers representing the dimensions of the BSpline grid used 
    #    to deform the moving (Source) image to the fixed (Target) image; 
    #    or None.  If not None, image registration will be skipped to save 
    #    computational time.
    # 
    #GridRes : list of floats or None (optional; None by default)
    #    List of floats representing the resolution of the BSpline grid used to 
    #    deform the moving (Source) image to the fixed (Target) image; or None.
    #    If not None, image registration will be skipped to save computational
    #    time.
    #
    #VectGridData : list of float strings or None (optional; None by default)
    #    List of floats representing the vector deformations that deform the
    #    moving (Source) image to the fixed (Target) image; or None.
    #    If not None, image registration will be skipped to save computational
    #    time.
    
    Dro : Pydicom Object (optional; None by default)
        A Spatial or Deformable Spatial Registration Object that contains the
        parameters required to register the moving (Source) image to the fixed
        (Target) image.
    
    OverwriteDro : boolean (optional; False by default)
        If True, Dro (if not None) will be overwritten.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    DictOfInputs : dictionary (optional; empty by default)
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings (optional; empty by default)
        A list containing the inputs that were called to CopyRoi().
    
    ListOfTimings : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopySeg().
          
                  
    Outputs:
    *******
        
    TrgSeg : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target SEG object.
    
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgSeg, None otherwise.
    
    DictOfInputs : dictionary
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
        
    ListOfTimings : list of strings
        A list (for each message) of timings for individual operations.
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target segmentations will appear displaced  
    w.r.t. the anatomical features highlighted in the Source image.
    """
    
    import importlib
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    import DroTools
    importlib.reload(DroTools)
    import RoiCopyTools
    importlib.reload(RoiCopyTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from ImageTools import ImportImage#, GetImageInfo
    #from SegTools import GetSegDataOfInterest, ProportionOfSegsInExtent
    #from SegTools import CreateSeg
    from GeneralTools import UniqueItems#, PrintIndsByRoi#, PrintTitle
    #from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    #from GeneralTools import ShiftFrame
    from DroTools import CreateSpaDro, CreateDefDro
    from RoiCopyTools import CheckValidityOfInputs
    
    """ Start timing. """
    times = []
    times.append(time.time())
    #""" Log timing messages. """
    #TimingMsgs = []
    
    
    SrcSeg = dcmread(SrcSegFpath)
    if TrgSegFpath:
        TrgSeg = dcmread(TrgSegFpath)
    else:
        TrgSeg = None
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print('Running of CopySeg():')
        print('\n\n', '-'*120)
        
        SrcSegLabels = GetRoiLabels(SrcSeg)
        
        print(f'\nSrcSegLabels = {SrcSegLabels}')
        
        if TrgSegFpath:
            TrgSegLabels = GetRoiLabels(TrgSeg)
        
            print(f'\nTrgSegLabels = {TrgSegLabels}')
    
    
    """ Check the validity of the combination of inputs. """
    CheckValidityOfInputs(SrcSeg, SrcSegLabel, SrcSliceNum, TrgSeg, TrgSegLabel,
                          TrgSliceNum, LogToConsole)
    
    #""" Determine which Use Case to apply. """
    #UseCaseThatApplies,\
    #UseCaseToApply = WhichUseCase(SrcSliceNum, FromSegLabel, TrgSliceNum, 
    #                              SrcDcmDir, TrgDcmDir, ForceRegistration,
    #                              LogToConsole=True)
    
    #""" Manipulate UseCase if ForceRegistration = True and UseCase was 3 or 4. """
    #if ForceRegistration and not '5' in UseCase:
    #    msg = f'\nUseCase was {UseCase} but will be treated as UseCase '
    #    
    #    UseCase.replace('3', '5').replace('4', '5')
    #    
    #    msg = f'{UseCase}.'
    #    
    #    print(msg)
        
    """ Get the data of interest from the Source SEG. """
    SrcPixArrBySeg,\
    SrcF2SindsBySeg = GetSegDataOfInterest(SrcSegFpath, SrcSliceNum, 
                                           SrcSegLabel, SrcDcmDir, 
                                           LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source SEG file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Load the Source SEG: """
    SrcSeg = dcmread(SrcSegFpath)
    
    
    """ Raise exception if there is no data of interest from the Source SEG."""
    if not UniqueItems(SrcF2SindsBySeg):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        #from SegTools import GetPFFGStoSliceIndsBySeg
        
        SrcSeg = dcmread(SrcSegFpath)
        
        Labels = GetRoiLabels(SrcSeg)
        L = len(Labels)
        
        SOPuids = GetDicomSOPuids(SrcDcmDir)
        
        F2SindsBySeg = GetPFFGStoSliceIndsBySeg(SrcSeg, SOPuids)
        
        msg = f"There are no segmentations on slice {SrcSliceNum} in any "\
              + f"segment in the Source SEG. There are {L} segments in the "\
              + "Source SEG:"
        
        for i in range(L):
            msg += f"\nSegment {i+1} with label '{Labels[i]}' has "\
                   + f"{len(F2SindsBySeg[i])} segmentations that correspond "\
                   + f"to slices {F2SindsBySeg[i]}" 
                   
        raise Exception(msg)
    
    """
    Determine whether the segmentations that make up the segment(s) intersect 
    with the Target image extent. 
    """
    FracProp = ProportionOfSegsInExtent(PixArrBySeg=SrcPixArrBySeg,
                                        F2SindsBySeg=SrcF2SindsBySeg, 
                                        SrcDicomDir=SrcDcmDir,
                                        TrgDicomDir=TrgDcmDir,
                                        LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to check whether the Source segment(s) intersect '\
          + 'the Target grid space.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the indices in the SEG lie',
              'within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the indices in the SEG lie within the physical extent'\
              + 'of the Target image.'
        raise Exception(msg)
        return None
    
    """ Import the 3D images. """
    SrcIm = ImportImage(SrcDcmDir)
    TrgIm = ImportImage(TrgDcmDir)
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    #if UseCase in ['1', '2a', '3a']: # 08/02
    if UseCaseToApply in ['1', '2a']: # 08/02
        """ Direct copy of a segmentation without resampling or registration
        transformation. """
        
        from GeneralTools import ReplaceIndInC2SindsByRoi
        from GeneralTools import ShiftFramesInPixArrBySeg
        
        """ The segmentation on SrcSliceNum is to be copied to TrgSliceNum. 
        Modify the frame-to-slice index (from SrcSliceNum) to TrgSliceNum. """
        SrcF2SindsBySeg = ReplaceIndInC2SindsByRoi(C2SindsByRoi=SrcF2SindsBySeg,
                                                   IndToReplace=SrcSliceNum, 
                                                   ReplacementInd=TrgSliceNum)
        
        if LogToConsole:
                print(f'\nModfied SrcF2SindsBySeg = {SrcF2SindsBySeg}.')
                
        if TrgSegFpath:
            """ Direct copy of a segmentation to an existing SEG. """
            
            """ An existing Target SEG is to be modified.  Since this is a  
            Direct copy, any existing segmentations in Target with segment 
            label matching FromSegLabel are to be preserved. """
            
            """ Get the data of interest from the Target RTS. """
            TrgPixArrBySeg,\
            TrgF2SindsBySeg = GetSegDataOfInterest(TrgSegFpath, SrcSliceNum,
                                                   #SrcSegLabel, TrgDcmDir,
                                                   TrgSegLabel, TrgDcmDir,
                                                   LogToConsole)
            """ 05/03/21: Previously SrcSegLabel was used as an input in 
            GetSegDataOfInterest for Target, hence requiring that the Target
            segment have the same label as the Source segment. Now TrgSegLabel 
            has been added to the list of inputs allowing for the Target segment
            to have a different label. """
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to parse the Target SEG file for '\
                  + 'existing segmentations within the ROI.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
            
            if LogToConsole:
                print(f'TrgF2SindsByRoi = {TrgF2SindsBySeg}.')
        else:
            TrgPixArrBySeg = None
            TrgF2SindsBySeg = None
        

    #if UseCase in ['1', '2a']:
        #""" Direct copy without resampling or registration transformation.
        """
        The segmentation will be copied to a different slice location, so the
        z-components of the indices will need to be modified. """
        
        """ Shift the in-plane (x & y) and out-of-plane (z) components of the 
        indices in SrcPixArrBySeg to account for the shift from SrcSliceNum to
        TrgSliceNum. """
        SrcPixArrBySeg, SrcF2SindsBySeg\
        = ShiftFramesInPixArrBySeg(PixArrBySeg=SrcPixArrBySeg, 
                                   F2SindsBySeg=SrcF2SindsBySeg, 
                                   SrcImage=SrcIm, 
                                   SrcSliceNum=SrcSliceNum,
                                   TrgImage=TrgIm, 
                                   TrgSliceNum=TrgSliceNum, 
                                   RefImage=TrgIm,
                                   ShiftInX=True,
                                   ShiftInY=True,
                                   ShiftInZ=True,
                                   Fractional=False,
                                   LogToConsole=LogToConsole)
        
        if LogToConsole:
                print(f'\nShifted SrcF2SindsBySeg = {SrcF2SindsBySeg}.')
                
        
        if TrgSegFpath:
            """ An existing Target SEG is to be modified. 
            Concatenate TrgF2SindsBySeg and SrcF2SindsBySeg, and 
            TrgPixArrBySeg and SrcPixArrBySeg. """
            
            TrgF2SindsBySeg = [TrgF2SindsBySeg[s] + SrcF2SindsBySeg[s] for s in range(len(TrgF2SindsBySeg))]
            
            TrgPixArrBySeg = [TrgPixArrBySeg[s] + SrcPixArrBySeg[s] for s in range(len(TrgPixArrBySeg))]
        
        else:
            """ A Target SEG was not provided. """
            
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            
            if LogToConsole:
                print(f'\nSrcF2SindsBySeg = {SrcF2SindsBySeg}')
                print(f'TrgF2SindsBySeg = {TrgF2SindsBySeg}')
    
    
    
    if UseCaseToApply == '2b':
        """ Relationship-preserving copy without resampling or registration
        transformation. """
        
        if False:
            """ TrgPixArrBySeg is simply SrcPixArrBySeg, and 
            TrgF2SindsBySeg is SrcF2SindsBySeg. """
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
        
        #from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1
        from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2
        
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        TrgF2SindsBySeg = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcF2SindsBySeg, 
                                                                SrcIm, TrgIm)
        
        
        """ TrgPixArrBySeg is simply SrcPixArrBySeg. """
        TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
        
        if LogToConsole:
            print(f'\nSrcF2SindsBySeg = {SrcF2SindsBySeg}')
            print('TrgF2SindsBySeg (= shifted SrcF2SindsBySeg) =',
                  f'{SrcF2SindsBySeg}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        #from ImageTools import ResampleLabImByRoi#, GetImageInfo
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PixArrByRoi2LabImByRoi
        #from GeneralTools import NumOfListsAtDepthTwo
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabImBySeg = PixArrByRoi2LabImByRoi(PixArrByRoi=SrcPixArrBySeg, 
                                               F2SindsByRoi=SrcF2SindsBySeg, 
                                               RefIm=SrcIm,
                                               LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to generate labelmap images from the Source '\
              'pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'\n*{msg}')
    
    
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        from ImageTools import ResampleLabImByRoi
        
        
        """ Resample the source labelmaps. """
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
        = ResampleLabImByRoi(LabImByRoi=SrcLabImBySeg, 
                             #F2SindsByRoi=SrcC2SindsByRoi,
                             F2SindsByRoi=SrcF2SindsBySeg,
                             SrcIm=SrcIm, 
                             TrgIm=TrgIm,
                             Interp=ResInterp,
                             PreResVariance=PreResVariance,
                             #PostResThresh=PostResThresh,
                             ApplyPostResBlur=ApplyPostResBlur,
                             PostResVariance=PostResVariance,
                             LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the interpolation that was used (the actual interpolation used 
        might have been different from the interpolation set) and the threshold
        used to re-binarise the resampled labelmap image (if applicable): """
        
        #print(f'\n\n\nMetadata keys in ResSrcLabImBySeg[0]:',
        #      ResSrcLabImBySeg[0].GetMetaDataKeys())
        
        for i in range(len(ResSrcLabImBySeg)):
            Interp = ResSrcLabImBySeg[i].GetMetaData("ResInterpUsed")
            
            ListOfInputs.append(f'ResInterpUsedForSeg{i} = {Interp}')
            DictOfInputs[f'ResInterpUsedForSeg{i}'] = Interp
            
            if Interp == 'BlurThenLinear':
                Thresh = ResSrcLabImBySeg[i].GetMetaData("PostResThreshUsed")
                
                ListOfInputs.append(f'PostResThreshUsedForSeg{i} = {Thresh}')
                DictOfInputs[f'PostResThreshUsedForSeg{i}'] = Thresh
    else:
        ResSrcLabImBySeg = None
        ResSrcPixArrBySeg = None
        ResSrcF2SindsBySeg = None
                
        
    
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        from importlib import reload
        import DroTools
        reload(DroTools)
        
        if SelxOrSitk == 'Selx':
            from SelxTools import RegisterImagesSelx, GetParamFromSelxImFilt
        else:
            from ImageTools import RegisterImagesSitk
        from ImageTools import TransformImWithListOfSitkTxs
        #from ImageTools import TransformLabImByRoi
        from ImageTools import TransformLabImByRoiWithListOfSitkTxs
        #from DroTools import CreateSitkTxFromTxMatrix
        from DroTools import CreateSitkTxFromSpaDro
        from DroTools import CreatePreDefSitkTxFromDefDro
        from DroTools import CreateSitkTxFromDefDro
        from GeneralTools import AreListsEqualToWithinEpsilon
        
        if not Transform in ['rigid', 'affine', 'bspline']:
            msg = f'The chosen transform, {Transform}, is not one of the '\
                  + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
            
            raise Exception(msg)
        
        #if (Transform in ['rigid', 'affine'] and TxMatrix == None) \
        #or (Transform == 'bspline' and VectGridData == None):
        #if Dro == None:
        #if Dro == None or (Dro != None and UseDroForTx): # 01/09/21
        if Dro == None or (Dro != None and not UseDroForTx): # 01/09/21
            """ The transform parameters required to register SrcIm to TrgIm 
            were not found from any matching subject assessors in XNAT, or
            a suitable DRO was found but UseDroForTx = False. 
            Image registration to be performed. """
            
            if Dro != None and UseDroForTx:
                print('Image registration will be performed despite a',
                      'suitable DRO found on XNAT.\n')
            
            if LogToConsole == False:
                print(f'Performing image registration using {Transform}',
                      'transform...\n')
            
            times.append(time.time())
            
            """ Register SrcIm to TrgIm. """
            if SelxOrSitk == 'Selx':
                RegIm, SelxImFiltOrSitkTx, TxParamMapDict\
                = RegisterImagesSelx(FixIm=TrgIm, MovIm=SrcIm, 
                                     Transform=Transform,
                                     MaxNumOfIters=MaxIters,
                                     InitMethod=InitMethod,
                                     FixFidsFpath=TrgFidsFpath, 
                                     MovFidsFpath=SrcFidsFpath,
                                     FlipK=False, LogToConsole=LogToConsole)
                
                #""" Get the registration transform parameters: """
                #TxParams = list(TxParamMapDict['TransformParameters'])
                
                """ Get the registration transform parameters for the final
                transform parameter map: (10/05/21) """
                TxParams = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                                  Param='TransformParameters')
                
            if SelxOrSitk == 'Sitk':
                RegIm, LandmarkTx, SelxImFiltOrSitkTx\
                = RegisterImagesSitk(FixIm=TrgIm, MovIm=SrcIm, 
                                     Transform=Transform, InitMethod=InitMethod,
                                     FixFidsFpath=TrgFidsFpath, 
                                     MovFidsFpath=SrcFidsFpath,
                                     LogToConsole=LogToConsole)
                
                """ Get the registration transform parameters: """
                TxParams = [str(item) for item in SelxImFiltOrSitkTx.GetParameters()]
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to register the Source image to the Target '\
                  + f'image using {Transform} transform.\n'
            ListOfTimings.append(msg)
            if LogToConsole == False:
                print(f'*{msg}')
            
            #""" Add the row vector [0, 0, 0, 1] to TxParams to complete the 
            #Frame of Reference Transformation Matrix 
            #(http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
            #"""
            ##TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
            #TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
            #            str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
            #            str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
            #            str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
            
            
            """ 
            04/06/2021:
            
            Due to issues with initialising a SimpleITK BSpline Transform
            using fiducials, the workflow for deformable registrations using
            fiducial pre-alignment consists of two resampling steps (as opposed 
            to one for deformable registrations not using fiducial pre-alignment
            or non-deformable registrations using fiducials or not). Hence
            create a list ListOfSitkTxs containing SelxImFiltOrSitkTx so
            that TransformLabImByRoi() can be looped for each item in the list.
            """
            ListOfSitkTxs = [SelxImFiltOrSitkTx]
            
        else:
            """ The transform parameters required to register SrcIm to TrgIm 
            are obtained from the transform matrix, found from the DRO (subject
            assessor) in XNAT. Create a SimpleITK transform from Dro: 
            """
            
            """ 04/06/21 Workaround described above. """
            ListOfSitkTxs = []
            
            if Transform in ['rigid', 'affine']:
                #SelxImFiltOrSitkTx = CreateSitkTxFromTxMatrix(TrgIm, TxMatrix, 
                #                                              Transform)
                #
                #""" TxMatrix has additional bottom row [0 0 0 1] (see 
                #http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
                #
                #Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. 
                #Hence:
                #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
                #"""
                #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 06/05/21
                
                SitkTx = CreateSitkTxFromSpaDro(Dro, LogToConsole)
                
                TxParams = list(SitkTx.GetParameters())
                
                ListOfSitkTxs.append(SitkTx)
            else:
                PreRegTx = CreatePreDefSitkTxFromDefDro(Dro, LogToConsole)
                
                TxParams = list(PreRegTx.GetParameters())
                
                Identity = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
                
                if not AreListsEqualToWithinEpsilon(TxParams, Identity):
                    """ The pre-registration matrix is not the identity matrix
                    so need to resample/transform the label image using 
                    PreRegTx. """
                
                    ListOfSitkTxs.append(PreRegTx)
                
                SitkTx = CreateSitkTxFromDefDro(Dro, RefIm=TrgIm,
                                                LogToConsole=LogToConsole)
                
                ListOfSitkTxs.append(SitkTx)
                
            SelxImFiltOrSitkTx = deepcopy(SitkTx) # 01/09/21 added to allow 
            # CreateSpaDroFromSitkTx() (line 3942) to run
            
            """ Since image registration was skipped, create RegIm by applying
            ListOfSitkTxs iteratively: """
            RegIm = TransformImWithListOfSitkTxs(MovIm=SrcIm, FixIm=TrgIm, 
                                                 ListOfSitkTxs=ListOfSitkTxs, 
                                                 Interp='Linear', 
                                                 LogToConsole=LogToConsole)
        
        
        """ Transform SrcLabImBySeg using the registration transformation. 
        
        Notes:
            
        Even though the data is being transformed the prefix 'Res' will be
        used since further operations below will need to be applied to
        either resampled or transformed data. 
            
        04/06/21: Workaround to issue with initialising a SimpleITK BSpline
        transform using fiducials involves a list (of potentially length 1)
        for the variable ListOfSitkTxs.  This is not ideal not only due to
        stacking of resampling errors but also due to use of a nearest
        neighbor interpolator.
        """
        
        #ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
        #= TransformLabImByRoi(LabImByRoi=SrcLabImBySeg,
        #                      F2SindsByRoi=SrcF2SindsBySeg,
        #                      TrgIm=TrgIm,
        #                      SelxImFiltOrSitkTx=SelxImFiltOrSitkTx,
        #                      Interp=TxInterp,
        #                      ApplyPostTxBlur=ApplyPostTxBlur,
        #                      PostTxVariance=PostTxVariance,
        #                      ApplyPostTxBin=ApplyPostTxBin,
        #                      #ThreshPostTx=ThreshPostTx, 
        #                      LogToConsole=LogToConsole)
        
        ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
        = TransformLabImByRoiWithListOfSitkTxs(LabImByRoi=SrcLabImBySeg,
                                               F2SindsByRoi=SrcF2SindsBySeg,
                                               TrgIm=TrgIm,
                                               ListOfSitkTxs=ListOfSitkTxs,
                                               Interp=TxInterp,
                                               ApplyPostTxBlur=ApplyPostTxBlur,
                                               PostTxVariance=PostTxVariance,
                                               ApplyPostTxBin=ApplyPostTxBin,
                                               LogToConsole=LogToConsole)
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
        
        if ApplyPostTxBlur:
            """ Get the threshold used to re-binarise the transformed label 
            image: """
            for i in range(len(ResSrcLabImBySeg)):
                Thresh = ResSrcLabImBySeg[i].GetMetaData("PostTxThreshUsed")
                
                ListOfInputs.append(f'PostTxThreshUsedForSeg{i} = {Thresh}')
                DictOfInputs[f'PostTxThreshUsedForSeg{i}'] = Thresh
    else:
        RegIm = None
        SelxImFiltOrSitkTx = None
        TxParams = None
        ListOfSitkTxs = None
              
                
    
    if UseCaseToApply in ['3a', '4a', '5a']:
        """ Direct copy with resampling/registration and averaging of multi-
        framed pixel arrays.
        
        Note:
        
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled/registered labelmap image).  So if after resampling there are 
        more than one frame, the frames will be averaged or a logical OR 
        operation will be performed.
        """
        
        from GeneralTools import ShiftFramesInPixArrBySeg
        from PlottingTools import PlotPixArrBySeg
        
        R = len(ResSrcPixArrBySeg) # number of segments
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
        MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
        
        if MaxNumOfFramesBySeg > 1:
            """ Reduce the number of frames in each pixel array to 1. """
            
            """ Perform averaging or OR operation? """
            #ReduceUsing = 'mean'
            ReduceUsing = 'OR'
            
            if LogToConsole:
                print(f'\nThere are {NumOfFramesBySeg} frames in each',
                      f'ROI/segment. Using {ReduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                PlotPixArrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                f'Prior to {ReduceUsing} operation')
                
            if ReduceUsing == 'mean':
                from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5 # <-- is this too high?
    
                ResSrcPixArrBySeg, ResF2SindsBySeg\
                = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                         F2SindsBySeg=ResSrcF2SindsBySeg,
                                         MakeBinary=True, 
                                         BinaryThresh=BinaryThresh,
                                         LogToConsole=LogToConsole)
            
            else:
                from GeneralTools import OrFrameOfPixArrBySeg
                
                ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
                = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                       F2SindsBySeg=ResSrcF2SindsBySeg,
                                       LogToConsole=LogToConsole)
        
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment, and',
                      f'the new ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}.')
                
                PlotPixArrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                'After {ReduceUsing} operation')
                
                
        
        """ Shift the out-of-plane elements in ResSrcPixArrBySeg to account for
        the shift from SrcSliceNumber to TrgSliceNum. 
        Note:
            The function below was incorrectly applying in-plane shifts as well
            up until the added inputs ShiftInX / Y / Z was added on 10/02. 
        
        14/02: I think SrcSliceNum shouldn't be SrcSliceNum but instead
        ResSrcF2SindsBySeg[0][0], since this is the frame-to-slice number of the
        resampled/transformed labelmap.  The slice number of the original
        (Source) labelmap is not relevant.
        
        Likewise, the inputs SrcImage and TrgImage should both be TrgImage 
        since the resampled labelmap shares the TrgImage grid, and it will be
        shifted to another location, still within the TrgImage grid.
        """
        
        """ SrcSliceNum --> ResSrcSliceNum following resampling or 
        transformation. """
        ResSrcSliceNum = ResSrcF2SindsBySeg[0][0]
        
        if LogToConsole:
            #print(f'\nSrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
            #      f'{ResSrcSliceNum}')
            #print(f'TrgSliceNum = {TrgSliceNum}')
            print(f'SrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
                  f'{ResSrcSliceNum} following resampling/transformation -->',
                  f'TrgSliceNum = {TrgSliceNum} for direct copy')
            print('\nResSrcF2SindsBySeg prior to shifting of out-of-plane',
                  f'elements: {ResSrcF2SindsBySeg}\n')
            
        ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
        = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg, 
                                   F2SindsBySeg=ResSrcF2SindsBySeg, 
                                   SrcImage=TrgIm, 
                                   SrcSliceNum=ResSrcSliceNum,
                                   TrgImage=TrgIm, 
                                   TrgSliceNum=TrgSliceNum, 
                                   RefImage=TrgIm,
                                   ShiftInX=False,
                                   ShiftInY=False,
                                   ShiftInZ=True,
                                   Fractional=False,
                                   LogToConsole=LogToConsole)
        
        if LogToConsole:
            print('After running ShiftFramesInPixArrBySeg():')
            print(f'ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}')
            print(f'len(ResSrcPixArrBySeg) = {len(ResSrcPixArrBySeg)}')
            for r in range(len(ResSrcPixArrBySeg)):
                print(f'type(ResSrcPixArrBySeg[{r}]) = {type(ResSrcPixArrBySeg[r])}')
                print(f'ResSrcPixArrBySeg[{r}].shape = {ResSrcPixArrBySeg[r].shape}')
            print('')
                
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copying with resampling or
        registration-transformation. """
    
        if UseCaseToApply in ['3a', '4a', '5a']:
            """ Direct copying with resampling/transformation. """
            
            if TrgSegFpath:
                """ An existing Target SEG is to be modified. 
                Concatenate TrgF2SindsBySeg and ResSrcF2SindsBySeg, and
                TrgPixArrBySeg and ResSrcPixArrBySeg. """
                
                TrgF2SindsBySeg = [TrgF2SindsBySeg[s] + ResSrcF2SindsBySeg[s] for s in range(len(TrgF2SindsBySeg))]
                
                TrgPixArrBySeg = [TrgPixArrBySeg[s] + ResSrcPixArrBySeg[s] for s in range(len(TrgPixArrBySeg))]
                
            else:
                """ A Target SEG was not provided. """
                
                TrgF2SindsBySeg = deepcopy(ResSrcF2SindsBySeg)
                
                TrgPixArrBySeg = deepcopy(ResSrcPixArrBySeg)
        
        
        else:
            """ UseCaseToApply is '3b', '4b' or '5b'. """
            
            TrgPixArrBySeg = deepcopy(ResSrcPixArrBySeg)
            
            TrgF2SindsBySeg = deepcopy(ResSrcF2SindsBySeg)
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    #print(f'\n\n\n\n\nTrgPixArrBySeg = {TrgPixArrBySeg}')
    
    
    """ Create/overwrite the Target SEG. """
    
    print('Creating the new Target SEG object...\n')
    
    #print(f'\nTrgDcmDir = {TrgDcmDir}\n')
    
    TrgSeg = CreateSeg(SrcSegFpath, TrgSegFpath, TrgSegLabel, TrgPixArrBySeg, 
                       TrgF2SindsBySeg, TrgDcmDir, SrcSegLabel,
                       TxtToAddToSegLabel, LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to create the SEG object.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Create a DICOM Registration Object (if applicable). """
    #if UseCaseToApply in ['5a', '5b'] and (Dro != None and OverwriteDro
    #                                       or Dro == None):
    if UseCaseToApply in ['5a', '5b']:
        from DroTools import CreateSpaDroFromSitkTx, CreateDefDroFromBsplineTx
        
        # ContentDescription:
        """ Consider adding this as an input variable for the user to specify.
        """
        ContentDesc = ''
        
        if Transform in ['rigid', 'affine']:
            """ Create a Spatial Registration DRO. """
            #Dro = CreateSpaDro(SrcDcmDir, TrgDcmDir, TxParams, ContentDesc, 
            #                   LogToConsole)
            Dro = CreateSpaDroFromSitkTx(SampleDroDir,
                                         SrcDcmDir, TrgDcmDir, 
                                         SelxImFiltOrSitkTx,
                                         ContentDesc, LogToConsole)
        else:
            #""" Get the bspline grid related registration transform parameters 
            #for the final transform parameter map: (10/05/21) """
            #GridOrig = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1, 
            #                                  'GridOrigin')
            #
            #GridDir = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                 'GridDirection')
            #
            #GridDims = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                  'GridSize')
            #GridDims = [int(item) for item in GridDims]
            #
            #GridRes = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                 'GridSpacing')
            #GridRes = [float(item) for item in GridRes]
            #
            #VectGridData = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
            #                                      'TransformParameters')
            #VectGridData = [float(item) for item in VectGridData]
            #
            #""" Get the non-deformable registration transform parameters 
            #for the second-last transform parameter map: (12/05/21) """
            #TxParams = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -2,
            #                                  'TransformParameters')
            #if TxParams != None:
            #    TxParams = [float(item) for item in TxParams]
            
            """ Create a Deformable Spatial Registration DRO. """
            #Dro = CreateDefDro(SrcDcmDir, TrgDcmDir, GridOrig, GridDir, 
            #                   GridDims, GridRes, VectGridData, TxParams,
            #                   ContentDesc, LogToConsole)
            #Dro = CreateDefDroFromBsplineTx(SrcDcmDir, TrgDcmDir, 
            #                                BsplineTx=SelxImFiltOrSitkTx, 
            #                                PreRegTx=LandmarkTx, 
            #                                Description=ContentDesc, 
            #                                LogToConsole=LogToConsole)
            """ Note:
                The first Transform in ListOfSitkTxs is the landmark-based
                Transform, the second/last one is the final BSpline Transform.
            """
            Dro = CreateDefDroFromBsplineTx(SampleDroDir,
                                            SrcDcmDir, TrgDcmDir, 
                                            BsplineTx=ListOfSitkTxs[-1], 
                                            PreRegTx=ListOfSitkTxs[0],
                                            Description=ContentDesc, 
                                            LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to create the DRO.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
    else:
        Dro = None
        
        
    if LogToConsole:
        print('\nEnd of CopySeg().')
        print('-'*120)
        
    #return TrgSeg, Dro, DictOfInputs, ListOfInputs, ListOfTimings
    ##return TrgSeg, TrgPixArrBySeg, TrgF2SindsBySeg, ListOfTimings
    return SrcDcmDir, TrgDcmDir, SrcIm, TrgIm, SrcSeg,\
           SrcPixArrBySeg, SrcF2SindsBySeg, SrcLabImBySeg,\
           TrgPixArrBySeg, TrgF2SindsBySeg,\
           RegIm, ListOfSitkTxs, SelxImFiltOrSitkTx, TxParams,\
           ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg,\
           TrgSeg, Dro, DictOfInputs, ListOfInputs, ListOfTimings





