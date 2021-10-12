# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 09:31:12 2021

@author: ctorti
"""


"""
Old code removed from various modules that were used in the OOP (object 
oriented programming) related code.
"""



""" From conversion_tools.inds_pts_cntdata.py """

def ind_to_pt_210921(index, refIm):
    """
    14/09/21: This is unnecessarily complex.
    
    Convert an index (in the Image Coordinate System) to a physical point (in 
    the Patient Coordinate System).
    
    Parameters
    ----------
    index : list of ints or floats
        A list (for each dimension) of an index, e.g. 
        [[i0, j0, k0], [i1, j1, k1], ...].
    refIm : SimpleITK image
        A 3D image that occupies the grid that index belongs to.
        
    Returns
    -------
    point : list of floats
        A list (for each dimension) of the coordinates of index, e.g.
        [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    #import numpy
    #from general_tools.general import get_list_of_dtypes
    
    dtypes = get_list_of_dtypes(index)
    
    if len(dtypes) == 1:
        dtype = dtypes[0]
        
        """All items in index are of a common data type."""
        if dtype == int:
            point = refIm.TransformIndexToPhysicalPoint(index)
        
        elif dtype == float:
            point = refIm.TransformContinuousIndexToPhysicalPoint(index)
        
        elif dtype in [numpy.ndarray, numpy.float64, numpy.int32]:
            """ Convert from numpy array to list: """
            index = index.tolist()
            
            dtypes = get_list_of_dtypes(index)
            
            if len(dtypes) == 1:
                dtype = dtypes[0]
                
                """All items in index are of a common data type."""
                if dtype == int:
                    point = refIm.TransformIndexToPhysicalPoint(index)
                
                elif dtype == float:
                    point = refIm.TransformContinuousIndexToPhysicalPoint(index)
                
                else:
                    msg = f"The data type of index is {dtype}. It must be "\
                          + "'int', 'float', 'numpy.ndarray' or 'numpy.float64'."
            
                    raise Exception(msg)
                    
            elif len(dtypes) == 2 and int in dtypes and float in dtypes:
                """Allow for the possibility of a mixture of integers and 
                floats."""
                
                """ Convert integers to float: """
                index = [float(item) for item in index]
                
                point = refIm.TransformContinuousIndexToPhysicalPoint(index)
                
            else:
                msg = f"The data type of index is {dtypes}. It must be 'int',"\
                      + " 'float' or 'numpy.ndarray'."
            
        else:
            msg = f"The data type of index is {dtypes}. It must be 'int',"\
                  + " 'float', 'numpy.ndarray' or 'numpy.float64'."
            raise Exception(msg)
            
    elif len(dtypes) == 2 and int in dtypes and float in dtypes:
        """Allow for the possibility of a mixture of integers and floats."""
        
        """ Convert integers to float: """
        index = [float(item) for item in index]
        
        point = refIm.TransformContinuousIndexToPhysicalPoint(index)
        
    else:
        msg = f"The data type of index is {dtypes}. It must be either "\
              + "'int' or 'float'."
        raise Exception(msg)
         
    #return point # 29/01/2021
    return list(point) # 29/01/2021


""" From dicom_tools.seg_data.py """

def get_pixarr_in_seg_OLD(seg, dicomDir, searchStr):
    """
    Get the frames in a SEG's pixel array that belong to a segment with
    specified label.
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    searchStr : str
        All or part of the Segment Label containing the segment of interest.
                       
    Returns
    -------
    pixarrInSeg : Numpy array
        Sub-array of the SEG's pixel array that contains the frames that belong
        to the segment of interest. 
    f2sIndsInSeg : list of ints
        List (for each segmentation) of the slice numbers that correspond to 
        each frame in pixarr.
    """
    
    pixarr_all = seg.pixel_array
    
    F_all, R, C = pixarr_all.shape
    
    # Get the frame numbers of the SEG's pixel array that correspond to the 
    # segment of interest, and the corresponding Per-frame Functional Groups
    # Sequence-to-slice indices:
    frameNumsInSeg, f2sIndsInSeg = get_frameNums(seg, searchStr, 
                                                           dicomDir)
    
    F = len(frameNumsInSeg)
    
    pixarrInSeg = np.zeros((F, R, C), dtype='uint')
    
    for i in range(F):  
        pixarrInSeg[i] = pixarr_all[frameNumsInSeg[i]]
    
    return pixarrInSeg, f2sIndsInSeg

def get_pixarrBySeg_OLD(seg, dicomDir, p2c=False):
    """
    Get a list of pixel arrays grouped by segment in a SEG.
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object.
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
    f2sIndsBySeg : list of a list of ints
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_pixarr_by_seg():')
        print('\n\n', '-'*120)
        
    pixarr_all = seg.pixel_array
    
    """ If the seg has a single frame it will have shape (NumOfRows, NumOfCols).
    For compatibility with multi-framed pixel arrays, convert single-framed
    pixel arrays to shape (1, NumOfRows, NumOfCols). """
    
    if len(pixarr_all.shape) < 3:
        R, C = pixarr_all.shape
        
        pixarr_all = np.reshape(pixarr_all, (1, R, C))
    
    F_all, R, C = pixarr_all.shape
    
    studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(dicomDir)
    
    f2sIndsBySeg = get_f2sIndsBySeg(seg, SOPUIDs)
    
    Nsegs = len(f2sIndsBySeg)
    
    if p2c:
        print(f'   pixarr_all.shape = {pixarr_all.shape}')
        #print('   f2sIndsBySeg =')
        print_indsByRoi(f2sIndsBySeg)
        
    pixarrBySeg = []
    
    j = 0 # total frame counter for all frames in pixarr_all
    
    for s in range(Nsegs):
        # The number of frames in pixarr_all that belong to this segment:
        F = len(f2sIndsBySeg[s])
        
        # Initialise the pixel array for this segment:
        pixarr = np.zeros((F, R, C), dtype='uint')
        
        if F > 0:
            for i in range(F):
                pixarr[i] = pixarr_all[j]
                
                j += 1 # increment the total frame counter
        
        pixarrBySeg.append(pixarr)
    
    if p2c:
        print_shape_of_pixarrBySeg(pixarrBySeg)
        print('-'*120)
        
    return pixarrBySeg, f2sIndsBySeg


""" From dicom_tools.seg_metadata.py """

def AddToPFFGS(seg):
    """
    Append the last item in the Per-FrameFunctionalGroupsSequence of the
    SEG to itself (i.e. increase the length of the sequence by one). 
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
        
    Returns
    -------
    newSeg : Pydicom Object
        Modified ROI object with Per-FrameFunctionalGroupsSequence lengthened
        by one.
    """
    
    from copy import deepcopy
    
    # Use seg as a template for newSeg: 
    newSeg = deepcopy(seg)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = deepcopy(seg.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    newSeg.PerFrameFunctionalGroupsSequence.append(last)
             
    return newSeg

def AddToSS(seg):
    """
    Append the last item in the SegmentSequence of the SEG Object to itself 
    (i.e. increase the length of the sequence by one). 
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
        
    Returns
    -------
    newSeg : Pydicom Object
        Modified seg with SegmentSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use seg as a template for newSeg: 
    newSeg = deepcopy(seg)
             
    # The last item in Segment Sequence:
    last = deepcopy(seg.SegmentSequence[-1])
    
    # Append to the Segment Sequence:
    newSeg.SegmentSequence.append(last)
             
    return newSeg

# Old functions below

def get_frame_nums_OLD(seg, searchStr, dicomDir):
    """
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Parameters
    ----------
    seg : Pydicom Object
        SEG object.
    searchStr : str
        All or part of the Segment Label containing the segmentation of
        interest.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    
    Returns
    -------
    frameNums : list of ints
        List of the frame numbers that correspond to all frames in the SEG's 
        pixel array that belongs to the segment of interest.  If only one frame 
        corresponds to the segment, the returned list will be of length one 
        (e.g. [3]).
    f2sInds : list of ints
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    from dicom_tools.metadata import get_dcm_uids, get_roicol_nums
    
    # Get the DimensionIndexValues (divs) and group the divs by segment.
    """ Note:  The divs are integers that start from 1. """
    divs = get_DIVs(seg)
    
    # Get the DICOM SOP UIDs:
    studyUID, seriesUID, FORUID, sopuids = get_dcm_uids(dicomDir)
    
    # Get the PerFrameFunctionalGroupsSequence-to-slice indices
    # (f2sInds) and f2sInds grouped by segment
    # (f2sIndsBySeg):
    """ Note:  The f2sInds are integers that start from 0. """
    f2sInds = get_f2sInds(seg, sopuids)
    
    f2sIndsBySeg = group_list_by_seg(listToGroup=f2sInds, 
                                           divs=divs)
    
    segNum = get_roicol_nums(seg, searchStr)
    
    # The f2sInds for the segment of interest:
    f2sInds = f2sIndsBySeg[segNum]
    
    # The number of segments:
    S = len(f2sIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(f2sInds)
    
    n = 0 # initialise frame number counter
    
    frameNums = []
    
    for s in range(S):
        if s == segNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of f2sInds[f] in f2sIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                frameNum = n + f2sIndsBySeg[s].index(f2sInds[f])
                
                frameNums.append(frameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(f2sIndsBySeg[s])
    
    return frameNums, f2sInds


""" From dro_tools.create_dro.py """

def create_dro_OLD_DELETE(srcDataset, trgDataset, newDataset, params):
    # TODO update docstrings
    """
    Create a spatial or deformable DICOM Registration Object if applicable.
    
    If a transform exists as a result of image registration a DRO will be
    created. Otherwise None will be returned. The DRO will be added to 
    propObj.
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    newDataset : Propagator Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    dro : Pydicom Object or None
        DICOM (spatial or deformable) Registration Object if applicable, None
        otherwise.
    """
    
    """ 
    Consider adding this as an input variable for the user to specify.
    Empty for now.
    """
    # ContentDescription:
    contentDesc = ''
        
    useCaseToApply = params.cfgDict['useCaseToApply']
    # Name of the registration transform:
    regTxName = params.cfgDict['regTxName']
    # The registration and pre-registration transforms:
    regTx = newDataset.finalTx
    preRegTx = newDataset.preRegTx
    sampleDroDir = params.cfgDict['sampleDroDir']
    p2c = params.cfgDict['p2c']
    
    if useCaseToApply in ['5a', '5b']:
        timingMsg = "Creating the DICOM Registration Object...\n"
        params.add_timestamp(timingMsg)
        
        if p2c:
            print('regTx =', print(regTx), '\n')
        
        if regTxName in ['rigid', 'affine']:
            # Create a Spatial Registration DRO.
            dro = create_spa_dro_from_tx(
                srcDataset, 
                trgDataset,
                newDataset, 
                params, 
                description=contentDesc
                )
            # TODO Consider adding preRegTx (but need to return it from create_spa_dro_from_tx..) (09/07/21)
        else:
            # Create a Deformable Spatial Registration DRO.
            dro = create_def_dro_from_tx(
                srcDataset, 
                trgDataset,
                newDataset, 
                params, 
                description=contentDesc
                )
        
        timingMsg = "Took [*] s to create the DICOM Registration Object.\n"
        params.add_timestamp(timingMsg)
    else:
        dro = None
    
    return dro


""" From general_tools.config.py """

def create_main_params_file_OLD(cfgDir):
    """
    Generate a JSON file containing parameters for defined configurations of
    different runs.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory to export the config file.
    
    Returns
    -------
    None.
    
    Notes
    -----
    XNAT config settings:
        If the values for the keys "url", "user" or "password" are left empty,    
        the user will be prompted to input it/them.
    
    Resampling / Registration / Transformation settings:
        If the values for the keys "srcFidsFpath" and "trgFidsFpath" are not
        empty, the files containing the fiducials will be used to initialise
        image registration.  If no fiducials are available leave empty.
        
        By default image registration will only be performed if the source and 
        target images have different FrameOfReferenceUIDs.  However it is possible
        for differing ImagePositionPatient and/or ImageOrientationPatient to exist be
        between two images with the same FoR (perhaps down to human error). 
        
        If results of an ROI copy between two datasets with the same FoR arise
        try setting the parameter 'forceReg' to True (it is False by default) to
        force registration to be performed despite equal FoRs.
        
        Image registration can be initialised using one of three methods:
            1. 'landmarks'
            2. 'moments' (or 'centerofgravity')
            3. 'geometry' (or 'geometricalcenter')
            
        If 'landmarks' is chosen but srcFidsFpath and trgFidsFpath are empty or
        not formatted properly, a 'moments'-based initialisation will be performed.
        
        If useDroForTx = False, the transformation parameters will be determined
        based on image registration rather than use of the parameters present in
        the DRO.
        
        Suitable limits for the maximum iterations ("maxIters") are 512 for 
        rigid/affine, and 100 for bspline registrations.
        
        The following three possible values are allowed for the key "resInterp"
        (resampling interpolator):
            1. 'NearestNeighbor'
            2. 'LabelGaussian'
            3. 'BlurThenLinear'
    """
    
    cfg = {}
    
    # Datasets to propagate a SEG using affine registration of soft tissue:
    cfg['Soft_tissue_seg_aff'] = {
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin", 
        'password': "",
        'projID' : 'Soft-tissue-Sarcoma', 
        'subjLab' : 'STS_004',
        'srcExpLab' : 'STS_004_PETCT',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'CT_Ser2_tumour',
        #'srcRoicolName' : 'Left_thigh_tumour',
        'srcRoiName' : None,
        'roicolMod' : 'SEG',
        'trgExpLab' : 'STS_004_MR',
        'trgScanID' : '501',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : os.path.join(cwd, r'STS_004_CT_fiducials.txt'),
        'trgFidsFpath' : os.path.join(cwd, r'STS_004_MR_fiducials.txt'), 
        'forceReg' : False,
        'regTxName' : 'affine',
        'initMethod' : 'landmarks',
        'maxIters' : 512,
        'useDroForTx' : True,
        'applyPreResBlur' : False,
        'preResVar' : (1,1,1),
        'resInterp' : 'BlurThenLinear',
        'applyPostResBlur' : True,
        'postResVar' : (1,1,1)
        }
    
    # Datasets to propagate a RTS using affine registration of soft tissue 
    # (start with copy of 'Soft_tissue_seg_aff' since nearly identical):
    cfg['Soft_tissue_rts_aff'] = dict(cfg['Soft_tissue_seg_aff'])
    cfg['Soft_tissue_rts_aff']['roicolMod'] = 'RTSTRUCT'
    #cfg['Soft_tissue_rts_aff']['srcRoicolName'] = 'CT_Ser2_tumour'
    
    # Datasets to propagate SEG using bspine registration of soft tissue:
    cfg['Soft_tissue_seg_def'] = dict(cfg['Soft_tissue_seg_aff'])
    cfg['Soft_tissue_seg_def']['regTxName'] = 'bspline'
    cfg['Soft_tissue_seg_def']['maxIters'] = 100
    
    # Datasets to propagate RTS using bspine registration of soft tissue:
    cfg['Soft_tissue_rts_def'] = dict(cfg['Soft_tissue_seg_aff'])
    cfg['Soft_tissue_rts_def']['roicolMod'] = 'RTSTRUCT'
    
    
    
    """
    Following are datasets for test runs using the NCITA_TEST collection (brain
    images) to either propagate ("relationship-preserving copy") or copy 
    ("direct copy") SEG and RTSTRUCT ROI Collections using affine registrations:
    
    https://docs.google.com/document/d/1_LxC1FUdQME4xuXZFGjjGBGYRb1cJX8YiEADjtZTn_g/edit#
    """
    
    # Datasets to propagate RTS within the same DICOM series:
    runID = 'NCITA_test_RR1'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin", 
        'password': "",
        'projID' : 'NCITA_TEST', 
        'subjLab' : 'TCGA-BB-A5HY',
        'srcExpLab' : 'Session2',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'Left eye CT',
        'srcRoiName' : None,
        'roicolMod' : 'RTSTRUCT',
        'trgExpLab' : 'Session2',
        'trgScanID' : '4',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : None,
        'trgFidsFpath' : None, 
        'forceReg' : False,
        'regTxName' : 'affine',
        'initMethod' : 'moments',
        'maxIters' : 512,
        'useDroForTx' : False,
        'applyPreResBlur' : False,
        'preResVar' : (1,1,1),
        'resInterp' : 'BlurThenLinear',
        'applyPostResBlur' : True,
        'postResVar' : (1,1,1)
        }
    
    # Datasets to propagate RTS across two DICOM series with the same slice
    # thickness:
    runID = 'NCITA_test_RR2'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session4'
    cfg[runID]['srcScanID'] = '13'
    cfg[runID]['srcRoicolName'] = 'Nasal cavity'
    cfg[runID]['trgExpLab'] = 'Session4'
    cfg[runID]['trgScanID'] = '10'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR3'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS from one DICOM series to another with smaller
    # slice thickness:
    runID = 'NCITA_test_RR4'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '11'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR5'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '14'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR6'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '11'
    cfg[runID]['srcRoicolName'] = 'Ventricles'
    cfg[runID]['trgScanID'] = '10'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR7'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '12'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR8'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '13'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR9'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '14'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR10'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session3'
    cfg[runID]['trgScanID'] = '6'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR11'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '13'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR12'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session5'
    cfg[runID]['srcScanID'] = '12'
    cfg[runID]['srcRoicolName'] = 'Left eye'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR13'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session1'
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR14'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '13'
    cfg[runID]['srcRoicolName'] = 'Nasal cavity'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR15'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session5'
    cfg[runID]['srcScanID'] = '12'
    cfg[runID]['srcRoicolName'] = 'Left eye'
    
    # Dataset to propagate RTS to an existing RTS:
    """ 
    Note:
        This was not one of the original test configurations so unknown whether
    this will work.
    """
    runID = 'NCITA_test_RR16'
    cfg[runID] = dict(cfg['NCITA_test_RR14'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session4'
    cfg[runID]['trgScanID'] = '11'
    cfg[runID]['trgRoicolName'] = 'Ventricles'
    
    # Datasets to copy RTS within the same DICOM series:
    runID = 'NCITA_test_RD1'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = 5
    cfg[runID]['srcRoicolName'] = 'Right eye CT (single slice)'
    cfg[runID]['trgScanID'] = '2'
    cfg[runID]['trgSlcNum'] = 6
    
    # Datasets to copy RTS within the same DICOM series:
    runID = 'NCITA_test_RD2'
    cfg[runID] = dict(cfg['NCITA_test_RD1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '4'
    
    # Datasets to copy RTS across two DICOM series with different IOPs:
    """
    Note:
        Source slice 5 maps to Target slices 55-64.
    """
    runID = 'NCITA_test_RD3'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    cfg[runID]['trgSlcNum'] = 65
    
    # Datasets to copy RTS from one DICOM series to another with equal slice
    # thickness:
    """
    Note:
        Source slice 5 maps to Target slice 7.
    """
    runID = 'NCITA_test_RD4'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session6'
    cfg[runID]['trgScanID'] = '3'
    cfg[runID]['trgSlcNum'] = 8
    
    # Datasets to copy RTS from one DICOM series to another with smaller
    # slice thickness:
    """
    Note:
        Source slice 5 maps to Target slices 67-77.
    """
    runID = 'NCITA_test_RD5'
    cfg[runID] = dict(cfg['NCITA_test_RD4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    cfg[runID]['trgSlcNum'] = 78
    
    # Datasets to copy RTS across two DICOM series with equal slice
    # thickness, demonstrating the displacement within the same series:
    runID = 'NCITA_test_RD6'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = 18
    cfg[runID]['trgSlcNum'] = 21
    
    # Datasets to propagate SEG within the same DICOM series:
    runID = 'NCITA_test_SR1'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with the same slice
    # thickness:
    runID = 'NCITA_test_SR2'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR3'
    cfg[runID] = dict(cfg['NCITA_test_RR3'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with smaller
    # slice thickness:
    runID = 'NCITA_test_SR4'
    cfg[runID] = dict(cfg['NCITA_test_RR4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR5'
    cfg[runID] = dict(cfg['NCITA_test_RR5'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR6'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR7'
    cfg[runID] = dict(cfg['NCITA_test_RR7'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR8'
    cfg[runID] = dict(cfg['NCITA_test_RR8'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR9'
    cfg[runID] = dict(cfg['NCITA_test_RR9'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR10'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR11'
    cfg[runID] = dict(cfg['NCITA_test_RR11'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR12'
    cfg[runID] = dict(cfg['NCITA_test_RR12'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR13'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR14'
    cfg[runID] = dict(cfg['NCITA_test_RR14'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR15'
    cfg[runID] = dict(cfg['NCITA_test_RR15'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG within the same DICOM series:
    runID = 'NCITA_test_SD1'
    cfg[runID] = dict(cfg['NCITA_test_RD1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG within the same DICOM series:
    runID = 'NCITA_test_SD2'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG across two DICOM series with different IOPs:
    """
    Note:
        Source slice 5 maps to Target slices 55-64.
    """
    runID = 'NCITA_test_SD3'
    cfg[runID] = dict(cfg['NCITA_test_RD3'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG from one DICOM series to another with equal slice
    # thickness:
    """
    Note:
        Source slice 5 maps to Target slice 7.
    """
    runID = 'NCITA_test_SD4'
    cfg[runID] = dict(cfg['NCITA_test_RD4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG from one DICOM series to another with smaller
    # slice thickness:
    """
    Note:
        Source slice 5 maps to Target slices 67-77.
    """
    runID = 'NCITA_test_SD5'
    cfg[runID] = dict(cfg['NCITA_test_RD5'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG across two DICOM series with equal slice
    # thickness, demonstrating the displacement within the same series:
    runID = 'NCITA_test_SD6'
    cfg[runID] = dict(cfg['NCITA_test_RD6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # TODO Verify the above NCITA_TEST configs and comments
    
    export_dict_to_json(dictionary=cfg, filename='main_params.json',
                        exportDir=cfgDir)
    
def add_password_to_main_params_file_OLD(cfgDir, password):
    
    fpath = os.path.join(cfgDir, r'main_params.json')
    
    cfg = import_dict_from_json(fpath)
    
    for runID in cfg.keys():
        cfg[runID]['password'] = password
    
    export_dict_to_json(dictionary=cfg, filename='main_params.json',
                        exportDir=cfgDir)
    
def create_general_params_file_OLD(cfgDir):
    """
    Create general parameters.
    
    It is assumed that the general parameters will not likely be changed 
    between runs with different datasets and/or resampling / registering / 
    transforming settings. The file general_params.json will be exported to a 
    directory named "config" within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory to export the config file.
    
    Returns
    -------
    None.
    
    Notes
    -----
    If uploadDro = True, the DRO will be uploaded to XNAT as a subject 
    assessor.
    
    If overwriteDro = True (and uploadDro = True) and there exists a DRO
    (as a subject assessor) in XNAT, the existing DRO will be over-written.
    """
    
    
    """ 
    Select one of the following suggested download directories for XNAT data:
    """
    #rootExportDir = ddd
    rootExportDir = os.path.join(ddd, currentDate)
    #rootExportDir = os.path.join(ddd, '2021-07-12') # change accordingly
    
    """
    Define the export directories for the (new) target ROI Collections, 
    DRO, plots, logs and configuration files:
    """
    rtsExportDir = os.path.join(rootExportDir, 'new_RTS')
    segExportDir = os.path.join(rootExportDir, 'new_SEG')
    
    droExportDir = os.path.join(rootExportDir, 'new_DRO')
    
    rtsPlotsExportDir = os.path.join(rootExportDir, 'plots_RTS')
    segPlotsExportDir = os.path.join(rootExportDir, 'plots_SEG')
    
    logsExportDir = os.path.join(rootExportDir, 'logs')
    
    #cfgExportDir = os.path.join(rootExportDir, 'configs')
    
    """
    Chose whether or not to export the (new) target ROI Collection (i.e.
    RTS or SEG), DRO, plots and logs:
    """
    #exportRoicol = False
    exportRoicol = True
    
    #exportDro = False
    exportDro = True
    
    #exportPlots = False
    exportPlots = True
    #exportLogs = False
    exportLogs = True
    
    
    """
    Chose whether or not to upload the (new) target DRO to XNAT:
    """
    #uploadDro = False
    uploadDro = True # see Notes
    
    """
    Chose whether or not to overwrite an existing DRO to XNAT:
    """
    overwriteDro = False
    #overwriteDro = True # see Notes
    
    """
    Chose whether or not to print results to console (e.g. for debugging):
    """
    #p2c = False
    p2c = True
    
    cfg = {'rootExportDir' : rootExportDir,
           #'cfgExportDir' : cfgExportDir,
           'rtsExportDir' : rtsExportDir,
           'segExportDir' : segExportDir,
           'droExportDir' : droExportDir,
           'rtsPlotsExportDir' : rtsPlotsExportDir,
           'segPlotsExportDir' : segPlotsExportDir,
           'logsExportDir' : logsExportDir,
           'exportTrgRoicol' : exportTrgRoicol,
           'exportDro' : exportDro,
           'exportPlots' : exportPlots,
           'exportLogs' : exportLogs,
           'uploadDro' : uploadDro,
           'overwriteDro' : overwriteDro,
           'p2c' : p2c}
    
    export_dict_to_json(dictionary=cfg, filename='general_params.json',
                        exportDir=cfgDir)

def create_xnat_cfg_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create a configuration file for XNAT connections.
    
    The file xnat_config.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Note
    ----
    If the values for keys "user" and "password" are left empty, i.e.
    
    cfg = {"url" : "http://10.1.1.20", 
           "user": "", 
           "password": ""}
    
    the user will be prompted to enter the credentials.
    """
    
    cfg = {"url" : "http://10.1.1.20", 
           "username": "admin", 
           "password": "admin"}
    
    export_dict_to_json(dictionary=cfg, filename='xnat_config.json',
                        exportDir=cfgDir)

def create_src_xnat_params_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create a parameters file for XNAT parameters for Source.
    
    The file xnat_params.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Notes
    -----
    roicol => ROI Collection (i.e. RTS or SEG)
    
    projID = XNAT project ID
    subjLab = XNAT subject label
    expLab = XNAT experiment label
    scanID = XNAT scan ID
    slcNum = Slice number within a colection of DICOM image stack
    roicolMod = Modality of the ROI Collection
    roicolName = XNAT name of the ROI Collection
    roiName = DICOM ROIName or SegmentLabel of an ROI Collection
    """
    
    cfg = {'projID' : 'Soft-tissue-Sarcoma', 
           'subjLab' : 'STS_004',
           'expLab' : 'STS_004_PETCT',
           'scanID' : '2',
           'slcNum' : None,
           'roicolMod' : 'RTSTRUCT',
           'roicolName' : 'CT_Ser2_tumour',
           'roiName' : None}
    
    export_dict_to_json(dictionary=cfg, filename='src_xnat_params.json',
                        exportDir=cfgDir)

def create_trg_xnat_params_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create a parameters file for XNAT parameters for Target.
    
    The file xnat_params.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Notes
    -----
    roicol => ROI Collection (i.e. RTS or SEG)
    
    projID = XNAT project ID
    subjLab = XNAT subject label
    expLab = XNAT experiment label
    scanID = XNAT scan ID
    slcNum = Slice number within a colection of DICOM image stack
    roicolMod = Modality of the ROI Collection
    roicolName = XNAT name of the ROI Collection
    roiName = DICOM ROIName or SegmentLabel of an ROI Collection
    """
    
    cfg = {'projID' : 'Soft-tissue-Sarcoma', 
           'subjLab' : 'STS_004',
           'expLab' : 'STS_004_MR',
           'scanID' : '501',
           'slcNum' : None,
           'roicolMod' : 'RTSTRUCT',
           'roicolName' : None,
           'roiName' : None}
           #'addToTrgRoicolName' : f' from {srcRoicolName}'}
    
    export_dict_to_json(dictionary=cfg, filename='trg_xnat_params.json',
                        exportDir=cfgDir)

def create_res_params_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create resampling / registration / transformation-related parameters.
    
    The file res_params.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Notes
    -----
    By default image registration will only be performed if the source and 
    target images have different FrameOfReferenceUIDs.  However it is possible
    for differing ImagePositionPatient and/or ImageOrientationPatient to exist be
    between two images with the same FoR (perhaps down to human error). 
    
    If results of an ROI copy between two datasets with the same FoR arise
    try setting the parameter 'forceReg' to True (it is False by default) to
    force registration to be performed despite equal FoRs.
    
    Image registration can be initialised using one of three methods:
        1. 'landmarks'
        2. 'moments' (or 'centerofgravity')
        3. 'geometry' (or 'geometricalcenter')
        
    If 'landmarks' is chosen but srcFidsFpath and trgFidsFpath are empty or
    not formatted properly, a 'moments'-based initialisation will be performed.
    
    If useDroForTx = False, the transformation parameters will be determined
    based on image registration rather than use of the parameters present in
    the DRO.
    """
    
    
    """ 
    Define the filepaths to the source and target fiducials txt files (if
    applicable):
    """
    #srcFidsFpath = '' # no fiducials
    #trgFidsFpath = '' # no fiducials
    srcFidsFpath = os.path.join(cwd, r'STS_004_CT_fiducials.txt')
    trgFidsFpath = os.path.join(cwd, r'STS_004_MR_fiducials.txt')
    
    
    """ 
    Define registration parameters:
    """
    forceReg = False
    #forceReg = True # see Notes
    
    #regTxName = 'rigid'
    regTxName = 'affine'
    #regTxName = 'bspline'
    
    initMethod = 'landmarks' # see Notes
    #initMethod = 'moments'
    #initMethod = 'geometry'
    
    maxIters = 512 # suitable for rigid/affine
    #maxIters = 100 # suitable for bspline
    
    
    """ 
    Define resampling / transforming parameters:
    """
    useDroForTx = True
    #useDroForTx = False # see Notes
    
    applyPreResBlur = False
    #applyPreResBlur = True
    
    preResVar = (1,1,1) # variance
    
    #resInterp = 'NearestNeighbor'
    #resInterp = 'LabelGaussian'
    resInterp = 'BlurThenLinear'
    
    #applyPostResBlur = False
    applyPostResBlur = True
    
    postResVar = (1,1,1) # variance
    
    cfg = {'srcFidsFpath' : srcFidsFpath,
           'trgFidsFpath' : trgFidsFpath, 
           'forceReg' : forceReg,
           'regTxName' : regTxName,
           'initMethod' : initMethod,
           'maxIters' : maxIters,
           'useDroForTx' : useDroForTx,
           'applyPreResBlur' : applyPreResBlur,
           'preResVar' : preResVar,
           'resInterp' : resInterp,
           'applyPostResBlur' : applyPostResBlur,
           'postResVar' : postResVar}
    
    export_dict_to_json(dictionary=cfg, filename='reg_res_tx_params.json',
                        exportDir=cfgDir)
    
    
""" From general_tools.general.py """

def are_items_equal_to_within_eps_REDUNDANT(item0, item1, epsilon=1e-06):
    """
    07/07/21:
        This is redundant - use are_lists_equal_to_within_eps(), which was
        latter named are_items_equal_to_within_eps.
        
    Are two items equal to within epsilon?
    
    Parameters
    ----------
    item0 : float or int
    item1 : float or int
    epsilon : float, optional
        Maximum allowed variation between items not considered to be unique. 
        The default is 1e-5.

    Returns
    -------
    is_equal : bool
        True if the items are equal to within epsilon.
    """
    
    abs_diff = abs(item0 - item1)
    
    if abs_diff < epsilon:
        is_equal = True
    else:
        is_equal = False
        
    return is_equal


""" From general_tools.shifting.py """

def shift_frames_in_pixarrBySeg_OBSOLETE(PixArrBySeg, F2SindsBySeg, srcIm, 
                                srcSlcNum, trgIm, trgSlcNum, refImage, 
                                shiftInX=True, shiftInY=True, shiftInZ=True, 
                                fractional=False, p2c=False):
    """
    See shift_frames_in_pixarrBySeg in propagate.propagate.py
    
    Shift the items in each frame of each pixel array in a list of 
    pixel-arrays-by-segment.  Example uses: To compensate for differences in 
    z-positions of a single-framed Source pixel array to be "direct" copied to
    a Target slice at a different stack position; or compensating for 
    differences in the origin of Source and Target images for relationship-
    preserving copies of images with equal voxel spacings.
    
    Parameters
    ----------
    PixArrBySeg : list of Numpy arrays
        A list (for each segment) of pixel arrays.
    F2SindsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of slice numbers that 
        correspond to each frame in the pixel arrays in PixArrBySeg.
    srcIm : SimpleITK image
        The Source 3D image.
    srcSlcNum : int
        The slice number within the Source image stack that the segment relates
        to.
    trgIm : SimpleITK image
        The Target 3D image.
    trgSlcNum : int
        The slice number within the Target image stack that the segment is to  
        be copied to following the shift in z-components.
    refImage : SimpleITK image
        The reference image whose voxel spacings will be used to convert apply
        the shift in pixels (e.g. trgIm).
    shiftInX / shiftInY / shiftInZ : bool, optional (True by default)
        If True the pixel shift between the origins of srcIm and trgIm
        will be applied along the x / y / z direction.
    fractional : bool, optional (False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    PixArrBySeg : list of Numpy arrays
        As above but with shifted frame with shape (1, R, C).
    F2SindsBySeg : list of a list of ints
        As above but with shifted slice indices that correspond to the shifted 
        PixArrBySeg.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of shift_frames_in_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    # Get pixel shift between FromSliceNum in srcIm and ToSliceNum in 
    # trgIm in the Target image domain:
    voxShift = get_voxel_shift_bt_slices(
        image0=srcIm, sliceNum0=srcSlcNum, image1=trgIm, 
        sliceNum1=trgSlcNum, refImage=trgIm
        )
    
    if p2c:
        print(f'   \nvoxShift between slices = {voxShift}')
        print(f'   \nF2SindsBySeg prior to shifting = {F2SindsBySeg}')
        
        #print(f'\n   Prior to shifting frames:')
        ##print(f'PixArrBySeg = {PixArrBySeg}')
        #print(f'   len(PixArrBySeg) = {len(PixArrBySeg)}')
        #for r in range(len(PixArrBySeg)):
        #    print(f'      type(PixArrBySeg[{r}]) = {type(PixArrBySeg[r])}')
        #    print(f'      PixArrBySeg[{r}].shape = {PixArrBySeg[r].shape}')
    
    if not shiftInX:
        voxShift[0] = 0
    if not shiftInY:
        voxShift[1] = 0
    if not shiftInZ:
        voxShift[2] = 0
    
    if p2c:
        print(f'   \nvoxShift that will be applied = {voxShift}')
        
    for s in range(len(PixArrBySeg)):
        #if PixArrBySeg[s]:
        if PixArrBySeg[s].shape[0]:
            #print('   \nBefore shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            # Replace PixArrBySeg[s] with the result of shifting the in-plane 
            # elements, and replace F2SindsBySeg[s] with [trgSlcNum]:
            PixArrBySeg[s] = shift_frame(
                frame=PixArrBySeg[s], voxShift=voxShift
                )
            
            #print('   \nAfter shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            #F2SindsBySeg[s] = [trgSlcNum] # 11/02
            
            # Shift the frame-to-slice indices by voxShift[2]:
            F2SindsBySeg[s] = [F2SindsBySeg[s][i] + voxShift[2] for i in range(len(F2SindsBySeg[s]))]
            
            if p2c:
                unique = get_unique_items(Items=PixArrBySeg[s], IgnoreZero=False)
                
                print(f'\nThere are {len(unique)} unique items in',
                      f'PixArrBySeg[{s}] after shifting the frame.')
    
    if p2c:
        print(f'   F2SindsBySeg after shifting = {F2SindsBySeg}')
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg

def shift_ptsByCntByRoi_210914(
        ptsByCntByRoi, c2sIndsByRoi, srcIm, srcSlcNum, 
        trgIm, trgSlcNum, refImage, 
        shiftInX=True, shiftInY=True, shiftInZ=True, 
        fractional=False, p2c=False
        ):
    """
    14/09/21: 
    
    Shift the points in each list of points in each list of contours in each
    list of ROIs.  
    
    Example uses: 
        - to compensate for differences in z-positions of a single-contour 
        Source points to be "direct" copied to a Target slice at a different  
        stack position
        - compensating for differences in the origin of Source and Target 
        images for relationship-preserving copies of images with equal voxel 
        spacings
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of indices that denote
        the slice number that the contour relates to.
    newZind : int
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    dicomDir : str
        Directory containing the corresponding DICOMs.
    
    Returns
    -------
    shiftedPtsByCntByRoi : list of list of a list of a list of floats
        As ptsByCntByRoi but with modified z-components.
    """
    
    #from ImageTools import ImportImage
    #from ConversionTools import Points2Indices, Indices2Points
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of shift_ptsByCntByRoi():')
        print('\n\n', '-'*120)
        
    """ Get voxel shift between FromSliceNum in srcIm and ToSliceNum in 
    trgIm in the Target image domain: """
    voxShift = get_voxel_shift_bt_slices(
        image0=srcIm, sliceNum0=srcSlcNum, 
        image1=trgIm, sliceNum1=trgSlcNum,
        refImage=trgIm
        )
    
    if p2c:
        print(f'   voxShift between slices = {voxShift}')
    
    if not shiftInX:
        voxShift[0] = 0
    if not shiftInY:
        voxShift[1] = 0
    if not shiftInZ:
        voxShift[2] = 0
    
    if p2c:
        print(f'   voxShift that will be applied = {voxShift}')
    
    
    """ Convert ptsByCntByRoi to indices, modify the indices, then convert back
    to physical points. """
    
    #indsByCntByRoi = []
    shiftedPtsByCntByRoi = []
    shiftedC2SindsByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(ptsByCntByRoi)): 
        #indsByCnt = []
        shiftedPtsByCnt = []
        shiftedC2Sinds = []
        
        if ptsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(ptsByCntByRoi[r])):
                pts = ptsByCntByRoi[r][c]
                
                """ inds are the list of relationship-preserving indices that
                correspond to pts: """
                inds = pts_to_inds(points=pts, refIm=trgIm, rounding=False)
                
                #print(f'\n\nvoxShift = {voxShift}, \ninds = {inds}')
                
                """ Apply the required pixel shifts: """
                inds = [inds[i] + voxShift[i] for i in range(len(voxShift))]
                
                #indsByCnt.append(inds)
                #print(f'\n\ntype(inds) = {type(inds)}')
                #print(f'\n\ntype(inds[0]) = {type(inds[0])}')
                
                """ Convert back to physical points: """
                pts = inds_to_pts(inds, refIm=trgIm)
                
                shiftedPtsByCnt.append(pts)
                
                """ Shift the contour-to-slice indices by voxShift[2]: """
                #c2sIndsByRoi[r][c] += voxShift[2]
                shiftedC2Sinds.append(c2sIndsByRoi[r][c] + voxShift[2])
        
        #indsByCntByRoi.append(indsByCnt)
        shiftedPtsByCntByRoi.append(shiftedPtsByCnt)
        shiftedC2SindsByRoi.append(shiftedC2Sinds)
        
    return shiftedPtsByCntByRoi, shiftedC2SindsByRoi


""" From image_tools.resampling.py """

def resample_labim_OLD(labim, f2sInds, srcIm, trgIm, interp, 
                   preResVar=(1,1,1), applyPostResBlur=True, 
                   postResVar=(1,1,1), p2c=False):
    """
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------  
    labim : SimpleITK image
        The 3D label image to be resampled.
    f2sInds : list of ints
        List (for each frame) of slice numbers that correspond to each frame in 
        labim.
    srcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of labim. 
    trgIm : SimpleITK image
        The Target 3D image whose gridspace the labim will be resampled to.
    interp : str
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    preResVar : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling.
    applyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    postResVar : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    resLabim : SimpleITK images
        The resampled 3D label image.
    resPixarr : Numpy arrays
        The resampled pixel array (converted from resLabim).
    resF2Sinds : list of ints
        The list (for each frame) of the frame-to-slice indices in resPixarr.
    
    Note
    ----
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (label) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled label,
    and hence an empty list of frame-to-slice indices, f2sInds = []) even if 
    there were contours (i.e. if SrcC2Sinds != []). If this is the case try
    suggestions from Ziv Yaniv (See Link below):
            
    1. Use the sitkLabelGaussian interpolator.
    
    2. Gaussian blur the segmentation, then resample using a linear 
    interpolator, and then threshold so that values in [a<1<b] are mapped to 1.
    You will need to select appropriate values for a, b.
    
    3. Compute the distance map of the segmentation, then resample using a
    linear interpolator, and then threshold the absolute value of the distance
    map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you need to 
    select an appropriate dist_threshold.
    
    What has been tried so far:
    
    1. Resulting Gaussian interpolated resampled image had zeros only 
    (the image was converted to float prior to resampling).
    
    2. Works.
    
    Comment from Bradley Lowecamp:
        
    The implementation of the sitkLabelGaussian interpolator in the 
    ResampleImageFilter sets the sigma to the pixel spacing of the input image 
    to the resample filter. This is why it does not do a better job in reducing 
    aliasing.
    
    Link: https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of resample_labim():')
        print('\n\n', '-'*120)
        
    if not interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f'The chosen interpolation, {interp}, is not one of the '\
              + 'accepted inputs: \'NearestNeighbor\', \'LabelGaussian\', or '\
              + '\'BlurThenLinear\'.'
        
        raise Exception(msg)
    
    # Store the interpolation set as metadata:
    labim.SetMetaData("resInterpSet", interp)
    
    #srcSpacings = srcIm.GetSpacing()
    #trgSpacings = trgIm.GetSpacing()
    
    #spacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(srcSpacings, 
    #                                                      trgSpacings))
    
    #srcSpacings = np.array(srcIm.GetSpacing())
    #trgSpacings = np.array(trgIm.GetSpacing())
        
    #spacingsRatio = np.divide(trgSpacings, srcSpacings)
    
    #volumeRatio = np.prod(spacingsRatio)
    
    #""" 
    #The FWHM of the Gaussian is:
    #    fwhm = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
    #    
    #Due to the Nyquist-Shannon sampling theorem, the FWHM should be at
    #least 2x the voxel spacings of trgIm. For symmetry it makes sense for
    #one Source voxel width to be blurred to three Target voxels widths:
    #    sigma = (3*trgSpacings)/2.355
    #"""
    
    #if p2c:
    #    print(f'\nsrcSpacings = {srcSpacings}')
    #    print(f'trgSpacings = {trgSpacings}')
    #    print(f'\nspacingsRatio = {spacingsRatio}')
    
    if interp in ['NearestNeighbor', 'LabelGaussian']:
        if p2c:
            print('Attempting to resample labim using the',
                  f'{interp} interpolator...\n')
        
        """ Attempt to resample labim to the Target image's grid using the
        chosen labelmap interpolator. """
        resLabim = resample_im(im=labim, refIm=trgIm, interp=interp)
        
        if p2c:
            print('Image info for resLabim after resampling using',
                  f'{interp} interpolator:')
        pixID, pixIDTypeAsStr, uniqueVals,\
        resF2Sinds = get_im_info(resLabim, p2c)
        if p2c:
            print('')
        
        if applyPostResBlur:
            # Gaussian blur the image.
            resLabim = gaussian_blur_im(im=resLabim, var=postResVar)
            
            if p2c:
                print('Image info for resLabim after Gaussian blurring:')
            pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
                resLabim, p2c
                )
            if p2c:
                print('')
            
            #""" The original binary pixel array: """
            #PixArr_B, f2sInds = im_to_pixarr(labim)
            #
            #""" The non-binary pixel array following blur: """
            #pixArr_NB, f2sInds = im_to_pixarr(resLabim)
            #
            #ratio = GetVoxelVolumeRatio(labim, resLabim)
            #
            #if p2c:
            #    print(f'pixArr_B.shape = {pixArr_B.shape}')
            #    print(f'pixArr_NB.shape = {pixArr_NB.shape}\n')
            
            # Binarise the resampled labelmap if required.
            if len(uniqueVals) != 2 or sum(uniqueVals) != 1:
                """Find suitable threshold value that approximately preserves
                the number of pre-blurred truth values scaled by volumeRatio: 
                """
                thresh = find_thresh(
                    binaryIm=labim, nonBinaryIm=resLabim, p2c=p2c
                    )
                
                # Binary threshold the image.
                resLabim = binarise_im(
                    im=resLabim, #thresh=postResThresh)
                    thresh=thresh) 
                
                if p2c:
                    print('\nImage info for resLabim after binary',
                          #f'thresholding at {postResThresh}:')
                          f'thresholding at {thresh}:')
                pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                    resLabim, p2c
                    )
                if p2c:
                    print('')
        
        #print(f'\n   resF2Sinds = {resF2Sinds}')
        
        """ Is resF2Sinds empty and not expected to be? If so, try the 
        "BlurThenLinear" approach.
        
        Note:
            If f2sInds isn't empty, resF2Sinds shouldn't be empty either.
        
            f2sInds will be empty if there were no segmentations/contours 
            of interest for the r^th ROI. In this case an empty resF2Sinds 
            is acceptable. """
        
        if resF2Sinds == []:
            print(f'There are {len(f2sInds)} non-empty masks in the input',
                  f'labelmap image but {len(resF2Sinds)} non-empty frames in',
                  f'the resampled labelmap image using {interp}.',
                  'Try using the \'BlurThenLinear\' approach.\n')
            
            interp = 'BlurThenLinear'

    if interp == 'BlurThenLinear':
        if p2c:
            print('\nResampling labim by applying Gaussian blur,',
                  'resampling using a linear interpolator, then binary',
                  'thresholding...')
        
        #""" The original binary pixel array: """
        #pixArr_B, f2sInds = im_to_pixarr(labim)
        
        """ Convert labim from 32-bit unsigned integer to float.
        Note:  Result of resampling results in empty labelmap unless
        image is converted to float prior to resampling. """
        labim = change_im_dtype(im=labim, newPixType='Float32')
        
        if p2c:
            print('\nImage info for labim after converting to float:')
            pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                labim, p2c
                )
              
        # Gaussian blur the image.
        labim = gaussian_blur_im(im=labim, variance=preResVar)
        
        # Use the RecursiveGaussian image filter:
        #resLabim = recursive_gaussian_blur_im(
        #    resLabim, sigma=3, direction=2
        #    )
        
        if p2c:
            print('\nImage info for labim after applying Gaussian blur:')
            pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                labim, p2c
                )
        
        if p2c:
            print('\nlabim prior to resampling:')
            print(f'   labim.GetSize() = {labim.GetSize()}')
            print(f'   labim.GetSpacing() = {labim.GetSpacing()}')
            print(f'   trgIm.GetSize() = {trgIm.GetSize()}')
            print(f'   trgIm.GetSpacing() = {trgIm.GetSpacing()}')
        
        # Linearly resample labim to the Target image's grid.
        resLabim = resample_im(im=labim, refIm=trgIm, interp='Linear')
        
        if p2c:
            print('\nImage info for resLabim after resampling using a',
                  'linear interpolator:')
            
        pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
            resLabim, p2c
            )
        
        if applyPostResBlur:
            # Gaussian blur the image.
            resLabim = gaussian_blur_im(
                im=resLabim, var=postResVar
                )
        
        #""" The non-binary pixel array following blur + linear resampling: """
        #pixArr_NB, f2sInds = im_to_pixarr(resLabim)
        #
        #if p2c:
        #        print(f'pixArr_B.shape = {pixArr_B.shape}')
        #        print(f'pixArr_NB.shape = {pixArr_NB.shape}\n')
            
        """Find suitable threshold value that approximately preserves the
        number of pre-blurred + linear resampled truth values scaled by
        volumeRatio: """
        thresh = find_thresh(
            binaryIm=labim, nonBinaryIm=resLabim, p2c=p2c
            )
        
        # Binary threshold the image. 
        resLabim = binarise_im(
            im=resLabim, #thresh=postResThresh)
            thresh=thresh
            ) 
        
        if p2c:
            print('\nImage info for resLabim after binary thresholding',
                  #f'at {postResThresh}:')
                  f'at {thresh}:')
        pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
            resLabim, p2c
            )
    
    # Ensure that resLabim is a 32-bit unsigned integer (pixID = 5).
    if pixID != 5: 
        if p2c:
            print(f'\nresLabim has PixelID = {pixID} ({pixIDTypeAsStr})).')
        
        # Convert resLabim from float to 32-bit unsigned integer. 
        resLabim = change_im_dtype(im=resLabim, newPixType='UInt32')
        
        if p2c:
            print('\nImage info for resLabim after converting to 32-bit',
                  'unsigned int:')
        #print(f'\nThe metadata keys are:', resLabim.GetMetaDataKeys())
        pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
            resLabim, p2c
            )
    
    # Convert resLabim to a pixel array. 
    resPixarr, resF2Sinds = im_to_pixarr(resLabim)
    
    # Store the interpolation used as metadata (which may be the same or 
    # different from the interpolation set): 
    resLabim.SetMetaData("resInterpUsed", interp)
    
    if p2c:
        print(f'\nThe key "resInterpUsed" with value "{interp}" was',
              'added to the metadata of resLabim. \nThe metadata keys are:',
              resLabim.GetMetaDataKeys())
    
    if interp == 'BlurThenLinear':
        """ Store the threshold used as metadata: """
        resLabim.SetMetaData("postResThreshUsed", f"{thresh}")
        
        if p2c:
            print(f'\nThe key "postResThreshUsed" with value "{thresh}" was',
                  'added to the metadata of resLabim. \nThe metadata keys:',
                  resLabim.GetMetaDataKeys())
            
    if p2c:
        print('\nAfter converting resLabim to a pixel array:')
        print(f'resPixarr.shape = {resPixarr.shape}')
        print(f'resF2Sinds = {resF2Sinds}')
        print('-'*120)
        
    return resLabim, resPixarr, resF2Sinds


""" From io_tools.download_data.py """

    def __init__210817(self, cfgDict, xnatSession=None):
        
        p2c = cfgDict['p2c']
        
        self.cfgDict = cfgDict
        
        if xnatSession == None:
            # Establish XNAT connection:
            xnatSession = create_session(xnatCfg=cfgDict)
            
        self.xnatSession = xnatSession
        
        # Download data and create pathsDict:
        self.get_pathsDict(xnatSession=xnatSession)
        #self.pathsDict = self.get_pathsDict(xnatSession)
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(config=self.cfgDict)
        
        if p2c:
            print(f"cfgDict['runID'] = {cfgDict['runID']}")
            print(f'useCaseThatApplies = {useCaseThatApplies}')
            print(f'useCaseToApply = {useCaseToApply}')
        
        # Update cfgDict:
        self.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        self.cfgDict['useCaseToApply'] = useCaseToApply
    
    def get_pathsDict_Pre_210809(self, xnatSession):
        """
        Downloads scans and ROI Collections and returns a dictionary containing
        filepaths.
        
        Parameters
        ----------
        params : Params Object
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        xnatCfg = self.xnatCfg
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        genParams = self.genParams
        
        pathsDict = download_scan(xnatCfg=xnatCfg, 
                                  xnatParams=srcXnatParams, 
                                  genParams=genParams,
                                  xnatSession=xnatSession)
    
        pathsDict = download_scan(xnatCfg=xnatCfg, 
                                  xnatParams=trgXnatParams, 
                                  genParams=genParams, 
                                  xnatSession=xnatSession,
                                  pathsDict=pathsDict)
    
        pathsDict = download_im_asr(xnatCfg=xnatCfg, 
                                    xnatParams=srcXnatParams,
                                    genParams=genParams, 
                                    xnatSession=xnatSession,
                                    pathsDict=pathsDict)
        
        trgRoicolName = trgXnatParams['roicolName']
        
        if trgRoicolName != None:
            pathsDict = download_im_asr(xnatCfg=xnatCfg, 
                                        xnatParams=trgXnatParams,
                                        genParams=genParams, 
                                        xnatSession=xnatSession,
                                        pathsDict=pathsDict)
        
        #return pathsDict
        self.pathsDict = pathsDict
        
    def add_dirs_and_fpaths_210810(self):
        """
        Add directories and filepaths to params.srcXnatParams and 
        params.trgXnatParams from pathsDict.

        Returns
        -------
        None.

        """
        pathsDict = self.pathsDict
        
        #print(f'pathsDict = {pathsDict}\n')
        
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        
        #print(f'srcXnatParams = {srcXnatParams}\n')
        
        projID = srcXnatParams['projID']
        subjLab = srcXnatParams['subjLab']
        srcExpLab = srcXnatParams['expLab']
        srcScanID = srcXnatParams['scanID']
        srcRoicolMod = srcXnatParams['roicolMod']
        srcRoicolName = srcXnatParams['roicolName']
        trgExpLab = trgXnatParams['expLab']
        trgScanID = trgXnatParams['scanID']
        trgRoicolMod = trgXnatParams['roicolMod']
        trgRoicolName = trgXnatParams['roicolName']
        
        srcDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['scans'][srcScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        trgDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][trgExpLab]['scans'][trgScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        srcRoicolFname, srcAsrID\
            = get_scan_asr_fname_and_id(pathsDict, projID, subjLab, srcExpLab,
                                        srcRoicolMod, srcRoicolName)
        
        srcRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['assessors'][srcAsrID]['resources']\
                [srcRoicolMod]['files'][srcRoicolFname]['roicolFpath']
        
        if trgRoicolName == None:
            trgRoicolFpath = None
        else:
            trgRoicolFname, trgAsrID\
                = get_scan_asr_fname_and_id(pathsDict, projID, subjLab, 
                                            trgExpLab, trgRoicolMod, 
                                            trgRoicolName)
        
            trgRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
                ['experiments'][trgExpLab]['assessors'][trgAsrID]['resources']\
                    [trgRoicolMod]['files'][trgRoicolFname]['roicolFpath']
        
        # Add filepaths and directories to srcXnatParams and trgXnatParams:
        srcXnatParams['asrID'] = srcAsrID
        srcXnatParams['roicolFname'] = srcRoicolFname
        srcXnatParams['roicolFpath'] = srcRoicolFpath
        srcXnatParams['dicomDir'] = srcDcmDir
        
        self.srcXnatParams = srcXnatParams
        
        trgXnatParams['dicomDir'] = trgDcmDir
        if trgRoicolName != None:
            trgXnatParams['asrID'] = trgAsrID
            trgXnatParams['roicolFname'] = trgRoicolFname
            trgXnatParams['roicolFpath'] = trgRoicolFpath
        else:
            trgXnatParams['asrID'] = None
            trgXnatParams['roicolFname'] = None
            trgXnatParams['roicolFpath'] = None
            
        self.trgXnatParams = trgXnatParams
    
    def download_210818(self):
        # TODO update docstrings
        """
        Download scans and ROI Collections, and create dictionary containing
        filepaths (pathsDict).
        
        Parameters
        ----------
        xnatSession : Requests Object
            Requests session of XNAT connection.
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        p2c = self.cfgDict['p2c']
        
        # Download data and create pathsDict:
        self.get_pathsDict()
        #self.pathsDict = self.get_pathsDict()
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(config=self.cfgDict)
        
        # Update cfgDict:
        self.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        self.cfgDict['useCaseToApply'] = useCaseToApply
        
        if p2c:
            print(f"cfgDict['runID'] = {self.cfgDict['runID']}")
            print(f'useCaseThatApplies = {useCaseThatApplies}')
            print(f'useCaseToApply = {useCaseToApply}')
            #print(f"cfgDict = {self.cfgDict}")
    

""" From io_tools.exports.py """

def export_newRoicol_MOVED(newRoicol, srcRoicolFpath, exportDir, fname=''):
    """
    16/09/21: Moved to dicom_tools.create_roicol.py
    
    Export ROI Collection (RTS/SEG) to disk.  
    
    Parameters
    ----------
    newRoicol : Pydicom object
        Target ROI Collection (RTS/SEG) to be exported.
    srcRoicolFpath : str
        Full path of the Source RTS/SEG file (used to generate the filename of
        the new RTS/SEG file).
    exportDir : str
        Directory where the new RTS/SEG is to be exported.
    fname : str
        File name to assign.     
    
    Returns
    -------
    newRoicolFpath : str
        Full path of the exported new Target RTS/SEG file.
    """
    
    # Get the filename of the original RTS/SEG file:
    #SrcRoiFname = os.path.split(SrcRoiFpath)[1]
    
    if not os.path.isdir(exportDir):
        #os.mkdir(ExportDir)
        Path(exportDir).mkdir(parents=True)
    
    currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #FnamePrefix = DateTime + '_' + NamePrefix + '_from_'
    #FnamePrefix = DateTime + '_' + TxtToAddToFname.replace(' ', '_')
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = FnamePrefix + SrcRoiFname
    #TrgRoiFname = FnamePrefix
    if fname == '':
        newRoicolFname = currentDateTime + '.dcm'
    else:
        #TrgRoiFname = Fname.replace(' ', '_') + '.dcm'
        #newRoicolFname = f'{fname}_{currentDateTime}.dcm'
        newRoicolFname = f'{fname}.dcm'
    
    newRoicolFpath = os.path.join(exportDir, newRoicolFname)
    
    newRoicol.save_as(newRoicolFpath)
        
    print(f'New Target {newRoicol.Modality} exported to:\n {newRoicolFpath}\n')
    
    return newRoicolFpath


""" From io_tools.fetch_params_and_data.py """

    def get_main_params_OBSOLETE(self, runID):
        """
        Fetches main configuration parameters from main_params.json.
        
        The parameters are returned as a dictionary.
        
        Parameters
        ----------
        runID : str
            The ID that determines the main configuration parameters to use for
            the run (i.e. the key in the base level dictionary contained in
            main_params.json).
        
        Returns
        -------
        dict
            A dictionary containing XNAT configuration parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, 'main_params.json')
        
        # Load all configuration settings:
        allCfgs = import_dict_from_json(filepath=fpath)
        #print(allCfgs.keys())
        
        try:
            mainCfg = allCfgs[runID]
        except KeyError:
            msg = f"No key '{runID}' in main_params.json"
            raise Exception(msg)
        
        # url and/or username and/or password may be blank (e.g. for security
        # reasons). If so prompt user for input and update xnatCfg:
        if mainCfg['url'] == '':
            url = getpass("Enter XNAT url: ")
            mainCfg['url'] = url
        
        if mainCfg['username'] == '':
            username = getpass("Enter user name: ")
            mainCfg['username'] = username
        
        if mainCfg['password'] == '':
            password = getpass("Enter password: ")
            mainCfg['password'] = password
        
        self.mainParams = mainCfg
    
    def get_gen_params_OBSOLETE(self):
        """
        Fetches general parameters from general_params.json.
        
        The parameters are returned as a dictionary.
        
        Returns
        -------
        dict
            A dictionary containing general parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')
        
        fpath = os.path.join(cfgDir, 'general_params.json')
        
        self.genParams = import_dict_from_json(filepath=fpath)
        #return import_dict_from_json(filepath=fpath)
    
    def get_xnat_cfg_OBSOLETE(self):
        """
        OBSOLETE (09/08/21)
        
        Fetches XNAT configuration parameters from xnat_config.json.
        
        The parameters are returned as a dictionary.
        
        Returns
        -------
        dict
            A dictionary containing XNAT configuration parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, 'xnat_config.json')
        
        self.xnatCfg = import_dict_from_json(filepath=fpath)
        #return import_dict_from_json(filepath=fpath)
        
        xnatCfg = import_dict_from_json(filepath=fpath)
        
        # url and/or username and/or password may be blank (e.g. for security
        # reasons). If so prompt user for input and update xnatCfg:
        if xnatCfg['url'] == '':
            url = getpass("Enter XNAT url: ")
            xnatCfg['url'] = url
        
        if xnatCfg['username'] == '':
            username = getpass("Enter user name: ")
            xnatCfg['username'] = username
        
        if xnatCfg['password'] == '':
            password = getpass("Enter password: ")
            xnatCfg['password'] = password
        
        self.xnatCfg = xnatCfg
        
    def get_xnat_params_OBSOLETE(self, srcORtrg):
        """
        OBSOLETE (09/08/21)
        
        Fetches XNAT parameters for Source (from src_xnat_params.json) or
        for Target (from trg_xnat_params.json).
        
        The parameters are returned as a dictionary.
        
        Parameters
        ----------
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        dict
            A dictionary containing XNAT parameters for Source or Target.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, f'{srcORtrg}_xnat_params.json')
        
        if srcORtrg == 'src':
            self.srcXnatParams = import_dict_from_json(filepath=fpath)
        else:
            self.trgXnatParams = import_dict_from_json(filepath=fpath)
        #return import_dict_from_json(filepath=fpath)
    
    def get_res_params_OBSOLETE(self):
        """
        OBSOLETE (09/08/21)
        
        Fetches resampling-/registration-/transforming-related parameters from 
        reg_res_tx_params.json.
        
        The parameters are returned as a dictionary.
        
        Returns
        -------
        dict
            A dictionary containing resampling, registration and transformation
            related parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, 'res_params.json')
        
        self.resParams = import_dict_from_json(filepath=fpath)
        
    def get_pathsDict_Pre_210809(self, xnatSession):
        """
        Downloads scans and ROI Collections and returns a dictionary containing
        filepaths.
        
        Parameters
        ----------
        params : Params Object
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        xnatCfg = self.xnatCfg
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        genParams = self.genParams
        
        pathsDict = download_scan(xnatCfg=xnatCfg, 
                                  xnatParams=srcXnatParams, 
                                  genParams=genParams,
                                  xnatSession=xnatSession)
    
        pathsDict = download_scan(xnatCfg=xnatCfg, 
                                  xnatParams=trgXnatParams, 
                                  genParams=genParams, 
                                  xnatSession=xnatSession,
                                  pathsDict=pathsDict)
    
        pathsDict = download_im_asr(xnatCfg=xnatCfg, 
                                    xnatParams=srcXnatParams,
                                    genParams=genParams, 
                                    xnatSession=xnatSession,
                                    pathsDict=pathsDict)
        
        trgRoicolName = trgXnatParams['roicolName']
        
        if trgRoicolName != None:
            pathsDict = download_im_asr(xnatCfg=xnatCfg, 
                                        xnatParams=trgXnatParams,
                                        genParams=genParams, 
                                        xnatSession=xnatSession,
                                        pathsDict=pathsDict)
        
        #return pathsDict
        self.pathsDict = pathsDict
    
    def add_dirs_and_fpaths_210810(self):
        """
        Add directories and filepaths to params.srcXnatParams and 
        params.trgXnatParams from pathsDict.

        Returns
        -------
        None.

        """
        pathsDict = self.pathsDict
        
        #print(f'pathsDict = {pathsDict}\n')
        
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        
        #print(f'srcXnatParams = {srcXnatParams}\n')
        
        projID = srcXnatParams['projID']
        subjLab = srcXnatParams['subjLab']
        srcExpLab = srcXnatParams['expLab']
        srcScanID = srcXnatParams['scanID']
        srcRoicolMod = srcXnatParams['roicolMod']
        srcRoicolName = srcXnatParams['roicolName']
        trgExpLab = trgXnatParams['expLab']
        trgScanID = trgXnatParams['scanID']
        trgRoicolMod = trgXnatParams['roicolMod']
        trgRoicolName = trgXnatParams['roicolName']
        
        srcDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['scans'][srcScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        trgDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][trgExpLab]['scans'][trgScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        srcRoicolFname, srcAsrID\
            = get_scan_asr_fname_and_id(pathsDict, projID, subjLab, srcExpLab,
                                        srcRoicolMod, srcRoicolName)
        
        srcRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['assessors'][srcAsrID]['resources']\
                [srcRoicolMod]['files'][srcRoicolFname]['roicolFpath']
        
        if trgRoicolName == None:
            trgRoicolFpath = None
        else:
            trgRoicolFname, trgAsrID\
                = get_scan_asr_fname_and_id(pathsDict, projID, subjLab, 
                                            trgExpLab, trgRoicolMod, 
                                            trgRoicolName)
        
            trgRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
                ['experiments'][trgExpLab]['assessors'][trgAsrID]['resources']\
                    [trgRoicolMod]['files'][trgRoicolFname]['roicolFpath']
        
        # Add filepaths and directories to srcXnatParams and trgXnatParams:
        srcXnatParams['asrID'] = srcAsrID
        srcXnatParams['roicolFname'] = srcRoicolFname
        srcXnatParams['roicolFpath'] = srcRoicolFpath
        srcXnatParams['dicomDir'] = srcDcmDir
        
        self.srcXnatParams = srcXnatParams
        
        trgXnatParams['dicomDir'] = trgDcmDir
        if trgRoicolName != None:
            trgXnatParams['asrID'] = trgAsrID
            trgXnatParams['roicolFname'] = trgRoicolFname
            trgXnatParams['roicolFpath'] = trgRoicolFpath
        else:
            trgXnatParams['asrID'] = None
            trgXnatParams['roicolFname'] = None
            trgXnatParams['roicolFpath'] = None
            
        self.trgXnatParams = trgXnatParams
        

""" From io_tools.import_data.py """

    def __init__210817(self, params, srcORtrg):
        
        cfgDict = params.cfgDict
        
        p2c = cfgDict['p2c']
        
        if srcORtrg == 'src':
            dicomDir = cfgDict['srcDicomDir']
            roicolFpath = cfgDict['srcRoicolFpath']
            # Store the slice number in self for latter use:
            self.slcNum = cfgDict['srcSlcNum']
        else:
            dicomDir = cfgDict['trgDicomDir']
            roicolFpath = cfgDict['trgRoicolFpath']
            # Store the slice number in self for latter use:
            self.slcNum = cfgDict['trgSlcNum']
        
        # TODO move DRO somewhere else since DroImporter requires params for
        # Source and Target (06/08/21)
        
        # Import an appropriate DRO:
        #self.droData = DroImporter(params)
        #print(f'type(self.droData.dro) = {type(self.droData.dro)}\n')
        
        # Get DICOM UIDs:
        self.get_dicom_uids(dicomDir)
        
        # Import the ROI Collection:
        self.import_roicol(roicolFpath)
        print(f'type(self.roicol) = {type(self.roicol)}\n')
        
        # Get ROI Names / segment labels (required elsewhere):
        """ 
        Note: roiNames may refer to either ROI Names (for RTS data) or SEG 
        labels (for SEG data). It does not refer to the name of the ROI 
        Collection itself.
        """
        self.roiNames = get_roicol_labels(self.roicol)
        
        # TODO Move this check somewhere else since it requires the Target 
        # roiNames:
        """
        # Check if combination of input parameters are valid:
        valid, errMsg = are_inputs_valid(params, trgRoiNames)
        
        if not valid:
            # Raise exception:
            raise Exception(errMsg)
        """
        
        # Import a DICOM series as a list of Pydicom Objects:
        self.import_dicoms(dicomDir, p2c)
        
        # Import a DICOM series as a SimpleITK Image:
        self.import_dicom_image(dicomDir)
        
        # Get the RTS/SEG data of interest:
        self.get_seg_metadata(cfgDict, srcORtrg)
    
    def __init__210805(self, params):
        
        #srcDcmDir = params.srcXnatParams['dicomDir']
        
        #trgDcmDir = params.trgXnatParams['dicomDir']
        
        #srcRoicolFpath = params.srcXnatParams['roicolFpath']
        
        # Import an appropriate DRO:
        self.droData = DroImporter(params)
        print(f'type(self.droData.dro) = {type(self.droData.dro)}\n')
        
        # Import the Source and Target (if applicable) ROI Collections:
        #self.RoicollectionImporter(params)
        self.import_roicols(self, params)
        print(f'type(self.srcRoicol) = {type(self.srcRoicol)}\n')
        
        # Get ROI Names / segment labels (required elsewhere):
        """ 
        Note: trgRoiNames may refer to either ROI Names or SEG labels, but
        does not relate to the name of the ROI Collection itself.
        """
        srcRoiNames = get_roicol_labels(self.srcRoicol)
        if self.trgRoicol:
            trgRoiNames = get_roicol_labels(self.trgRoicol)
        else:
            trgRoiNames = None
        
        self.srcRoiNames = srcRoiNames
        self.trgRoiNames = trgRoiNames
        
        # Check if combination of input parameters are valid:
        valid, errMsg = are_inputs_valid(params, trgRoiNames)
        
        if not valid:
            # Raise exception:
            raise Exception(errMsg)
        
        # Get DICOM UIDs for Source and Target:
        #self.get_dicom_uids(params, 'src')
        #self.get_dicom_uids(params, 'trg')
        self.get_dicom_uids(params)
        
        # Get the RTS/SEG data of interest:
        self.get_seg_metadata(params, 'src')
        if self.trgRoicol:
            self.get_seg_metadata(params, 'trg')
        
        # Import the Source and Target DICOM series as a list of Pydicom 
        # Objects:
        self.get_dicoms(self, params)
        
        # Import the Source and Target DICOM series as a SimpleITK Images:
        self.get_images(self, params)
    
    def add_timestamp_OLD(self, timingMsg):
        # TODO update docstrings
        """
        Add timestamp and timing message. Replace the keyword 'dTime' in
        timingMsg with the actual dTime calculated below.
        """
        self.timings.append(time.time())
        dTime = round(self.timings[-1] - self.timings[-2], 1)
        timingMsg.replace('dTime', f'{dTime}')
    #    self.timingMsgs.append(timingMsg)
        print(f'*{timingMsg}')
    
    def get_dicom_uids_210805(self, params):
        """
        Get StudyUID, SeriesUID, FrameOfReferenceUID and SOPInstanceUIDs for
        both Source and Target DICOM series.
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
        
        Returns
        -------
        self.srcStudyUID or self.trgStudyUID : str
            DICOM Study UID.
        self.srcSeriesUID or self.trgSeriesUID : str
            DICOM Series UID.
        self.srcFORUID or self.trgFORUID : str
            DICOM Frame of reference UID.
        self.srcSOPUID or self.trgSOPUID : str
            DICOM SOP UID.
        """
        srcDcmDir = params.srcXnatParams['dicomDir']
        trgDcmDir = params.trgXnatParams['dicomDir']
            
        self.srcStudyUID, self.srcSeriesUID, self.srcFORUID, self.srcSOPUIDs\
            = get_dcm_uids(srcDcmDir)
        
        self.trgStudyUID, self.trgSeriesUID, self.trgFORUID, self.trgSOPUIDs\
            = get_dcm_uids(trgDcmDir)
    
    def get_dicom_uids_210721(self, params, srcORtrg):
        """
        Get StudyUID, SeriesUID, FrameOfReferenceUID and SOPInstanceUIDs.
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        None.
        
        """
        
        if srcORtrg == 'src':
            dicomDir = params.srcXnatParams['dicomDir']
            
            studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(dicomDir)
            
            self.srcStudyUID = studyUID
            self.srcSeriesUID = seriesUID
            self.srcFORUID = FORUID
            self.srcSOPUIDs = SOPUIDs
        else:
            dicomDir = params.trgXnatParams['dicomDir']
            
            studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(dicomDir)
            
            self.trgStudyUID = studyUID
            self.trgSeriesUID = seriesUID
            self.trgFORUID = FORUID
            self.trgSOPUIDs = SOPUIDs
    
    def import_roicols_210805(self, params):
        """
        Imports ROI Collections as Pydicom Objects.
        
        Since SEG data may not exist (e.g. for Target) this method may assign 
        NoneType to the attribute (e.g. self.trgRoiCol).
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
            
        Returns
        -------
        self.srcRoicol or self.trgRoicol : Pydicom Object or None
            Pydicom representation of the ROI Collection (RTS/SEG) if not None.
        """
        
        srcFpath = params.srcXnatParams['roicolFpath']
        trgFpath = params.trgXnatParams['roicolFpath']
        
        if srcFpath != None:
            self.srcRoicol = dcmread(srcFpath)
        else:
            self.srcRoicol = None
        
        if trgFpath != None:
            self.trgRoicol = dcmread(trgFpath)
        else:
            self.trgRoicol = None
    
    def import_roicol_210721(self, params, srcORtrg):
        """
        Imports an ROI Collection as a Pydicom Object.
        
        Since SEG data may not be provided for Target, this function will only
        get SEG data for Source or Target (via the parameter srcORtrg).
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
        srcORtrg : str
            'src' or 'trg' for Source or Target.
            
        Returns
        -------
        None.
        """
        
        if srcORtrg == 'src':
            filepath = params.srcXnatParams['roicolFpath']
            
            if filepath != None:
                self.srcRoicol = dcmread(filepath)
            else:
                self.srcRoicol = None
        else:
            # i.e.: srcORtrg == 'trg'
            filepath = params.trgXnatParams['roicolFpath']
            
            if filepath != None:
                self.trgRoicol = dcmread(filepath)
            else:
                self.trgRoicol = None
    
    def get_seg_metadata_210809(self, xnatParams, p2c=False):
        """
        Get metadata for a SEG. 
        
        Since SEG data may not be provided for Target, this function will only
        get SEG data for Source or Target (via the parameter srcORtrg).
        
        Parameters
        ----------
        xnatParams : dict
            Contains various XNAT-associated parameters that relate to the Source 
        or Target DICOM series. An instance variable of a Params Object.
        p2c : bool (optional)
            If True some results will be printed to the console. False by 
            default.
        
        Returns
        -------
        self.divs : list of ints or None
            Dimensional Index Values
        self.allF2Sinds : list of ints or None
            Frame-to-slice indices (i.e. slice indices corresponding to each
            frame for the ROI Collection).
        self.allF2SindsBySeg : list of a list of ints or None
            A list of frame-to-slice indices for each ROI/segment in the ROI
            Collection.
        self.roiNums : list of ints or None
            List of ints for each ROI/segment in the ROI Collection (e.g.
            [0, 1] => 2 ROIs/segments).
        self.pixArrBySeg : list of Numpy data array Objects or None
            List of pixel arrays - one for each ROI/segment.
        self.f2sIndsBySeg : list of a list of ints or None
            List for each ROI/segment of a list of frame-to-slice indices.  
        """
        
        # Get metadata:
        if self.roicol == None:
            self.divs = None
            self.allF2Sinds = None
            self.allF2SindsBySeg = None
            self.roiNums = None
            self.pixArrBySeg = None
            self.f2sIndsBySeg = None
        else:
            # Get list of DimensionIndexValues:
            print(self.roicol)
            self.divs = get_divs(self.roicol)
            
            # Get list of PerFrameFunctionalGroupsSequence-to-slice indices
            # and the same grouped by segment:
            self.allF2Sinds = get_p2sInds(
                seg=self.roicol, 
                sopuids=self.sopuids)

            self.allF2SindsBySeg = group_list_by_seg(
                listToGroup=self.allF2Sinds, 
                divs=self.divs)
            
            # Get the segment number(s) that match the segment label of
            # interest (roiName):
            """ 
            e.g. 0 if the segment of interest is the only segment, or if it
            is the first of several.
            
            Using the attribute srcRoiNums rather than srcSegNums for 
            interchangability for ROI data.
            """
            self.roiNums = get_roicol_nums(
                seg=self.roicol, 
                searchStr=xnatParams['roiName'])
            
            # Get the pixel array by segment and frame-to-slice indices by
            # segment for the data of interest (to be copied):
            self.pixArrBySeg, self.f2sIndsBySeg\
                = get_seg_data_of_interest(
                    seg=self.roicol, 
                    allF2SindsBySeg=self.allF2SindsBySeg, 
                    segNums=self.roiNums, 
                    segLabs=self.roiNames,
                    segLab=xnatParams['roiName'],
                    slcNum=xnatParams['slcNum'],
                    p2c=p2c)
            
            # Raise exception if self.srcF2SindsBySeg is empty (i.e. if 
            # there were no segmentations that matched the input parameters)
            raise_error_if_no_seg_data_of_interest(
                f2sIndsBySeg=self.f2sIndsBySeg, 
                segLabs=self.roiNames, 
                slcNum=xnatParams['slcNum'])
    
    def get_seg_metadata_210805(self, params, srcORtrg):
        """
        Get metadata for a SEG. 
        
        Since SEG data may not be provided for Target, this function will only
        get SEG data for Source or Target (via the parameter srcORtrg).
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        self.srcDIVs or self.trgDIVs : list of ints or None
            Dimensional Index Values
        self.srcAllF2Sinds or self.trgAllF2Sinds : list of ints or None
            Frame-to-slice indices (i.e. slice indices corresponding to each
            frame for the ROI Collection).
        self.srcAllF2SindsBySeg or self.trgAllF2SindsBySeg : list of a list of
            ints or None
            A list of frame-to-slice indices for each ROI/segment in the ROI
            Collection.
        self.srcRoiNums or self.trgRoiNums : list of ints or None
            List of ints for each ROI/segment in the ROI Collection (e.g.
            [0, 1] => 2 ROIs/segments).
        self.srcPixArrBySeg or self.trgPixArrBySeg : list of Numpy data array 
            Objects or None
            List of pixel arrays - one for each ROI/segment.
        self.srcF2SindsBySeg or self.trgF2SindsBySeg : list of a list of ints
            or None
            List for each ROI/segment of a list of frame-to-slice indices.  
        """
        
        # Get required inputs specific to source or target:
        if srcORtrg == 'src':
            roiName = params.srcXnatParams['roiName']
            slcNum = params.srcXnatParams['slcNum']
            #seg = self.srcSeg
            roicol = self.srcRoicol
            SOPUIDs = self.srcSOPUIDs
            f2sInds = self.srcF2Sinds
            DIVs = self.srcDIVs
            roiNums = self.srcRoiNums 
            roiNames = self.srcRoiNames
        else:
            # i.e.: srcORtrg == 'trg'
            roiName = params.trgXnatParams['roiName']
            slcNum = params.trgXnatParams['slcNum']
            #seg = self.trgSeg
            roicol = self.trgRoicol
            SOPUIDs = self.trgSOPUIDs
            f2sInds = self.trgF2Sinds
            DIVs = self.trgDIVs
            roiNums = self.trgRoiNums 
            roiNames = self.trgRoiNames
        
        # Get metadata:
        if roicol == None:
            DIVs = None
            allF2Sinds = None
            allF2SindsBySeg = None
            roiNums = None
            pixArrBySeg = None
            f2sIndsBySeg = None
        else:
            # Get list of DimensionIndexValues:
            DIVs = get_divs(roicol)
            
            # Get list of PerFrameFunctionalGroupsSequence-to-slice indices
            # and the same grouped by segment:
            allF2Sinds = get_p2sInds(
                seg=roicol, 
                SOPUIDs=SOPUIDs)

            allF2SindsBySeg = group_list_by_seg(
                listToGroup=f2sInds, 
                DIVs=DIVs)
            
            # Get the segment number(s) that match the segment label of
            # interest (roiName):
            """ 
            e.g. 0 if the segment of interest is the only segment, or if it
            is the first of several.
            
            Using the attribute srcRoiNums rather than srcSegNums for 
            interchangability for ROI data.
            """
            roiNums = get_roicol_nums(seg=roicol, searchStr=roiName)
            
            # Get the pixel array by segment and frame-to-slice indices by
            # segment for the data of interest (to be copied):
            pixArrBySeg, f2sIndsBySeg\
                = get_seg_data_of_interest(
                    seg=roicol, 
                    allF2SindsBySeg=allF2SindsBySeg, 
                    segNums=roiNums, 
                    segLabs=roiNames,
                    segLab=roiName,
                    slcNum=slcNum,
                    p2c=params.genParams['p2c'])
            
            # Raise exception if self.srcF2SindsBySeg is empty (i.e. if 
            # there were no segmentations that matched the input parameters)
            raise_error_if_no_seg_data_of_interest(
                f2sIndsBySeg=f2sIndsBySeg, 
                segLabs=roiNames, 
                slcNum=slcNum)
            
        # Assign metadata to self for src or trg:
        if srcORtrg == 'src':
            self.srcDIVs = DIVs
            self.srcAllF2Sinds = allF2Sinds
            self.srcAllF2SindsBySeg = allF2SindsBySeg
            self.srcRoiNums = roiNums
            self.srcPixArrBySeg = pixArrBySeg
            self.srcF2SindsBySeg = f2sIndsBySeg
        else:
            # i.e.: srcORtrg == 'trg'
            self.trgDIVs = DIVs
            self.trgAllF2Sinds = allF2Sinds
            self.trgAllF2SindsBySeg = allF2SindsBySeg
            self.trgRoiNums = roiNums
            self.trgPixArrBySeg = pixArrBySeg
            self.trgF2SindsBySeg = f2sIndsBySeg
    
    def get_seg_metadata_210721(self, params, srcORtrg):
        """
        Get metadata for a SEG. 
        
        Since SEG data may not be provided for Target, this function will only
        get SEG data for Source or Target (via the parameter srcORtrg).
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        None.
        
        """
        
        if srcORtrg == 'src':
            #seg = self.srcSeg
            
            
            if self.srcSeg == None:
                self.srcDIVs = None
            else:
                # Get list of DimensionIndexValues:
                self.srcDIVs = get_divs(self.srcSeg)
                
                # Get list of PerFrameFunctionalGroupsSequence-to-slice indices
                # and the same grouped by segment:
                self.srcAllF2Sinds = get_p2sInds(
                    seg=self.srcSeg, 
                    SOPUIDs=self.srcSOPUIDs)
    
                self.srcAllF2SindsBySeg = group_list_by_seg(
                    listToGroup=self.srcF2Sinds, 
                    DIVs=self.srcDIVs)
                
                # Get the segment number(s) that match the segment label of
                # interest (roiName):
                """ 
                e.g. 0 if the segment of interest is the only segment, or if it
                is the first of several.
                
                Using the attribute srcRoiNums rather than srcSegNums for 
                interchangability for ROI data.
                """
                self.srcRoiNums = get_roicol_nums(
                    seg=self.srcSeg, 
                    searchStr=params.srcXnatParams['roiName'])
                
                # Get the pixel array by segment and frame-to-slice indices by
                # segment for the data of interest (to be copied):
                self.srcPixArrBySeg, self.srcF2SindsBySeg\
                    = get_seg_data_of_interest(
                        seg=self.srcSeg, 
                        allF2SindsBySeg=self.srcAllF2SindsBySeg, 
                        segNums=self.srcRoiNums, 
                        segLabs=self.srcRoiNames,
                        segLab=params.srcXnatParams['roiName'],
                        slcNum=params.srcXnatParams['slcNum'],
                        p2c=params.genParams['p2c'])
                
                # Raise exception if self.srcF2SindsBySeg is empty (i.e. if 
                # there were no segmentations that matched the input parameters)
                raise_error_if_no_seg_data_of_interest(
                    f2sIndsBySeg=self.srcF2SindsBySeg, 
                    segLabs=self.srcRoiNames, 
                    slcNum=params.srcXnatParams['slcNum'])
                
                """ Instead of checking the proportion of segs in extent 
                within this importing module, do this after importing... """
        
        else:
            #seg = self.trgSeg
            
            """ Repeat the above lines but for trg? """
    
    def import_dicoms_210805(self, params):
        """
        Imports DICOM series as a list of Pydicom Objects.
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
            
        Returns
        -------
        self.srcDcms or self.trgDcms : list of Pydicom Objects
            A list (for each DICOM) of Pydicom representations of the DICOM
            series.
        """
        srcDcmDir = params.srcXnatParams['dicomDir']
        trgDcmDir = params.trgXnatParams['dicomDir']
        
        p2c = params.genParams['p2c']
        
        self.srcDcms = import_dcms(dicomDir=srcDcmDir, p2c=p2c)
        self.trgDcms = import_dcms(dicomDir=trgDcmDir, p2c=p2c)
    
    def import_images_210805(self, params):
        """
        Imports DICOM series as a SimpleITK Image.
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
            
        Returns
        -------
        self.srcIm or self.trgIm : SimpleITK Image
            SimpleITK Image representation of the DICOM series.
        """
        srcDcmDir = params.srcXnatParams['dicomDir']
        trgDcmDir = params.trgXnatParams['dicomDir']
        
        self.srcIm = import_dicoms_as_im(dicomDir=srcDcmDir)
        self.trgIm = import_dicoms_as_im(dicomDir=trgDcmDir)
    

""" From io_tools.inputs_checker.py """

def are_inputs_valid_210810(params, trgRoiNames):
    """
    Check if the combination of input parameters are valid.

    Parameters
    ----------
    params : Params Object
        Object containing various parameters.
    trgRoiNames : list of strs or None
        List of ROI Names / segment labels of the Target ROI Collection 
        (RTS/SEG), or None if there is no Target ROI Collection.
    
    Returns
    -------
    valid : bool
        True if the combinations of inputs are valid.
    msg : str
        The error message. msg should be '' if valid = True, and non-empty
        otherwise.
    """
    
    mod = params.srcXnatParams['roicolMod']
    srcRoiName = params.srcXnatParams['roiName']
    srcSlcNum = params.srcXnatParams['slcNum']
    
    trgRoiName = params.trgXnatParams['roiName']
    trgSlcNum = params.trgXnatParams['slcNum']
    
    """ Useful variables:
    RorS = 'ROI' or 'segment',  
    CorS = 'contours' or 'segmentations', and
    NorL = 'name' or 'label'
    """
    if mod == 'RTSTRUCT':
        RorS = 'ROI'
        CorS = 'contours'
        NorL = 'name'
    else:
        RorS = 'segment'
        CorS = 'segmentations'
        NorL = 'label'
    
    valid = False
    
    msg = ''
    
    if not srcRoiName:
        if srcSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All {CorS} in all {RorS}s "\
                + f"in the Source {mod} are to be copied. \nBut srcSlcNum (= "\
                + f"{srcSlcNum}) is specified (i.e. != None). "
        
        if trgRoiName:
            msg += f"srcRoiName = {srcRoiName} => All {CorS} in all {RorS}s "\
                + f"in the Source {mod} are to be copied.\nBut trgRoiName (= "\
                + f"{trgRoiName}) is specified(i.e. != None). "
        
        if trgSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All {CorS} in all {RorS}s "\
                + f"in the Source {mod} are to be copied. \nBut trgSlcNum (= "\
                + f"{trgSlcNum}) is specified (i.e. != None). "
    else:
        if not srcSlcNum and trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => All {CorS} in the {RorS} with"\
                + f" {NorL} '{srcRoiName}' are to be copied.\nBut trgSlcNum "\
                + f"(= {trgSlcNum}) is specified. "
        
        if srcSlcNum and not trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => The {CorS}(s) for the "\
                + f"specified DICOM slice in the {RorS} with {NorL} "\
                + f"'{srcRoiName}' are to be copied.\nBut trgSlcNum (= "\
                + f"{trgSlcNum}) is not specified. "
        
        if trgRoiNames != None:
            T = len(trgRoiNames)
            
            if T > 1 and not trgRoiName:
                msg += f"srcRoiName (= {srcRoiName}) has been specified but "\
                    + f"trgRoiName (= {trgRoiName}) has not, and there are "\
                    + f"{T} {RorS}s in the Target ROI Collection. "
    
    if msg:
        msg += "This is not a valid combination of inputs."
        #raise Exception(msg)
    else:
        valid = True
    
    return valid, msg

def are_inputs_valid_OLD(srcRoi, srcRoiName, srcSlcNum, trgRoi, trgRoiName, 
                     trgSlcNum, p2c=False):
    """
    Check if the combination of inputs are valid.

    Parameters
    ----------
    srcRoi : Pydicom object
        The Source RTS/SEG object.
    srcRoiName : string
        The Source ROIName or SegmentLabel.
    srcSlcNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If srcSlcNum = None, a relationship-preserving copy will be made.
    trgRoi : Pydicom object or None
        The Target RTS/SEG object if applicable; None otherwise.
    trgRoiName : string or None
        The Target ROIName or SegmentLabel if applicable; None otherwise.
    trgSlcNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If trgSlcNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    p2c : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    IsValid : boolean
        True if the combinations of inputs are valid.
    """
    
    msg = ''
              
    #if srcRoiName == None:
    if not srcRoiName:
        #if srcSlcNum != None:
        if srcSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut srcSlcNum (= {srcSlcNum}) is "\
                   + "specified (i.e. != None). "
        
        #if trgRoiName != None:
        if trgRoiName:
            msg += f"srcRoiName = {srcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut trgRoiName (= {trgRoiName}) is specified "\
                   + "(i.e. != None). "
        
        #if trgSlcNum != None:
        if trgSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut trgSlcNum (= {trgSlcNum}) is "\
                   + "specified (i.e. != None). "
        
    else:
        #if srcSlcNum == None and trgSlcNum != None:
        if not srcSlcNum and trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => All contours/segmentations"\
                   + f"in the ROI/segment with name/label '{srcRoiName}' are "\
                   + f"to be copied. \nBut trgSlcNum (= {trgSlcNum}) is "\
                   + "specified. "
        
        #if srcSlcNum != None and trgSlcNum == None:
        if srcSlcNum and not trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => The contour(s)/"\
                   + "segmentation(s) for the specified DICOM slice in the "\
                   + f"ROI/segment with name/label '{srcRoiName}' are to be "\
                   + f"copied. \nBut trgSlcNum (= {trgSlcNum}) is not "\
                   + "specified. "
        
        if trgRoi != None:
            """ Get the names/labels of all ROIs/segments in trgRoi: """
            trgRoiNames = get_roicol_labels(trgRoi)
            
            T = len(trgRoiNames)
            
            #if T > 1 and trgRoiName == None:
            if T > 1 and not trgRoiName:
                msg += f"srcRoiName (= {srcRoiName}) has been specified but "\
                       + f"trgRoiName (= {trgRoiName}) has not been specified"\
                       + f" and there are {T} ROIs/segments in trgRoi. "
        
    if msg:
        msg += "This is not a valid combination of inputs."
            
        raise Exception(msg)
    
    return True

def which_use_case_210810(srcXnatParams, trgXnatParams, resParams, genParams):
    """
    Determine which of 5 Use Cases applies.
    
    Parameters
    ----------
    srcXnatParams : dict
        Dictionary of XNAT parameters for Source.
    trgXnatParams : dict
        Dictionary of XNAT parameters for Target.
    resParams : dict
        Dictionary of general parameters.
    genParams : dict
        Dictionary containing resampling, registration and transformation
        related parameters.
    
    Returns
    -------
    useCaseThatApplies : str
        String that denotes the use case that applies.
    useCaseToApply : str
        String that denotes the use case that will be applied.
    
    Note
    ----
    useCaseToApply will be the same as useCaseThatApplies unless forceReg is 
    True and useCaseThatApplies is '3a', '3b', '4a' or '4b'.
    """
    
    srcDcmDir = srcXnatParams['dicomDir']
    
    trgDcmDir = trgXnatParams['dicomDir']
    trgSlcNum = trgXnatParams['slcNum']
    
    forceReg = resParams['forceReg']
    
    p2c = genParams['p2c']
    
    # Get the DICOM UIDs:
    srcStudyUID, srcSeriesUID, srcFORUID,\
    srcSOPUIDs = get_dcm_uids(dicomDir=srcDcmDir)
    
    trgStudyUID, trgSeriesUID, trgFORUID,\
    trgSOPUIDs = get_dcm_uids(dicomDir=trgDcmDir)
        
    if trgSlcNum:
        if trgSlcNum > len(trgSOPUIDs) - 1:
            msg = f'The input argument "trgSlcNum" = {trgSlcNum} exceeds '\
                  + f'the number of Target DICOMs ({len(trgSOPUIDs)}).'
            raise Exception(msg)
    
    # Get the Image Attributes for Source and Target:
    srcSize, srcSpacing, srcST, srcIPPs, srcDirs,\
    warnings = get_im_attrs(dicomDir=srcDcmDir, package='sitk')
    
    trgSize, trgSpacing, trgST, trgIPPs, trgDirs,\
    warnings = get_im_attrs(dicomDir=trgDcmDir, package='sitk')
    
    
    # Are the Source and Target part of the same study, series and do they 
    # have the same SOP UIDs?
    if srcStudyUID == trgStudyUID:
        sameStudy = True
    else:
        sameStudy = False
        
    if srcSeriesUID == trgSeriesUID:
        sameSeries = True
    else:
        sameSeries = False
        
    if srcFORUID == trgFORUID:
        sameFOR = True
    else:
        sameFOR = False
        
    if srcSOPUIDs == trgSOPUIDs:
        sameSOPs = True
    else:
        sameSOPs = False
    
    # Are the Source and Target spacings, directions and origins equal to  
    # within 1e-06?
    sameSpacings = are_items_equal_to_within_eps(srcSpacing, trgSpacing)
    sameDirs = are_items_equal_to_within_eps(srcDirs, trgDirs)
    sameOrigins = are_items_equal_to_within_eps(srcIPPs[0], trgIPPs[0])
    sameST = are_items_equal_to_within_eps([srcST], [trgST])
    
    if srcSize == trgSize:
        sameSize = True
    else:
        sameSize = False
    
    # Check which case follows:
    if sameFOR:
        if sameDirs:
            """ 
            Comment 11/02:
                
            It's not entirely necessary to resample the Source image to the 
            Target domain just because their origins differ. It depends on
            whether the slice in Source coincide with slices in Target.
            
            Even if the slices don't coincide, an approximation can be made as
            follows.  When dealing contours, the in-plane coordinates don't 
            need to be changed, and the destination slice can be determined by 
            matching on the out-of-plane (z) coordinates as best as possible. 
            When dealing with segmentations, the destination slice can be 
            determined by matching on the z components of the IPPs as best as
            possible.  Once the destination slice number is known, the 
            necessary shift in pixels for the in-plane components can be 
            applied accordingly. 
            
            Perhaps the approximation can be justified if the difference 
            between the planes in Source and Target are close to an integer.
            If the difference is large, the Source image should be resampled to
            the Target image.  For now I will treat UseCase 2 data sets as 
            UseCase 3 if their origins are not the same. Further complexity
            (taking into account the error in the above approximation) can be
            considered for future work. 
            """
            if sameSpacings and sameOrigins:
            #if sameSpacings:
                if sameSeries:
                    useCaseThatApplies = '1' # Direct copy
                else:
                    if isinstance(trgSlcNum, int):
                        useCaseThatApplies = '2a' # Direct copy
                    else:
                        useCaseThatApplies = '2b' # Relationship-preserving copy
            else:
                """ 
                Comment 11/02:
                    
                Arguably resampling is only required if the spacings along any
                direction are different when dealing with segmentations.  When
                dealing with contours, resampling is probably unnecessary if
                the spacings are different only in the out-of-plane (z)
                direction. 
                """
                if isinstance(trgSlcNum, int):
                    useCaseThatApplies = '3a' # Direct copy
                else:
                    useCaseThatApplies = '3b' # Relationship-preserving copy
        else:
            if isinstance(trgSlcNum, int):
                useCaseThatApplies = '4a' # Direct copy
            else:
                useCaseThatApplies = '4b' # Relationship-preserving copy
    else:
        if isinstance(trgSlcNum, int):
            useCaseThatApplies = '5a' # Direct copy
        else:
            useCaseThatApplies = '5b' # Relationship-preserving copy
    
    if p2c:
        # Print comparisons to the console:
        srcIm = import_im(srcDcmDir)
        trgIm = import_im(trgDcmDir)
        
        # Get the image extents
        srcImExtent = get_im_extent(srcIm)
        trgImExtent = get_im_extent(trgIm)
        
        # Are the image extents equal to within 1e-06?
        sameExtents = are_items_equal_to_within_eps(srcImExtent, trgImExtent)
    
        #CompareImageAttributes(srcDcmDir, trgDcmDir)
        
        print('\n\n', '-'*120)
        print('Results of running which_use_case():')
        print('\n\n', '-'*120)
        print('The Source and Target are in the/have:')
        
        if sameStudy:
            print('   * same study')
        else:
            print('   * different studies')
        
        if sameFOR:
            print('   * same Frame of Reference')
        else:
            print('   * different Frames of Reference')
            
        if sameSeries:
            print('   * same series')
        else:
            print('   * different series')
            
        if sameSOPs:
            print('   * same DICOMs')
        else:
            print('   * different DICOMs')
        
        if sameSize:
            print(f'   * same image size\n      Source/Target: {srcSize}')
        else:
            print(f'   * different image sizes:\n      Source: {srcSize}',
                  f'\n      Target: {trgSize}')
        
        if sameSpacings:
            print(f'   * same voxel sizes\n      Source/Target: {srcSpacing}')
        else:
            print(f'   * different voxel sizes:\n      Source: {srcSpacing}',
                  f'\n      Target: {trgSpacing}')
        
        if sameST:
            print(f'   * same slice thickness\n      Source/Target: {srcST}')
        else:
            print(f'   * different slice thickness:\n      Source: {srcST}',
                  f'\n      Target: {trgST}')
        
        if sameExtents:
            print(f'   * same image extents\n      Source/Target: {srcImExtent}')
        else:
            print(f'   * different image extents:\n      Source: {srcImExtent}',
                  f'\n      Target: {trgImExtent}')
        
        if sameOrigins:
            print(f'   * same origin\n      Source/Target: {srcIPPs[0]}')
        else:
            print(f'   * different origins:\n      Source: {srcIPPs[0]}',
                  f'\n      Target: {trgIPPs[0]}')
            
        if sameDirs:
            print('   * same patient orientation\n',
                  f'      Source/Target: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                  f'                      {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                  f'                      {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]')
        else:
            print('   * different patient orientations:\n',
                  f'      Source: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                  f'               {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                  f'               {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]\n',
                  f'      Target: [{trgDirs[0]}, {trgDirs[1]}, {trgDirs[2]},\n',
                  f'               {trgDirs[3]}, {trgDirs[4]}, {trgDirs[5]},\n',
                  f'               {trgDirs[6]}, {trgDirs[7]}, {trgDirs[8]}]')
    
    if True:#p2c:
        #PrintTitle(f'Case {UseCase} applies.')
        if forceReg and useCaseThatApplies in ['3a', '3b', '4a', '4b']:
            useCaseToApply = useCaseThatApplies.replace('3', '5').replace('4', '5')
            
            msg = f'\n*UseCase {useCaseThatApplies} applies but since '\
                  + f'forceReg = {forceReg}, UseCase {useCaseToApply} will '\
                  + 'be applied.'
            print(msg)
        else:
            useCaseToApply = useCaseThatApplies
            
            print(f'\n*UseCase {useCaseThatApplies} applies.')
            
        print('-'*120)
        
    return useCaseThatApplies, useCaseToApply

def which_use_case_Pre_210810(srcSlcNum, trgSlcNum, srcDcmDir, trgDcmDir, forceReg, 
                   p2c=False):
    """
    Determine which of 5 Use Cases applies.
    
    Parameters
    ----------
    srcSlcNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If srcSlcNum = None, a relationship-preserving copy will be made.
    trgSlcNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If trgSlcNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    srcDcmDir : string
        Directory containing the Source DICOMs.
    trgDcmDir : string
        Directory containing the Target DICOMs.
    forceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    p2c : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    useCaseThatApplies : string
        String that denotes the Use Case that applies.
    useCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as useCaseThatApplies unless forceRegistration is True and
        useCaseThatApplies is '3a', '3b', '4a' or '4b'.
    """
    
    from dicom_tools.metadata import get_dcm_uids
    from image_tools.attrs_info import get_im_attrs
    from general_tools.general import items_unique_to_within
    
    
    """ Get the DICOM UIDs. """
    srcStudyuid, srcSeriesuid, srcFORuid,\
    srcSOPuids = get_dcm_uids(DicomDir=srcDcmDir)
    
    trgStudyuid, trgSeriesuid, trgFORuid,\
    trgSOPuids = get_dcm_uids(DicomDir=trgDcmDir)
        
    if trgSlcNum:
        if trgSlcNum > len(trgSOPuids) - 1:
            msg = f'The input argument "trgSlcNum" = {trgSlcNum} exceeds '\
                  + f'the number of Target DICOMs ({len(trgSOPuids)}).'
            
            raise Exception(msg)
    
    
    """ Get the Image Attributes for Source and Target using SimpleITK. """
    srcSize, srcSpacing, srcST, srcIPPs, srcDirs,\
    warnings = get_im_attrs(DicomDir=srcDcmDir, Package='sitk')
    
    trgSize, trgSpacing, trgST, trgIPPs, trgDirs,\
    warnings = get_im_attrs(DicomDir=trgDcmDir, Package='sitk')
    
    
    """ Are the Source and Target part of the same study, series and do they 
    have the same SOP UIDs? """
    if srcStudyuid == trgStudyuid:
        sameStudy = True
    else:
        sameStudy = False
        
    if srcSeriesuid == trgSeriesuid:
        SameSeries = True
    else:
        SameSeries = False
        
    if srcFORuid == trgFORuid:
        sameFOR = True
    else:
        sameFOR = False
        
    if srcSOPuids == trgSOPuids:
        sameSOPs = True
    else:
        sameSOPs = False
    
    """ Are the Source and Target spacings, directions and origins equal to  
    within 1e-06? """
    sameSpacings = items_unique_to_within(srcSpacing, trgSpacing)
    sameDirs = items_unique_to_within(srcDirs, trgDirs)
    sameOrigins = items_unique_to_within(srcIPPs[0], trgIPPs[0])
    sameST = items_unique_to_within([srcST], [trgST])
    
    if srcSize == trgSize:
        sameSize = True
    else:
        sameSize = False
    
        
    """ Check which case follows. """
    if sameFOR:
        if sameDirs:
            """ 
            Comment 11/02:
                
            It's not entirely necessary to resample the Source image to the 
            Target domain just because their origins differ. It depends on
            whether the slice in Source coincide with slices in Target.
            
            Even if the slices don't coincide, an approximation can be made as
            follows.  When dealing contours, the in-plane coordinates don't 
            need to be changed, and the destination slice can be determined by 
            matching on the out-of-plane (z) coordinates as best as possible. 
            When dealing with segmentations, the destination slice can be 
            determined by matching on the z components of the IPPs as best as
            possible.  Once the destination slice number is known, the 
            necessary shift in pixels for the in-plane components can be 
            applied accordingly. 
            
            Perhaps the approximation can be justified if the difference 
            between the planes in Source and Target are close to an integer.
            If the difference is large, the Source image should be resampled to
            the Target image.  For now I will treat UseCase 2 data sets as 
            UseCase 3 if their origins are not the same. Further complexity
            (taking into account the error in the above approximation) can be
            considered for future work. """
            if sameSpacings and sameOrigins:
            #if sameSpacings:
                if SameSeries:
                    useCaseThatApplies = '1' # Direct copy
                else:
                    if isinstance(trgSlcNum, int):
                        useCaseThatApplies = '2a' # Direct copy
                    else:
                        useCaseThatApplies = '2b' # Relationship-preserving copy
            else:
                """ Comment 11/02:
                    
                Arguably resampling is only required if the spacings along any
                direction are different when dealing with segmentations.  When
                dealing with contours, resampling is probably unnecessary if
                the spacings are different only in the out-of-plane (z)
                direction. """
                if isinstance(trgSlcNum, int):
                    useCaseThatApplies = '3a' # Direct copy
                else:
                    useCaseThatApplies = '3b' # Relationship-preserving copy
        else:
            if isinstance(trgSlcNum, int):
                useCaseThatApplies = '4a' # Direct copy
            else:
                useCaseThatApplies = '4b' # Relationship-preserving copy
    else:
        if isinstance(trgSlcNum, int):
            useCaseThatApplies = '5a' # Direct copy
        else:
            useCaseThatApplies = '5b' # Relationship-preserving copy
    
        
    """ Print comparisons to the console. """
    if p2c:
        from ImageTools import import_im, GetImageExtent
        
        srcIm = import_im(srcDcmDir)
        trgIm = import_im(trgDcmDir)
        
        """ Get the image extents. """
        srcImExtent = GetImageExtent(srcIm)
        trgImExtent = GetImageExtent(trgIm)
        
        """ Are the image extents equal to within 1e-06? """
        sameExtents = items_unique_to_within(srcImExtent, trgImExtent)
    
        #CompareImageAttributes(srcDcmDir, trgDcmDir)
        
        print('\n\n', '-'*120)
        print('Results of running WhichUseCase():')
        print('\n\n', '-'*120)
        print('The Source and Target are in the/have:')
        
        if sameStudy:
            print('   * same study')
        else:
            print('   * different studies')
        
        if sameFOR:
            print('   * same Frame of Reference')
        else:
            print('   * different Frames of Reference')
            
        if SameSeries:
            print('   * same series')
        else:
            print('   * different series')
            
        if sameSOPs:
            print('   * same DICOMs')
        else:
            print('   * different DICOMs')
        
        if sameSize:
            print(f'   * same image size\n      Source/Target: {srcSize}')
        else:
            print(f'   * different image sizes:\n      Source: {srcSize}',
                  f'\n      Target: {trgSize}')
        
        if sameSpacings:
            print(f'   * same voxel sizes\n      Source/Target: {srcSpacing}')
        else:
            print(f'   * different voxel sizes:\n      Source: {srcSpacing}',
                  f'\n      Target: {trgSpacing}')
        
        if sameST:
            print(f'   * same slice thickness\n      Source/Target: {srcST}')
        else:
            print(f'   * different slice thickness:\n      Source: {srcST}',
                  f'\n      Target: {trgST}')
        
        if sameExtents:
            print(f'   * same image extents\n      Source/Target: {srcImExtent}')
        else:
            print(f'   * different image extents:\n      Source: {srcImExtent}',
                  f'\n      Target: {trgImExtent}')
        
        if sameOrigins:
            print(f'   * same origin\n      Source/Target: {srcIPPs[0]}')
        else:
            print(f'   * different origins:\n      Source: {srcIPPs[0]}',
                  f'\n      Target: {trgIPPs[0]}')
            
        if sameDirs:
            print('   * same patient orientation\n',
                  f'      Source/Target: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                  f'                      {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                  f'                      {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]')
        else:
            print('   * different patient orientations:\n',
                  f'      Source: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                  f'               {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                  f'               {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]\n',
                  f'      Target: [{trgDirs[0]}, {trgDirs[1]}, {trgDirs[2]},\n',
                  f'               {trgDirs[3]}, {trgDirs[4]}, {trgDirs[5]},\n',
                  f'               {trgDirs[6]}, {trgDirs[7]}, {trgDirs[8]}]')
        
    
    if True:#p2c:
        #PrintTitle(f'Case {UseCase} applies.')
        if forceReg and useCaseThatApplies in ['3a', '3b', '4a', '4b']:
            useCaseToApply = useCaseThatApplies.replace('3', '5').replace('4', '5')
            
            msg = f'\n*UseCase {useCaseThatApplies} applies but since '\
                  + f'forceReg = {forceReg}, UseCase {useCaseToApply} will '\
                  + 'be applied.'
            
            print(msg)
        else:
            useCaseToApply = useCaseThatApplies
            
            print(f'\n*UseCase {useCaseThatApplies} applies.')
            
        print('-'*120)
        
    return useCaseThatApplies, useCaseToApply


""" From io_tools.propagate.py """

    def replace_indices_USE_OTHER(self, srcDataset, trgDataset):
        # TODO modify the docstrings
        """
        # 14/09/21: Use replace_ind_in_C2SindsByRoi in general_tools.shifting.py
        
        #13/09/21: This method is obsolete. shift_frames() adjusts the indices
        #appropriately.
        
        06/09/21: This was called replace_ind_in_f2sIndsByRoi prior to 06/09.
        
        Replace an index in a list of frame-to-slice indices for making
        (non-relationship-preserving) copies of ROIs.
        
        Parameters
        ----------
        srcDataset : Importer Object
            Importer Object for the source DICOM series.
        trgDataset : Importer Object
            Importer Object for the target DICOM series.
        
        Returns
        -------
        self.replacedF2SindsByRoi : list of list of ints
            List (for each ROI/segment) of a list (for each contour/frame) of
            the slice indices.
        
        Note
        ----
        The segmentation on srcDataset.slcNum is to be copied to 
        trgDataset.slcNum. The purpose of this method is to replace the
        source's slcNum in srcDataset.f2sIndsByRoi with that of the target.
        The modified F2SindsByRoi will be added to the new dataset.
        """
        
        srcF2SindsByRoi = srcDataset.f2sIndsByRoi
        indToReplace = srcDataset.slcNum
        replacementInd = trgDataset.slcNum
        
        # Start with copy of srcF2SindsByRoi for newF2SindsByRoi:
        #replacedF2SindsByRoi = deepcopy(srcF2SindsByRoi)
        newF2SindsByRoi = list(srcF2SindsByRoi)
        
        #print(f'type(replacedF2SindsByRoi) = {type(replacedF2SindsByRoi)}\n')
        
        for r in range(len(srcF2SindsByRoi)):
            for c in range(len(srcF2SindsByRoi[r])):
                if srcF2SindsByRoi[r][c] == indToReplace:
                    #replacedF2SindsByRoi[r][c] = replacementInd
                    newF2SindsByRoi[r][c] = replacementInd
        
        #self.replacedF2SindsByRoi = replacedF2SindsByRoi
        #self.f2sIndsByRoi = replacedF2SindsByRoi
        self.f2sIndsByRoi = newF2SindsByRoi
    
    def shift_frames_210913(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        06/09/21: This was called shift_frames_in_pixarrBySeg prior to 06/09.
            
        Shift the pixels in frames for (non-relationship-preserving) copies of
        ROIs. 
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        #isSrcResampled : bool
        #    True if the source pixel arrays have been resampled (i.e.
        #    srcPixarrByRoi --> resSrcPixarrByRoi), and False otherwise.
        
        Returns
        -------
        self.pixarrBySeg : list of Numpy Arrays
            The result of shifting the voxels in srcPixarrByRoi to make a non-
            relationship-preserving copy for target.
        self.f2sIndsBySeg : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.shiftedPixarrBySeg.
        """
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of shift_frames():')
            print('-'*120, '\n')
            
        #print(f'dir(self):\n{dir(self)}\n')
        
        # Have the source pixel arrays been resampled or not?
        #if hasattr(self, 'resPixarrByRoi'):
        if useCaseToApply in ['1', '2a', '2b']:
            slcNum = srcDataset.slcNum # initial value
            f2sIndsBySeg = srcDataset.f2sIndsBySeg # initial value
            pixarrBySeg = srcDataset.pixarrBySeg # initial value
        
        elif useCaseToApply in ['3a', '3b', '4a', '4b']:
            # Use the resampled source slice number for srcSlcNum, and the
            # resampled pixel arrays for srcPixarrByRoi:
            #srcSlcNum = self.resSlcNum # does this still exist? (06/09/21)
            #srcPixarrByRoi = self.resPixarrByRoi # does this still exist? (06/09/21)
            slcNum = self.slcNum # (06/09/21)
            f2sIndsBySeg = self.f2sIndsBySeg # (06/09/21)
            pixarrBySeg = self.pixarrBySeg # (06/09/21)
        
        """ Note:
        self.dcmIm was initialised in the Propagator class as trgIm, but 
        depending on the use case might have been replaced with the result of
        resampling srcIm to the target domain.
        """
        im = self.dcmIm
        
        trgIm = trgDataset.image
        trgSlcNum = trgDataset.slcNum
        
        #print(f'slcNum = {slcNum}')
        #print(f'trgSlcNum = {trgSlcNum}')
        
        """
        if not hasattr(self, 'resPixarrByRoi'):
            # Replace all indices in source with srcDataset.slcNum with
            # trgDataset.slcNum:
            #newF2SindsByRoi = deepcopy(srcDataset.f2sIndsBySeg)
            #self.replace_indices(srcDataset, trgDataset)
            self.replace_indices(srcDataset, trgDataset)
            replacedF2SindsByRoi = self.replacedF2SindsByRoi
            #f2sIndsBySeg = self.f2sIndsBySeg
        else:
            # Pass self.resF2SindsByRoi instead:
            replacedF2SindsByRoi = self.resF2SindsByRoi
        """
        
        """
        06/09/21: This is unnecessary since the index shift along z will 
        already be done below using voxShift!
        
        if useCaseToApply in ['1', '2a']:
            # Replace all indices in source with srcDataset.slcNum with
            # trgDataset.slcNum to get newF2SindsByRoi (self.f2sIndsBySeg):
            self.replace_indices(srcDataset, trgDataset)
        """
        
        # Although the f2sInds will be shifted, the pixel data will not, so
        # simply make a copy of the source's pixarrBySeg:
        #shiftedPixarrBySeg = deepcopy(srcDataset.pixarrBySeg)
        
        """
        # Get voxel shift between srcSlcNum in srcIm and trgSlcNum in trgIm
        # in the trgIm domain:
        voxShift = get_voxel_shift_bt_slices(
            image0=srcIm, sliceNum0=srcSlcNum,
            image1=trgIm, sliceNum1=trgSlcNum,
            refImage=trgIm
            )
        """
        
        # Get voxel shift (in the trgIm domain) between slice slcNum in im and
        # slice trgSlcNum in trgIm:
        voxShift = get_voxel_shift_bt_slices(
            image0=im, sliceNum0=slcNum,
            image1=trgIm, sliceNum1=trgSlcNum,
            refImage=trgIm
            )
        
        
        if p2c:
            print(f'   f2sIndsBySeg prior to shifting = {self.f2sIndsBySeg}')
            print(f'   voxShift that will be applied = {voxShift}')
        
        shiftedF2SindsBySeg = [] # initial value
        shiftedPixarrBySeg = [] # initial value
        
        # Loop through each pixel array:
        #for s in range(len(srcPixarrByRoi)):
        for s in range(len(pixarrBySeg)):
            # Proceed only if there is at least one frame in this pixel array:
            #if srcPixarrByRoi[s].shape[0]:
            if pixarrBySeg[s].shape[0]:
                # Replace srcPixarrByRoi[s] with the result of shifting the 
                # in-plane elements and add the pixel shift along z to
                # replacedF2SindsByRoi:
                """
                shiftedPixarrBySeg.append(
                    shift_frame(
                        frame=srcPixarrByRoi[s], voxShift=voxShift
                        )
                    )
                """
                shiftedPixarrBySeg.append(
                    shift_frame(
                        frame=pixarrBySeg[s], voxShift=voxShift
                        )
                    )
                
                #F = len(replacedF2SindsByRoi[s])
                F = len(f2sIndsBySeg[s])
                
                # Shift the frame-to-slice indices by voxShift[2] to account
                # for the z-shift:
                """
                shiftedF2SindsBySeg.append(
                    [replacedF2SindsByRoi[s][i] + voxShift[2] for i in range(F)]
                    #[f2sIndsBySeg[s][i] + voxShift[2] for i in range(F)]
                    )
                """
                shiftedF2SindsBySeg.append(
                    [f2sIndsBySeg[s][i] + voxShift[2] for i in range(F)]
                    )
                
                if p2c:
                    unique = get_unique_items(shiftedPixarrBySeg[s])
                    
                    print(f'   There are {len(unique)} unique items in',
                          f'shiftedPixarrBySeg[{s}] after shifting the frame')
        
        if p2c:
            print('   shiftedF2SindsBySeg after shifting =',
                  f'{shiftedF2SindsBySeg}')
            print('end of shift_frames\n')
            print('-'*120)
        
        """
        self.shiftedPixarrBySeg = shiftedPixarrBySeg
        self.shiftedF2SindsBySeg = shiftedF2SindsBySeg
        #self.pixarrBySeg = shiftedPixarrBySeg
        #self.f2sIndsBySeg = shiftedF2SindsBySeg
        """
        self.pixarrBySeg = shiftedPixarrBySeg
        self.f2sIndsBySeg = shiftedF2SindsBySeg
    
    def add_modified_data_to_existing_trgDataset_pre270921(self, trgDataset, params):
        """
        Concatenate a list (by ROI/segment) of points/pixel arrays and of 
        frame-to-slice indices (f2sIndsByRoi) to the corresponding instance
        variables in trgDataset.
        
        The relevant instance variables in trgDataset may be None, e.g. if 
        there was no ROI Collection imported for the target DICOM series
        i.e. trgDataset.roicol = None). If so f2sIndsByRoi and ptsByCntByRoi/
        pixarrByRoi will be shiftedF2SindsByRoi and shiftedPtsByCntByRoi/
        shiftedPixarrByRoi.
        
        Parameters
        ----------
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        self.pixarrByRoi : list of Numpy Arrays
            The list (for each ROI) of the pixel arrays of the non-relationship
            -preserving copy of the source ROI.
        self.f2sIndsByRoi : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.pixarrByRoi.
        """
        
        roicolMod = params.cfgDict['roicolMod']
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of add_modified_data_to_existing_trgDataset():')
            print('-'*120, '\n')
        
        if trgDataset.roicol == None:
            if p2c:
                print('   There is no ROI data to add from trgDataset')
                if roicolMod == 'RTSTRUCT':
                    print(f'   self.c2sIndsByRoi = {self.c2sIndsByRoi}')
                else:
                    print(f'   self.f2sIndsBySeg = {self.f2sIndsBySeg}')
                    #print('   self.pixarrBySeg = shiftedPixarrBySeg')
                print('-'*120)
        else:
            
            #shiftedF2SindsByRoi = self.shiftedF2SindsByRoi
            #shiftedPixarrByRoi = self.shiftedPixarrByRoi
            shiftedF2SindsBySeg = self.f2sIndsBySeg
            shiftedPixarrBySeg = self.pixarrBySeg
            shiftedC2SindsByRoi = self.c2sIndsByRoi
            shiftedPtsByCntByRoi = self.ptsByCntByRoi
            
            # Concatenate the existing data from target with the shifted data:
            trgF2SindsBySeg = trgDataset.f2sIndsBySeg
            trgPixarrBySeg = trgDataset.pixarrBySeg
            trgC2SindsByRoi = trgDataset.c2sIndsByRoi
            trgPtsByCntByRoi = trgDataset.ptsByCntByRoi
            
            if roicolMod == 'RTSTRUCT':
                R = len(trgC2SindsByRoi)
            else:
                R = len(trgF2SindsBySeg)
                
            newF2SindsBySeg = [
                shiftedF2SindsBySeg[r] + trgF2SindsBySeg[r] for r in range(R)
                ]
            
            newPixarrBySeg = [
                shiftedPixarrBySeg[r] + trgPixarrBySeg[r] for r in range(R)
                ]
            
            newC2SindsByRoi = [
                shiftedC2SindsByRoi[r] + trgC2SindsByRoi[r] for r in range(R)
                ]
            
            newPtsByCntByRoi = [
                shiftedPtsByCntByRoi[r] + trgPtsByCntByRoi[r] for r in range(R)
                ]
            
            self.f2sIndsBySeg = newF2SindsBySeg
            self.pixarrBySeg = newPixarrBySeg
            self.c2sIndsByRoi = newC2SindsByRoi
            self.ptsByCntByRoi = newPtsByCntByRoi
            
            if p2c:
                if roicolMod == 'RTSTRUCT':
                    print('   There was ROI data to add from trgDataset')
                    print(f'   shiftedC2SindsByRoi = {shiftedC2SindsByRoi}')
                    print('   after adding existing c2sIndsByRoi from trgDataset',
                          f'= {newC2SindsByRoi}')
                else:
                    print('   There was SEG data to add from trgDataset')
                    print(f'   shiftedF2SindsBySeg = {shiftedF2SindsBySeg}')
                    print('   after adding existing f2sIndsBySeg from trgDataset',
                          f'= {newF2SindsBySeg}')
                print('-'*120)
    
    def get_labimByRoi_MOVED(self, params, dataset):
        # TODO update docstrings
        """
        This has been moved to get_seg_metadata() in io_tools.import_data.
        
        Get a list of label images-by-ROI from a list of pixel arrays-by-ROI in
        a Dataset.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains various parameters.
        dataset : DataImporter Object
            DataImporter Object for a DICOM series.
            
        Returns
        -------
        None.
        """
        
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of get_labimByRoi():')
            print('-'*120, '\n')
        
        pixarrByRoi = dataset.pixarrByRoi
        f2sIndsByRoi = dataset.f2sIndsByRoi
        refIm = dataset.dcmIm
        
        labimByRoi = []
        f2sIndsByRoi_new = []
    
        for r in range(len(pixarrByRoi)):
            labim = pixarr_to_im(
                pixarr=pixarrByRoi[r], f2sInds=f2sIndsByRoi[r], refIm=refIm
                )
            
            labimByRoi.append(labim)
        
            if p2c:
                print(f'\n   Image info for labimByRoi[{r}]:')
                
            pixID, pixIDtypeAsStr, uniqueVals, f2sInds = \
                get_im_info(labim, p2c)
            
            f2sIndsByRoi_new.append(f2sInds)
        
        if p2c:
            print('-'*120)
            
        self.labimByRoi
        
    def export_labims(self, srcDataset, trgDataset, params):
        
        print('* Exporting label images..\n')
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        useDroForTx = cfgDict['useDroForTx']
        labimExportDir = cfgDict['labimExportDir']
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        # List of datasets:
        dsets = [srcDataset, trgDataset, self]
        
        dsetNames = ['src', 'trg', 'new']
        
        # List of labimByRoi to export:
        listOfLabimByRoi = [dset.labimByRoi for dset in dsets]
        
        # List of file names (will be modified in loop below):
        if useDroForTx:
            fnames = [
                f'{runID}_DRO_{name}Labim_{currentDateTime}' for name in dsetNames
                ]
        else:
            fnames = [
                f'{runID}_{name}Labim_{currentDateTime}' for name in dsetNames
                ]
        
        for d in range(len(listOfLabimByRoi)):
            labimByRoi = listOfLabimByRoi[d]
            fname = fnames[d]
            
            if labimByRoi: # if not None
                # Loop through each labim (for each ROI):
                for r in range(len(labimByRoi)):
                    export_im(
                        labimByRoi[r], filename=f'{fname}{r}',
                        fileFormat='HDF5ImageIO', exportDir=labimExportDir
                    )
    
    def export_tx_and_params(self, params):
        
        print('* Exporting transforms and transform parameters..\n')
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        useDroForTx = cfgDict['useDroForTx']
        txExportDir = cfgDict['txExportDir']
        
        # Create the export directory if it doesn't already exist:
        if not os.path.isdir(txExportDir):
            #os.mkdir(exportDir)
            Path(txExportDir).mkdir(parents=True)
        
        resTx = self.resTx
        initRegTx = self.initRegTx
        preRegTx = self.preRegTx
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        # List of transforms to export and their names:
        txs = [resTx, initRegTx, preRegTx]
        txNames = ['resTx', 'initRegTx', 'preRegTx']
        
        for tx, txName in zip(txs, txNames):
            if tx: # is not None
                if useDroForTx:
                    fname = f'{runID}_DRO_{txName}_{currentDateTime}'
                else:
                    fname = f'{runID}_{txName}_{currentDateTime}'
                
                # Export tx as a TFM:
                fpath = os.path.join(txExportDir, fname + '.tfm')
                sitk.WriteTransform(tx, fpath)
                
                # Export the tx parameters as a TXT:
                export_list_to_txt(
                    items=tx.GetParameters(),
                    filename=fname,
                    exportDir=txExportDir
                    )


""" From testing.test_download_data.py """

    def test_fetch_params_and_data_210812(runID):
        """
        Test fetch_params_and_data.
        
        Parameters
        ----------
        runID : str
            The ID that determines the main configuration parameters to use for the
            run (i.e. the key in the base level dictionary contained in 
            main_params.json).
        
        Returns
        -------
        Fetcher Object
            A Fetcher Object containing various parameters for the specified runID.
        
        Notes
        -----
        To run unit test:
            from testing.test_fetch_params_and_data import test_fetch_params_and_data
            
            params = test_fetch_params_and_data(runID)
        """
        
        cfgDir = r'C:\Code\WP1.3_multiple_modalities\Python_code\configs'
        
        print('Fetching parameters (locally) and data (from XNAT) without pre-',
              'existing XNAT session...\n')
        params = Fetcher(cfgDir, runID)
        
        #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
        #params = Fetcher(cfgDir, params.xnatSession)
        
        return params


""" From testing.test_fetch_config.py """

def test_fetch_params_and_data_210812(runID):
    """
    Test fetch_params_and_data.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    
    Returns
    -------
    Fetcher Object
        A Fetcher Object containing various parameters for the specified runID.
    
    Notes
    -----
    To run unit test:
        from testing.test_fetch_params_and_data import test_fetch_params_and_data
        
        params = test_fetch_params_and_data(runID)
    """
    
    cfgDir = r'C:\Code\WP1.3_multiple_modalities\Python_code\configs'
    
    print('Fetching parameters (locally) and data (from XNAT) without pre-',
          'existing XNAT session...\n')
    params = Fetcher(cfgDir, runID)
    
    #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
    #params = Fetcher(cfgDir, params.xnatSession)
    
    return params


""" From testing.test_fetch_params_and_data.py """

def test_fetch_params_and_data_210812(runID):
    """
    Test fetch_params_and_data.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    
    Returns
    -------
    Fetcher Object
        A Fetcher Object containing various parameters for the specified runID.
    
    Notes
    -----
    To run unit test:
        from testing.test_fetch_params_and_data import test_fetch_params_and_data
        
        params = test_fetch_params_and_data(runID)
    """
    
    cfgDir = r'C:\Code\WP1.3_multiple_modalities\Python_code\configs'
    
    print('Fetching parameters (locally) and data (from XNAT) without pre-',
          'existing XNAT session...\n')
    params = Fetcher(cfgDir, runID)
    
    #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
    #params = Fetcher(cfgDir, params.xnatSession)
    
    return params


""" From testing.test_import_data.py """

def test_import_data_210812(runID, params=None):
    """
    Run unit or integrated test of import_data.
    
    To run unit test, provide Fetcher Object. To run integrated test, leave 
    input argument blank.
    
    Parameters
    ----------
    params : Fetcher Object, optional
        A Fetcher object. The default is None.
        
    Returns
    -------
    srcDataset : Dataset Object
        Dataset Object for Source DICOM series.
    trgDataset : Dataset Object
        Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data(runID, params)
    
    To integrated test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data(runID)
    """
        
    if params == None:
        # Fetch the parameters:
        print('Fetching parameters (locally) and data (from XNAT)...\n')
        params = test_fetch_params_and_data(runID)
    
    # Import data for source and target:
    print('Importing source dataset...\n')
    srcDataset = Importer(params=params, srcORtrg='src')
    
    print('Importing target dataset...\n')
    trgDataset = Importer(params=params, srcORtrg='trg')
    
    return srcDataset, trgDataset

def test_import_data_210809(params=None):
    """
    Run unit or integrated test of import_data.
    
    To run unit test, provide Params Object. To run integrated test, leave 
    input argument blank.
    
    Parameters
    ----------
    params : Params Object, optional
        A Params object. The default is None.
        
    Returns
    -------
    srcDataset : Dataset Object
        Dataset Object for Source DICOM series.
    trgDataset : Dataset Object
        Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data(Params)
    
    To integrated test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data()
    """
    
    if params == None:
        # Fetch the parameters:
        params = test_fetch_params_and_data()
    
    # Import data for source and target:
    srcDataset = Dataset(xnatParams=params.srcXnatParams,
                         genParams=params.genParams)
    
    trgDataset = Dataset(xnatParams=params.trgXnatParams,
                         genParams=params.genParams)
    
    return srcDataset, trgDataset



""" From testing.test_propagate_data.py """

    def export_newRoiCol_MOVED(self):
        """ See dicom_tools.create_roicol.py """
        
        print('* Exporting the new Target ROI Collection..\n')
        
        srcDataset = self.srcDataset
        #trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        cfgDict = params.cfgDict
        roicolMod = cfgDict['roicolMod']
        runID = cfgDict['runID']
        useDroForTx = cfgDict['useDroForTx']
        
        if roicolMod == 'RTSTRUCT': 
            exportDir = cfgDict['rtsExportDir']
        else:
            exportDir = cfgDict['segExportDir']
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        if useDroForTx:
            fname = f'{runID}_DRO_newRoicol_{currentDateTime}'
        else:
            fname = f'{runID}_newRoicol_{currentDateTime}'
        
        newRoicolFpath = export_newRoicol(
            newRoicol=newDataset.roicol, 
            srcRoicolFpath=srcDataset.roicolFpath, 
            exportDir=exportDir, fname=fname
            )
        
        self.newDataset.roicolFpath = newRoicolFpath
    
    def plot_metric_v_iters_MOVED(self):
        """ 
        See io_tools.propagate.py 
        
        Plot the metric values v iterations from the registeration.
        
        Note that metric values will only exist if registration was performed,
        so this plot will only be produced if useDroForTx = False.
        """
        
        params = self.params
        cfgDict = params.cfgDict
        useDroForTx = cfgDict['useDroForTx']
        
        if not useDroForTx:
            runID = cfgDict['runID']
            exportPlot = cfgDict['exportPlots']
            #p2c = cfgDict['p2c']
            useCaseToApply = cfgDict['useCaseToApply']
            #forceReg = cfgDict['forceReg']
            
            regTxName = cfgDict['regTxName']
            initMethod = cfgDict['initMethod']
            resInterp = cfgDict['resInterp']
            
            resExportDir = cfgDict['resPlotsExportDir']
            
            #print(useCaseToApply)
            
            metricValues = self.newDataset.metricValues
            multiresIters = self.newDataset.multiresIters
            
            print('* Plotting metric v iterations from registeration..\n')
            
            # Prepare plot title:
            if useCaseToApply in ['5a', '5b'] and not useDroForTx:
                resTitle = f'Metric v iters for Src reg to Trg ({regTxName}, '\
                    + f'{initMethod}, {resInterp})'
                
                # Prepare filename for exported plot:
                fname = f'{runID}_' + resTitle.replace(' ', '_').replace(',', '')
                fname = fname.replace('(', '').replace(')', '')
                
                plot_metricValues_v_iters(
                    metricValues=metricValues, multiresIters=multiresIters, 
                    exportPlot=exportPlot, exportDir=resExportDir,
                    fname=fname
                    )
    
    def plot_res_results_MOVED(self):
        """ 
        See io_tools.propagate.py
        
        Plot a single slice from trgIm and resIm and compare to assess result
        of resampling/registering.
        """
        
        print('* Plotting resampled/registered results..\n')
        
        #srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        #roicolMod = cfgDict['roicolMod']
        exportPlot = cfgDict['exportPlots']
        #p2c = cfgDict['p2c']
        useCaseToApply = cfgDict['useCaseToApply']
        #forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        trgIm = trgDataset.image
        resIm = newDataset.image # resampled source or source-registered-to-target 
        alignedIm = newDataset.alignedIm # pre-registration aligned image
        # (which may be None, e.g. if registration wasn't performed)
        
        resExportDir = cfgDict['resPlotsExportDir']
        
        #print(useCaseToApply)
        
        # Prepare plot title for the aligned image:
        aliTitle = f'Src aligned to Trg ({initMethod})'
        
        # Prepare plot title for the resampled/registered image:
        if useCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
            resTitle = 'Src '
            if useCaseToApply in ['3a', '3b', '4a', '4b']:
                resTitle += f'res to Trg ({resInterp})'
            elif useCaseToApply in ['5a', '5b']:
                if useDroForTx:
                    resTitle += f'tx to Trg ({regTxName} DRO, {resInterp})'
                else:
                    resTitle += f'reg to Trg ({regTxName}, {initMethod}, '\
                        + f'{resInterp})'
            
            midInd = trgIm.GetSize()[2] // 2
            
            # List of images to plot and the plot titles:
            images = [alignedIm, resIm]
            titles = [aliTitle, resTitle]
            
            for image, title in zip(images, titles):
                if image: # i.e. not None
                    # Prepare filename for exported plot:
                    fname = f'{runID}_' + \
                        title.replace(' ', '_').replace(',', '')
                    fname = fname.replace('(', '').replace(')', '')
                    
                    compare_res_results(
                        resIm0=trgIm, resIm1=image, resInd=midInd,
                        resTitle0='Target image', resTitle1=title,
                        exportPlot=exportPlot, exportDir=resExportDir,
                        fname=fname
                    )
    
    def plot_roi_over_dicom_ims_MOVED(self):
        """ 
        See io_tools.propagate.py
        
        Plot copied/propagated ROI overlaid on DICOM images 
        """
        
        print('* Plotting ROI over DICOM images..\n')
        
        srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        roicolMod = cfgDict['roicolMod']
        exportPlot = cfgDict['exportPlots']
        p2c = cfgDict['p2c']
        useCaseToApply = cfgDict['useCaseToApply']
        forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        #trgIm = trgDataset.image
        #resIm = newDataset.image # resampled source or source-registered-to-target 
        
        #resExportDir = cfgDict['resPlotsExportDir']
        
        if roicolMod == 'RTSTRUCT':
            roiExportDir = cfgDict['rtsPlotsExportDir']
        else:
            roiExportDir = cfgDict['segPlotsExportDir']
        
        # Prepare plot title for new dataset:
        method = 'Src '
        if useCaseToApply in ['1', '2a']:
            method += 'copied to Trg'
        elif useCaseToApply == '2b':
            method += 'propagated to Trg'
        elif useCaseToApply in ['3a', '4a']:
            method += f'copied to Trg (res, {resInterp})'
        elif useCaseToApply in ['3b', '4b']:
            method += f'propagated to Trg (res, {resInterp})'
        elif useCaseToApply == '5a':
            if useDroForTx:
                method += f'copied to Trg ({regTxName} DRO, {resInterp})'
            else:
                method += f'copied to Trg ({regTxName} reg, {initMethod}, '\
                    + f'{resInterp})'
        elif useCaseToApply == '5b':
            if useDroForTx:
                method += f'propagated to Trg ({regTxName} DRO, {resInterp})'
            else:
                method += f'propagated to Trg ({regTxName} reg, {initMethod},'\
                    + f' {resInterp})'
        
        #print(method)
        
        listOfSegs = [srcDataset.roicol, newDataset.roicol]
        """
        Note: newDataset does not have a dicomDir since DICOMs are not
        generated for the resampled/registered source dataset, but for the
        purposes of pulling metadata (sizes, spacings, IPPs, directions, etc)
        relevant to the resampled/registered dataset, getting it from the
        target dataset is equivalent.
        """
        listOfDicomDirs = [srcDataset.dicomDir, trgDataset.dicomDir]
        
        listOfImages = [srcDataset.image, newDataset.image]
        listOfPlotTitles = ['Src', method]
        
        """
        plot_pixarrs_from_list_of_segs_v3(
        listOfSegs, listOfDicomDirs, listOfPlotTitles, exportPlot=exportPlot,
            exportDir=exportDir, runID=runID, useCaseToApply=useCaseToApply,
            forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
            initMethod=initMethod, resInterp=resInterp, p2c=p2c
        )
        """
        
        # Prepare filename for exported plot:
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        if useDroForTx:
            fname = f'{runID}_DRO'
        else:
            fname = f'{runID}_'
        fname += method.replace(' ', '_').replace(',', '')
        fname = fname.replace('(', '').replace(')', '')
        fname += f'_{currentDateTime}'
            
        plot_pixarrs_from_list_of_segs_and_images(
            listOfSegs, listOfImages, listOfDicomDirs, listOfPlotTitles, 
            exportPlot=exportPlot, exportDir=roiExportDir,
            runID=runID, useCaseToApply=useCaseToApply,
            forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
            initMethod=initMethod, resInterp=resInterp, 
            fname=fname, p2c=p2c
            #fname='', p2c=p2c
        )

    def export_labims_MOVED(self):
        """ See io_tools.propagate.py """
        
        print('* Exporting label images..\n')
        
        srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        useDroForTx = cfgDict['useDroForTx']
        labimExportDir = cfgDict['labimExportDir']
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        # List of datasets:
        dsets = [srcDataset, trgDataset, newDataset]
        
        dsetNames = ['src', 'trg', 'new']
        
        # List of labimByRoi to export:
        listOfLabimByRoi = [dset.labimByRoi for dset in dsets]
        
        # List of file names (will be modified in loop below):
        if useDroForTx:
            fnames = [
                f'{runID}_DRO_{name}Labim_{currentDateTime}' for name in dsetNames
                ]
        else:
            fnames = [
                f'{runID}_{name}Labim_{currentDateTime}' for name in dsetNames
                ]
        
        for d in range(len(listOfLabimByRoi)):
            labimByRoi = listOfLabimByRoi[d]
            fname = fnames[d]
            
            if labimByRoi: # if not None
                # Loop through each labim (for each ROI):
                for r in range(len(labimByRoi)):
                    export_im(
                        labimByRoi[r], filename=f'{fname}{r}',
                        fileFormat='HDF5ImageIO', exportDir=labimExportDir
                    )
    
    def export_tx_and_params_MOVED(self):
        """ See io_tools.propagate.py """
        
        print('* Exporting transforms and transform parameters..\n')
        
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        useDroForTx = cfgDict['useDroForTx']
        txExportDir = cfgDict['txExportDir']
        
        # Create the export directory if it doesn't already exist:
        if not os.path.isdir(txExportDir):
            #os.mkdir(exportDir)
            Path(txExportDir).mkdir(parents=True)
        
        resTx = self.newDataset.resTx
        initRegTx = self.newDataset.initRegTx
        preRegTx = self.newDataset.preRegTx
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        # List of transforms to export and their names:
        txs = [resTx, initRegTx, preRegTx]
        txNames = ['resTx', 'initRegTx', 'preRegTx']
        
        for tx, txName in zip(txs, txNames):
            if tx: # is not None
                if useDroForTx:
                    fname = f'{runID}_DRO_{txName}_{currentDateTime}'
                else:
                    fname = f'{runID}_{txName}_{currentDateTime}'
                
                # Export tx as a TFM:
                fpath = os.path.join(txExportDir, fname + '.tfm')
                sitk.WriteTransform(tx, fpath)
                
                # Export the tx parameters as a TXT:
                export_list_to_txt(
                    items=tx.GetParameters(),
                    filename=fname,
                    exportDir=txExportDir
                    )
        
def test_propagator_210812(runID, params=None, srcDataset=None, trgDataset=None):
    """
    Run unit or integrated test of propagate.
    
    To run unit test, provide Importer Objects for source and target. To run 
    integrated test, leave optional input arguments blank.
    
    Parameters
    ----------
    params : Fetcher Object, optional
        A Fetcher object containing various parameters. The default is None.
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
        
    Returns
    -------
    newDataset : Dataset Object
        New Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_propagate_data import test_propagator
        
        newDataset = test_propagator(
            runID, params, srcDataset, trgDataset
            )
    
    To integrated test:
        from testing.test_propagate_data import test_propagator
        
        newDataset = test_propagator(runID)
    """
    
    if params == None:
        # Fetch the parameters:
        print('Fetching parameters (locally) and data (from XNAT)...\n')
        params = test_fetch_params_and_data(runID)
    
    if srcDataset == None:
        # Import data for source:
        print('Importing source dataset...\n')
        srcDataset = Importer(params=params, srcORtrg='src')
    
    if trgDataset == None:
        # Import data for target:
        print('Importing target dataset...\n')
        trgDataset = Importer(params=params, srcORtrg='trg')
    
    return srcDataset, trgDataset



""" From xnat_tools.im_assessors.py """

def download_im_asr_210810(xnatCfg, xnatParams, genParams, xnatSession=None,
                    pathsDict=None):
    """
    Download an image assessor (ROI Collection) of interest for a particular 
    scan from XNAT.
    
    Parameters
    ----------
    xnatCfg : dict
        Dictionary of XNAT configuration settings.
    xnatParams : dict
        Dictionary of XNAT parameters.
    genParams : dict
        Dictionary of general parameters.
    xnatSession : requests.models.Response, optional
        If provided a new session request will be avoided.
    pathsDict : dict, optional (None by default)
        Dictionary containing paths of data downloaded. 
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    
    Notes
    -----
    Although the structure of the REST call below:
    
    assessor = session.get(f'{url}/data/projects/{projID}/'
                           f'subjects/{subjLab}/'
                           f'experiments/{expLab}/'
                           f'assessors/{AsrID}/resources/{AsrMod}/files')
    
    follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Also, Jason suggested adding the file name after 'files', which worked. 
    """
    
    if xnatSession == None:
        xnatSession = create_session(xnatCfg)
    
    url = xnatSession.url
    
    rootExportDir = genParams['rootExportDir']
    
    projID = xnatParams['projID']
    subjLab = xnatParams['subjLab']
    expLab = xnatParams['expLab']
    scanID = xnatParams['scanID']
    roicolMod = xnatParams['roicolMod']
    roicolName = xnatParams['roicolName']
    
    p2c = genParams['p2c']
    
    if p2c:
        print('Inputs to download_im_asr():')
        print(f'  projID = {projID}')
        print(f'  subjLab = {subjLab}')
        print(f'  expLab = {expLab}')
        print(f'  scanID = {scanID}')
        print(f'  roicolMod = {roicolMod}')
        print(f'  roicolName = {roicolName}\n')
    
    """ The 'collectionType' key in 'data_fields' is either 'AIM' or 'SEG'. """
    if roicolMod == 'RTSTRUCT':
        collType = 'AIM'
    
    """ Get the experiment of interest: """
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}?format=json'
                              
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    experiment = request.json()
    
    """ There are likely more than one item in the list 
    experiment['items'][0]['children'] (e.g. one for scans and one for 
    image assessors). Get the index that corresponds to assessors, asr_ind, and
    the index that corresponds to scans, scan_ind: """
    asrInd = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'assessor' in experiment['items'][0]['children'][i]['field']:
            asrInd = i
    
    if asrInd == None:
        msg = f"There are no assessors for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    scanInd = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'scan' in experiment['items'][0]['children'][i]['field']:
            scanInd = i
    
    if scanInd == None:
        msg = f"There are no scans for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    if p2c:
        print(f'asrInd = {asrInd}')
        print(f'scanInd = {scanInd}\n')
    
    """ Get the SeriesUID for the scan of interest: """
    seriesUID = '' # initial value
    
    if pathsDict != None:
        try:
            seriesUID\
                = pathsDict['projects'][projID]['subjects'][subjLab]\
                    ['experiments'][expLab]['scans'][scanID]['seriesUID']
            
            excRaised = False # exception raised
            
        except Exception:
            excRaised = True
    
    if pathsDict == None or excRaised:
        for scan in experiment['items'][0]['children'][scanInd]['items']:
            if p2c:
                print(f"ID = {scan['data_fields']['ID']}")
                print(f"UID = {scan['data_fields']['UID']}\n")
            
            if scan['data_fields']['ID'] == scanID:
                seriesUID = scan['data_fields']['UID']
                
                if p2c:
                    print(f'seriesUID = {seriesUID}\n')
    
    if seriesUID == '':
        msg = f"No image assessors were found for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    """ Get the assessor ID and label of interest: """
    roicolInd = None # initial value
    asrID = None # initial value
    label = None # initial value
    asrDate = None # initial value
    asrTime = None # initial value
    
    for i in range(len(experiment['items'][0]['children'][asrInd]['items'])):
        if p2c:
            print(f"Item no {i}:")
        
        children = experiment['items'][0]['children'][asrInd]['items'][i]\
            ['children']
        
        if p2c:
            print(f"  children has {len(children)} items")
        
        fields = [child['field'] for child in children]
        
        # Get the index corresponding to the item containing the fields
        # 'out/file' and 'references/seriesUID':
        if 'out/file' in fields:
            fileInd = fields.index('out/file')
        else:
            fileInd = None
        
        if 'references/seriesUID' in fields:
            refInd = fields.index('references/seriesUID')
        else:
            refInd = None
        
        if p2c:
            print(f"  fields = {fields}")
            print(f"  fileInd = {fileInd}")
            print(f"  refInd = {refInd}\n")
        
        """ The following check is probably not required in most cases, but 
        covers cases where a ROI collection has been deleted and although the
        item with field 'out/file' is removed from the JSON, the item with
        field 'references/seriesUID' can persist. """
        # Check that items exist with both 'out/file' and 
        # 'references/seriesUID' fields:
        #if 'out/file' in fields and 'references/seriesUID' in fields:
        if isinstance(fileInd, int) and isinstance(refInd, int):
            if p2c:
                print("  This item's children contains items with fields",
                      "'out/file' and 'references/seriesUID'.\n")
            
            uid = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['children'][refInd]['items'][0]['data_fields']['seriesUID']
            
            type = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['data_fields']['collectionType']
            
            name = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['data_fields']['name']
            
            if p2c:
                print(f"  uid = {uid}")
                print(f"  type = {type}")
                print(f"  name = {name}\n")
            
                #if uid == seriesUID:
                #    print(f"   * Matched collType = {type}")
                #    print(f"   * Matched name = {name}\n")
            
            if (uid == seriesUID and type == collType and name == roicolName):
                if p2c:
                    print(f"   * Matched seriesUID = {seriesUID}")
                    print(f"   * Matched collType = {collType}")
                    print(f"   * Matched name = {name}\n")
                
                roicolInd = i
                
                asrID = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['id']
                
                label = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['label']
                
                asrDate = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['date']
                
                asrTime = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['time']
        else:
            if p2c:
                print("  This item's children DOES NOT contain items with",
                      "both field 'out/file' and 'references/seriesUID'.\n")
    
    if p2c:
        print(f"Index of the ROI Collection of interest = {roicolInd}")
        print(f"ID of the ROI Collection of interest = {asrID}")
        print(f"Label of the ROI Collection of interest = {label}")
        print(f"Date and time of the ROI Collection of interest = {asrDate}",
              f"{asrTime}")
    
    if roicolInd == None:
        msg = f"No image assessors were found for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'\ncontaining matches on:"\
              + f"\n  roicolMod = '{roicolMod}'\n  roicolName = '{roicolName}'."
        raise Exception(msg)
    
    """ Get the assessor of interest: """
    fname = label + '.dcm'
    
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/assessors/{asrID}/'\
          + f'resources/{roicolMod}/files/{fname}'
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    exportDir = os.path.join(rootExportDir, 'projects', projID, 'subjects',
                             subjLab, 'experiments', expLab, 'assessors',
                             asrID, 'resources', roicolMod, 'files')

    if not os.path.isdir(exportDir):
        Path(exportDir).mkdir(parents=True)
        print(f'Created directory:\n {exportDir}\n')
    
    fpath = os.path.join(exportDir, fname)
    
    with open(fpath, 'wb') as file:
        file.write(request.content)
    
    print(f'Image assessor downloaded to:\n {fpath}\n')
    
    """ Store info in a dictionary: """
    pathsDict, keys\
        = create_pathsDict_for_im_asr(projID, subjLab, expLab, asrID,
                                      roicolMod, fname, pathsDict)
    
    pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
    [expLab]['assessors'][asrID]['resources'][roicolMod]['files'][fname]\
    .update({'asrDate' : asrDate,
             'asrTime' : asrTime,
             'roicolName' : roicolName,
             'label' : label,
             'asrID' : asrID,
             'roicolMod' : roicolMod,
             'roicolDir' : exportDir,
             'roicolFpath' : fpath})
    
    return pathsDict

def download_im_asr_Pre_210810(url, projId, subjLab, expLab, scanId, mod, 
                    asrName, session=None, rootExportDir='', 
                    pathsDict=None, username=None, password=None, 
                    p2c=False):
    """
    Download an image assessor (ROI Collection) of interest for a particular 
    scan from XNAT.
    
    Parameters
    ----------
    url : str
        Address of XNAT (e.g. 'http://10.1.1.20').
    projId : str
        The project ID of interest.
    subjLab : str
        The subject label of interest. 
    expLab : str
        The DICOM study / XNAT experiment label of interest. 
    scanId : str
        The DICOM series label / XNAT scan ID of interest.
    mod : str
        The modality of the image assessor (ROI Collection) of interest.
    asrName : str
        The name of the image assessor (ROI Collection) of interest.
    session : requests session, optional
        If provided a new session request will be avoided.
    rootExportDir : str, optional
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, subjLab, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "xnat_downloads" within the default downloads directory.
    pathsDict : dict, optional
        Dictionary containing paths of data downloaded.
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = None)
        the user will be prompted to enter a user name.
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = None)
        the user will be prompted to enter a password.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    session : requests session
    
    Notes
    -----
    
    Although the structure of the REST call below:
    
    assessor = session.get(f'{url}/data/projects/{projId}/'
                           f'subjects/{subjLab}/'
                           f'experiments/{expLab}/'
                           f'assessors/{AsrId}/resources/{AsrMod}/files')
    
    follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Also, Jason suggested adding the file name after 'files', which worked. 
    """
    
    #p2c = True
    
    if p2c:
        print('Inputs to download_im_asr():')
        print(f'  projId = {projId}')
        print(f'  subjLab = {subjLab}')
        print(f'  expLab = {expLab}')
        print(f'  scanId = {scanId}')
        print(f'  mod = {mod}')
        print(f'  asrName = {asrName}\n')
    
    """ The 'collectionType' key in 'data_fields' is either 'AIM' or 'SEG'. """
    if mod == 'RTSTRUCT':
        collType = 'AIM'
    
    if session == None:
        session = create_session(url, username, password)
    
    if rootExportDir == '':
        rootExportDir = os.path.join(Path.home(), "Downloads", 
                                     "xnat_downloads")
    
    """ Get the experiment of interest: """
    uri = f'{url}/data/projects/{projId}/subjects/{subjLab}/'\
          + f'experiments/{expLab}?format=json'
                              
    request = session.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    experiment = request.json()
    
    """ There are likely more than one item in the list 
    experiment['items'][0]['children'] (e.g. one for scans and one for 
    image assessors). Get the index that corresponds to assessors, asr_ind, and
    the index that corresponds to scans, scan_ind: """
    asr_ind = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'assessor' in experiment['items'][0]['children'][i]['field']:
            asr_ind = i
    
    if asr_ind == None:
        msg = f"There are no assessors for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
        raise Exception(msg)
    
    scan_ind = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'scan' in experiment['items'][0]['children'][i]['field']:
            scan_ind = i
    
    if scan_ind == None:
        msg = f"There are no scans for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
        raise Exception(msg)
    
    if p2c:
        print(f'asr_ind = {asr_ind}')
        print(f'scan_ind = {scan_ind}\n')
    
    """ Get the seriesUid for the scan of interest: """
    seriesUid = '' # initial value
    
    if pathsDict != None:
        try:
            seriesUid = pathsDict['projects'][projId]['subjects'][subjLab]\
                        ['experiments'][expLab]['scans'][scanId]['seriesUid']
            
            exception_raised = False
            
        except Exception:
            exception_raised = True
    
    if pathsDict == None or exception_raised:
        for scan in experiment['items'][0]['children'][scan_ind]['items']:
            if p2c:
                print(f"ID = {scan['data_fields']['ID']}")
                print(f"UID = {scan['data_fields']['UID']}\n")
            
            if scan['data_fields']['ID'] == scanId:
                seriesUid = scan['data_fields']['UID']
                
                if p2c:
                    print(f'seriesUid = {seriesUid}\n')
    
    
    if seriesUid == '':
        msg = f"No image assessors were found for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
        raise Exception(msg)
    
    
    """ Get the assessor ID and label of interest: """
    roicolInd = None # initial value
    asrId = None # initial value
    label = None # initial value
    date = None # initial value
    time = None # initial value
    
    for i in range(len(experiment['items'][0]['children'][asr_ind]['items'])):
        if p2c:
            print(f"Item no {i}:")
        
        children = experiment['items'][0]['children'][asr_ind]['items'][i]\
                   ['children']
        
        if p2c:
            print(f"  children has {len(children)} items")
        
        fields = [child['field'] for child in children]
        
        # Get the index corresponding to the item containing the fields
        # 'out/file' and 'references/seriesUid':
        if 'out/file' in fields:
            fileInd = fields.index('out/file')
        else:
            fileInd = None
        
        if 'references/seriesUid' in fields:
            refInd = fields.index('references/seriesUid')
        else:
            refInd = None
        
        if p2c:
            print(f"  fields = {fields}\n")
        
        """ The following check is probably not required in most cases, but 
        covers cases where a ROI collection has been deleted and although the
        item with field 'out/file' is removed from the JSON, the item with
        field 'references/seriesUid' can persist. """
        # Check that items exist with both 'out/file' and 
        # 'references/seriesUid' fields:
        #if 'out/file' in fields and 'references/seriesUid' in fields:
        if isinstance(fileInd, int) and isinstance(refInd, int):
            if p2c:
                print("  This item's children contains an items with fields",
                      "'out/file' and 'references/seriesUid'.\n")
            
            seriesUid = experiment['items'][0]['children'][asr_ind]['items'][i]\
                        ['children'][refInd]['items'][0]['data_fields']\
                        ['seriesUid']
        
            collType = experiment['items'][0]['children'][asr_ind]['items'][i]\
                       ['data_fields']['collectionType']
            
            name = experiment['items'][0]['children'][asr_ind]['items'][i]\
                   ['data_fields']['name']
            
            if p2c:
                print(f"seriesUid = {seriesUid}")
                print(f"collType = {collType}")
                print(f"name = {name}\n")
            
                #if seriesUid == seriesUid:
                #    print(f"   * Matched collectionType = {collectionType}")
                #    print(f"   * Matched name = {name}\n")
            
            if seriesUid == seriesUid and collType == collType and name == asrName:
                if p2c:
                    print(f"   * Matched seriesUid = {seriesUid}")
                    print(f"   * Matched collType = {collType}")
                    print(f"   * Matched name = {name}\n")
                
                roicolInd = i
                
                asrId = experiment['items'][0]['children'][asr_ind]['items'][i]\
                         ['data_fields']['id']
                
                label = experiment['items'][0]['children'][asr_ind]['items'][i]\
                        ['data_fields']['label']
                
                date = experiment['items'][0]['children'][asr_ind]['items'][i]\
                       ['data_fields']['date']
                
                time = experiment['items'][0]['children'][asr_ind]['items'][i]\
                       ['data_fields']['time']
        else:
            if p2c:
                print("  This item's children DOES NOT contain an items with",
                      "both field 'out/file' and 'references/seriesUid'.\n")
    
    
    if p2c:
        print(f"Index of the ROI Collection of interest = {roicolInd}")
        print(f"ID of the ROI Collection of interest = {asrId}")
        print(f"Label of the ROI Collection of interest = {label}")
        print(f"Date and time of the ROI Collection of interest = {date} {time}")
    
    if roicolInd == None:
        msg = f"No image assessors were found for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'\ncontaining matches on:"\
              + f"\n  mod = '{mod}'\n  asrName = '{asrName}'."
        raise Exception(msg)
    
    
    """ Get the assessor of interest: """
    fname = label + '.dcm'
    
    uri = f'{url}/data/projects/{projId}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/assessors/{asrId}/'\
          + f'resources/{mod}/files/{fname}'
    
    request = session.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    export_dir = os.path.join(rootExportDir, 'projects', projId, 'subjects',
                              subjLab, 'experiments', expLab, 'assessors',
                              asrId, 'resources', mod, 'files')

    if not os.path.isdir(export_dir):
        Path(export_dir).mkdir(parents=True)
        print(f'Created directory:\n {export_dir}\n')
    
    fpath = os.path.join(export_dir, fname)
    
    with open(fpath, 'wb') as file:
        file.write(request.content)
    
    print(f'Image assessor downloaded to:\n {fpath}\n')
    
    
    """ Store info in a dictionary: """
    pathsDict, keys\
        = create_paths_dict_for_im_asr(projId, subjLab, expLab, asrId, mod,
                                       fname, pathsDict)
    
    pathsDict['projects'][projId]['subjects'][subjLab]['experiments']\
    [expLab]['assessors'][asrId]['resources'][mod]['files'][fname]\
    .update({'date' : date,
             'time' : time,
             'asrName' : asrName,
             'label' : label,
             'asrId' : asrId,
             'mod' : mod,
             'dir' : export_dir,
             'fpath' : fpath})
    
    return pathsDict, session



""" From xnat_tools.scans.py """

def download_scan_210809(xnatCfg, xnatParams, genParams, xnatSession=None,
                  pathsDict=None):
    """
    Download scan from XNAT.
    
    Returns a dictionary containing paths mimacking XNAT file hierarchy and the 
    XNAT requests session.
    
    Parameters
    ----------
    xnatCfg : dict
        Dictionary of XNAT configuration settings.
    xnatParams : dict
        Dictionary of XNAT parameters.
    genParams : dict
        Dictionary of general parameters
    xnatSession : requests.models.Response, optional
        If provided a new session request will be avoided.
    pathsDict : dict, optional (None by default)
        Dictionary containing paths of data downloaded. 
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    
    Notes
    -----
    GET - /data/projects/{project-id}/subjects/{subject-id | subject-label}/
           experiments/{experiment-id | experiment-label}/scans/{scan-id}/files
           
    GET - /data/experiments/{experiment-id}/scans/{scan-id}/files
    
    https://wiki.xnat.org/display/XAPI/Image+Scan+Resource+API
    """
    
    
    """ Set to True if issue with zip file (e.g. duplication of files): """
    OVERWRITE_ZIP = False
    
    if xnatSession == None:
        xnatSession = create_session(xnatCfg)
    
    url = xnatSession.url
    
    rootExportDir = genParams['rootExportDir']
    
    projID = xnatParams['projID']
    subjLab = xnatParams['subjLab']
    expLab = xnatParams['expLab']
    scanID = xnatParams['scanID']
    
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/scans/{scanID}/'\
          + 'resources/DICOM/files?format=zip'
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    exportDir = os.path.join(rootExportDir, 'projects', projID, 
                             'subjects', subjLab, 'experiments')

    if not os.path.isdir(exportDir):
        Path(exportDir).mkdir(parents=True)
        print(f'Created directory:\n {exportDir}\n')
    
    filepath = os.path.join(exportDir,
                            f'Experiment_{expLab}__Scan_{scanID}.zip')
    
    if not os.path.exists(filepath) or OVERWRITE_ZIP:
        with open(filepath, 'wb') as file:
            file.write(request.content)
            
        print(f'Zipped file downloaded to: \n{filepath}\n')
        
        #zip_file = ZipFile(filepath, mode='r')
        #zip_file.extractall(exportDir)
        
        with ZipFile(filepath, 'r') as file:
            file.extractall(exportDir)
        
        print(f'Zipped file extracted to: \n{exportDir}\n')
    else:
        if not OVERWRITE_ZIP:
            print(f'Zipped file already downloaded to: \n{filepath}\n')
    
    
    """ Get the directory name of the exported DICOMs: """
    # Get the filepath of the first file in the zip file:
    with ZipFile(filepath, 'r') as file:
        filepath = file.infolist()[0].filename
    
    fullFpath = os.path.join(exportDir, filepath)
    
    #print(f'fullFpath = {fullFpath}\n')
    
    dicomDir = os.path.split(fullFpath)[0]
    
    print(f"DICOMs downloaded to:\n{dicomDir}\n")
    
    """ Store info in a dictionary: """
    pathsDict, keys\
        = create_pathsDict_for_scan(projID, subjLab, expLab, scanID, pathsDict)
    
    print(f'pathsDict = {pathsDict}\n')
    
    if not 'dir' in keys:
        """ Get the experiment of interest: """
        uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
              + f'experiments/{expLab}?format=json'
        
        request = xnatSession.get(uri)
        
        # Raise status error if not None:
        if request.raise_for_status() != None:
            print(request.raise_for_status())
            
        experiment = request.json()
        
        #print(experiment['items'], '\n')
        #print(experiment['items'][0], '\n')
        #print(experiment['items'][0]['children'], '\n')
        #print(experiment['items'][0]['children'][1], '\n')
        #print(experiment['items'][0]['children'][1]['items'], '\n')
        #print(experiment['items'][0]['children'][1]['items'][0], '\n')
        #print(experiment['items'][0]['children'][1]['items'][0]['data_fields'], '\n')

        #studyUID = experiment['items'][0]['children'][1]['items'][0]['data_fields']['type']['UID']
        #studyUID = experiment['items'][0]['children'][1]['items'][0]['data_fields']['UID'] # 27/05/2021
        studyUID = experiment['items'][0]['data_fields']['UID'] # 27/05/2021
        
        """ There may be more than one item in the list 
        experiment['items'][0]['children'] (e.g. depending on whether there are 
        image assessors or not). Get the index that corresponds to scans: """
        scanInd = None # initial value
        
        for i in range(len(experiment['items'][0]['children'])):
            if 'scan' in experiment['items'][0]['children'][i]['field']:
                scanInd = i
        
        if scanInd == None:
            msg = f"There are no scans for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
            raise Exception(msg)
            
        
        """ Get the list of series labels (scan IDs), seriesUIDs and series 
        descriptions for this series/scan: """
        scanIDs = []
        seriesUIDs = []
        seriesDescs = []
        
        #for exp in experiment['items'][0]['children'][1]['items']: # 27/05/21
        for exp in experiment['items'][0]['children'][scanInd]['items']: # 27/05/21
            scanIDs.append(exp['data_fields']['ID'])
            seriesUIDs.append(exp['data_fields']['UID'])
            seriesDescs.append(exp['data_fields']['type'])
            
        
        seriesUID = seriesUIDs[scanIDs.index(scanID)]
        seriesDesc = seriesDescs[scanIDs.index(scanID)]
        
        """ Replacing code below with something more reliable (i.e. rather than
        replacing spaces, dashes and full stops with underscores) (01/06/21).
        """
            
        #""" The folder name between the 'scans' folder and 'resources' folder
        #has the form 'scanID-seriesDesc', but it appears that spaces, dashes and 
        #full stops within seriesDesc become underscores in the folder structure 
        #of the extracted zip file. """
        ##print(f'seriesUID = {seriesUID}')
        ##print(f'seriesDesc = {seriesDesc}\n')
        #
        ##DirName = f'{scanID}-{seriesDesc}'.replace(' ', '_').replace('.', '_') # 01/06/21
        #DirName = f'{scanID}-{seriesDesc}'.replace(' ', '_').replace('.', '_').replace('-', '_') # 01/06/21
        #
        ##ResourceDir = os.path.join(rootExportDir, expLab, 'scans', DirName,
        ##                           'resources')
        #
        #ResourceDir = os.path.join(rootExportDir, 'projects', projID, 
        #                           'subjects', subjLab, 'experiments', 
        #                           expLab, 'scans', DirName, 'resources')
        #
        #dicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        #
        #print(f"DICOMs downloaded to:\n{dicomDir}\n")
        
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID]['resources']['DICOM']['files']\
                .update({'studyUID' : studyUID,
                         'seriesDesc' : seriesDesc,
                         'seriesUID' : seriesUID,
                         'dicomDir' : dicomDir})
                
    #print(f'pathsDict = {pathsDict}\n')
    
    return pathsDict

def download_scan_Pre_210809(url, projId, subjLab, expLab, scanId,
                  xnatSession=None, rootExportDir='', pathsDict=None,
                  username=None, password=None):
    """
    Download scan from XNAT.
    
    Returns a dictionary containing paths mimacking XNAT file hierarchy and the 
    XNAT requests session.
    
    Parameters
    ----------
    url : str
        URL (e.g. 'http://10.1.1.20').
    projId : str
        The project ID of interest.
    subjLab : str
        The subject label of interest. 
    expLab : str
        The DICOM study / XNAT experiment label of interest. 
    scanId : str
        The DICOM series label / XNAT scan ID of interest.
    xnatSession : requests.models.Response, optional
        If provided a new session request will be avoided.
    rootExportDir : str, optional
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for proj_label, subjLab, scans or
        assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "xnat_downloads" within the default downloads directory.
    pathsDict : dict, optional
        Dictionary containing paths of data downloaded.
    username : str, optional
        If not provided user will be prompted to enter a user name.
    password : str, optional
        If not provided user will be prompted to enter a password.
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded. 
    xnatSession : requests.models.Response
    
    Notes
    -----
    GET - /data/projects/{project-id}/subjects/{subject-id | subject-label}/
           experiments/{experiment-id | experiment-label}/scans/{scan-id}/files
           
    GET - /data/experiments/{experiment-id}/scans/{scan-id}/files
    
    https://wiki.xnat.org/display/XAPI/Image+Scan+Resource+API
    """
    
    
    """ Set to True if issue with zip file (e.g. duplication of files): """
    OVERWRITE_ZIP = False
    
    if xnatSession == None:
        xnatSession = create_session(url, username, password)
    
    if rootExportDir == '':
        rootExportDir = os.path.join(Path.home(), "Downloads", 
                                      "xnat_downloads")
    
    uri = f'{url}/data/projects/{projId}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/scans/{scanId}/'\
          + 'resources/DICOM/files?format=zip'
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    exportDir = os.path.join(rootExportDir, 'projects', projId, 
                             'subjects', subjLab, 'experiments')

    if not os.path.isdir(exportDir):
        Path(exportDir).mkdir(parents=True)
        print(f'Created directory:\n {exportDir}\n')
    
    filepath = os.path.join(exportDir,
                            f'Experiment_{expLab}__Scan_{scanId}.zip')
    
    if not os.path.exists(filepath) or OVERWRITE_ZIP:
        with open(filepath, 'wb') as file:
            file.write(request.content)
            
        print(f'Zipped file downloaded to: \n{filepath}\n')
        
        #zip_file = ZipFile(filepath, mode='r')
        #zip_file.extractall(exportDir)
        
        with ZipFile(filepath, 'r') as file:
            file.extractall(exportDir)
        
        print(f'Zipped file extracted to: \n{exportDir}\n')
    else:
        if not OVERWRITE_ZIP:
            print(f'Zipped file already downloaded to: \n{filepath}\n')
    
    
    """ Get the directory name of the exported DICOMs: """
    # Get the filepath of the first file in the zip file:
    with ZipFile(filepath, 'r') as file:
        filepath = file.infolist()[0].filename
    
    full_filepath = os.path.join(exportDir, filepath)
    
    #print(f'full_filepath = {full_filepath}\n')
    
    dicomDir = os.path.split(full_filepath)[0]
    
    print(f"DICOMs downloaded to:\n{dicomDir}\n")
    
    """ Store info in a dictionary: """
    pathsDict, keys = create_pathsDict_for_scan(projId, subjLab, 
                                                expLab, scanId,
                                                pathsDict)
    
    if not 'Dir' in keys:
        """ Get the experiment of interest: """
        uri = f'{url}/data/projects/{projId}/subjects/{subjLab}/'\
              + f'experiments/{expLab}?format=json'
        
        request = xnatSession.get(uri)
        
        # Raise status error if not None:
        if request.raise_for_status() != None:
            print(request.raise_for_status())
            
        experiment = request.json()
        
        #print(experiment['items'], '\n')
        #print(experiment['items'][0], '\n')
        #print(experiment['items'][0]['children'], '\n')
        #print(experiment['items'][0]['children'][1], '\n')
        #print(experiment['items'][0]['children'][1]['items'], '\n')
        #print(experiment['items'][0]['children'][1]['items'][0], '\n')
        #print(experiment['items'][0]['children'][1]['items'][0]['data_fields'], '\n')

        #study_uid = experiment['items'][0]['children'][1]['items'][0]['data_fields']['type']['UID']
        #study_uid = experiment['items'][0]['children'][1]['items'][0]['data_fields']['UID'] # 27/05/2021
        study_uid = experiment['items'][0]['data_fields']['UID'] # 27/05/2021
        
        """ There may be more than one item in the list 
        experiment['items'][0]['children'] (e.g. depending on whether there are 
        image assessors or not). Get the index that corresponds to scans: """
        scanInd = None # initial value
        
        for i in range(len(experiment['items'][0]['children'])):
            if 'scan' in experiment['items'][0]['children'][i]['field']:
                scanInd = i
        
        if scanInd == None:
            msg = f"There are no scans for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
            raise Exception(msg)
            
        
        """ Get the list of series labels (scan IDs), seriesUids and series 
        descriptions for this series/scan: """
        scanIds = []
        seriesUids = []
        seriesDescs = []
        
        #for exp in experiment['items'][0]['children'][1]['items']: # 27/05/21
        for exp in experiment['items'][0]['children'][scanInd]['items']: # 27/05/21
            scanIds.append(exp['data_fields']['ID'])
            seriesUids.append(exp['data_fields']['UID'])
            seriesDescs.append(exp['data_fields']['type'])
            
        
        seriesUid = seriesUids[scanIds.index(scanId)]
        seriesDesc = seriesDescs[scanIds.index(scanId)]
        
        """ Replacing code below with something more reliable (i.e. rather than
        replacing spaces, dashes and full stops with underscores) (01/06/21).
        """
            
        #""" The folder name between the 'scans' folder and 'resources' folder
        #has the form 'scanId-seriesDesc', but it appears that spaces, dashes and 
        #full stops within seriesDesc become underscores in the folder structure 
        #of the extracted zip file. """
        ##print(f'seriesUid = {seriesUid}')
        ##print(f'seriesDesc = {seriesDesc}\n')
        #
        ##DirName = f'{scanId}-{seriesDesc}'.replace(' ', '_').replace('.', '_') # 01/06/21
        #DirName = f'{scanId}-{seriesDesc}'.replace(' ', '_').replace('.', '_').replace('-', '_') # 01/06/21
        #
        ##ResourceDir = os.path.join(rootExportDir, expLab, 'scans', DirName,
        ##                           'resources')
        #
        #ResourceDir = os.path.join(rootExportDir, 'projects', projId, 
        #                           'subjects', subjLab, 'experiments', 
        #                           expLab, 'scans', DirName, 'resources')
        #
        #dicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        #
        #print(f"DICOMs downloaded to:\n{dicomDir}\n")
        
        pathsDict['projects'][projId]['subjects'][subjLab]['experiments']\
        [expLab]['scans'][scanId]['resources']['DICOM']['files']\
        .update({'study_uid' : study_uid,
                 'seriesDesc' : seriesDesc,
                 'seriesUid' : seriesUid,
                 'dicomDir' : dicomDir})
    
    return pathsDict, xnatSession



""" From xnat_tools.sessions.py """

def create_session_OLD(url, username=None, password=None):
    """Create a requests session.
    
    Parameters
    ----------
    url : str
        URL (e.g. 'http://10.1.1.20').
    username : str, optional
        If not provided user will be prompted to enter a user name.
    password : str, optional
        If not provided user will be prompted to enter a password.
    
    Returns
    -------
    xnatSession : requests.models.Response
    """
    
    from getpass import getpass
    import requests
    #from pathlib import Path
    
    if username == None:
        username = getpass("Enter user name: ")
    
    if password == None:
        password = getpass("Enter password:")
        
    with requests.Session() as xnatSession:
        xnatSession.url = url
        xnatSession.auth = (f'{username}', f'{password}')
    
    # Test connection:
    request = requests.get(xnatSession.url, auth=xnatSession.auth)
    
    # Raise status error if not None:
    if request.raise_for_status() == None:
        print(f'Connection established to XNAT at {url}\n')
    else:
        print(request.raise_for_status())
    
    return xnatSession

