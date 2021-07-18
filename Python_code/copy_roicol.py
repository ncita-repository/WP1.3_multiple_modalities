# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 15:24:22 2021

@author: ctorti
"""

import time
from datetime import datetime
import os
#from pypref import Preferences
#from pydicom import dcmread
#from importlib import reload

from RtsTools import CopyRts
from SegTools import CopySeg

from io_tools.imports import import_dict_from_json
from io_tools.inputs_checker import which_use_case
from io_tools.exports import export_dict_to_json, export_list_to_txt
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
#from xnat_tools.scan_assessors import get_fname_and_id
from xnat_tools.format_paths_dict import get_scan_asr_fname_and_id
from xnat_tools.dro import search_dro

def copy_roicol(cfgDir, xnatSession=None):
    """
                XnatUrl, xnatSession, ProjId, SubjLabel, 
                SrcExpLabel, SrcScanId, SrcSliceNum,
                SrcRoiColMod, SrcRoiColName, SrcRoiName, 
                TrgExpLabel, TrgScanId, TrgSliceNum=None, 
                TrgRoiColMod=None, TrgRoiColName=None, TrgRoiName=None, 
                TxtToAddToTrgRoiColName='',
                pathsDict=None, XnatDownloadDir='default', 
                LogToConsole=False, ExportLogFiles=False,
                ResInterp='BlurThenLinear', PreResVar=(1,1,1), 
                ApplyPostResBlur=False, PostResVar=(1,1,1),
                ForceReg=False, UseDroForTx=True, 
                Transform='affine', MaxIters='512', 
                InitMethod='centerofgravity',
                SrcFidsFpath=None, TrgFidsFpath=None,
                TxInterp='NearestNeighbor', ApplyPostTxBin=True, 
                ApplyPostTxBlur=True, PostTxVar=(1,1,1), DictOfInputs=None):
    """
    """
    Copy a contour/ROI/RTSTRUCT/segmentation/segment/SEG from data downloaded
    from XNAT.
    
    Parameters
    ----------
    XnatUrl : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    xnatSession : requests session or None
    ProjId : str
        The project ID of interest.
    SubjLabel : str
        The subject label of interest. 
    SrcExpLabel : str
        The Source DICOM study / XNAT experiment label. 
    SrcScanId : str
        The Source DICOM series label / XNAT scan ID.
    SrcSliceNum : int
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied. The integer is zero-indexed.
    SrcRoiColMod : str
        The modality of the Source ROI Collection.
    SrcRoiColName : str
        The Source ROI Collection name (StructureSetLabel or SeriesDescription).
    SrcRoiName : str
        The Source ROIName or SegmentLabel.
    TrgExpLabel : str
        The Target DICOM study / XNAT experiment label. 
    TrgScanId : str
        The Target DICOM series label / XNAT scan ID.
    TrgSliceNum : int, optional unless making a direct copy (None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to.  This only applies for direct copies, 
        hence the default value None. The integer is zero-indexed.
    TrgRoiColMod : str, optional unless TrgRoiColName != None (None by default)
        The Target image assessor modality.
    TrgRoiColName : str, optional unless TrgRoiColMod != None (None by default)
        The Target ROI Collection name (StructureSetLabel or SeriesDescription). 
        If provided and if a direct copy is to be made, existing contours/
        segmentations for the ROI/segment will be preserved.
    TrgRoiName : str, optional (None by default)
        The Target ROIName or SegmentLabel.
    TxtToAddToTrgRoiColName : str, optional ('' by default)
        If provided the string of text will be appended to the ROI Collection 
        name (StructureSetLabel or SeriesDescription), and to the filename of 
        the exported (new) Target RTS/SEG ROI Collection.
    pathsDict : dict, optional (None by default)
        Dictionary containing paths of data downloaded. If provided, paths of
        newly downloaded data will be added to pathsDict. If not provided, a
        new dictionary will be initialised.
    XnatDownloadDir : str, optional ('default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjectLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    LogToConsole : bool, optional (False by default)
        If True, some results will be logged to the console.
    ExportLogFiles : bool, optional (False by default)
        If True, log files will be exported.
    ResInterp : str, optional ('BlurThenLinear' by default)
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding)
    PreResVar : tuple of floats, optional ((1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    ApplyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    PostResVar : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
    ForceReg : bool, optional (False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    UseDroForTx : bool, optional (True by default)
        If True the transformation parameters from any suitable DRO (for the
        required Transform) will be used instead of performing image
        registration.
    Transform : str, optional ('affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    MaxIters : str, optional ('512' by default)
        The maximum number of iterations used for the optimiser during image 
        registration (if applicable).
    InitMethod : str, optional ('centerofgravity' by default)
        The initialisation method used for initial alignment prior to image
        registration (if applicable).
    SrcFidsFpath : str, optional ('moving_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Source/Moving image.
    TrgFidsFpath : str, optional ('fixed_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Target/Fixed image.
    TxInterp : str, optional ('NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    ApplyPostTxBin : bool, optional (True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
    ApplyPostTxBlur : bool, optional (True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
    PostTxVar : tuple of floats, optional ((1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
    
    Returns
    -------
    TrgRoiCol : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing RTS/SEG) 
        Target ROI Collection (RTS/SEG) object.
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgSeg, None otherwise.
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    xnatSession : requests session
    DictOfInputs : dict
        A dictionary containing the inputs that were called to CopyRoi().
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
    TimingMsgs : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopyRoi() and subsequent tasks called within CopyRts() or CopySeg().
    
    Notes
    -----
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours/segmentations will appear 
    displaced w.r.t. the anatomical features highlighted in the Source image.
    """
    
    #runDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    currentTime = datetime.now()
    runDateTime = currentTime.strftime('%Y%m%d_%H%M%S.%f')
    
    timings = [f'runDateTime: {runDateTime}\n']
    
    times = []
    
    
    """ Read in the config settings: """
    fpath = os.path.join(cfgDir, 'general_params.json')
    genParams = import_dict_from_json(filepath=fpath)
    
    fpath = os.path.join(cfgDir, 'xnat_config.json')
    xnatCfg = import_dict_from_json(filepath=fpath)
    
    fpath = os.path.join(cfgDir, 'xnat_src_params.json')
    srcXnatParams = import_dict_from_json(filepath=fpath)
    
    fpath = os.path.join(cfgDir, 'xnat_trg_params.json')
    trgXnatParams = import_dict_from_json(filepath=fpath)
    
    
    """ Create XNAT session. """
    if xnatSession == None:
        xnatSession = create_session(xnatCfg)
    
    # Start timing
    times.append(time.time())
    
    """ Download data from XNAT. """
    pathsDict = download_scan(xnatCfg, srcXnatParams, genParams, xnatSession)
    
    pathsDict = download_scan(xnatCfg, trgXnatParams, genParams, xnatSession,
                              pathsDict)

    pathsDict = download_im_asr(xnatCfg, srcXnatParams, genParams, xnatSession,
                                pathsDict)
    
    trgRoicolName = trgXnatParams['roicolName']
    
    if trgRoicolName != None:
        pathsDict = download_im_asr(xnatCfg, trgXnatParams, genParams, 
                                    xnatSession, pathsDict)
        
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to download data from XNAT.\n'
        timings.append(msg)
        print(f'*{msg}')
    
    #print('\n', pathsDict, '\n')
    
    """ Get the DICOM directories and ROI Collection file paths: """
    srcDcmDir = pathsDict['projects'][projId]['subjects'][subjLab]\
                ['experiments'][srcExpLab]['scans'][srcScanId]\
                ['resources']['DICOM']['files']['Dir']
    
    trgDcmDir = pathsDict['projects'][projId]['subjects'][subjLab]\
                ['experiments'][trgExpLab]['scans'][trgScanId]\
                ['resources']['DICOM']['files']['Dir']
    
    srcRoicolFname, srcAsrId\
        = get_fname_and_id(pathsDict, projId, subjLab, srcExpLab,
                           srcRoicolMod, srcRoicolName)
    
    srcRoicolFpath = pathsDict['projects'][projId]['subjects'][subjLab]\
                     ['experiments'][srcExpLab]['assessors'][srcAsrId]\
                     ['resources'][srcRoicolMod]['files'][srcRoicolFname]\
                     ['Fpath']
    
    if trgRoicolName == None:
        trgRoicolFpath = None
    else:
        trgRoicolFname, trgAsrId\
            = get_fname_and_id(pathsDict, projId, subjLab, trgExpLab,
                               trgRoicolMod, trgRoicolName)
    
        trgRoicolFpath = pathsDict['projects'][projId]['subjects'][subjLab]\
                         ['experiments'][trgExpLab]['assessors'][trgAsrId]\
                         ['resources'][trgRoicolMod]['files'][trgRoicolFname]\
                         ['Fpath']
    
    
    print(f'srcDcmDir = {srcDcmDir}\n')
    print(f'trgDcmDir = {trgDcmDir}\n')
    
    #""" Establish whether the inputs are valid. """
    #CheckValidityOfInputs(SrcRoiName, srcSlcNum, TrgRoiName, trgSlcNum, 
    #                      LogToConsole)
    
    """ Determine which use case to apply. """
    useCaseThatApplies, useCaseToApply\
        = which_use_case(srcSlcNum, trgSlcNum, srcDcmDir, trgDcmDir, ForceReg,
                         LogToConsole=True)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        #print(f'\n\nDone.  Took {Dtime} s to run.')
        msg = f'Took {Dtime} s to determine which UseCase applies.\n'
        timings.append(msg)
        print(f'*{msg}')
        
    
    DictOfInputs['useCaseThatApplies'] = useCaseThatApplies
    DictOfInputs['useCaseToApply'] = useCaseToApply
    
    
    """ If image registration will be performed (UseCase = '5a' or '5b'), check
    whether there exists a DRO for this set of Source and Target images: """
    if useCaseToApply in ['5a', '5b']:
        times.append(time.time())
        
        Dro, TxMatrix, GridDims, GridRes,\
        VectGridData = search_dro(XnatUrl=XnatUrl, projId=projId, 
                                  subjLab=subjLab, 
                                  srcExpLab=srcExpLab, srcScanId=srcScanId,
                                  trgExpLab=trgExpLab, trgScanId=trgScanId,
                                  Transform=Transform, xnatSession=xnatSession, 
                                  LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to search XNAT for an existing DRO.\n'
        timings.append(msg)
        print(f'*{msg}')
        
        #print(f"type(Dro) after SearchXnatForDro() = {type(Dro)}\n")
    else:
        Dro = None
        #TxMatrix = None
        #GridDims = None
        #GridRes = None
        #VectGridData = None
        
    #DictOfInputs['TxMatrix'] = TxMatrix
    #DictOfInputs['GridDims'] = GridDims
    #DictOfInputs['GridRes'] = GridRes
    #DictOfInputs['VectGridData'] = VectGridData
    
    #print(type(TxMatrix))
    #print(f'TxMatrix = {TxMatrix}\n')
    
    
    
    if srcRoicolMod == 'RTSTRUCT':
        #print(f"type(Dro) before CopyRts() = {type(Dro)}\n")
        
        srcDcmDir, trgDcmDir, SrcIm, TrgIm, SrcRoiCol,\
        SrcPtsByCntByRoi, SrcC2SindsByRoi, TrgPtsByCntByRoi, TrgC2SindsByRoi,\
        SrcPixArrByRoi, SrcLabImByRoi, SrcF2SindsByRoi,\
        RegIm, ListOfSitkTxs, SelxImFiltOrSitkTx, TxParams,\
        ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi,\
        ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi, ResSrcC2SindsByRoi,\
        TrgRoiCol, Dro, DictOfInputs, ListOfInputs, timings\
        = CopyRts(SrcRtsFpath=srcRoicolFpath, srcSlcNum=srcSlcNum, 
                  SrcRoiName=SrcRoiName, srcDcmDir=srcDcmDir, trgDcmDir=trgDcmDir, 
                  useCaseToApply=useCaseToApply, TrgRtsFpath=trgRoicolFpath, 
                  TrgRoiName=TrgRoiName, trgSlcNum=trgSlcNum, 
                  ResInterp=ResInterp, PreResVariance=PreResVar, 
                  ApplyPostResBlur=ApplyPostResBlur, PostResVariance=PostResVar, 
                  ForceReg=ForceReg, UseDroForTx=UseDroForTx, 
                  Transform=Transform, MaxIters=MaxIters, InitMethod=InitMethod,
                  SrcFidsFpath=SrcFidsFpath, TrgFidsFpath=TrgFidsFpath, 
                  TxInterp=TxInterp, ApplyPostTxBin=ApplyPostTxBin, 
                  ApplyPostTxBlur=ApplyPostTxBlur, PostTxVariance=PostTxVar, 
                  #TxMatrix, GridDims, GridRes, VectGridData, 
                  Dro=Dro, TxtToAddToTrgRoiName=TxtToAddToTrgRoiColName, 
                  LogToConsole=LogToConsole,
                  DictOfInputs=DictOfInputs, #ListOfInputs=ListOfInputs, 
                  timings=timings)
        
        SrcPixArrBySeg = None
        SrcF2SindsBySeg = None
        SrcLabImBySeg = None
        TrgPixArrBySeg = None
        TrgF2SindsBySeg = None
        ResSrcLabImBySeg = None
        ResSrcPixArrBySeg = None
        ResSrcF2SindsBySeg = None
        
        #print(f"type(Dro) after CopyRts() = {type(Dro)}\n")
    
    elif srcRoicolMod == 'SEG':
        srcDcmDir, trgDcmDir, SrcIm, TrgIm, SrcRoiCol,\
        SrcPixArrBySeg, SrcF2SindsBySeg, SrcLabImBySeg,\
        TrgPixArrBySeg, TrgF2SindsBySeg,\
        RegIm, ListOfSitkTxs, SelxImFiltOrSitkTx, TxParams,\
        ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg,\
        TrgRoiCol, Dro, DictOfInputs, ListOfInputs, timings\
        = CopySeg(SrcSegFpath=srcRoicolFpath, srcSlcNum=srcSlcNum, 
                  SrcSegLabel=SrcRoiName, srcDcmDir=srcDcmDir, trgDcmDir=trgDcmDir, 
                  useCaseToApply=useCaseToApply, TrgSegFpath=trgRoicolFpath, 
                  TrgSegLabel=TrgRoiName, trgSlcNum=trgSlcNum, 
                  ResInterp=ResInterp, PreResVariance=PreResVar, 
                  ApplyPostResBlur=ApplyPostResBlur, PostResVariance=PostResVar, 
                  ForceReg=ForceReg, UseDroForTx=UseDroForTx,
                  Transform=Transform, MaxIters=MaxIters, InitMethod=InitMethod,
                  SrcFidsFpath=SrcFidsFpath, TrgFidsFpath=TrgFidsFpath, 
                  TxInterp=TxInterp, ApplyPostTxBin=ApplyPostTxBin, 
                  ApplyPostTxBlur=ApplyPostTxBlur, PostTxVariance=PostTxVar, 
                  #TxMatrix, GridDims, GridRes, VectGridData, 
                  Dro=Dro, TxtToAddToTrgRoiName=TxtToAddToTrgRoiColName, 
                  LogToConsole=LogToConsole,
                  DictOfInputs=DictOfInputs, #ListOfInputs=ListOfInputs, 
                  timings=timings)
        
        SrcPtsByCntByRoi = None
        SrcC2SindsByRoi = None
        TrgPtsByCntByRoi = None
        TrgC2SindsByRoi = None
        SrcPixArrByRoi = None
        SrcLabImByRoi = None
        SrcF2SindsByRoi = None
        ResSrcLabImByRoi = None
        ResSrcPixArrByRoi = None
        ResSrcF2SindsByRoi = None
        
    else:
        msg = f'The modality of the Source ROI Collection ({srcRoicolMod}) '\
              + 'must be "RTSTRUCT" or "SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to copy the ROI(s) and create the '\
              + f'{srcRoicolMod} object.\n'  
        timings.append(msg)
        print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export DictOfInputs as a json, and ListOfInputs as a txt, to a sub-
        directory "logs" in the current working directory. """
        JsonFname = runDateTime + '_DictOfInputs.json'
        TxtFname = runDateTime + '_ListOfInputs.txt'
        
        CWD = os.getcwd()
        
        ExportDir = os.path.join(CWD, 'logs')
        
        JsonFpath = os.path.join(ExportDir, JsonFname)
        TxtFpath = os.path.join(ExportDir, TxtFname)
        
        export_dict_to_json(DictOfInputs, JsonFname, ExportDir)
        export_list_to_txt(ListOfInputs, TxtFname, ExportDir)
            
        print('Dictionary of inputs to CopyRoi() saved to:\n', JsonFpath, 
              '\nand list of inputs saved to:\n', TxtFpath, '\n')
        
    #return TrgRoiCol, Dro, pathsDict, xnatSession, DictOfInputs,\
    #       ListOfInputs, TimingMsgs
    return srcDcmDir, trgDcmDir, SrcIm, TrgIm, SrcRoiCol,\
           SrcPtsByCntByRoi, SrcC2SindsByRoi, TrgPtsByCntByRoi, TrgC2SindsByRoi,\
           SrcPixArrByRoi, SrcLabImByRoi, SrcF2SindsByRoi,\
           SrcPixArrBySeg, SrcF2SindsBySeg, SrcLabImBySeg,\
           TrgPixArrBySeg, TrgF2SindsBySeg,\
           RegIm, ListOfSitkTxs, SelxImFiltOrSitkTx, TxParams,\
           ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi,\
           ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg,\
           ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi, ResSrcC2SindsByRoi,\
           TrgRoiCol, Dro, pathsDict, xnatSession, DictOfInputs, timings