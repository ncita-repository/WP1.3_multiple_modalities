# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:31:59 2021

@author: ctorti
"""



#from importlib import reload
import time
from copy import deepcopy
from pydicom import dcmread
import SimpleITK as sitk

from general_tools.general import get_unique_items, are_items_equal_to_within_eps
from general_tools.geometry import prop_of_segs_in_extent
from general_tools.shifting import replace_ind_in_C2SindsByRoi
from general_tools.shifting import shift_frames_in_pixarrBySeg
from general_tools.shifting import get_src_C2SindsByRoi_to_trg_C2SindsByRoi
from general_tools.pixarr_ops import mean_frame_in_pixarrBySeg
from general_tools.pixarr_ops import or_frame_of_pixarrBySeg
from inputs_checker import are_inputs_valid
from conversion_tools.pixarrs_ims import pixarrByRoi_to_imByRoi
from dicom_tools.metadata import get_roicol_labels, get_dcm_uids
from image_tools.imports import import_im
from image_tools.resampling import resample_im, resample_labimByRoi
from image_tools.registering import register_ims
from seg_tools.metadata import get_PFFGS_to_slice_inds_by_seg
from seg_tools.get_seg_data_of_interest import get_seg_data_of_interest
from seg_tools.create import create_seg
from dro_tools.dro_to_tx import create_tx_from_spa_dro
from dro_tools.dro_to_tx import create_pre_tx_from_def_dro
from dro_tools.dro_to_tx import create_tx_from_def_dro
from dro_tools.tx_to_spa_dro import create_spa_dro_from_tx
from dro_tools.tx_to_def_dro import create_def_dro_from_tx
from plotting_tools.plotting import plot_pixarrBySeg

#from DicomTools import GetRoiLabels, GetDicomSOPuids
#from DroTools import CreateSpaDro, CreateDefDro

def copy_seg(params, ):
    
    
    """
    Check prop_of_segs_in_extent
    """


def copy_seg_OLD(SrcSegFpath, SrcSliceNum, SrcSegLabel, SrcDcmDir, TrgDcmDir,
             UseCaseToApply, TrgSegFpath=None, TrgSegLabel=None, 
             TrgSliceNum=None, ResInterp='BlurThenLinear', 
             ApplyPreResBlur=False, PreResVariance=(1,1,1),
             ApplyPostResBlur=False, PostResVariance=(1,1,1),
             ForceReg=False, UseDroForTx=True, Transform='affine',
             MaxIters='512', InitMethod='landmarks', 
             SrcFidsFpath='', TrgFidsFpath='', TxInterp='NearestNeighbor', 
             ApplyPreTxBlur=False, PreTxVariance=(1,1,1),
             ApplyPostTxBlur=True, PostTxVariance=(1,1,1),
             Dro=None, TxtToAddToSegLabel='', #OverwriteDro=False, 
             LogToConsole=False, 
             DictOfInputs={}, ListOfInputs=[], ListOfTimings=[]):
    """
    Copy a DICOM-SEG from Source to Target dataset.
    
    This function will, depending on the inputs, copy either:
        1. A specific segmentation
        2. All segmentations in a specific segment
        3. All segmentations in all segments
     
    Parameters
    ----------
    SrcSegFpath : str
        Filepath of the Source SEG file.
    SrcSliceNum : int or None
        The slice indeces within the Source DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a  
        single segmentation). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
    SrcSegLabel : str
        All or part of the Source SegmentLabel for the segment containing the
        segmentation(s) to be copied.
    SrcDcmDir : str
        Directory containing the Source DICOMs.
    TrgDcmDir : str
        Directory containing the Target DICOMs.
    UseCaseToApply : str
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
    TrgRtsFpath : str or None, optional (None by default)
        Filepath to the Target SEG file that the segmentation(s) is/are to be 
        copied to. TrgSegFpath != None only if an existing Target SEG exists 
        and is to added to. An existing SEG will be added to only for Direct 
        copy operations.     
    TrgSegLabel : str or None, optional (None by default) 
        All or part of the SegmentLabel of the destination segment.
    TrgSliceNum : int, optional (None by default)
        The slice index within the Target DICOM stack where the segmentation 
        will be copied to (applies only for the case of direct copies of single 
        segmentations). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the segmentation(s) will be copied to 
        will not depend on user input.
    ResInterp : str, optional ('BlurThenLinear' by default)
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    PreResVariance : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    ApplyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    PostResVariance : tuple of floats, optional ((1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
    ForceReg : bool, optional (False by default)
        If True the Source image will be registered to the Target image, even 
        if the Source and Target images have the same FrameOfReferenceUID.  
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
    PostTxVariance : tuple of floats, optional ((1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
    TxtToAddToSegLabel : str, optional ('' by default)
        String of text to add to the segment label of the new Target SEG.
    Dro : Pydicom Object, optional (None by default)
        A Spatial or Deformable Spatial Registration Object that contains the
        parameters required to register the moving (Source) image to the fixed
        (Target) image.
    OverwriteDro : bool, optional (False by default)
        If True, Dro (if not None) will be overwritten.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    DictOfInputs : dict, optional ({} by default)
        A dictionary containing the inputs that were called to CopyRoi().
    ListOfInputs : list of str, optional ([] by default)
        A list containing the inputs that were called to CopyRoi().
    ListOfTimings : list of str, optional ([] by default)
        A list of the time to execute certain tasks during the calling of 
        CopySeg().
          
    Returns
    -------
    TrgSeg : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target SEG object.
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgSeg, None otherwise.
    DictOfInputs : dict
        A dictionary containing the inputs that were called to CopyRoi().
    ListOfInputs : list of str
        A list containing the inputs that were called to CopyRoi().
    ListOfTimings : list of str
        A list (for each message) of timings for individual operations.
        
    Notes
    -----
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target segmentations will appear displaced  
    w.r.t. the anatomical features highlighted in the Source image.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
    #""" Log timing messages. """
    #TimingMsgs = []
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of copy_seg():')
        print('\n\n', '-'*120)
    
    """ Load the Source SEG (and Target SEG if applicable): """
    SrcSeg = dcmread(SrcSegFpath)
    if TrgSegFpath:
        TrgSeg = dcmread(TrgSegFpath)
    else:
        TrgSeg = None
    
    """ Check the validity of the combination of inputs. """
    are_inputs_valid(SrcSeg, SrcSegLabel, SrcSliceNum, TrgSeg, TrgSegLabel,
                     TrgSliceNum, LogToConsole)
    
    
    SrcSegLabels = get_roicol_labels(SrcSeg)
    if TrgSegFpath:
        TrgSegLabels = get_roicol_labels(TrgSeg)
        
    if LogToConsole:
        print(f'\nSrcSegLabels = {SrcSegLabels}')
        if TrgSegFpath:
            print(f'\nTrgSegLabels = {TrgSegLabels}')
    
    """ Get the data of interest from the Source SEG. """
    SrcPixArrBySeg,\
    SrcF2SindsBySeg = get_seg_data_of_interest(SrcSegFpath, SrcSliceNum, 
                                               SrcSegLabel, SrcDcmDir, 
                                               LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source SEG file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    """ Raise exception if there is no data of interest from the Source SEG."""
    if not get_unique_items(SrcF2SindsBySeg):
        Studyuid, Seriesuid, FORuid, SOPuids = get_dcm_uids(SrcDcmDir)
        
        F2SindsBySeg = get_PFFGS_to_slice_inds_by_seg(SrcSeg, SOPuids)
        
        msg = f"There are no segmentations on slice {SrcSliceNum} in any "\
              + f"segment in the Source SEG. There are {len(SrcSegLabels)} "\
              + "segments in the Source SEG:"
        
        for i in range(len(SrcSegLabels)):
            msg += f"\nSegment {i+1} with label '{SrcSegLabels[i]}' has "\
                   + f"{len(F2SindsBySeg[i])} segmentations that correspond "\
                   + f"to slices {F2SindsBySeg[i]}" 
        raise Exception(msg)
    
    """
    Determine whether the segmentations that make up the segment(s) intersect 
    with the Target image extent. 
    """
    FracProp = prop_of_segs_in_extent(PixArrBySeg=SrcPixArrBySeg,
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
    SrcIm = import_im(SrcDcmDir)
    TrgIm = import_im(TrgDcmDir)
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    if UseCaseToApply in ['1', '2a']: # 08/02
        """ Direct copy of a segmentation without resampling or registration
        transformation. """
        
        #from general_tools.shifting import replace_ind_in_C2SindsByRoi
        #from general_tools.shifting import shift_frames_in_pixarrBySeg
        
        """ The segmentation on SrcSliceNum is to be copied to TrgSliceNum. 
        Modify the frame-to-slice index (from SrcSliceNum) to TrgSliceNum. """
        SrcF2SindsBySeg\
            = replace_ind_in_C2SindsByRoi(C2SindsByRoi=SrcF2SindsBySeg,
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
            TrgPixArrBySeg, TrgF2SindsBySeg\
            = get_seg_data_of_interest(TrgSegFpath, SrcSliceNum,
                                       #SrcSegLabel, TrgDcmDir,
                                       TrgSegLabel, TrgDcmDir,
                                       LogToConsole)
            """ 05/03/21: Previously SrcSegLabel was used as an input in 
            get_seg_data_of_interest for Target, hence requiring that the Target
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
        
        """
        The segmentation will be copied to a different slice location, so the
        z-components of the indices will need to be modified. """
        
        """ Shift the in-plane (x & y) and out-of-plane (z) components of the 
        indices in SrcPixArrBySeg to account for the shift from SrcSliceNum to
        TrgSliceNum. """
        SrcPixArrBySeg, SrcF2SindsBySeg\
            = shift_frames_in_pixarrBySeg(PixArrBySeg=SrcPixArrBySeg, 
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
        
        #from general_tools.shifting import get_src_C2SindsByRoi_to_trg_C2SindsByRoi
        
        if False:
            """ TrgPixArrBySeg is simply SrcPixArrBySeg, and 
            TrgF2SindsBySeg is SrcF2SindsBySeg. """
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        TrgF2SindsBySeg\
            = get_src_C2SindsByRoi_to_trg_C2SindsByRoi(SrcF2SindsBySeg, 
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
        
        ##import numpy as np
        ##from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        ##from ImageTools import ResampleLabImByRoi#, GetImageInfo
        ##from ConversionTools import ConvertImagePixelType
        ##from ConversionTools import PixArr2Image, Image2PixArr
        #from conversion_tools.pixarrs_ims import pixarrByRoi_to_imByRoi
        ##from GeneralTools import NumOfListsAtDepthTwo
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabImBySeg = pixarrByRoi_to_imByRoi(PixArrByRoi=SrcPixArrBySeg, 
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
        
        #from image_tools.resampling import resample_labimByRoi
        
        """ Resample the source label images. """
        ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
            = resample_labimByRoi(LabImByRoi=SrcLabImBySeg, 
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
        
        #from image_tools.registering import register_ims
        #from ImageTools import TransformImWithListOfSitkTxs
        #from ImageTools import TransformLabImByRoi
        #from ImageTools import TransformLabImByRoiWithListOfSitkTxs
        #from DroTools import CreateSitkTxFromTxMatrix
        #from dro_tools.dro_to_tx import create_tx_from_spa_dro
        #from dro_tools.dro_to_tx import create_pre_tx_from_def_dro
        #from dro_tools.dro_to_tx import create_tx_from_def_dro
        #from general_tools.general import are_items_equal_to_within_eps
        
        if Dro == None or (Dro != None and UseDroForTx):
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
            RegIm, InitialTx, FinalTx\
                = register_ims(FixIm=TrgIm, MovIm=SrcIm, 
                               RegTxName=Transform, InitMethod=InitMethod,
                               FixFidsFpath=TrgFidsFpath, 
                               MovFidsFpath=SrcFidsFpath,
                               #SamplingPercentage=SamplingPercentage,
                               LogToConsole=LogToConsole)
            
            """ Get the registration transform parameters: """
            TxParams = [str(item) for item in FinalTx.GetParameters()]
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to register the Source image to the Target '\
                  + f'image using {Transform} transform.\n'
            ListOfTimings.append(msg)
            if LogToConsole == False:
                print(f'*{msg}')
        else:
            """ The transform parameters required to register SrcIm to TrgIm 
            are obtained from the transform matrix, found from the DRO (subject
            assessor) in XNAT. Create a SimpleITK transform from Dro: 
            """
            #""" 04/06/21 Workaround described above. """
            #ListOfSitkTxs = []
            
            if Transform in ['rigid', 'affine']:
                SitkTx = create_tx_from_spa_dro(Dro, LogToConsole)
                
                TxParams = list(SitkTx.GetParameters())
                
                #ListOfSitkTxs.append(SitkTx)
                
                """ 
                09/07/21:
                --> Consider returning the pre-reg tx if it's not the identity.
                PreRegTx = None or identity tx?
                """
            else:
                # Create the deformable SimpleITK Transform:
                SitkTx = create_tx_from_def_dro(Dro, RefIm=TrgIm,
                                                LogToConsole=LogToConsole)
                
                # Create the pre-deformable SimpleITK Transform:
                PreRegTx = create_pre_tx_from_def_dro(Dro, LogToConsole)
                
                TxParams = list(PreRegTx.GetParameters())
                
                Identity = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
                
                if are_items_equal_to_within_eps(TxParams, Identity):
                    """ 
                    The pre-registration matrix is simply the identity matrix,
                    so no need to resample/transform the label image using
                    PreRegTx. 
                    """
                    #SitkTx = create_tx_from_def_dro(Dro, RefIm=TrgIm,
                    #                                LogToConsole=LogToConsole)
                    
                    RegIm = resample_im(Im=SrcIm, RefIm=TrgIm, sitkTx=SitkTx, 
                                        Interp='Linear', LogToConsole=False)
                else:
                    """ 
                    The pre-registration matrix is not the identity matrix so
                    need to resample/transform the label image using PreRegTx. 
                    """
                    #ListOfSitkTxs.append(PreRegTx)
                    
                    #SitkTx = create_tx_from_def_dro(Dro, RefIm=TrgIm,
                    #                                LogToConsole=LogToConsole)
                    
                    # Compose transforms: (09/07/21)
                    comp_tx = sitk.CompositeTransform(PreRegTx)
                    comp_tx.AddTransform(SitkTx)
                    
                    RegIm = resample_im(Im=SrcIm, RefIm=TrgIm, sitkTx=comp_tx, 
                                        Interp='Linear', LogToConsole=False)
                
                #ListOfSitkTxs.append(SitkTx)
            
            #""" Since image registration was skipped, create RegIm by applying
            #ListOfSitkTxs iteratively: """
            #RegIm = TransformImWithListOfSitkTxs(MovIm=SrcIm, FixIm=TrgIm, 
            #                                     ListOfSitkTxs=ListOfSitkTxs, 
            #                                     Interp='Linear', 
            #                                     LogToConsole=LogToConsole)
            
            #RegIm = resample_im(Im=SrcIm, RefIm=TrgIm, sitkTx=SitkTx, 
            #                    Interp='Linear', LogToConsole=False)
        
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
        
        #ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
        #= TransformLabImByRoiWithListOfSitkTxs(LabImByRoi=SrcLabImBySeg,
        #                                       F2SindsByRoi=SrcF2SindsBySeg,
        #                                       TrgIm=TrgIm,
        #                                       ListOfSitkTxs=ListOfSitkTxs,
        #                                       Interp=TxInterp,
        #                                       ApplyPostTxBlur=ApplyPostTxBlur,
        #                                       PostTxVariance=PostTxVariance,
        #                                       ApplyPostTxBin=ApplyPostTxBin,
        #                                       LogToConsole=LogToConsole)
        
        ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
            = resample_labimByRoi(LabImByRoi=SrcLabImBySeg,
                                  F2SindsByRoi=SrcF2SindsBySeg,
                                  SrcIm=SrcIm, TrgIm=TrgIm,
                                  sitkTx=SitkTx, Interp=TxInterp,
                                  ApplyPreResBlur=ApplyPreTxBlur,
                                  PreResVariance=PreResVariance,
                                  ApplyPostResBlur=ApplyPostTxBlur,
                                  PostResVariance=PostTxVariance,
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
        #SelxImFiltOrSitkTx = None
        SitkTx = None
        TxParams = None
        #ListOfSitkTxs = None
    
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
        #from general_tools.shifting import shift_frames_in_pixarrBySeg
        #from PlottingTools import PlotPixArrBySeg
        
        R = len(ResSrcPixArrBySeg) # number of segments
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
        MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
        
        if MaxNumOfFramesBySeg > 1:
            """ Reduce the number of frames in each pixel array to 1. """
            
            # Perform averaging or OR operation?
            #ReduceUsing = 'mean'
            ReduceUsing = 'OR'
            
            if LogToConsole:
                print(f'\nThere are {NumOfFramesBySeg} frames in each',
                      f'ROI/segment. Using {ReduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                plot_pixarrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                 f'Prior to {ReduceUsing} operation')
                
            if ReduceUsing == 'mean':
                #from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5 # <-- is this too high?
    
                ResSrcPixArrBySeg, ResF2SindsBySeg\
                    = mean_frame_in_pixarrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                                F2SindsBySeg=ResSrcF2SindsBySeg,
                                                MakeBinary=True, 
                                                BinaryThresh=BinaryThresh,
                                                LogToConsole=LogToConsole)
            else:
                #from GeneralTools import OrFrameOfPixArrBySeg
                
                ResSrcPixArrBySeg, ResSrcF2SindsBySeg\
                    = or_frame_of_pixarrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                              F2SindsBySeg=ResSrcF2SindsBySeg,
                                              LogToConsole=LogToConsole)
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment, and',
                      f'the new ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}.')
                
                plot_pixarrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
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
            = shift_frames_in_pixarrBySeg(PixArrBySeg=ResSrcPixArrBySeg, 
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
            print('After running shift_frames_in_pixarrBySeg():')
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
    TrgSeg = create_seg(SrcSegFpath, TrgSegFpath, TrgSegLabel, TrgPixArrBySeg, 
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
        #from DroTools import CreateSpaDroFromSitkTx, CreateDefDroFromBsplineTx
        
        # ContentDescription:
        """ Consider adding this as an input variable for the user to specify.
        """
        ContentDesc = ''
        
        if Transform in ['rigid', 'affine']:
            """ Create a Spatial Registration DRO. """
            #Dro = CreateSpaDro(SrcDcmDir, TrgDcmDir, TxParams, ContentDesc, 
            #                   LogToConsole)
            Dro = create_spa_dro_from_tx(SrcDcmDir, TrgDcmDir, SitkTx,
                                         ContentDesc, LogToConsole)
            """
            09/07/21:
                --> Consider adding PreRegTx (but need to return it from
                    create_spa_dro_from_tx..)
            """
        else:
            """ Create a Deformable Spatial Registration DRO. """
            #Dro = CreateDefDro(SrcDcmDir, TrgDcmDir, GridOrig, GridDir, 
            #                   GridDims, GridRes, VectGridData, TxParams,
            #                   ContentDesc, LogToConsole)
            #Dro = CreateDefDroFromBsplineTx(SrcDcmDir, TrgDcmDir, 
            #                                BsplineTx=SelxImFiltOrSitkTx, 
            #                                PreRegTx=LandmarkTx, 
            #                                Description=ContentDesc, 
            #                                LogToConsole=LogToConsole)
            """ 
            Note:
                The first Transform in ListOfSitkTxs is the landmark-based
                Transform, the second/last one is the final BSpline Transform.
            """
            #Dro = CreateDefDroFromBsplineTx(SrcDcmDir, TrgDcmDir, 
            #                                BsplineTx=ListOfSitkTxs[-1], 
            #                                PreRegTx=ListOfSitkTxs[0], 
            #                                Description=ContentDesc, 
            #                                LogToConsole=LogToConsole)
            
            """ 
            09/07/21:
                Since work-around was no longer required there may be no need
                to store a pre-registration transform (PreRegTx), so omit that
                parameter (argument defaults to None):
            """
            Dro = create_def_dro_from_tx(SrcDcmDir, TrgDcmDir, 
                                         BsplineTx=SitkTx,
                                         PreRegTx=PreRegTx,
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
           SrcPixArrBySeg, SrcLabImBySeg, SrcF2SindsBySeg,\
           RegIm, SitkTx, TxParams,\
           ResSrcLabImBySeg, ResSrcPixArrBySeg, ResSrcF2SindsBySeg,\
           TrgSeg, Dro, DictOfInputs, ListOfInputs, ListOfTimings