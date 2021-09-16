# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 17:05:12 2021

@author: ctorti
"""

from importlib import reload

import general_tools.shifting
reload(general_tools.shifting)
import general_tools.pixarr_ops
reload(general_tools.pixarr_ops)
import general_tools.general
reload(general_tools.general)
import image_tools.resampling
reload(image_tools.resampling)
import image_tools.registering
reload(image_tools.registering)
import dro_tools.create_tx_from_dro
reload(dro_tools.create_tx_from_dro)
import plotting_tools.general
reload(plotting_tools.general)
import general_tools.geometry
reload(general_tools.geometry)
import general_tools.console_printing
reload(general_tools.console_printing)

import time
import os
from pathlib import Path
#from copy import deepcopy
import SimpleITK as sitk

from general_tools.shifting import (
    get_voxel_shift_bt_slices, shift_frame, replace_ind_in_C2SindsByRoi,
    z_shift_ptsByCntByRoi, shift_ptsByCntByRoi
    )
from general_tools.general import (
    get_unique_items, are_items_equal_to_within_eps, generate_reg_fname
    #does_instance_variable_exist
    )
from conversion_tools.pixarrs_ims import pixarr_to_im
from image_tools.attrs_info import get_im_info
from image_tools.resampling import resample_im, resample_labimBySeg
from image_tools.registering import register_im
from dro_tools.create_tx_from_dro import (
    create_tx_from_spa_dro, create_pre_tx_from_def_dro, create_tx_from_def_dro
    )
from general_tools.pixarr_ops import (
    mean_frame_in_pixarrBySeg, or_frame_of_pixarrBySeg
    )
from plotting_tools.general import plot_pixarrBySeg
from io_tools.exports import export_list_to_txt
from general_tools.geometry import (
    prop_of_segs_in_extent, prop_of_rois_in_extent
    )
from general_tools.console_printing import print_indsByRoi, print_ptsByCntByRoi
from plotting_tools.res_reg_results import compare_res_results

class Propagator:
    # TODO modify the docstrings
    """
    This class copies or propagates DataImporter Objects based on the settings
    in the DataDownloader Object.
    
    Parameters
    ----------
    srcDataset : DataImporter Object
        Contains various data for the source DICOM series.
    trgDataset : DataImporter Object
        Contains various data for the target DICOM series.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    
    
    Note
    ----
    There are two types of copying:
        1. Non-relationship-preserving (NRP)
        2. Relationship-preserving (RP) or propagation
    
    An ROI Collection can consist of one of three entity types:
        A. A single contour or segmentation
        B. A single collection of contours (ROI) or segmentations (segment)
        C. A collection of one or more ROIs or segments
    
    A RP copy/propagation of an ROI Collection will involve the entire 
    source entity (type A/B/C), and for types B/C, will return a the same
    entity type. For entity type A a single contour/segmentation may propagate
    to multiple contours/segmentations (depending on the relative frames of 
    reference and image sampling of the source and target images).
    
    By definition, a NRP copy is always from a single contour/segmentation to
    a single contour/segmentation (i.e. only entities of type A can be NRP 
    copied).  Hence the resulting entity is always a single contour/segmentation
    irrespective of how many frames are in the resampled/registration 
    transformed label image, for example.
    """
    
    def __init__(self, srcDataset, trgDataset, params, dro=None):
        
        # Copy selected instance attributes from srcDataset that will be
        # modified for useCases 3-5:
        self.dcmIm = srcDataset.dcmIm
        self.dcmPixarr = srcDataset.dcmPixarr
        
        # Copy selected instance attributes from trgDataset:
        self.imSize = trgDataset.imSize
        self.imSpacings = trgDataset.imSpacings
        self.imSlcThick = trgDataset.imSlcThick
        self.imPositions = trgDataset.imPositions
        self.imDirections = trgDataset.imDirections
        #self.dicomDir = trgDataset.dicomDir
        #self.dicoms = trgDataset.dicoms
        #self.divs = trgDataset.divs
        #self.f2sIndsByRoi = trgDataset.f2sIndsByRoi
        self.foruid = trgDataset.foruid
        #self.dcmIm = trgDataset.dcmIm # this will be changed for useCases 3-5; 14/09/21
        #self.pixarrByRoi = trgDataset.pixarrByRoi
        #self.roiName = trgDataset.roiName
        #self.roiNames = trgDataset.roiNames
        #self.roiNums = trgDataset.roiNums
        #self.roicol = trgDataset.roicol
        #self.seriesuid = trgDataset.seriesuid
        #self.slcNum = trgDataset.slcNum
        #self.sopuids = trgDataset.sopuids
        self.studyuid = trgDataset.studyuid
        
        # Copy srcDataset and trgDataset:
        #self.srcDataset = srcDataset
        #self.trgDataset = trgDataset
        
        # Initialise instance attributes that will result from making a
        # non-relationship-preserving copy, relationship-preserving copy or
        # propagation of the source ROI Collection to the target domain:
        #self.pixarrByRoi = None # 14/09/21
        #self.f2sIndsByRoi = None # 14/09/21
        #self.labimByRoi = None # 14/09/21
        #self.dcmIm = None
        self.f2sIndsBySeg = srcDataset.f2sIndsBySeg # 15/09/21
        self.pixarrBySeg = srcDataset.pixarrBySeg # 15/09/21
        self.labimBySeg = srcDataset.labimBySeg # 15/09/21
        self.c2sIndsByRoi = srcDataset.c2sIndsByRoi # 15/09/21
        self.ptsByCntByRoi = srcDataset.ptsByCntByRoi # 14/09/21
        self.cntdataByCntByRoi = srcDataset.cntdataByCntByRoi # 14/09/21
        self.labimByRoi = srcDataset.labimByRoi # 14/09/21
        
        # Initialise the SimpleITK Transform that will be used for resampling 
        # source label images to be the identity transform:
        """
        Note:
        
        If image registration is performed, resTx will be updated to be the 
        final registration transform. If a suitable DRO is found, resTx will be
        updated to be the transform created from the DRO and used to transform/
        resample the source label image.
        """
        #self.sitkTx = sitk.Transform()
        self.resTx = sitk.Transform(3, sitk.sitkIdentity)
        self.resTxParams = self.resTx.GetParameters()
        
        # Initialise instance attributes that will be updated following
        # registration or parsing of a DRO:
        self.initRegTx = None # initial registration transform
        self.initRegTxParams = None
        self.alignedIm = None # from resampling usinig initRegTx
        #self.finalTx = None # from registration or DRO
        #self.sitkTx = None
        #self.regIm = None
        
        # Initialise instance attributes that will be updated following
        # parsing of a DRO:
        #self.sitkTx = None
        #self.txParams = None # from registration or DRO
        #self.regIm = None
        self.preRegTx = None # will be updated for deformable DROs only
        self.preRegTxParams = None
        
        # TODO are inputs valid below
        """
        # Check if combination of input parameters are valid:
        valid, errMsg = are_inputs_valid(params, trgRoiNames)
        
        if not valid:
            # Raise exception:
            raise Exception(errMsg)
        """
        
        # Check whether the RTS/SEG of interest intersects the target image's
        # extent:
        self.get_intersection_of_roi_and_trgIm(srcDataset, trgDataset, params)
        
        # Determine which use case applies (and will be applied):
        self.which_use_case(srcDataset, trgDataset, params)
        
        if params.cfgDict['p2c']:
            print(f"runID = {params.cfgDict['runID']}")
            print(f"useCaseThatApplies = {params.cfgDict['useCaseThatApplies']}")
            print(f"useCaseToApply = {params.cfgDict['useCaseToApply']}\n")
        
        # Initialise list of timestamps and timing messages to be stored:
        self.timings = [time.time()]
        self.timingMsgs = []
        
        #print(f'srcIm direction:\n{srcDataset.dcmIm.GetDirection()}\n')
        #print(f'trgIm direction:\n{trgDataset.dcmIm.GetDirection()}\n')
        
        #breakpoint()
        
        #self.propagate(params, srcDataset, trgDataset)
    
    def get_intersection_of_roi_and_trgIm(
            self, srcDataset, trgDataset, params
            ):
        """
        Get the intersection of the ROI and target image domain as a fraction.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            Contains various data for the source DICOM series.
        trgDataset : DataImporter Object
            Contains various data for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.fracProp : float
            The fractional proportion (normalised to 1) of the points or voxels
            that define the contour(s)/segmentation(s) to be copied/propagated
            intersect the extent of trgIm.
        """
    
        cfgDict = params.cfgDict
        roicolMod = cfgDict['roicolMod']
        p2c = cfgDict['p2c']
        
        if roicolMod == 'SEG':
            fracProp = prop_of_segs_in_extent(
                pixarrBySeg=srcDataset.pixarrBySeg, 
                f2sIndsBySeg=srcDataset.f2sIndsBySeg, 
                srcIm=srcDataset.dcmIm, 
                trgIm=trgDataset.dcmIm, 
                p2c=p2c
                )
        else:
            fracProp = prop_of_rois_in_extent(
                ptsByCntByRoi=srcDataset.ptsByCntByRoi,
                trgIm=trgDataset.dcmIm, 
                p2c=p2c
                )
        
        self.fracProp = fracProp
            
    def which_use_case(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Determine which Use Cases applies.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        params.cfgDict['useCaseThatApplies'] : str
            The use case that applies.
        params.cfgDict['useCaseToApply'] : str
            The use case that will be applied.
        
        Note
        ----
        This is similar to the function which_use_case in 
        io_tools.inputs_checker.py, except that the functional form takes a
        cfgDict (params.cfgDict) and imports various data. This method uses
        data already imported via the DataImporter Objects srcDataset and
        trgDataset.
        """
        
        trgSlcNum = params.cfgDict['trgSlcNum']
        forceReg = params.cfgDict['forceReg']
        p2c = params.cfgDict['p2c']
        
        srcStudyUID = srcDataset.studyuid
        srcSeriesUID = srcDataset.seriesuid
        srcSOPuids = srcDataset.sopuids
        srcFORuid = srcDataset.foruid
        srcIPPs = srcDataset.imPositions
        srcDirs = srcDataset.imDirections
        srcSpacings = srcDataset.imSpacings
        srcST = srcDataset.imSlcThick
        srcSize = srcDataset.imSize
        #srcIm = srcDataset.image
        srcImExtent = srcDataset.imExtent
        
        trgStudyUID = trgDataset.studyuid
        trgSeriesUID = trgDataset.seriesuid
        trgSOPuids = trgDataset.sopuids
        trgFORuid = trgDataset.foruid
        trgIPPs = trgDataset.imPositions
        trgDirs = trgDataset.imDirections
        trgSpacings = trgDataset.imSpacings
        trgST = trgDataset.imSlcThick
        trgSize = trgDataset.imSize
        #trgIm = trgDataset.image
        trgImExtent = trgDataset.imExtent
        
        if trgSlcNum:
            if trgSlcNum > len(trgSOPuids) - 1:
                msg = f'The input argument "trgSlcNum" = {trgSlcNum} exceeds '\
                    + f'the number of Target DICOMs ({len(trgSOPuids)}).'
                raise Exception(msg)
        
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
            
        if srcFORuid == trgFORuid:
            sameFOR = True
        else:
            sameFOR = False
            
        if srcSOPuids == trgSOPuids:
            sameSOPs = True
        else:
            sameSOPs = False
        
        # Are the Source and Target spacings, directions and origins equal to  
        # within 1e-06?
        sameSpacings = are_items_equal_to_within_eps(srcSpacings, trgSpacings)
        sameDirs = are_items_equal_to_within_eps(srcDirs, trgDirs)
        sameOrigins = are_items_equal_to_within_eps(srcIPPs[0], trgIPPs[0])
        #sameST = are_items_equal_to_within_eps([srcST], [trgST])
        sameST = are_items_equal_to_within_eps(srcST, trgST)
        
        if srcSize == trgSize:
            sameSize = True
        else:
            sameSize = False
        
        # Check which case follows:
        if sameFOR:
            if sameDirs:
                """ 
                See Comment 11/02 #1 in which_use_case() in 
                io_tools.inputs_checker.py
                """
                if sameSpacings and sameOrigins:
                    if sameSeries:
                        # Direct copy
                        useCaseThatApplies = '1'
                    else:
                        if isinstance(trgSlcNum, int):
                            useCaseThatApplies = '2a' # Direct copy
                        else:
                            useCaseThatApplies = '2b' # Relationship-preserving copy
                else:
                    """ 
                    See Comment 11/02 #2 in which_use_case() in 
                    io_tools.inputs_checker.py
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
            
            # Are the image extents equal to within 1e-06?
            sameExtents = are_items_equal_to_within_eps(
                srcImExtent, trgImExtent
                )
            
            print('\n\n', '-'*120)
            print('Results of running which_use_case():')
            print('\n\n', '-'*120)
            print('Source and Target are in the/have:')
            
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
                print(f'   * same voxel sizes\n      Source/Target: {srcSpacings}')
            else:
                print(f'   * different voxel sizes:\n      Source: {srcSpacings}',
                      f'\n      Target: {trgSpacings}')
            
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
            if forceReg and useCaseThatApplies in ['3a', '3b', '4a', '4b']:
                useCaseToApply = useCaseThatApplies.replace('3', '5').replace('4', '5')
                
                msg = f'\n*UseCase {useCaseThatApplies} applies but since '\
                      + f'forceReg = {forceReg}, UseCase {useCaseToApply} '\
                      + 'will be applied.'
                print(msg)
            else:
                useCaseToApply = useCaseThatApplies
                
                print(f'\n*UseCase {useCaseThatApplies} applies.')
                
            print('-'*120)
        
        params.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        params.cfgDict['useCaseToApply'] = useCaseToApply
    
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
    
    def get_voxel_shift(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Determine the voxel shift required to make a non-relationship
        -preserving copy of an ROI. 
        
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
        self.voxShift : list of ints
            A list (for each dimension) of the required voxel shift.
        
        Note
        ----
        The shifting required will depend on whether the source pixel array(s)
        have been resampled or not. The source pixel array(s) for cases 1/2a
        have not been resampled, so to account for the shift from srcSlcNum to
        trgSlcNumpixel shifts are required in x, y and z directions for 
        relationship-preserving copies, and in only z for non-relationship
        -preserving copies.
        
        For cases 3a/4a/5a the source pixel array(s) have been resampled (or 
        transformed using a registration transform) so that the in-plane
        (x and y) components do not require further modification. However due 
        to the change from srcSlcNum to trgSlcNum, a shift will be required for
        the out-of-plane (z) components.
        
        Example uses: 
        - To compensate for differences in z-positions of a single-framed
        Source pixel array to be "direct" copied to a Target slice at a 
        different stack position; 
        - Compensating for differences in the origin of Source and Target 
        images for relationship-preserving copies of images with equal voxel
        spacings.
        
        At the moment there's no obvious way of knowing whether the pixel 
        arrays are from resampled label images or not. So will rely on the 
        value of params.cfgDict['useCaseToApply']. If useCaseToApply is '1' or
        '2a', shifting will be performed in all directions, else along the
        z-direction only.
        """
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of get_voxel_shift():')
            print('-'*120, '\n')
        
        # Have the source pixel arrays been resampled or not?
        #if hasattr(self, 'resPixarrByRoi'):
        #if useCaseToApply in ['1', '2a', '2b']:
        if useCaseToApply in ['1', '2a']:
            im = srcDataset.dcmIm
            #slcNum = srcDataset.slcNum # initial value
            #f2sIndsByRoi = srcDataset.f2sIndsByRoi # initial value
            #pixarrByRoi = srcDataset.pixarrByRoi # initial value
            
            #if useCaseToApply in ['1', '2a']:
            if p2c:
                print('Source pixel arrays have not been resampled, and',
                      'non-relationship-preserving copy to be made.',
                      '\nApplying pixel shifts along z direction only.\n')
            
            # Apply pixel shift in z direction only:
            shiftInX = False
            shiftInY = False
            shiftInZ = True
            """
            else:
                if p2c:
                    print('Source pixel arrays have not been resampled, and',
                          'relationship-preserving copy to be made.',
                          '\nApplying pixel shifts along all directions.\n')
                
                # Apply pixel shift in all directions:
                shiftInX = True
                shiftInY = True
                shiftInZ = True
            """
        
        elif useCaseToApply in ['3a', '4a']:
            """ 
            Determining the voxel shift is only relevant for non-relationship
            -preserving copies for use cases 3-4.
            """
            # Use the resampled source slice number for srcSlcNum, and the
            # resampled pixel arrays for srcPixarrByRoi:
            #im = self.dcmIm # 14/09/21
            im = srcDataset.dcmIm # 14/09/21
            #srcSlcNum = self.resSlcNum # does this still exist? (06/09/21)
            #srcPixarrByRoi = self.resPixarrByRoi # does this still exist? (06/09/21)
            #slcNum = self.slcNum # (06/09/21)
            #f2sIndsByRoi = self.f2sIndsByRoi # (06/09/21)
            #pixarrByRoi = self.pixarrByRoi # (06/09/21)
            
            if p2c:
                print('Source pixel arrays have been resampled so ',
                      'applying pixel shifts along z direction only.\n')
            
            # Apply pixel shift along z dimension only:
            shiftInX = False
            shiftInY = False
            shiftInZ = True
        
        """ Note:
        self.dcmIm was initialised in the Propagator class as trgIm, but 
        depending on the use case might have been replaced with the result of
        resampling srcIm to the target domain.
        """
        #srcIm = srcDataset.dcmIm
        slcNum = srcDataset.slcNum
        
        trgIm = trgDataset.dcmIm
        trgSlcNum = trgDataset.slcNum
        
        #print(f'slcNum = {slcNum}')
        #print(f'trgSlcNum = {trgSlcNum}')
        
        """
        if not hasattr(self, 'resPixarrByRoi'):
            # Replace all indices in source with srcDataset.slcNum with
            # trgDataset.slcNum:
            #newF2SindsByRoi = deepcopy(srcDataset.f2sIndsByRoi)
            #self.replace_indices(srcDataset, trgDataset)
            self.replace_indices(srcDataset, trgDataset)
            replacedF2SindsByRoi = self.replacedF2SindsByRoi
            #f2sIndsByRoi = self.f2sIndsByRoi
        else:
            # Pass self.resF2SindsByRoi instead:
            replacedF2SindsByRoi = self.resF2SindsByRoi
        """
        
        # Although the f2sInds will be shifted, the pixel data will not, so
        # simply make a copy of the source's pixarrByRoi:
        #shiftedPixarrByRoi = deepcopy(srcDataset.pixarrByRoi)
        
        """
        # Get voxel shift between srcSlcNum in srcIm and trgSlcNum in trgIm
        # in the trgIm domain:
        voxShift = get_voxel_shift_bt_slices(
            image0=srcIm, sliceNum0=srcSlcNum,
            image1=trgIm, sliceNum1=trgSlcNum,
            refImage=trgIm
            )
        """
        
        print(f'\n\nslcNum = {slcNum}')
        print(f'trgSlcNum = {trgSlcNum}')
        
        
        # Get voxel shift (in the trgIm domain) between slice slcNum in im and
        # slice trgSlcNum in trgIm:
        voxShift = get_voxel_shift_bt_slices(
            image0=im, sliceNum0=slcNum,
            image1=trgIm, sliceNum1=trgSlcNum,
            refImage=trgIm
            )
        
        if p2c:
            print(f'   voxShift between slices = {voxShift}')
            #print(f'   f2sIndsByRoi prior to shifting = {replacedF2SindsByRoi}')
            #print(f'   f2sIndsByRoi prior to shifting = {f2sIndsByRoi}')
            #print(f'   f2sIndsByRoi prior to shifting = {self.f2sIndsByRoi}')
        
        if not shiftInX:
            voxShift[0] = 0
        if not shiftInY:
            voxShift[1] = 0
        if not shiftInZ:
            voxShift[2] = 0
        
        if p2c:
            print(f'   voxShift that will be applied = {voxShift}')
            print('\nEnd of get_voxel_shift\n')
            print('-'*120)
        
        self.voxShift = voxShift
            
    def shift_frames(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Shift the pixels in frames based on the required voxel shift.
        
        For non-relationship-preserving copies of ROIs.
        
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
        self.pixarrByRoi : list of Numpy Arrays
            The result of shifting the voxels in srcPixarrByRoi to make a non-
            relationship-preserving copy for target.
        self.f2sIndsBySeg : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.shiftedPixarrByRoi.
        
        Note
        ----
        06/09/21: This was called shift_frames_in_pixarrByRoi prior to 06/09.
        """
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        roicolMod = params.cfgDict['roicolMod']
        p2c = params.cfgDict['p2c']
        voxShift = self.voxShift
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of shift_frames():')
            print('-'*120, '\n')
            
        #print(f'dir(self):\n{dir(self)}\n')
        
        # Have the source pixel arrays been resampled or not?
        #if hasattr(self, 'resPixarrByRoi'):
        if useCaseToApply in ['1', '2a', '2b']:
            #slcNum = srcDataset.slcNum # initial value
            if roicolMod == 'RTSTRUCT':
                _2sIndsBy_ = srcDataset.c2sIndsByRoi # initial value
                pixarrBy_ = srcDataset.pixarrByRoi # initial value
            else:
                _2sIndsBy_ = srcDataset.f2sIndsBySeg # initial value
                pixarrBy_ = srcDataset.pixarrBySeg # initial value
        
        elif useCaseToApply in ['3a', '3b', '4a', '4b']:
            # Use the resampled source slice number for srcSlcNum, and the
            # resampled pixel arrays for srcPixarrByRoi:
            #srcSlcNum = self.resSlcNum # does this still exist? (06/09/21)
            #srcPixarrByRoi = self.resPixarrByRoi # does this still exist? (06/09/21)
            #slcNum = self.slcNum # (06/09/21)
            if roicolMod == 'RTSTRUCT':
                _2sIndsBy_ = self.c2sIndsByRoi # (15/09/21)
                pixarrBy_ = self.pixarrByRoi # (15/09/21)
            else:
                _2sIndsBy_ = self.f2sIndsBySeg # (06/09/21)
                pixarrBy_ = self.pixarrBySeg # (06/09/21)
        
        if p2c:
            print(f'   _2sIndsBy_ prior to shifting = {_2sIndsBy_}')
            print(f'   voxShift that will be applied = {voxShift}')
        
        shifted_2SindsBy_ = [] # initial value
        shiftedPixarrBy_ = [] # initial value
        
        # Loop through each pixel array:
        #for s in range(len(srcPixarrByRoi)):
        for s in range(len(pixarrBy_)):
            # Proceed only if there is at least one frame in this pixel array:
            #if srcPixarrByRoi[s].shape[0]:
            if pixarrBy_[s].shape[0]:
                # Replace pixarrBy_[s] with the result of shifting the 
                # in-plane elements and add the voxel shift along z to:
                """
                shiftedPixarrBySeg.append(
                    shift_frame(
                        frame=srcPixarrByRoi[s], voxShift=voxShift
                        )
                    )
                """
                shiftedPixarrBy_.append(
                    shift_frame(
                        frame=pixarrBy_[s], voxShift=voxShift
                        )
                    )
                
                #F = len(replacedF2SindsByRoi[s])
                F = len(_2sIndsBy_[s])
                
                # Shift the contour-/frame-to-slice indices by voxShift[2] to 
                # account for the z-shift:
                """
                shiftedF2SindsBySeg.append(
                    [replacedF2SindsByRoi[s][i] + voxShift[2] for i in range(F)]
                    #[f2sIndsBySeg[s][i] + voxShift[2] for i in range(F)]
                    )
                """
                shifted_2SindsBy_.append(
                    [_2sIndsBy_[s][i] + voxShift[2] for i in range(F)]
                    )
                
                if p2c:
                    unique = get_unique_items(shiftedPixarrBy_[s])
                    
                    print(f'   There are {len(unique)} unique items in',
                          f'shiftedPixarrBy_[{s}] after shifting the frame')
        
        if p2c:
            print('   shifted_2SindsBy_ after shifting =',
                  f'{shifted_2SindsBy_}')
            print('end of shift_frames\n')
            print('-'*120)
        
        """
        self.shiftedPixarrBySeg = shiftedPixarrBySeg
        self.shiftedF2SindsBySeg = shiftedF2SindsBySeg
        #self.pixarrBySeg = shiftedPixarrBySeg
        #self.f2sIndsBySeg = shiftedF2SindsBySeg
        """
        if roicolMod == 'RTSTRUCT':
            self.c2sIndsByRoi = shifted_2SindsBy_
            self.pixarrByRoi = shiftedPixarrBy_
        else:
            self.f2sIndsBySeg = shifted_2SindsBy_
            self.pixarrBySeg = shiftedPixarrBy_
    
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
    
    def shift_points(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Shift the coordinates in contour data based on required voxel shift.
        
        For non-relationship-preserving copies of ROIs.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.shiftedPtsByCntByRoi : list of list of a list of a list of floats
            List (for each ROI) of a list (for each contour) of a list (for 
            each point) of a list (for each dimension) of coordinates.
        self.shiftedC2SindsByRoi : list of a list of ints
            The list (for each ROI) of a list (for each contour) of the
            contour-to-slice indices in self.shiftedPtsByCntByRoi.
        
        Note
        ----
        13/09/21: This was formally done in two separate functions called 
        ZshiftPtsByCntByRoi() and ShiftPtsByCntByRoi().
        
        The shifting required will depend on whether the source pixel array(s)
        have been resampled or not. The source pixel array(s) for Use cases 
        1/2a have not been resampled, so pixel shifts are required in x, y and
        z directions to account for the shift from srcSlcNum to trgSlcNum.
        
        For Use cases 3a/4a/5a the source pixel array(s) have been resampled
        (or transformed using a registration transform) so that the in-plane
        (x and y) components do not require further modification. However due 
        to the change from srcSlcNum to trgSlcNum, a shift will be required for
        the out-of-plane (z) components.
        
        Example uses: 
        - To compensate for differences in z-positions of a single-framed
        Source pixel array to be "direct" copied to a Target slice at a 
        different stack position; 
        - Compensating for differences in the origin of Source and Target 
        images for relationship-preserving copies of images with equal voxel
        spacings.
        
        At the moment there's no obvious way of knowing whether the pixel 
        arrays are from resampled label images or not. So will rely on the 
        value of params.cfgDict['useCaseToApply']. If useCaseToApply is '1' or
        '2a', shifting will be performed in all directions, else along the
        z-direction only.
        """
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        p2c = params.cfgDict['p2c']
        voxShift = self.voxShift
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of shift_points():')
            print('-'*120, '\n')
        
        #print(f'dir(self):\n{dir(self)}\n')
        
        # Have the source pixel arrays been resampled or not?
        #if hasattr(self, 'resPixarrByRoi'):
        if useCaseToApply in ['1', '2a', '2b']:
            #slcNum = srcDataset.slcNum # initial value
            c2sIndsByRoi = srcDataset.c2sIndsByRoi # initial value
            ptsByCntByRoi = srcDataset.ptsByCntByRoi # initial value
            refIm = trgDataset.dcmIm
            
            #print('\nPrior to running shift_ptsByCntByRoi:')
            #print_indsByRoi(c2sIndsByRoi)
            #print_ptsByCntByRoi(ptsByCntByRoi)
            
            shiftedPtsByCntByRoi, shiftedC2SindsByRoi = shift_ptsByCntByRoi(
                ptsByCntByRoi, c2sIndsByRoi, voxShift, refIm, p2c
                )
        
        if p2c:
            print(f'   voxShift between slices = {voxShift}')
            #print(f'   f2sIndsByRoi prior to shifting = {replacedF2SindsByRoi}')
            print(f'   c2sIndsByRoi prior to shifting = {c2sIndsByRoi}')
            #print(f'   f2sIndsByRoi prior to shifting = {self.f2sIndsByRoi}')
            #print('   shiftedF2SindsByRoi after shifting =',
            #      f'{shiftedF2SindsByRoi}')
            print('   shiftedC2SindsByRoi after shifting =',
                  f'{shiftedC2SindsByRoi}')
            #print('-'*120)
            print('\nEnd of shift_points\n')
            print('-'*120)
            
        self.c2sIndsByRoi = shiftedC2SindsByRoi
        self.ptsByCntByRoi = shiftedPtsByCntByRoi
    
    def add_modified_data_to_existing_trgDataset(self, trgDataset, params):
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
    
    def make_non_relationship_preserving_copy(
            self, srcDataset, trgDataset, params
            ):
        # TODO update docstrings
        """
        Perform a non-relationship-preserving copy of the segmentation(s)/
        contour(s) in srcDataset to trgDataset.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        self.shiftedPixarrByRoi : list of Numpy Arrays
            The result of shifting the voxels in srcPixarrByRoi to make a non-
            relationship-preserving copy for target.
        self.shiftedF2SindsByRoi : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.shiftedPixarrByRoi.
        self.pixarrByRoi : list of Numpy Arrays
            The list (for each ROI) of the pixel arrays of the non-relationship
            -preserving copy of the source ROI.
        self.f2sIndsByRoi : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.pixarrByRoi.
        params.timings : list of Time timestamps
            Additional timestamp appended.
        params.timingMsgs : list of strs
            Additional message appended.
        
        Note
        ----
        The copy will be done by manipulation of the source segmentation/contour
        data so that the ROI Collection will overlay as desired onto the target
        image without preserving the spatial relationships between the ROI and 
        its original image domain.
        """
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        roicolMod = params.cfgDict['roicolMod']
        
        timingMsg = "Making a non-relationship-preserving copy of the "\
                + "source ROI Collection.\n"
        params.add_timestamp(timingMsg)
        
        # Replace all indices in source with srcDataset.slcNum with
        # trgDataset.slcNum:
        """ This was moved to inside shift_frames: """
        #self.replace_indices(srcDataset, trgDataset)
        
        # Get the required voxel shift to account for change from srcSlcNum to 
        # trgSlcNum:
        self.get_voxel_shift(srcDataset, trgDataset, params)
        
        # Shift the points/frames based on the required voxel shift:
        if roicolMod == 'RTSTRUCT':
            if useCaseToApply in ['1', '2a']:
                # Replace indices in frame-to-slice indices:
                #self.replace_indices(srcDataset, trgDataset)
                self.c2sIndsByRoi = replace_ind_in_C2SindsByRoi(
                    c2sIndsByRoi=srcDataset.c2sIndsByRoi,
                    indToReplace=srcDataset.slcNum,
                    replacementInd=trgDataset.slcNum
                    )
                
                # Shift the out-of-plane (z) components of the points:
                self.ptsByCntByRoi = z_shift_ptsByCntByRoi(
                    ptsByCntByRoi=srcDataset.ptsByCntByRoi,
                    newZind=trgDataset.slcNum,
                    dcmIm=srcDataset.dcmIm
                    )
            elif useCaseToApply == '2b':
                # Shift points:
                self.shift_points(srcDataset, trgDataset, params)
            else:
                # Contour data was converted to pixel data, so shift frames:
                self.shift_frames(srcDataset, trgDataset, params)
        else:
            self.shift_frames(srcDataset, trgDataset, params)
        
        # Any pixel array(s) and corresponding F2Sinds for any ROI(s)/segment(s)
        # that may be present in trgDataset:
        self.add_modified_data_to_existing_trgDataset(trgDataset, params)
        
        timingMsg = "Took [*] s to make a non-relationship-preserving "\
                + "copy of the source ROI Collection.\n"
        params.add_timestamp(timingMsg)
    
    def shift_indices(
            self, srcDataset, trgDataset, params
            ):
        """
        Shift the indices to account for any differences in the origins of the
        source and target image domains.
        
        For relationship-preserving copies of ROIs.
    
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
    
        Returns
        -------
        self.f2sIndsBySeg : list of a list of ints
            List (for each segment/ROI) of a list (for each segmentation/
            contour) of slice numbers that correspond to each Source 
            segmentation/contour in the Target image domain.
        
        Note
        ----
        06/09/21 This was called get_trgF2SindsByRoi_from_srcF2SindsByRoi.
        """
        
        roicolMod = params.cfgDict['roicolMod']
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of shift_indices():')
            print('-'*120, '\n')
        
        srcIm = srcDataset.dcmIm
        trgIm = trgDataset.dcmIm
        
        #srcF2SindsBySeg = srcDataset.f2sIndsBySeg # 15/09/21
        
        """
        Note: In this function srcF2SindsByRoi will refer to either frame-to
        -slice indices by segment, or contour-to-slice indices by ROI.
        """
        if roicolMod == 'RTSTRUCT':
            srcF2SindsByRoi = srcDataset.c2sIndsByRoi  # 15/09/21
        else:
            srcF2SindsByRoi = srcDataset.f2sIndsBySeg  # 15/09/21
        
        srcOrigin = srcIm.GetOrigin()
        trgOrigin = trgIm.GetOrigin()
        
        #dim = len(srcOrigin)
        #mmDiff = [trgOrigin[i] - srcOrigin[i] for i in range(dim)]
        
        dz_mm = trgOrigin[2] - srcOrigin[2]
        
        dk = trgIm.GetSpacing()[2]
        
        dz_pix = round(dz_mm/dk)
        
        #trgF2SindsByRoi = []
        newF2SindsByRoi = []
        
        for r in range(len(srcF2SindsByRoi)):
            #trgF2Sinds = []
            newF2Sinds = []
            
            for c in range(len(srcF2SindsByRoi[r])):
                ##trgF2Sinds.append(srcF2SindsByRoi[r][c] + dz_pix)
                #trgF2Sinds.append(srcF2SindsByRoi[r][c] - dz_pix)
                newF2Sinds.append(srcF2SindsByRoi[r][c] - dz_pix)
                
            #trgF2SindsByRoi.append(trgF2Sinds)
            newF2SindsByRoi.append(newF2Sinds)
        
        #self.f2sIndsByRoi = trgF2SindsByRoi
        self.f2sIndsByRoi = newF2SindsByRoi
        
        if p2c:
                print(f'   srcF2SindsByRoi = {srcF2SindsByRoi}')
                #print(f'   trgF2SindsByRoi = {trgF2SindsByRoi}')
                print(f'   newF2SindsByRoi = {newF2SindsByRoi}')
                print('-'*120)

    def make_relationship_preserving_copy(
            self, params, srcDataset, trgDataset
            ):
        # TODO update docstrings
        """
        Perform a relationship-preserving copy of the segmentation(s)/
        contour(s) in srcDataset to trgDataset.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
            
        Returns
        -------
        self.pixarrByRoi : list of Numpy Arrays
            List (for each segment/ROI) of the relationship-preserved copy of 
            the pixel array representations of each source segmentation/contour.
        self.f2sIndsByRoi : list of a list of ints
            List (for each segment/ROI) of a list (for each segmentation/
            contour) in self.pixarrByRoi.
        params.timings : list of Time timestamps
            Additional timestamp appended.
        params.timingMsgs : list of strs
            Additional message appended.
        
        Note
        ----
        The copy will be done by manipulating the source segmentation/contour
        data such that spatial relationships are preserved between the ROI and
        the original image domain.
        """
        
        timingMsg = "Making a relationship-preserving "\
                + "copy of the source ROI Collection.\n"
        params.add_timestamp(timingMsg)
        
        # Shift the frame-to-slice indices to account for differences in the 
        # origins of the images:
        self.shift_indices(srcDataset, trgDataset, params)
        
        """
        Since the spatial relationships are not to be altered, 
        newPixarrByRoi = srcPixarrBySeg, ...:
        """
        #self.pixarrByRoi = deepcopy(srcDataset.pixarrByRoi)
        #self.pixarrBySeg = list(srcDataset.pixarrBySeg) # this is in __init__
        
        timingMsg = "Took [*] s to make a relationship-preserving "\
                + "copy of the source ROI Collection.\n"
        params.add_timestamp(timingMsg)
    
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
        
    def resample_src_labims(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Resample the source labimByRoi using the target image grid.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        resTx : SimpleITK Transform
            The transform to use when resampling the source image to the target
            domain.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        self.labimByRoi :  list of SimpleITK Images
            The list (for each ROI) of resampled 3D label image.
        self.pixarrByRoi : list of Numpy Arrays
            The list (for each ROI) of the pixel array representations of
            self.resLabimByRoi.
        self.f2sIndsByRoi : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.resPixarrByRoi.
        params.timings : list of Time timestamps
            Additional timestamp appended.
        params.timingMsgs : list of strs
            Additional message appended.
        """
        
        roicolMod = params.cfgDict['roicolMod']
        
        if roicolMod == 'RTSTRUCT':
            _2sIndsBy_ = srcDataset.c2sIndsByRoi
            labimBy_ = srcDataset.labimByRoi
        else:
            _2sIndsBy_ = srcDataset.f2sIndsBySeg
            labimBy_ = srcDataset.labimBySeg
            
        labimBy_, pixarrBy_, _2sIndsBy_ = resample_labimBySeg(
            labimBySeg=labimBy_,
            f2sIndsBySeg=_2sIndsBy_,
            im=srcDataset.dcmIm,
            refIm=trgDataset.dcmIm,
            sitkTx=self.resTx, # 03/09/21
            #sitkTx=self.sitkTx, # 01/09/21
            #sitkTx=self.finalTx, # 01/09/21
            interp=params.cfgDict['resInterp'], 
            applyPreResBlur=params.cfgDict['applyPreResBlur'],
            preResVar=params.cfgDict['preResVar'],
            applyPostResBlur=params.cfgDict['applyPostResBlur'],
            postResVar=params.cfgDict['postResVar'],
            p2c=params.cfgDict['p2c']
            )
        
        """
        Although the desired interpolation (resInterp) was provided, the actual
        interpolation that was used might have been different. 
        
        Metadata was written to the resampled label images (see 
        image_tools.resampling.resample_labim for details). Pull this metadata
        from the first label image and add to params.cfgDict.
        """
        
        params.cfgDict['resInterpUsed']\
            = labimBy_[0].GetMetaData('resInterpUsed')
        
        if params.cfgDict['p2c']:
            print(f"resInterp = {params.cfgDict['resInterp']}")
            print(f"resInterpUsed = {params.cfgDict['resInterpUsed']}")
            print('\n\n\nMetadata keys in labimBy_[0]:',
                  labimBy_[0].GetMetaDataKeys())
        
        timingMsg = "Took [*] s to resample source label images.\n"
        params.add_timestamp(timingMsg)
        
        if roicolMod == 'RTSTRUCT':
            self.c2sIndsByRoi = _2sIndsBy_
            self.pixarrByRoi = pixarrBy_
            self.labimByRoi = labimBy_
        else:
            self.f2sIndsBySeg = _2sIndsBy_
            self.pixarrBySeg = pixarrBy_
            self.labimBySeg = labimBy_
    
    def register_image(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Register the source labimByRoi using the target image grid.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        self.dcmIm : SimpleITK Image
            The source image registered to the target image.
        self.initialTx : SimpleITK Transform
            The transform used during initialisation of the registration.
        self.resTx : list of floats
            The transform that registers source image to target image.
        params.timings : list of Time timestamps
            Additional timestamp appended.
        params.timingMsgs : list of strs
            Additional message appended.
        """
        
        fixIm = trgDataset.dcmIm
        movIm = srcDataset.dcmIm
        regTxName = params.cfgDict['regTxName']
        initMethod = params.cfgDict['initMethod']
        fixFidsFpath = params.cfgDict['trgFidsFpath']
        movFidsFpath = params.cfgDict['srcFidsFpath']
        p2c = params.cfgDict['p2c']
        
        # Required for generating a file path for exported registeration plot:
        resPlotsExportDir = params.cfgDict['resPlotsExportDir']
        srcExpLab = params.cfgDict['srcExpLab']
        srcScanID = params.cfgDict['srcScanID']
        trgExpLab = params.cfgDict['trgExpLab']
        trgScanID = params.cfgDict['trgScanID']
        
        #initMethod = 'geometry'
        #print('\n* On 06/09 changed initMethod from "moments" to "geometry"\n\n')
        
        resPlotFname = generate_reg_fname(
                srcExpLab, srcScanID, trgExpLab, trgScanID, regTxName
                )
        
        resPlotFpath = os.path.join(resPlotsExportDir, resPlotFname)
        
        # Create directory if it doesn't already exist:
        if not os.path.isdir(resPlotsExportDir):
            #os.mkdir(exportDir)
            Path(resPlotsExportDir).mkdir(parents=True)
        
        timingMsg = "* Registering the source to target DICOM scans using "\
            + f"{regTxName} registration with {initMethod} initialisation...\n"
        params.add_timestamp(timingMsg)
        
        print(f'fixFidsFpath = {fixFidsFpath}')
        print(f'movFidsFpath = {movFidsFpath}\n')
        
        self.initRegTx, self.alignedIm, self.resTx, self.dcmIm,\
            self.metricValues, self.multiresIters = register_im(
                fixIm=fixIm, movIm=movIm, 
                regTxName=regTxName, initMethod=initMethod,
                fixFidsFpath=fixFidsFpath, 
                movFidsFpath=movFidsFpath,
                p2c=p2c, regPlotFpath=resPlotFpath
                )
        
        self.dcmPixarr = sitk.GetArrayViewFromImage(self.dcmIm)
        
        print(f'\nmetricValues = {self.metricValues}')
        print(f'\nmultiresIters = {self.multiresIters}\n')
        
        #raise Exception('quitting')
        
        # Plot registration result:
        midInd = fixIm.GetSize()[2] // 2
        
        compare_res_results(
            resIm0=fixIm, resIm1=self.dcmIm, resInd=midInd,
            resTitle0='Target image', resTitle1='Registered'
        )
        
        # Get the registration transform parameters:
        #self.txParams = [str(item) for item in self.finalTx.GetParameters()]  # 01/09/21
        self.resTxParams = self.resTx.GetParameters() # 03/09/21
        self.initRegTxParams = self.initRegTx.GetParameters() # 03/09/21
        
        timingMsg = "Took [*] s to register the source to target DICOM scans.\n"
        params.add_timestamp(timingMsg)
    
    def create_tx_from_dro(self, srcDataset, trgDataset, dro, params):
        # TODO update docstrings
        """
        Create a SimpleITK Transform from a DRO.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        dro : Pydicom Object
            The DICOM Registration Object as a Pydicom Object.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
         
        Returns
        -------
        self.resTx : SimpleITK Transform
            The transform parsed from the DRO.
        self.txParams : list of floats
            The transform parameters (i.e. flat matrix) corresponding to
            self.sitkTx.
        self.dcmIm : SimpleITK Image
            The resampled source image based on self.sitkTx (analogue of the
            registered image - see Note).
        
        Note
        ----
        Although not required, the source image is resampled below so that the
        result is accessible from self.regIm as would have been the case
        following image registration. Note the use of regIm rather than resIm
        even though it did not result from image registration. This was because
        the resampling using the transform parsed from the DRO replaced 
        registration and hence consistency in attributes is maintained.
        """
        
        srcIm = srcDataset.dcmIm
        trgIm = trgDataset.dcmIm
        p2c = params.cfgDict['p2c']
        
        if params.cfgDict['regTxName'] in ['rigid', 'affine']:
            resTx = create_tx_from_spa_dro(dro, p2c)
            self.resTx = resTx # update resTx
            
            resTxParams = list(resTx.GetParameters())
            self.resTxParams = resTxParams # update resTxParams
            
            #listOfSitkTxs.append(sitkTx)
            
            # Resample srcIm:
            resIm = resample_im(
                im=srcIm, refIm=trgIm, sitkTx=resTx, interp='Linear', 
                p2c=False
                )
        else:
            # Create the deformable SimpleITK Transform:
            resTx = create_tx_from_def_dro(dro, refIm=trgIm, p2c=p2c)
            
            # Create the pre-deformable SimpleITK Transform:
            preRegTx = create_pre_tx_from_def_dro(dro, p2c)
            self.preRegTx = preRegTx
            self.preRegTxParams = list(preRegTx.GetParameters())
            
            IDENTITY = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
            
            if are_items_equal_to_within_eps(self.preRegTxParams, IDENTITY):
                """ 
                The pre-registration matrix is simply the identity matrix,
                so no need to resample/transform the label image using
                preRegTx - only need to use sitkTx. 
                """
                # Resample srcIm using sitkTx:
                resIm = resample_im(
                    im=srcIm, refIm=trgIm, sitkTx=resTx, 
                    interp='Linear', p2c=False
                    )
                
                # Update resTx and resTxParams:
                self.resTx = resTx
                self.resTxParams = list(resTx.GetParameters())
            else:
                """ 
                The pre-registration matrix is not the identity matrix so
                need to resample/transform the label image using the composite
                transform of preRegTx and sitkTx. 
                """
                # Compose transforms: (09/07/21)
                compTx = sitk.CompositeTransform(preRegTx)
                compTx.AddTransform(resTx)
            
                # Resample srcIm usig compTx:
                resIm = resample_im(
                    im=srcIm, refIm=trgIm, sitkTx=compTx, 
                    interp='Linear', p2c=False
                    )
                
                # Update resTx and resTxParams:
                self.resTx = compTx
                self.resTxParams = list(compTx.GetParameters())
        
        self.dcmIm = resIm
    
    def reduce_frames(self, params):
        # TODO update docstrings
        """
        Reduce the number of frames in multi-framed pixel arrays to a
        single frame using averaging or a logical or operation.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        self.pixarrBySeg : list of Numpy Arrays
            List (for each ROI/segment) of resampled pixel arrays for source.
        self.f2sIndsBySeg : list of floats
            List (for each ROI/segment) of resampled pixel arrays for target.
        
        Returns
        -------
        self.pixarrBySeg : list of Numpy Arrays
            List (for each ROI/segment) of single-framed pixel arrays for 
            source by frame averaging or a logical or operation.
        self.f2sIndsBySeg : list of floats
            List (for each ROI/segment) of single-framed pixel arrays for 
            target by frame averaging or a logical or operation.
        self.slcNum : int
            The slice index corresponding to the reduced resampled pixel array
            for source.
        """
        
        roicolMod = params.cfgDict['roicolMod']
        p2c = params.cfgDict['p2c']
        
        if roicolMod == 'RTSTRUCT':
            pixarrBy_ = self.pixarrByRoi
            _2sIndsBy_ = self.c2sIndsByRoi
        else:
            pixarrBy_ = self.pixarrBySeg
            _2sIndsBy_ = self.f2sIndsBySeg
        
        R = len(pixarrBy_) # number of segments
        
        # Get the number of frames in each pixel array:
        numOfFramesBy_ = [pixarrBy_[r].shape[0] for r in range(R)]
        
        maxNumOfFramesBy_ = max(numOfFramesBy_)
        
        if maxNumOfFramesBy_ > 1:
            # Reduce the number of frames in each pixel array to 1:
            
            # Perform averaging or or operation?
            #reduceUsing = 'mean'
            reduceUsing = 'or'
            
            if p2c:
                print(f'\nThere are {numOfFramesBy_} frames in each',
                      f'ROI/segment. Using {reduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                plot_pixarrBySeg(
                    pixarrBy_, _2sIndsBy_, 
                    f'Prior to {reduceUsing} operation'
                    )
                
            if reduceUsing == 'mean':
                """ 14/02:  
                There seems to be a problem with the averaging operation that 
                results in many more than 1 f2sInd and many more than 1 contour
                """
                
                # The threshold used when converting the averaged (non-binary)  
                # pixel array to a binary pixel array:
                binaryThresh = 0.5 # <-- is this too high?
    
                pixarrBy_, _2sIndsBy_ = mean_frame_in_pixarrBySeg(
                    pixarrBySeg=pixarrBy_,
                    f2sIndsBySeg=_2sIndsBy_,
                    makeBinary=True, 
                    binaryThresh=binaryThresh,
                    p2c=p2c
                    )
            else:
                pixarrBy_, _2sIndsBy_ = or_frame_of_pixarrBySeg(
                    pixarrBySeg=pixarrBy_,
                    f2sIndsBySeg=_2sIndsBy_,
                    p2c=p2c
                    )
            
            if p2c:
                numOfFramesBy_ = [
                    pixarrBy_[r].shape[0] for r in range(R)
                    ]
        
                maxNumOfFramesBy_ = max(numOfFramesBy_)
                
                print(f'\nFollowing {reduceUsing} operation there are',
                      f'{numOfFramesBy_} frames in each ROI/segment, and',
                      f'the new _2sIndsBy_ = {_2sIndsBy_}.')
                
                plot_pixarrBySeg(
                    pixarrBy_, _2sIndsBy_, 
                    'After {reduceUsing} operation'
                    )
        
        """ 
        Note:
            srcSlcNum --> resSlcNum following resampling or transformation. 
            Since the pixel arrays in Source's _2sIndsBy_ are single framed, 
            resSlcNum is simply the first index in resampled _2sIndsBy_:
        """
        resSlcNum = _2sIndsBy_[0][0]
        
        if roicolMod == 'RTSTRUCT':
            self.pixarrByRoi = pixarrBy_
            self.c2sIndsByRoi = _2sIndsBy_
        else:
            self.pixarrBySeg = pixarrBy_
            self.f2sIndsBySeg = _2sIndsBy_
        self.slcNum = resSlcNum
    
    def make_non_relationship_preserving_propagation(
            self, srcDataset, trgDataset, params
            ):
        # TODO update docstrings
        """
        Perform a non-relationship-preserving propagation of the 
        segmentation(s)/contour(s) in srcDataset to trgDataset.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        
        Note
        ----
        
        """
        
        """
        May need to reduce multi-framed NRP copies to a single framed pixel 
        array (using frame averaging or logical OR operation).
        """
        self.reduce_frames(params)
        
        # Get the required voxel shift to account for change from srcSlcNum 
        # to trgSlcNum:
        self.get_voxel_shift(srcDataset, trgDataset, params)
    
        # Shift the frames to account for change from srcSlcNum to 
        # trgSlcNum:
        self.shift_frames(srcDataset, trgDataset, params)
        
        # Any pixel array(s) and corresponding F2Sinds for any ROI(s)/
        # segment(s) that may be present in trgDataset:
        self.add_modified_data_to_existing_trgDataset(trgDataset, params)
    
    def convert_pixarr_to_cntdata(
            self, srcDataset, trgDataset, newDataset, params
            ):
        # TODO update docstrings
        """
        Convert pixel array (segmentation) data to contour data, and update the
        relevant object attributes.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        self.pixarrBySeg : list of Numpy Arrays
            List (for each ROI/segment) of resampled pixel arrays for source.
        self.f2sIndsBySeg : list of floats
            List (for each ROI/segment) of resampled pixel arrays for target.
        
        Returns
        -------
        self.ptsByCntByRoi : list of list of a list of a list of floats
            List (for each ROI) of a list (for all contours) of a list (for 
            each point) of a list (for each dimension) of coordinates.
        self.cntdataByCntByRoi : list of a list of a list of strs
            List (for each ROI) of a list (for all contours) of a flat list of 
            coordinates in ptsByCntByRoi converted from floats to strings.
        """
        
        pass
    
    def execute(self, srcDataset, trgDataset, params, dro):
        # TODO update docstrings
        """
        Execute the methods that copy/propagate the source ROI Collection to
        the target domain.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        dro : Pydicom Object
            The DICOM Registration Object as a Pydicom Object.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self....
        
        Notes
        -----
        Propagation will be either by resampling the source label images to the
        target domain, or by transforming them using the transform that 
        registers the source image to the target image, or parsed from an
        appropriate spatial or deformable DRO.
        
        Source segmentation (pixel) data will be converted to SimpleITK Image(s)
        (i.e. label images), and back to the pixel array representation 
        following the propagation, i.e.:
            Source Pixel data/array(s) --> Label image(s) 
            Propagated to the Target domain
            Target Label image(s) --> Pixel array(s)/data
        
        Source contour data will be converted to SimpleITK Image(s) via an 
        intermediate conversion to pixel data, and back to the contour data
        representation following propagation via an intermediate conversion to
        pixel data, i.e.:
            Source Contour data --> Pixel array(s) --> Label image(s) 
            Propagated to the Target domain
            Target Label image(s) --> Pixel array(s) --> Contour data
            
        Even when label images are transformed (i.e. resampled rather than 
        registred) the prefix 'res' will be used since further operations will
        need to be applied to either resampled or transformed label images, and
        hence it's useful for instance attribute names to be consistent. 
        """
        
        cfgDict = params.cfgDict
        
        useCase = cfgDict['useCaseToApply']
        roicolMod = cfgDict['roicolMod']
        #resInterp = cfgDict['resInterp']
        #applyPreResBlur = cfgDict['applyPreResBlur']
        #preResVar = cfgDict['preResVar']
        #applyPostResBlur = cfgDict['applyPostResBlur']
        #postResVar = cfgDict['postResVar']
        useDroForTx = cfgDict['useDroForTx']
        #regTxName = cfgDict['regTxName']
        p2c = cfgDict['p2c']
        
        #srcLabimByRoi = srcDataset.labimByRoi
        #srcF2SindsByRoi = srcDataset.f2sIndsBySeg
        #srcIm = srcDataset.image
        #trgIm = trgDataset.image
        
        
        # TODO where to import DRO?
        #dro = None # for now
        
        if p2c:
            print(f"params.cfgDict['runID'] = {params.cfgDict['runID']}")
            print(f'useCase = {useCase}')
        
        if useCase in ['1', '2a']:
            """ 
            NRP copy without resampling or registration transformation.
            """
            if p2c:
                print('Running useCase in ["1", "2a"]\n')
            # Make a non-relationship-preserving copy:
            self.make_non_relationship_preserving_copy(
                srcDataset, trgDataset, params
                )
            
        if useCase == '2b':
            """ 
            RP copy/propagation without resampling or registration 
            transformation.
            """
            if p2c:
                print('Running useCase = "2b"\n')
            # Make a relationship-preserving copy:
            self.make_relationship_preserving_copy(
                params, srcDataset, trgDataset
                )
        
        #if useCase in ['3a', '3b', '4a', '4b', '5a', '5b']:
        #    doesExist = does_instance_variable_exist(
        #        instanceObj=srcDataset, varName='labimByRoi', p2c=p2c
        #        )
        #    
        #    # Get a list of label images-by-ROI for source:
        #    #srcDataset.get_labimByRoi(params, srcDataset)
        #    
        #    """ 
        #    Above line doesn't work since srcDataset -> DataImporter object,
        #    does not have method get)labimByRoi. This has been moved to
        #    get_seg_metadata() in io_tools.import_data, even though it's not
        #    needed for use cases 1-2. 
        #    """
        #    
        #    doesExist = does_instance_variable_exist(
        #        instanceObj=srcDataset, varName='labimByRoi', p2c=p2c
        #        )
        
        if useCase in ['3a', '3b', '4a', '4b']:
            """
            NRP copy (3a, 4a) or RP propagation (3b, 4b) by resampling using
            the identity transform.
            """
            if p2c:
                print('Running useCase in ["3a", "3b", "4a" and "4b"]\n')
            self.resample_src_labims(srcDataset, trgDataset, params)
            
            """ 
            Although resampling of the source image to the target domain is 
            not required for copying/propagating the source ROI Collection to  
            the target domain, it is useful to have the resampled source image 
            (e.g. for overlays of the resampled ROI Collection on the resampled
            DICOM image).
            """ 
            self.dcmIm = resample_im(
                srcDataset.dcmIm, refIm=trgDataset.dcmIm,
                sitkTx=self.resTx, p2c=params.cfgDict['p2c']
                )
            
            self.dcmPixarr = sitk.GetArrayViewFromImage(self.dcmIm)
        
        if useCase in ['5a', '5b']:
            """
            NRP copy (5a) or RP propagation (5b) with registration 
            transformation.
            """
            if p2c:
                print('Running useCase in ["5a", "5b"]\n')
            
            #if dro == None or (dro != None and useDroForTx): # 01/09/21
            if dro == None or (dro != None and not useDroForTx): # 01/09/21
                """
                No suitable DRO was found on XNAT or user wishes to 
                override the DRO and perform image registration.
                """
                
                if (p2c and dro != None and useDroForTx):
                    print('Image registeration will be performed instead',
                          'of using the DRO from XNAT.\n')
                
                self.register_image(srcDataset, trgDataset, params)
            else:
                # TODO register_image creates self.dcmIm and self.dcmPixarr,
                # but those attributes aren't created in this case...
                """
                A suitable DRO was found on XNAT to bypass image 
                registration.
                """
                
                self.create_tx_from_dro(srcDataset, trgDataset, dro, params)
            
            # Resample the source label images using the registration 
            # transform (i.e. transform the source label images):
            self.resample_src_labims(srcDataset, trgDataset, params)
            
        if useCase in ['3a', '4a', '5a']:
            self.make_non_relationship_preserving_propagation(
                srcDataset, trgDataset, params
                )
        
        """
        if useCase in ['3b', '4b', '5b']:
            self.pixarrBySeg = self.resPixarrByRoi
            self.f2sIndsBySeg = self.resF2SindsByRoi
        """
        
        if (useCase in ['3a', '3b', '4a', '4b', '5a', '5b'] and
                roicolMod == 'RTSTRUCT'):
            # Convert from pixel array data to contour data:
            print('\nNeed to convert from pixel array to contour data...')
            
    
    def export_data(self, srcDataset, trgDataset, params):
        
        #exportLabim = params.cfgDict['exportLabim']
        #exportTx = params.cfgDict['exportTx']
        
        labimExportDir = params.cfgDict['labimExportDir']
        txExportDir = params.cfgDict['txExportDir']
        
        if not os.path.isdir(txExportDir):
            Path(txExportDir).mkdir(parents=True)
        
        srcLabimByRoi = srcDataset.labimByRoi
        trgLabimByRoi = trgDataset.labimByRoi
        newLabimByRoi = self.labimByRoi
        
        resTx = self.resTx
        initRegTx = self.initRegTx
        preRegTx = self.preRegTx
        
        #resTxParams = self.resTxParams
        #initRegTxParams = self.initRegTxParams
        #preRegTxParams = self.preRegTxParams
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        useDroForTx = params.cfgDict['useDroForTx']
        regTxName = params.cfgDict['regTxName']
        
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        """ Export label images """
        for labimByRoi, ds in zip(
                [srcLabimByRoi, trgLabimByRoi, newLabimByRoi],
                ['src', 'trg', 'new']
                ):
            if labimByRoi:
                for r in range(len(labimByRoi)):
                    fname = f'{currentDateTime}_{ds}Labim_{r}'
                    
                    fpath = os.path.join(labimExportDir, fname)
                    
                    
                
        
        """ Export transforms """
        
        # Beginning of file name will be current datetime:
        resFname = f'{currentDateTime}_'
        initRegFname = f'{currentDateTime}_'
        preRegFname = f'{currentDateTime}_'
        
        if useCaseToApply in ['3a', '3b', '4a', '4b']:
            resFname += 'resample_tx'
            
        elif useCaseToApply in ['5a', '5a']:
            if useDroForTx:
                resFname += f'res_tx_{regTxName}_from_DRO'
                
                if regTxName == 'bspline':
                    preRegFname += 'pre_reg_tx_from_DRO'
                
            else:
                resFname += f'reg_tx_{regTxName}'
                initRegFname += 'reg_init_tx'
        
        for tx, fname in zip(
                [resTx, initRegTx, preRegTx], 
                [resFname, initRegFname, preRegFname]
                ):
            if tx:
                # Write transform in TFM and HDF formats:
                fpath = os.path.join(txExportDir, fname + '.tfm')
                sitk.WriteTransform(tx, fpath)
                
                fpath = os.path.join(txExportDir, fname + '.hdf')
                sitk.WriteTransform(tx, fpath)
                
                # Write transform parameters in TXT format:
                export_list_to_txt(
                    items=tx.GetParameters(), 
                    filename=fname+'_params.txt', 
                    exportDir=txExportDir
                    )
        