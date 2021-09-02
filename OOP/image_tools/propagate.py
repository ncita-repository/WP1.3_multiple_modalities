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
import plotting_tools.plotting
reload(plotting_tools.plotting)

import time
import os
from pathlib import Path
from copy import deepcopy
import SimpleITK as sitk

from general_tools.shifting import get_voxel_shift_bt_slices, shift_frame
from general_tools.general import (
    get_unique_items, are_items_equal_to_within_eps, generate_reg_fname
    #does_instance_variable_exist
    )
from conversion_tools.pixarrs_ims import pixarr_to_im
from image_tools.attrs_info import get_im_info
from image_tools.resampling import resample_im, resample_labimByRoi
from image_tools.registering import register_im
from dro_tools.create_tx_from_dro import (
    create_tx_from_spa_dro, create_pre_tx_from_def_dro, create_tx_from_def_dro
    )
from general_tools.pixarr_ops import (
    mean_frame_in_pixarrBySeg, or_frame_of_pixarrBySeg
    )
from plotting_tools.plotting import plot_pixarrBySeg

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
        
        #cfgDict = params.cfgDict
        #useCase = cfgDict['useCaseToApply']
        #useDroForTx = cfgDict['useDroForTx']
        #regTxName = cfgDict['regTxName']
        #p2c = cfgDict['p2c']
        
        # TODO where to import DRO?
        #dro = None # for now
        
        ## Create a new target dataset:
        #self.newDataset = deepcopy(trgDataset)
        #self = deepcopy(trgDataset)
        
        # Copy selected instance attributes from trgDataset:
        #self.dicomDir = trgDataset.dicomDir
        #self.dicoms = trgDataset.dicoms
        #self.divs = trgDataset.divs
        #self.f2sIndsByRoi = trgDataset.f2sIndsByRoi
        self.foruid = trgDataset.foruid
        #self.image = trgDataset.image
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
        self.pixarrByRoi = None
        self.f2sIndsByRoi = None
        self.labimByRoi = None
        self.image = None
        
        # Initialise sitkTx that will be used for resampling source label
        # images to be the identity transform:
        #self.sitkTx = sitk.Transform()
        self.sitkTx = sitk.Transform(3, sitk.sitkIdentity)
        
        # Initialise instance attributes that will be updated following
        # registration:
        self.initialTx = None
        self.finalTx = None
        #self.sitkTx = None
        #self.regIm = None
        
        # Initialise instance attributes that will be updated following
        # parsing of a DRO:
        #self.sitkTx = None
        self.txParams = None
        #self.regIm = None
        """
        Note:
            self.finalTx and self.regIm will be updated after parsing a DRO.
        """
        
        # Determine which use case applies (and will be applied):
        self.which_use_case(srcDataset, trgDataset, params)
        
        if params.cfgDict['p2c']:
            print(f"runID = {params.cfgDict['runID']}")
            print(f"useCaseThatApplies = {params.cfgDict['useCaseThatApplies']}")
            print(f"useCaseToApply = {params.cfgDict['useCaseToApply']}\n")
        
        # Initialise list of timestamps and timing messages to be stored:
        self.timings = [time.time()]
        self.timingMsgs = []
        
        #self.propagate(params, srcDataset, trgDataset)
    
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
    
    def replace_ind_in_f2sIndsByRoi(self, srcDataset, trgDataset):
        # TODO modfy the docstrings
        """
        Replace an index in a list of frame-to-slice indices.
        
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
        replacedF2SindsByRoi = deepcopy(srcF2SindsByRoi)
        
        #print(f'type(replacedF2SindsByRoi) = {type(replacedF2SindsByRoi)}\n')
        
        for r in range(len(srcF2SindsByRoi)):
            for c in range(len(srcF2SindsByRoi[r])):
                if srcF2SindsByRoi[r][c] == indToReplace:
                    replacedF2SindsByRoi[r][c] = replacementInd
        
        self.replacedF2SindsByRoi = replacedF2SindsByRoi
        #self.f2sIndsByRoi = replacedF2SindsByRoi
    
    def shift_frames_in_pixarrByRoi(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        19/08/21:
            The problem with wrapping this in a method in this way is that
            the arguments required for use case 3a/4a/5a is diffent than that
            required for use case 1/2a.
            
        Shift the pixels in each frame of srcPixarrByRoi. 
        
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
        self.shiftedPixarrByRoi : list of Numpy Arrays
            The result of shifting the voxels in srcPixarrByRoi to make a non-
            relationship-preserving copy for target.
        self.shiftedF2SindsByRoi : list of a list of ints
            The list (for each ROI) of a list (for each frame) of the
            frame-to-slice indices in self.shiftedPixarrByRoi.
        
        Note
        ----
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
            images for relationship-preserving copies of images with equal 
            voxel spacings.
        """
        
        p2c = params.cfgDict['p2c']
        
        # Are we working with resampled source pixel arrays or not?
        """ Could check value of self.useCase but this might be better. """ 
        if hasattr(self, 'resPixarrByRoi'):
            if p2c:
                print('Source pixel arrays have been resampled, so applying',
                      'pixel shifts along z direction only.\n')
            
            # Apply pixel shift along z dimension only:
            shiftInX = False
            shiftInY = False
            shiftInZ = True
            
            # Since the source pixel array has been resampled it shares the 
            # same image grid as the target image, so use the target image for
            # srcIm: 
            srcIm = trgDataset.image
            
            # Use the resampled source slice number for srcSlcNum, and the
            # resampled pixel arrays for srcPixarrByRoi:
            srcSlcNum = self.resSlcNum
            srcPixarrByRoi = self.resPixarrByRoi
        else:
            if p2c:
                print('Source pixel arrays have not been resampled, so',
                      'applying pixel shifts along all directions.\n')
            
            # Apply pixel shift in all dimensions:
            shiftInX = True
            shiftInY = True
            shiftInZ = True
                                    
            srcIm = srcDataset.image
            srcSlcNum = srcDataset.slcNum
            srcPixarrByRoi = srcDataset.pixarrByRoi
        
        trgIm = trgDataset.image
        trgSlcNum = trgDataset.slcNum
        
        if not hasattr(self, 'resPixarrByRoi'):
            # Replace all indices in source with srcDataset.slcNum with
            # trgDataset.slcNum:
            #newF2SindsByRoi = deepcopy(srcDataset.f2sIndsByRoi)
            #self.replace_ind_in_f2sIndsByRoi(srcDataset, trgDataset)
            self.replace_ind_in_f2sIndsByRoi(srcDataset, trgDataset)
            replacedF2SindsByRoi = self.replacedF2SindsByRoi
            #f2sIndsByRoi = self.f2sIndsByRoi
        else:
            # Pass self.resF2SindsByRoi instead:
            replacedF2SindsByRoi = self.resF2SindsByRoi
        
        # Although the f2sInds will be shifted, the pixel data will not, so
        # simply make a copy of the source's pixarrByRoi:
        #shiftedPixarrByRoi = deepcopy(srcDataset.pixarrByRoi)
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of shift_frames_in_pixarrByRoi():')
            print('-'*120, '\n')
            
        # Get voxel shift between srcSlcNum in srcIm and trgSlcNum in trgIm
        # in the trgIm domain:
        voxShift = get_voxel_shift_bt_slices(
            image0=srcIm, sliceNum0=srcSlcNum, image1=trgIm, 
            sliceNum1=trgSlcNum, refImage=trgIm
            )
        
        if p2c:
            print(f'   voxShift between slices = {voxShift}')
            print(f'   f2sIndsByRoi prior to shifting = {replacedF2SindsByRoi}')
            #print(f'   f2sIndsByRoi prior to shifting = {f2sIndsByRoi}')
        
        if not shiftInX:
            voxShift[0] = 0
        if not shiftInY:
            voxShift[1] = 0
        if not shiftInZ:
            voxShift[2] = 0
        
        if p2c:
            print(f'   voxShift that will be applied = {voxShift}')
        
        shiftedF2SindsByRoi = [] # initial value
        shiftedPixarrByRoi = [] # initial value
        
        # Loop through each pixel array:
        for s in range(len(srcPixarrByRoi)):
            # Proceed only if there is at least one frame in this pixel array:
            if srcPixarrByRoi[s].shape[0]:
                # Replace srcPixarrByRoi[s] with the result of shifting the 
                # in-plane elements and add the pixel shift along z to
                # replacedF2SindsByRoi:
                shiftedPixarrByRoi.append(
                    shift_frame(
                        frame=srcPixarrByRoi[s], voxShift=voxShift
                        )
                    )
                
                F = len(replacedF2SindsByRoi[s])
                #F = len(f2sIndsByRoi[s])
                
                # Shift the frame-to-slice indices by voxShift[2] to account
                # for the z-shift:
                shiftedF2SindsByRoi.append(
                    [replacedF2SindsByRoi[s][i] + voxShift[2] for i in range(F)]
                    #[f2sIndsByRoi[s][i] + voxShift[2] for i in range(F)]
                    )
                
                if p2c:
                    unique = get_unique_items(shiftedPixarrByRoi[s])
                    
                    print(f'   There are {len(unique)} unique items in',
                          f'shiftedPixarrByRoi[{s}] after shifting the frame')
        
        if p2c:
            print('   shiftedF2SindsByRoi after shifting =',
                  f'{shiftedF2SindsByRoi}')
            print('-'*120)
        
        self.shiftedPixarrByRoi = shiftedPixarrByRoi
        self.shiftedF2SindsByRoi = shiftedF2SindsByRoi
        #self.pixarrByRoi = shiftedPixarrByRoi
        #self.f2sIndsByRoi = shiftedF2SindsByRoi
    
    def add_modified_data_to_existing_trgDataset(self, trgDataset, params):
        """
        Concatenate a list (by segment) of pixel arrays (pixarrByRoi) and of 
        frame-to-slice indices (f2sIndsByRoi) to the corresponding instance
        variables in trgDataset.
        
        The relevant instance variables in trgDataset may be None, e.g. if 
        there was no ROI Collection imported for the target DICOM series
        i.e. trgDataset.roicol = None). If so f2sIndsByRoi and pixarrByRoi will 
        will be shiftedF2SindsByRoi and shiftedPixarrByRoi.
        
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
        
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of add_modified_data_to_existing_trgDataset():')
            print('-'*120, '\n')
        
        shiftedF2SindsByRoi = self.shiftedF2SindsByRoi
        shiftedPixarrByRoi = self.shiftedPixarrByRoi
        
        if trgDataset.roicol == None:
            self.f2sIndsByRoi = shiftedF2SindsByRoi
            self.pixarrByRoi = shiftedPixarrByRoi
            
            if p2c:
                print('   There are no pixel arrays to add from trgDataset')
                print('   self.f2sIndsByRoi = shiftedF2SindsByRoi =',
                      f'{shiftedF2SindsByRoi}')
                #print('   self.pixarrByRoi = shiftedPixarrByRoi')
                print('-'*120)
        else:
            # Concatenate the existing data from target with the shifted data:
            trgF2SindsByRoi = trgDataset.f2sIndsByRoi
            trgPixarrByRoi = trgDataset.pixarrByRoi
            
            S = len(trgF2SindsByRoi)
            
            newF2SindsByRoi = [
                shiftedF2SindsByRoi[s] + trgF2SindsByRoi[s] for s in range(S)
                ]
            
            newPixarrByRoi = [
                shiftedPixarrByRoi[s] + trgPixarrByRoi[s] for s in range(S)
                ]
            
            self.f2sIndsByRoi = newF2SindsByRoi
            self.pixarrByRoi = newPixarrByRoi
            
            if p2c:
                print('   There were pixel arrays to add from trgDataset')
                print('   f2sIndsByRoi from shiftedF2SindsByRoi =',
                      f'{shiftedF2SindsByRoi}')
                print('   f2sIndsByRoi after adding existing f2sIndsByRoi from',
                      f'trgDataset = {newF2SindsByRoi}')
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
        
        timingMsg = "Making a non-relationship-preserving copy of the "\
                + "source ROI Collection.\n"
        params.add_timestamp(timingMsg)
        
        # Replace all indices in source with srcDataset.slcNum with
        # trgDataset.slcNum:
        """ This was moved to inside shift_frames_in_pixarrByRoi: """
        #self.replace_ind_in_f2sIndsByRoi(srcDataset, trgDataset)
        
        # Shift the frames to account for change from srcSlcNum to trgSlcNum:
        self.shift_frames_in_pixarrByRoi(srcDataset, trgDataset, params)
        
        # Any pixel array(s) and corresponding F2Sinds for any ROI(s)/segment(s)
        # that may be present in trgDataset:
        self.add_modified_data_to_existing_trgDataset(trgDataset, params)
        
        timingMsg = "Took [*] s to make a non-relationship-preserving "\
                + "copy of the source ROI Collection.\n"
        params.add_timestamp(timingMsg)
    
    def get_trgF2SindsByRoi_from_srcF2SindsByRoi(
            self, params, srcDataset, trgDataset
            ):
        """
        Shift the indices in srcF2SindsByRoi to account for any differences in 
        the origins of the source and target image domains.
    
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
        self.f2sIndsByRoi : list of a list of ints
            List (for each segment/ROI) of a list (for each segmentation/
            contour) of slice numbers that correspond to each Source 
            segmentation/contour in the Target image domain.
        """
        
        p2c = params.cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print(' Running of get_trgF2SindsByRoi_from_srcF2SindsByRoi():')
            print('-'*120, '\n')
        
        srcIm = srcDataset.image
        srcF2SindsByRoi = srcDataset.f2sIndsByRoi
        trgIm = trgDataset.image
        
        srcOrigin = srcIm.GetOrigin()
        trgOrigin = trgIm.GetOrigin()
        
        #dim = len(srcOrigin)
        #mmDiff = [trgOrigin[i] - srcOrigin[i] for i in range(dim)]
        
        dz_mm = trgOrigin[2] - srcOrigin[2]
        
        dk = trgIm.GetSpacing()[2]
        
        dz_pix = round(dz_mm/dk)
        
        trgF2SindsByRoi = []
        
        for r in range(len(srcF2SindsByRoi)):
            trgF2Sinds = []
            
            for c in range(len(srcF2SindsByRoi[r])):
                #trgF2Sinds.append(srcF2SindsByRoi[r][c] + dz_pix)
                trgF2Sinds.append(srcF2SindsByRoi[r][c] - dz_pix)
                
            trgF2SindsByRoi.append(trgF2Sinds)
        
        self.f2sIndsByRoi = trgF2SindsByRoi
        
        if p2c:
                print(f'   srcF2SindsByRoi = {srcF2SindsByRoi}')
                print(f'   trgF2SindsByRoi = {trgF2SindsByRoi}')
                print('-'*120)

    def make_relationship_preserving_copy(self, params, srcDataset, trgDataset):
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
        
        # Shift the indices is srcF2SindsByRoi to get trgF2SindsByRoi:
        self.get_trgF2SindsByRoi_from_srcF2SindsByRoi(
            params, srcDataset, trgDataset
            )
        
        # Since the spatial relationships are not to be modified, 
        # trgPixarrByRoi = srcPixarrBySeg:
        self.pixarrByRoi = deepcopy(srcDataset.pixarrByRoi)
        
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
        refIm = dataset.image
        
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
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        #sitkTx : SimpleITK Transform, optional
        #    The transform to use when resampling the source image to the target
        #    domain. Default value is sitk.Transform() (identity transform).
            
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
        
        self.labimByRoi, self.pixarrByRoi, self.f2sIndsByRoi\
                = resample_labimByRoi(
                    labimByRoi=srcDataset.labimByRoi,
                    f2sIndsByRoi=srcDataset.f2sIndsByRoi,
                    im=srcDataset.image,
                    refIm=trgDataset.image,
                    #sitkTx=sitkTx,
                    #sitkTx=self.sitkTx, # 01/09/21
                    sitkTx=self.finalTx, # 01/09/21
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
            = self.labimByRoi[0].GetMetaData('resInterpUsed')
        
        if params.cfgDict['p2c']:
            print(f"resInterp = {params.cfgDict['resInterp']}")
            print(f"resInterpUsed = {params.cfgDict['resInterpUsed']}")
            print('\n\n\nMetadata keys in resSrcLabimByRoi[0]:',
                  self.labimByRoi[0].GetMetaDataKeys())
        
        timingMsg = "Took [*] s to resample source label images.\n"
        params.add_timestamp(timingMsg)
    
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
        self.image : SimpleITK Image
            The source image registered to the target image.
        self.initialTx : SimpleITK Transform
            The transform used during initialisation of the registration.
        self.finalTx : list of floats
            The transform that registers source image to target image.
        params.timings : list of Time timestamps
            Additional timestamp appended.
        params.timingMsgs : list of strs
            Additional message appended.
        """
        
        fixIm = trgDataset.image
        movIm = srcDataset.image
        regTxName = params.cfgDict['regTxName']
        initMethod = params.cfgDict['initMethod']
        fixFidsFpath = params.cfgDict['trgFidsFpath']
        movFidsFpath = params.cfgDict['srcFidsFpath']
        p2c = params.cfgDict['p2c']
        
        # Required for generating a file path for exported registeration plot:
        regPlotsExportDir = params.cfgDict['regPlotsExportDir']
        srcExpLab = params.cfgDict['srcExpLab']
        srcScanID = params.cfgDict['srcScanID']
        trgExpLab = params.cfgDict['trgExpLab']
        trgScanID = params.cfgDict['trgScanID']
        
        regPlotFname = generate_reg_fname(
                srcExpLab, srcScanID, trgExpLab, trgScanID, regTxName
                )
        regPlotFpath = os.path.join(regPlotsExportDir, regPlotFname)
        
        # Create directory if it doesn't already exist:
        if not os.path.isdir(regPlotsExportDir):
            #os.mkdir(exportDir)
            Path(regPlotsExportDir).mkdir(parents=True)
        
        timingMsg = "* Registering the source to target DICOM scans using "\
            + f"{regTxName} registration with {initMethod} initialisation...\n"
        params.add_timestamp(timingMsg)
        
        print(f'fixFidsFpath = {fixFidsFpath}')
        print(f'movFidsFpath = {movFidsFpath}\n')
        
        self.image, self.initialTx, self.finalTx, self.metricValues,\
            self.multiresIters = register_im(
                fixIm=fixIm, movIm=movIm, 
                regTxName=regTxName, initMethod=initMethod,
                fixFidsFpath=fixFidsFpath, 
                movFidsFpath=movFidsFpath,
                p2c=p2c, regPlotFpath=regPlotFpath
                )
        
        # Get the registration transform parameters:
        #self.txParams = [str(item) for item in self.finalTx.GetParameters()]  # 01/09/21
        self.txParams = self.finalTx.GetParameters() # 01/09/21
        
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
        self.finalTx : SimpleITK Transform
            The transform parsed from the DRO.
        self.txParams : list of floats
            The transform parameters (i.e. flat matrix) corresponding to
            self.sitkTx.
        self.image : SimpleITK Image
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
        
        srcIm = srcDataset.image
        trgIm = trgDataset.image
        p2c = params.cfgDict['p2c']
        
        if params.cfgDict['regTxName'] in ['rigid', 'affine']:
            finalTx = create_tx_from_spa_dro(dro, p2c)
            self.finalTx = finalTx
            
            txParams = list(finalTx.GetParameters())
            self.txParams = txParams
            
            #listOfSitkTxs.append(sitkTx)
            
            # Resample srcIm:
            resIm = resample_im(
                    im=srcIm, refIm=trgIm, sitkTx=finalTx, 
                    interp='Linear', p2c=False
                    )
        else:
            # Create the deformable SimpleITK Transform:
            finalTx = create_tx_from_def_dro(dro, refIm=trgIm, p2c=p2c)
            self.finalTx = finalTx
            
            # Create the pre-deformable SimpleITK Transform:
            preRegTx = create_pre_tx_from_def_dro(dro, p2c)
            self.preRegTx = preRegTx
            
            txParams = list(preRegTx.GetParameters())
            self.txParams = txParams
            
            IDENTITY = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
            
            if are_items_equal_to_within_eps(txParams, IDENTITY):
                """ 
                The pre-registration matrix is simply the identity matrix,
                so no need to resample/transform the label image using
                preRegTx - only need to use sitkTx. 
                """
                #finalTx = create_tx_from_def_dro(Dro, refIm=trgIm, p2c=p2c)
                
                # Resample srcIm using sitkTx:
                resIm = resample_im(
                    im=srcIm, refIm=trgIm, sitkTx=finalTx, 
                    interp='Linear', p2c=False
                    )
            else:
                """ 
                The pre-registration matrix is not the identity matrix so
                need to resample/transform the label image using the composite
                transform of preRegTx and sitkTx. 
                """
                #listOfSitkTxs.append(preRegTx)
                
                #finalTx = create_tx_from_def_dro(Dro, refIm=trgIm, p2c=p2c)
                
                # Compose transforms: (09/07/21)
                compTx = sitk.CompositeTransform(preRegTx)
                compTx.AddTransform(finalTx)
                # Update self.sitkTx:
                self.finalTx = compTx
                
                # Resample srcIm usig compTx:
                resIm = resample_im(
                    im=srcIm, refIm=trgIm, sitkTx=compTx, 
                    interp='Linear', p2c=False
                    )
        
        self.image = resIm
    
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
        self.pixarrByRoi : list of Numpy Arrays
            List (for each ROI/segment) of resampled pixel arrays for source.
        self.f2sIndsByRoi : list of floats
            List (for each ROI/segment) of resampled pixel arrays for target.
        
        Returns
        -------
        self.pixarrByRoi : list of Numpy Arrays
            List (for each ROI/segment) of single-framed pixel arrays for 
            source by frame averaging or a logical or operation.
        self.f2sIndsByRoi : list of floats
            List (for each ROI/segment) of single-framed pixel arrays for 
            target by frame averaging or a logical or operation.
        self.slcNum : int
            The slice index corresponding to the reduced resampled pixel array
            for source.
        """
        
        pixarrByRoi = self.pixarrByRoi
        f2sIndsByRoi = self.f2sIndsByRoi
        p2c = params.cfgDict['p2c']
        
        R = len(pixarrByRoi) # number of segments
        
        # Get the number of frames in each pixel array:
        numOfFramesByRoi = [pixarrByRoi[r].shape[0] for r in range(R)]
        
        maxNumOfFramesByRoi = max(numOfFramesByRoi)
        
        if maxNumOfFramesByRoi > 1:
            # Reduce the number of frames in each pixel array to 1:
            
            # Perform averaging or or operation?
            #reduceUsing = 'mean'
            reduceUsing = 'or'
            
            if p2c:
                print(f'\nThere are {numOfFramesByRoi} frames in each',
                      f'ROI/segment. Using {reduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                plot_pixarrBySeg(
                    pixarrByRoi, f2sIndsByRoi, 
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
    
                pixarrByRoi, f2sIndsBySeg = mean_frame_in_pixarrBySeg(
                    pixarrBySeg=pixarrByRoi,
                    f2sIndsBySeg=f2sIndsByRoi,
                    makeBinary=True, 
                    binaryThresh=binaryThresh,
                    p2c=p2c
                    )
            else:
                pixarrByRoi, f2sIndsByRoi = or_frame_of_pixarrBySeg(
                    pixarrBySeg=pixarrByRoi,
                    f2sIndsBySeg=f2sIndsByRoi,
                    p2c=p2c
                    )
            
            if p2c:
                numOfFramesByRoi = [
                    pixarrByRoi[r].shape[0] for r in range(R)
                    ]
        
                maxNumOfFramesByRoi = max(numOfFramesByRoi)
                
                print(f'\nFollowing {reduceUsing} operation there are',
                      f'{numOfFramesByRoi} frames in each ROI/segment, and',
                      f'the new f2sIndsByRoi = {f2sIndsByRoi}.')
                
                plot_pixarrBySeg(
                    pixarrByRoi, f2sIndsByRoi, 
                    'After {reduceUsing} operation'
                    )
        
        """ 
        Note:
            srcSlcNum --> resSlcNum following resampling or transformation. 
            Since the pixel arrays in srcF2SindsBySeg are single framed,
            resSlcNum is simply the first frame-to-slice index in
            resF2SindsBySeg:
        """
        resSlcNum = f2sIndsByRoi[0][0]
        
        self.pixarrByRoi = pixarrByRoi
        self.f2sIndsByRoi = f2sIndsByRoi
        self.slcNum = resSlcNum
    
    #def create_roicol(self, ):
    
    def propagate(self, srcDataset, trgDataset, params, dro):
        # TODO update docstrings
        """
        Propagate the source ROI Collection to the target domain.
        
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
        #resInterp = cfgDict['resInterp']
        #applyPreResBlur = cfgDict['applyPreResBlur']
        #preResVar = cfgDict['preResVar']
        #applyPostResBlur = cfgDict['applyPostResBlur']
        #postResVar = cfgDict['postResVar']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        p2c = cfgDict['p2c']
        
        #srcLabimByRoi = srcDataset.labimByRoi
        #srcF2SindsByRoi = srcDataset.f2sIndsByRoi
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
                print('Running useCase in [1, 2a]\n')
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
                print('Running useCase = 2b\n')
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
            NRP copy (3a, 4a) or RP propagation (3b, 4b) with resampling (not 
            registration transformation).
            """
            if p2c:
                print('Running useCase in [3a, 3b, 4a and 4b]\n')
            self.resample_src_labims(srcDataset, trgDataset, params)
        
        if useCase in ['5a', '5b']:
            """
            NRP copy (5a) or RP propagation (5b) with registration 
            transformation.
            """
            if p2c:
                print('Running useCase in [5a, 5b]\n')
            
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
                """
                A suitable DRO was found on XNAT to bypass image 
                registration.
                """
                
                self.create_tx_from_dro(srcDataset, trgDataset, dro, params)
            
            # Resample the source label images using the registration 
            # transform (i.e. transform the source label images):
            self.resample_src_labims(srcDataset, trgDataset, params)
            
        if useCase in ['3a', '4a', '5a']:
            """
            May need to reduce multi-framed NRP copies to a single framed pixel 
            array (using frame averaging or logical OR operation).
            """
            
            self.reduce_frames(params)
            
            # Shift the frames to account for change from srcSlcNum to trgSlcNum:
            self.shift_frames_in_pixarrByRoi(srcDataset, trgDataset, params)
            
            # Any pixel array(s) and corresponding F2Sinds for any ROI(s)/segment(s)
            # that may be present in trgDataset:
            self.add_modified_data_to_existing_trgDataset(trgDataset, params)
            
        """
        if useCase in ['3b', '4b', '5b']:
            self.pixarrByRoi = self.resPixarrByRoi
            self.f2sIndsByRoi = self.resF2SindsByRoi
        """
        