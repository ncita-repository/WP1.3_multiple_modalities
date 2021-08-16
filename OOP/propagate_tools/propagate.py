# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 17:05:12 2021

@author: ctorti
"""

from importlib import reload

import general_tools.shifting
reload(general_tools.shifting)
import image_tools.resampling
reload(image_tools.resampling)

from copy import deepcopy

from general_tools.shifting import get_voxel_shift_bt_slices, shift_frame
from general_tools.general import get_unique_items, are_items_equal_to_within_eps
from general_tools.general import does_instance_variable_exist
from conversion_tools.pixarrs_ims import pixarr_to_im
from image_tools.attrs_info import get_im_info
from image_tools.resampling import resample_labimByRoi

class Propagator:
    # TODO modify the docstrings
    """
    This class propagates Importer Objects based on the settings in the Fetcher
    Object.
    
    Parameters
    ----------
    params : DataDownloader Object
        Contains various parameters.
    srcDataset : DataImporter Object
        Contains various data for the source DICOM series.
    trgDataset : DataImporter Object
        Contains various data for the target DICOM series.
    
    Returns
    -------
    
    """
    
    def __init__(self, params, srcDataset, trgDataset, dro=None):
        
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
        self.dicoms = trgDataset.dicoms
        self.divs = trgDataset.divs
        self.f2sIndsByRoi = trgDataset.f2sIndsByRoi
        self.foruid = trgDataset.foruid
        self.image = trgDataset.image
        self.pixarrByRoi = trgDataset.pixarrByRoi
        self.roiName = trgDataset.roiName
        self.roiNames = trgDataset.roiNames
        self.roiNums = trgDataset.roiNums
        self.roicol = trgDataset.roicol
        self.seriesuid = trgDataset.seriesuid
        self.slcNum = trgDataset.slcNum
        self.sopuids = trgDataset.sopuids
        self.studyuid = trgDataset.studyuid
        
        # Copy srcDataset and trgDataset:
        #self.srcDataset = srcDataset
        #self.trgDataset = trgDataset
        
        self.propagate(params, srcDataset, trgDataset)
    
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
        replacedF2SindsByRoi : list of list of ints
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
    
    def shift_frames_in_pixarrByRoi(self, params, srcDataset, trgDataset):
        # TODO update docstrings
        """
        Shift the items in each frame of each pixel array in a list of 
        pixel-arrays-by-segment.  
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains various parameters.
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        
        Returns
        -------
        self. : DataImporter Object
            Copy of srcDataset with changes to .f2sIndsByRoi and .pixarrByRoi
            attributes.
        
        Notes
        -----
        
        
        Example uses: 
            - To compensate for differences in z-positions of a single-framed 
            Source pixel array to be "direct" copied to a Target slice at a 
            different stack position; 
            - Compensating for differences in the origin of Source and Target 
            images for relationship-preserving copies of images with equal 
            voxel spacings.
        """
        
        p2c = params.cfgDict['p2c']
        
        # Apply pixel shift in all dimensions:
        shiftInX = True
        shiftInY = True
        shiftInZ = True
                                
        srcIm = srcDataset.image
        srcSlcNum = srcDataset.slcNum
        srcPixarrByRoi = srcDataset.pixarrByRoi
        trgIm = trgDataset.image
        trgSlcNum = trgDataset.slcNum
        
        # Replace all indices in source with srcDataset.slcNum with
        # trgDataset.slcNum:
        #newF2SindsByRoi = deepcopy(srcDataset.f2sIndsByRoi)
        #self.replace_ind_in_f2sIndsByRoi(srcDataset, trgDataset)
        self.replace_ind_in_f2sIndsByRoi(srcDataset, trgDataset)
        replacedF2SindsByRoi = self.replacedF2SindsByRoi
        
        ## Although the f2sInds will be shifted, the pixel data will not, so
        ## simply make a copy of the source's pixarrByRoi:
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
                
                # Shift the frame-to-slice indices by voxShift[2] to account
                # for the z-shift:
                shiftedF2SindsByRoi.append(
                    [replacedF2SindsByRoi[s][i] + voxShift[2] for i in range(F)]
                    )
                
                if p2c:
                    unique = get_unique_items(shiftedPixarrByRoi[s])
                    
                    print(f'   There are {len(unique)} unique items in',
                          f'shiftedPixarrByRoi[{s}] after shifting the frame')
        
        if p2c:
            print('   shiftedF2SindsByRoi after shifting =',
                  f'{shiftedF2SindsByRoi}')
            print('-'*120)
        
        self.shiftedF2SindsByRoi = shiftedF2SindsByRoi
        self.shiftedPixarrByRoi = shiftedPixarrByRoi
    
    def add_modified_data_to_existing_trgDataset(self, params, trgDataset):
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
        params : DataDownloader Object
            Contains various parameters.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
            
        Returns
        -------
        None.
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
    
    def make_non_relationship_preserving_copy(self, params, srcDataset, trgDataset):
        # TODO update docstrings
        """
        Perform a non-relationship-preserving copy of the segmentation(s)/
        contour(s) in srcDataset to trgDataset.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains various parameters.
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
            
        Returns
        -------
        None.
        """
        
        # Replace all indices in source with srcDataset.slcNum with
        # trgDataset.slcNum:
        """ This was moved to inside shift_frames_in_pixarrByRoi: """
        #self.replace_ind_in_f2sIndsByRoi(srcDataset, trgDataset)
        
        # Shift the frames in srcDataset.pixarrByRoi accordingly:
        self.shift_frames_in_pixarrByRoi(params, srcDataset, trgDataset)
        
        # Any pixel array(s) and corresponding F2Sinds for any ROI(s)/segment(s)
        # that may be present in trgDataset:
        self.add_modified_data_to_existing_trgDataset(params, trgDataset)
    
    def get_trgF2SindsByRoi_from_srcF2SindsByRoi(self, params, srcDataset, 
                                                 trgDataset):
        """
        Shift the indices in srcF2SindsByRoi to account for any differences in 
        the origins of the source and target image domains.
    
        Parameters
        ----------
        params : DataDownloader Object
            Contains various parameters.
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
    
        Returns
        -------
        None
        
        Note
        ----
        trgF2SindsByRoi : list of a list of ints
            List (for each segment/ROI) of a list (for each segmentation/contour) 
            of slice numbers that correspond to each Source segmentation/contour in 
            the Target image domain.
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
            Contains various parameters.
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
            
        Returns
        -------
        None.
        """
        
        # Shift the indices is srcF2SindsByRoi to get trgF2SindsByRoi:
        self.get_trgF2SindsByRoi_from_srcF2SindsByRoi(
            params, srcDataset, trgDataset
            )
        
        # Since the spatial relationships are not to be modified, 
        # trgPixarrByRoi = srcPixarrBySeg:
        self.pixarrByRoi = deepcopy(srcDataset.pixarrByRoi)
    
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
        
    #def resample_images(self, ):
    
    
    #def register_images(self, ):
        
    
    #def reduce_frames(self, ):
        
    
    #def create_roicol(self, ):
    
    def propagate(self, params, srcDataset, trgDataset, dro):
        
        cfgDict = params.cfgDict
        
        useCase = cfgDict['useCaseToApply']
        resInterp = cfgDict['resInterp']
        applyPreResBlur = cfgDict['applyPreResBlur']
        preResVariance = cfgDict['preResVariance']
        applyPostResBlur = cfgDict['applyPostResBlur']
        postResVariance = cfgDict['postResVariance']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        p2c = cfgDict['p2c']
        
        srcLabimByRoi = srcDataset.labimByRoi
        srcF2SindsByRoi = srcDataset.f2sIndsByRoi
        srcIm = srcDataset.image
        trgIm = trgDataset.image
        
        
        # TODO where to import DRO?
        #dro = None # for now
        
        if p2c:
            print(f"params.cfgDict['runID'] = {params.cfgDict['runID']}")
            print(f'useCase = {useCase}')
        
        if useCase in ['1', '2a']:
            # Make a non-relationship-preserving copy:
            self.make_non_relationship_preserving_copy(
                params, srcDataset, trgDataset
                )
            
        if useCase == '2b':
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
            This used to be under:
                if useCase in ['3a', '3b', '4a', '4b', '5a', '5b']
            but the label images by ROI is fetched along with other data.
            """
            print('Need to do useCase 3a, 3b, 4a and 4b')
            resSrcLabimByRoi, resSrcPixarrByRoi, resSrcF2SindsByRoi\
                = resample_labimByRoi(
                    labimByRoi=srcLabimByRoi, 
                    f2sIndsByRoi=srcF2SindsByRoi, 
                    im=srcIm, 
                    refIm=trgIm,  
                    interp=resInterp, 
                    applyPreResBlur=applyPreResBlur, 
                    preResVariance=preResVariance, 
                    applyPostResBlur=applyPostResBlur, 
                    postResVariance=postResVariance, 
                    p2c=p2c
                    )
        
        else:
            if dro == None or (dro != None and useDroForTx):
                """
                No suitable DRO was found on XNAT or user wishes to 
                override the DRO and perform image registration.
                """
                
                if (dro != None and useDroForTx):
                    print('Image registeration will be performed instead',
                          'of using the DRO from XNAT.\n')
                
                print('Need to write register_images()')
                #self.register_images()
                
                # Get the reg tx params
                
                # more...
            else:
                """
                A suitable DRO was found on XNAT to bypass image 
                registration.
                """
                
                if regTxName in ['rigid', 'affine']:
                    print('Need to write create_tx_from_spa_dro()')
                    #self.create_tx_from_spa_dro()
                    
                    # do other stuff
                    
                else:
                    print('Need to write create_tx_from_def_dro()')
                    #self.create_tx_from_def_dro()
                    
                    print('Need to write create_pre_tx_from_def_dro()')
                    #self.create_pre_tx_from_def_dro()
                    
                    # do other stuff
                    
                    if are_items_equal_to_within_eps():
                        print('Need to write resample_image()')
                        #self.resample_image()
                    
                    else:
                        # Composite transforms
                        print('Need to write resample_image()')
                        #self.resample_image()
                        
                        # TODO consider using only one call of resample_image
                        # and instead modify self.sitkTx accordingly
                
                print('Need to write resample_labimByRoi()')
                #self.resample_labimByRoi()
                    
                # TODO what is the else with RegIm = None, SitkTx = None...
                # at the bottom of page 10 for?
            
        if useCase in ['3a', '4a', '5a']:
            print('Need to write reduce_frames()')
            #self.reduce_frames()
            
            # shift_frames_in_pixarrByRoi (pg 12)
        
            # concatenate...
            
        #if useCase in ['3b', '4b', '5b']:
            # deepcopying may not be required..
            
        
        """ Create ROI and DRO in separate classes... """
        # Create/overwrite the Target ROI Collection:
        print('Need to write create_roicol()')
        #trgDataset.create_roicol()