# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:47:21 2021

@author: ctorti
"""

from importlib import reload

#import xnat_tools.dros
#reload(xnat_tools.dros)

#import io_tools.fetch_params_and_data
#reload(io_tools.fetch_params_and_data)

import io_tools.import_dro
reload(io_tools.import_dro)
import io_tools.import_roicol
reload(io_tools.import_roicol)
import io_tools.inputs_checker
reload(io_tools.inputs_checker)
import seg_tools.metadata
reload(seg_tools.metadata)
import seg_tools.get_seg_data_of_interest
reload(seg_tools.get_data_of_interest)


#from xnat_tools.sessions import create_session
#from io_tools.fetch_params_and_data import Params
#from xnat_tools.dros import search_dro
#from xnat_tools.scan_assessors import get_scan_asr_fname_and_id
#from xnat_tools.format_paths_dict import get_scan_asr_fname_and_id
 
from io_tools.import_dro import DroImporter
from io_tools.import_roicol import RoicollectionImporter
from io_tools.inputs_checker import are_inputs_valid
from dicom_tools.metadata import (
    get_dcm_uids, get_roicol_labels, get_roicol_nums
    )
from seg_tools.metadata import (
    get_DIVs, group_list_by_seg, get_PFFGS_to_sliceinds
    )
from seg_tools.get_seg_data_of_interest import (
    get_seg_data_of_interest, raise_error_if_no_seg_data_of_interest
    )

class Dataset:
    """
    This class imports data downloaded from XNAT.
    
    Parameters
    ----------
    params : Params Object
        Object containing various parameters.
    
    Returns
    -------
    
    """
    
    
    def __init__(self, params):
        
        
        srcDcmDir = params.srcXnatParams['dicomDir']
        
        trgDcmDir = params.trgXnatParams['dicomDir']
        
        #srcRoicolFpath = params.srcXnatParams['roicolFpath']
        
        # Import an appropriate DRO:
        self.droData = DroImporter(params)
        print(f'type(self.droData.dro) = {type(self.droData.dro)}\n')
        
        # Import the Source and Target (if applicable) ROI Collections:
        self.RoicollectionImporter(params)
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
        self.get_dicom_uids(params, 'src')
        self.get_dicom_uids(params, 'trg')
        
        
        # Get the RTS/SEG data of interest:
        
        
        
        #self.import_dicoms_and_roicols(self)
    
    
    def get_dicom_uids(self, params, srcORtrg):
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
        
    def get_seg_metadata(self, params, srcORtrg):
        """
        Get metadata for a SEG.
        
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
                self.srcDIVs = get_DIVs(self.srcSeg)
                
                # Get list of PerFrameFunctionalGroupsSequence-to-slice indices
                # and the same grouped by segment:
                self.srcAllF2Sinds = get_PFFGS_to_sliceinds(
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
                        segNames=self.srcRoiNames, 
                        p2c=params.genParams['p2c'])
                
                # Raise exception if self.srcF2SindsBySeg is empty (i.e. if 
                # there were no segmentations that matched the input parameters)
                raise_error_if_no_seg_data_of_interest(
                    f2sIndsBySeg=self.srcF2SindsBySeg, 
                    segLabs=self.srcRoiNames, 
                    slcNum=params.srcXnatParams['slcNum'])
        else:
            #seg = self.trgSeg
        
        
        
        
    
    def import_dicoms_and_roicols(self):
        """
        Imports scans and ROI Collections.

        Returns
        -------
        None.

        """
        
        pathsDict = self.pathsDict
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        
        
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
        
        
        
        