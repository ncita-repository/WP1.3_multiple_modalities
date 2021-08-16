# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:47:21 2021

@author: ctorti
"""

"""
21/07:
    I'm confused by the incomplete code here to import ROI Collections 
    and import_roicol in other tab. Seems I've changed my mind and there 
    seems to be some redundancy.
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
import seg_tools.seg_data
reload(seg_tools.seg_data)
import dicom_tools.imports
reload(dicom_tools.imports)
import image_tools.imports
reload(image_tools.imports)
import conversion_tools.pixarrs_ims
reload(conversion_tools.pixarrs_ims)

#from xnat_tools.sessions import create_session
#from io_tools.fetch_params_and_data import Params
#from xnat_tools.dros import search_dro
#from xnat_tools.scan_assessors import get_scan_asr_fname_and_id
#from xnat_tools.format_paths_dict import get_scan_asr_fname_and_id

from pydicom import dcmread

from io_tools.import_dro import DroImporter
#from io_tools.import_roicol import RoicollectionImporter
from io_tools.inputs_checker import are_inputs_valid
from dicom_tools.metadata import (
    get_dcm_uids, get_roicol_labels, get_roicol_nums
    )
from seg_tools.metadata import (
    get_divs, group_list_by_seg, get_p2sInds
    )
from seg_tools.seg_data import (
    get_seg_data_of_interest, raise_error_if_no_seg_data_of_interest
    )
from conversion_tools.pixarrs_ims import pixarrByRoi_to_imByRoi
from dicom_tools.imports import import_dcms
from image_tools.imports import import_im

class DataImporter:
    # TODO modify the docstrings
    """
    This class imports data downloaded from XNAT and assigns them as instance
    variables.
    
    Parameters
    ----------
    params : DataDownloader Object
        Contains parameters and file paths in cfgDict and pathsDict.
    srcORtrg : str
        'src' or 'trg' for source or target dataset to be imported.
    
    Returns
    -------
    self.droData : droData Object
        Object contains sub-attributes (see fetch_dro in 
        io_tools/import_dro.py).
    self.srcRoicol or self.trgRoicol : Pydicom Object or None
            Pydicom representation of the ROI Collection (RTS/SEG) if not None.
    self.srcDIVs or self.trgDIVs : list of ints or None
            Dimensional Index Values
    self.srcAllF2Sinds or self.trgAllF2Sinds : list of ints or None
        Frame-to-slice indices (i.e. slice indices corresponding to each
        frame for the ROI Collection).
    self.srcAllF2SindsByRoi or self.trgAllF2SindsByRoi : list of a list of
        ints or None
        A list of frame-to-slice indices for each ROI/segment in the ROI
        Collection.
    self.srcRoiNums or self.trgRoiNums : list of ints or None
        List of ints for each ROI/segment in the ROI Collection (e.g.
        [0, 1] => 2 ROIs/segments).
    self.srcPixArrByRoi or self.trgPixArrByRoi : list of Numpy data array 
        Objects or None
        List of pixel arrays - one for each ROI/segment.
    self.srcF2SindsByRoi or self.trgF2SindsByRoi : list of a list of ints
        or None
        List for each ROI/segment of a list of frame-to-slice indices.
    self.srcIm or self.trgIm : SimpleITK Image
        SimpleITK Image representation of the DICOM series.
    self.srcDcms or self.trgDcms : list of Pydicom Objects
        A list (for each DICOM) of Pydicom representations of the DICOM series.
    self.srcStudyUID or self.trgStudyUID : str
        DICOM Study UID.
    self.srcSeriesUID or self.trgSeriesUID : str
        DICOM Series UID.
    self.srcFORUID or self.trgFORUID : str
        DICOM Frame of reference UID.
    self.srcSOPUID or self.trgSOPUID : str
        DICOM SOP UID.
    """
    
    def __init__(self, params, srcORtrg):
        
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
        self.import_image(dicomDir)
        
        # Get the RTS/SEG data of interest:
        self.get_seg_metadata(cfgDict, srcORtrg)
        
        #print(f'type(self.labimByRoi) = {type(self.labimByRoi)}')
        #if not self.labimByRoi is None:
        #    print(f'type(self.labimByRoi[0]) = {type(self.labimByRoi[0])}')
    
    def __init__210809(self, xnatParams, genParams):
        
        # TODO move DRO somewhere else since DroImporter requires params for
        # Source and Target (06/08/21)
        
        # Import an appropriate DRO:
        #self.droData = DroImporter(params)
        #print(f'type(self.droData.dro) = {type(self.droData.dro)}\n')
        
        # Get DICOM UIDs:
        self.get_dicom_uids(xnatParams['dicomDir'])
        
        # Import the ROI Collection:
        self.import_roicol(xnatParams['roicolFpath'])
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
        
        # TODO continue modifying code below to deal with only one dataset 
        # (source or target), and necessary modification to functions to allow
        # for input = None:
        
        # Get the RTS/SEG data of interest:
        self.get_seg_metadata(xnatParams, genParams['p2c'])
        
        # Import a DICOM series as a list of Pydicom Objects:
        self.get_dicoms(self, xnatParams['dicomDir'], genParams['p2c'])
        
        # Import a DICOM series as a SimpleITK Image:
        self.get_image(self, xnatParams['dicomDir'])
    
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
    
    def get_dicom_uids(self, fdir):
        """
        Get StudyUID, SeriesUID, FrameOfReferenceUID and SOPInstanceUIDs for
        a DICOM series.
        
        Parameters
        ----------
        fdir : str
            Directory containing a series of DICOM files.
        
        Returns
        -------
        self.studyuid : str
            DICOM Study UID of the DICOM series.
        self.seriesuid : str
            DICOM Series UID of the DICOM series.
        self.foruid : str
            DICOM Frame of reference UID of the DICOM series.
        self.sopuids : list of str
            List of DICOM SOP UIDs for each DICOM in the series.
        """
        self.studyuid, self.seriesuid, self.foruid, self.sopuids\
            = get_dcm_uids(fdir)
    
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
    
    def import_roicol(self, fpath):
        """
        Imports an ROI Collection as a Pydicom Object.
        
        If fpath = None (e.g. SEG data may not exist for Target, hence fpath
        will be None), the roicol attribute will be also be None.
        
        Parameters
        ----------
        fpath : str or None
            File path to DICOM-RTSTRUCT or DICOM-SEG file. None if a ROI
            collection doesn't exist.
        
        Returns
        -------
        self.roicol : Pydicom Object or None
        """
        if fpath != None:
            self.roicol = dcmread(fpath)
        else:
            self.roicol = None
        
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
    
    def get_seg_metadata(self, cfgDict, srcORtrg):
        """
        Get metadata for a SEG. 
        
        Since SEG data may not be provided for Target, this function will only
        get SEG data for Source or Target (via the parameter srcORtrg).
        
        Parameters
        ----------
        cfgDict : dict
            Dictionary containing various parameters.
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        self.divs : list of ints or None
            Dimensional Index Values
        self.allF2Sinds : list of ints or None
            Frame-to-slice indices (i.e. slice indices corresponding to each
            frame for the ROI Collection).
        self.allF2SindsByRoi : list of a list of ints or None
            A list of frame-to-slice indices for each ROI/segment in the ROI
            Collection.
        self.roiNums : list of ints or None
            List of ints for each ROI/segment in the ROI Collection (e.g.
            [0, 1] => 2 ROIs/segments).
        self.pixarrByRoi : list of Numpy data array Objects or None
            List of pixel arrays - one for each ROI/segment.
        self.f2sIndsByRoi : list of a list of ints or None
            List for each ROI/segment of a list of frame-to-slice indices.  
        """
        
        if srcORtrg == 'src':
            roiName = cfgDict['srcRoiName']
            #slcNum = cfgDict['srcSlcNum']
        else:
            roiName = cfgDict['trgRoiName']
            #slcNum = cfgDict['trgSlcNum']
        # Store roiName for latter use:
        self.roiName = roiName
        
        # Get metadata:
        if self.roicol == None:
            self.divs = None
            self.allF2Sinds = None
            self.allF2SindsByRoi = None
            self.roiNums = None
            self.pixarrByRoi = None
            self.f2sIndsByRoi = None
            self.labimByRoi = None
        else:
            #print(self.roicol)
            print(type(self.roicol), '\n')
            
            # Get list of DimensionIndexValues:
            self.divs = get_divs(self.roicol)
            
            # Get list of PerFrameFunctionalGroupsSequence-to-slice indices
            # and the same grouped by segment:
            self.allF2Sinds = get_p2sInds(
                seg=self.roicol, 
                sopuids=self.sopuids)

            self.allF2SindsByRoi = group_list_by_seg(
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
                roicol=self.roicol, 
                searchStr=self.roiName)
            
            self.roiNames = get_roicol_labels(roicol=self.roicol)
            
            # Get the pixel array by segment and frame-to-slice indices by
            # segment for the data of interest (to be copied):
            self.pixarrByRoi, self.f2sIndsByRoi = get_seg_data_of_interest(
                seg=self.roicol, 
                allF2SindsBySeg=self.allF2SindsByRoi, 
                segNums=self.roiNums, 
                segLabs=self.roiNames,
                segLab=self.roiName,
                slcNum=self.slcNum,
                p2c=cfgDict['p2c']
                )
            
            # Raise exception if self.f2sIndsByRoi is empty (i.e. if 
            # there were no segmentations that matched the input parameters)
            raise_error_if_no_seg_data_of_interest(
                f2sIndsBySeg=self.f2sIndsByRoi, 
                segLabs=self.roiNames, 
                slcNum=self.slcNum
                )
            
            # Get list of label images-by-ROI:
            """ 
            Note this won't be used for use cases 1-2 but is needed for 3-5.
            """
            self.labimByRoi = pixarrByRoi_to_imByRoi(
                pixarrByRoi=self.pixarrByRoi, 
                f2sIndsByRoi=self.f2sIndsByRoi, 
                refIm=self.image,
                p2c=cfgDict['p2c']
                )
    
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
    
    def import_dicoms(self, dpath, p2c=False):
        """
        Imports DICOM series as a list of Pydicom Objects.
        
        Parameters
        dpath : str
            Path to a directory containing a DICOM series.
        p2c : bool (optional)
            If True some results will be printed to the console. False by 
            default.
        
        Returns
        -------
        self.dicoms : list of Pydicom Objects
            A list (for each DICOM) of Pydicom representations of the DICOM
            series.
        """
        
        self.dicoms = import_dcms(dicomDir=dpath, p2c=p2c)
    
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
    
    def import_image(self, dpath):
        """
        Imports a DICOM series as a SimpleITK Image.
        
        Parameters
        ----------
        dpath : str
            Path to a directory containing a DICOM series.
            
        Returns
        -------
        self.image : SimpleITK Image
            SimpleITK Image representation of the DICOM series.
        """
        
        self.image = import_im(dicomDir=dpath)
    
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
        
        self.srcIm = import_im(dicomDir=srcDcmDir)
        self.trgIm = import_im(dicomDir=trgDcmDir)
    