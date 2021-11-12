# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:47:21 2021

@author: ctorti
"""

"""
from importlib import reload
#import xnat_tools.dros
#reload(xnat_tools.dros)
#import io_tools.fetch_params_and_data
#reload(io_tools.fetch_params_and_data)
import io_tools.download_data
reload(io_tools.download_data)
import io_tools.imports
reload(io_tools.imports)
import io_tools.import_dro
reload(io_tools.import_dro)
#import io_tools.import_roicol
#reload(io_tools.import_roicol)
import io_tools.inputs_checker
reload(io_tools.inputs_checker)
import dicom_tools.seg_metadata
reload(dicom_tools.seg_metadata)
import dicom_tools.seg_data
reload(dicom_tools.seg_data)
import dicom_tools.rts_metadata
reload(dicom_tools.rts_metadata)
import dicom_tools.rts_data
reload(dicom_tools.rts_data)
import dicom_tools.dcm_metadata
reload(dicom_tools.dcm_metadata)
import conversion_tools.pixarrs_ims
reload(conversion_tools.pixarrs_ims)
import conversion_tools.inds_pts_pixarrs
reload(conversion_tools.inds_pts_pixarrs)
#import conversion_tools.inds_pts_cntdata
#reload(conversion_tools.inds_pts_cntdata)
"""

#from xnat_tools.sessions import create_session
#from io_tools.fetch_params_and_data import Params
#from xnat_tools.dros import search_dro
#from xnat_tools.scan_assessors import get_scan_asr_fname_and_id
#from xnat_tools.format_paths_dict import get_scan_asr_fname_and_id

from pydicom import dcmread

#from io_tools.import_dro import DroImporter
#from io_tools.import_roicol import RoicollectionImporter
from io_tools.imports import import_dcms, import_dicoms_as_im
#from io_tools.inputs_checker import are_inputs_valid, which_use_case
from dicom_tools.dcm_metadata import (
    get_dcm_uids, get_roicol_labels, get_roicol_nums
    )
from dicom_tools.seg_metadata import (
    get_DIVs, group_list_by_seg, get_f2sInds
    )
from dicom_tools.seg_data import (
    get_seg_data_of_interest, raise_error_if_no_seg_data_of_interest
    )
from dicom_tools.rts_metadata import (
    get_ptsByCntByRoi
    )
from dicom_tools.rts_data import (
    get_rts_data_of_interest, raise_error_if_no_rts_data_of_interest
    )
from conversion_tools.pixarrs_ims import pixarrBySeg_to_labimBySeg
from conversion_tools.inds_pts_pixarrs import ptsByCntByRoi_to_pixarrByRoi
#from conversion_tools.inds_pts_cntdata import ptsByCntByRoi_to_cntdataByCntByRoi
from image_tools.attrs_info import get_im_attrs
from general_tools.geometry import get_im_extent
from io_tools.exports import export_im

class DataImporter:
    # TODO modify the docstrings
    """
    This class imports data downloaded from XNAT and assigns them as instance
    variables.
    
    Parameters
    ----------
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
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
        self.cfgDict = cfgDict
        self.srcORtrg = srcORtrg
        
        if srcORtrg == 'src':
            self.dicomDir = cfgDict['srcDicomDir']
            self.roicolFpath = cfgDict['srcRoicolFpath']
            self.slcNum = cfgDict['srcSlcNum']
        else:
            self.dicomDir = cfgDict['trgDicomDir']
            self.roicolFpath = cfgDict['trgRoicolFpath']
            self.slcNum = cfgDict['trgSlcNum']
        
        # Initialise attributes that will be updated:
        self.divs = []
        self.allF2Sinds = []
        self.allF2SindsBySeg = []
        self.roiNums = []
        self.f2sIndsBySeg = []
        self.pixarrBySeg = []
        self.labimBySeg = []
        self.allC2Sinds = []
        self.allC2SindsByRoi = []
        self.c2sIndsByRoi = []
        self.ptsByCntByRoi = []
        self.cntdataByCntByRoi = []
        self.f2sIndsByRoi = []
        self.pixarrByRoi = []
        self.labimByRoi = []
        
        #self.timings = params.timings
        #self.timingMsgs = params.timingMsgs
    
    def get_dicom_uids(self):
        """
        Get StudyUID, SeriesUID, FrameOfReferenceUID and SOPInstanceUIDs for
        a DICOM series.
        
        Parameters
        ----------
        self.dicomDir : str
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
            = get_dcm_uids(self.dicomDir)
    
    def import_roicol(self):
        """
        Imports an ROI Collection as a Pydicom Object.
        
        If self.roicolFpath = None (e.g. SEG data may not exist for Target, 
        hence self.roicolFpath will be None), the roicol attribute will be also
        be None.
        
        Parameters
        ----------
        self.roicolFpath : str or None
            File path to DICOM-RTSTRUCT or DICOM-SEG file. None if a ROI
            collection doesn't exist.
        
        Returns
        -------
        self.roicol : Pydicom Object or None
        """
        if self.roicolFpath != None:
            self.roicol = dcmread(self.roicolFpath)
        else:
            self.roicol = None
    
    def get_seg_metadata(self):
        """
        Get metadata for a SEG. 
        
        Since SEG data may not be provided for Target, this function will only
        get SEG data for Source or Target (via the parameter srcORtrg).
        
        Parameters
        ----------
        self.cfgDict : dict
            Dictionary containing various parameters.
        self.srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        self.divs : list of ints or None
            Dimensional Index Values
        self.allF2Sinds : list of ints or None
            Frame-to-slice indices (i.e. slice indices corresponding to each
            frame for the ROI Collection).
        self.allF2SindsBySeg : list of a list of ints or None
            A list of frame-to-slice indices for each segment in the ROI
            Collection.
        self.roiNums : list of ints or None
            List of ints for each segment in the ROI Collection (e.g.
            [0, 1] => 2 segments).
        self.pixarrBySeg : list of Numpy data array Objects or None
            List of pixel arrays - one for each segment.
        self.f2sIndsBySeg : list of a list of ints or None
            List for each segment of a list of frame-to-slice indices.  
        """
        
        cfgDict = self.cfgDict
        
        if self.srcORtrg == 'src':
            roiName = cfgDict['srcRoiName']
            #slcNum = cfgDict['srcSlcNum']
        else:
            roiName = cfgDict['trgRoiName']
            #slcNum = cfgDict['trgSlcNum']
        # Store roiName for latter use:
        self.roiName = roiName
        
        # Get metadata:
        if self.roicol != None:
            #print(self.roicol)
            #print(type(self.roicol), '\n')
            
            # Get list of DimensionIndexValues:
            self.divs = get_DIVs(self.roicol)
            
            # Get list of PerFrameFunctionalGroupsSequence-to-slice indices
            # (referred to as frame-to-slice indices, or f2sInds), and the same
            # grouped by segment:
            self.allF2Sinds = get_f2sInds(
                seg=self.roicol, 
                sopuids=self.sopuids
                )

            self.allF2SindsBySeg = group_list_by_seg(
                listToGroup=self.allF2Sinds, 
                divs=self.divs
                )
            
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
            self.pixarrBySeg, self.f2sIndsBySeg = get_seg_data_of_interest(
                seg=self.roicol, 
                allF2SindsBySeg=self.allF2SindsBySeg, 
                segNums=self.roiNums, 
                allSegLabs=self.roiNames,
                segLab=self.roiName,
                slcNum=self.slcNum,
                p2c=cfgDict['p2c']
                )
            
            # Raise exception if self.f2sIndsBySeg is empty (i.e. if 
            # there were no segmentations that matched the input parameters)
            raise_error_if_no_seg_data_of_interest(
                f2sIndsBySeg=self.f2sIndsBySeg, 
                allSegLabs=self.roiNames, 
                slcNum=self.slcNum
                )
            
            # Get list of label images-by-ROI:
            """ 
            Note this won't be used for use cases 1-2 but is needed for 3-5.
            """
            self.labimBySeg, self.f2sIndsBySeg = pixarrBySeg_to_labimBySeg(
                pixarrBySeg=self.pixarrBySeg, 
                f2sIndsBySeg=self.f2sIndsBySeg, 
                refIm=self.dcmIm,
                p2c=cfgDict['p2c']
                )
    
    def get_rts_metadata(self):
        # TODO update docstrings
        """
        Get metadata for a RTSTRUCT. 
        
        Since RTSTRUCT data may not be provided for Target, this function will
        only get RTSTRUCT data for Source or Target (via the parameter 
        srcORtrg).
        
        Parameters
        ----------
        self.cfgDict : dict
            Dictionary containing various parameters.
        self.srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        self.divs : list of ints or None
            Dimensional Index Values
        self.allC2Sinds : list of ints or None
            Contour-to-slice indices (i.e. slice indices corresponding to each
            contour for the ROI Collection).
        self.allC2SindsByRoi : list of a list of ints or None
            A list of contour-to-slice indices for each ROI in the ROI 
            Collection.
        self.roiNums : list of ints or None
            List of ints for each ROI in the ROI Collection (e.g.
            [0, 1] => 2 ROIs).
        self.c2sIndsByRoi : list of a list of ints or None
            List for each ROI of a list of contour-to-slice indices.  
        """
        
        cfgDict = self.cfgDict
        p2c = cfgDict['p2c']
        
        if self.srcORtrg == 'src':
            roiName = cfgDict['srcRoiName']
            #slcNum = cfgDict['srcSlcNum']
        else:
            roiName = cfgDict['trgRoiName']
            #slcNum = cfgDict['trgSlcNum']
        # Store roiName for latter use:
        self.roiName = roiName
        
        # Get metadata:
        if self.roicol != None:
            #print(self.roicol)
            #print(type(self.roicol), '\n')
            
            """
            # Get list of DimensionIndexValues:
            self.divs = get_DIVs(self.roicol)
            
            # Get list of ContourSequence-to-slice indices
            # (referred to as frame-to-slice indices, or f2sInds for
            # interchangability with SEG data), and the same
            # grouped by segment:
            self.allF2Sinds = get_f2sInds(
                seg=self.roicol, 
                sopuids=self.sopuids
                )

            self.allF2SindsByRoi = group_list_by_seg(
                listToGroup=self.allF2Sinds, 
                divs=self.divs
                )
            """
            
            #self.allF2SindsByRoi = get_c2sIndsByRoi(
            #    rts=self.roicol, 
            #    sopuids=self.sopuids
            #    )
            
            self.allPtsByCntByRoi, self.allCntdataByCntByRoi,\
                self.allC2SindsByRoi, self.allC2Sinds = \
                    get_ptsByCntByRoi(
                        rts=self.roicol, 
                        sopuids=self.sopuids,
                        p2c=p2c
                        )
            
            # Get the segment number(s) that match the ROI name of
            # interest (roiName):
            """ 
            e.g. 0 if the ROI of interest is the only ROI, or if it
            is the first of several.
            
            Using the attribute srcRoiNums rather than srcSegNums for 
            interchangability for ROI data.
            """
            self.roiNums = get_roicol_nums(
                roicol=self.roicol,
                searchStr=self.roiName)
            
            self.roiNames = get_roicol_labels(roicol=self.roicol)
            
            # Get the list of points by contour by ROI and the list of 
            # contour-to-slice indices by ROI for the data of interest (to be 
            # copied):
            self.ptsByCntByRoi, self.cntdataByCntByRoi, self.c2sIndsByRoi =\
                get_rts_data_of_interest(
                    rts=self.roicol,
                    allPtsByCntByRoi=self.allPtsByCntByRoi,
                    allC2SindsByRoi=self.allC2SindsByRoi,
                    roiNums=self.roiNums, 
                    allRoiNames=self.roiNames,
                    roiName=self.roiName,
                    slcNum=self.slcNum,
                    p2c=p2c
                    )
            
            # Raise exception if self.c2sIndsByRoi is empty (i.e. if 
            # there were no contours that matched the input parameters)
            raise_error_if_no_rts_data_of_interest(
                c2sIndsByRoi=self.c2sIndsByRoi, 
                allRoiNames=self.roiNames, 
                slcNum=self.slcNum
                )
    
    def convert_rts_to_seg_metadata(self):
        # TODO update docstrings
        """
        Get SEG metadata from RTSTRUCT metadata. 
        
        Parameters
        ----------
        self.cfgDict : dict
            Dictionary containing various parameters.
        self.srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        self.pixarrBySeg : list of Numpy data array Objects or None
            List of pixel arrays - one for each ROI/segment.
        self.labimBySeg : list of SimpleITK Images or None
            List for each ROI/segment of the SimpleITK Image representation of
            the 3D pixel array.
        self.f2sIndsBySeg : list of list of ints
            A list (for each ROI/segment) of a list (for each frame) of the 
            slice numbers that correspond to each frame in the label image
            representation of the pixel arrays. 
        """
        
        if self.roicol != None:
            # Convert the RTS data of interest to a list of 3D pixel arrays for
            # each ROI:
            """
            Note: 
            
            Conversion of contour data to pixel data and label images are
            not required for use cases 1, 2a or 2b but for consistency with
            SEG data those class attributes will be populated.
            
            The frames in pixarrByRoi are conversions of contours in the same
            order of the points in ptsByCntByRoi, and they correspond to the
            DICOM slice indices in c2sIndsByRoi (which are not necessarily in
            a sorted order!). Hence the frames in pixarrByRoi are also not
            necessarily in a sorted (e.g. ascending) order.
            
            The conversion from ptsByCnt to pixarr preserves info on th
            presence of multiple contours for any given slice, since each
            frame in pixarr corresponds to a contour in ptsByCnt.
            """
            self.pixarrByRoi = ptsByCntByRoi_to_pixarrByRoi(
                ptsByCntByRoi=self.ptsByCntByRoi, 
                refIm=self.dcmIm, 
                p2c=self.cfgDict['p2c']
                )
            
            # Get list of label images-by-ROI:
            """
            Note:
            
            This won't be used for use cases 1-2 but is needed for 3-5.
            
            The conversion from pixarr to labim results in "loss" of knowledge
            of the pixel data that makes up each frame in labim. There may be
            cases where multiple segmentations (or contours) existed for the 
            same slice index (repeated indices in f2sIndsBySeg or c2sIndsByRoi).
            Hence if converting from labim back to pixarr using 
            imBySeg_to_pixarrBySeg(), the list of f2sIndsBySeg returned from
            that function may have fewer items than the original list of
            f2sIndsBySeg (or c2sIndsByRoi) that relate to the pixel data parsed
            from the SEG (or converted from contour data parsed from the RTS).
            """
            self.labimByRoi, self.f2sIndsByRoi = pixarrBySeg_to_labimBySeg(
                pixarrBySeg=self.pixarrByRoi, 
                f2sIndsBySeg=self.c2sIndsByRoi, 
                refIm=self.dcmIm,
                p2c=self.cfgDict['p2c']
                )
    
    def import_dicoms(self):
        """
        Imports DICOM series as a list of Pydicom Objects.
        
        Parameters
        self.dicomDir : str
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
        
        dicomDir = self.dicomDir
        p2c = self.cfgDict['p2c']
        
        self.dicoms = import_dcms(dicomDir, p2c)
    
    def get_image_attributes_from_dicoms(self, params):
        """
        Get image size, voxel spacings, image positions (IPPs), and direction 
        cosines from a DICOM series.
        
        Parameters
        ----------
        dicomDir : str
            Directory containing DICOMs.
        package : str, optional
            Package to use; acceptable inputs are:
             - 'pydicom' (default)
             - 'sitk' (SimpleITK)
            Default value is 'pydicom'.
        p2c : bool, optional
            If True results will be printed to the console. Default value is
            False.
        
        Returns
        -------
        self.imSize : list of int
            The size/dimensions of the 3D image along x, y and z, 
            e.g. [c, r, s] 
            where c = number of columns
                  r = number of rows
                  s = number of slices
        self.imSpacings : list of float
            The pixel spacings along x, y and z (= SliceThickness appended to 
            PixelSpacing), e.g. [di, dj, dk]
        self.imSlcThick : float
            The slice thickness.
        self.imPositions : list of list of float
            The ImagePositionPatient of all slices in the DICOM series, 
            e.g. [[x0_0, y0_0, z0_0], [x1_0, y1_0, z1_0], ...]
        self.imDirections : list of float
            The direction cosine along x (rows), y (columns) and z (slices).
        self.imWarnings : list of str
            List of any warnings encountered during parsing of attributes.
              
        Note
        ----
        See image_tools.attrs_info.py.
        """
        
        dicomDir = self.dicomDir
        p2c = self.cfgDict['p2c']
        
        self.imSize, self.imSpacings, self.imSlcThick, self.imPositions,\
            self.imDirections, self.imWarnings = get_im_attrs(
                dicomDir, package='pydicom', p2c=p2c
                )
    
    def import_dicom_image(self):
        """
        Imports a DICOM series as a 3D SimpleITK Image and a 3D Numpy pixel 
        array.
        
        Parameters
        ----------
        self.dicomDir : str
            Path to a directory containing a DICOM series.
            
        Returns
        -------
        self.dcmIm : SimpleITK Image
            SimpleITK Image representation of the DICOM series.
        self.dcmPixarr : SimpleITK Image
            SimpleITK Image representation of the DICOM series.
        """
        
        self.dcmIm, self.dcmPixarr = import_dicoms_as_im(self.dicomDir)
    
    def export_dicom_image(self, params):
        """
        Exports a SimpleITK Image.
        
        This is useful for exporting of a 3D DICOM image without importing of
        the full dataset, and prior to any further processing.
        
        Parameters
        ----------
        self.dcmIm : SimpleITK Image
            SimpleITK Image representation of the DICOM series.
        params : Params Object
            Object containing various parameters.
            
        Returns
        -------
        None
        
        Note
        ----
        The image can be exported in HDF, NIFTI or NRRD formats.
        """
        
        srcORtrg = self.srcORtrg
        image = self.dcmIm
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        imExportDir = cfgDict['imExportDir']
        
        print(f'* Exporting {srcORtrg} image..\n')
        
        fname = f'{runID}_{srcORtrg}Im'
        
        export_im(
            image, filename=fname,
            fileFormat='HDF5ImageIO', exportDir=imExportDir
        )
        
    def import_data(self, params):
        # TODO update docstrings
        """
        Import various data.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        Returns
        -------
        
        """
        
        cfgDict = self.cfgDict
        roicolMod = cfgDict['roicolMod']
        p2c = cfgDict['p2c']
        
        if self.srcORtrg == 'src':
            toImport = 'source'
        else:
            toImport = 'target'
            
        timingMsg = f"* Importing the {toImport} DICOM scans and ROI "\
            + "Collection...\n"
        params.add_timestamp(timingMsg)
        
        # Import an appropriate DRO:
        #self.droData = DroImporter(params)
        #print(f'type(self.droData.dro) = {type(self.droData.dro)}\n')
        
        # TODO move DRO somewhere else since DroImporter requires params for
        # Source and Target (06/08/21)
        
        # Get DICOM UIDs:
        self.get_dicom_uids()
        
        # Import the ROI Collection:
        self.import_roicol()
        if p2c:
            print(f'type(self.roicol) = {type(self.roicol)}\n')
        
        """ This is already done in get_seg_metadata (23/08/21): """
        # Get ROI Names / segment labels (required elsewhere):
        """ 
        Note: roiNames may refer to either ROI Names (for RTS data) or SEG 
        labels (for SEG data). It does not refer to the name of the ROI 
        Collection itself.
        """
        #self.roiNames = get_roicol_labels(self.roicol)
        
        """ This is already done in get_seg_metadata (23/08/21): """
        # The list of ROI/segment numbers for the ROI(s)/segment(s) of interest:
        #self.roiNums = get_roicol_nums(self.roicol, self.roiName)
        
        # Import a DICOM series as a list of Pydicom Objects:
        self.import_dicoms()
        
        # Import a DICOM series as a SimpleITK Image:
        self.import_dicom_image()
        
        # Get the image attributes:
        self.get_image_attributes_from_dicoms(params)
        
        # Get the image extent:
        self.imExtent = get_im_extent(self.dcmIm)
        
        # Get the RTS/SEG data of interest:
        if roicolMod == 'SEG':
            self.get_seg_metadata()
        else:
            self.get_rts_metadata()
            
            self.convert_rts_to_seg_metadata()
        
        #print(f'type(self.labimByRoi) = {type(self.labimByRoi)}')
        #if not self.labimByRoi is None:
        #    print(f'type(self.labimByRoi[0]) = {type(self.labimByRoi[0])}')
        
        timingMsg = f"Took [*] to import the {toImport} DICOM scans and "\
                + "ROI Collection.\n"
        params.add_timestamp(timingMsg)