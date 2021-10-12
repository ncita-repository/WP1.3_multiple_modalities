# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 10:39:55 2021

@author: ctorti
"""


from importlib import reload



class Resampler:
    # TODO modify the docstrings
    """
    This class resamples label images in srcDataset (Importer Object) based on 
    user settings (DataDownloader Object) and trgDataset, i.e. reference,
    (Importer Object).
    
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
    
    def __init__(self, params, srcDataset, trgDataset):
        
        #cfgDict = params.cfgDict
        #regTxName = cfgDict['regTxName']
        #p2c = cfgDict['p2c']
        
        # Copy selected instance attributes from trgDataset:
        """
        Resampled srcDataset (resSrcDataset) will acquire some of trgDataset's
        attributes. 
        """
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
        self.seriesuid = trgDataset.seriesuid
        self.slcNum = trgDataset.slcNum
        #self.sopuids = trgDataset.sopuids
        self.studyuid = trgDataset.studyuid
        
        # resSrcDataset will have same/similar ROI Names as srcDataset:
        # TODO consider modifying name to reflect the fact that it has been
        # resampled?
        self.roiName = srcDataset.roiName
        self.roiNames = srcDataset.roiNames
        self.roiNums = srcDataset.roiNums
        
        
        # Copy srcDataset and trgDataset:
        #self.srcDataset = srcDataset
        #self.trgDataset = trgDataset
        
        self.resample_labimByRoi(params, srcDataset, trgDataset)