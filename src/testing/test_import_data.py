# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:38:46 2021

@author: ctorti
"""

from importlib import reload

#import testing.test_fetch_config
#reload(testing.test_fetch_config)

import testing.test_download_data
reload(testing.test_download_data)

import io_tools.import_data
reload(io_tools.import_data)

import io_tools.exports
reload(io_tools.exports)

import plotting_tools.conversion_results
reload(plotting_tools.conversion_results)

#from testing.test_fetch_config import ConfigFetcherTester
from testing.test_download_data import DataDownloaderTester
from io_tools.import_data import DataImporter
from io_tools.exports import export_im
from plotting_tools.conversion_results import plot_mask_to_cnt_conversion


class DataImporterTester:
    """
    Run unit or integrated test of import_data.
    
    To run unit test, provide Fetcher Object. To run integrated test, leave 
    input argument blank.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    params : ConfigFetcher Object, optional
        A ConfigFetcher object. The default is None.
        
    Returns
    -------
    srcDataset : Dataset Object
        Dataset Object for Source DICOM series.
    trgDataset : Dataset Object
        Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_import_data import DataImporterTester
        
        testResult = ImporterTester(runID, params)
    
    To integrated test:
        from testing.test_import_data import DataImporterTester
        
        testResult = ImporterTester(runID)
    """
    
    def __init__(self, cfgDir, runID, params=None):
        
        if params == None:
            # Run DataDownloaderTester to get the config parameters and 
            # download the data from XNAT:
            testResult = DataDownloaderTester(cfgDir=cfgDir, runID=runID)
            params = testResult.params
            
            # Or run the internal operations within DataDownloaderTester:
            """
            # Run ConfigFetcherTester:
            cfgObj = ConfigFetcherTester(cfgDir=cfgDir, runID=runID)
            
            # Instantiate a DataDownloader object:
            params = DataDownloader(cfgDict=cfgObj.cfgDict)
            """
        
        # Initialise a DataImporter objects for source and target:
        #self.DataImporter(params=params, srcORtrg='src')
        srcDataset = DataImporter(params=params, srcORtrg='src')
        #self.DataImporter(params=params, srcORtrg='trg')
        trgDataset = DataImporter(params=params, srcORtrg='trg')
        
        # Import source and target data:
        print('Importing source dataset...\n')
        srcDataset.import_data(params)
        print('Importing target dataset...\n')
        trgDataset.import_data(params)
        
        self.params = params
        self.srcDataset = srcDataset
        self.trgDataset = trgDataset
        
        p2c = params.cfgDict['p2c']
        roicolMod = params.cfgDict['roicolMod']
        
        if p2c and roicolMod == 'RTSTRUCT':
            # Check the result of conversion from contour to pixel data:
            listOfIms = [srcDataset.dcmIm]
            listOfDcmPixarr = [srcDataset.dcmPixarr]
            listOfF2SindsByRoi = [srcDataset.f2sIndsByRoi]
            listOfPtsByCntByRoi = [srcDataset.ptsByCntByRoi]
            listOfPixarrByRoi = [srcDataset.pixarrByRoi]
            listOfDcmDirs = [srcDataset.dicomDir]
            listOfPlotTitles = ['Conversion for srcDataset']
            
            if trgDataset.roicol != None:
                listOfIms.append(trgDataset.dcmIm)
                listOfDcmPixarr.append(trgDataset.dcmPixarr)
                listOfF2SindsByRoi.append(trgDataset.f2sIndsByRoi)
                listOfPtsByCntByRoi.append(trgDataset.ptsByCntByRoi)
                listOfPixarrByRoi.append(trgDataset.pixarrByRoi)
                listOfDcmDirs.append(trgDataset.dicomDir) 
                listOfPlotTitles.append('Conversion for srcDataset')
            
            plot_mask_to_cnt_conversion(
                listOfIms=listOfIms, 
                listOfDcmPixarr=listOfDcmPixarr, 
                listOfF2SindsByRoi=listOfF2SindsByRoi,
                listOfPtsByCntByRoi=listOfPtsByCntByRoi, 
                listOfPixarrByRoi=listOfPixarrByRoi,
                listOfDcmDirs=listOfDcmDirs, 
                listOfPlotTitles=listOfPlotTitles,
                p2c=p2c
            )

    def export_ims(self):
            
        print('* Exporting images..\n')
        
        srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        imExportDir = cfgDict['imExportDir']
        
        # List of datasets:
        dsets = [srcDataset, trgDataset]
        
        dsetNames = ['src', 'trg']
        
        # List of images to export:
        images = [dset.image for dset in dsets]
        
        # List of file names (will be modified in loop below):
        fnames = [f'{runID}_{name}Im' for name in dsetNames]
        
        for d in range(len(images)):
            export_im(
                images[d], filename=fnames[d],
                fileFormat='HDF5ImageIO', exportDir=imExportDir
            )