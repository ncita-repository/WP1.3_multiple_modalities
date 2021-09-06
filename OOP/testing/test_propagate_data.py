# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:52:32 2021

@author: ctorti
"""


from importlib import reload

#import testing.test_import_data
#reload(testing.test_import_data)
#from testing.test_import_data import FetcherTester

import testing.test_import_data
reload(testing.test_import_data)

import io_tools.import_dro
reload(io_tools.import_dro)

import image_tools.propagate
reload(image_tools.propagate)

import dro_tools.create_dro
reload(dro_tools.create_dro)

import seg_tools.create_seg
reload(seg_tools.create_seg)

import plotting_tools.plotting
reload(plotting_tools.plotting)

import io_tools.exports
reload(io_tools.exports)

import os
import time
from pathlib import Path
import SimpleITK as sitk
from testing.test_import_data import DataImporterTester
from io_tools.import_dro import DroImporter
from image_tools.propagate import Propagator
from dro_tools.create_dro import DroCreator
from seg_tools.create_seg import create_seg
#from plotting_tools.plotting import plot_pixarrs_from_list_of_segs_v3
from plotting_tools.plotting import (
    plot_pixarrs_from_list_of_segs_and_images, compare_res_results,
    plot_metricValues_v_iters
    )
from io_tools.exports import (export_im, export_list_to_txt, export_newRoicol)

class PropagatorTester:
    # TODO update docstrings
    """
    Run unit or integrated test of propagate.
    
    To run unit test, provide Importer Objects for source and target. To run 
    integrated test, leave optional input arguments blank.
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
        
    Returns
    -------
    newDataset : Dataset Object
        New Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_propagate_data import PropagatorTester
        
        testResult = PropagatorTester(
            runID, params, srcDataset, trgDataset
            )
    
    To integrated test:
        from testing.test_propagate_data import PropagatorTester
        
        testResult = PropagatorTester(runID)
    """
    
    def __init__(
            self, cfgDir, runID, srcDataset=None, trgDataset=None, params=None
            ):
        if (params == None or srcDataset == None or trgDataset == None):
            # Import data for source:
            print('Running ImporterTester...\n')
            testResult = DataImporterTester(cfgDir=cfgDir, runID=runID)
            
            params = testResult.params
            srcDataset = testResult.srcDataset
            trgDataset = testResult.trgDataset
        
        self.params = params
        self.srcDataset = srcDataset
        self.trgDataset = trgDataset
        """ 
        useCaseToApply isn't an attribute until the Propagator class is
        implemented
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        """
        
        
        #if '5' in params.cfgDict['useCaseToApply']:
        #    """
        #    Image registration will be required unless a suitable DRO can be
        #    found on XNAT.  Search XNAT for an appropriate DRO:
        #    """
        # The following was under the if statement:
        
        if params.cfgDict['useDroForTx']:
            # Instantiate a droData Object:
            droData = DroImporter()
            droData.fetch_dro(params)
            self.droData = droData
            dro = droData.dro
            #print(f'type(droData.dro) = {type(droData.dro)}\n')
        else:
            dro = None
        
        #print(f'dir(srcDataset) = {dir(srcDataset)}\n')
        
        # Instantiate a Propagator object:
        newDataset = Propagator(srcDataset, trgDataset, params)
        
        # Propagate/copy the source ROI Collection:
        newDataset.propagate(srcDataset, trgDataset, params, dro)
        
        self.newDataset = newDataset
        # Update params:
        params = self.params
        useCaseToApply = params.cfgDict['useCaseToApply']
        
        """ Creating a new ROI Collection and DRO (if applicable) is not
        required but adding it to more fully test (23/08/2021). """
        newRoicol = create_seg(srcDataset, trgDataset, newDataset, params)
        
        self.newDataset.roicol = newRoicol
        
        """ Note create_seg will need to be wrapped in another function along
        with the yet-to-be-refactored create_rts. """
        # TODO refactor create_rts and wrap it and create_seg into a new func (create_new_roicol)
        
        """ Export the new ROI Collection """
        self.export_newRoiCol()
        
        """ Create and export a new DRO (if applicable) """
        if '5' in useCaseToApply:
            # Instantiate a DroCreator Object:
            newDroData = DroCreator(newDataset, params)
            
            # Create the new DRO:
            newDroData.create_dro(srcDataset, trgDataset, newDataset, params)
            
            # Export the new DRO:
            newDroData.export_dro(params)
            
            self.newDroData = newDroData
        
    
        """ Export label images, tranforms and transform parameters """
        self.export_labims()
        self.export_tx_and_params()
        
        """ Plot results """
        self.plot_metric_v_iters()
        self.plot_res_results()
        self.plot_roi_over_dicom_ims()
    
    def export_newRoiCol(self):
        
        print('* Exporting the new Target ROI Collection..\n')
        
        srcDataset = self.srcDataset
        #trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        roicolMod = params.cfgDict['roicolMod']
        runID = params.cfgDict['runID']
        
        if roicolMod == 'RTSTRUCT': 
            exportDir = params.cfgDict['rtsExportDir']
        else:
            exportDir = params.cfgDict['segExportDir']
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        fname = f'{runID}_newRoicol_{currentDateTime}'
        
        newRoicolFpath = export_newRoicol(
            newRoicol=newDataset.roicol, 
            srcRoicolFpath=srcDataset.roicolFpath, 
            exportDir=exportDir, fname=fname
            )
        
        self.newDataset.roicolFpath = newRoicolFpath
    
    def plot_metric_v_iters(self):
        """ 
        Plot the metric values v iterations from the registeration.
        """
        
        print('* Plotting metric v iterations from registeration..\n')
        
        metricValues = self.newDataset.metricValues
        multiresIters = self.newDataset.multiresIters
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        exportPlot = cfgDict['exportPlots']
        #p2c = cfgDict['p2c']
        useCaseToApply = cfgDict['useCaseToApply']
        #forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        
        resExportDir = cfgDict['resPlotsExportDir']
        
        #print(useCaseToApply)
        
        # Prepare plot title:
        if useCaseToApply in ['5a', '5b'] and not useDroForTx:
            resTitle = f'Metric v iters for Src reg to Trg ({regTxName}, '\
                + f'{initMethod}, {resInterp})'
            
            # Prepare filename for exported plot:
            fname = f'{runID}_' + resTitle.replace(' ', '_').replace(',', '')
            fname = fname.replace('(', '').replace(')', '')
            
            plot_metricValues_v_iters(
                metricValues=metricValues, multiresIters=multiresIters, 
                exportPlot=exportPlot, exportDir=resExportDir,
                fname=fname
            )
    
    def plot_res_results(self):
        """ 
        Plot a single slice from trgIm and resIm and compare to assess result
        of resampling/registering.
        """
        
        print('* Plotting resampled/registered results..\n')
        
        #srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        #roicolMod = cfgDict['roicolMod']
        exportPlot = cfgDict['exportPlots']
        #p2c = cfgDict['p2c']
        useCaseToApply = cfgDict['useCaseToApply']
        #forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        trgIm = trgDataset.image
        resIm = newDataset.image # resampled source or source-registered-to-target 
        
        resExportDir = cfgDict['resPlotsExportDir']
        
        #print(useCaseToApply)
        
        # Prepare plot title for new dataset:
        if useCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
            resTitle = 'Src '
            if useCaseToApply in ['3a', '3b', '4a', '4b']:
                resTitle += f'res to Trg ({resInterp})'
            elif useCaseToApply in ['5a', '5b']:
                if useDroForTx:
                    resTitle += f'tx to Trg ({regTxName} DRO, {resInterp})'
                else:
                    resTitle += f'reg to Trg ({regTxName}, {initMethod}, '\
                        + f'{resInterp})'
            
            midInd = trgIm.GetSize()[2] // 2
            
            # Prepare filename for exported plot:
            fname = f'{runID}_' + resTitle.replace(' ', '_').replace(',', '')
            fname = fname.replace('(', '').replace(')', '')
            
            compare_res_results(
                resIm0=trgIm, resIm1=resIm, resInd=midInd,
                resTitle0='Target image', resTitle1=resTitle,
                exportPlot=exportPlot, exportDir=resExportDir,
                fname=fname
            )
    
    def plot_roi_over_dicom_ims(self):
        """ 
        Plot copied/propagated ROI overlaid on DICOM images 
        """
        
        print('* Plotting ROI over DICOM images..\n')
        
        srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        roicolMod = cfgDict['roicolMod']
        exportPlot = cfgDict['exportPlots']
        p2c = cfgDict['p2c']
        useCaseToApply = cfgDict['useCaseToApply']
        forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        #trgIm = trgDataset.image
        #resIm = newDataset.image # resampled source or source-registered-to-target 
        
        #resExportDir = cfgDict['resPlotsExportDir']
        
        if roicolMod == 'RTSTRUCT':
            roiExportDir = cfgDict['rtsPlotsExportDir']
        else:
            roiExportDir = cfgDict['segPlotsExportDir']
        
        # Prepare plot title for new dataset:
        method = 'Src '
        if useCaseToApply in ['1', '2a']:
            method += 'copied to Trg'
        elif useCaseToApply == '2b':
            method += 'propagated to Trg'
        elif useCaseToApply in ['3a', '4a']:
            method += f'copied to Trg (res, {resInterp})'
        elif useCaseToApply in ['3b', '4b']:
            method += f'propagated to Trg (res, {resInterp})'
        elif useCaseToApply == '5a':
            if useDroForTx:
                method += f'copied to Trg ({regTxName} DRO, {resInterp})'
            else:
                method += f'copied to Trg ({regTxName} reg, {initMethod}, '\
                    + f'{resInterp})'
        elif useCaseToApply == '5b':
            if useDroForTx:
                method += f'propagated to Trg ({regTxName} DRO, {resInterp})'
            else:
                method += f'propagated to Trg ({regTxName} reg, {initMethod},'\
                    + f' {resInterp})'
        
        #print(method)
        
        listOfSegs = [srcDataset.roicol, newDataset.roicol]
        """
        Note: newDataset does not have a dicomDir since DICOMs are not
        generated for the resampled/registered source dataset, but for the
        purposes of pulling metadata (sizes, spacings, IPPs, directions, etc)
        relevant to the resampled/registered dataset, getting it from the
        target dataset is equivalent.
        """
        listOfDicomDirs = [srcDataset.dicomDir, trgDataset.dicomDir]
        
        listOfImages = [srcDataset.image, newDataset.image]
        listOfPlotTitles = ['Src', method]
        
        """
        plot_pixarrs_from_list_of_segs_v3(
        listOfSegs, listOfDicomDirs, listOfPlotTitles, exportPlot=exportPlot,
            exportDir=exportDir, runID=runID, useCaseToApply=useCaseToApply,
            forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
            initMethod=initMethod, resInterp=resInterp, p2c=p2c
        )
        """
        
        # Prepare filename for exported plot:
        fname = f'{runID}_' + method.replace(' ', '_').replace(',', '')
        fname = fname.replace('(', '').replace(')', '')
            
        plot_pixarrs_from_list_of_segs_and_images(
            listOfSegs, listOfImages, listOfDicomDirs, listOfPlotTitles, 
            exportPlot=exportPlot, exportDir=roiExportDir,
            runID=runID, useCaseToApply=useCaseToApply,
            forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
            initMethod=initMethod, resInterp=resInterp, 
            #txtToAddToFname=fname, p2c=p2c
            txtToAddToFname='', p2c=p2c
        )

    def export_labims(self):
        
        print('* Exporting label images..\n')
        
        srcDataset = self.srcDataset
        trgDataset = self.trgDataset
        newDataset = self.newDataset
        params = self.params
        
        labimExportDir = params.cfgDict['labimExportDir']
        runID = params.cfgDict['runID']
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        # List of datasets:
        dsets = [srcDataset, trgDataset, newDataset]
        
        dsetNames = ['src', 'trg', 'new']
        
        # List of labimByRoi to export:
        listOfLabimByRoi = [dset.labimByRoi for dset in dsets]
        
        # List of file names (will be modified in loop below):
        fnames = [
            f'{runID}_{name}Labim_{currentDateTime}' for name in dsetNames
            ]
        
        for d in range(len(listOfLabimByRoi)):
            labimByRoi = listOfLabimByRoi[d]
            fname = fnames[d]
            
            if labimByRoi: # if not None
                # Loop through each labim (for each ROI):
                for r in range(len(labimByRoi)):
                    export_im(
                        labimByRoi[r], filename=f'{fname}{r}',
                        fileFormat='HDF5ImageIO', exportDir=labimExportDir
                    )
    
    def export_tx_and_params(self):
        
        print('* Exporting transforms and transform parameters..\n')
        
        params = self.params
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        
        txExportDir = params.cfgDict['txExportDir']
        
        # Create the export directory if it doesn't already exist:
        if not os.path.isdir(txExportDir):
            #os.mkdir(exportDir)
            Path(txExportDir).mkdir(parents=True)
        
        resTx = self.newDataset.resTx
        initRegTx = self.newDataset.initRegTx
        preRegTx = self.newDataset.preRegTx
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        # List of transforms to export and their names:
        txs = [resTx, initRegTx, preRegTx]
        txNames = ['resTx', 'initRegTx', 'preRegTx']
        
        for tx, txName in zip(txs, txNames):
            if tx: # is not None
                fname = f'{runID}_{txName}_{currentDateTime}'
                
                # Export tx as a TFM:
                fpath = os.path.join(txExportDir, fname + '.tfm')
                sitk.WriteTransform(tx, fpath)
                
                # Export the tx parameters as a TXT:
                export_list_to_txt(
                    items=tx.GetParameters(),
                    filename=fname,
                    exportDir=txExportDir
                    )
        
def test_propagator_210812(runID, params=None, srcDataset=None, trgDataset=None):
    """
    Run unit or integrated test of propagate.
    
    To run unit test, provide Importer Objects for source and target. To run 
    integrated test, leave optional input arguments blank.
    
    Parameters
    ----------
    params : Fetcher Object, optional
        A Fetcher object containing various parameters. The default is None.
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
        
    Returns
    -------
    newDataset : Dataset Object
        New Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_propagate_data import test_propagator
        
        newDataset = test_propagator(
            runID, params, srcDataset, trgDataset
            )
    
    To integrated test:
        from testing.test_propagate_data import test_propagator
        
        newDataset = test_propagator(runID)
    """
    
    if params == None:
        # Fetch the parameters:
        print('Fetching parameters (locally) and data (from XNAT)...\n')
        params = test_fetch_params_and_data(runID)
    
    if srcDataset == None:
        # Import data for source:
        print('Importing source dataset...\n')
        srcDataset = Importer(params=params, srcORtrg='src')
    
    if trgDataset == None:
        # Import data for target:
        print('Importing target dataset...\n')
        trgDataset = Importer(params=params, srcORtrg='trg')
    
    return srcDataset, trgDataset