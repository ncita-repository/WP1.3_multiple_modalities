# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:24:09 2021

@author: ctorti
"""



""" 
Batch script for running the ROI Collection propagation.

Note:
This is intended for running in a Python console/IDE or in a command shell,
and relies on parameters stored in a configuration file. Based on the 
parameters the necessary data is fetched from XNAT using REST APIs and 
downloaded locally.  If image registration is required, a search for a 
suitable Spatial Registration Object (for rigid including affine registrations)
or Deformable Registration Object (for deformable/bspline registrations) and if
found it will be downloaded locally, unless the variable useDroForTx in the
config file is set to False.

Once the copy/propagation is complete, the new ROI Collection will be uploaded
to XNAT (via the REST API), as too will a DRO (if image registration was 
performed) unless the variable uploadDro is set to False.

Other outputs (useful for debugging) will be exported to the local drive (e.g.
transforms, images, label images, log files, plots, etc.) unless the relevant
variables (exportTx, exportIm, exportLabim, etc.) are set to False.

This implementation should cover all use case scenarios defined.

See app_ohif.py for a newer version intended for use as a backend script for
the OHIF-Viewer (still in development as of Nov 4, 2021).
""" 

import os
import sys

#code_root = r'C:\Code\WP1.3_multiple_modalities\src'
code_root = os.getcwd()
#print(f'code_root = {code_root}\n')

# Add code_root to the system path so packages can be imported from it:
sys.path.append(code_root)


from importlib import reload
import io_tools.fetch_config
reload(io_tools.fetch_config)
import io_tools.download_data
reload(io_tools.download_data)
import io_tools.import_data
reload(io_tools.import_data)
import io_tools.import_dro
reload(io_tools.import_dro)
import io_tools.propagate
reload(io_tools.propagate)
import dicom_tools.create_roicol
reload(dicom_tools.create_roicol)
import dro_tools.create_dro
reload(dro_tools.create_dro)


import time
import argparse
from io_tools.fetch_config import ConfigFetcher
from io_tools.download_data import DataDownloader
from io_tools.import_data import DataImporter
from io_tools.import_dro import DroImporter
from io_tools.propagate import Propagator
from dicom_tools.create_roicol import RoicolCreator
from dro_tools.create_dro import DroCreator


def main(
        cfgDir, runID, printSummary=False, plotResults=False):
    """
    Main script for fetching the config settings, downloading data from XNAT,
    importing of source ROI Collection and source and target DICOM series, 
    copy/propagation of source ROI Collection to the target dataset, creation
    of new ROI Collection, creation of a DRO and uploading to XNAT (if 
    applicable).
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    printSummary : bool, optional
        If True, summarising results will be printed. The default is False.
    plotResults : bool, optional
        If True, results will be printed. The default is False.
    
    Returns
    -------
    None.
    """
    
    #cfgDir = r'' + cfgDir # convert to an r path
    #print(f'\ncfgDir = {cfgDir}\n')
    
    # Store time stamps for various steps:
    times = [time.time()]
    
    # Instanstantiate a ConfigFetcher object, get the config settings, and 
    # get the alias token (which may or may not exist):
    cfgObj = ConfigFetcher(cfgDir, runID)
    cfgObj.get_config()
    cfgObj.update_cfgDict()
    #cfgObj.get_alias_token()
    
    # Instantiate a DataDownloader object, download the data and create
    # pathsDict:
    params = DataDownloader(cfgObj)
    params.download_and_get_pathsDict()
    
    # Instantiate a DataImporter object for source and import the data:
    srcDataset = DataImporter(params, 'src')
    srcDataset.import_data(params)
    
    # Instantiate a DataImporter object for target and import the data:
    trgDataset = DataImporter(params, 'trg')
    trgDataset.import_data(params)
    
    # Instantiate a DROImporter object and fetch the DRO (if applicable):
    droObj = DroImporter()
    droObj.fetch_dro(params)
    
    times.append(time.time())
    dTime = times[-1] - times[-2]
    
    # Instantiate a Propagator object and copy/propagate the source ROI
    # Collection to the target dataset:
    newDataset = Propagator(srcDataset, trgDataset, params)
    newDataset.execute(srcDataset, trgDataset, params, droObj.dro)
    
    times.append(time.time())
    dTime = times[-1] - times[-2]
    
    print('\n\n\n*** TIMINGS ***')
    [print(msg) for msg in params.timingMsgs if 'Took' in msg]
    print(f'Total time {dTime:.1f} s ({dTime/60:.1f} min) to',
          f'execute runID {runID}.\n\n\n')
    
    if printSummary:
        newDataset.print_summary_of_results(srcDataset, trgDataset)
    
    if plotResults:
        newDataset.plot_metric_v_iters(params)
        newDataset.plot_res_results(srcDataset, trgDataset, params)
        newDataset.plot_roi_over_dicom_im(
            srcDataset, trgDataset, params
            )
    
    # Instantiate RoicolCreator, create the new ROI Collection, check
    # for errors, export, upload to XNAT, and plot results 
    # (conditional):
    roicolObj = RoicolCreator()
    roicolObj.create_roicol(srcDataset, trgDataset, newDataset, params)
    roicolObj.error_check_roicol(srcDataset, trgDataset, newDataset, params)
    roicolObj.export_roicol(params)
    roicolObj.upload_roicol(params)
    if plotResults:
        roicolObj.plot_roi_over_dicoms(
            srcDataset, trgDataset, newDataset, params
            )
    
    # Instantiate a DroCreator object, create a new DRO, export it to
    # disk, and upload to XNAT:
    newDroObj = DroCreator(newDataset, params)
    newDroObj.create_dro(srcDataset, trgDataset, newDataset, params)
    newDroObj.export_dro(params)
    newDroObj.upload_dro(params)

if __name__ == '__main__':
    """
    Run app.py as a script.
    
    Example usage in a console:
    
    python app.py configs NCITA_TEST_RR2
    """
    
    parser = argparse.ArgumentParser(description='Arguments for main()')
    
    parser.add_argument(
        "cfgDir", 
        help="The directory containing config files"
        )
    
    parser.add_argument(
        "runID", 
        help="The run ID to execute"
        )
    
    parser.add_argument(
        "--printSummary", 
        action="store_true",
        #type=bool,
        help="Print summary if True"
        )
    
    parser.add_argument(
        "--plotResults", 
        action="store_true",
        #type=bool,
        help="Plot results if True"
        )
    
    args = parser.parse_args()
    
    # Run main():
    main(args.cfgDir, args.runID, args.printSummary, args.plotResults)