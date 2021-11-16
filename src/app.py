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
import io_tools.fetch_configs
reload(io_tools.fetch_configs)
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


#import time
import argparse
from io_tools.fetch_configs import ConfigFetcher
#from io_tools.imports import import_dict_from_json
from io_tools.download_data import DataDownloader
from io_tools.import_data import DataImporter
from io_tools.import_dro import DroImporter
from io_tools.propagate import Propagator
from dicom_tools.create_roicol import RoicolCreator
from dro_tools.create_dro import DroCreator


def main(
        xnatCfgFname='xnatCfg', printSummary=False, plotResults=False):
    """
    Main script for fetching the config settings, downloading data from XNAT,
    importing of source ROI Collection and source and target DICOM series, 
    copy/propagation of source ROI Collection to the target dataset, creation
    of new ROI Collection, creation of a DRO and uploading to XNAT (if 
    applicable).
    
    Parameters
    ----------
    xnatCfgFname : str, optional
        The file name of the XNAT config JSON file (in src/) containing the 
        parameters to be run. The default value is 'xnatCfg'.
    printSummary : bool, optional
        If True, summarising results will be printed. The default is False.
    plotResults : bool, optional
        If True, results will be printed. The default is False.
    
    Returns
    -------
    None.
    """
    
    #print(f'\ncfgDir = {cfgDir}\n')
    
    # Store time stamps for various steps:
    #times = [time.time()]
    
    # Instanstantiate a ConfigFetcher object, get the config settings and
    # export it to src/cfgDict.json:
    cfgObj = ConfigFetcher(xnatCfgFname)
    
    # Overwrite p2c parameter to True for verbose output to the console:
    #cfgObj.cfgDict['p2c'] = True
    
    # Instantiate a DataDownloader object, establish a connection to XNAT
    # (use or creating an XNAT Alias Token), download the data and create
    # pathsDict:
    params = DataDownloader(cfgObj)
    params.download_and_get_pathsDict()
    
    # Instantiate a DataImporter object for source and import the data:
    srcDataset = DataImporter(params, 'src')
    srcDataset.import_data(params)
    
    # Instantiate a DataImporter object for target and import the data:
    trgDataset = DataImporter(params, 'trg')
    trgDataset.import_data(params)
    
    # Check whether the RTS/SEG of interest intersects the target image's
    # extent:
    cfgObj.get_intersection_of_roi_and_trgIm(srcDataset, trgDataset, params)
    
    # Determine which use case applies (and will be applied):
    cfgObj.which_use_case(srcDataset, trgDataset, params)
    
    if params.cfgDict['p2c']:
        print(f"runID = {params.cfgDict['runID']}")
        print(f"useCaseThatApplies = {params.cfgDict['useCaseThatApplies']}")
        print(f"useCaseToApply = {params.cfgDict['useCaseToApply']}\n")
    
    # Instantiate a DROImporter object and fetch the DRO (if applicable):
    droObj = DroImporter(params)
    
    #times.append(time.time())
    #dTime = times[-1] - times[-2]
    
    # Instantiate a Propagator object and copy/propagate the source ROI
    # Collection to the target dataset:
    newDataset = Propagator(srcDataset, trgDataset, params)
    newDataset.execute(srcDataset, trgDataset, params, droObj.dro)
    
    #times.append(time.time())
    #dTime = times[-1] - times[-2]
    
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
    
    timingMsg = "Took total of [*] to execute the run.\n"
    params.add_timestamp(timingMsg)
    
    print('\n\nSUMMARY\n*******')
    [print(msg) for msg in params.timingMsgs if 'Took' in msg]
    
    #dTime = params.timings[-1] - params.timings[0]
    #timingMsg = f"Took {dTime:.1f} s ({dTime/60:.1f} min) to execute runID " +\
    #    f"{cfgObj.cfgDict['runID']}.\n"
    #print(timingMsg)

if __name__ == '__main__':
    """
    Run app.py as a script.
    
    Example usage in a console:
    
    python app.py
    
    or 
    
    python app.py --xnatCfgFname=xnatCfg_34j2cf
    """
    
    parser = argparse.ArgumentParser(description='Arguments for main()')
    
    parser.add_argument(
        "--xnatCfgFname",
        nargs='?', default='xnatCfg', const='xnatCfg',
        help="Optional file name of XNAT config file (default is xnatCfg)"
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
    
    #main(args.cfgDir, args.runID, args.printSummary, args.plotResults)
    main(args.xnatCfgFname, args.printSummary, args.plotResults)