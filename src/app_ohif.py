# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:24:09 2021

@author: ctorti
"""



""" 
Batch script for running the ROI Collection propagation.

Note:
This is intended as a backend script working alongside frontend work in the
OHIF-Viewer.

As opposed to the stand-alone version (app.py) this will not fetch data from
XNAT but instead will work on the basis that the data has been stored in a
temporary folder on the localhost. Likewise it will not search for a suitable
DRO in the subject assessors.

Hence for development the ROI Collections and scans have been stored in the
inputs directory.

This implementation will only need to perform ROI propagations (i.e. those
involving relationship-preserving that require resampling or image registration
to account for different voxel resolutions and/or frames of reference), since
simple ("direct") copying of ROIs will be handled exclusively in the front-end.

Note that none of the code that handles the simple cases (use cases 1-2) has
been removed.
""" 


import sys

code_root = r'C:\Code\WP1.3_multiple_modalities\src'

# Add code_root to the system path so packages can be imported from it:
sys.path.append(code_root)


from importlib import reload
import io_tools.fetch_config_ohif
reload(io_tools.fetch_config_ohif)
#import io_tools.download_data
#reload(io_tools.download_data)
import io_tools.import_data
reload(io_tools.import_data)
import io_tools.import_dro_ohif
reload(io_tools.import_dro_ohif)
import io_tools.propagate
reload(io_tools.propagate)
import dicom_tools.create_roicol
reload(dicom_tools.create_roicol)
import dro_tools.create_dro
reload(dro_tools.create_dro)


import time
import argparse
from io_tools.fetch_config_ohif import ConfigFromOhif
#from io_tools.download_data import DataDownloader
from io_tools.import_data import DataImporter
from io_tools.import_dro_ohif import DroImporterOhif
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
    
    cfgDir = r'' + cfgDir # convert to an r path
    #print(f'\ncfgDir = {cfgDir}\n')
    
    # Store time stamps for various steps:
    times = [time.time()]
    
    # Import the dictionary stored in cfgDict_ohif.json, which is a manually
    # created dictionary that contains metadata, directory paths of src and trg
    # scans, and filepaths to src and trg ROI Collections, all of which will be
    # assumed to be made available from the OHIF-Viewer front-end (hence by-
    # passing the steps taking in app.py involving the ConfigFetcher and 
    # DataDownloader classes):
    cfgObj = ConfigFromOhif(cfgDir, runID)
    cfgObj.get_config()
    
    #print(f"cfgDict = {cfgObj.cfgDict}")
    #print(f"cfgDict['cwd'] = {cfgObj.cfgDict['cwd']}")
    
    # Instantiate a DataImporter object for source and import the data:
    """
    Since data doesn't need to be downloaded, pass cfgObj to DataImporter
    rather than a DataDownloader class object.
    """
    srcDataset = DataImporter(cfgObj, 'src')
    srcDataset.import_data(cfgObj)
    
    # Instantiate a DataImporter object for target and import the data:
    trgDataset = DataImporter(cfgObj, 'trg')
    trgDataset.import_data(cfgObj)
    
    # Check whether the RTS/SEG of interest intersects the target image's
    # extent:
    """ Note:
        Instead of providing a params object (as in app.py), provide cfgObj,
        since it is the cfgDict attribute that is used.
    """
    cfgObj.get_intersection_of_roi_and_trgIm(srcDataset, trgDataset, cfgObj)
    
    # Determine which use case applies (and will be applied):
    """ Note:
        Instead of providing a params object (as in app.py), provide cfgObj,
        since it is the cfgDict attribute that is used.
    """
    cfgObj.which_use_case(srcDataset, trgDataset, cfgObj)
    
    # Instantiate a DROImporter object and fetch the DRO (if applicable):
    droObj = DroImporterOhif()
    droObj.import_dro(cfgObj)
    
    times.append(time.time())
    dTime = times[-1] - times[-2]
    
    # Instantiate a Propagator object and copy/propagate the source ROI
    # Collection to the target dataset:
    newDataset = Propagator(srcDataset, trgDataset, cfgObj)
    newDataset.execute(srcDataset, trgDataset, cfgObj, droObj.dro)
    
    times.append(time.time())
    dTime = times[-1] - times[-2]
    
    print('\n\n\n*** TIMINGS ***')
    [print(msg) for msg in cfgObj.timingMsgs if 'Took' in msg]
    print(f'Total time {dTime:.1f} s ({dTime/60:.1f} min) to',
          f'execute runID {runID}.\n\n\n')
    
    if printSummary:
        newDataset.print_summary_of_results(srcDataset, trgDataset)
    
    if plotResults:
        newDataset.plot_metric_v_iters(cfgObj)
        newDataset.plot_res_results(srcDataset, trgDataset, cfgObj)
        newDataset.plot_roi_over_dicom_im(
            srcDataset, trgDataset, cfgObj
            )
    
    # Instantiate RoicolCreator, create the new ROI Collection, check
    # for errors, export, (DO NOT*) upload to XNAT, and plot results 
    # (conditional) (*for this version):
    roicolObj = RoicolCreator()
    roicolObj.create_roicol(srcDataset, trgDataset, newDataset, cfgObj)
    roicolObj.error_check_roicol(srcDataset, trgDataset, newDataset, cfgObj)
    roicolObj.export_roicol(cfgObj)
    #roicolObj.upload_roicol(cfgObj)
    if plotResults:
        roicolObj.plot_roi_over_dicoms(
            srcDataset, trgDataset, newDataset, cfgObj
            )
    
    # Instantiate a DroCreator object, create a new DRO, export it to
    # disk, and DO NOT* upload to XNAT (*for this version):
    newDroObj = DroCreator(newDataset, cfgObj)
    newDroObj.create_dro(srcDataset, trgDataset, newDataset, cfgObj)
    newDroObj.export_dro(cfgObj)
    #newDroObj.upload_dro(cfgObj)

if __name__ == '__main__':
    """
    Run app.py as a script.
    
    Example usage in a console:
    
    python app.py C:\Code\WP1.3_multiple_modalities\src\configs NCITA_TEST_RR2
    """
    
    parser = argparse.ArgumentParser(description='Arguments for main()')
    
    msg = "The directory containing cfgDict_ohif.json (should be the src " +\
        "directory)."
    
    parser.add_argument(
        "cfgDir", 
        help=msg
        )
    
    msg = "The run ID to execute. Acceptable runIDs include:" +\
        "\n- 'RR3_contour'" +\
        "\n- 'RR3_roi'" +\
        "\n- 'RR3_roicol'" +\
        "\n- 'SR3_contour'" +\
        "\n- 'SR3_roi'" +\
        "\n- 'SR3_roicol'" +\
        "\n- 'RR11_contour'" +\
        "\n- 'RR11_roi'" +\
        "\n- 'RR11_roicol'" +\
        "\n- 'SR11_contour'" +\
        "\n- 'SR11_roi'" +\
        "\n- 'SR11_roicol'"
    
    parser.add_argument(
        "runID", 
        help=msg
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