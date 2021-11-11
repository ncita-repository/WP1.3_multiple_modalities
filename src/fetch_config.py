# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 10:31:20 2021

@author: ctorti
"""

from importlib import reload
import io_tools.general
reload(io_tools.general)

import os
import argparse
from io_tools.general import get_fpaths, get_user_input_as_int
from io_tools.imports import import_dict_from_json
from io_tools.exports import export_dict_to_json
from general_tools.geometry import (
    prop_of_segs_in_extent, prop_of_rois_in_extent
    )
from general_tools.general import are_items_equal_to_within_eps

class ConfigFetcher:
    # TODO! Update docstrings
    """
    This class fetches parameters from a local file and returns an
    object containing various parameters for a specified runID (stored in the 
    dictionary self.cfgDict).
    
    The parameters are imported from a JSON file for the desired run.
    
    Parameters
    ----------
    #cfgDir : str
    #    The path to the directory containing config files.
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    
    Returns
    -------
    #self.cfgDir : str
    #    The path to the directory containing config files.
    self.runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    self.globalVars : dict
        Dictionary containing the global variables.
    self.cfgDict : dict
        Dictionary containing the parameters for the desired run.
    """
    
    def __init__(self, runID):
        self.runID = runID
        #self.cfgDir = cfgDir
        #self.cfgDict = {} # will update later
        
        self.get_global_vars()
        #print(f"Global variables:\n\n {self.globalVars}\n\n")
        # Get XNAT config for chosen runID:
        self.get_xnat_config()
        #print(f"XNAT config:\n\n {self.cfgDict}\n\n")
        self.add_global_vars()
        #print(f"XNAT config with added global vars:\n\n {self.cfgDict}\n\n")
        self.update_paths()
        #print(f"XNAT config with updated paths:\n\n {self.cfgDict}\n\n")
        self.export_cfgDict()
        
    def get_global_vars(self):
        """
        Fetches the global variables from src/global_variables.json.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        self.globalVars : dict
            Dictionary containing the global variables.
        self.xnatCfgDir : str
            The XNAT configurations directory.
        """
        
        fname = 'global_variables.json'
        
        print(f'Fetching global variables from {fname}\n')
        
        try:
            self.globalVars = import_dict_from_json(fname)
            self.xnatCfgDir = self.globalVars['xnatCfgDir']
        except FileNotFoundError:
            print(f'File {fname} not found')
        
    def get_xnat_config(self):
        """
        Fetches the XNAT configuration parameters from a JSON in 
        src/xnat_configs with filename that matches self.runID.
        
        Parameters
        ----------
        self.runID : str
            The ID that determines the main configuration parameters to use for
            the run, and should match with a file name of a JSON.
        
        Returns
        -------
        self.cfgDict : dict
            Dictionary containing the XNAT parameters for the desired run.
        """
        
        #cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')
        runID = self.runID
        xnatCfgDir = self.xnatCfgDir

        cfgFpath = os.path.join(xnatCfgDir, f'{runID}.json')
        
        print(f"Fetching parameters for runID '{runID}' from {cfgFpath}\n")
        
        try:
            cfgDict = import_dict_from_json(cfgFpath)
        
        except FileNotFoundError:
            # Get a list of JSONs in xnatCfgDir:
            fpaths = get_fpaths(dpath=xnatCfgDir, ext='json')
            
            F = len(fpaths)
            
            if F == 0:
                msg = f"There were no JSON files found in {xnatCfgDir}.\n"\
                      + "Please check the input 'cfgDir' and try again."
            else:
                msg = f"There were no JSON files found in {xnatCfgDir} with "\
                    + f"file name matching '{runID}'.\nPlease select from one"\
                    + f" of the {F} files below:"
                
                for i in range(len(fpaths)):
                    msg += f'\n{i + 1}.  {fpaths[i]}'
                
                msg += "Enter an integer"
            
            #raise Exception(msg)
            choice = get_user_input_as_int(message=msg, minVal=1, maxVal=F)
            
            fpath = fpaths[choice - 1] # choice is 1-indexed
            
            cfgDict = import_dict_from_json(filepath=fpath)
            
        self.cfgDict = cfgDict
    
    def add_global_vars(self):
        """
        Add global variables to cfgDict that are not already defined in the
        dictionary.
        
        Any variables already defined in cfgDict will take priority over the 
        values in globalVars.
        
        Parameters
        ----------
        self.globalVars : dict
            Dictionary containing the global variables.
        self.cfgDict : dict
            Dictionary containing the parameters for the desired run.
        
        Returns
        -------
        self.cfgDict : dict
            Updated dictionary containing additional parameters from 
            self.globalVars.
        """
        
        globalVars = self.globalVars
        cfgDict = self.cfgDict
        
        for key, value in globalVars.items():
            if not key in cfgDict.keys():
                cfgDict[key] = value
        
        print('Global variables added to cfgDict\n')
        
        self.cfgDict = cfgDict
        
    def update_paths(self):
        """
        Update the current working directory (and subsequent child directories)
        in the configuration dictionary to reflect the 'workdir' set as an 
        environmental variable (if applicable).
        
        Parameters
        ----------
        self.cfgDict : dict
            Dictionary containing the parameters for the desired run.
        
        Returns
        -------
        self.cfgDict : dict
            Updated dictionary.
        """
        
        cfgDict = self.cfgDict
        
        # Try to get the current working directory from an environmnt variable:
        cwd = os.getenv('workdir')
        #print(f'os.getenv (in fetch_config.py) = {cwd}')
        if cwd == None:
            cwd = os.getcwd()
        #print(f'cwd (in fetch_config.py) = {cwd}')
        
        # Lines below added 03/11/21
        inputsDir = os.path.join(cwd, r'inputs')
        outputsDir = os.path.join(cwd, r'outputs')
        sampleDroDir = os.path.join(inputsDir, r'sample_DROs')
        fidsDir = os.path.join(inputsDir, r'fiducials')
        rtsExportDir = os.path.join(outputsDir, r'new_RTS')
        segExportDir = os.path.join(outputsDir, r'new_SEG')
        droExportDir = os.path.join(outputsDir, r'new_DRO')
        rtsPlotsExportDir = os.path.join(outputsDir, 'plots_RTS')
        segPlotsExportDir = os.path.join(outputsDir, 'plots_SEG')
        resPlotsExportDir = os.path.join(outputsDir, r'plots_res')
        txExportDir = os.path.join(outputsDir, r'transforms')
        imExportDir = os.path.join(outputsDir, r'images')
        labimExportDir = os.path.join(outputsDir, r'label_images')
        logsExportDir = os.path.join(outputsDir, r'logs')
    
        
        cfgDict['cwd'] = cwd
        # Lines below added 03/11/21
        cfgDict['inputsDir'] = inputsDir
        cfgDict['outputsDir'] = outputsDir
        cfgDict['sampleDroDir'] = sampleDroDir
        cfgDict['fidsDir'] = fidsDir
        cfgDict['rtsExportDir'] = rtsExportDir
        cfgDict['segExportDir'] = segExportDir
        cfgDict['droExportDir'] = droExportDir
        cfgDict['rtsPlotsExportDir'] = rtsPlotsExportDir
        cfgDict['segPlotsExportDir'] = segPlotsExportDir
        cfgDict['resPlotsExportDir'] = resPlotsExportDir
        cfgDict['txExportDir'] = txExportDir
        cfgDict['imExportDir'] = imExportDir
        cfgDict['labimExportDir'] = labimExportDir
        cfgDict['logsExportDir'] = logsExportDir
        
        print('Paths updated in cfgDict\n')
        
        self.cfgDict = cfgDict
    
    def export_cfgDict(self):
        """
        Export the final configuration dictionary to disk.
        
        Parameters
        ----------
        self.cfgDict : dict
            Dictionary containing the parameters for the desired run.
        
        Returns
        -------
        None.
        """
        
        cfgDict = self.cfgDict
        cwd = cfgDict['cwd']
        
        fname = 'cfgDict.json'
        
        export_dict_to_json(
            dictionary=cfgDict, filename=fname, exportDir=cwd
            )
        
        print(f'cfgDict exported to {os.path.join(cwd, fname)}\n')
    
    def get_intersection_of_roi_and_trgIm(
            self, srcDataset, trgDataset, params
            ):
        """
        Get the intersection of the ROI and target image domain as a fraction.
        
        This method requires objects that won't exist at the time the class is
        instantiated, but has been placed in this class because it builds up
        parameters containted within cfgDict.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            Contains various data for the source DICOM series.
        trgDataset : DataImporter Object
            Contains various data for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.fracProp : float
            The fractional proportion (normalised to 1) of the points or voxels
            that define the contour(s)/segmentation(s) to be copied/propagated
            intersect the extent of trgIm.
        """
    
        cfgDict = params.cfgDict
        roicolMod = cfgDict['roicolMod']
        p2c = cfgDict['p2c']
        
        if roicolMod == 'SEG':
            fracProp = prop_of_segs_in_extent(
                pixarrBySeg=srcDataset.pixarrBySeg, 
                f2sIndsBySeg=srcDataset.f2sIndsBySeg, 
                srcIm=srcDataset.dcmIm, 
                trgIm=trgDataset.dcmIm, 
                p2c=p2c
                )
        else:
            fracProp = prop_of_rois_in_extent(
                ptsByCntByRoi=srcDataset.ptsByCntByRoi,
                trgIm=trgDataset.dcmIm, 
                p2c=p2c
                )
        
        if p2c:
            print(f'{100*fracProp:.2f}% of the Source ROI Collection '
                  'intersects the Target image domain.\n')
            
        self.fracProp = fracProp
    
    def which_use_case(self, srcDataset, trgDataset, params):
        # TODO update docstrings
        """
        Determine which Use Cases applies.
        
        This method requires objects that won't exist at the time the class is
        instantiated, but has been placed in this class because it builds up
        parameters containted within cfgDict.
        
        Parameters
        ----------
        srcDataset : DataImporter Object
            DataImporter Object for the source DICOM series.
        trgDataset : DataImporter Object
            DataImporter Object for the target DICOM series.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
            
        Returns
        -------
        params.cfgDict['useCaseThatApplies'] : str
            The use case that applies.
        params.cfgDict['useCaseToApply'] : str
            The use case that will be applied.
        
        Note
        ----
        This is similar to the function which_use_case in 
        io_tools.inputs_checker.py, except that the functional form takes a
        cfgDict (params.cfgDict) and imports various data. This method uses
        data already imported via the DataImporter Objects srcDataset and
        trgDataset.
        """
        
        trgSlcNum = params.cfgDict['trgSlcNum']
        forceReg = params.cfgDict['forceReg']
        p2c = params.cfgDict['p2c']
        
        srcStudyUID = srcDataset.studyuid
        srcSeriesUID = srcDataset.seriesuid
        srcSOPuids = srcDataset.sopuids
        srcFORuid = srcDataset.foruid
        srcIPPs = srcDataset.imPositions
        srcDirs = srcDataset.imDirections
        srcSpacings = srcDataset.imSpacings
        srcST = srcDataset.imSlcThick
        srcSize = srcDataset.imSize
        #srcIm = srcDataset.image
        srcImExtent = srcDataset.imExtent
        
        trgStudyUID = trgDataset.studyuid
        trgSeriesUID = trgDataset.seriesuid
        trgSOPuids = trgDataset.sopuids
        trgFORuid = trgDataset.foruid
        trgIPPs = trgDataset.imPositions
        trgDirs = trgDataset.imDirections
        trgSpacings = trgDataset.imSpacings
        trgST = trgDataset.imSlcThick
        trgSize = trgDataset.imSize
        #trgIm = trgDataset.image
        trgImExtent = trgDataset.imExtent
        
        if trgSlcNum:
            if trgSlcNum > len(trgSOPuids) - 1:
                msg = f'The input argument "trgSlcNum" = {trgSlcNum} exceeds '\
                    + f'the number of Target DICOMs ({len(trgSOPuids)}).'
                raise Exception(msg)
        
        # Are the Source and Target part of the same study, series and do they 
        # have the same SOP UIDs?
        if srcStudyUID == trgStudyUID:
            sameStudy = True
        else:
            sameStudy = False
            
        if srcSeriesUID == trgSeriesUID:
            sameSeries = True
        else:
            sameSeries = False
            
        if srcFORuid == trgFORuid:
            sameFOR = True
        else:
            sameFOR = False
            
        if srcSOPuids == trgSOPuids:
            sameSOPs = True
        else:
            sameSOPs = False
        
        # Are the Source and Target spacings, directions and origins equal to  
        # within 1e-06?
        sameSpacings = are_items_equal_to_within_eps(srcSpacings, trgSpacings)
        sameDirs = are_items_equal_to_within_eps(srcDirs, trgDirs)
        sameOrigins = are_items_equal_to_within_eps(srcIPPs[0], trgIPPs[0])
        #sameST = are_items_equal_to_within_eps([srcST], [trgST])
        sameST = are_items_equal_to_within_eps(srcST, trgST)
        
        if srcSize == trgSize:
            sameSize = True
        else:
            sameSize = False
        
        # Check which case follows:
        if sameFOR:
            if sameDirs:
                """ 
                See Comment 11/02 #1 in which_use_case() in 
                io_tools.inputs_checker.py
                """
                if sameSpacings and sameOrigins:
                    if sameSeries:
                        # Direct copy
                        useCaseThatApplies = '1'
                    else:
                        if isinstance(trgSlcNum, int):
                            useCaseThatApplies = '2a' # Direct copy
                        else:
                            useCaseThatApplies = '2b' # Relationship-preserving copy
                else:
                    """ 
                    See Comment 11/02 #2 in which_use_case() in 
                    io_tools.inputs_checker.py
                    """
                    if isinstance(trgSlcNum, int):
                        useCaseThatApplies = '3a' # Direct copy
                    else:
                        useCaseThatApplies = '3b' # Relationship-preserving copy
            else:
                if isinstance(trgSlcNum, int):
                    useCaseThatApplies = '4a' # Direct copy
                else:
                    useCaseThatApplies = '4b' # Relationship-preserving copy
        else:
            if isinstance(trgSlcNum, int):
                useCaseThatApplies = '5a' # Direct copy 
            else:
                useCaseThatApplies = '5b' # Relationship-preserving copy
        
        if p2c:
            # Print comparisons to the console:
            
            # Are the image extents equal to within 1e-06?
            sameExtents = are_items_equal_to_within_eps(
                srcImExtent, trgImExtent
                )
            
            print('\n\n', '-'*120)
            print('Results of running which_use_case():')
            print('\n\n', '-'*120)
            print('Source and Target are in the/have:')
            
            if sameStudy:
                print('   * same study')
            else:
                print('   * different studies')
            
            if sameFOR:
                print('   * same Frame of Reference')
            else:
                print('   * different Frames of Reference')
                
            if sameSeries:
                print('   * same series')
            else:
                print('   * different series')
                
            if sameSOPs:
                print('   * same DICOMs')
            else:
                print('   * different DICOMs')
            
            if sameSize:
                print(f'   * same image size\n      Source/Target: {srcSize}')
            else:
                print(f'   * different image sizes:\n      Source: {srcSize}',
                      f'\n      Target: {trgSize}')
            
            if sameSpacings:
                print(f'   * same voxel sizes\n      Source/Target: {srcSpacings}')
            else:
                print(f'   * different voxel sizes:\n      Source: {srcSpacings}',
                      f'\n      Target: {trgSpacings}')
            
            if sameST:
                print(f'   * same slice thickness\n      Source/Target: {srcST}')
            else:
                print(f'   * different slice thickness:\n      Source: {srcST}',
                      f'\n      Target: {trgST}')
            
            if sameExtents:
                print(f'   * same image extents\n      Source/Target: {srcImExtent}')
            else:
                print(f'   * different image extents:\n      Source: {srcImExtent}',
                      f'\n      Target: {trgImExtent}')
            
            if sameOrigins:
                print(f'   * same origin\n      Source/Target: {srcIPPs[0]}')
            else:
                print(f'   * different origins:\n      Source: {srcIPPs[0]}',
                      f'\n      Target: {trgIPPs[0]}')
                
            if sameDirs:
                print('   * same patient orientation\n',
                      f'      Source/Target: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                      f'                      {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                      f'                      {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]')
            else:
                print('   * different patient orientations:\n',
                      f'      Source: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                      f'               {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                      f'               {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]\n',
                      f'      Target: [{trgDirs[0]}, {trgDirs[1]}, {trgDirs[2]},\n',
                      f'               {trgDirs[3]}, {trgDirs[4]}, {trgDirs[5]},\n',
                      f'               {trgDirs[6]}, {trgDirs[7]}, {trgDirs[8]}]')
        
        if True:#p2c:
            if forceReg and useCaseThatApplies in ['3a', '3b', '4a', '4b']:
                useCaseToApply = useCaseThatApplies.replace('3', '5').replace('4', '5')
                
                msg = f'\n*UseCase {useCaseThatApplies} applies but since '\
                      + f'forceReg = {forceReg}, UseCase {useCaseToApply} '\
                      + 'will be applied.'
                print(msg)
            else:
                useCaseToApply = useCaseThatApplies
                
                print(f'\n*UseCase {useCaseThatApplies} applies.')
                
            print('-'*120)
        
        params.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        params.cfgDict['useCaseToApply'] = useCaseToApply


if __name__ == '__main__':
    """
    Run fetch_config.py as a script.
    
    Example usage in a console:
    
    cd C:\Code\WP1.3_multiple_modalities\src
    python fetch_config.py runID
    """
    
    
    parser = argparse.ArgumentParser(
        description='Arguments for fetch_config()'
        )
    
    parser.add_argument(
        "runID", 
        help="The run ID to execute (name of file in src/xnat_configs/)"
        )
    
    args = parser.parse_args()
    
    # Instantiate ConfigFetcher object:
    cfgObj = ConfigFetcher(args.runID)