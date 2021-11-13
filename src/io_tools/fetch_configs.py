# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:04:23 2021

@author: ctorti
"""

import os
from io_tools.imports import import_dict_from_json
#from select_xnat_config import get_global_vars, get_xnat_config
from general_tools.geometry import (
    prop_of_segs_in_extent, prop_of_rois_in_extent
    )
from general_tools.general import are_items_equal_to_within_eps


def get_global_vars():
    """
    Fetches the global variables from global_variables.json in the current
    working directory.
    
    Parameters
    ----------
    None.
    
    Returns
    -------
    globalVars : dict
        Dictionary containing the global variables.
    """
    
    globalVars = {}
    
    fname = 'global_variables.json'
    
    print(f'\nFetching global variables from {fname}\n')
    
    try:
        globalVars = import_dict_from_json(fname)
        #print(f"Global variables:\n\n {globalVars}\n\n")
    except FileNotFoundError:
        print(f'File {fname} not found')
    
    return globalVars


class ConfigFetcher:
    # TODO! Update docstrings
    """
    This class imports two dictionaries from JSON files: global_variables.json
    (containing global variables) and xnatCfg.json (containing XNAT parameters
    specific to the run) from src/, and combines them to a new dictionary
    (cfgDict).
    
    The class also contains methods that add parameters to cfgDict once other 
    data are available.
    
    Parameters
    ----------
    xnatCfgFname : str, optional
        The file name of the XNAT config file (in src/configs) containing the 
        XNAT-related parameters for the run. The default value is 'xnatCfg'.
    
    Returns
    -------
    self.globalVars : dict
        Dictionary containing the global variables.
    self.xnatCfgFname : str
        Same as input argument xnatCfgFname.
    self.xnatCfg : dict
        Dictionary containing the XNAT config parameters for the desired run.
    self.cfgDict : dict
        Dictionary containing the XNAT config parameters for the desired run
        and the global variables.
    self.runID : str
        The ID for the run.
    self.fracProp : float
        The fractional proportion (normalised to 1) of the points or voxels
        that define the contour(s)/segmentation(s) to be copied/propagated
        intersect the extent of trgIm.
    
    Note
    ----
    The optional argument xnatCfgFname allows for the possibility of scaling
    up for multiple runs concurrently. Each run can import a unique XNAT config
    file (default name 'xnatCfg'), e.g. 'xnatCfg_34j2cf'.
    """
    
    def __init__(self, xnatCfgFname='xnatCfg'):
        self.globalVars = get_global_vars()
        
        self.xnatCfgFname = xnatCfgFname
        self.import_xnatCfg()
        
        self.add_global_vars()
        self.update_paths()
        
        # runID is already in xnatCfg but copy value to self.runID for 
        # convenience:
        self.runID = self.xnatCfg['runID']
        
    def import_xnatCfg(self):
        """
        Imports the JSON file containing XNAT config parameters, stored in 
        src/ with filename self.xnatCfgFname.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        self.xnatCfg : dict
            Dictionary containing the XNAT config parameters for the desired 
            run.
        """
        
        cwd = self.globalVars['cwd']
        fname = self.xnatCfgFname + '.json'
        
        print(f'Fetching XNAT config from {os.path.join(cwd, fname)}\n')
        
        try:
            self.xnatCfg = import_dict_from_json(fname)
        except FileNotFoundError:
            print(f'File {fname} not found')
    
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
        self.xnatCfg : dict
            Dictionary containing the XNAT config parameters for the desired 
            run.
        
        Returns
        -------
        self.cfgDict : dict
            Dictionary containing the items in self.xnatCfg and all items in
            self.globalVars not already present in self.xnatCfg.
        """
        
        globalVars = self.globalVars
        xnatCfg = self.xnatCfg
        
        # Initialise cfgDict by copying xnatCfg:
        cfgDict = dict(xnatCfg)
        
        # Add the key-values in globalVars for keys not already in cfgDict:
        for key, value in globalVars.items():
            if not key in cfgDict.keys():
                cfgDict[key] = value
        
        print('Global variables added to XNAT config\n')
        
        self.cfgDict = cfgDict
    
    def update_paths(self):
        """
        Update the current working directory (and subsequent child directories)
        in the configuration dictionary to reflect the 'workdir' set as an 
        environmental variable (if applicable).
        
        Parameters
        ----------
        self.cfgDict : dict
            Dictionary containing the XNAT config parameters for the desired 
            run and the global variables.
        
        Returns
        -------
        self.cfgDict : dict
            As above but with updated paths.
        """
        
        cfgDict = self.cfgDict
        globalVars = self.globalVars
        
        # Try to get the current working directory from an environmnt variable:
        cwd = os.getenv('workdir')
        #print(f'os.getenv (in fetch_config.py) = {cwd}')
        if cwd == None:
            cwd = os.getcwd()
        #print(f'cwd (in fetch_config.py) = {cwd}')
        
        cfgDict['cwd'] = cwd
        
        print(f"cfgDict['cwd'] = {cfgDict['cwd']}")
        
        # Change all directory paths (other than cwd) from relative (to cwd)
        # to absolute paths:
        for key in cfgDict.keys():
            if 'Dir' in key:
                cfgDict[key] = os.path.join(cwd, globalVars[key])
        
        print('Paths updated in cfgDict\n')
        
        self.cfgDict = cfgDict
    
    def get_intersection_of_roi_and_trgIm(
            self, srcDataset, trgDataset, params
            ):
        """
        Get the intersection of the entity (contour/segmentation/ROI/segment/
        ROI Collection) to be copied and target image domain as a fraction.
        
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
        
        timingMsg = "* Calculating the intersection between the entity to " +\
            "be copied and the target image extent...\n"
        params.add_timestamp(timingMsg)
        
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
        
        timingMsg = "Took [*] to calculate the intersection between the " +\
            "entity to be copied and the target image extent.\n"
        params.add_timestamp(timingMsg)
        
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
        self.cfgDict : dict
            Dictionary containing the parameters for the desired run, updated
            with additional keys 'useCaseThatApplies' and 'useCaseToApply'.
        
        Note
        ----
        This is similar to the function which_use_case in 
        io_tools.inputs_checker.py, except that the functional form takes a
        cfgDict (params.cfgDict) and imports various data. This method uses
        data already imported via the DataImporter Objects srcDataset and
        trgDataset.
        """
        
        timingMsg = "* Determining which use case applies and is to be " +\
            "applied...\n"
        params.add_timestamp(timingMsg)
        
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
        
        timingMsg = "Took [*] to determine which use case applies and is " +\
            "to be applied.\n"
        params.add_timestamp(timingMsg)
        
        params.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        params.cfgDict['useCaseToApply'] = useCaseToApply