# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 10:31:20 2021

@author: ctorti
"""

from importlib import reload
import io_tools.general
reload(io_tools.general)

import os
#from getpass import getpass
from io_tools.general import get_fpaths, get_user_input_as_int
from io_tools.imports import import_dict_from_json


class ConfigFetcher:
    # TODO! Update docstrings
    """
    This class fetches parameters from a local file and returns an
    object containing various parameters for a specified runID (stored in the 
    dictionary self.cfgDict), and an XNAT alias token (stored in the 
    dictionary self.aliasToken).
    
    The parameters are imported from a JSON file for the desired run.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    
    Returns
    -------
    self.cfgDir : str
        The path to the directory containing config files.
    self.runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    self.cfgDict : dict
        Dictionary containing the parameters for the desired run.
    #self.aliasToken : dict
    #    Dictionary containing an alias token, which may be valid, invalid,
    #    or empty (if no alias_token.json file was found).
    """
    
    def __init__(self, cfgDir, runID):
        self.cfgDir = cfgDir
        self.runID = runID
        self.cfgDict = {} # will update later
    
    def get_config(self):
        """
        Fetches the configuration parameters from an appropriate JSON.
        
        The parameters are in a dictionary.
        
        Parameters
        ----------
        self.cfgDir : str
            The path to the directory containing config files.
        self.runID : str
            The ID that determines the main configuration parameters to use for
            the run, and should match with a file name of a JSON.
        
        Returns
        -------
        cfgDict : dict
            Dictionary containing the parameters for the desired run.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')
        runID = self.runID

        fpath = os.path.join(cfgDir, f'{runID}.json')
        
        print('Fetching parameters from local drive...\n')
        
        try:
            cfgDict = import_dict_from_json(filepath=fpath)
        
        except FileNotFoundError:
            # Get a list of JSONs in cfgDir:
            fpaths = get_fpaths(dpath=cfgDir, ext='json')
            
            F = len(fpaths)
            
            if F == 0:
                msg = f"There were no JSON files found in {cfgDir}.\nPlease "\
                      + "check the input 'cfgDir' and try again."
            else:
                msg = f"There were no JSON files found in {cfgDir} whose file"\
                    + f" name matches '{runID}'.\nPlease select from one of "\
                    + f"the {F} files below or try a different 'cfgDir':"
                
                for i in range(len(fpaths)):
                    msg += f'\n{i + 1}.  {fpaths[i]}'
                
                msg += "Enter an integer"
            
            #raise Exception(msg)
            choice = get_user_input_as_int(message=msg, minVal=1, maxVal=F)
            
            fpath = fpaths[choice - 1] # choice is 1-indexed
            
            cfgDict = import_dict_from_json(filepath=fpath)
            
        self.cfgDict = cfgDict
    
    def update_cfgDict(self):
        """
        Update the current working directory (and subsequent child directories)
        in the configuration dictionary to reflect the 'workdir' set as an 
        environmental variable (if applicable).
        
        Parameters
        ----------
        cfgDict : dict
            Dictionary containing the parameters for the desired run.
        
        Returns
        -------
        cfgDict : dict
            Updated dictionary.
        """
        
        cfgDict = self.cfgDict
        
        # Try to get the current working directory from an environmnt variable:
        cwd = os.getenv('workdir')
        print(f'os.getenv (in fetch_config.py) = {cwd}')
        if cwd == None:
            cwd = os.getcwd()
        print(f'cwd (in fetch_config.py) = {cwd}')
        
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
        
        self.cfgDict = cfgDict