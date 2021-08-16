# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 10:31:20 2021

@author: ctorti
"""

from importlib import reload
import io_tools.general
reload(io_tools.general)

import os
from getpass import getpass
from io_tools.general import get_fpaths, get_user_input_as_int
from io_tools.imports import import_dict_from_json


class ConfigFetcher:
    """
    This class fetches parameters from a local file and downloads data from
    XNAT.
    
    The parameters are imported from a JSON file for the desired run.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    xnatSession : Requests Object, optional
        A Requests Object for an existing XNAT session. The default is None.
    
    Returns
    -------
    Config Object
        An Object containing various parameters (as a dictionary) for the 
        specified runID.
    """
    
    def __init__(self, cfgDir, runID, xnatSession=None):
        self.cfgDir = cfgDir
        
        self.get_config(runID)
    
    def get_config(self, runID):
        """
        Fetches the configuration parameters from an appropriate JSON.
        
        The parameters are in a dictionary.
        
        Parameters
        ----------
        runID : str
            The ID that determines the main configuration parameters to use for
            the run, and should match with a file name of a JSON.
        
        Returns
        -------
        cfgDict : dict
            Dictionary containing the parameters for the desired run.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, f'{runID}.json')
        
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
        
        # url and/or username and/or password may be blank (e.g. for security
        # reasons). If so prompt user for input and update cfgDict:
        if cfgDict['url'] == '':
            url = getpass("Enter XNAT url: ")
            cfgDict['url'] = url
        
        if cfgDict['username'] == '':
            username = getpass("Enter user name: ")
            cfgDict['username'] = username
        
        if cfgDict['password'] == '':
            password = getpass("Enter password: ")
            cfgDict['password'] = password
        
        self.cfgDict = cfgDict