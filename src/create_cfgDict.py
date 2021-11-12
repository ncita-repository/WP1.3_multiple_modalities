# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:04:47 2021

@author: ctorti
"""


import os
import argparse
from io_tools.general import get_fpaths, get_user_input_as_int
from io_tools.imports import import_dict_from_json
from io_tools.exports import export_dict_to_json

class ConfigCreator:
    # TODO! Update docstrings
    """
    This class creates cfgDict.json, containing parameters for a specified 
    run to execute.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    cfgFname : str, optional
        The file name to assign to the config file (in src/) containing the 
        parameters to be run. The default value is 'cfgDict'.
    
    Returns
    -------
    self.runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    self.globalVars : dict
        Dictionary containing the global variables.
    self.cfgDict : dict
        Dictionary containing the parameters for the desired run.
    
    Note
    ----
    These are the parameters that will be required from the XNAT container
    service, and they are a combination of global variables 
    (from src/global_variables.json) and a JSON config file from src/configs
    that are specific to the run to be executed.
    
    The optional argument cfgFname allows for the possibility of scaling
    up for multiple runs concurrently. Each run can generate a unique cfgDict
    with unique file name, e.g. cfgDict_34j2cf.json, that can be passed to
    ConfigFetcher.
    """
    
    def __init__(self, runID, cfgFname='cfgDict'):
        if not '.json' in cfgFname:
            cfgFname += '.json'
        
        self.runID = runID
        self.cfgFname = cfgFname
        
        # Get the global variables:
        self.get_global_vars()
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
        
        print(f'\nFetching global variables from {fname}\n')
        
        try:
            self.globalVars = import_dict_from_json(fname)
            self.xnatCfgDir = self.globalVars['xnatCfgDir']
            #print(f"Global variables:\n\n {self.globalVars}\n\n")
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
        
        fname = self.cfgFname
        
        export_dict_to_json(
            dictionary=cfgDict, filename=fname, exportDir=cwd
            )
        
        print(f'cfgDict exported to {os.path.join(cwd, fname)}\n')
    
if __name__ == '__main__':
    """
    Run create_cfgDict.py as a script.
    
    Example usage in a console:
    
    python create_cfgDict.py runID
    
    or
    
    python create_cfgDict.py runID --cfgFname=cfgDict_42md31
    """
    
    parser = argparse.ArgumentParser(
        description='Arguments for ConfigCreator()'
        )
    
    parser.add_argument(
        "runID", 
        help="The run ID to execute (name of file in src/xnat_configs/)"
        )
    
    parser.add_argument(
        "--cfgFname",
        nargs='?', default='cfgDict', const='cfgDict',
        help="Optional file name of cfgDict file (default is cfgDict)"
        )
    
    args = parser.parse_args()
    
    cfgObj = ConfigCreator(args.runID, args.cfgFname)