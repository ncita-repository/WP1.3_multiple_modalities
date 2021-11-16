# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:10:00 2021

@author: ctorti
"""


"""
The purpose of this module is to pull a JSON configuation file from src/configs 
and export the dictionary to src/. The parameters stored in the file (with 
default filename xnatCfg) are those that will be required from the XNAT 
container service. 

Since the code has been developed to work as a stand-alone tool, this module 
creates the file by copying a pre-populated config file from src/configs
that are specific to the run to be executed.

Hence the only function of this module is to provide a temporary means of 
creating a config file that will ultimately be provided by the Container
Service.

The optional argument xnatCfgFname allows for the possibility of scaling
up for multiple runs concurrently. Each run can generate a new xnatCfg file
with unique file name, e.g. 'xnatCfg_34j2cf'.
"""


import os
import argparse
from io_tools.general import get_fpaths, get_user_input_as_int
from io_tools.imports import import_dict_from_json
from io_tools.exports import export_dict_to_json
from io_tools.fetch_configs import get_global_vars


# Import the global variables:
globalVars = get_global_vars()

#print(f"globalVars['xnatCfgDir'] = {globalVars['xnatCfgDir']}")
    
def get_xnat_config(runID, xnatCfgDir=globalVars['xnatCfgDir']):
    """
    Fetches the XNAT configuration parameters from a JSON with file name runID
    from a directory containing XNAT config files whose path is specified by 
    key 'xnatCfgDir' in globalVars.
    
    The filename of the XNAT config file (= runID) should match the key 'runID'
    in the config file.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for
        the run, and should match with a file name of a JSON.
    xnatCfgDir : str, optional
        The directory containing the XNAT configuration files. The default
        value is globalVars['xnatCfgDir'].
    
    Returns
    -------
    xnatCfg : dict
        Dictionary containing the XNAT parameters for the desired run.
    """
    
    xnatCfg = {}
    
    xnatCfgFpath = os.path.join(xnatCfgDir, f'{runID}.json')
    
    print(f"Fetching parameters for runID '{runID}' from {xnatCfgFpath}\n")
    
    try:
        xnatCfg = import_dict_from_json(xnatCfgFpath)
        
        runID_from_file = xnatCfg['runID']
        
        if not runID == runID_from_file:
            msg = f"\nWARNING: The XNAT config file name '{runID}' does not "+\
                "match the runID key value '{runID_from_file}'.\n"
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
        
        xnatCfg = import_dict_from_json(filepath=fpath)
        
    return xnatCfg

def export_xnatCfg(xnatCfg, xnatCfgFname='xnatCfg'):
    """
    Export the final configuration dictionary to src/.
    
    Parameters
    ----------
    xnatCfg : dict
        Dictionary containing the parameters for the desired run.
    xnatCfgFname : str, optional
        The filename to assign to the exported JSON file. The default value
        is 'xnatCfg'.
    
    Returns
    -------
    None.
    """
    
    cwd = globalVars['cwd']
    
    export_dict_to_json(
        dictionary=xnatCfg, filename=xnatCfgFname, exportDir=cwd
        )
    
    print(f'xnatCfg exported to {os.path.join(cwd, xnatCfgFname)}\n')

if __name__ == '__main__':
    """
    Run select_xnat_config.py as a script.
    
    Example usage in a console:
    
    python select_xnat_config.py runID
    
    or
    
    python create_cfgDict.py runID --xnatCfgFname=xnatCfg_34j2cf
    """
    
    parser = argparse.ArgumentParser(
        description='Arguments for select_xnat_config.py'
        )
    
    parser.add_argument(
        "runID", 
        help="The run ID to execute (name of file in src/xnat_configs/)"
        )
    
    parser.add_argument(
        "--xnatCfgFname",
        nargs='?', default='xnatCfg', const='xnatCfg',
        help="Optional file name of the XNAT config file (default is xnatCfg)"
        )
    
    args = parser.parse_args()
    
    xnatCfg = get_xnat_config(args.runID, xnatCfgDir=globalVars['xnatCfgDir'])
    
    export_xnatCfg(xnatCfg, xnatCfgFname=args.xnatCfgFname)