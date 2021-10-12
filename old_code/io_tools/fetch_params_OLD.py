# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:12:11 2021

@author: ctorti
"""

import os
from io_tools.imports import import_dict_from_json

class Params:
    """
    This class fetches parameters from configuration files.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    dict
        A dictionary of XNAT config, or general or Source XNAT or Target XNAT
        parameters.
    """
    
    def __init__(self, cfgDir):
        self.cfgDir = cfgDir
        self.get_gen_params()
        self.get_xnat_cfg()
        self.get_xnat_params(srcORtrg='src')
        self.get_xnat_params(srcORtrg='trg')
    
    def get_gen_params(self):
        """
        Fetches general parameters from general_params.json.
        
        The parameters are returned as a dictionary.
        
        Returns
        -------
        dict
            A dictionary containing general parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')
        
        fpath = os.path.join(cfgDir, 'general_params.json')
        
        self.genParams = import_dict_from_json(filepath=fpath)
        #return import_dict_from_json(filepath=fpath)
    
    def get_xnat_cfg(self):
        """
        Fetches XNAT configuration parameters from xnat_config.json.
        
        The parameters are returned as a dictionary.
        
        Returns
        -------
        dict
            A dictionary containing XNAT configuration parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, 'xnat_config.json')
        
        self.xnatCfg = import_dict_from_json(filepath=fpath)
        #return import_dict_from_json(filepath=fpath)
        
    def get_xnat_params(self, srcORtrg):
        """
        Fetches XNAT parameters for Source (from src_xnat_params.json) or
        for Target (from trg_xnat_params.json).
        
        The parameters are returned as a dictionary.
        
        Parameters
        ----------
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        dict
            A dictionary containing XNAT parameters for Source or Target.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, f'{srcORtrg}_xnat_params.json')
        
        if srcORtrg == 'src':
            self.srcXnatParams = import_dict_from_json(filepath=fpath)
        else:
            self.trgXnatParams = import_dict_from_json(filepath=fpath)
        #return import_dict_from_json(filepath=fpath)