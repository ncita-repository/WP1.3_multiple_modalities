# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 10:31:20 2021

@author: ctorti
"""

from importlib import reload

import xnat_tools.scans
reload(xnat_tools.scans)
import xnat_tools.im_assessors
reload(xnat_tools.im_assessors)
import xnat_tools.format_pathsDict
reload(xnat_tools.format_pathsDict)
import io_tools.inputs_checker
reload(io_tools.inputs_checker)
import xnat_tools.dros
reload(xnat_tools.dros)

import os
from io_tools.imports import import_dict_from_json
from io_tools.inputs_checker import which_use_case
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
from xnat_tools.format_pathsDict import get_scan_asr_fname_and_id
#from xnat_tools.dros import search_dro


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
    
    def __init__(self, cfgDir, xnatSession=None):
        self.cfgDir = cfgDir
        
        # Get XNAT config, Source XNAT, Target XNAT, general and resampling/
        # registration/transformation parameters:
        self.get_xnat_cfg()
        self.get_xnat_params(srcORtrg='src')
        self.get_xnat_params(srcORtrg='trg')
        self.get_gen_params()
        self.get_res_params()
        #self.xnatCfg = self.get_xnat_cfg()
        #self.srcXnatParams = self.get_xnat_params(srcORtrg='src')
        #self.trgXnatParams = self.get_xnat_params(srcORtrg='trg')
        #self.genParams = self.get_gen_params()
        #self.resParams = self.get_res_params()
        
        if xnatSession == None:
            # Establish XNAT connection:
            xnatSession = create_session(self.xnatCfg)
            
        self.xnatSession = xnatSession
        
        # Download data and create pathsDict:
        self.get_pathsDict(xnatSession)
        #self.pathsDict = self.get_pathsDict(xnatSession)
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(self.srcXnatParams, self.trgXnatParams, 
                             self.resParams, self.genParams)
        
        self.useCaseThatApplies = useCaseThatApplies
        self.useCaseToApply = useCaseToApply
    
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
    
    def get_res_params(self):
        """
        Fetches resampling-/registration-/transforming-related parameters from 
        reg_res_tx_params.json.
        
        The parameters are returned as a dictionary.
        
        Returns
        -------
        dict
            A dictionary containing resampling, registration and transformation
            related parameters.
        """
        
        cfgDir = self.cfgDir
        #print(f'cfgDir = {cfgDir}\n')

        fpath = os.path.join(cfgDir, 'res_params.json')
        
        self.resParams = import_dict_from_json(filepath=fpath)
        
    def get_pathsDict(self, xnatSession):
        """
        Downloads scans and ROI Collections and returns a dictionary containing
        filepaths.
        
        Parameters
        ----------
        params : Params Object
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        xnatCfg = self.xnatCfg
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        genParams = self.genParams
        
        pathsDict = download_scan(xnatCfg=xnatCfg, 
                                  xnatParams=srcXnatParams, 
                                  genParams=genParams,
                                  xnatSession=xnatSession)
    
        pathsDict = download_scan(xnatCfg=xnatCfg, 
                                  xnatParams=trgXnatParams, 
                                  genParams=genParams, 
                                  xnatSession=xnatSession,
                                  pathsDict=pathsDict)
    
        pathsDict = download_im_asr(xnatCfg=xnatCfg, 
                                    xnatParams=srcXnatParams,
                                    genParams=genParams, 
                                    xnatSession=xnatSession,
                                    pathsDict=pathsDict)
        
        trgRoicolName = trgXnatParams['roicolName']
        
        if trgRoicolName != None:
            pathsDict = download_im_asr(xnatCfg=xnatCfg, 
                                        xnatParams=trgXnatParams,
                                        genParams=genParams, 
                                        xnatSession=xnatSession,
                                        pathsDict=pathsDict)
        
        #return pathsDict
        self.pathsDict = pathsDict
    
    def add_dirs_and_fpaths(self):
        """
        Add directories and filepaths to params.srcXnatParams and 
        params.trgXnatParams from pathsDict.

        Returns
        -------
        None.

        """
        pathsDict = self.pathsDict
        
        #print(f'pathsDict = {pathsDict}\n')
        
        srcXnatParams = self.srcXnatParams
        trgXnatParams = self.trgXnatParams
        
        #print(f'srcXnatParams = {srcXnatParams}\n')
        
        projID = srcXnatParams['projID']
        subjLab = srcXnatParams['subjLab']
        srcExpLab = srcXnatParams['expLab']
        srcScanID = srcXnatParams['scanID']
        srcRoicolMod = srcXnatParams['roicolMod']
        srcRoicolName = srcXnatParams['roicolName']
        trgExpLab = trgXnatParams['expLab']
        trgScanID = trgXnatParams['scanID']
        trgRoicolMod = trgXnatParams['roicolMod']
        trgRoicolName = trgXnatParams['roicolName']
        
        srcDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['scans'][srcScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        trgDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][trgExpLab]['scans'][trgScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        srcRoicolFname, srcAsrID\
            = get_scan_asr_fname_and_id(pathsDict, projID, subjLab, srcExpLab,
                                        srcRoicolMod, srcRoicolName)
        
        srcRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['assessors'][srcAsrID]['resources']\
                [srcRoicolMod]['files'][srcRoicolFname]['roicolFpath']
        
        if trgRoicolName == None:
            trgRoicolFpath = None
        else:
            trgRoicolFname, trgAsrID\
                = get_scan_asr_fname_and_id(pathsDict, projID, subjLab, 
                                            trgExpLab, trgRoicolMod, 
                                            trgRoicolName)
        
            trgRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
                ['experiments'][trgExpLab]['assessors'][trgAsrID]['resources']\
                    [trgRoicolMod]['files'][trgRoicolFname]['roicolFpath']
        
        # Add filepaths and directories to srcXnatParams and trgXnatParams:
        srcXnatParams['asrID'] = srcAsrID
        srcXnatParams['roicolFname'] = srcRoicolFname
        srcXnatParams['roicolFpath'] = srcRoicolFpath
        srcXnatParams['dicomDir'] = srcDcmDir
        
        self.srcXnatParams = srcXnatParams
        
        trgXnatParams['dicomDir'] = trgDcmDir
        if trgRoicolName != None:
            trgXnatParams['asrID'] = trgAsrID
            trgXnatParams['roicolFname'] = trgRoicolFname
            trgXnatParams['roicolFpath'] = trgRoicolFpath
        else:
            trgXnatParams['asrID'] = None
            trgXnatParams['roicolFname'] = None
            trgXnatParams['roicolFpath'] = None
            
        self.trgXnatParams = trgXnatParams
        
        