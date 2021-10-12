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
import io_tools.general
reload(io_tools.general)
import io_tools.inputs_checker
reload(io_tools.inputs_checker)
import xnat_tools.dros
reload(xnat_tools.dros)

import os
from getpass import getpass
from io_tools.general import get_fpaths, get_user_input_as_int
from io_tools.imports import import_dict_from_json
from io_tools.inputs_checker import which_use_case
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
from xnat_tools.format_pathsDict import get_scan_asr_fname_and_id
#from xnat_tools.dros import search_dro


class Fetcher:
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
    Params Object
        A Params Object containing various parameters for the specified runID.
    """
    
    def __init__(self, cfgDir, runID, xnatSession=None):
        self.cfgDir = cfgDir
        
        # Get XNAT config, Source XNAT, Target XNAT, general and resampling/
        # registration/transformation parameters:
        #self.get_xnat_cfg()
        #self.get_xnat_params(srcORtrg='src')
        #self.get_xnat_params(srcORtrg='trg')
        #self.get_gen_params()
        #self.get_res_params()
        #self.get_main_params(runID)
        #self.xnatCfg = self.get_xnat_cfg()
        #self.srcXnatParams = self.get_xnat_params(srcORtrg='src')
        #self.trgXnatParams = self.get_xnat_params(srcORtrg='trg')
        #self.genParams = self.get_gen_params()
        #self.resParams = self.get_res_params()
        #
        #self.get_gen_params()
        #self.get_main_params(runID)
        #
        self.get_config(runID)
        
        if xnatSession == None:
            # Establish XNAT connection:
            #xnatSession = create_session(self.xnatCfg)
            #xnatSession = create_session(self.mainParams)
            xnatSession = create_session(xnatCfg=self.cfgDict)
            
        self.xnatSession = xnatSession
        
        # Download data and create pathsDict:
        self.get_pathsDict(xnatSession=xnatSession)
        #self.pathsDict = self.get_pathsDict(xnatSession)
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(config=self.cfgDict)
        
        #self.useCaseThatApplies = useCaseThatApplies
        #self.useCaseToApply = useCaseToApply
        
        # Update cfgDict:
        self.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        self.cfgDict['useCaseToApply'] = useCaseToApply
        #self.mainParams = mainParams
    
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
        None
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
    
    def get_pathsDict(self, xnatSession):
        """
        Downloads scans and ROI Collections and returns a dictionary containing
        filepaths.
        
        Parameters
        ----------
        xnatSession : Requests Object
            Requests session of XNAT connection.
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        cfgDict = self.cfgDict
        
        print('*** Fetching source DICOM scan from XNAT...\n')
        
        pathsDict = download_scan(config=cfgDict, 
                                  srcORtrg='src',
                                  xnatSession=xnatSession)
        
        print('*** Fetching target DICOM scan from XNAT...\n')
        
        pathsDict = download_scan(config=cfgDict,
                                  srcORtrg='trg',
                                  xnatSession=xnatSession,
                                  pathsDict=pathsDict)
        
        print('*** Fetching source ROI Collection from XNAT...\n')
        
        pathsDict = download_im_asr(config=cfgDict,
                                    srcORtrg='src',
                                    xnatSession=xnatSession,
                                    pathsDict=pathsDict)
        
        trgRoicolName = cfgDict['trgRoicolName']
        
        if trgRoicolName != None:
            print('*** Fetching target ROI Collection from XNAT...\n')
            
            pathsDict = download_im_asr(config=cfgDict,
                                        srcORtrg='trg',
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
        
        print(f'pathsDict = {pathsDict}\n')
        
        cfgDict = self.cfgDict
        #print(f'cfgDict = {cfgDict}\n')
        
        projID = cfgDict['projID']
        subjLab = cfgDict['subjLab']
        roicolMod = cfgDict['roicolMod']
        srcExpLab = cfgDict['srcExpLab']
        srcScanID = cfgDict['srcScanID']
        srcRoicolName = cfgDict['srcRoicolName']
        trgExpLab = cfgDict['trgExpLab']
        trgScanID = cfgDict['trgScanID']
        trgRoicolName = cfgDict['trgRoicolName']
        
        srcDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['scans'][srcScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        trgDcmDir = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][trgExpLab]['scans'][trgScanID]['resources']\
                ['DICOM']['files']['dicomDir']
        
        srcRoicolFname, srcAsrID\
            = get_scan_asr_fname_and_id(
                pathsDict, projID, subjLab, srcExpLab, roicolMod, srcRoicolName
                )
        
        srcRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['assessors'][srcAsrID]['resources']\
                [roicolMod]['files'][srcRoicolFname]['roicolFpath']
        
        if trgRoicolName == None:
            trgRoicolFpath = None
        else:
            trgRoicolFname, trgAsrID\
                = get_scan_asr_fname_and_id(
                    pathsDict, projID, subjLab, trgExpLab, roicolMod, 
                    trgRoicolName
                    )
            
            trgRoicolFpath = pathsDict['projects'][projID]['subjects'][subjLab]\
                ['experiments'][trgExpLab]['assessors'][trgAsrID]['resources']\
                    [roicolMod]['files'][trgRoicolFname]['roicolFpath']
        
        # Add filepaths and directories to cfgDict:
        cfgDict['srcAsrID'] = srcAsrID
        cfgDict['srcRoicolFname'] = srcRoicolFname
        cfgDict['srcRoicolFpath'] = srcRoicolFpath
        cfgDict['srcDicomDir'] = srcDcmDir
        
        cfgDict['trgDicomDir'] = trgDcmDir
        if trgRoicolName != None:
            cfgDict['trgAsrID'] = trgAsrID
            cfgDict['trgRoicolFname'] = trgRoicolFname
            cfgDict['trgRoicolFpath'] = trgRoicolFpath
        else:
            cfgDict['trgAsrID'] = None
            cfgDict['trgRoicolFname'] = None
            cfgDict['trgRoicolFpath'] = None
            
        self.cfgDict = cfgDict