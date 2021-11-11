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
import xnat_tools.alias_tokens
reload(xnat_tools.alias_tokens)
#import io_tools.inputs_checker
#reload(io_tools.inputs_checker)
#import xnat_tools.dros
#reload(xnat_tools.dros)
#import io_tools.fetch_config
#reload(io_tools.fetch_config)

import os
import time
#import requests
#from io_tools.inputs_checker import which_use_case
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
from xnat_tools.format_pathsDict import get_scan_asr_fname_and_id
from xnat_tools.alias_tokens import (
    import_alias_token, is_alias_token_valid, generate_alias_token,
    export_alias_token
    )
#from xnat_tools.dros import search_dro


class DataDownloader:
    # TODO! Update docstrings
    """
    This class downloads data from XNAT and creates stores the file paths in a 
    dictionary.
    
    Parameters
    ----------
    #cfgObj : ConfigFetcher Object
    #    An object containing various parameters for a specified runID (stored 
    #    in the dictionary self.cfgDict), and an XNAT alias token (stored in the 
    #    dictionary self.aliasToken).
    self.cfgDict : dict
        Dictionary containing the parameters for the desired run.
    xnatSession : Requests Object, optional
        A Requests Object for an existing XNAT session. The default is None.
    
    Returns
    -------
    cfgDict : dict
        Updated dictionary of configuration parameters.
    pathsDict : dict
        A dictionary containing file paths of the downloaded data.
    
    Note
    ----
    The DataDownloader object used in subsequent classes has the object 
    variable name 'params'.
    """
    
    #def __init__(self, cfgObj, xnatSession=None):
    #def __init__(self, cfgObj):
    def __init__(self, cfgDict):
        
        #self.cfgDir = cfgObj.cfgDir
        #self.cfgDict = cfgObj.cfgDict
        self.cfgDict = cfgDict
        #self.aliasToken = cfgObj.aliasToken
        self.aliasToken = {} # initial value
        
        """
        if xnatSession == None:
            # Establish XNAT connection:
            xnatSession = create_session(xnatCfg=cfgDict)
            
        self.xnatSession = xnatSession
        """
        
        self.establish_xnat_connection()
        
        # Initialise list of timestamps and timing messages to be stored:
        self.timings = [time.time()]
        self.timingMsgs = []
        
    def establish_xnat_connection(self):
        """
        Establish XNAT connection.
        
        Parameters
        ----------
        self.cfgDir : str
            The path to the directory containing config files.
        self.cfgDict : dict
            Dictionary containing the parameters for the desired run.
        self.aliasToken : dict
            Dictionary containing an XNAT alias token. The dictionary may be
            empty, and if not empty, the token may be valid or invalid.
        
        Returns
        -------
        self.xnatSession : requests.models.Response
            A requests session containing the url (xnatSession.url) and 
            authentication details (xnatSession.auth = ('username', 'password')).
        """
        
        #cfgDir = self.cfgDir
        cwd = self.cfgDict['cwd']
        url = self.cfgDict['url']
        username = self.cfgDict['username']
        aliasToken = self.aliasToken
        
        # Assumed directory that may contain an XNAT alias token:
        tokenDir = os.path.join(cwd, 'tokens')
        
        if aliasToken == {}:
            #print(f'Searching for an XNAT alias token in {cfgDir}')
            print(f'Searching for an XNAT alias token in {tokenDir}')
            
            # Import an XNAT alias token:
            #aliasToken = import_alias_token(cfgDir)
            aliasToken = import_alias_token(tokenDir)
        
        # Is the alias token valid?
        if is_alias_token_valid(aliasToken, url):
            print('Using XNAT alias token to establish XNAT connection.')
            xnatSession = create_session(url=url, aliasToken=aliasToken)
        else:
            print('Establishing XNAT connection without XNAT alias token.')
            xnatSession = create_session(url=url, username=username)
        
        # Generate an XNAT alias token to avoid making a new authenticated
        # user session:
        aliasToken = generate_alias_token(xnatSession)
        
        #export_alias_token(aliasToken, cfgDir)
        export_alias_token(aliasToken, tokenDir)
        
        self.xnatSession = xnatSession
        self.aliasToken = aliasToken
        
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
    
    def add_timestamp(self, timingMsg):
        # TODO update docstrings
        """
        Add timestamp and timing message. Replace the "keyword" '[*]' in
        timingMsg with the time difference dTime calculated below. If the
        keywork doesn't exist it's not a time-related message but just info.
        """
        self.timings.append(time.time())
        
        dTime = self.timings[-1] - self.timings[-2]
        
        if '[*]' in timingMsg:
            timingMsg = timingMsg.replace('[*]', f'{dTime:.2f}')
        self.timingMsgs.append(timingMsg)
        print(f'*{timingMsg}')
    
    def print_cfg_params_to_console(self):
        """ Print some useful info to the console. """
        
        cfgDict = self.cfgDict
        runID = cfgDict['runID']
        projID = cfgDict['projID']
        subjLab = cfgDict['subjLab']
        srcExpLab = cfgDict['srcExpLab']
        srcScanID = cfgDict['srcScanID']
        srcSlcNum = cfgDict['srcSlcNum']
        srcRoicolName = cfgDict['srcRoicolName']
        srcRoiName = cfgDict['srcRoiName']
        roicolMod = cfgDict['roicolMod']
        trgExpLab = cfgDict['trgExpLab']
        trgScanID = cfgDict['trgScanID']
        trgSlcNum = cfgDict['trgSlcNum']
        trgRoicolName = cfgDict['trgRoicolName']
        trgRoiName = cfgDict['trgRoiName']
        
        items = [
            projID, subjLab, srcExpLab, srcScanID, srcSlcNum,\
            srcRoicolName, srcRoiName, roicolMod, trgExpLab,\
            trgScanID, trgSlcNum, trgRoicolName, trgRoiName
            ]
        
        names = [
            'projID', 'subjLab', 'srcExpLab', 'srcScanID', 'srcSlcNum',\
            'srcRoicolName', 'srcRoiName', 'roicolMod', 'trgExpLab',\
            'trgScanID', 'trgSlcNum', 'trgRoicolName', 'trgRoiName'
            ]
        
        maxL = max([len(name) for name in names])
        
        print(f'* Fetching data from XNAT for runID = {runID}')
        """
        print(f'  projID = {projID}\n  subjLab = {subjLab}\n  srcExpLab ='
              f'{srcExpLab}\n  srcScanID = {srcScanID}\n  srcSlcNum =',
              f'{srcSlcNum}\n  srcRoicolName = {srcRoicolName}\n  srcRoiName',
              f'= {srcRoiName}\n  roicolMod = {roicolMod}\n  trgExpLab =',
              f'{trgExpLab}\n  trgScanID = {trgScanID}\n  trgSlcNum =',
              f'{trgSlcNum}\n  trgRoicolName = {trgRoicolName}\n  trgRoiName',
              f'= {trgRoiName}\n')
        """
        for name, item in zip(names, items):
            L = len(name)
            d = maxL - L
            print(' '*d + f'{name} = {item}')
        print('')
    
    def download_and_get_pathsDict(self):
        # TODO update docstrings
        """
        17/08: xnatSession was input.
        
        Downloads scans and ROI Collections and creates a dictionary containing
        filepaths.
        
        Returns
        -------
        self.pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        #xnatSession = self.xnatSession
        #cfgDict = self.cfgDict
        #p2c = cfgDict['p2c']
        
        self.print_cfg_params_to_console()
        
        print('*** Fetching source DICOM scan from XNAT...\n')
        
        self.pathsDict = download_scan(
            config=self.cfgDict, srcORtrg='src', xnatSession=self.xnatSession
            )
        
        print('*** Fetching target DICOM scan from XNAT...\n')
        
        self.pathsDict = download_scan(
            config=self.cfgDict, srcORtrg='trg', xnatSession=self.xnatSession,
            pathsDict=self.pathsDict
            )
        
        print('*** Fetching source ROI Collection from XNAT...\n')
        
        self.pathsDict = download_im_asr(
            config=self.cfgDict, srcORtrg='src', xnatSession=self.xnatSession,
            pathsDict=self.pathsDict
            )
        
        trgRoicolName = self.cfgDict['trgRoicolName']
        
        if trgRoicolName != None:
            print('*** Fetching target ROI Collection from XNAT...\n')
            
            self.pathsDict = download_im_asr(
                config=self.cfgDict, srcORtrg='trg', 
                xnatSession=self.xnatSession, pathsDict=self.pathsDict
                )
        
        #self.pathsDict = pathsDict
        #self.cfgDict = cfgDict
        
        if trgRoicolName == None:
            timingMsg = "Took [*] s to download the source and target DICOM"\
                + " scans and the source ROI Collection.\n"
        else:
            timingMsg = "Took [*] s to download the source and target DICOM"\
                + " scans and the source and target ROI Collections.\n"
        self.add_timestamp(timingMsg)
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        """
        Moved to image_tools.propagate.py:
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(config=self.cfgDict)
        
        # Update cfgDict:
        self.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        self.cfgDict['useCaseToApply'] = useCaseToApply
        
        
        if self.cfgDict['p2c']:
            print(f"cfgDict['runID'] = {self.cfgDict['runID']}")
            print(f'useCaseThatApplies = {useCaseThatApplies}')
            print(f'useCaseToApply = {useCaseToApply}')
            #print(f"cfgDict = {self.cfgDict}")
        """