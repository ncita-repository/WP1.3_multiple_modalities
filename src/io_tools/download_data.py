# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 10:31:20 2021

@author: ctorti
"""

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
"""

import os
import time
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
from xnat_tools.format_pathsDict import get_scan_asr_fname_and_id
#from xnat_tools.alias_tokens import (
#    import_alias_token, is_alias_token_valid, generate_alias_token,
#    export_alias_token
#    )
from xnat_tools.alias_tokens import (
    import_valid_alias_token, generate_alias_token, export_alias_token
    )


class DataDownloader:
    # TODO! Update docstrings
    """
    This class downloads data from XNAT and creates stores the file paths in a 
    dictionary.
    
    Parameters
    ----------
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
    
    def __init__(self, cfgObj):
        self.cfgDict = cfgObj.cfgDict
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
        tokenDir = os.path.join(cwd, 'xnat_tokens')
        
        if aliasToken == {}:
            print(f'Searching for a valid XNAT alias token in {tokenDir}')
            
            # Import a valid XNAT alias token (or return {}):
            aliasToken = import_valid_alias_token(tokenDir, url)
        
        if aliasToken:
            print('Using XNAT alias token to establish XNAT connection.')
            xnatSession = create_session(url=url, aliasToken=aliasToken)
        else:
            print('No valid XNAT alias token was found.')
            print('Establishing XNAT connection without an XNAT alias token.')
            xnatSession = create_session(url=url, username=username)
        
            # Generate an XNAT alias token to avoid making a new authenticated
            # user session and export to tokenDir:
            aliasToken = generate_alias_token(xnatSession)
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
        
        if 'total' in timingMsg:
            dTime = self.timings[-1] - self.timings[0]
        else:
            dTime = self.timings[-1] - self.timings[-2]
        
        if '[*]' in timingMsg:
            if 'total' in timingMsg:
                timingMsg = timingMsg.replace(
                    '[*]', f'{dTime:.2f} s ({dTime/60:.1f} min)'
                    )
            else:
                timingMsg = timingMsg.replace('[*]', f'{dTime:.2f} s')
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
        Downloads scans and ROI Collections and creates a dictionary containing
        filepaths.
        
        Returns
        -------
        self.pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
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
        
        if trgRoicolName == None:
            timingMsg = "Took [*] to download the source and target DICOM"\
                + " scans and the source ROI Collection.\n"
        else:
            timingMsg = "Took [*] to download the source and target DICOM"\
                + " scans and the source and target ROI Collections.\n"
        self.add_timestamp(timingMsg)
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()