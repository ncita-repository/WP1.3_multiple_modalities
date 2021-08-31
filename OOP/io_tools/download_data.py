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
#import io_tools.inputs_checker
#reload(io_tools.inputs_checker)
#import xnat_tools.dros
#reload(xnat_tools.dros)
#import io_tools.fetch_config
#reload(io_tools.fetch_config)

import time
#from io_tools.inputs_checker import which_use_case
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
from xnat_tools.format_pathsDict import get_scan_asr_fname_and_id
#from xnat_tools.dros import search_dro


class DataDownloader:
    """
    This class downloads data from XNAT and creates stores the file paths in a 
    dictionary.
    
    Parameters
    ----------
    cfgDict : dict
        A dictionary containing configuration settings (parameters) for the run
        to be performed.
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
    
    def __init__(self, cfgDict, xnatSession=None):
        
        self.cfgDict = cfgDict
        
        if xnatSession == None:
            # Establish XNAT connection:
            xnatSession = create_session(xnatCfg=cfgDict)
            
        self.xnatSession = xnatSession
        
        # Initialise list of timestamps and timing messages to be stored:
        self.timings = [time.time()]
        self.timingMsgs = []
    
    def __init__210817(self, cfgDict, xnatSession=None):
        
        p2c = cfgDict['p2c']
        
        self.cfgDict = cfgDict
        
        if xnatSession == None:
            # Establish XNAT connection:
            xnatSession = create_session(xnatCfg=cfgDict)
            
        self.xnatSession = xnatSession
        
        # Download data and create pathsDict:
        self.get_pathsDict(xnatSession=xnatSession)
        #self.pathsDict = self.get_pathsDict(xnatSession)
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(config=self.cfgDict)
        
        if p2c:
            print(f"cfgDict['runID'] = {cfgDict['runID']}")
            print(f'useCaseThatApplies = {useCaseThatApplies}')
            print(f'useCaseToApply = {useCaseToApply}')
        
        # Update cfgDict:
        self.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        self.cfgDict['useCaseToApply'] = useCaseToApply
    
    
    
    def get_pathsDict_Pre_210809(self, xnatSession):
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
        
    def add_dirs_and_fpaths_210810(self):
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
        
        print(f'* Fetching data from XNAT for runID = {runID}')
        print(f'  projID = {projID}\n  subjLab = {subjLab}\n  srcExpLab ='
              f'{srcExpLab}\n  srcScanID = {srcScanID}\n  srcSlcNum =',
              f'{srcSlcNum}\n  srcRoicolName = {srcRoicolName}\n  srcRoiName',
              f'= {srcRoiName}\n  roicolMod = {roicolMod}\n  trgExpLab =',
              f'{trgExpLab}\n  trgScanID = {trgScanID}\n  trgSlcNum =',
              f'{trgSlcNum}\n  trgRoicolName = {trgRoicolName}\n  trgRoiName',
              f'= {trgRoiName}\n')
        
    def download_210818(self):
        # TODO update docstrings
        """
        Download scans and ROI Collections, and create dictionary containing
        filepaths (pathsDict).
        
        Parameters
        ----------
        xnatSession : Requests Object
            Requests session of XNAT connection.
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        p2c = self.cfgDict['p2c']
        
        # Download data and create pathsDict:
        self.get_pathsDict()
        #self.pathsDict = self.get_pathsDict()
        
        # Modify Source and Target XNAT parameters:
        self.add_dirs_and_fpaths()
        
        # Determine which use case applies and which to apply:
        useCaseThatApplies, useCaseToApply\
            = which_use_case(config=self.cfgDict)
        
        # Update cfgDict:
        self.cfgDict['useCaseThatApplies'] = useCaseThatApplies
        self.cfgDict['useCaseToApply'] = useCaseToApply
        
        if p2c:
            print(f"cfgDict['runID'] = {self.cfgDict['runID']}")
            print(f'useCaseThatApplies = {useCaseThatApplies}')
            print(f'useCaseToApply = {useCaseToApply}')
            #print(f"cfgDict = {self.cfgDict}")
    
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