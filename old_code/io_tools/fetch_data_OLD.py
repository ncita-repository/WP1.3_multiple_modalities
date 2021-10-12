# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:22:27 2021

@author: ctorti
"""

from importlib import reload
import io_tools.fetch_params
reload(io_tools.fetch_params)

import xnat_tools.scans
reload(xnat_tools.scans)
import xnat_tools.im_assessors
reload(xnat_tools.im_assessors)

from io_tools.fetch_params import Params
from xnat_tools.sessions import create_session
from xnat_tools.scans import download_scan
from xnat_tools.im_assessors import download_im_asr
from xnat_tools.format_paths_dict import get_scan_asr_fname_and_id
#from xnat_tools.dro import search_dro


class PathsDict:
    """
    This class downloads data from XNAT and stores filepaths into a dictionary.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    xnatSession : requests.models.Response, optional
        If provided a new session request will be avoided.
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    """
    
    def __init__(self, cfgDir, xnatSession=None):
        #self.params = Params(cfgDir)
        ##self.Params(cfgDir)
        #
        #print(f'dir(self) = {dir(self)}\n')
        #print(f'self.params = {self.params}\n')
        #print(f'self.params.xnatCfg = {self.params.xnatCfg}\n')
        
        params = Params(cfgDir)
        xnatCfg = params.xnatCfg
        
        """ I expected these variables to be accessible to get_pathsDict()
        but they aren't for some reason...
        xnatCfg = params.xnatCfg
        srcXnatParams = params.srcXnatParams
        trgXnatParams = params.trgXnatParams
        genParams = params.genParams
        """
        
        if xnatSession == None:
            ##self.xnatSession = create_session(self.xnatCfg)
            #self.xnatSession = create_session(self.params.xnatCfg)
            
            #params = Params(cfgDir)
            #self.xnatSession = create_session(params.get_xnat_cfg())
            #self.xnatSession = create_session(xnatCfg)
            xnatSession = create_session(xnatCfg)
        #else:
        #    self.xnatSession = xnatSession
        
        ###self.download_scans_and_roicols(self)
        ##self.download_scans_and_roicols()
        #self.download_scans_and_roicols(params=params) # pathsDict -> 
        ## io_tools.fetch_data.PathsDict object
        
        #self.params = params
        self.get_pathsDict(params, xnatSession)
    
    def get_pathsDict(self, params, xnatSession):
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
        
        #params = self.params
        
        #xnatSession = self.xnatSession
        #xnatCfg = params.get_xnat_cfg()
        #srcXnatParams = params.get_xnat_params('src')
        #trgXnatParams = params.get_xnat_params('trg')
        #genParams = params.get_gen_params()
        
        xnatCfg = params.xnatCfg
        srcXnatParams = params.srcXnatParams
        trgXnatParams = params.trgXnatParams
        genParams = params.genParams
        
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
    
    def add_dirs_and_fpaths(self, params):
        """
        Add directories and filepaths to params.srcXnatParams and 
        params.trgXnatParams from pathsDict.

        Returns
        -------
        None.

        """
        pathsDict = self.pathsDict
        
        srcXnatParams = params.srcXnatParams
        trgXnatParams = params.trgXnatParams
        
        projId = srcXnatParams['projId']
        subjLab = srcXnatParams['subjLab']
        srcExpLab = srcXnatParams['srcExpLab']
        srcScanId = srcXnatParams['srcScanId']
        srcRoicolMod = srcXnatParams['srcRoicolMod']
        srcRoicolName = srcXnatParams['srcRoicolName']
        trgExpLab = trgXnatParams['trgExpLab']
        trgScanId = trgXnatParams['trgScanId']
        trgRoicolMod = trgXnatParams['trgRoicolMod']
        trgRoicolName = trgXnatParams['trgRoicolName']
        
        srcDcmDir = pathsDict['projects'][projId]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['scans'][srcScanId]['resources']\
                ['DICOM']['files']['dir']
        
        trgDcmDir = pathsDict['projects'][projId]['subjects'][subjLab]\
            ['experiments'][trgExpLab]['scans'][trgScanId]['resources']\
                ['DICOM']['files']['dir']
        
        srcRoicolFname, srcAsrId\
            = get_scan_asr_fname_and_id(pathsDict, projId, subjLab, srcExpLab,
                                        srcRoicolMod, srcRoicolName)
        
        srcRoicolFpath = pathsDict['projects'][projId]['subjects'][subjLab]\
            ['experiments'][srcExpLab]['assessors'][srcAsrId]['resources']\
                [srcRoicolMod]['files'][srcRoicolFname]['fpath']
        
        if trgRoicolName == None:
            trgRoicolFpath = None
        else:
            trgRoicolFname, trgAsrId\
                = get_scan_asr_fname_and_id(pathsDict, projId, subjLab, 
                                            trgExpLab, trgRoicolMod, 
                                            trgRoicolName)
        
            trgRoicolFpath = pathsDict['projects'][projId]['subjects'][subjLab]\
                ['experiments'][trgExpLab]['assessors'][trgAsrId]['resources']\
                    [trgRoicolMod]['files'][trgRoicolFname]['fpath']
        
        # Add filepaths and directories to srcXnatParams and trgXnatParams:
        srcXnatParams['srcAsrId'] = srcAsrId
        srcXnatParams['srcRoicolFname'] = srcRoicolFname
        srcXnatParams['srcRoicolFpath'] = srcRoicolFpath
        srcXnatParams['srcDcmDir'] = srcDcmDir
        
        params.srcXnatParams = srcXnatParams
        
        if trgRoicolName != None:
            trgXnatParams['trgAsrId'] = trgAsrId
            trgXnatParams['trgRoicolFname'] = trgRoicolFname
            trgXnatParams['trgRoicolFpath'] = trgRoicolFpath
            trgXnatParams['trgDcmDir'] = trgDcmDir
            
            params.trgXnatParams = trgXnatParams
        
        self.pathsDict = pathsDict
    
    """
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            raise AttributeError("No such attribute: " + name)
    
    def __setattr__(self, name, value):
        #self[name] = value
        # https://stackoverflow.com/questions/17020115/how-to-use-setattr-correctly-avoiding-infinite-recursion
        self.__dict__[name] = value

    def __delattr__(self, name):
        if name in self.__dict__:
            del self.__dict__[name]
        else:
            raise AttributeError("No such attribute: " + name)
    """
    
    def download_scans_and_roicols(self, params):
        """
        14/07: Delete this
        
        Downloads scans and ROI Collections and returns a dictionary containing
        filepaths.
        
        Parameters
        ----------
        params : Params Object
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        #"""
        
        params = self.params
        
        xnatSession = self.xnatSession
        xnatCfg = params.get_xnat_cfg()
        srcXnatParams = params.get_xnat_params('src')
        trgXnatParams = params.get_xnat_params('trg')
        genParams = params.get_gen_params()
        
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
        
        return pathsDict
    
    def download_scans_and_roicols_OLD(self):
        """
        14/07: Delete this
        
        Downloads scans and ROI Collections and returns a dictionary containing
        filepaths.
        
        Parameters
        ----------
        
        
        Returns
        -------
        pathsDict : dict
            Dictionary containing paths of data downloaded.
        """
        
        params = self.params
        xnatSession = self.xnatSession
        
        xnatCfg = params.xnatCfg
        srcXnatParams = params.srcXnatParams
        trgXnatParams = params.trgXnatParams
        genParams = params.genParams
        
        pathsDict = download_scan(xnatCfg=xnatCfg, xnatParams=srcXnatParams, 
                                  genParams=genParams, xnatSession=xnatSession)
    
        pathsDict = download_scan(xnatCfg=xnatCfg, xnatParams=trgXnatParams, 
                                  genParams=genParams, xnatSession=xnatSession,
                                  pathsDict=pathsDict)
    
        pathsDict = download_im_asr(xnatCfg=xnatCfg, xnatParams=srcXnatParams,
                                    genParams=genParams, xnatSession=xnatSession,
                                    pathsDict=pathsDict)
        
        trgRoicolName = trgXnatParams['roicolName']
        
        if trgRoicolName != None:
            pathsDict = download_im_asr(xnatCfg=xnatCfg, xnatParams=trgXnatParams,
                                        genParams=genParams, 
                                        xnatSession=xnatSession,
                                        pathsDict=pathsDict)
            
        self.pathsDict = pathsDict
        
