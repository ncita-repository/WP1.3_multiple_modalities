# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 15:54:12 2021

@author: ctorti
"""


def get_scan_asr_fname_and_id(pathsDict, projID, subjLab, expLab, asrMod, 
                              asrName):
    """
    Get the filename and scan ID corresponding to a scan assessor given a
    name and modality.

    Parameters
    ----------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    projID : str
        The project ID of interest.
    subjLab : str
        The subject label of interest.
    expLab : str
        The DICOM study / XNAT experiment label of interest.
    asrMod : str
        The modality of the scan assessor.
    asrName : str
        The name of the scan assessor.

    Returns
    -------
    asrFname : str
        File name of the scan assessor.
    asrID : str
        ID of the scan assessor.
    """
    
    # Get all IDs:
    IDs = list(pathsDict['projects'][projID]['subjects'][subjLab]\
               ['experiments'][expLab]['assessors'].keys()
               )

    asrFname = None
    asrID = None

    for ID in IDs:
        fnames = list(pathsDict['projects'][projID]['subjects'][subjLab]\
                      ['experiments'][expLab]['assessors'][ID]\
                      ['resources'][asrMod]['files'].keys()
                      )

        for fname in fnames:
            name = pathsDict['projects'][projID]['subjects'][subjLab]\
                            ['experiments'][expLab]['assessors'][ID]\
                            ['resources'][asrMod]['files'][fname]['name']
            
            if name == asrName:
                asrFname = fname
                asrID = ID
        
    #if ScanAsrID == None:
    #    msg = f"There was no match on ScanAsrName = {ScanAsrName}."
    #    
    #    raise Exception(msg)
    
    return asrFname, asrID
