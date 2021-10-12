# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 16:28:41 2021

@author: ctorti
"""


def create_pathsDict_for_exp(projID, subjLab, expLab, pathsDict=None):
    """ 
    Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level. 
    """
    
    if pathsDict == None:
        pathsDict = {}
        keys = []
    else:
        keys = list(pathsDict.keys())
    
    
    if not 'projects' in keys:
        pathsDict.update({'projects' : {}})
    
    keys = list(pathsDict['projects'].keys())
    
    if not projID in keys:
        pathsDict['projects'].update({projID : {}})
    
    keys = list(pathsDict['projects'][projID].keys())
    
    if not 'subjects' in keys:
        pathsDict['projects'][projID].update({'subjects' : {}})
    
    keys = list(pathsDict['projects'][projID]['subjects'].keys())
    
    if not subjLab in keys:
        pathsDict['projects'][projID]['subjects'].update({subjLab : {}})
    
    keys = list(pathsDict['projects'][projID]['subjects'][subjLab].keys())
    
    if not 'experiments' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]\
            .update({'experiments' : {}})
        
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            .keys()
                )
    
    if not expLab in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            .update({expLab : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab].keys())
    
    return pathsDict, keys

def create_pathsDict_for_scan(projID, subjLab, expLab, scanID, pathsDict=None):
    """ 
    Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level. 
    """
    
    pathsDict, keys\
        = create_pathsDict_for_exp(projID, subjLab, expLab, pathsDict)
    
    if not 'scans' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab].update({'scans' : {}})
        
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'].keys()
                )
    
    if not scanID in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'].update({scanID : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID].keys()
                )
    
    if not 'resources' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID].update({'resources' : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID]['resources'].keys()
                )
    
    if not 'DICOM' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID]['resources'].update({'DICOM' : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID]['resources']['DICOM'].keys()
                )
                
    
    if not 'files' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID]['resources']['DICOM']\
                .update({'files' : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['scans'][scanID]['resources']['DICOM']['files'].keys()
                )
    
    return pathsDict, keys

def create_pathsDict_for_im_asr(
        projID, subjLab, expLab, asrID, mod, fname, pathsDict=None
        ):
    """ 
    Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level.
    """
    
    pathsDict, keys\
        = create_pathsDict_for_exp(projID, subjLab, expLab, pathsDict)
    
    if not 'assessors' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab].update({'assessors' : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'].keys())
    
    if not asrID in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'].update({asrID : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID].keys())
    
    if not 'resources' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID].update({'resources' : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID]['resources'].keys())
    
    if not mod in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID]['resources'].update({mod : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID]['resources'][mod].keys()
                )
    
    if not 'files' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID]['resources'][mod]\
                .update({'files' : {}})
    
    keys = list(
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID]['resources'][mod]['files'].keys()
                )
    
    if not fname in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
            [expLab]['assessors'][asrID]['resources'][mod]['files']\
                .update({fname : {}})
    
    return pathsDict, keys

def get_scan_asr_fname_and_id(
        pathsDict, projID, subjLab, expLab, asrMod, asrName
        ):
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
               ['experiments'][expLab]['assessors'].keys())

    asrFname = None
    asrID = None

    for ID in IDs:
        fnames = list(
            pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
                [expLab]['assessors'][ID]['resources'][asrMod]['files'].keys()
                      )

        for fname in fnames:
            name = pathsDict['projects'][projID]['subjects'][subjLab]\
                ['experiments'][expLab]['assessors'][ID]['resources'][asrMod]\
                    ['files'][fname]['roicolName']
            
            if name == asrName:
                asrFname = fname
                asrID = ID
        
    #if ScanAsrID == None:
    #    msg = f"There was no match on ScanAsrName = {ScanAsrName}."
    #    
    #    raise Exception(msg)
    
    return asrFname, asrID

def create_pathsDict_for_subj_asr(
        projID, subjLab, resID, name, pathsDict=None
        ):
    """ 
    Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level. 
    """
    
    if pathsDict == None:
        pathsDict = {}
        keys = []
    else:
        keys = list(pathsDict.keys())
    
    if not 'projects' in keys:
        pathsDict.update({'projects' : {}})
    
    keys = list(pathsDict['projects'].keys())
    
    if not projID in keys:
        pathsDict['projects'].update({projID : {}})
    
    keys = list(pathsDict['projects'][projID].keys())
    
    if not 'subjects' in keys:
        pathsDict['projects'][projID].update({'subjects' : {}})
    
    keys = list(pathsDict['projects'][projID]['subjects'].keys())
    
    if not subjLab in keys:
        pathsDict['projects'][projID]['subjects'].update({subjLab : {}})
    
    keys = list(pathsDict['projects'][projID]['subjects'][subjLab].keys())
    
    if not 'resources' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]\
        .update({'resources' : {}})
        
    keys = list(pathsDict['projects'][projID]['subjects'][subjLab]\
                ['resources'].keys())
    
    if not resID in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['resources']\
        .update({resID : {}})
    
    keys = list(pathsDict['projects'][projID]['subjects'][subjLab]\
                ['resources'][resID].keys())
    
    if not 'files' in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['resources']\
        [resID].update({'files' : {}})
        
    keys = list(pathsDict['projects'][projID]['subjects'][subjLab]\
                ['resources'][resID]['files'].keys())
    
    if not name in keys:
        pathsDict['projects'][projID]['subjects'][subjLab]['resources']\
        [resID]['files'].update({name : {}})
    
    keys = list(pathsDict['projects'][projID]['subjects'][subjLab]\
                ['resources'][resID]['files'][name].keys())
    
    return pathsDict, keys

def reorder_keys_and_fill_zeros(data_by_proj):
    """
    Reorder the keys in the dictionary data_by_proj and zero-fill missing keys.
    
    E.g. if there are no PET scans for a particular project, there will be no
    'PET' key, so add it and assign the number 0.
    
    Parameters
    ----------
    data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Returns
    -------
    new_data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    Ensures that PI comes before CIs.
    Orders image sessions as follows: MR, CT, PET, SC, RTIMAGE.
    
    It isn't clear what are the image session types available within XNAT. At
    present that info is not available, e.g. see:
    https://wiki.xnat.org/documentation/xnat-administration/managing-data-types-in-xnat/adding-and-managing-image-scan-data-types
    
    so will define a list containing keys in the desired order. Any items not
    included in the list will follow whatever order they were included in the
    original dictionary.
    """
    
    new_order = [
        'PI', 'Project description', 'Co-investigators', 
        'No. of users', 'Users', 
        'First upload', 'Last upload', 
        'No. of subjects', 'No. of experiments', 
        'No. of MR sessions', 'No. of CT sessions', 'No. of PET sessions',
        'Ave. image scans/session', 
        'No. of MR scans', 'No. of CT scans', 'No. of PET scans',
        'No. of OT scans',
        'No. of DICOM image files [k]', 'Size of DICOM image files [MB]',
        'No. of RTSTRUCT scans', #'Size of RTSTRUCT scans [MB]',
        'No. of RTSTRUCT ROI files', #'Size of RTSTRUCT ROI files [MB]',
        'No. of SEG ROI files', #'Size of SEG ROI files [MB]',
        'No. of DICOM files [k]', 'Size of DICOM files [MB]'
        ]
    
    #print(f'new_order = {new_order}\n')
    
    new_data_by_proj = {}
    
    for project in data_by_proj.keys():
        new_data_by_proj[project] = {}
        
        old_order = list(data_by_proj[project].keys())
        
        #print(f'old_order = {old_order}\n')
        
        # First add the keys in the order they appear in new_order:
        for key in new_order:
            if key in old_order:
                new_data_by_proj[project][key] = data_by_proj[project][key]
            else:
                #print(f"key '{key}' not in new_order\n")
                new_data_by_proj[project][key] = 0
                #new_data_by_proj[project].update({key : 0})
        
        # Add any keys in old_order that were not in new_order:
        for key in old_order:
            if not key in new_order:
                new_data_by_proj[project][key] = data_by_proj[project][key]
    
    return new_data_by_proj