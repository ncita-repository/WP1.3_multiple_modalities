# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:30:55 2021

@author: ctorti
"""

from importlib import reload

import io_tools.general
reload(io_tools.general)

import general_tools.general
reload(general_tools.general)

import os
from pathlib import Path
from pydicom import dcmread
from io import BytesIO
from xnat_tools.sessions import create_session
from xnat_tools.format_pathsDict import create_pathsDict_for_im_asr
from io_tools.general import get_user_input_as_int
from general_tools.general import (
    combine_dates_and_times, get_ind_of_newest_and_oldest_datetime
    )

def download_im_asr(
        config, srcORtrg, xnatSession=None, pathsDict=None#, 
        #whichSrcRoicol='user'
        ):
    """
    Download an image assessor (ROI Collection) of interest for a particular 
    scan from XNAT.
    
    Parameters
    ----------
    config : dict
        Dictionary of configuration settings.
    srcORtrg : str
            'src' or 'trg' for Source or Target.
    xnatSession : requests.models.Response, optional
        If provided a new session request will be avoided. The default value
        is None.
    pathsDict : dict, optional
        Dictionary containing paths of data downloaded. The default value is
        None.
    whichSrcRoicol : str, optional
        Determine how to proceed if multiple image assessors (ROI Collections)
        are found that meet the criteria.  Allowed values are 'oldest' -> the
        oldest file, 'newest' -> The newest file, and 'user' -> Prompt the user
        to decide. The default value is 'user'.
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    
    Notes
    -----
    Although the structure of the REST call below:
    
    assessor = session.get(f'{url}/data/projects/{projID}/'
                           f'subjects/{subjLab}/'
                           f'experiments/{expLab}/'
                           f'assessors/{AsrID}/resources/{AsrMod}/files')
    
    follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Also, Jason suggested adding the file name after 'files', which worked. 
    """
    
    if xnatSession == None:
        xnatSession = create_session(xnatCfg=config)
    
    url = xnatSession.url
    
    rootExportDir = config['rootExportDir']
    
    projID = config['projID']
    subjLab = config['subjLab']
    roicolMod = config['roicolMod']
    if srcORtrg == 'src':
        expLab = config['srcExpLab']
        scanID = config['srcScanID']
        roicolName_req = config['srcRoicolName']
    else:
        expLab = config['trgExpLab']
        scanID = config['trgScanID']
        roicolName_req = config['trgRoicolName']
    
    whichSrcRoicol = config['whichSrcRoicol']
    p2c = config['p2c']
    
    if p2c:
        print('Inputs to download_im_asr():')
        print(f'  projID = {projID}')
        print(f'  subjLab = {subjLab}')
        print(f'  expLab = {expLab}')
        print(f'  scanID = {scanID}')
        print(f'  roicolMod = {roicolMod}')
        print(f'  roicolName_req = {roicolName_req}')
        print(f'  whichSrcRoicol = {whichSrcRoicol}\n')
    
    """ 
    The 'collectionType' key in 'data_fields' is either 'AIM' or 'SEG'. The 
    required collection type, collType_req, is 'AIM' for roicolMod = 'RTSTRUCT'
    and 'SEG' for roicolMod = 'SEG'.
    """
    if roicolMod == 'RTSTRUCT':
        #collType_req = 'AIM' # 10/09/21 I can't remember why I made it AIM
        collType_req = 'RTSTRUCT' # 10/09/21
    else:
        collType_req = 'SEG'
    
    """ Get the experiment of interest: """
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}?format=json'
                              
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    experiment = request.json()
    
    """ 
    There are likely more than one item in the list 
    experiment['items'][0]['children'] (e.g. one for scans and one for 
    image assessors). Get the index that corresponds to assessors, asr_ind, and
    the index that corresponds to scans, scan_ind:
    """
    asrInd = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'assessor' in experiment['items'][0]['children'][i]['field']:
            asrInd = i
    
    if asrInd == None:
        msg = f"There are no assessors for:\n  projID = '{projID}'"\
            + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
            + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    scanInd = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'scan' in experiment['items'][0]['children'][i]['field']:
            scanInd = i
    
    if scanInd == None:
        msg = f"There are no scans for:\n  projID = '{projID}'"\
            + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
            + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    if p2c:
        print(f'asrInd = {asrInd}')
        print(f'scanInd = {scanInd}\n')
    
    """ Get the SeriesUID for the scan of interest: """
    seriesUID_req = None # initial value
    
    if pathsDict != None:
        try:
            seriesUID_req = \
                pathsDict['projects'][projID]['subjects'][subjLab]\
                    ['experiments'][expLab]['scans'][scanID]['seriesUID']
            
            excRaised = False # exception raised
            
        except Exception:
            excRaised = True
    
    if pathsDict == None or excRaised:
        for scan in experiment['items'][0]['children'][scanInd]['items']:
            if p2c:
                print(f"ID = {scan['data_fields']['ID']}")
                print(f"UID = {scan['data_fields']['UID']}\n")
            
            if scan['data_fields']['ID'] == scanID:
                seriesUID_req = scan['data_fields']['UID']
                
                if p2c:
                    print(f'seriesUID_req = {seriesUID_req}\n')
    
    if seriesUID_req == None:
        msg = f"No image assessors were found for:\n  projID = '{projID}'"\
            + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
            + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    if p2c:
        print('Items to match:')
        print(f'   seriesUID_req = {seriesUID_req}')
        print(f'   collType_req = {collType_req}')
        print(f'   roicolName_req = {roicolName_req}\n')
    
    """ Get the assessor ID and label of interest: """
    #roicolInd = None # initial value
    #asrID = None # initial value
    #label = None # initial value
    #asrDate = None # initial value
    #asrTime = None # initial value
    matchingInds = [] # initial value
    matchingIDs = [] # initial value
    matchingLabels = [] # initial value
    matchingNames = [] # initial value
    matchingDates = [] # initial value
    matchingTimes = [] # initial value
    
    for i in range(len(experiment['items'][0]['children'][asrInd]['items'])):
        if p2c:
            print(f"Item no {i}:")
        
        children = experiment['items'][0]['children'][asrInd]['items'][i]\
            ['children']
        
        if p2c:
            print(f"  children has {len(children)} items")
        
        fields = [child['field'] for child in children]
        
        # Get the index corresponding to the item containing the fields
        # 'out/file' and 'references/seriesUID':
        if 'out/file' in fields:
            fileInd = fields.index('out/file')
        else:
            fileInd = None
        
        if 'references/seriesUID' in fields:
            refInd = fields.index('references/seriesUID')
        else:
            refInd = None
        
        if p2c:
            print(f"  fields = {fields}")
            print(f"  fileInd = {fileInd}")
            print(f"  refInd = {refInd}\n")
        
        """
        The following check is probably not required in most cases, but 
        covers cases where a ROI collection has been deleted and although the
        item with field 'out/file' is removed from the JSON, the item with
        field 'references/seriesUID' can persist.
        """
        # Check that items exist with both 'out/file' and 
        # 'references/seriesUID' fields:
        #if 'out/file' in fields and 'references/seriesUID' in fields:
        if isinstance(fileInd, int) and isinstance(refInd, int):
            if p2c:
                print("  This item's children contains items with fields",
                      "'out/file' and 'references/seriesUID'.\n")
            
            seriesUid = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['children'][refInd]['items'][0]['data_fields']['seriesUID']
            
            collType = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['data_fields']['collectionType']
            
            roicolName = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['data_fields']['name']
            
            if p2c:
                print(f"  seriesUid = {seriesUid}")
                print(f"  collType = {collType}")
                print(f"  roicolName = {roicolName}\n")
            
                #if seriesUid == seriesUID_req:
                #    print(f"   * Matched collType = {collType}")
                #    print(f"   * Matched roicolName = {roicolName}\n")
            
            if (seriesUid == seriesUID_req and collType == collType_req and 
                    roicolName == roicolName_req):
                
                asrID = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['id']
                
                label = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['label']
                
                asrDate = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['date']
                
                asrTime = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['time']
                
                #roicolInd = i
                matchingInds.append(i)
                matchingIDs.append(asrID)
                matchingLabels.append(label)
                matchingNames.append(roicolName)
                matchingDates.append(asrDate)
                matchingTimes.append(asrTime)
                
                if p2c:
                    print(f"   * Match for ROI Collection {i}")
                    print(f"   * Matched seriesUid = {seriesUid}")
                    print(f"   * Matched collType = {collType}")
                    print(f"   * Matched roicolName = {roicolName}\n")
            else:
                if p2c:
                    print('   * Not a match.\n')
        else:
            if p2c:
                print("  This item's children DOES NOT contain items with",
                      "both field 'out/file' and 'references/seriesUID'.\n")
    
    if len(matchingInds) == 0:
        msg = f"No image assessors were found for:\n  projID = '{projID}'" +\
            f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'" +\
            f"\n  scanID = '{scanID}'\ncontaining matches on:" +\
            f"\n  roicolMod = '{roicolMod}'\n  roicolName = " +\
            f"{roicolName_req}."
        raise Exception(msg)
    
    if p2c:
        print(f'Matching indices = {matchingInds}')
    
    if len(matchingInds) > 1:
        msg = "More than one matching ROI Collection was found. "\
            + "Following are their timestamps:"
        
        for i in range(len(matchingInds)):
            msg += f"\n{i+1}. {matchingDates[i]} {matchingTimes[i]}"
        
        if whichSrcRoicol == 'user':
            msg += "\nEnter an integer corresponding to the ROI Collection to be "\
                + "imported"
            
            choice = get_user_input_as_int(
                message=msg, minVal=1, maxVal=len(matchingInds)
                )
            
            roicolInd = choice - 1 # userInput is 1-indexed
            
            #print(f'choice = {choice}, roicolInd = {roicolInd}\n')
        else:
            matchingDateTimes = combine_dates_and_times(
                matchingDates, matchingTimes
                )
            
            newestInd, oldestInd = get_ind_of_newest_and_oldest_datetime(
                matchingDateTimes
                )
            
            if whichSrcRoicol == 'oldest':
                roicolInd = int(oldestInd)
            else:
                roicolInd = int(newestInd)
             
        asrID = matchingIDs[roicolInd]
        label = matchingLabels[roicolInd]
        roicolName = matchingNames[roicolInd]
        asrDate = matchingDates[roicolInd]
        asrTime = matchingTimes[roicolInd]
    else:
        roicolInd = matchingInds[0]
        
        asrID = matchingIDs[0]
        label = matchingLabels[0]
        roicolName = matchingNames[0]
        asrDate = matchingDates[0]
        asrTime = matchingTimes[0]
    
    #if p2c:
    #    print(f'roicolInd = {roicolInd}\n')
    
    #asrID = matchingIDs[roicolInd]
    #label = matchingLabels[roicolInd]
    #roicolName = matchingNames[roicolInd]
    #asrDate = matchingDates[roicolInd]
    #asrTime = matchingTimes[roicolInd]
    
    if p2c:
        print(f"Index of the ROI Collection of interest = {roicolInd}")
        print(f"Name of the ROI Collection of interest = {roicolName}")
        print(f"ID of the ROI Collection of interest = {asrID}")
        print(f"Label of the ROI Collection of interest = {label}")
        print(f"Date and time of the ROI Collection of interest = {asrDate}",
              f"{asrTime}")
    
    """ Moved above:
    
    if roicolInd == None:
        msg = f"No image assessors were found for:\n  projID = '{projID}'"\
            + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
            + f"\n  scanID = '{scanID}'\ncontaining matches on:"\
            + f"\n  roicolMod = '{roicolMod}'\n  roicolName_req ="\
            + f"{roicolName_req}."
        raise Exception(msg)
    """
    
    """ Get the assessor of interest: """
    fname = label + '.dcm'
    
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/assessors/{asrID}/'\
          + f'resources/{roicolMod}/files/{fname}'
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    exportDir = os.path.join(
        rootExportDir, 'projects', projID, 'subjects', subjLab,
        'experiments', expLab, 'assessors', asrID, 'resources', roicolMod,
        'files'
        )

    if not os.path.isdir(exportDir):
        Path(exportDir).mkdir(parents=True)
        print(f'Created directory:\n {exportDir}\n')
    
    fpath = os.path.join(exportDir, fname)
    
    with open(fpath, 'wb') as file:
        file.write(request.content)
    
    print(f'Image assessor downloaded to:\n {fpath}\n')
    
    """ Store info in a dictionary: """
    pathsDict, keys\
        = create_pathsDict_for_im_asr(
            projID, subjLab, expLab, asrID, roicolMod, fname, pathsDict
            )
    
    pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
    [expLab]['assessors'][asrID]['resources'][roicolMod]['files'][fname]\
    .update({'asrDate' : asrDate,
             'asrTime' : asrTime,
             'roicolName' : roicolName,
             'label' : label,
             'asrID' : asrID,
             'roicolMod' : roicolMod,
             'roicolDir' : exportDir,
             'roicolFpath' : fpath})
    
    return pathsDict

def download_im_asr_210810(xnatCfg, xnatParams, genParams, xnatSession=None,
                    pathsDict=None):
    """
    Download an image assessor (ROI Collection) of interest for a particular 
    scan from XNAT.
    
    Parameters
    ----------
    xnatCfg : dict
        Dictionary of XNAT configuration settings.
    xnatParams : dict
        Dictionary of XNAT parameters.
    genParams : dict
        Dictionary of general parameters.
    xnatSession : requests.models.Response, optional
        If provided a new session request will be avoided.
    pathsDict : dict, optional (None by default)
        Dictionary containing paths of data downloaded. 
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    
    Notes
    -----
    Although the structure of the REST call below:
    
    assessor = session.get(f'{url}/data/projects/{projID}/'
                           f'subjects/{subjLab}/'
                           f'experiments/{expLab}/'
                           f'assessors/{AsrID}/resources/{AsrMod}/files')
    
    follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Also, Jason suggested adding the file name after 'files', which worked. 
    """
    
    if xnatSession == None:
        xnatSession = create_session(xnatCfg)
    
    url = xnatSession.url
    
    rootExportDir = genParams['rootExportDir']
    
    projID = xnatParams['projID']
    subjLab = xnatParams['subjLab']
    expLab = xnatParams['expLab']
    scanID = xnatParams['scanID']
    roicolMod = xnatParams['roicolMod']
    roicolName = xnatParams['roicolName']
    
    p2c = genParams['p2c']
    
    if p2c:
        print('Inputs to download_im_asr():')
        print(f'  projID = {projID}')
        print(f'  subjLab = {subjLab}')
        print(f'  expLab = {expLab}')
        print(f'  scanID = {scanID}')
        print(f'  roicolMod = {roicolMod}')
        print(f'  roicolName = {roicolName}\n')
    
    """ The 'collectionType' key in 'data_fields' is either 'AIM' or 'SEG'. """
    if roicolMod == 'RTSTRUCT':
        collType = 'AIM'
    
    """ Get the experiment of interest: """
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}?format=json'
                              
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    experiment = request.json()
    
    """ There are likely more than one item in the list 
    experiment['items'][0]['children'] (e.g. one for scans and one for 
    image assessors). Get the index that corresponds to assessors, asr_ind, and
    the index that corresponds to scans, scan_ind: """
    asrInd = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'assessor' in experiment['items'][0]['children'][i]['field']:
            asrInd = i
    
    if asrInd == None:
        msg = f"There are no assessors for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    scanInd = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'scan' in experiment['items'][0]['children'][i]['field']:
            scanInd = i
    
    if scanInd == None:
        msg = f"There are no scans for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    if p2c:
        print(f'asrInd = {asrInd}')
        print(f'scanInd = {scanInd}\n')
    
    """ Get the SeriesUID for the scan of interest: """
    seriesUID = '' # initial value
    
    if pathsDict != None:
        try:
            seriesUID\
                = pathsDict['projects'][projID]['subjects'][subjLab]\
                    ['experiments'][expLab]['scans'][scanID]['seriesUID']
            
            excRaised = False # exception raised
            
        except Exception:
            excRaised = True
    
    if pathsDict == None or excRaised:
        for scan in experiment['items'][0]['children'][scanInd]['items']:
            if p2c:
                print(f"ID = {scan['data_fields']['ID']}")
                print(f"UID = {scan['data_fields']['UID']}\n")
            
            if scan['data_fields']['ID'] == scanID:
                seriesUID = scan['data_fields']['UID']
                
                if p2c:
                    print(f'seriesUID = {seriesUID}\n')
    
    if seriesUID == '':
        msg = f"No image assessors were found for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'."
        raise Exception(msg)
    
    """ Get the assessor ID and label of interest: """
    roicolInd = None # initial value
    asrID = None # initial value
    label = None # initial value
    asrDate = None # initial value
    asrTime = None # initial value
    
    for i in range(len(experiment['items'][0]['children'][asrInd]['items'])):
        if p2c:
            print(f"Item no {i}:")
        
        children = experiment['items'][0]['children'][asrInd]['items'][i]\
            ['children']
        
        if p2c:
            print(f"  children has {len(children)} items")
        
        fields = [child['field'] for child in children]
        
        # Get the index corresponding to the item containing the fields
        # 'out/file' and 'references/seriesUID':
        if 'out/file' in fields:
            fileInd = fields.index('out/file')
        else:
            fileInd = None
        
        if 'references/seriesUID' in fields:
            refInd = fields.index('references/seriesUID')
        else:
            refInd = None
        
        if p2c:
            print(f"  fields = {fields}")
            print(f"  fileInd = {fileInd}")
            print(f"  refInd = {refInd}\n")
        
        """ The following check is probably not required in most cases, but 
        covers cases where a ROI collection has been deleted and although the
        item with field 'out/file' is removed from the JSON, the item with
        field 'references/seriesUID' can persist. """
        # Check that items exist with both 'out/file' and 
        # 'references/seriesUID' fields:
        #if 'out/file' in fields and 'references/seriesUID' in fields:
        if isinstance(fileInd, int) and isinstance(refInd, int):
            if p2c:
                print("  This item's children contains items with fields",
                      "'out/file' and 'references/seriesUID'.\n")
            
            uid = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['children'][refInd]['items'][0]['data_fields']['seriesUID']
            
            type = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['data_fields']['collectionType']
            
            name = experiment['items'][0]['children'][asrInd]['items'][i]\
                ['data_fields']['name']
            
            if p2c:
                print(f"  uid = {uid}")
                print(f"  type = {type}")
                print(f"  name = {name}\n")
            
                #if uid == seriesUID:
                #    print(f"   * Matched collType = {type}")
                #    print(f"   * Matched name = {name}\n")
            
            if (uid == seriesUID and type == collType and name == roicolName):
                if p2c:
                    print(f"   * Matched seriesUID = {seriesUID}")
                    print(f"   * Matched collType = {collType}")
                    print(f"   * Matched name = {name}\n")
                
                roicolInd = i
                
                asrID = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['id']
                
                label = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['label']
                
                asrDate = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['date']
                
                asrTime = experiment['items'][0]['children'][asrInd]['items'][i]\
                    ['data_fields']['time']
        else:
            if p2c:
                print("  This item's children DOES NOT contain items with",
                      "both field 'out/file' and 'references/seriesUID'.\n")
    
    if p2c:
        print(f"Index of the ROI Collection of interest = {roicolInd}")
        print(f"ID of the ROI Collection of interest = {asrID}")
        print(f"Label of the ROI Collection of interest = {label}")
        print(f"Date and time of the ROI Collection of interest = {asrDate}",
              f"{asrTime}")
    
    if roicolInd == None:
        msg = f"No image assessors were found for:\n  projID = '{projID}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanID = '{scanID}'\ncontaining matches on:"\
              + f"\n  roicolMod = '{roicolMod}'\n  roicolName = '{roicolName}'."
        raise Exception(msg)
    
    """ Get the assessor of interest: """
    fname = label + '.dcm'
    
    uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/assessors/{asrID}/'\
          + f'resources/{roicolMod}/files/{fname}'
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    exportDir = os.path.join(rootExportDir, 'projects', projID, 'subjects',
                             subjLab, 'experiments', expLab, 'assessors',
                             asrID, 'resources', roicolMod, 'files')

    if not os.path.isdir(exportDir):
        Path(exportDir).mkdir(parents=True)
        print(f'Created directory:\n {exportDir}\n')
    
    fpath = os.path.join(exportDir, fname)
    
    with open(fpath, 'wb') as file:
        file.write(request.content)
    
    print(f'Image assessor downloaded to:\n {fpath}\n')
    
    """ Store info in a dictionary: """
    pathsDict, keys\
        = create_pathsDict_for_im_asr(projID, subjLab, expLab, asrID,
                                      roicolMod, fname, pathsDict)
    
    pathsDict['projects'][projID]['subjects'][subjLab]['experiments']\
    [expLab]['assessors'][asrID]['resources'][roicolMod]['files'][fname]\
    .update({'asrDate' : asrDate,
             'asrTime' : asrTime,
             'roicolName' : roicolName,
             'label' : label,
             'asrID' : asrID,
             'roicolMod' : roicolMod,
             'roicolDir' : exportDir,
             'roicolFpath' : fpath})
    
    return pathsDict

def download_im_asr_Pre_210810(url, projId, subjLab, expLab, scanId, mod, 
                    asrName, session=None, rootExportDir='', 
                    pathsDict=None, username=None, password=None, 
                    p2c=False):
    """
    Download an image assessor (ROI Collection) of interest for a particular 
    scan from XNAT.
    
    Parameters
    ----------
    url : str
        Address of XNAT (e.g. 'http://10.1.1.20').
    projId : str
        The project ID of interest.
    subjLab : str
        The subject label of interest. 
    expLab : str
        The DICOM study / XNAT experiment label of interest. 
    scanId : str
        The DICOM series label / XNAT scan ID of interest.
    mod : str
        The modality of the image assessor (ROI Collection) of interest.
    asrName : str
        The name of the image assessor (ROI Collection) of interest.
    session : requests session, optional
        If provided a new session request will be avoided.
    rootExportDir : str, optional
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, subjLab, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "xnat_downloads" within the default downloads directory.
    pathsDict : dict, optional
        Dictionary containing paths of data downloaded.
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = None)
        the user will be prompted to enter a user name.
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = None)
        the user will be prompted to enter a password.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    pathsDict : dict
        Dictionary containing paths of data downloaded.
    session : requests session
    
    Notes
    -----
    
    Although the structure of the REST call below:
    
    assessor = session.get(f'{url}/data/projects/{projId}/'
                           f'subjects/{subjLab}/'
                           f'experiments/{expLab}/'
                           f'assessors/{AsrId}/resources/{AsrMod}/files')
    
    follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Also, Jason suggested adding the file name after 'files', which worked. 
    """
    
    #p2c = True
    
    if p2c:
        print('Inputs to download_im_asr():')
        print(f'  projId = {projId}')
        print(f'  subjLab = {subjLab}')
        print(f'  expLab = {expLab}')
        print(f'  scanId = {scanId}')
        print(f'  mod = {mod}')
        print(f'  asrName = {asrName}\n')
    
    """ The 'collectionType' key in 'data_fields' is either 'AIM' or 'SEG'. """
    if mod == 'RTSTRUCT':
        collType = 'AIM'
    
    if session == None:
        session = create_session(url, username, password)
    
    if rootExportDir == '':
        rootExportDir = os.path.join(Path.home(), "Downloads", 
                                     "xnat_downloads")
    
    """ Get the experiment of interest: """
    uri = f'{url}/data/projects/{projId}/subjects/{subjLab}/'\
          + f'experiments/{expLab}?format=json'
                              
    request = session.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    experiment = request.json()
    
    """ There are likely more than one item in the list 
    experiment['items'][0]['children'] (e.g. one for scans and one for 
    image assessors). Get the index that corresponds to assessors, asr_ind, and
    the index that corresponds to scans, scan_ind: """
    asr_ind = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'assessor' in experiment['items'][0]['children'][i]['field']:
            asr_ind = i
    
    if asr_ind == None:
        msg = f"There are no assessors for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
        raise Exception(msg)
    
    scan_ind = None # initial value
    
    for i in range(len(experiment['items'][0]['children'])):
        if 'scan' in experiment['items'][0]['children'][i]['field']:
            scan_ind = i
    
    if scan_ind == None:
        msg = f"There are no scans for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
        raise Exception(msg)
    
    if p2c:
        print(f'asr_ind = {asr_ind}')
        print(f'scan_ind = {scan_ind}\n')
    
    """ Get the seriesUid for the scan of interest: """
    seriesUid = '' # initial value
    
    if pathsDict != None:
        try:
            seriesUid = pathsDict['projects'][projId]['subjects'][subjLab]\
                        ['experiments'][expLab]['scans'][scanId]['seriesUid']
            
            exception_raised = False
            
        except Exception:
            exception_raised = True
    
    if pathsDict == None or exception_raised:
        for scan in experiment['items'][0]['children'][scan_ind]['items']:
            if p2c:
                print(f"ID = {scan['data_fields']['ID']}")
                print(f"UID = {scan['data_fields']['UID']}\n")
            
            if scan['data_fields']['ID'] == scanId:
                seriesUid = scan['data_fields']['UID']
                
                if p2c:
                    print(f'seriesUid = {seriesUid}\n')
    
    
    if seriesUid == '':
        msg = f"No image assessors were found for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'."
        raise Exception(msg)
    
    
    """ Get the assessor ID and label of interest: """
    roicolInd = None # initial value
    asrId = None # initial value
    label = None # initial value
    date = None # initial value
    time = None # initial value
    
    for i in range(len(experiment['items'][0]['children'][asr_ind]['items'])):
        if p2c:
            print(f"Item no {i}:")
        
        children = experiment['items'][0]['children'][asr_ind]['items'][i]\
                   ['children']
        
        if p2c:
            print(f"  children has {len(children)} items")
        
        fields = [child['field'] for child in children]
        
        # Get the index corresponding to the item containing the fields
        # 'out/file' and 'references/seriesUid':
        if 'out/file' in fields:
            fileInd = fields.index('out/file')
        else:
            fileInd = None
        
        if 'references/seriesUid' in fields:
            refInd = fields.index('references/seriesUid')
        else:
            refInd = None
        
        if p2c:
            print(f"  fields = {fields}\n")
        
        """ The following check is probably not required in most cases, but 
        covers cases where a ROI collection has been deleted and although the
        item with field 'out/file' is removed from the JSON, the item with
        field 'references/seriesUid' can persist. """
        # Check that items exist with both 'out/file' and 
        # 'references/seriesUid' fields:
        #if 'out/file' in fields and 'references/seriesUid' in fields:
        if isinstance(fileInd, int) and isinstance(refInd, int):
            if p2c:
                print("  This item's children contains an items with fields",
                      "'out/file' and 'references/seriesUid'.\n")
            
            seriesUid = experiment['items'][0]['children'][asr_ind]['items'][i]\
                        ['children'][refInd]['items'][0]['data_fields']\
                        ['seriesUid']
        
            collType = experiment['items'][0]['children'][asr_ind]['items'][i]\
                       ['data_fields']['collectionType']
            
            name = experiment['items'][0]['children'][asr_ind]['items'][i]\
                   ['data_fields']['name']
            
            if p2c:
                print(f"seriesUid = {seriesUid}")
                print(f"collType = {collType}")
                print(f"name = {name}\n")
            
                #if seriesUid == seriesUid:
                #    print(f"   * Matched collectionType = {collectionType}")
                #    print(f"   * Matched name = {name}\n")
            
            if seriesUid == seriesUid and collType == collType and name == asrName:
                if p2c:
                    print(f"   * Matched seriesUid = {seriesUid}")
                    print(f"   * Matched collType = {collType}")
                    print(f"   * Matched name = {name}\n")
                
                roicolInd = i
                
                asrId = experiment['items'][0]['children'][asr_ind]['items'][i]\
                         ['data_fields']['id']
                
                label = experiment['items'][0]['children'][asr_ind]['items'][i]\
                        ['data_fields']['label']
                
                date = experiment['items'][0]['children'][asr_ind]['items'][i]\
                       ['data_fields']['date']
                
                time = experiment['items'][0]['children'][asr_ind]['items'][i]\
                       ['data_fields']['time']
        else:
            if p2c:
                print("  This item's children DOES NOT contain an items with",
                      "both field 'out/file' and 'references/seriesUid'.\n")
    
    
    if p2c:
        print(f"Index of the ROI Collection of interest = {roicolInd}")
        print(f"ID of the ROI Collection of interest = {asrId}")
        print(f"Label of the ROI Collection of interest = {label}")
        print(f"Date and time of the ROI Collection of interest = {date} {time}")
    
    if roicolInd == None:
        msg = f"No image assessors were found for:\n  projId = '{projId}'"\
              + f"\n  subjLab = '{subjLab}'\n  expLab = '{expLab}'"\
              + f"\n  scanId = '{scanId}'\ncontaining matches on:"\
              + f"\n  mod = '{mod}'\n  asrName = '{asrName}'."
        raise Exception(msg)
    
    
    """ Get the assessor of interest: """
    fname = label + '.dcm'
    
    uri = f'{url}/data/projects/{projId}/subjects/{subjLab}/'\
          + f'experiments/{expLab}/assessors/{asrId}/'\
          + f'resources/{mod}/files/{fname}'
    
    request = session.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    export_dir = os.path.join(rootExportDir, 'projects', projId, 'subjects',
                              subjLab, 'experiments', expLab, 'assessors',
                              asrId, 'resources', mod, 'files')

    if not os.path.isdir(export_dir):
        Path(export_dir).mkdir(parents=True)
        print(f'Created directory:\n {export_dir}\n')
    
    fpath = os.path.join(export_dir, fname)
    
    with open(fpath, 'wb') as file:
        file.write(request.content)
    
    print(f'Image assessor downloaded to:\n {fpath}\n')
    
    
    """ Store info in a dictionary: """
    pathsDict, keys\
        = create_paths_dict_for_im_asr(projId, subjLab, expLab, asrId, mod,
                                       fname, pathsDict)
    
    pathsDict['projects'][projId]['subjects'][subjLab]['experiments']\
    [expLab]['assessors'][asrId]['resources'][mod]['files'][fname]\
    .update({'date' : date,
             'time' : time,
             'asrName' : asrName,
             'label' : label,
             'asrId' : asrId,
             'mod' : mod,
             'dir' : export_dir,
             'fpath' : fpath})
    
    return pathsDict, session

def upload_roicol(
        roicol_fpath, url, proj_id, session_id, coll_label='',
        session=None, username=None, password=None
    ):
    """
    Upload a ROI Collection to a particular image session in XNAT.
    
    Parameters
    ----------
    roicol_fpath : str
        Full filepath of the ROI Collection to be uploaded to XNAT.
    url : str
        Address of XNAT (e.g. 'http://10.1.1.20').
    proj_id : str
        The project ID of the project destination.
    session_id : str
        The session ID label of the image session destination. 
    coll_label : str, optional
        String to assign to the collection label. If not provided the file name
        from roicol_fpath will be used. The default value is ''.
    session : requests session, optional
        If provided a new session request will be avoided.
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = None)
        the user will be prompted to enter a user name.
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = None)
        the user will be prompted to enter a password.
    
    Returns
    -------
    session : requests session
    """
    
    overwrite = 'false'
    #overwrite = 'true'
    
    if session == None:
        session = create_session(url, username, password)
    
    if coll_label == '':
        fname_with_ext = os.path.split(roicol_fpath)[1]

        fname, ext = os.path.splitext(
            os.path.split(fname_with_ext)[1]
            )

        #fname += f'_{coll_label}'
        #fname_with_ext = fname + ext
        
        coll_label = str(fname)
    
    #print(f'coll_label = {coll_label}')
    
    # Upload the ROI Collection:
    with open(roicol_fpath, 'rb') as file:
        ds = dcmread(file)
        
        mod = ds.Modality
        
        buf = BytesIO(file.read())
        
        uri = f"{url}/xapi/roi/projects/{proj_id}/sessions/{session_id}/" +\
            f"collections/{coll_label}?overwrite={overwrite}&type={mod}"
        
        #print(f'\nuri = {uri}\n')
        
        request = session.put(uri, data=buf)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    print(f'The ROI Collection was uploaded to: \n{uri}\n')
    
    return session