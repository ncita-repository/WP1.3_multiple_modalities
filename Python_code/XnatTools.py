# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:10:47 2021

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
GENERAL XNAT FUNCTIONS
******************************************************************************
******************************************************************************
"""


"""
Note:

How to read DICOM from a 'requests' (in memory) rather than from a file:

e.g.

from pydicom import dcmread
from pydicom.filebase import DicomBytesIO

# Get a listing of all resource files for the subject:
request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                          f'subjects/{SubjLabel}/'
                          f'files')

Fnames = []

for result in request.json()['ResultSet']['Result']:
    Fnames.append(result['Name'])

# Get the contents of a resource file:
request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                          f'subjects/{SubjLabel}/'
                          f'files/{Fnames[0]}')

raw = DicomBytesIO(request.content)

dro = dcmread(raw)
"""


def CreateXnatSession(XnatUrl, Username=None, Password=None):
    """
    Create an XNAT session.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
        
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    XnatSession : requests session
    """
    
    from getpass import getpass
    import requests
    #from pathlib import Path
    
    if Username == None:
        Username = getpass("Enter user name: ")
    
    if Password == None:
        Password = getpass("Enter password:")
        
    
    with requests.Session() as XnatSession:
        XnatSession.auth = (f'{Username}', f'{Password}')
        
    
    print(f'Connection established to XNAT at {XnatUrl}\n')
    
    return XnatSession



    """
    scan_id_list = []
    
    scans_request = XnatSession.get(f"{self.domain}/data/experiments/{self.experiment_id}/scans")
    
    scan_list_json = scans_request.json()
    
    [scan_id_list.append(scan['ID']) for scan in scan_list_json['ResultSet']['Result']]
    
    scan_full_download = XnatSession.get(f"{self.domain}/data/experiments/{self.experiment_id}/scans/{scan}/"
                                         f"resources/DICOM/files?format=zip")
    
    with open("path/to/save/data/at/", 'wb') as file:
        file.write(scan_full_download.content)
    
    """




def CreatePathsDictForExp(ProjId, SubjLabel, ExpLabel, PathsDict=None):
    """ Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level. """
    
    
    if PathsDict == None:
        PathsDict = {}
        keys = []
    else:
        keys = list(PathsDict.keys())
    
    
    if not 'projects' in keys:
        PathsDict.update({'projects' : {}})
    
    keys = list(PathsDict['projects'].keys())
    
    if not ProjId in keys:
        PathsDict['projects'].update({ProjId : {}})
    
    keys = list(PathsDict['projects'][ProjId].keys())
    
    if not 'subjects' in keys:
        PathsDict['projects'][ProjId].update({'subjects' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'].keys())
    
    if not SubjLabel in keys:
        PathsDict['projects'][ProjId]['subjects'].update({SubjLabel : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel].keys())
    
    if not 'experiments' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
        .update({'experiments' : {}})
        
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'].keys())
    
    if not ExpLabel in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        .update({ExpLabel : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel].keys())
    
    return PathsDict, keys







def CreatePathsDictForScan(ProjId, SubjLabel, ExpLabel, ScanId, 
                           PathsDict=None):
    """ Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level. """
    
    
    PathsDict, keys = CreatePathsDictForExp(ProjId, SubjLabel, ExpLabel, 
                                            PathsDict)
    
    if not 'scans' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel].update({'scans' : {}})
        
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['scans'].keys())
    
    if not ScanId in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['scans'].update({ScanId : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['scans'][ScanId].keys())
    
    if not 'resources' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['scans'][ScanId].update({'resources' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['scans'][ScanId]['resources']\
                .keys())
    
    if not 'DICOM' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['scans'][ScanId]['resources'].update({'DICOM' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['scans'][ScanId]['resources']\
                ['DICOM'].keys())
    
    if not 'files' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['scans'][ScanId]['resources']['DICOM']\
        .update({'files' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['scans'][ScanId]['resources']\
                ['DICOM']['files'].keys())
    
    return PathsDict, keys







def DownloadScan(XnatUrl, ProjId, SubjLabel, ExpLabel, ScanId,
                 XnatSession=None, ExportRootDir='default', PathsDict=None,
                 Username=None, Password=None):
    """
    Download scan from XNAT and the director of the downloaded DICOMs in a 
    dictionary.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    ExpLabel : string
        The DICOM study / XNAT experiment label of interest. 
    
    ScanId : string
        The DICOM series label / XNAT scan ID of interest.
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    ExportRootDir : string (optional; 'default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
        
    
        
    Outputs:
    *******
    
    PathsDict : dictionary
        Dictionary containing paths of data downloaded.
        
    XnatSession : requests session
    
    
    Notes:
    *****
    
    GET - /data/projects/{project-id}/subjects/{subject-id | subject-label}/
           experiments/{experiment-id | experiment-label}/scans/{scan-id}/files
           
    GET - /data/experiments/{experiment-id}/scans/{scan-id}/files
    
    https://wiki.xnat.org/display/XAPI/Image+Scan+Resource+API
    """
    
    import os
    from pathlib import Path
    #import requests
    import zipfile
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatData")
    
    
    
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}/scans/{ScanId}/'\
          + f'resources/DICOM/files?format=zip'
    
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
        else:
            msg = str(request)
            
        raise Exception(msg)
    
    
    
    ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 
                             'subjects', SubjLabel, 'experiments')

    if not os.path.isdir(ExportDir):
        Path(ExportDir).mkdir(parents=True)
        print(f'Created directory:\n {ExportDir}\n')
    
    Fpath = os.path.join(ExportDir,
                         f'Experiment_{ExpLabel}__Scan_{ScanId}.zip')
    
    if not os.path.exists(Fpath):
    
        with open(Fpath, 'wb') as file:
            file.write(request.content)
            
        print(f'Zipped file downloaded to: \n{Fpath}\n')
    
    
        zip_file = zipfile.ZipFile(Fpath, mode='r')
    
        zip_file.extractall(ExportDir)
        
        print(f'Zipped file extracted to: \n{ExportDir}\n')
    else:
        print(f'Zipped file already downloaded to: \n{Fpath}\n')
    
    
    
    """ Store info in a dictionary: """
    PathsDict, keys = CreatePathsDictForScan(ProjId, SubjLabel, ExpLabel, 
                                             ScanId, PathsDict)
    
    if not 'Dir' in keys:
        """ Get the experiment of interest: """
        uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
              + f'experiments/{ExpLabel}?format=json'
        
        request = XnatSession.get(uri)
        
        if not '200' in str(request):
            if '403' in str(request):
                msg = f'Response 403: Forbidden'
            elif '404' in str(request):
                msg = f'Response 403: Not found'
            elif '409' in str(request):
                msg = f'Response 409: Conflict'
            elif '422' in str(request):
                msg = f'Response 422: Unprocessable entity'
                
            raise Exception(msg)
            
        Experiment = request.json()
        
        #print(Experiment['items'], '\n')
        #print(Experiment['items'][0], '\n')
        #print(Experiment['items'][0]['children'], '\n')
        #print(Experiment['items'][0]['children'][1], '\n')
        #print(Experiment['items'][0]['children'][1]['items'], '\n')
        #print(Experiment['items'][0]['children'][1]['items'][0], '\n')
        #print(Experiment['items'][0]['children'][1]['items'][0]['data_fields'], '\n')
        
        #return XnatSession, Experiment

        #StudyUid = Experiment['items'][0]['children'][1]['items'][0]['data_fields']['type']['UID']
        #StudyUid = Experiment['items'][0]['children'][1]['items'][0]['data_fields']['UID'] # 27/05/2021
        StudyUid = Experiment['items'][0]['data_fields']['UID'] # 27/05/2021
        
        """ There may be more than one item in the list 
        Experiment['items'][0]['children'] (e.g. depending on whether there are 
        image assessors or not). Get the index that corresponds to scans: """
        ScanInd = None # initial value
        
        for i in range(len(Experiment['items'][0]['children'])):
            if 'scan' in Experiment['items'][0]['children'][i]['field']:
                ScanInd = i
        
        if ScanInd == None:
            msg = f"There are no scans for:\n  ProjId = '{ProjId}'"\
              + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
              + f"\n  ScanId = '{ScanId}'."
            raise Exception(msg)
            
        
        """ Get the list of series labels (scan IDs), SeriesUIDs and series 
        descriptions for this series/scan: """
        ScanIds = []
        SeriesUids = []
        SeriesDescs = []
        
        #for experiment in Experiment['items'][0]['children'][1]['items']: # 27/05/21
        for experiment in Experiment['items'][0]['children'][ScanInd]['items']: # 27/05/21
            ScanIds.append(experiment['data_fields']['ID'])
            SeriesUids.append(experiment['data_fields']['UID'])
            SeriesDescs.append(experiment['data_fields']['type'])
            
        
        SeriesUid = SeriesUids[ScanIds.index(ScanId)]
        SeriesDesc = SeriesDescs[ScanIds.index(ScanId)]
        
        """ The folder name between the 'scans' folder and 'resources' folder
        has the form 'ScanId-SeriesDesc', but it appears that spaces and 
        full stops within SeriesDesc become underscores in the folder structure 
        of the extracted zip file. """
        DirName = f'{ScanId}-{SeriesDesc}'.replace(' ', '_').replace('.', '_')
        
        #ResourceDir = os.path.join(ExportRootDir, ExpLabel, 'scans', DirName,
        #                           'resources')
        
        ResourceDir = os.path.join(ExportRootDir, 'projects', ProjId, 
                                   'subjects', SubjLabel, 'experiments', 
                                   ExpLabel, 'scans', DirName, 'resources')
        
        DicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['scans'][ScanId]['resources']['DICOM']['files']\
        .update({'StudyUID' : StudyUid,
                 'SeriesDesc' : SeriesDesc,
                 'SeriesUID' : SeriesUid,
                 'Dir' : DicomDir})
    
    
    return PathsDict, XnatSession






def CreatePathsDictForImAsr(ProjId, SubjLabel, ExpLabel, Id, Mod, Fname,
                            PathsDict=None):
    """ Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level.
    """
    
    PathsDict, keys = CreatePathsDictForExp(ProjId, SubjLabel, ExpLabel, 
                                            PathsDict)
    
    if not 'assessors' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel].update({'assessors' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['assessors'].keys())
    
    if not Id in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['assessors'].update({Id : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['assessors'][Id].keys())
    
    if not 'resources' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['assessors'][Id].update({'resources' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['assessors'][Id]\
                ['resources'].keys())
    
    if not Mod in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['assessors'][Id]['resources'].update({Mod : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['assessors'][Id]\
                ['resources'][Mod].keys())
    
    if not 'files' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['assessors'][Id]['resources'][Mod].update({'files' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['experiments'][ExpLabel]['assessors'][Id]\
                ['resources'][Mod]['files'].keys())
    
    if not Fname in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
        [ExpLabel]['assessors'][Id]['resources'][Mod]['files']\
        .update({Fname : {}})
    
    return PathsDict, keys







def DownloadImAsrs_INCOMPLETE(XnatUrl, ProjId, SubjLabel, ExpLabel, 
                     XnatSession=None, ExportRootDir='default', AsZip=True,
                     PathsDict=None, Username=None, Password=None):
    """
    INCOMPLETE - NOTES: 
        1. DownloadImAsr() downloads the assessor that corresponds to a scan
        of interest based on a search through the SeriesUID of all assessors.
        2. Hence DownloadImAsr() is more useful than this one.
        3. Currently only works as expected if AsZip=True (not completed due to
        #1-2.)
    
    Download all image assessors for a particular scan from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    ExpLabel : string
        The DICOM study / XNAT experiment label of interest. 
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    ExportRootDir : string (optional; 'default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    
    AsZip : boolean (optional; True by default)
        If True, all scan assessors will be downloaded as a zip file and 
        extracted in the desired directory. If False, the files will be 
        downloaded one-by-one.
    
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    XnatSession : requests session
    
    """
    
    import os
    from pathlib import Path
    #import requests
    import zipfile
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatData")
    
    """ Get a listing of all resource files for the study: """
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}/assessors/ALL/files?file_format=DICOM'
              
    if AsZip:
        uri += f'&format=zip'
                                   
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    
    
    ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 'subjects', 
                             SubjLabel, 'experiments', ExpLabel)

    if not os.path.isdir(ExportDir):
        print('Creating directory', ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    if AsZip:
        Fpath = os.path.join(ExportDir, 'assessors.zip')
        
        with open(Fpath, 'wb') as file:
            file.write(request.content)
        
        print(f'Zipped file downloaded to: \n{Fpath}\n')
        
        
        zip_file = zipfile.ZipFile(Fpath, mode='r')
        
        ExtractDir = os.path.join(ExportDir, 'assessors')
        
        if not os.path.isdir(ExtractDir):
            print('Creating directory', ExtractDir)
            Path(ExtractDir).mkdir(parents=True)
        
        zip_file.extractall(ExtractDir)
        
        print(f'Zipped file extracted to: \n{ExportDir}\n.')
    #else:
    #    for result in request.json()['ResultSet']['Result']:
    #        fulluri = f"{XnatUrl}{result['URI']}"
    #        
    #        Id = result['cat_ID']
    #        
    #        Name = result['Name']
    #        
    #        ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 'subjects',
    #                                 SubjLabel, 'resources', Id, 'files',
    #                                 Name)
    #        
    #        SubjAsrFpath = os.path.join(ExportDir, Name)
    #        
    #        request = XnatSession.get(fulluri)
    #        
    #        with open(SubjAsrFpath, 'wb') as file:
    #            file.write(request.content)
    
    return XnatSession







def DownloadImAsr(XnatUrl, ProjId, SubjLabel, ExpLabel, ScanId, Mod, Name,
                  XnatSession=None, ExportRootDir='default', PathsDict=None, 
                  Username=None, Password=None):
    """
    Download an image assessor of interest for a particular scan from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    ExpLabel : string
        The DICOM study / XNAT experiment label of interest. 
    
    ScanId : string
        The DICOM series label / XNAT scan ID of interest.
    
    Mod : string
        The modality of the image assessor of interest.
        
    Name : string
        The name of the image assessor of interest.
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    ExportRootDir : string (optional; 'default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    PathsDict : dictionary
        Dictionary containing paths of data downloaded.
        
    XnatSession : requests session
    
    """
    
    import os
    from pathlib import Path
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatData")
    
    
    
    """ Get the experiment of interest: """
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}?format=json'
                              
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
        
    Experiment = request.json()
    
    
    """ There are likely more than one item in the list 
    Experiment['items'][0]['children'] (e.g. one for scans and one for 
    image assessors). Get the index that corresponds to assessors, AsrInd, and
    the index that corresponds to scans, ScanInd: """
    AsrInd = None # initial value
    
    for i in range(len(Experiment['items'][0]['children'])):
        if 'assessor' in Experiment['items'][0]['children'][i]['field']:
            AsrInd = i
    
    if AsrInd == None:
        msg = f"There are no assessors for:\n  ProjId = '{ProjId}'"\
          + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
          + f"\n  ScanId = '{ScanId}'."
        raise Exception(msg)
    
    ScanInd = None # initial value
    
    for i in range(len(Experiment['items'][0]['children'])):
        if 'scan' in Experiment['items'][0]['children'][i]['field']:
            ScanInd = i
    
    if ScanInd == None:
        msg = f"There are no scans for:\n  ProjId = '{ProjId}'"\
          + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
          + f"\n  ScanId = '{ScanId}'."
        raise Exception(msg)
    
    #""" Get the SeriesUID for the scan of interest: """
    #if PathsDict != None:
    #    try:
    #        SeriesUid = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
    #        ['experiments'][ExpLabel]['scans'][ScanId]['SeriesUID']
    #        
    #    except Exception:
    #        """ Get the list of scan IDs and SeriesUIDs for this scan: """
    #        ScanIds = []
    #        SeriesUids = []
    #        
    #        #for experiment in Experiment['items'][0]['children'][1]['items']: # 27/05/21
    #        for experiment in Experiment['items'][0]['children'][ScanInd]['items']: # 27/05/21
    #            ScanIds.append(experiment['data_fields']['ID'])
    #            SeriesUids.append(experiment['data_fields']['UID'])
    #            
    #        """ The SeriesUID and series description for the scan of interest 
    #        (ScanId): """
    #        SeriesUid = SeriesUids[ScanIds.index(ScanId)]
    
    """ Above modified on 27/05/21. """
    
    """ Get the SeriesUID for the scan of interest: """
    if PathsDict != None:
        try:
            SeriesUid = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
            ['experiments'][ExpLabel]['scans'][ScanId]['SeriesUID']
            
            ExceptionRaised = False
            
        except Exception:
            ExceptionRaised = True
    
    if PathsDict == None or ExceptionRaised:
        """ Get the list of scan IDs and SeriesUIDs for this scan: """
        ScanIds = []
        SeriesUids = []
        
        #for experiment in Experiment['items'][0]['children'][1]['items']: # 27/05/21
        for experiment in Experiment['items'][0]['children'][ScanInd]['items']: # 27/05/21
            ScanIds.append(experiment['data_fields']['ID'])
            SeriesUids.append(experiment['data_fields']['UID'])
            
        """ The SeriesUID and series description for the scan of interest 
        (ScanId): """
        SeriesUid = SeriesUids[ScanIds.index(ScanId)]

    
    """ Get the (referenced) SeriesUIDs of all assessors: """
    RefSeriesUids = []

    #for experiment in Experiment['items'][0]['children'][0]['items']: # 27/05/21
    for experiment in Experiment['items'][0]['children'][AsrInd]['items']: # 27/05/21
        RefSeriesUids.append(experiment['children'][1]['items'][0]\
                             ['data_fields']['seriesUID'])
    
    """ The indices of the RefSeriesUids that match SeriesUid: """
    inds = [i for i, uid in enumerate(RefSeriesUids) if uid == SeriesUid]
    
    
    if inds == []:
        #msg = f"There are no scan assessors for ScanId '{ScanId}' "\
        #      + f"in ExpLabel '{ExpLabel}' for SubjLabel '{SubjLabel}' "\
        #      + f"in ProjId '{ProjId}'."
        msg = f"There are no scan assessors for:\n  ProjId = '{ProjId}'"\
              + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
              + f"\n  ScanId = '{ScanId}'\n  ScanLabel = '{Name}'."
        raise Exception(msg)
    
    
    """ Get the names and other meta data for the assessors for the scan of
    interest: """
    Dates = []
    Times = []
    Names = []
    Labels = []
    Mods = []
    Ids = []
    
    for ind in inds:
        #experiment = Experiment['items'][0]['children'][0]['items'][ind] # 27/05/21
        experiment = Experiment['items'][0]['children'][AsrInd]['items'][ind] # 27/05/21
        
        Dates.append(experiment['data_fields']['date'])
        Times.append(experiment['data_fields']['time'])
        Names.append(experiment['data_fields']['name'])
        Labels.append(experiment['data_fields']['label'])
        Mods.append(experiment['data_fields']['collectionType'])
        Ids.append(experiment['data_fields']['id'])
    
    
    """ The indices of assessors that match Mod.: """
    inds = [i for i, mod in enumerate(Mods) if mod == Mod]
    
    """ Use inds to filter the lists. """
    Dates = [Dates[i] for i in inds]
    Times = [Times[i] for i in inds]
    Names = [Names[i] for i in inds]
    Labels = [Labels[i] for i in inds]
    Mods = [Mods[i] for i in inds]
    Ids = [Ids[i] for i in inds]
    
    """ The index of the assessor matching Name: """
    try:
        ind = Names.index(Name)
    
    except Exception:
        #msg = f"There are no assessors for Scan '{ScanId}' with name "\
        #      + f"'{Name}' for ScanId '{ScanId}' in ExpLabel "\
        #      + f"'{ExpLabel}' for SubjLabel '{SubjLabel}' in ProjId "\
        #      + f"'{ProjId}'."
        msg = f"There are no scan assessors for:\n  ProjId = '{ProjId}'"\
              + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
              + f"\n  ScanId = '{ScanId}'\n  ScanLabel = '{Name}'."
        raise Exception(msg)
    
    Date = Dates[ind]
    Time = Times[ind]
    Label = Labels[ind]
    #Mod = Mods[ind]
    Id = Ids[ind]
    
    
    """ Get the assessor of interest: """
    
    """ 
    Note:
        
    Although the structure of the REST call below follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Assessor = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                               f'subjects/{SubjLabel}/'
                               f'experiments/{ExpLabel}/'
                               f'assessors/{AsrId}/resources/{AsrMod}/files')
    
    Jason suggested adding the file name after 'files', which worked. 
    """
    Fname = Label + '.dcm'
    
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}/assessors/{Id}/'\
          + f'resources/{Mod}/files/{Fname}'
    
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    
    
    #ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 'subjects',
    #                         SubjLabel, 'experiments', ExpLabel, 'assessors')
    ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 'subjects',
                             SubjLabel, 'experiments', ExpLabel, 'assessors',
                             Id, 'resources', Mod, 'files')

    if not os.path.isdir(ExportDir):
        Path(ExportDir).mkdir(parents=True)
        print(f'Created directory:\n {ExportDir}\n')
    
    #Fname = AsrLabel + '.dcm'
    Fpath = os.path.join(ExportDir, Fname)
    
    with open(Fpath, 'wb') as file:
        file.write(request.content)
    
    print(f'Image assessor downloaded to:\n {Fpath}\n')
    
    
    """ Store info in a dictionary: """
    PathsDict, keys\
    = CreatePathsDictForImAsr(ProjId, SubjLabel, ExpLabel, Id, Mod, Fname, 
                                PathsDict)
    
    PathsDict['projects'][ProjId]['subjects'][SubjLabel]['experiments']\
    [ExpLabel]['assessors'][Id]['resources'][Mod]['files'][Fname]\
    .update({'Date' : Date,
             'Time' : Time,
             'Name' : Name,
             'Label' : Label,
             'Id' : Id,
             'Mod' : Mod,
             'Dir' : ExportDir,
             'Fpath' : Fpath})
    
    return PathsDict, XnatSession






def GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, ExpLabel, Mod, Name):

    """ Get the scan assessor filename and scan ID corresponding to a scan
    assessor name of interest with given modality. """
    
    # Get all Ids:
    Ids = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
               ['experiments'][ExpLabel]['assessors'].keys()
               )

    Fname = None
    ID = None

    for Id in Ids:
        fnames = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                      ['experiments'][ExpLabel]['assessors'][Id]\
                      ['resources'][Mod]['files'].keys()
                      )

        for fname in fnames:
            name = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                            ['experiments'][ExpLabel]['assessors'][Id]\
                            ['resources'][Mod]['files'][fname]['Name']
            
            if name == Name:
                Fname = fname
                ID = Id
        
    
    #if ScanAsrId == None:
    #    msg = f"There was no match on ScanAsrName = {ScanAsrName}."
    #    
    #    raise Exception(msg)
    
    return Fname, ID






def CreatePathsDictForSubjAsr(ProjId, SubjLabel, Id, Name, PathsDict=None):
    """ Populate a dictionary mimacking XNAT file hierarchy, and return the 
    dictionary and the keys at the highest level. """
    
    
    if PathsDict == None:
        PathsDict = {}
        keys = []
    else:
        keys = list(PathsDict.keys())
    
    
    if not 'projects' in keys:
        PathsDict.update({'projects' : {}})
    
    keys = list(PathsDict['projects'].keys())
    
    if not ProjId in keys:
        PathsDict['projects'].update({ProjId : {}})
    
    keys = list(PathsDict['projects'][ProjId].keys())
    
    if not 'subjects' in keys:
        PathsDict['projects'][ProjId].update({'subjects' : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'].keys())
    
    if not SubjLabel in keys:
        PathsDict['projects'][ProjId]['subjects'].update({SubjLabel : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel].keys())
    
    if not 'resources' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
        .update({'resources' : {}})
        
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['resources'].keys())
    
    if not Id in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['resources']\
        .update({Id : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['resources'][Id].keys())
    
    if not 'files' in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['resources']\
        [Id].update({'files' : {}})
        
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['resources'][Id]['files'].keys())
    
    if not Name in keys:
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['resources']\
        [Id]['files'].update({Name : {}})
    
    keys = list(PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                ['resources'][Id]['files'][Name].keys())
    
    return PathsDict, keys







def DownloadSubjAsrs(XnatUrl, ProjId, SubjLabel, XnatSession=None, 
                     ExportRootDir='default', PathsDict=None, 
                     Username=None, Password=None):
    """
    Download all assessors for a subject of interest from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    ExportRootDir : string (optional; 'default' by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    PathsDict : dictionary
        Dictionary containing paths of data downloaded.
    
    XnatSession : requests session
    
    """
    
    import os
    from pathlib import Path
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatData")
    
    
    
    """ Get a listing of all resource files for the subject: """
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/files'
                              
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    for result in request.json()['ResultSet']['Result']:
        fulluri = f"{XnatUrl}{result['URI']}"
        
        Id = result['cat_ID']
        
        Name = result['Name']
        
        ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 'subjects',
                                 SubjLabel, 'resources', Id, 'files',
                                 Name)
        
        SubjAsrFpath = os.path.join(ExportDir, Name)
        
        request = XnatSession.get(fulluri)
        
        with open(SubjAsrFpath, 'wb') as file:
            file.write(request.content)
        
        """ Store info in a dictionary: """
        PathsDict, keys = CreatePathsDictForSubjAsr(ProjId, SubjLabel, 
                                                    Id, Name, 
                                                    PathsDict)
        
        PathsDict['projects'][ProjId]['subjects'][SubjLabel]['resources']\
        [Id]['files'][Name].update({'Name' : Name,
                                                  'Id' : Id,
                                                  'SubjAsrDir' : ExportDir,
                                                  'SubjAsrFpath' : SubjAsrFpath})
    
    
    return PathsDict, XnatSession







def UploadSubjAsr(SubjAsrFpath, XnatUrl, ProjId, SubjLabel, 
                  ContentLabel='SRO-DRO', XnatSession=None, 
                  Username=None, Password=None):
    """
    Upload a subject assessor to a particular subject in XNAT.
    
    Inputs:
    ******
    
    Fpath : string
        Full filepath of assessor to be uploaded to XNAT.
    
    XnatUrl : string
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    ContentLabel : string (optional; 'SRO-DRO' by default)
        String to assign to content label. Suggested values are:
            - 'SRO-DRO' for Spatial Registration Object type of DICOM 
            Registration Object
            - 'DSRO-DRO' for Deformable Spatial Registration Object type of 
            DICOM Registration Object
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
        
    XnatSession : requests session
    
    """
    
    import os
    #from pathlib import Path
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    SubjAsrFname = os.path.split(SubjAsrFpath)[1]
    
    
    """ Upload the resource: """
    with open(SubjAsrFpath, 'rb') as file:
        uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
              + f'resources/DICOM/files/{SubjAsrFname}'
        
        request = XnatSession.put(f'{uri}?inbody=true&file_format=DICOM'
                                  '&overwrite=true', data=file)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
        
    print(f'The DRO was uploaded as a subject assessor to: \n{uri}\n')
    
    return XnatSession







def GetFORuid(XnatUrl, ProjId, SubjLabel, ExpLabel, ScanId,
              XnatSession=None, Username=None, Password=None):
    """
    Get the FrameOfReferenceUID of a single DICOM within a specified scan.
    
    XnatUrl : string
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    ExpLabel : string
        The DICOM study / XNAT experiment label of interest. 
    
    ScanId : string
        The DICOM series label / XNAT scan ID of interest.
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    PathsDict : dictionary (optional; None by default)
        Dictionary containing paths of data downloaded.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    FORuid : string
        The FrameOfReferenceUID of a DICOM within the DICOM series / XNAT scan
        of interest.
    """
    
    from pydicom.filebase import DicomBytesIO
    from pydicom import dcmread
    
    # Get the Source scan files:
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}/scans/{ScanId}/files'
    
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    files = request.json()
    
    # Get the URI of the first file:
    uri = files['ResultSet']['Result'][0]['URI']
    
    fulluri = f"{XnatUrl}{uri}"
    
    # Get the file:
    request = XnatSession.get(fulluri)
    
    raw = DicomBytesIO(request.content)
            
    dcm = dcmread(raw)
    
    FORuid = dcm.FrameOfReferenceUID
    
    return FORuid






def SearchXnatForDro(XnatUrl, ProjId, SubjLabel, 
                     SrcExpLabel, SrcScanId, TrgExpLabel, TrgScanId,
                     Transform,
                     XnatSession=None, Username=None, Password=None):
    """
    Search XNAT for a DRO (as a subject assessor) whose FrameOfReferenceUIDs   
    match the Source (moving) and Target (fixed) FrameOfReferenceUIDs, and 
    whose SOP Class UID matches the transformation type specified by Transform. 
    If one (or more) matches, obtain (from the newest DRO for multiple matches)
    the FrameOfReferenceTransformation (for non-deformable transforms), or the 
    GridDimensions, GridResolution and VectorGridData (for deformable 
    transformation) of the moving (Source) sequence.
    
    Inputs:
    ******
    
    SrcDcmDir : string
        Directory containing the Source (moving) DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target (fixed) DICOMs.
    
    XnatUrl : string
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    SrcExpLabel : string
        The Source DICOM study / XNAT experiment label of interest. 
    
    SrcScanId : string
        The Source DICOM series label / XNAT scan ID of interest.
    
    TrgExpLabel : string
        The Target DICOM study / XNAT experiment label of interest. 
    
    TrgScanId : string
        The Target DICOM series label / XNAT scan ID of interest.
    
    Transform : string (optional; 'affine' by default)
        Denotes type of registration transformation. Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    TxMatrix : list of float strings or None
        List of float strings representing the non-deformable transformation 
        that transforms the moving (Source) image to the fixed (Target) image 
        (if there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    
    GridDims : list of integers or None
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image (if 
        there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    
    GridRes : list of floats or None
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image (if there 
        exists a DRO (as a subject assessor) that matches the requirements); 
        None if not.
    
    VectGridData : list of float strings or None
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image (if there exists a
        DRO (as a subject assessor) that matches the requirements); None if not.
    """
    
    from pydicom.filebase import DicomBytesIO
    from pydicom import dcmread
    import time
    
    #print(XnatSession)
    
    SrcFORuid = GetFORuid(XnatUrl, ProjId, SubjLabel, SrcExpLabel, SrcScanId,
                          XnatSession)
    
    TrgFORuid = GetFORuid(XnatUrl, ProjId, SubjLabel, TrgExpLabel, TrgScanId,
                          XnatSession)
    
    # Get a listing of all resource files for the subject:
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/files'
    
    request = XnatSession.get(uri)
    
    if not '200' in str(request):
        if '403' in str(request):
            msg = f'Response 403: Forbidden'
        elif '404' in str(request):
            msg = f'Response 403: Not found'
        elif '409' in str(request):
            msg = f'Response 409: Conflict'
        elif '422' in str(request):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    
    Fnames = []

    for result in request.json()['ResultSet']['Result']:
        #Fnames.append(result['Name']) # 07/05/21
        fname = result['Name'] # 07/05/21
        if '.dcm' in fname: # 07/05/21
            Fnames.append(fname) # 07/05/21
    
    # Get a list of the FrameOfReferenceUID as a tuple of pairs, 
    # e.g. (SrcForUid, TrgForUid), for all items:
    #ForUidPairs = []
    
    # Initialise a list of all FrameOfReferenceTransformations (from Spatial
    # Registration Objects), GridDimensions, GridResolutions and VectorGridData
    # (for Deformable Spatial Registration Objects):
    TxMatrices = []
    #TxMatrixTypes = []
    GridDimss = []
    GridRess = []
    VectGridDatas = []
    
    
    #ContentDateTimes = []
    
    # The indices of any DROs whose FORuids match the Source and Target FORuids: 
    Inds = []
    
    # Find the time difference between now and each ContentDateTime:
    Now = time.time()

    TimeDiffs = []
    
    for i in range(len(Fnames)):
        uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
              + f'files/{Fnames[i]}'
        
        request = XnatSession.get(uri)
        
        raw = DicomBytesIO(request.content)
        
        dro = dcmread(raw)
        
        Mod = f'{dro.Modality}'
        
        if Mod == 'REG': # 07/05/21
            DroType = f'{dro.SOPClassUID}'
            
            if 'Deformable' in DroType:
                if Transform == 'bspline':
                    MatchOnDroType = True
                else:
                    MatchOnDroType = False
            else:
                if Transform == 'bspline':
                    MatchOnDroType = False
                else:
                    MatchOnDroType = True
            
            #ForUidPairs.append((dro.RegistrationSequence[0].FrameOfReferenceUID,
            #                    dro.RegistrationSequence[1].FrameOfReferenceUID))
                       
            # Use f-strings instead of deepcopy:
            if 'Deformable' in DroType:
                FORuid0 = f'{dro.DeformableRegistrationSequence[0].SourceFrameOfReferenceUID}'
                FORuid1 = f'{dro.DeformableRegistrationSequence[1].SourceFrameOfReferenceUID}'
            else:
                FORuid0 = f'{dro.RegistrationSequence[0].FrameOfReferenceUID}'
                FORuid1 = f'{dro.RegistrationSequence[1].FrameOfReferenceUID}'
            
            #Type = f'{dro.RegistrationSequence[1].MatrixRegistrationSequence[0].MatrixSequence[0].FrameOfReferenceTransformationMatrixType}'
            #Type = Type.lower()
            
            #Mod = f'{dro.Modality}' # 07/05/21
                        
            #print(Type)
            
            #if FORuid0 == TrgFORuid and FORuid1 == SrcFORuid and Type == Transform\
            #and Mod == 'REG': # 07/05/21
            if FORuid0 == TrgFORuid and FORuid1 == SrcFORuid and MatchOnDroType:
                Inds.append(i)
                
                if 'Deformable' in DroType:
                    GridDimss.append(dro.DeformableRegistrationSequence[1]\
                                        .DeformableRegistrationGridSequence[0]\
                                        .GridDimensions)
                    
                    GridRess.append(dro.DeformableRegistrationSequence[1]\
                                       .DeformableRegistrationGridSequence[0]\
                                       .GridResolution)
                    
                    VectGridDatas.append(dro.DeformableRegistrationSequence[1]\
                                            .DeformableRegistrationGridSequence[0]\
                                            .VectorGridData)
                else:
                    TxMatrices.append(dro.RegistrationSequence[1]\
                                         .MatrixRegistrationSequence[0]\
                                         .MatrixSequence[0]\
                                         .FrameOfReferenceTransformationMatrix)
            
                #TxMatrixTypes.append(dro.RegistrationSequence[1]\
                #                        .MatrixRegistrationSequence[0]\
                #                        .MatrixSequence[0]\
                #                        .FrameOfReferenceTransformationMatrixType)
    
                ContentDateTime = f'{dro.ContentDate}{dro.ContentTime}'
    
                #ContentDateTimes.append(ContentDateTime)
    
                TimeDiffs.append(Now - time.mktime(time.strptime(ContentDateTime, 
                                                                 '%Y%m%d%H%M%S')))

    
    if TimeDiffs:
        # The index of the most recent ContentDateTime:
        ind = TimeDiffs.index(min(TimeDiffs))
        
        if 'Deformable' in DroType:
            TxMatrix = None
            GridDims = GridDimss[ind]
            GridRes = GridRess[ind]
            VectGridData = VectGridDatas[ind]
            
            msg = f"There were {len(Fnames)} subject assessors found, of which "\
                  + f"{len(TimeDiffs)} contained matching FrameOfReferenceUIDs "\
                  + "and transform type, the most \nrecent of which contains the "\
                  + f"grid dimensions (GridDims):\n{GridDims}\ngrid resolution:"\
                  + f"{GridRes}\nand vector grid data containing {len(VectGridData)}"\
                  + " elements"
        else:
            TxMatrix = TxMatrices[ind]
            GridDims = None
            GridRes = None
            VectGridData = None
        
            msg = f"There were {len(Fnames)} subject assessors found, of which "\
                  + f"{len(TimeDiffs)} contained matching FrameOfReferenceUIDs "\
                  + "and transform type, the most \nrecent of which contains the "\
                  + f"transformation matrix (TxMatrix):\n{TxMatrix}\n"
    else:
        msg = "There were no DROs with FrameOfReferenceUIDs matching those of"\
              + f" the Source ({SrcFORuid}) and Target ({TrgFORuid}) image "\
              + "series, and whose SOP Class UID matched the desired "\
              + f"transform type ({Transform}).\n"
        
        TxMatrix = None
        GridDims = None
        GridRes = None
        VectGridData = None
    
    print(msg)
    
    return TxMatrix, GridDims, GridRes, VectGridData