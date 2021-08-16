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







def CreateXnatSession(XnatUrl, Username=None, Password=None):
    """Create a requests session.
    
    Parameters
    ----------
    
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
        
    Username : string (optional; None by default)
        The username for XNAT log-in.  If not provided (i.e. Username = None)
        the user will be prompted to enter a user name.
    
    Password : string (optional; None by default)
        The password for XNAT log-in.  If not provided (i.e. Password = None)
        the user will be prompted to enter a password.
    
    Returns
    -------
    
    XnatSession : requests session
    """
    
    import time
    from getpass import getpass
    import requests
    #from pathlib import Path
    
    if Username == None:
        Username = getpass("Enter user name: ")
    
    if Password == None:
        Password = getpass("Enter password:")
        
    
    with requests.Session() as XnatSession:
        XnatSession.url = XnatUrl
        XnatSession.auth = (f'{Username}', f'{Password}')
    
    # Test connection:
    request = requests.get(XnatSession.url, auth=XnatSession.auth)
    
    #if request.status_code == 200:
    #    print(f'Connection established to XNAT at {XnatUrl}\n')
    #elif request.status_code == 401:
    #    print(f'Connection not established to XNAT at {XnatUrl}.',
    #          f'\nAuthentication error ({request.status_code}). Check your',
    #          'log-in details and try again.\n')
    #else:
    #    print(f'Connection not established to XNAT at {XnatUrl}. \nError code',
    #          f'{request.status_code}.\n')
    
    # Raise status error if not None:
    if request.raise_for_status() == None:
        print(f'Connection established to XNAT at {XnatUrl}\n')
    else:
        print(request.raise_for_status())
    
    return XnatSession



"""
01/07/2021

Considered writing a simple function to check whether an existing requests
session (e.g. XnatSession) has expired, and if so, create a new session.

The idea was to not create multiple log-in sessions, and simply use on if it
exists, i.e.:
    
if not 'XnatSession' in locals():    
    from XnatTools import CreateXnatSession
    XnatSession = CreateXnatSession(XnatUrl)
    
However it's possible for sessions to expire.  And during some testing I ran
into errors when making requests since an existing XnatSession had expired
(which wasn't obvious). 

The plan was to make a simple request, and check the headers for an expiration 
date-time:
    
import requests

request = requests.get(XnatSession.url, auth=XnatSession.auth)

request.headers will return, for example:

{'Server': 'nginx/1.14.0 (Ubuntu)', 
 'Date': 'Wed, 30 Jun 2021 03:16:32 GMT', 
 'Content-Type': 'text/html;charset=ISO-8859-1', 
 'Transfer-Encoding': 'chunked', 
 'Connection': 'keep-alive', 
 'Set-Cookie': 'JSESSIONID=4F6F9FB75FCD730761EA33F4B6BFBEB4; Path=/; HttpOnly, SESSION_EXPIRATION_TIME="1625022991874,432000000"; Version=1; Path=/', 
 'Expires': 'Tue, 03 Jul 2001 06:00:00 GMT', 
 'Last-Modified': 'Wed Jun 30 03:16:31 UTC 2021', 
 'Cache-Control': 'no-store, no-cache, must-revalidate, max-age=0, post-check=0, pre-check=0', 
 'Pragma': 'no-cache', 
 'X-Content-Type-Options': 'nosniff', 
 'X-XSS-Protection': '1; mode=block', 
 'X-Frame-Options': 'SAMEORIGIN', 
 'Content-Security-Policy': "frame-ancestors 'self'", 
 'Content-Language': 'en-US', 
 'Content-Encoding': 'gzip'}

Although the value for 'Expires' had a sensible day and month, the year was
2001!! 

I hoped that the values in 'Set-Cookie' following 'SESSION_EXPIRATION_TIME' 
might be useable:
    
session_exp_time = request.headers['Set-Cookie'].split('SESSION_EXPIRATION_TIME="')[1].split(',')[0]

datetime.datetime.fromtimestamp(int(session_exp_time)/1000) # convert from ms to s

> datetime.datetime(2021, 6, 30, 3, 14, 34, 266000)

which was in the past.  Since I was unable to identify a sensible expiration 
date-time from the headers I gave up.
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
    from zipfile import ZipFile
    
    """ Set to True if issue with zip file (e.g. duplication of files): """
    OverwriteZip = False
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatData")
    
    
    
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}/scans/{ScanId}/'\
          + 'resources/DICOM/files?format=zip'
    
    request = XnatSession.get(uri)
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #    else:
    #        msg = str(request)
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    
    ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 
                             'subjects', SubjLabel, 'experiments')

    if not os.path.isdir(ExportDir):
        Path(ExportDir).mkdir(parents=True)
        print(f'Created directory:\n {ExportDir}\n')
    
    Fpath = os.path.join(ExportDir,
                         f'Experiment_{ExpLabel}__Scan_{ScanId}.zip')
    
    if not os.path.exists(Fpath) or OverwriteZip:
        with open(Fpath, 'wb') as file:
            file.write(request.content)
            
        print(f'Zipped file downloaded to: \n{Fpath}\n')
        
        #zip_file = ZipFile(Fpath, mode='r')
        #zip_file.extractall(ExportDir)
        
        with ZipFile(Fpath, 'r') as file:
            file.extractall(ExportDir)
        
        print(f'Zipped file extracted to: \n{ExportDir}\n')
    else:
        if not OverwriteZip:
            print(f'Zipped file already downloaded to: \n{Fpath}\n')
    
    
    """ Get the directory name of the exported DICOMs: """
    # Get the filepath of the first file in the zip file:
    with ZipFile(Fpath, 'r') as file:
        fpath = file.infolist()[0].filename
    
    fullfpath = os.path.join(ExportDir, fpath)
    
    #print(f'fullfpath = {fullfpath}\n')
    
    DicomDir = os.path.split(fullfpath)[0]
    
    print(f"DICOMs downloaded to:\n{DicomDir}\n")
    
    """ Store info in a dictionary: """
    PathsDict, keys = CreatePathsDictForScan(ProjId, SubjLabel, ExpLabel, 
                                             ScanId, PathsDict)
    
    if not 'Dir' in keys:
        """ Get the experiment of interest: """
        uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
              + f'experiments/{ExpLabel}?format=json'
        
        request = XnatSession.get(uri)
        
        #if not '200' in str(request):
        #    if '403' in str(request):
        #        msg = 'Response 403: Forbidden'
        #    elif '404' in str(request):
        #        msg = 'Response 403: Not found'
        #    elif '409' in str(request):
        #        msg = 'Response 409: Conflict'
        #    elif '422' in str(request):
        #        msg = 'Response 422: Unprocessable entity'
        #        
        #    raise Exception(msg)
        
        # Raise status error if not None:
        if request.raise_for_status() != None:
            print(request.raise_for_status())
            
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
        
        """ Replacing code below with something more reliable (i.e. rather than
        replacing spaces, dashes and full stops with underscores) (01/06/21).
        """
            
        #""" The folder name between the 'scans' folder and 'resources' folder
        #has the form 'ScanId-SeriesDesc', but it appears that spaces, dashes and 
        #full stops within SeriesDesc become underscores in the folder structure 
        #of the extracted zip file. """
        ##print(f'SeriesUid = {SeriesUid}')
        ##print(f'SeriesDesc = {SeriesDesc}\n')
        #
        ##DirName = f'{ScanId}-{SeriesDesc}'.replace(' ', '_').replace('.', '_') # 01/06/21
        #DirName = f'{ScanId}-{SeriesDesc}'.replace(' ', '_').replace('.', '_').replace('-', '_') # 01/06/21
        #
        ##ResourceDir = os.path.join(ExportRootDir, ExpLabel, 'scans', DirName,
        ##                           'resources')
        #
        #ResourceDir = os.path.join(ExportRootDir, 'projects', ProjId, 
        #                           'subjects', SubjLabel, 'experiments', 
        #                           ExpLabel, 'scans', DirName, 'resources')
        #
        #DicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        #
        #print(f"DICOMs downloaded to:\n{DicomDir}\n")
        
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
        uri += '&format=zip'
                                   
    request = XnatSession.get(uri)
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    
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







def DownloadImAsr_OLD(XnatUrl, ProjId, SubjLabel, ExpLabel, ScanId, Mod, Name,
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
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
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
    
    """ Below commented section was modified on 27/05/21 - see below. """
    
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
    for assessor in Experiment['items'][0]['children'][AsrInd]['items']: # 27/05/21
        RefSeriesUids.append(assessor['children'][1]['items'][0]\
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
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    
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






def DownloadImAsr(XnatUrl, ProjId, SubjLabel, ExpLabel, ScanId, Mod, Name,
                  XnatSession=None, ExportRootDir='default', PathsDict=None, 
                  Username=None, Password=None, LogToConsole=False):
    """
    Download an image assessor (ROI Collection) of interest for a particular 
    scan from XNAT.
    
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
        The modality of the image assessor (ROI Collection) of interest.
        
    Name : string
        The name of the image assessor (ROI Collection) of interest.
    
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
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    PathsDict : dictionary
        Dictionary containing paths of data downloaded.
        
    XnatSession : requests session
    
    
    Notes:
    ******
    
    Although the structure of the REST call below:
    
    Assessor = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                               f'subjects/{SubjLabel}/'
                               f'experiments/{ExpLabel}/'
                               f'assessors/{AsrId}/resources/{AsrMod}/files')
    
    follows the how-to guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Also, Jason suggested adding the file name after 'files', which worked. 
    """
    
    import os
    from pathlib import Path
    
    #LogToConsole = True
    
    if LogToConsole:
        print('Inputs to DownloadImAsr():')
        print(f'  ProjId = {ProjId}')
        print(f'  SubjLabel = {SubjLabel}')
        print(f'  ExpLabel = {ExpLabel}')
        print(f'  ScanId = {ScanId}')
        print(f'  Mod = {Mod}')
        print(f'  Name = {Name}\n')
    
    """ The 'collectionType' key in 'data_fields' is either 'AIM' or 'SEG'. """
    if Mod == 'RTSTRUCT':
        CollType = 'AIM'
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatData")
    
    
    
    """ Get the experiment of interest: """
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}?format=json'
                              
    request = XnatSession.get(uri)
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
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
    
    if LogToConsole:
        print(f'AsrInd = {AsrInd}')
        print(f'ScanInd = {ScanInd}\n')
    
    
    """ Get the SeriesUID for the scan of interest: """
    SeriesUid = '' # initial value
    
    if PathsDict != None:
        try:
            SeriesUid = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                        ['experiments'][ExpLabel]['scans'][ScanId]['SeriesUID']
            
            ExceptionRaised = False
            
        except Exception:
            ExceptionRaised = True
    
    if PathsDict == None or ExceptionRaised:
        for scan in Experiment['items'][0]['children'][ScanInd]['items']:
            if LogToConsole:
                print(f"ID = {scan['data_fields']['ID']}")
                print(f"UID = {scan['data_fields']['UID']}\n")
            
            if scan['data_fields']['ID'] == ScanId:
                SeriesUid = scan['data_fields']['UID']
                
                if LogToConsole:
                    print(f'SeriesUid = {SeriesUid}\n')
    
    
    if SeriesUid == '':
        msg = f"No image assessors were found for:\n  ProjId = '{ProjId}'"\
              + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
              + f"\n  ScanId = '{ScanId}'."
        raise Exception(msg)
    
    
    """ Get the assessor ID and label of interest: """
    RoiColInd = None # initial value
    Id = None # initial value
    Label = None # initial value
    Date = None # initial value
    Time = None # initial value
    
    for i in range(len(Experiment['items'][0]['children'][AsrInd]['items'])):
        if LogToConsole:
            print(f"Item no {i}:")
        
        children = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                   ['children']
        
        if LogToConsole:
            print(f"  children has {len(children)} items")
        
        fields = [child['field'] for child in children]
        
        # Get the index corresponding to the item containing the fields
        # 'out/file' and 'references/seriesUID':
        if 'out/file' in fields:
            FileInd = fields.index('out/file')
        else:
            FileInd = None
        
        if 'references/seriesUID' in fields:
            RefInd = fields.index('references/seriesUID')
        else:
            RefInd = None
        
        if LogToConsole:
            print(f"  fields = {fields}\n")
        
        """ The following check is probably not required in most cases, but 
        covers cases where a ROI collection has been deleted and although the
        item with field 'out/file' is removed from the JSON, the item with
        field 'references/seriesUID' can persist. """
        # Check that items exist with both 'out/file' and 
        # 'references/seriesUID' fields:
        #if 'out/file' in fields and 'references/seriesUID' in fields:
        if isinstance(FileInd, int) and isinstance(RefInd, int):
            if LogToConsole:
                print("  This item's children contains an items with fields",
                      "'out/file' and 'references/seriesUID'.\n")
            
            seriesUid = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                        ['children'][RefInd]['items'][0]['data_fields']\
                        ['seriesUID']
        
            collType = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                       ['data_fields']['collectionType']
            
            name = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                   ['data_fields']['name']
            
            if LogToConsole:
                print(f"  seriesUID = {seriesUid}")
                print(f"  collType = {collType}")
                print(f"  name = {name}\n")
            
                #if seriesUid == SeriesUid:
                #    print(f"   * Matched collectionType = {collectionType}")
                #    print(f"   * Matched name = {name}\n")
            
            if seriesUid == SeriesUid and collType == CollType and name == Name:
                if LogToConsole:
                    print(f"   * Matched seriesUID = {seriesUid}")
                    print(f"   * Matched collType = {collType}")
                    print(f"   * Matched name = {name}\n")
                
                RoiColInd = i
                
                Id = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                     ['data_fields']['id']
                
                Label = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                        ['data_fields']['label']
                
                Date = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                       ['data_fields']['date']
                
                Time = Experiment['items'][0]['children'][AsrInd]['items'][i]\
                       ['data_fields']['time']
        else:
            if LogToConsole:
                print("  This item's children DOES NOT contain an items with",
                      "both field 'out/file' and 'references/seriesUID'.\n")
    
    
    if LogToConsole:
        print(f"Index of the ROI Collection of interest = {RoiColInd}")
        print(f"ID of the ROI Collection of interest = {Id}")
        print(f"Label of the ROI Collection of interest = {Label}")
        print(f"Date and time of the ROI Collection of interest = {Date} {Time}")
    
    if RoiColInd == None:
        msg = f"No image assessors were found for:\n  ProjId = '{ProjId}'"\
              + f"\n  SubjLabel = '{SubjLabel}'\n  ExpLabel = '{ExpLabel}'"\
              + f"\n  ScanId = '{ScanId}'\ncontaining matches on:"\
              + f"\n  Mod = '{Mod}'\n  Name = '{Name}'."
        raise Exception(msg)
    
    
    """ Get the assessor of interest: """
    Fname = Label + '.dcm'
    
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
          + f'experiments/{ExpLabel}/assessors/{Id}/'\
          + f'resources/{Mod}/files/{Fname}'
    
    request = XnatSession.get(uri)
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    ExportDir = os.path.join(ExportRootDir, 'projects', ProjId, 'subjects',
                             SubjLabel, 'experiments', ExpLabel, 'assessors',
                             Id, 'resources', Mod, 'files')

    if not os.path.isdir(ExportDir):
        Path(ExportDir).mkdir(parents=True)
        print(f'Created directory:\n {ExportDir}\n')
    
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
    """ 
    06/07/21:
        See get_scan_asr_fname_and_id in xnat_tools.format_paths_dict.py
        (also xnat_tools.scan_assessors.py).
        
    Get the scan assessor filename and scan ID corresponding to a scan
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
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
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
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
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
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = f'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = f'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = f'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = f'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    files = request.json()
    
    # Get the URI of the first DICOM file:
    for file in files['ResultSet']['Result']:
        if file['collection'] == 'DICOM':
            uri = file['URI']
            break
    
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
                     XnatSession=None, Username=None, Password=None,
                     LogToConsole=False):
    """
    Search XNAT for a DRO (as a subject assessor) whose FrameOfReferenceUIDs   
    match the Source (moving) and Target (fixed) FrameOfReferenceUIDs, and 
    whose SOP Class UID matches the transformation type specified by Transform. 
    If one (or more) matches, return the (newest for multiple matches) DRO,
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
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom Object or None
        A DICOM Registration Object (if there exists a DRO, as a subject 
        assessor) that matches the requirements); None if not.
    
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
    #import time
    from datetime import datetime
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of SearchXnatForDro():')
        print('\n\n', '-'*120)
        
    #print(XnatSession)
    
    SrcFORuid = GetFORuid(XnatUrl, ProjId, SubjLabel, SrcExpLabel, SrcScanId,
                          XnatSession)
    
    TrgFORuid = GetFORuid(XnatUrl, ProjId, SubjLabel, TrgExpLabel, TrgScanId,
                          XnatSession)
    
    # Get a listing of all resource files for the subject:
    uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/files'
    
    request = XnatSession.get(uri)
    
    #if not '200' in str(request):
    #    if '403' in str(request):
    #        msg = 'Response 403: Forbidden'
    #    elif '404' in str(request):
    #        msg = 'Response 403: Not found'
    #    elif '409' in str(request):
    #        msg = 'Response 409: Conflict'
    #    elif '422' in str(request):
    #        msg = 'Response 422: Unprocessable entity'
    #        
    #    raise Exception(msg)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    
    Fnames = []

    for result in request.json()['ResultSet']['Result']:
        #Fnames.append(result['Name']) # 07/05/21
        fname = result['Name'] # 07/05/21
        if '.dcm' in fname: # 07/05/21
            Fnames.append(fname) # 07/05/21
    
    if LogToConsole:
        print(f'  There were {len(Fnames)} DICOM subject assessors found.')
                
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
    
    # The indices of any DROs whose FORuids match the Source and Target FORuids
    # and whose SOPClassUID is compatible with the Transform: 
    Inds = []
    
    # Find the time difference between now and each ContentDateTime:
    #Now = time.time()
    #Now = datetime.now()
    #Now = datetime.strftime(Now, '%Y%m%d%H%M%S.%f')
    
    #if LogToConsole:
    #    print(f'  Now = {Now}')

    #TimeDiffs = []
    ContentDateTimes = []
    Dros = []
    
    for i in range(len(Fnames)):
        uri = f'{XnatUrl}/data/projects/{ProjId}/subjects/{SubjLabel}/'\
              + f'files/{Fnames[i]}'
        
        request = XnatSession.get(uri)
        
        raw = DicomBytesIO(request.content)
        
        # Read in the file (which may or may not be a DRO):
        Dro = dcmread(raw)
        
        Dros.append(Dro)
        
        #print('Dro =', Dro, '\n')
        #print(f'Dro.SOPClassUID = {Dro.SOPClassUID}\n')
        #print(f'Dro[0x0008, 0x0016].name = {Dro[0x0008, 0x0016].name}\n')
        #print(f'Dro[0x0008, 0x0016].value = {Dro[0x0008, 0x0016].value}\n')
        
        Mod = f'{Dro.Modality}'
        
        # Only proceed if the DICOM file is a registration object:
        if Mod == 'REG': # 07/05/21
            #DroType = f'{Dro.SOPClassUID}' # this is a number (e.g. 1.2.840.10008.5.1.4.1.1.66.3)
            DroType = f'{Dro[0x0008, 0x0016].repval}' # the representation of the element's value
            
            if 'Deformable' in DroType:
                if Transform == 'bspline':
                    MatchOnDroType = True
                else:
                    MatchOnDroType = False
            else:
                if Transform in ['rigid', 'affine']:
                    MatchOnDroType = True
                else:
                    MatchOnDroType = False
            
            if LogToConsole:
                print(f'  File {i+1} of {len(Fnames)}')
                print(f'    DroType = {DroType}\n')
                print(f'    MatchOnDroType = {MatchOnDroType}\n')
            
            #ForUidPairs.append((Dro.RegistrationSequence[0].FrameOfReferenceUID,
            #                    Dro.RegistrationSequence[1].FrameOfReferenceUID))
                       
            # Use f-strings instead of deepcopy:
            if 'Deformable' in DroType:
                FORuid0 = f'{Dro.DeformableRegistrationSequence[0].SourceFrameOfReferenceUID}'
                FORuid1 = f'{Dro.DeformableRegistrationSequence[1].SourceFrameOfReferenceUID}'
            else:
                FORuid0 = f'{Dro.RegistrationSequence[0].FrameOfReferenceUID}'
                FORuid1 = f'{Dro.RegistrationSequence[1].FrameOfReferenceUID}'
            
            #Type = f'{Dro.RegistrationSequence[1].MatrixRegistrationSequence[0].MatrixSequence[0].FrameOfReferenceTransformationMatrixType}'
            #Type = Type.lower()
            
            #Mod = f'{Dro.Modality}' # 07/05/21
                        
            #print(Type)
            
            #if FORuid0 == TrgFORuid and FORuid1 == SrcFORuid and Type == Transform\
            #and Mod == 'REG': # 07/05/21
            if FORuid0 == TrgFORuid and FORuid1 == SrcFORuid and MatchOnDroType:
                Inds.append(i)
                
                if 'Deformable' in DroType:
                    GridDimss.append(Dro.DeformableRegistrationSequence[1]\
                                        .DeformableRegistrationGridSequence[0]\
                                        .GridDimensions)
                    
                    GridRess.append(Dro.DeformableRegistrationSequence[1]\
                                       .DeformableRegistrationGridSequence[0]\
                                       .GridResolution)
                    
                    VectGridDatas.append(Dro.DeformableRegistrationSequence[1]\
                                            .DeformableRegistrationGridSequence[0]\
                                            .VectorGridData)
                else:
                    TxMatrices.append(Dro.RegistrationSequence[1]\
                                         .MatrixRegistrationSequence[0]\
                                         .MatrixSequence[0]\
                                         .FrameOfReferenceTransformationMatrix)
            
                #TxMatrixTypes.append(Dro.RegistrationSequence[1]\
                #                        .MatrixRegistrationSequence[0]\
                #                        .MatrixSequence[0]\
                #                        .FrameOfReferenceTransformationMatrixType)
    
                ContentDateTime = f'{Dro.ContentDate}{Dro.ContentTime}'
                
                if LogToConsole:
                    print(f'  ContentDateTime = {ContentDateTime}')
                
                # Convert to datetime object:
                ContentDateTime = datetime.strptime(ContentDateTime, 
                                                    '%Y%m%d%H%M%S.%f')
                ContentDateTimes.append(ContentDateTime)
                
                #print(f'ContentDateTime = {ContentDateTime}')
                #print(f'Now = {Now}')
                
                #t0 = time.mktime(time.strptime(ContentDateTime, 
                #                               '%Y%m%d%H%M%S.%f'))
                #t0 = datetime.strptime(ContentDateTime, '%Y%m%d%H%M%S.%f')
                
                #if LogToConsole:
                #    print(f'  t0 = {t0}')
                
                #diff = Now - t0
                
                #print(f't0 = {t0}')
                #print(f'diff = {diff}\n')
    
                #TimeDiffs.append(diff)
    
    if LogToConsole:
        print(f'  The indices of the files that match = {Inds}')
        #print(f'  TimeDiffs = {TimeDiffs}')
        print(f'  ContentDateTimes = {ContentDateTimes}')
    
    if ContentDateTimes:#TimeDiffs:
        # The index of the most recent ContentDateTime:
        #NewInd = TimeDiffs.index(min(TimeDiffs))
        NewInd = ContentDateTimes.index(max(ContentDateTimes))
        
        # The index of the most recent file that matches:
        DroInd = Inds[NewInd]
        
        if LogToConsole:
            print(f'  The index of the most recent file that matches = {DroInd}')
        
        Dro = Dros[DroInd]
        
        #if 'Deformable' in DroType:
        if Transform == 'bspline':
            TxMatrix = None
            GridDims = GridDimss[NewInd]
            GridRes = GridRess[NewInd]
            VectGridData = VectGridDatas[NewInd]
            
            GridRes_rounded = [round(item, 2) for item in GridRes]
            
            msg = f"  There were {len(Fnames)} subject assessors found, of "\
                  + f"which {len(ContentDateTimes)} was a {DroType}, with "\
                  + "\n  FrameOfReferenceUIDs matching those of the Source "\
                  + "and Target image series, matching the transform type, "\
                  + "\n  the most recent of which contains the grid "\
                  + f"dimensions {GridDims}, grid resolution {GridRes_rounded},"\
                  + f" and vector grid data containing {len(VectGridData)} "\
                  + "elements.\n"
        #else:
        if Transform in ['rigid', 'affine']:
            TxMatrix = TxMatrices[NewInd]
            GridDims = None
            GridRes = None
            VectGridData = None
            
            TxMatrix_rounded = [round(item, 3) for item in TxMatrix]
            
            msg = f"  There were {len(Fnames)} subject assessors found, of "\
                  + f"which {len(ContentDateTimes)} was a {DroType}, with "\
                  + "\n  FrameOfReferenceUIDs matching those of the Source and"\
                  + " Target image series, matching the transform type, "\
                  + "\n  the most recent of which contains the transformation"\
                  + f" matrix \n  {TxMatrix_rounded}\n"
    else:
        msg = "  There were no DROs with FrameOfReferenceUIDs matching those "\
              + "of the Source ({SrcFORuid}) \nand Target ({TrgFORuid}) "\
              + "image series, and whose SOP Class matched the desired "\
              + f"transform type ({Transform}).\n"
        
        TxMatrix = None
        GridDims = None
        GridRes = None
        VectGridData = None
    
    print(msg)
    
    if LogToConsole:
        print('-'*120)
    
    return Dro, TxMatrix, GridDims, GridRes, VectGridData






def GetProjectIds(XnatUrl, XnatSession):
    """
    Get a list of project IDs for all projects in an XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    XnatSession : Requests Session
        A valid requests session.
    
    
    Outputs:
    *******
    
    ProjectIds : list of strings
        A list (for each project) of the project name.
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/projects')
    
    Projects = request.json()['ResultSet']['Result']
    
    ProjectIds = [project['ID'] for project in Projects]

    return ProjectIds




def GetProjectNames(XnatUrl, XnatSession):
    """
    Get a list of project names for all projects in an XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    XnatSession : Requests Session
        A valid requests session.
    
    
    Outputs:
    *******
    
    ProjectNames : list of strings
        A list (for each project) of the project name.
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/projects')
    
    Projects = request.json()['ResultSet']['Result']
    
    ProjectNames = [project['name'] for project in Projects]

    return ProjectNames





def GetInvestigatorsByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing the PIs and co-investigators organised by
    project.  If passed a dictionary (with projects as keys), the investigator
    details will be added to existing data (e.g. subjects organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    Can get list of investigators from:
      request = session.get(f'{XnatAddress}/data/investigators')
      investigators = request.json()
      investigators['ResultSet']['Result']
      
     but resulting list of dictionaries don't include projects.
     
     Instead:
        request = session.get(f'{XnatAddress}/xapi/investigators')
        investigators = request.json()
    
    returns list of dictionaries containing details of project, and whether the
    investigator is a PI or co-investigator for the project.
    
    A flat list (not organised by project) of investigators (PIs and CIs):
        Investigators = [f"{item['firstname']} {item['lastname']}" for item in investigators]
    
    Similarly a flat list of PIs:
        request = XnatSession.get(f'{XnatUrl}/data/projects')
        projects = request.json()
        PIs = [f"{item['pi_firstname']} {item['pi_lastname']}" for item in projects['ResultSet']['Result']]
    """
    
    request = XnatSession.get(f'{XnatUrl}/xapi/investigators')
    
    investigators = request.json()
    
    if DataByProj == None:
        DataByProj = {}
    
    for investigator in investigators:
        firstname = investigator['firstname']
        lastname = investigator['lastname']
        priProjects = investigator['primaryProjects']
        coProjects = investigator['investigatorProjects']
        
        if priProjects:
            for project in priProjects:
                if project in list(DataByProj.keys()):
                    DataByProj[project].update({'PI' : f'{firstname} {lastname}'})
                else:
                    DataByProj[project] = {'PI' : f'{firstname} {lastname}'}
        
        if coProjects:
            for project in coProjects:
                if project in list(DataByProj.keys()):
                    if 'Co-investigators' in list(DataByProj[project].keys()):
                        DataByProj[project]['Co-investigators'].append(f'{firstname} {lastname}')
                    else:
                        DataByProj[project].update({'Co-investigators' : [f'{firstname} {lastname}']})
                else:
                    DataByProj[project] = {'Co-investigators' : [f'{firstname} {lastname}']}
    
    # Add the keys 'PI' and 'CIs' to any project that did not have a PI or CI:
    for project in list(DataByProj.keys()):
        if not 'PI' in list(DataByProj[project].keys()):
            DataByProj[project].update({'PI' : ''})
        
        if not 'Co-investigators' in list(DataByProj[project].keys()):
            DataByProj[project].update({'Co-investigators' : []})
    
    return DataByProj





def GetSubjectsByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing the number of subjects organised by project.
    If passed a dictionary (with projects as keys), the subject details will be
    added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    A flat list (not organised by project) of subjects:
        Subjects = [item['label'] for item in subjects['ResultSet']['Result']]
    
    The list of subjects will be deleted once the number of subjects has been
    computed.
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/subjects')
    
    subjects = request.json()
    
    if DataByProj == None:
        DataByProj = {}
        
    for subject in subjects['ResultSet']['Result']:
        project = subject['project']
        label = subject['label']
        
        if project in list(DataByProj.keys()):
            if 'Subjects' in list(DataByProj[project].keys()):
                DataByProj[project]['Subjects'].append(label)
            else:
                DataByProj[project].update({'Subjects' : [label]})
        else:
            DataByProj[project] = {'Subjects' : [label]}
    
    # Add the key 'NumOfSubjects' to all projects:
    for project in list(DataByProj.keys()):
        subjects = DataByProj[project]['Subjects']
        
        DataByProj[project].update({'No. of subjects' : len(subjects)})
        
        del DataByProj[project]['Subjects']
    
    return DataByProj




def GetImageSessionsByType_OLD(XnatUrl, XnatSession):
    """
    Return a dictionary of image sessions organised by type (e.g. MR, PET).
    

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.

    Returns
    -------
    ScanByScanType : dictionary
        A dictionary with scan type as keys, containing data organised by 
        scan type.
    
    Notes
    -----
    The keys within each sub-dictionary will be `xnat_imagescandata_id', e.g.:
    ImSessionByType['CT']['202'] = {'Project' : 'project_name',
                                    'ExpType' : 'experiment_type',
                                    'ExpLabel' : 'experiment_label',
                                    'ExpId' : 'experiment_id',
                                    'Date' : 'experiment_date',
                                    'SeriesDesc' : 'series_description'
                                    }
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()
    
    ImSessionByType = {}
    
    for experiment in experiments['ResultSet']['Result']:
        expType = experiment['xsiType']
        project = experiment['project']
        expLabel = experiment['label']
        expId = experiment['ID']
        date = experiment['date']
        
        if False:#True:
            print(f' expType = {expType}')
            print(f' project = {project}')
            print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            print(f' date = {date}')
            print('')
        
        if 'SessionData' in expType:
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}/scans')
    
            #print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' expLabel = {expLabel}')
                print(f' expId = {expId}')
                print(f' date = {date}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                ScanDict = {}
    
                scanType = scan['xsiType']
                scanMod = expType.split('xnat:')[1].split('ScanData')[0].upper()
    
                imageScanDataId = scan['xnat_imagescandata_id']
                seriesDesc = scan['series_description']
                scanId = scan['ID']
    
                ScanDict = {'Project' : project,
                            'ExpType' : expType,
                            'ExpLabel' : expLabel,
                            'ExpId' : expId,
                            'Date' : date,
                            'ScanId' : scanId,
                            'SeriesDesc' : seriesDesc}
                
                if scanType in list(ImSessionByType.keys()):
                    ImSessionByType[scanType].update({imageScanDataId : ScanDict})
                else:
                    ImSessionByType.update({scanType : {imageScanDataId : ScanDict}
                                           })
    
    return ImSessionByType





def GetImageSessionsByType(XnatUrl, XnatSession):
    """
    Return a dictionary of image sessions organised by type (e.g. MR, PET).
    

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.

    Returns
    -------
    ScanByScanType : dictionary
        A dictionary with scan type as keys, containing data organised by 
        scan type.
    
    Notes
    -----
    The keys within each sub-dictionary will be `xnat_imagescandata_id', e.g.:
    ImSessionByType['CT'] = {'Project' : 'project_name',
                             'ExpType' : 'experiment_type',
                             'ExpLabel' : 'experiment_label',
                             'ExpId' : 'experiment_id',
                             'Date' : 'experiment_date',
                             'ScanType' : scanType,
                             'ScanId' : imageScanDataId,
                             'SeriesDesc' : 'series_description'
                             }
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()
    
    ImSessionByType = {}
    
    for experiment in experiments['ResultSet']['Result']:
        expType = experiment['xsiType']
        project = experiment['project']
        expLabel = experiment['label']
        expId = experiment['ID']
        date = experiment['date']
        
        if True:
            print(f' expType = {expType}')
            print(f' project = {project}')
            print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            print(f' date = {date}')
            print('')
        
        if 'SessionData' in expType:
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}/scans')
    
            print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' expLabel = {expLabel}')
                print(f' expId = {expId}')
                print(f' date = {date}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                ScanDict = {}
    
                scanType = scan['xsiType']
                scanMod = scanType.split('xnat:')[1].split('ScanData')[0].upper()
                
                print(f'   scanType = {scanType}')
                print(f'   scanMod = {scanMod}')
    
                imageScanDataId = scan['xnat_imagescandata_id']
                seriesDesc = scan['series_description']
    
                ScanDict = {'Project' : project,
                            'ExpType' : expType,
                            'ExpLabel' : expLabel,
                            'ExpId' : expId,
                            'Date' : date,
                            'ScanType' : scanType,
                            'ScanId' : imageScanDataId,
                            'SeriesDesc' : seriesDesc}
    
                #if scanType in list(ImSessionByType.keys()):
                if scanMod in list(ImSessionByType.keys()):
                    #ImSessionByType[scanType].append(ScanDict)
                    ImSessionByType[scanMod].append(ScanDict)
                    
                else:
                    #ImSessionByType.update({scanType : [ScanDict]})
                    #ImSessionByType.update({scanMod : [ScanDict]}) # equivalent to below
                    ImSessionByType[scanMod] = [ScanDict]
    
    return ImSessionByType





def GetExperimentsByType_INCOMPLETE(XnatUrl, XnatSession):
    """
    Return a dictionary of experiments organised by type (e.g. MR, PET, Image
    Scan Data, ROI Collection, etc.).
    

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.

    Returns
    -------
    ExperimentsByType : dictionary
        A dictionary with scan type as keys, containing data organised by 
        experiment type.
    
    Notes
    -----
    The keys within each sub-dictionary will be `xnat_imagescandata_id', e.g.:
    ExperimentsByType['CT'] = {'Project' : 'project_name',
                               'ExpType' : 'experiment_type',
                               'ExpLabel' : 'experiment_label',
                               'ExpId' : 'experiment_id',
                               'Date' : 'experiment_date',
                               'ScanType' : scanType,
                               'ScanId' : imageScanDataId,
                               'SeriesDesc' : 'series_description'
                               }
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()
    
    ExperimentsByType = {}
    
    for experiment in experiments['ResultSet']['Result']:
        expType = experiment['xsiType']
        project = experiment['project']
        expLabel = experiment['label']
        expId = experiment['ID']
        date = experiment['date']
        uri = experiment['URI']
        
        if False:#True:
            print(f' expType = {expType}')
            print(f' project = {project}')
            print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            print(f' date = {date}')
            print('')
        
        """ The following will return detailed info about the experiment that
        might be useful later:
        
        request = XnatSession.get(f'{XnatUrl}{uri}?format=json')
        
        experiment = request.json()

        print(f"keys in experiment = {list(experiment.keys())}\n")
        
        print(f"len(experiment['items'] = {len(experiment['items'])}\n")
        
        print(f"keys in experiment['items'][0] = {list(experiment['items'][0].keys())}\n")
        
        experiment['items'][0]
        """

        # continue...
    
    return ExperimentsByType



def GetNumOfExpsByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary of the number of experiments organised by project.
    If passed a dictionary (with projects as keys), the number of experiments
    will be added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    ExperimentsByType : dictionary
        A dictionary with scan type as keys, containing data organised by 
        experiment type.
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()
    
    if DataByProj == None:
        DataByProj = {}
    
    for experiment in experiments['ResultSet']['Result']:
        #expType = experiment['xsiType']
        project = experiment['project']
        #expLabel = experiment['label']
        #expId = experiment['ID']
        #date = experiment['date']
        #uri = experiment['URI']
        
        #if False:#True:
        #    print(f' expType = {expType}')
        #    print(f' project = {project}')
        #    print(f' expLabel = {expLabel}')
        #    print(f' expId = {expId}')
        #    print(f' date = {date}')
        #    print('')
        
        """ The following will return detailed info about the experiment that
        might be useful later:
        
        request = XnatSession.get(f'{XnatUrl}{uri}?format=json')
        
        experiment = request.json()

        print(f"keys in experiment = {list(experiment.keys())}\n")
        
        print(f"len(experiment['items'] = {len(experiment['items'])}\n")
        
        print(f"keys in experiment['items'][0] = {list(experiment['items'][0].keys())}\n")
        
        experiment['items'][0]
        """

        if project in list(DataByProj.keys()):
            if 'No. of Experiments' in list(DataByProj[project].keys()):
                #print(f'     before = {DataByProj[project]['No. of Experiments']}')
                DataByProj[project]['No. of Experiments'] += 1
                #print(f'     after = {DataByProj[project]['No. of Experiments']}')
                
            else:
                #print('     before = 0')
                DataByProj[project]['No. of Experiments'] = 1
                #print(f'     after = {DataByProj[project]['No. of Experiments']}')
        else:
            #print('     before = 0')
            #DataByProj[project]['No. of Experiments'] = 1
            DataByProj[project] = {'No. of Experiments' : 1}
            #print(f'     after = {DataByProj[project]['No. of Experiments']}')
    
    return DataByProj




def GetImSessionByTypeByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing image sessions organised by image type by
    project. If passed a dictionary (with projects as keys), the image session 
    details will be added to existing data (e.g. investigators organised by 
    project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    The keys within each sub-dictionary will be e.g.:
    DataByProj['ProjectName]['CT'] = {'ExpType' : 'experiment_type',
                                      'ExpLabel' : 'experiment_label',
                                      'ExpId' : 'experiment_id',
                                      'Date' : 'experiment_date',
                                      'ScanType' : scanType,
                                      'ScanId' : imageScanDataId,
                                      'SeriesDesc' : 'series_description'
                                      }
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()['ResultSet']['Result']
    
    if DataByProj == None:
        DataByProj = {}
    
    for experiment in experiments:
        expType = experiment['xsiType']
        project = experiment['project']
        expLabel = experiment['label']
        expId = experiment['ID']
        expDate = experiment['date']
        insertDate = experiment['insert_date']
        
        if True:
            print(f' expType = {expType}')
            print(f' project = {project}')
            print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            print(f' expDate = {expDate}')
            print(f' insertDate = {insertDate}')
            print('')
        
        if 'SessionData' in expType:
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}/scans')
    
            print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' expLabel = {expLabel}')
                print(f' expId = {expId}')
                print(f' expDate = {expDate}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                ScanDict = {}
    
                scanType = scan['xsiType']
                scanMod = scanType.split('xnat:')[1].split('ScanData')[0].upper()
                
                print(f'   scanType = {scanType}')
                print(f'   scanMod = {scanMod}')
    
                imageScanDataId = scan['xnat_imagescandata_id']
                seriesDesc = scan['series_description']
    
                ScanDict = {'ExpType' : expType,
                            'ExpLabel' : expLabel,
                            'ExpId' : expId,
                            'ExpDate' : expDate,
                            'ScanType' : scanType,
                            'ScanId' : imageScanDataId,
                            'UploadDate' : insertDate,
                            'SeriesDesc' : seriesDesc}
                
                if project in list(DataByProj.keys()):
                    if scanMod in list(DataByProj[project].keys()):
                        DataByProj[project][scanMod].append(ScanDict)
                        
                    else:
                        #DataByProj[project].update({scanMod : [ScanDict]}) # equivalent to below
                        DataByProj[project][scanMod] = [ScanDict]
                else:
                    DataByProj[project] = {scanMod : [ScanDict]}
        
    
    return DataByProj





def GetNumOfImSessionByTypeByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing the number of image sessions organised by 
    image type by project. If passed a dictionary (with projects as keys), the 
    image session details will be added to existing data (e.g. investigators 
    organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    It isn't clear what are the image session types available within XNAT. At
    present that info is not available, e.g. see:
    https://wiki.xnat.org/documentation/xnat-administration/managing-data-types-in-xnat/adding-and-managing-image-scan-data-types
    """
    
    #AllDataByProj = GetImSessionByTypeByProject(XnatUrl, XnatSession, DataByProj)
    
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()
    
    if DataByProj == None:
        DataByProj = {}
    
    for experiment in experiments['ResultSet']['Result']:
        expType = experiment['xsiType']
        project = experiment['project']
        expLabel = experiment['label']
        expId = experiment['ID']
        expDate = experiment['date']
        insertDate = experiment['insert_date']
        
        if False:#True:
            print(f' expType = {expType}')
            print(f' project = {project}')
            print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            print(f' expDate = {expDate}')
            print(f' insertDate = {insertDate}')
            print('')
        
        if 'SessionData' in expType:
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}/scans')
    
            #print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' expLabel = {expLabel}')
                print(f' expId = {expId}')
                print(f' expDate = {expDate}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                scanType = scan['xsiType']
                scanMod = scanType.split('xnat:')[1].split('ScanData')[0].upper()
                
                key = f"No. of {scanMod} scans"
                
                #print(f'   scanType = {scanType}')
                #print(f'   scanMod = {scanMod}')
                #print(f'   key = {key}')
                
                if project in list(DataByProj.keys()):
                    if key in list(DataByProj[project].keys()):
                        #print(f'     before = {DataByProj[project][key]}')
                        DataByProj[project][key] += 1
                        #print(f'     after = {DataByProj[project][key]}')
                        
                    else:
                        #print('     before = 0')
                        DataByProj[project][key] = 1
                        #print(f'     after = {DataByProj[project][key]}')
                else:
                    #print('     before = 0')
                    #DataByProj[project][key] = 1
                    DataByProj[project] = {key : 1}
                    #print(f'     after = {DataByProj[project][key]}')
        
    
    return DataByProj





def GetExpsSessionsScansByProject_OLD(XnatUrl, XnatSession, DataByProj=None):
    """
    18/06/21: Expansion of combination of GetNumOfExpsByProject and 
    GetNumOfImSessionByTypeByProject.
    
    21/06/21: Error in calculation of the number of experiments.  Numbers seem
    too high. Also this function is redundant due to the newer function
    GetFileCountAndSizeByTypeByProject, which iterates through all projects,
    then through all experiments, rather than in this function, which iterates
    through all experiments, then gathers info by common project (probably more
    error prone).
    
    
    Return a dictionary containing the number of experiments, the number of
    image sessions per modality, the number of image scans per modality
    and the average number of image scans by image session, all organised by 
    project. If passed a dictionary (with projects as keys), the above details
    will be added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    It isn't clear what are the image session types available within XNAT. At
    present that info is not available, e.g. see:
    https://wiki.xnat.org/documentation/xnat-administration/managing-data-types-in-xnat/adding-and-managing-image-scan-data-types
    
    Edit: Check here:
    https://wiki.xnat.org/documentation/how-to-use-xnat/understanding-the-xnat-data-model
    
    The following will return detailed info about the experiment that might be
    useful later:
        
    request = XnatSession.get(f'{XnatUrl}{uri}?format=json')
    
    experiment = request.json()

    print(f"keys in experiment = {list(experiment.keys())}\n")
    
    print(f"len(experiment['items'] = {len(experiment['items'])}\n")
    
    print(f"keys in experiment['items'][0] = {list(experiment['items'][0].keys())}\n")
    
    experiment['items'][0]
    """
    
    #AllDataByProj = GetImSessionByTypeByProject(XnatUrl, XnatSession, DataByProj)
    
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()['ResultSet']['Result']
    
    if DataByProj == None:
        DataByProj = {}
    
    for experiment in experiments:
        expType = experiment['xsiType']
        project = experiment['project']
        expLabel = experiment['label']
        expId = experiment['ID']
        expDate = experiment['date']
        insertDate = experiment['insert_date']
        
        if False:#True:
            print(f' expType = {expType}')
            print(f' project = {project}')
            print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            print(f' expDate = {expDate}')
            print(f' insertDate = {insertDate}')
            print('')
        
        
        """ Proceed only if an image session: """
        if 'SessionData' in expType:
            
            """ Add the number of image session experiments: 
            (i.e. excludes icr:roiCollectionData) """
            if project in list(DataByProj.keys()):
                if 'No. of Experiments' in list(DataByProj[project].keys()):
                    #print(f'     before = {DataByProj[project]['No. of Experiments']}')
                    DataByProj[project]['No. of Experiments'] += 1
                    #print(f'     after = {DataByProj[project]['No. of Experiments']}')
                    
                else:
                    #print('     before = 0')
                    DataByProj[project]['No. of Experiments'] = 1
                    #print(f'     after = {DataByProj[project]['No. of Experiments']}')
            else:
                #print('     before = 0')
                DataByProj[project] = {'No. of Experiments' : 1}
                #print(f'     after = {DataByProj[project]['No. of Experiments']}')
            
            
            
            if False:#True:#False:#'500' in str(request):
                print(f' project = {project}')
                print(f' expLabel = {expLabel}')
                print(f' expId = {expId}')
                print(f' expDate = {expDate}')
    
            
            sessionMod = expType.split('xnat:')[1].split('SessionData')[0].upper()
            
            key = f"No. of {sessionMod} sessions"
            
            """ Get the number of this session type: """
            if project in list(DataByProj.keys()):
                if f'No. of {sessionMod} sessions' in list(DataByProj[project].keys()):
                    #print(f'     before = {DataByProj[project][key]}')
                    DataByProj[project][key] += 1
                    #print(f'     after = {DataByProj[project][key]}')
                    
                else:
                    #print('     before = 0')
                    DataByProj[project][key] = 1
                    #print(f'     after = {DataByProj[project][key]}')
            else:
                #print('     before = 0')
                DataByProj[project] = {key : 1}
                #print(f'     after = {DataByProj[project][key]}')
            
            
            
            """ Get info on scans: """
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}/scans')
    
            #print(request, '\n')
            
            scans = request.json()
            
            """ Get the number of scans by type: """
            for scan in scans['ResultSet']['Result']:
                scanType = scan['xsiType']
                scanMod = scanType.split('xnat:')[1].split('ScanData')[0].upper()
                
                
                """ Exclude RT image scan data: """
                if not scanMod == 'RTIMAGE':
                    if False:
                        print(f'   scanType = {scanType}')
                        print(f'   scanMod = {scanMod}')
                    
                    key = f"No. of {scanMod} scans"
                    
                    #print(f'   scanType = {scanType}')
                    #print(f'   scanMod = {scanMod}')
                    #print(f'   key = {key}')
                    
                    if project in list(DataByProj.keys()):
                        if key in list(DataByProj[project].keys()):
                            #print(f'     before = {DataByProj[project][key]}')
                            DataByProj[project][key] += 1
                            #print(f'     after = {DataByProj[project][key]}')
                            
                        else:
                            #print('     before = 0')
                            DataByProj[project][key] = 1
                            #print(f'     after = {DataByProj[project][key]}')
                    else:
                        #print('     before = 0')
                        DataByProj[project] = {key : 1}
                        #print(f'     after = {DataByProj[project][key]}')
        
    
    """ Get the average number of scans by project: """
    for project in list(DataByProj.keys()):
        NumSessions = 0
        NumScans = 0
        for key in list(DataByProj[project].keys()):
            if 'sessions' in key:
                NumSessions += DataByProj[project][key]
            
            if 'scans' in key:
                NumScans += DataByProj[project][key]
        
        AveScans = round(NumScans/NumSessions, 1)
        
        DataByProj[project]['Ave. scans/session'] = AveScans
        
    
    return DataByProj





def GetFileCountAndSizeByTypeByProject_OLD(XnatUrl, XnatSession, DataByProj=None,
                                       LogToConsole=False):
    """
    18/06/21: Expansion of combination of GetNumOfExpsByProject and 
    GetNumOfImSessionByTypeByProject.
    
    Return a dictionary containing the number of files and the total size of
    files by modality, all organised by project.  If passed a dictionary (with 
    projects as keys), the above details will be added to existing data 
    (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.
    LogToConsole : boolean (optional; False by default)
        If True some info will be printed to the console.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
        
    request = XnatSession.get(f'{XnatUrl}{uri}?format=json')
    
    or
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}?format=json')
    
    Exp = request.json()

    
    Exp is a dictionary with one key: 'items'.
        
    Exp['items'] is a list of length 1.
    
    Exp['items'][0] is a dictionary with keys: 'children', 'meta'
    and 'data_fields'.
    
    Exp['items'][0]['meta'] and Exp['items'][0]['data_fields']
    are simple dictionaries.
    
    Exp['items'][0]['children'] is a list of variable length. It
    seems that if there are image assessors, 
    Exp['items'][0]['children'][0] relate to assessors, and
    Exp['items'][0]['children'][1] relate to scans.
    
    Exp['items'][0]['children'][i] is a dictionary with keys:
    'field' and 'items'.
    
    Exp['items'][0]['children'][i]['field'] is either 
    'assessors/assessor' or 'scans/scan'.
    
    Exp['items'][0]['children'][i]['items'] is a list for each
    assessor/scan.
    
    Each item in Exp['items'][0]['children'][i]['items'] is a 
    dictionary with keys: 'children', 'meta' and 'data_fields'.
    
    Exp['items'][0]['children'][i]['items'][j]['meta'] and 
    Exp['items'][0]['children'][i]['items'][j]['data_fields'] are
    simple dictionaries.
    
    Exp['items'][0]['children'][i]['items'][j]['children'] is a
    list of length 2.
    
    Exp['items'][0]['children'][i]['items'][j]['children'][0] is a
    dictionary with keys 'field' and 'items'.
    
    Exp['items'][0]['children'][i]['items'][j]['children'][0]['field']
    seems to always be 'out/file' (for assessors) or 'file' (for scans).
    
    If the item is an assessor:
        If the assessor is RTSTRUCT:
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items']
            is a list of length 2.
             
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][0]
            relates to the AIM, and 
            
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][1]
            relates to the RTSTRUCT.
    
        If the assessor is SEG:
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items']
            is a list of length 1.
             
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][0]
            relates to the SEG.
        
    
    Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][k] 
    is a dictionary with keys 'children' (= []), 'meta' and 'data_fields'.
    
    The file_size is located in
    Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][k]['data_fields'].
    
    Exp['items'][0]['children'][i]['items'][j]['children'][1] is a
    simple dictionary that references the SeriesUID, and has keys 'field' 
    and 'items'.  
    
    Exp['items'][0]['children'][i]['items'][j]['children'][1]['field']
    seems to always be 'references/seriesUID'.
    
    expeExpriment['items'][0]['children'][i]['items'][j]['children'][1]['items']
    is a list of length 1.
    
    Exp['items'][0]['children'][i]['items'][j]['children'][1]['items'][0]
    is a dictionary with keys 'children' (= []), 'meta' and 'data_fields'.
    """
    
    if DataByProj == None:
        DataByProj = {}
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    Exps = request.json()['ResultSet']['Result']
    
    if LogToConsole:
        print(f'There are {len(Exps)} experiments.')
    
    for ExpNum in range(len(Exps)):
        if LogToConsole:
            print('--------------------------------------------------------')
        
        ExpType = Exps[ExpNum]['xsiType']
        Project = Exps[ExpNum]['project']
        ExpLabel = Exps[ExpNum]['label']
        ExpId = Exps[ExpNum]['ID']
        ExpDate = Exps[ExpNum]['date']
        InsertDate = Exps[ExpNum]['insert_date']
        
        if LogToConsole:
            print(f'Experiment {ExpNum}:')
            print(f'  ExpType = {ExpType}')
            print(f'  Project = {Project}')
            print(f'  ExpLabel = {ExpLabel}')
            print(f'  ExpId = {ExpId}')
            print(f'  ExpDate = {ExpDate}')
            print(f'  InsertDate = {InsertDate}')
            print('')
        
        """ It seems that an image session experiment references ROI 
        Collections that correspond to that imaging session, and an ROI 
        Collection experiment references image sessions that correspond to that
        ROI Collection. To avoid redundancy, proceed only for experiments that 
        are not ROI Collections: """
        if ExpType == 'icr:roiCollectionData':
            if LogToConsole:
                print(f"Skipping experiment {ExpNum} since it is type {ExpType}.\n")
        else:
            """ Get more detailed info for this experiment: """
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{ExpId}?format=json')

            Exp = request.json()['items'][0]['children']
        
            #print(f'  Experiment {ExpNum} has {len(Exp)} categories (i.e. scans and/or assessors).\n')
    
            # Loop through each experiment "category" ("category" here means  
            # scans or assessor (if applicable)):
            for CatNum in range(len(Exp)):
                if 'assessor' in Exp[CatNum]['field']:
                    ExpCat = 'assessor'
                elif 'scan' in Exp[CatNum]['field']:
                    ExpCat = 'scan'
                else:
                    msg = f"Unexpected field value {Exp[CatNum]['field']}. Was expecting "\
                          + "assessors/assessors or scans/scan."
                    raise Exception(msg)
                
                if LogToConsole:
                    #print(f'    ExpCat = {ExpCat}\n')
                    print(f"  There are {len(Exp[CatNum]['items'])} {ExpCat}s")
    
                # Loop through each assessor/scan:
                for ItemNum in range(len(Exp[CatNum]['items'])):
    
                    if 'file' in Exp[CatNum]['items'][ItemNum]['children'][0]['field']:
                        # NumFiles = 1 for scans or SEG, 2 for RTSTRUCT
                        NumFiles = len(Exp[CatNum]['items'][ItemNum]['children'][0]['items'])
                        
                        if LogToConsole:
                            print(f"    {ExpCat} {ItemNum} has {NumFiles} items")
    
                        for FileNo in range(NumFiles):
                            # FileLabel = AIM, RTSTRUCT, SEG or DICOM:
                            FileLabel = Exp[CatNum]['items'][ItemNum]['children'][0]['items'][FileNo]['data_fields']['label']
                            
                            if LogToConsole:
                                print(f'      File {FileNo} has {FileLabel} label')
                                
                            """ It seems that file size is not listed for
                            SNAPSHOTS so skip them: """
                            
                            if FileLabel == 'SNAPSHOTS':
                                if LogToConsole:
                                    print(f"Skipping file no. {FileNo} since",
                                          f"it is type {FileLabel}.\n")
                            else:
                                
                                FileCount = Exp[CatNum]['items'][ItemNum]['children'][0]['items'][FileNo]['data_fields']['file_count']
                                FileFormat = Exp[CatNum]['items'][ItemNum]['children'][0]['items'][FileNo]['data_fields']['format']
                                FileSize = Exp[CatNum]['items'][ItemNum]['children'][0]['items'][FileNo]['data_fields']['file_size']
                                
                                if LogToConsole:
                                    print(f'      File {FileNo} has {FileCount} files')
                                    print(f'      File {FileNo} has {FileFormat} format')
                                    print(f'      File {FileNo} has {FileSize} file size\n')
                                
                                FileCountKey = f"No. of {FileLabel} files"
                                FileSizeKey = f"Total {FileLabel} file size"
                                
                                if Project in DataByProj.keys():
                                    if FileCountKey in DataByProj[Project].keys():
                                        DataByProj[Project][FileCountKey] += FileCount
                                    else:
                                        DataByProj[Project][FileCountKey] = FileCount
                                    
                                    if FileSizeKey in DataByProj[Project].keys():
                                        DataByProj[Project][FileSizeKey] += FileSize/1000000
                                    else:
                                        DataByProj[Project][FileSizeKey] = FileSize/1000000
                                else:
                                    DataByProj[Project] = {FileCountKey : FileCount,
                                                        FileSizeKey : FileSize/1000000}
        
        
        
        if LogToConsole:
            print('--------------------------------------------------------\n')
    
    
    """ Round the file sizes to 1 decimal: """
    for project in DataByProj.keys():
        for key, val in DataByProj[project].items():
            if 'file size' in key:
                DataByProj[project][key] = round(val, 1)
        
    return DataByProj





def GetFileCountAndSizeByTypeByProject(XnatUrl, XnatSession, DataByProj=None,
                                       LogToConsole=False):
    """
    18/06/21: Expansion of combination of GetNumOfExpsByProject and 
    GetNumOfImSessionByTypeByProject, and replaces the latter.
    
    Return a dictionary containing the number of files and the total size of
    files by modality, all organised by project.  If passed a dictionary (with 
    projects as keys), the above details will be added to existing data 
    (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.
    LogToConsole : boolean (optional; False by default)
        If True some info will be printed to the console.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
        
    request = XnatSession.get(f'{XnatUrl}{uri}?format=json')
    
    or
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}?format=json')
    
    Exp = request.json()

    
    Exp is a dictionary with one key: 'items'.
        
    Exp['items'] is a list of length 1.
    
    Exp['items'][0] is a dictionary with keys: 'children', 'meta'
    and 'data_fields'.
    
    Exp['items'][0]['meta'] and Exp['items'][0]['data_fields']
    are simple dictionaries.
    
    Exp['items'][0]['children'] is a list of variable length. It
    seems that if there are image assessors, 
    Exp['items'][0]['children'][0] relate to assessors, and
    Exp['items'][0]['children'][1] relate to scans.
    
    Exp['items'][0]['children'][i] is a dictionary with keys:
    'field' and 'items'.
    
    Exp['items'][0]['children'][i]['field'] is either 
    'assessors/assessor' or 'scans/scan'.
    
    Exp['items'][0]['children'][i]['items'] is a list for each
    assessor/scan.
    
    Each item in Exp['items'][0]['children'][i]['items'] is a 
    dictionary with keys: 'children', 'meta' and 'data_fields'.
    
    Exp['items'][0]['children'][i]['items'][j]['meta'] and 
    Exp['items'][0]['children'][i]['items'][j]['data_fields'] are
    simple dictionaries.
    
    Exp['items'][0]['children'][i]['items'][j]['children'] is a
    list of length 2.
    
    Exp['items'][0]['children'][i]['items'][j]['children'][0] is a
    dictionary with keys 'field' and 'items'.
    
    Exp['items'][0]['children'][i]['items'][j]['children'][0]['field']
    seems to always be 'out/file' (for assessors) or 'file' (for scans).
    
    If the item is an assessor:
        If the assessor is RTSTRUCT:
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items']
            is a list of length 2.
             
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][0]
            relates to the AIM, and 
            
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][1]
            relates to the RTSTRUCT.
    
        If the assessor is SEG:
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items']
            is a list of length 1.
             
            Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][0]
            relates to the SEG.
        
    
    Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][k] 
    is a dictionary with keys 'children' (= []), 'meta' and 'data_fields'.
    
    The file_size is located in
    Exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][k]['data_fields'].
    
    Exp['items'][0]['children'][i]['items'][j]['children'][1] is a
    simple dictionary that references the SeriesUID, and has keys 'field' 
    and 'items'.  
    
    Exp['items'][0]['children'][i]['items'][j]['children'][1]['field']
    seems to always be 'references/seriesUID'.
    
    expeExpriment['items'][0]['children'][i]['items'][j]['children'][1]['items']
    is a list of length 1.
    
    Exp['items'][0]['children'][i]['items'][j]['children'][1]['items'][0]
    is a dictionary with keys 'children' (= []), 'meta' and 'data_fields'.
    """
    
    #from copy import deepcopy
    
    if DataByProj == None:
        DataByProj = {}
    
    EmptyDict = {'No. of experiments' : 0,
                 'No. of MR sessions' : 0,
                 'No. of CT sessions' : 0,
                 'No. of PET sessions' : 0,
                 'Ave. image scans/session' : 0,
                 #'Ave. scans/session' : 0,
                 'No. of MR scans' : 0,
                 'No. of CT scans' : 0,
                 'No. of PET scans' : 0,
                 'No. of OT scans' : 0,
                 'No. of DICOM image files [k]' : 0,
                 'Size of DICOM image files [MB]' : 0,
                 'No. of RTSTRUCT scans' : 0,
                 #'Size of RTSTRUCT scans [MB]' : 0,
                 #'No. of AIM ROI files' : 0,
                 #'Size of AIM ROI files [MB]' : 0,
                 'No. of RTSTRUCT ROI files' : 0,
                 #'Size of RTSTRUCT ROI files [MB]' : 0,
                 'No. of SEG ROI files' : 0,
                 #'Size of SEG ROI files [MB]' : 0,
                 'No. of DICOM files [k]' : 0,
                 'Size of DICOM files [MB]' : 0
                 }
    
    ProjIds = GetProjectIds(XnatUrl, XnatSession)
    
    for ProjId in ProjIds:
        if ProjId in DataByProj.keys():
            DataByProj[ProjId].update(EmptyDict)
        else:
            #DataByProj[ProjId] = {}
            #DataByProj[ProjId] = EmptyDict
            #DataByProj[ProjId] = deepcopy(EmptyDict)
            DataByProj[ProjId] = dict(EmptyDict) # or EmptyDict.copy()
            
            #print('EmptyDict = ', EmptyDict, '\n')
            
        
        # Get all experiments for this ProjId:
        request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/experiments')
        
        Exps = request.json()['ResultSet']['Result']
        
        if LogToConsole:
            print('****************************************************')
            print(f"Project {ProjId} has {len(Exps)} experiments.\n")
        
        key = 'No. of experiments'
        DataByProj[ProjId][key] = len(Exps)
            
        URIs = [Exp['URI'] for Exp in Exps]
        
        NumSessionsThisProj = 0 # initial value of no. of image sessions for this project
        NumScansThisProj = 0 # initial value of the no. of scans for this project
        NumImScansThisProj = 0 # initial value of the no. of image scans for this project
        
        # Loop through each experiment:
        for ExpNo in range(len(Exps)):
            ExpType = Exps[ExpNo]['xsiType']
            ExpMod = ExpType.split('xnat:')[1].split('SessionData')[0].upper()
            #Project = Exps[ExpNo]['project']
            ExpLabel = Exps[ExpNo]['label']
            ExpId = Exps[ExpNo]['ID']
            ExpDate = Exps[ExpNo]['date']
            InsertDate = Exps[ExpNo]['insert_date']
        
            if LogToConsole:
                print('----------------------------------------------------')
                print(f"Experiment {ExpNo}:\n")
                print(f' ExpType = {ExpType}')
                print(f' ExpMod = {ExpMod}')
                print(f' Project = {ProjId}')
                print(f' ExpLabel = {ExpLabel}')
                print(f' ExpId = {ExpId}')
                print(f' ExpDate = {ExpDate}')
                print(f' InsertDate = {InsertDate}\n')
                
            
            key = f'No. of {ExpMod} sessions'
            if key in DataByProj[ProjId].keys():
                DataByProj[ProjId][key] += 1
            else:
                DataByProj[ProjId][key] = 1
                
            if 'Session' in ExpType:
                NumSessionsThisProj += 1
            
            #URI = URIs[ExpNo]
            
            # Get more details for this experiment: 
            request = XnatSession.get(f'{XnatUrl}{URIs[ExpNo]}?format=json')
            
            Exp = request.json()['items'][0]['children']
            
            if LogToConsole:
                print(f" Experiment {ExpNo} has {len(Exp)} items.\n")
                
            
            for ExpItemNo in range(len(Exp)):
                #ExpItemType = Exp[ExpItemNo]['field'] # assessors/assessor or scans/scan
                
                if 'assessor' in Exp[ExpItemNo]['field']:
                    ExpItemType = 'assessor'
                elif 'scan' in Exp[ExpItemNo]['field']:
                    ExpItemType = 'scan'
                else:
                    msg = f"Unexpected field value {Exp[ExpItemNo]['field']}."\
                          + "Was expecting 'assessors/assessors' or 'scans/scan'."
                    raise Exception(msg)
                
                #print(f"    Experiment item {ExpItemNo} is a {ExpItemType}.\n")
                
                if LogToConsole:
                    print(f"    Item {ExpItemNo} is a {ExpItemType} and has",
                          f"{len(Exp[ExpItemNo]['items'])} items.\n")
                
                # Loop through each assessor/scan sub-item:
                for SubItemNo in range(len(Exp[ExpItemNo]['items'])):
                    # Proceed only for items with field 'out/file' (assessors)
                    # or 'file' (scans) (there are also items with field 
                    # 'references/seriesUID':
                    field = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['field']
                    
                    if LogToConsole:
                        print(f"      Sub-item {SubItemNo} has field {field}.")
                    
                    if 'file' in field:
                        
                        #Label = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][0]['data_fields']['label']
                        #Format = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][0]['data_fields']['format']
                        #
                        #print(f"      Sub-item {SubItemNo} has label {Label}.")
                        #print(f"      Sub-item {SubItemNo} has format {Format}.")
                        
                        # Modality doesn't exist for AIM:
                        #Modality = Exp[ExpItemNo]['items'][SubItemNo]['data_fields']['modality']
                        #try:
                        #    Modality = Exp[ExpItemNo]['items'][SubItemNo]['data_fields']['modality']
                        #except KeyError:
                        #    #return Exp[ExpItemNo]['items'][SubItemNo]
                        #    Modality = 'N/A'
                        
                        """ Seems that modality is not accessible from Exp[ExpItemNo]['items'][SubItemNo]['children']...
                        as are label, format, file_size, etc. (as obtained by looping through each sub-sub-item below). """
                        if ExpItemType == 'assessor':
                            Modality = Exp[ExpItemNo]['items'][SubItemNo]['data_fields']['collectionType']
                        else: # ExpItemType = 'scan'
                            Modality = Exp[ExpItemNo]['items'][SubItemNo]['data_fields']['modality']
                        
                        
                        """ Change 'PT' to 'PET' for consistency: """
                        if Modality == 'PT':
                            Modality = 'PET'
                        
                        N = len(Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'])
                        
                        if LogToConsole:
                            print(f"\n      Sub-item {SubItemNo} has {N} items.")
                        
                        # Loop through each sub-sub-item:
                        for SubSubItemNo in range(N):
                            #if LogToConsole:
                            #    print(Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo], '\n')
                            
                            Label = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo]['data_fields']['label']
                            Format = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo]['data_fields']['format']
                            
                            if LogToConsole:
                                print(f"        Sub-sub-item {SubSubItemNo} has modality {Modality}.")
                                print(f"        Sub-sub-item {SubSubItemNo} has label {Label}.")
                                print(f"        Sub-sub-item {SubSubItemNo} has format {Format}.")
                            
                            
                            # File size is not listed for SNAPSHOTS so skip them:
                            if Label == 'SNAPSHOTS':
                                #return Exp[ExpItemNo]['items'][SubItemNo]
                            
                                if LogToConsole:
                                    print(f"      Skipping SubSubItemNo {SubSubItemNo} since it is type {Label}.\n")
                            
                            else:
                                NumScansThisProj += 1
                                
                                if Label == 'DICOM':
                                    # This is an image scan:
                                    NumImScansThisProj += 1
                                
                                FileCount = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo]['data_fields']['file_count']
                                FileSize = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo]['data_fields']['file_size']
                                
                                if LogToConsole:
                                    print(f"        Sub-sub-item {SubSubItemNo} has FileCount {FileCount}.")
                                    print(f"        Sub-sub-item {SubSubItemNo} has FileSize {FileSize}.\n")
                                    
                                    if Label == 'secondary':
                                        print(Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo], '\n')
                                
                                # Proceed only for sub-sub-items with field 'out/file' (ignore 'references/seriesUID'):
                                if 'file' in Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['field']:
                                    FileCount = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo]['data_fields']['file_count']
                                    FileSize = Exp[ExpItemNo]['items'][SubItemNo]['children'][0]['items'][SubSubItemNo]['data_fields']['file_size']
                                    FileSize = FileSize/1000000 # convert to MB
                                    
                                    #print(f"        Sub-sub-item {SubSubItemNo} has modality {Modality}.")
                                    #print(f"        Sub-sub-item {SubSubItemNo} has FileCount {FileCount}.")
                                    #print(f"        Sub-sub-item {SubSubItemNo} has FileSize {FileSize}.\n")
                                    
                                    #if Label == 'AIM':
                                    #    return Exp[ExpItemNo]['items']#[SubItemNo]
                                    
                                    
                                    
                                    #if Label == 'SEG':
                                    #    return Exp[ExpItemNo]['items'][SubItemNo]
                                    #    print('')
                                    
                                    #if Label == 'secondary':
                                    #    return Exp[ExpItemNo]['items'][SubItemNo]#['children'][0]['items'][0]
                                    #    return Exp[ExpItemNo]['items'][SubItemNo - 1]
                                    #    print('')
                                    
                                    #if Modality == 'OT':
                                    #    return Exp[ExpItemNo]['items'][SubItemNo]
                                    
                                    
                                    
                                    #if Label in ['AIM', 'RTSTRUCT', 'SEG']:
                                    if Label in ['RTSTRUCT', 'SEG']:
                                        key = f'No. of {Label} ROI files'
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += FileCount
                                        else:
                                            DataByProj[ProjId][key] = FileCount
                                        
                                        """ Leave out file size for secondary
                                        ROI Collections for now:
                                        key = f'Size of {Label} ROI files [MB]'
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += FileSize
                                        else:
                                            DataByProj[ProjId][key] = FileSize
                                        """
                                    
                                    elif Label == 'DICOM':
                                        # This is a DICOM scan (images)
                                        
                                        key = 'No. of DICOM image files [k]'
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += FileCount/1000
                                        else:
                                            DataByProj[ProjId][key] = FileCount/1000
                                        
                                        key = 'Size of DICOM image files [MB]'
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += FileSize
                                        else:
                                            DataByProj[ProjId][key] = FileSize
                                        
                                        key = f"No. of {Modality} scans"
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += 1
                                        else:
                                            DataByProj[ProjId][key] = 1        
                                    
                                    if Format == 'DICOM':
                                        # This is a DICOM scan or RTSTRUCT scan or
                                        # ROI Collection of type AIM, RTSTRUCT or SEG:
                                        
                                        key = 'No. of DICOM files [k]'
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += FileCount/1000
                                        else:
                                            DataByProj[ProjId][key] = FileCount/1000
                                        
                                        key = 'Size of DICOM files [MB]'
                                        if key in DataByProj[ProjId].keys():
                                            DataByProj[ProjId][key] += FileSize
                                        else:
                                            DataByProj[ProjId][key] = FileSize
                                        
                                        
                                        
                                        if Label == 'secondary': 
                                            # Secondary DICOMs, e.g. RTSTRUCT or OT scan
                                            """ Replaced FileCount with 1 counting scans
                                            not number of files within scans. """
                                            key = f'No. of {Modality} scans'
                                            if key in DataByProj[ProjId].keys():
                                                #DataByProj[ProjId][key] += FileCount
                                                DataByProj[ProjId][key] += 1
                                            else:
                                                #DataByProj[ProjId][key] = FileCount
                                                DataByProj[ProjId][key] = 1
                                            
                                            """ Leave out file size for secondary
                                            DICOMs for now: 
                                            key = f'Size of {Modality} scans [MB]'
                                            if key in DataByProj[ProjId].keys():
                                                DataByProj[ProjId][key] += FileSize
                                            else:
                                                DataByProj[ProjId][key] = FileSize
                                            """
             
                #key = 'No. of experiments'
                ##if key in DataByProj[ProjId].keys():
                ##    DataByProj[ProjId][key] += 1
                ##else:
                ##    DataByProj[ProjId][key] = 1
                #DataByProj[ProjId][key] += 1
            
            
            if LogToConsole:
                print('----------------------------------------------------\n')
        
        """ Get the average number of scans by project: """
        #AveScans = round(NumScansThisProj/NumSessionsThisProj, 1)
        AveImScans = round(NumImScansThisProj/NumSessionsThisProj, 1)
        
        #DataByProj[ProjId]['Ave. scans/session'] = AveScans
        DataByProj[ProjId]['Ave. image scans/session'] = AveImScans
        
        if LogToConsole:
            print('****************************************************\n')
    
    
    """ Round the file sizes and file counts expressed in k to 2 decimals: """
    for project in list(DataByProj.keys()):
        for key, value in DataByProj[project].items():
            if '[MB]' in key or '[k]' in key:
                DataByProj[project][key] = round(value, 2)
        
    return DataByProj






def GetInsertDatesByProject(XnatUrl, XnatSession):
    """
    Return a dictionary containing a list of insert dates for all image session
    uploads organised by project.

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.

    Returns
    -------
    DatesByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    request = XnatSession.get(f'{XnatUrl}/data/experiments')
    
    experiments = request.json()['ResultSet']['Result']
    
    DatesByProject = {}
    
    for experiment in experiments:
        expType = experiment['xsiType']
        project = experiment['project']
        #expLabel = experiment['label']
        expId = experiment['ID']
        #expDate = experiment['date']
        insertDate = experiment['insert_date']
    
        if False:
            print(f' expType = {expType}')
            print(f' project = {project}')
            #print(f' expLabel = {expLabel}')
            print(f' expId = {expId}')
            #print(f' expDate = {expDate}')
            print(f' insertDate = {insertDate}')
            print('')
    
        if 'SessionData' in expType:
            request = XnatSession.get(f'{XnatUrl}/data/experiments/{expId}/scans')
    
            #print(request, '\n')
            
            if False:
                if project in list(DatesByProject.keys()):
                    if 'InsertDates' in list(DatesByProject[project].keys()):
                        DatesByProject[project]['InsertDates'].append(insertDate)
    
                    else:
                        #DatesByProject[project].update({'InsertDates' : [insertDate]}) # equivalent to below
                        DatesByProject[project]['InsertDates'] = [insertDate]
                else:
                    DatesByProject[project] = {'InsertDates' : [insertDate]}
            
            if project in list(DatesByProject.keys()):
                DatesByProject[project].append(insertDate)
            else:
                DatesByProject[project] = [insertDate]
    
    return DatesByProject





def GetFirstAndLastImSessionUploadsByProject(XnatUrl, XnatSession, 
                                             DataByProj=None):
    """
    Return a dictionary containing the first and last dates of image session
    uploads organised by project. If passed a dictionary (with projects as 
    keys), the upload details will be added to existing data (e.g. 
    investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    #from datetime import datetime
    
    # Get insert dates by project:
    DatesByProject = GetInsertDatesByProject(XnatUrl, XnatSession)
    
    #print(f'DatesByProject = {DatesByProject}')
    
    if DataByProj == None:
        DataByProj = {}
    
    for project, InsertDates in DatesByProject.items():
        #print(f'There are {len(InsertDates)} insert dates for project {project}')
        
        """
        if len(InsertDates) == 1:
            firstDate = InsertDates[0]
            lastDate = InsertDates[0]
        else:
            # Arbitrarily set initial values:
            firstDate = InsertDates[0]
            lastDate = InsertDates[0]
            
            for date in InsertDates:
                if date < firstDate:
                    firstDate = date
                if date > lastDate:
                    lastDate = date
        """
        
        # Arbitrarily set initial values (or only value if list of length 1):
        firstDate = InsertDates[0]
        lastDate = InsertDates[0]
        
        # Update firstDate and lastDate:
        for date in InsertDates:
            if date < firstDate:
                firstDate = date
            if date > lastDate:
                lastDate = date
            
        if project in list(DataByProj.keys()):
            if 'First upload' in list(DataByProj[project].keys()):
                DataByProj[project]['First upload'] = firstDate
            else:
                #DataByProj[project].update({'First upload' : firstDate}) # equivalent to below
                DataByProj[project]['First upload'] = firstDate
            
            if 'Last upload' in list(DataByProj[project].keys()):
                DataByProj[project]['Last upload'] = lastDate
            else:
                #DataByProj[project].update({'Last upload' : lastDate}) # equivalent to below
                DataByProj[project]['Last upload'] = lastDate
        else:
            DataByProj[project] = {'First upload' : firstDate}
            #DataByProj.update({project : {'First upload' : firstDate}})
            #DataByProj[project] = {'Last upload' : lastDate}
            DataByProj[project].update({'Last upload' : lastDate})
        
                    
    return DataByProj 





def GetProjectDescriptionByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing the project description organised by 
    project. If passed a dictionary (with projects as keys), the description
    will be added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    if DataByProj == None:
        DataByProj = {}
    
    # Get a list of project IDs:
    ProjIds = GetProjectIds(XnatUrl, XnatSession)
    
    for projId in ProjIds:
        project = XnatSession.get(f'{XnatUrl}/data/projects/{projId}?format=json').json()
        
        try:
            # Try to get the description (if it exists):
            desc = project['items'][0]['data_fields']['description']
        except KeyError:
            # There is no description for this project:
            desc = 'None'
        
        DataByProj[projId]['Project description'] = desc
    
    return DataByProj






def GetStartDateByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing the date when project metadata was added to
    XNAT, organised by project. If passed a dictionary (with projects as keys),
    the date will be added to existing data (e.g. investigators organised by 
    project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    
    Notes
    -----
    It seems that start date within children is when the investigator details
    were addded to XNAT, and the start date within meta is when project-related
    meta data was added (e.g. project description).  So neither dates have
    anything to do with first upload.
    """
    from datetime import datetime
    
    if DataByProj == None:
        DataByProj = {}
    
    # Get a list of project IDs:
    ProjIds = GetProjectIds(XnatUrl, XnatSession)
    
    for projId in ProjIds:
        project = XnatSession.get(f'{XnatUrl}/data/projects/{projId}?format=json').json()
        
        date = project['items'][0]['meta']['start_date']
    
        # Convert from ctime (e.g. Fri May 21 22:14:17 UTC 2021) to datetime 
        # (e.g. 2021-05-21 22:14:17):
        # Note: Both below give same results:
        #date = = datetime.strftime(date, '%Y-%m-%d %H:%M:%S')
        date = str(datetime.strptime(date, '%a %b %d %H:%M:%S UTC %Y'))
        
        DataByProj[projId]['Start date'] = date
    
    return DataByProj






def GetUsersByProject(XnatUrl, XnatSession, DataByProj=None):
    """
    Return a dictionary containing the project users organised by project. If 
    passed a dictionary (with projects as keys), the date will be added to 
    existing data (e.g. investigators organised by project).

    Parameters
    ----------
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    XnatSession : Requests Session
        A valid Requests Session for the XNAT of interest.
    DataByProj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    if DataByProj == None:
        DataByProj = {}
    
    # Get a list of project IDs:
    ProjIds = GetProjectIds(XnatUrl, XnatSession)
    
    for projId in ProjIds:
        request = XnatSession.get(f'{XnatUrl}/data/projects/{projId}/users')
        
        ProjectUsers = request.json()['ResultSet']['Result']
        
        NumOfUsers = len(ProjectUsers)
        
        Users = []
        
        for i in range(NumOfUsers):
            fname = ProjectUsers[i]['firstname']
            lname = ProjectUsers[i]['lastname']
            
            Users.append(f"{fname} {lname}")
        
        if projId in list(DataByProj.keys()):
            DataByProj[projId]['Users'] = Users
            DataByProj[projId]['No. of users'] = NumOfUsers
        else:
            DataByProj[projId] = {'Users' : Users}
            DataByProj[projId] = {'No. of users' : NumOfUsers}
    
    return DataByProj







def ReorderKeysAndFillZeros(DataByProj):
    """
    Reorder the keys in the dictionary DataByProj and zero-fill missing keys
    (i.e. if there are no PET scans for a particular project, there will be no
     'PET' key, so add it and assign the number 0).
    
    Inputs:
    ******
    
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Outputs:
    *******
    
    NewDataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    
    Notes:
    *****
    
    Ensures that PI comes before CIs.
    Orders image sessions as follows: MR, CT, PET, SC, RTIMAGE.
    
    It isn't clear what are the image session types available within XNAT. At
    present that info is not available, e.g. see:
    https://wiki.xnat.org/documentation/xnat-administration/managing-data-types-in-xnat/adding-and-managing-image-scan-data-types
    
    so will define a list containing keys in the desired order. Any items not
    included in the list will follow whatever order they were included in the
    original dictionary.
    """
    
    newOrder = ['PI', 'Project description', 'Co-investigators', 
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
    
    #print(f'newOrder = {newOrder}\n')
    
    NewDataByProj = {}
    
    for project in DataByProj.keys():
        NewDataByProj[project] = {}
        
        oldOrder = list(DataByProj[project].keys())
        
        #print(f'oldOrder = {oldOrder}\n')
        
        # First add the keys in the order they appear in newOrder:
        for key in newOrder:
            if key in oldOrder:
                NewDataByProj[project][key] = DataByProj[project][key]
            else:
                #print(f"key '{key}' not in newOrder\n")
                NewDataByProj[project][key] = 0
                #NewDataByProj[project].update({key : 0})
        
        # Add any keys in oldOrder that were not in newOrder:
        for key in oldOrder:
            if not key in newOrder:
                NewDataByProj[project][key] = DataByProj[project][key]
    
    return NewDataByProj






def GetXnatSnapshot(XnatUrl, Username=None, Password=None, XnatSession=None,
                    ExportToXlsx=False, LogToConsole=False):
    """
    Get a "snapshot" of info from an XNAT broken down by projects, and 
    XNAT-wide.
    
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
    
    XnatSession : requests session (optional; None by default)
        If provided a multiple XNAT session requests will be avoided.
    
    ExportToXlsx : boolean (optional; False by default)
    
    Outputs:
    *******
    
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    DataXnatWide : dictionary
        A dictionary containing XNAT-wide analytics.
    
    XnatSession : requests session
    
    
    Notes:
    *****
    
    List of items of interest (message from Simon on 25/05/21):

    Project name
    PI and co-investigators
    Number of subjects
    Number of imaging sessions
    Number of imaging sessions broken down by type (MR, PET, CT, etc.)
    First data upload
    Last data upload/activity
    Number of authorised users
    User list
    Project description
    Average number of scans per session - needs calculating
    Total data storage for project - needs calculating
    
    If ExportToXlsx = True, the data dictionary will be exported to the
    current working directory with an automatically generated filename of the
    form "YYYYMMDD_XNAT_snapshot.xlsx".
    
    Times for my XNAT:
        *Took 0.09 s to fetch investigators.
        *Took 0.01 s to fetch subjects.
        *Took 2.78 s to fetch file counts and sizes.
        *Took 0.04 s to fetch project descriptions.
        *Took 0.04 s to fetch users.
        *Took 0.29 s to fetch image session uploads.
        *Took 0.0 s to reorder keys.
        *Took 3.26 s to generate data-by-project snapshot.
        *Took 0.0 s to generate XNAT-wide snapshot.
        *Took 3.26 s to get all snapshots.
    """
    
    import time
    #from getpass import getpass
    #import requests
    #import os
    #from pathlib import Path
    #import datetime
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    #""" Get a list of projects: """
    #request = XnatSession.get(f'{XnatUrl}/data/projects')
    #projects = request.json()
    #Projects = [item['name'] for item in projects['ResultSet']['Result']]
    
    
    #""" Get a list of PIs: """
    #PIs = [f"{item['pi_firstname']} {item['pi_lastname']}" for item in projects['ResultSet']['Result']]
    
    
    
    """ Get a list of investigators grouped by project: """
    DataByProj = GetInvestigatorsByProject(XnatUrl, XnatSession)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to fetch investigators.\n')
        
    """ Get a list of subjects grouped by project: """
    DataByProj = GetSubjectsByProject(XnatUrl, XnatSession, DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to fetch subjects.\n')
    
    """ Get the number of experiments by project: """
    #DataByProj = GetNumOfExpsByProject(XnatUrl, XnatSession, DataByProj)
    
    
    """ Get the number of scans by image session type by project: """
    #DataByProj = GetNumOfImSessionByTypeByProject(XnatUrl, XnatSession, 
    #                                              DataByProj)
    
    """ Get the number of experiments, image sessions by modality, image scans 
    by modality, and average number of scans per session, all organised by 
    project: """
    #DataByProj = GetExpsSessionsScansByProject(XnatUrl, XnatSession, DataByProj)
    
    """ Get the no. of experiments, no. of image sessions by modality, no. of
    image scans by modality, ave. scans/session, ave. image scans/session,
    no. of DICOM/RTSTRUCT/AIM/SEG files, total size of DICOM/RTSTRUCT/SEG files,
    etc, all organised by project: """
    DataByProj = GetFileCountAndSizeByTypeByProject(XnatUrl, XnatSession,
                                                    DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to fetch file counts and sizes.\n')
    
    """ Get the project description by project: """
    DataByProj = GetProjectDescriptionByProject(XnatUrl, XnatSession, 
                                                DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to fetch project descriptions.\n')
        
    """ Start date might not be the same as first upload.. """
    #""" Get the start date by project: """
    #DataByProj = GetStartDateByProject(XnatUrl, XnatSession, DataByProj)
    
    
    """ Get the project users by project: """
    DataByProj = GetUsersByProject(XnatUrl, XnatSession, DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to fetch users.\n')
        
    """ Get the first and last upload dates by project: """
    DataByProj = GetFirstAndLastImSessionUploadsByProject(XnatUrl, XnatSession, 
                                                          DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to fetch image session uploads.\n')
        
    """ Re-order the keys at the second tier: """
    DataByProj = ReorderKeysAndFillZeros(DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to reorder keys.\n')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 2)
    print(f'*Took {Dtime} s to generate data-by-project snapshot.\n')
    
    """ Get XNAT-wide analytics (not broken down by project, i.e.:
        - no. of unique users
        - average no. of users per project
    """
    DataXnatwide = GetXnatWideSnapshot(DataByProj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if LogToConsole:
        print(f'*Took {Dtime} s to generate XNAT-wide snapshot.\n')
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 2)
    print(f'*Took {Dtime} s to get all snapshots.\n')
    
    
    """ Export the dictionary: """
    if ExportToXlsx:
        import datetime
        from importlib import reload
        import GeneralTools
        reload(GeneralTools)
        from GeneralTools import ExportDictionaryToXlsx
        from GeneralTools import ExportDictionaryOfDictionariesToXlsx
        
        TimeNow = datetime.datetime.now()
        DateTime = TimeNow.strftime('%Y%m%d_%H%M%S')
        
        FileName = DateTime + '_XNAT_by_project_snapshot.xlsx'
        
        # Export the by-project snapshot:
        ExportDictionaryOfDictionariesToXlsx(DataByProj, FileName, 
                                             ExportDir='cwd',
                                             NameOfFirstColumn='Project')
        
        FileName = DateTime + '_XNAT_wide_snapshot.xlsx'
        
        # Export the XNAT-wide snapshot:
        ExportDictionaryToXlsx(DataXnatwide, FileName, ExportDir='cwd')
        
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 2)
        print(f'*Took {Dtime} s to export the snapshots to XLSX files.\n')
    
    return DataByProj, DataXnatwide, XnatSession




def GetXnatWideSnapshot(DataByProj):
    """ 
    Get XNAT-wide analytics (not broken down by project).
    
    Inputs:
    ******
    
    DataByProj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Outputs:
    *******
    
    DataXnatWide : dictionary
        A dictionary containing XNAT-wide analytics.
    """
    
    # Number of decimals for rounding:
    #Decimals = 2
    
    NumOfProjects = len(list(DataByProj.keys()))
    
    # Initialise:
    UniqueUsers = []
    NumUsers = 0
    #NumUsersByProj = []
    #NumSubjectsByProj = []
    NumSubjects = 0
    NumExps = 0
    NumMrSessions = 0
    NumCtSessions = 0
    NumPetSessions = 0
    NumMrScans = 0
    NumCtScans = 0
    NumPetScans = 0
    NumOtScans = 0
    NumDicomImFiles = 0
    SizeDicomImFiles = 0
    
    
    for ProjId in DataByProj.keys():
        users = DataByProj[ProjId]['Users']
        
        NumUsers += len(users)
        
        #NumUsersByProj.append(len(users))
        
        for user in users:
            if not user in UniqueUsers:
                UniqueUsers.append(user)
        
        NumSubjects += DataByProj[ProjId]['No. of subjects']
        NumExps += DataByProj[ProjId]['No. of experiments']
        NumMrSessions += DataByProj[ProjId]['No. of MR sessions']
        NumCtSessions += DataByProj[ProjId]['No. of CT sessions']
        NumPetSessions += DataByProj[ProjId]['No. of PET sessions']
        NumMrScans += DataByProj[ProjId]['No. of MR scans']
        NumCtScans += DataByProj[ProjId]['No. of CT scans']
        NumPetScans += DataByProj[ProjId]['No. of PET scans']
        NumOtScans += DataByProj[ProjId]['No. of OT scans']
        NumDicomImFiles += DataByProj[ProjId]['No. of DICOM image files [k]']
        SizeDicomImFiles += DataByProj[ProjId]['Size of DICOM image files [MB]']
        
        #NumSubjectsByProj.append(DataByProj[ProjId]['No. of subjects'])
    
    AveNumUsersPerProj = round(NumUsers/NumOfProjects, 1)
    #AveNumUsersPerProj = round(NumUsers/NumOfProjects, 1)
    AveNumSubjectsPerProj = round(NumSubjects/NumOfProjects, 1)
    AveNumExpsPerProj = round(NumExps/NumOfProjects, 1)
    AveNumMrSessionsPerProj = round(NumMrSessions/NumOfProjects, 1)
    AveNumCtSessionsPerProj = round(NumCtSessions/NumOfProjects, 1)
    AveNumPetSessionsPerProj = round(NumPetSessions/NumOfProjects, 1)
    AveNumMrScansPerProj = round(NumMrScans/NumOfProjects, 1)
    AveNumCtScansPerProj = round(NumCtScans/NumOfProjects, 1)
    AveNumPetScansPerProj = round(NumPetScans/NumOfProjects, 1)
    AveNumOtScansPerProj = round(NumOtScans/NumOfProjects, 1)
    AveNumDicomImFilesPerProj = round(NumDicomImFiles/NumOfProjects, 2)
    AveSizeDicomImFilesPerProj = round(SizeDicomImFiles/NumOfProjects, 2)
    
    
    DataXnatWide = {'No. of unique users' : len(UniqueUsers),
                    'Ave. users/project' : AveNumUsersPerProj,
                    'No. of projects' : NumOfProjects,
                    'No. of subjects' : NumSubjects,
                    'Ave. subjects/project' : AveNumSubjectsPerProj,
                    'No. of experiments' : NumExps,
                    'Ave. experiments/project' : AveNumExpsPerProj,
                    'No. of MR sessions' : NumMrSessions,
                    'Ave. MR sessions/project' : AveNumMrSessionsPerProj,
                    'No. of CT sessions' : NumCtSessions,
                    'Ave. CT sessions/project' : AveNumCtSessionsPerProj,
                    'No. of PET sessions' : NumPetSessions,
                    'Ave. PET sessions/project' : AveNumPetSessionsPerProj,
                    'No. of MR scans' : NumMrScans,
                    'Ave. MR scans/project' : AveNumMrScansPerProj,
                    'No. of CT scans' : NumCtScans,
                    'Ave. CT scans/project' : AveNumCtScansPerProj,
                    'No. of PET scans' : NumPetScans,
                    'Ave. PET scans/project' : AveNumPetScansPerProj,
                    'No. of OT scans' : NumOtScans,
                    'Ave. OT scans/project' : AveNumOtScansPerProj,
                    'No. of DICOM image files [M]' : round(NumDicomImFiles/1000, 2),
                    'Ave. DICOM image files/project [k]' : AveNumDicomImFilesPerProj,
                    'Size of DICOM image files [GB]' : round(SizeDicomImFiles/1000, 2),
                    'Ave. size of DICOM image files/project [MB]' : AveSizeDicomImFilesPerProj
                    }
    
    return DataXnatWide




def UploadXnatSnapshot_WIP(XnatUrl, Username=None, Password=None, XnatSession=None):
    """
    need to edit below...
    
    
    Get a "snapshot" of info from an XNAT.
    
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
    
    
    Notes:
    *****
    
    List of items of interest (message from Simon on 25/05/21):

    Project name
    PI and co-investigators
    Number of subjects
    Number of imaging sessions
    Number of imaging sessions broken down by type (MR, PET, CT, etc.)
    First data upload
    Last data upload/activity
    Number of authorised users
    User list
    Project description
    Average number of scans per session - needs calculating
    Total data storage for project - needs calculating
    """
    
    #from getpass import getpass
    #import requests
    import os
    from pathlib import Path
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    #""" Get a list of projects: """
    #request = XnatSession.get(f'{XnatUrl}/data/projects')
    #projects = request.json()
    #Projects = [item['name'] for item in projects['ResultSet']['Result']]
    
    
    #""" Get a list of PIs: """
    #PIs = [f"{item['pi_firstname']} {item['pi_lastname']}" for item in projects['ResultSet']['Result']]
    
    
    
    """ Get a list of investigators grouped by project: """
    DataByProj = GetInvestigatorsByProject(XnatUrl, XnatSession)
    
    
    """ Get a list of subjects grouped by project: """
    DataByProj = GetSubjectsByProject(XnatUrl, XnatSession, DataByProj)
    
    
    """ Get the number of experiments by project: """
    DataByProj = GetNumOfExpsByProject(XnatUrl, XnatSession, DataByProj)
    
    
    """ Get the number of scans by image session type by project: """
    DataByProj = GetNumOfImSessionByTypeByProject(XnatUrl, XnatSession, 
                                                  DataByProj)
    
    
    return DataByProj




