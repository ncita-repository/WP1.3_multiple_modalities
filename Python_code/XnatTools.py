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


def DownloadXnatData_OLD(XnatUrl, Username=None, Password=None,
                     ProjectLabel=None, SubjLabel=None, 
                     ExperimentLabel=None, ExportDir=None):
    """
    Download data from XNAT.
    
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
    
    ProjectLabel : string (optional; None by default)
        The project label of interest. The label doesn't need to be exact but 
        can be a sub-set of the strings contained in the label (i.e. a key 
        word).  If not provided (i.e. ProjectLabel = None), all projects will 
        be considered of interest.
    
    SubjLabel : string (optional; None by default)
        The subject label of interest. The label doesn't need to be exact but 
        can be a sub-set of the strings contained in the label (i.e. a key 
        word).  If not provided (i.e. SubjLabel = None), all subjects will 
        be considered of interest within ProjectLabel.
    
    ExperimentLabel : string (optional; None by default)
        The experiment label of interest. The label doesn't need to be exact 
        but can be a sub-set of the strings contained in the label (i.e. a key 
        word).  If not provided (i.e. ExperimentLabel = None), all sessions 
        will be considered of interest within SubjLabel.
    
    ExportDir : string (optional; None by default)
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, SubjLabel, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    """
    
    import os
    from pathlib import Path
    import pyxnat
    from getpass import getpass
    
    if Username == None:
        Username = getpass("Enter user name: ")
    
    if Password == None:
        Password = getpass("Enter password:")
    
    if ExportDir == None:
        ExportDir = os.path.join(Path.home(), "Downloads", "XnatDownloads")
    
    
    """ Store info in a dictionary. """
    XnatDict = {}
    
    interface = pyxnat.Interface(server=XnatUrl, user=Username, 
                                 password=Password)
    
    
    if ProjectLabel:
        ProjectLabels = interface.select.projects().get()
        
        ind = None
        
        for i in range(len(ProjectLabels)):
            if ProjectLabel in ProjectLabels[i]:
                ind = i
                
                Projects = interface.select.projects(ProjectLabels[i])
        
        if ind == None:
            msg = 'There are no projects that contain the keyword'\
                  f'{ProjectLabel}'
            raise Exception(msg)
    
    else:
        Projects = interface.select.projects()
    
    
    
    for project in Projects:
        
        ProjLabel = project.label()
        
        print(f'Project label = {ProjLabel}')
        
        XnatDict.update({ProjLabel : {}})
        
        if SubjLabel:
            SubjLabels = project.subjects().get()
            
            ind = None
            
            for i in range(len(SubjLabels)):
                if SubjLabel in SubjLabels[i]:
                    ind = i
                    
                    Subjects = project.subjects(SubjLabels[i])
            
            if ind == None:
                msg = 'There are no subjects that contain the keyword'\
                      f'{SubjLabel}'
                raise Exception(msg)
            
            
        else:
            Subjects = project.subjects()
    
    
        
        for subject in Subjects:
            
            SubjLabel = subject.label()
            
            print(f'Subject label = {SubjLabel}')
            
            XnatDict[ProjLabel].update({SubjLabel : {}})
            
            if ExperimentLabel:
                ExperimentLabels = project.subjects().experiments().get()
                
                ind = None
                
                for i in range(len(ExperimentLabels)):
                    if ExperimentLabel in ExperimentLabels[i]:
                        ind = i
                        
                        Experiments = subject.experiments(ExperimentLabels[i])
                
                if ind == None:
                    msg = 'There are no experiments that contain the keyword'\
                          f'{ExperimentLabel}'
                    raise Exception(msg)
                
                
            else:
                Experiments = subject.experiments()
    
    
            
            for experiment in Experiments:
                
                ExpLabel = experiment.label()
                
                print(f'Experiment label = {ExpLabel}')
                
                XnatDict[ProjLabel][SubjLabel].update({ExpLabel : {}})
                
                Scans = experiment.scans()
                
                DownloadDir = os.path.join(ExportDir, ProjLabel, SubjLabel)
                
                print(f'\nDownloadDir = {DownloadDir}')
                
                if not os.path.isdir(DownloadDir):
                    print(False)
                    #os.mkdir(DownloadDir)
                    Path(DownloadDir).mkdir(parents=True)
                    #Path(DownloadDir).mkdir(parents=True, exist_ok=True)
                    
                """ Download all scans for this experiment: """
                #Scans.download(DownloadDir, type='ALL', extract=True)
                
                
                ScanLabels= []
                ScanUids = []
                
                for scan in Scans:
                    #ScanDict = {}
                    
                    #Files = scan.resources().files()
                    ScanLabel = scan.label()
                    ScanLabels.append(ScanLabel)
                    #ScanDesc = scan.attrs.get('type')
                    #ScanUid = scan.attrs.get('UID')
                    ScanDesc_Uid = scan.attrs.mget(['type', 'UID'])
                    ScanDesc = ScanDesc_Uid[0]
                    ScanUid = ScanDesc_Uid[1]
                    
                    ScanUids.append(ScanUid)
                    
                    ScanLabelDesc = f'{ScanLabel}-{ScanDesc}'
                    
                    ScanDir = os.path.join(DownloadDir, ExpLabel, 'scans',
                                           ScanLabelDesc, 'resources', 'DICOM',
                                           'files')
                    
                    
                    #XnatDict[ProjLabel][SubjLabel][ExpLabel].update({'ScanLabel' : ScanLabel,
                    #                                                 'ScanDesc' : ScanDesc,
                    #                                                 'ScanDirName' : ScanLabelDesc,
                    #                                                 'ScanDir' : ScanDir,
                    #                                                 'ScanUid' : ScanUid})
                    
                    ScanDict = {'ScanDesc' : ScanDesc,
                                'ScanDirName' : ScanLabelDesc,
                                'ScanDir' : ScanDir,
                                'ScanUid' : ScanUid}
                    
                    XnatDict[ProjLabel][SubjLabel][ExpLabel].update({ScanLabel : ScanDict})
                    
                    
                    
                
                
                
                Assessors = experiment.assessors().resources().files()
                
                for assessor in Assessors:
                    print(str(assessor.get()))
                    
                    #AsrLabel = assessor.label()
                    #AsrLabel = assessor.attrs.get('label')
                    #AsrDate = assessor.attrs.get('date')
                    AsrLabel_Date_Time = assessor.attrs.mget(['label', 'date', 'time'])
                    AsrLabel = AsrLabel_Date_Time[0]
                    AsrDate = AsrLabel_Date_Time[1]
                    AsrTime = AsrLabel_Date_Time[2]
                    
                    #content = str(assessor.get())
                    """ Download the assessor: """
                    content = assessor.get()
                    #content = content.content()
                    
                    print(f"\n{content}")
                    #print(f"\n{content.split('icr:name>')}")
                    
                    AsrName = content.split('icr:name>')[1].split('<')[0]
                    AsrMod = content.split('xnat:file label="')[1].split('"')[0]
                    ScanUid = content.split('seriesUID>')[1].split('<')[0]
                    AsrDate = content.split('xnat:date>')[1].split('<')[0]
                    AsrTime = content.split('xnat:time>')[1].split('<')[0]
                    
                    """ Index of the scan with this ScanUid: """
                    ind = ScanUids.index(ScanUid)
                    
                    ScanLabelDesc = XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['ScanDirName']
                    
                    AssessorDict = {'AsrLabel' : AsrLabel,
                                    'AsrDate' : AsrDate,
                                    'AsrTime' : AsrTime,
                                    #'AsrName' : AsrName,
                                    'AsrMod' : AsrMod,
                                    'ScanUid' : ScanUid,
                                    'ScanLabel' : ScanLabels[ind]}
                    
                    keys = XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel].keys()
                    
                    if 'Assessors' in keys:
                        keys = XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'].keys()
                        
                        if AsrMod in keys:
                            XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'][AsrMod].update({'AsrName' : AsrName})
                            #XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'][AsrMod][AsrName] = AssessorDict
                        else:
                            XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'].update({'AsrMod' : AsrMod})
                            XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'][AsrMod].update({'AsrName' : AsrName})
                    else:
                        XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel].update({'Assessors' : {}})
                        XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'].update({'AsrMod' : AsrMod})
                        XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'][AsrMod].update({'AsrName' : AsrName})
                    
                    XnatDict[ProjLabel][SubjLabel][ExpLabel][ScanLabel]['Assessors'][AsrMod][AsrName] = AssessorDict
                    
                    AssessorDir = os.path.join(DownloadDir, ExpLabel, 'scans',
                                               ScanLabelDesc, 'resources', 
                                               AsrMod, AsrName)
                
                    if not os.path.isdir(AssessorDir):
                        os.mkdir(AssessorDir)
            
            
            
                
        
    return XnatDict




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





def DownloadScan(XnatUrl, ProjId, SubjLabel, StudyLabel, SeriesLabel,
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
    
    StudyLabel : string
        The study / experiment label of interest. 
    
    SeriesLabel : string
        The series label / scan ID of interest.
    
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
    #import requests
    import zipfile
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatDownloads")
    
    
    
    
    
    Scan = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                           f'subjects/{SubjLabel}/'
                           f'experiments/{StudyLabel}/'
                           f'scans/{SeriesLabel}/'
                           f'resources/DICOM/files?format=zip')
    
    if not '200' in str(Scan):
        if '403' in str(Scan):
            msg = f'Response 403: Forbidden'
        elif '404' in str(Scan):
            msg = f'Response 403: Not found'
        elif '409' in str(Scan):
            msg = f'Response 409: Conflict'
        elif '422' in str(Scan):
            msg = f'Response 422: Unprocessable entity'
        else:
            msg = str(Scan)
            
        raise Exception(msg)
    
    
    
    ExportDir = os.path.join(ExportRootDir, ProjId, SubjLabel)

    if not os.path.isdir(ExportDir):
        print('Creating directory', ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    Fpath = os.path.join(ExportDir, 
                         f'Experiment_{StudyLabel}__Scan_{SeriesLabel}.zip')
    
    if not os.path.exists(Fpath):
    
        with open(Fpath, 'wb') as file:
            file.write(Scan.content)
            
        print(f'Zipped file downloaded to: \n{Fpath}\n')
    
    
        zip_file = zipfile.ZipFile(Fpath, mode='r')
    
        zip_file.extractall(ExportDir)
        
        print(f'Zipped file extracted to: \n{ExportDir}\n')
    else:
        print(f'Zipped file already downloaded to: \n{Fpath}\n')
    
    
    
    """ Store info in a dictionary: """
    
    if PathsDict == None:
        PathsDict = {}
        keys = []
    else:
        keys = list(PathsDict.keys())
    
    
    if not ProjId in keys:
        PathsDict.update({ProjId : {}})
    
    keys = list(PathsDict[ProjId].keys())
    
    if not SubjLabel in keys:
        PathsDict[ProjId].update({SubjLabel : {}})
    
    keys = list(PathsDict[ProjId][SubjLabel].keys())
    
    if not StudyLabel in keys:
        PathsDict[ProjId][SubjLabel].update({StudyLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel].keys())
    
    if not SeriesLabel in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel]\
        .update({SeriesLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
                .keys())
    
    if not 'DicomDir' in keys:
        """ Get the experiment of interest: """
        request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                                  f'subjects/{SubjLabel}/'
                                  f'experiments/{StudyLabel}?format=json')
        
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
        
        #print(f"{Experiment['items']}")
        #print(f"{Experiment['items'][0]}")
        #print(f"{Experiment['items'][0]['children']}")
        #print(f"{Experiment['items'][0]['children'][1]}")
        #print(f"{Experiment['items'][0]['children'][1]['items']}")
        #print(f"{Experiment['items'][0]['children'][1]['items'][0]}")
        #print(f"{Experiment['items'][0]['children'][1]['items'][0]['data_fields']}")
        
        #StudyUid = Experiment['items'][0]['children'][1]['items'][0]['data_fields']['type']['UID']
        StudyUid = Experiment['items'][0]['children'][1]['items'][0]['data_fields']['UID']
        
        """ Get the list of series labels (scan IDs), SeriesUIDs and series 
        descriptions for this series/scan: """
        SeriesLabels = []
        SeriesUids = []
        SeriesDescs = []
        
        for experiment in Experiment['items'][0]['children'][1]['items']:
            SeriesLabels.append(experiment['data_fields']['ID'])
            SeriesUids.append(experiment['data_fields']['UID'])
            SeriesDescs.append(experiment['data_fields']['type'])
            
        
        SeriesUid = SeriesUids[SeriesLabels.index(SeriesLabel)]
        SeriesDesc = SeriesDescs[SeriesLabels.index(SeriesLabel)]
        
        """ The folder name between the 'scans' folder and 'resources' folder
        has the form 'SeriesLabel-SeriesDesc', but it appears that spaces and 
        full stops within SeriesDesc become underscores in the folder structure 
        of the extracted zip file. """
        DirName = f'{SeriesLabel}-{SeriesDesc}'.replace(' ', '_').replace('.', '_')
        
        #ResourceDir = os.path.join(ExportRootDir, StudyLabel, 'scans', DirName,
        #                           'resources')
        
        ResourceDir = os.path.join(ExportRootDir, ProjId, SubjLabel, StudyLabel, 
                                   'scans', DirName, 'resources')
        
        DicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
        .update({'StudyUID' : StudyUid,
                 'SeriesDesc' : SeriesDesc,
                 'SeriesUID' : SeriesUid,
                 'DicomDir' : DicomDir})
    
    
    return PathsDict, XnatSession





def AddScanPathToDict(XnatUrl, ProjId, SubjLabel, StudyLabel, 
                      SeriesLabel, XnatSession=None, ExportRootDir='default', 
                      PathsDict=None, Username=None, Password=None):
    """
    Add the directory for a downloaded scan from XNAT and the SeriesUID of the
    scan to a dictionary.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    StudyLabel : string
        The study / experiment label of interest. 
    
    SeriesLabel : string
        The series label / scan ID of interest.
    
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
        Dictionary containing directory and file paths to data downloaded from
        XNAT.
    """
    
    import os
    from pathlib import Path
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatDownloads")
    
    
    if PathsDict == None:
        PathsDict = {}
        keys = []
    else:
        keys = list(PathsDict.keys())
    
    
    if not ProjId in keys:
        PathsDict.update({ProjId : {}})
    
    keys = list(PathsDict[ProjId].keys())
    
    if not SubjLabel in keys:
        PathsDict[ProjId].update({SubjLabel : {}})
    
    keys = list(PathsDict[ProjId][SubjLabel].keys())
    
    if not StudyLabel in keys:
        PathsDict[ProjId][SubjLabel].update({StudyLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel].keys())
    
    if not SeriesLabel in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel].update({SeriesLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel].keys())
    
    if not 'DicomDir' in keys:
        """
        # Get the type value (i.e. description) for the chosen scan:
        request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                                  f'subjects/{SubjLabel}/'
                                  f'experiments/{StudyLabel}/scans')
        
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
        
        Scans = request.json()
        
        SeriesLabels = []
    
        for scan in Scans['ResultSet']['Result']:
            SeriesLabels.append(scan['ID'])
        
        ind = SeriesLabels.index(SeriesLabel)
            
        ScanDesc = Scans['ResultSet']['Result'][ind]['type']
        """
        
        
        """ Get the experiment of interest: """
        request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                                  f'subjects/{SubjLabel}/'
                                  f'experiments/{StudyLabel}?format=json')
        
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
        
        StudyUid = Experiment['items'][0]['children'][1]['items'][0]\
                   ['data_fields']['type']['UID']
                   
        """ Get the list of scan IDs and SeriesUIDs for this scan: """
        SeriesLabels = []
        SeriesUids = []
        SeriesDescs = []
        
        for experiment in Experiment['items'][0]['children'][1]['items']:
            SeriesLabels.append(experiment['data_fields']['ID'])
            SeriesUids.append(experiment['data_fields']['UID'])
            SeriesDescs.append(experiment['data_fields']['type'])
            
        """ The SeriesUID and series description for the scan of interest 
        (SeriesLabel): """
        SeriesUid = SeriesUids[SeriesLabels.index(SeriesLabel)]
        SeriesDesc = SeriesDescs[SeriesLabels.index(SeriesLabel)]
        
        
        ResourceDir = os.path.join(ExportRootDir, StudyLabel, 'scans',
                                f'{SeriesLabel}-{SeriesDesc}', 'resources')
        
        DicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        
        #AsrDir = os.path.join(ResourceDir, '')
        
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel].update({'StudyUID' : StudyUid,
                                                                            'SeriesDesc' : SeriesDesc,
                                                                            'SeriesUID' : SeriesUid,
                                                                            'DicomDir' : DicomDir})
    
    
    return PathsDict







def DownloadAssessors(XnatUrl, ProjId, SubjLabel, StudyLabel, 
                      XnatSession=None, ExportRootDir='default', 
                      PathsDict=None, Username=None, Password=None):
    """
    Download all assessors for a particular scan from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    StudyLabel : string
        The study / experiment label of interest. 
    
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
    
    XnatSession : requests session
    
    """
    
    import os
    from pathlib import Path
    #import requests
    import zipfile
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == 'default':
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatDownloads")
    
    
    
    Assessors = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                                f'subjects/{SubjLabel}/'
                                f'experiments/{StudyLabel}/'
                                f'assessors/ALL/files?file_format=DICOM&format=zip')
    
    
    
    if not '200' in str(Assessors):
        if '403' in str(Assessors):
            msg = f'Response 403: Forbidden'
        elif '404' in str(Assessors):
            msg = f'Response 403: Not found'
        elif '409' in str(Assessors):
            msg = f'Response 409: Conflict'
        elif '422' in str(Assessors):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    
    
    ExportDir = os.path.join(ExportRootDir, ProjId, SubjLabel, 
                             StudyLabel)

    if not os.path.isdir(ExportDir):
        print('Creating directory', ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    Fpath = os.path.join(ExportDir, 'assessors.zip')
    
    with open(Fpath, 'wb') as file:
        file.write(Assessors.content)
    
    print(f'Zipped file downloaded to: \n{Fpath}\n')
    
    
    zip_file = zipfile.ZipFile(Fpath, mode='r')
    
    ExtractDir = os.path.join(ExportDir, 'assessors')
    
    if not os.path.isdir(ExtractDir):
        print('Creating directory', ExtractDir)
        Path(ExtractDir).mkdir(parents=True)
    
    zip_file.extractall(ExtractDir)
    
    print(f'Zipped file extracted to: \n{ExportDir}\n.')
    
    return XnatSession







def DownloadAssessor(XnatUrl, ProjId, SubjLabel, StudyLabel, SeriesLabel, 
                     AsrMod, AsrName, XnatSession=None, ExportRootDir='default', 
                     PathsDict=None, Username=None, Password=None):
    """
    Download an assessor of interest for a particular scan from XNAT.
    
    Inputs:
    ******
    
    XnatUrl : string
        Address of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    StudyLabel : string
        The study / experiment label of interest. 
    
    SeriesLabel : string
        The scan ID of interest.
    
    AsrMod : string
        The modality of the assessor of interest.
        
    AsrName : string
        The assessor name of interest.
    
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
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatDownloads")
    
    
    
    """ Get the experiment of interest: """
    request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                              f'subjects/{SubjLabel}/'
                              f'experiments/{StudyLabel}?format=json')
    
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
    
    
    """ Get the SeriesUID for the scan of interest: """
    if PathsDict != None:
        try:
            SeriesUid = PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
                        ['SeriesUID']
            
        except Exception:
            """ Get the list of scan IDs and SeriesUIDs for this scan: """
            SeriesLabels = []
            SeriesUids = []
            
            for experiment in Experiment['items'][0]['children'][1]['items']:
                SeriesLabels.append(experiment['data_fields']['ID'])
                SeriesUids.append(experiment['data_fields']['UID'])
                
            """ The SeriesUID and series description for the scan of interest 
            (SeriesLabel): """
            SeriesUid = SeriesUids[SeriesLabels.index(SeriesLabel)]

    
    """ Get the (referenced) SeriesUIDs of all assessors: """
    RefSeriesUids = []

    for experiment in Experiment['items'][0]['children'][0]['items']:
        RefSeriesUids.append(experiment['children'][1]['items'][0]\
                             ['data_fields']['seriesUID'])
    
    """ The indices of the RefSeriesUids that match SeriesUid: """
    inds = [i for i, uid in enumerate(RefSeriesUids) if uid == SeriesUid]
    
    
    if inds == []:
        msg = f"There are no scan assessors for SeriesLabel '{SeriesLabel}' "\
              + f"in StudyLabel '{StudyLabel}' for SubjLabel '{SubjLabel}' "\
              + f"in ProjId '{ProjId}'."
        raise Exception(msg)
    
    
    """ Get the names and other meta data for the assessors for the scan of
    interest: """
    AsrDates = []
    AsrTimes = []
    AsrNames = []
    AsrLabels = []
    AsrMods = []
    AsrIds = []
    
    for ind in inds:
        experiment = Experiment['items'][0]['children'][0]['items'][ind]
        
        AsrDates.append(experiment['data_fields']['date'])
        AsrTimes.append(experiment['data_fields']['time'])
        AsrNames.append(experiment['data_fields']['name'])
        AsrLabels.append(experiment['data_fields']['label'])
        AsrMods.append(experiment['data_fields']['collectionType'])
        AsrIds.append(experiment['data_fields']['id'])
    
    
    """ The indices of assessors that match AsrMod.: """
    inds = [i for i, mod in enumerate(AsrMods) if mod == AsrMod]
    
    """ Use inds to filter the lists. """
    AsrDates = [AsrDates[i] for i in inds]
    AsrTimes = [AsrTimes[i] for i in inds]
    AsrNames = [AsrNames[i] for i in inds]
    AsrLabels = [AsrLabels[i] for i in inds]
    AsrMods = [AsrMods[i] for i in inds]
    AsrIds = [AsrIds[i] for i in inds]
    
    """ The index of the assessor matching AsrName: """
    try:
        ind = AsrNames.index(AsrName)
    
    except Exception:
        msg = f"There are no assessors for Scan '{SeriesLabel}' with name "\
              + f"'{AsrName}' for SeriesLabel '{SeriesLabel}' in StudyLabel "\
              + f"'{StudyLabel}' for SubjLabel '{SubjLabel}' in ProjId "\
              + f"'{ProjId}'."
        raise Exception(msg)
    
    AsrDate = AsrDates[ind]
    AsrTime = AsrTimes[ind]
    AsrLabel = AsrLabels[ind]
    #AsrMod = AsrMods[ind]
    AsrId = AsrIds[ind]
    
    
    """ Get the assessor of interest: """
    
    """ 
    Note:
        
    Although the structure of the REST call below is as in the how-to 
    guide:
    
    https://wiki.xnat.org/display/XAPI/How+To+Download+Files+via+the+XNAT+REST+API
    
    the downloaded assessor was much smaller (hundreds of B rather than tens of
    kB as in the same file downloaded via the XNAT web UI) and was nearly empty
    of any tags.
    
    Assessor = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                               f'subjects/{SubjLabel}/'
                               f'experiments/{StudyLabel}/'
                               f'assessors/{AsrId}/resources/{AsrMod}/files')
    
    Jason suggested adding the file name after 'files', which worked. 
    """
    Fname = AsrLabel + '.dcm'
    
    Assessor = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                               f'subjects/{SubjLabel}/'
                               f'experiments/{StudyLabel}/'
                               f'assessors/{AsrId}/resources/{AsrMod}/'
                               f'files/{Fname}')
    
    if not '200' in str(Assessor):
        if '403' in str(Assessor):
            msg = f'Response 403: Forbidden'
        elif '404' in str(Assessor):
            msg = f'Response 403: Not found'
        elif '409' in str(Assessor):
            msg = f'Response 409: Conflict'
        elif '422' in str(Assessor):
            msg = f'Response 422: Unprocessable entity'
            
        raise Exception(msg)
    
    
    
    ExportDir = os.path.join(ExportRootDir, ProjId, SubjLabel, 
                             StudyLabel, 'assessors')

    if not os.path.isdir(ExportDir):
        print('Creating directory', ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    #Fname = AsrLabel + '.dcm'
    Fpath = os.path.join(ExportDir, Fname)
    
    with open(Fpath, 'wb') as file:
        file.write(Assessor.content)
    
    print(f'Assessor file downloaded to: \n{Fpath}\n')
    
    
    """ Store info in a dictionary: """
    
    if PathsDict == None:
        PathsDict = {}
        keys = []
    else:
        keys = list(PathsDict.keys())
    
    
    if not ProjId in keys:
        PathsDict.update({ProjId : {}})
    
    keys = list(PathsDict[ProjId].keys())
    
    if not SubjLabel in keys:
        PathsDict[ProjId].update({SubjLabel : {}})
    
    keys = list(PathsDict[ProjId][SubjLabel].keys())
    
    if not StudyLabel in keys:
        PathsDict[ProjId][SubjLabel].update({StudyLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel].keys())
    
    if not SeriesLabel in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel].update({SeriesLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel].keys())
    
    if not 'assessors' in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
        .update({'assessors' : {}})
    
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
                ['assessors'].keys())
    
    if not AsrMod in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
        ['assessors'].update({AsrMod : {}})
    
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
                ['assessors'][AsrMod].keys())
    
    if not AsrName in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]\
        ['assessors'][AsrMod].update({AsrName : {}})
    
    PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel]['assessors']\
    [AsrMod][AsrName].update({'AsrDate' : AsrDate,
                               'AsrTime' : AsrTime,
                               'AsrName' : AsrName,
                               'AsrLabel' : AsrLabel,
                               'AsrMod' : AsrMod,
                               'AsrId' : AsrId,
                               'AsrDir' : ExportDir,
                               'AsrFpath' : Fpath})
    
    
    return PathsDict, XnatSession





    """
    # Reduce list to those containing DICOMs:
    FileFormats = []
    
    for assessor in Assessors.json()['ResultSet']['Result']:
        FileFormats.append(assessor['file_format'])
    
    inds = [i for i, FileFormat in enumerate(FileFormats) if FileFormat == 'DICOM']
    
    DicomAssessors = [Assessors.json()['ResultSet']['Result'][i] for i in inds]
    
    AssessorUris = []
    AssessorFnames = []
    Fpaths = []
    
    for assessor in DicomAssessors:
        uri = assessor['URI']
        
        AssessorUris.append(uri)
        
        Fname = assessor['Name']
        
        AssessorFnames.append(Fname)
        
        Fpath = os.path.join(ExportDir, Fname)
        Fpaths.append(Fpath)
    """




def AddAssessorPathToDict(XnatUrl, ProjId, SubjLabel, StudyLabel, 
                          SeriesLabel, AsrName, XnatSession=None, ExportRootDir=None, 
                          PathsDict=None, Username=None, Password=None):
    """
    Add the file path for a downloaded assessor to a dictionary.
    
    Inputs:
    ******
    
    XnatUrl : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    
    ProjId : string
        The project ID of interest.
    
    SubjLabel : string
        The subject label of interest. 
    
    StudyLabel : string
        The study / experiment label of interest. 
    
    SeriesLabel : string
        The series label / scan ID of interest.
    
    AsrName : string
        The ROI name of interest.
    
    XnatSession : requests session (optional; None by default)
        If provided a new session request will be avoided.
    
    ExportRootDir : string (optional; None by default)
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
        Dictionary containing directory and file paths to data downloaded from
        XNAT.
    """
    
    import os
    from pathlib import Path
    
    if XnatSession == None:
        XnatSession = CreateXnatSession(XnatUrl, Username, Password)
    
    if ExportRootDir == None:
        ExportRootDir = os.path.join(Path.home(), "Downloads", "XnatDownloads")
    
    
    
    
    
    
    
    
    if PathsDict == None:
        PathsDict = {}
        keys = []
    else:
        keys = list(PathsDict.keys())
    
    
    if not ProjId in keys:
        PathsDict.update({ProjId : {}})
    
    keys = list(PathsDict[ProjId].keys())
    
    if not SubjLabel in keys:
        PathsDict[ProjId].update({SubjLabel : {}})
    
    keys = list(PathsDict[ProjId][SubjLabel].keys())
    
    if not StudyLabel in keys:
        PathsDict[ProjId][SubjLabel].update({StudyLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel].keys())
    
    if not SeriesLabel in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel].update({SeriesLabel : {}})
        
    keys = list(PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel].keys())
    
    if not 'assessors' in keys:
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel].update({'assessors' : {}})
        
        """ Get the type value (i.e. description) for the chosen scan: """
        
        """ 
        The following should work:
            
        request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                                  f'subjects/{SubjLabel}/'
                                  f'experiments/{StudyLabel}/
                                  f'scans/{SeriesLabel}')
        
        but I get the error:
        
        JSONDecodeError: Expecting value: line 2 column 1 (char 1)
        
        So use more involved work-around: Get a list of scan IDs, get the index
        of the scan ID of interest, get the corresponding scan description, and
        create the directory path:
        """

        request = XnatSession.get(f'{XnatUrl}/data/projects/{ProjId}/'
                                  f'subjects/{SubjLabel}/'
                                  f'experiments/{StudyLabel}/scans')
        
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
        
        Scans = request.json()
        
        SeriesLabels = []
    
        for scan in Scans['ResultSet']['Result']:
            SeriesLabels.append(scan['ID'])
        
        ind = SeriesLabels.index(SeriesLabel)
            
        ScanDesc = Scans['ResultSet']['Result'][ind]['type']
        
        ResourceDir = os.path.join(ExportRootDir, StudyLabel, 'scans',
                                f'{SeriesLabel}-{ScanDesc}', 'resources')
        
        DicomDir = os.path.join(ResourceDir, 'DICOM', 'files')
        
        #AsrDir = os.path.join(ResourceDir, '')
        
        PathsDict[ProjId][SubjLabel][StudyLabel][SeriesLabel].update({'DicomDir' : DicomDir})
    
    
    return PathsDict