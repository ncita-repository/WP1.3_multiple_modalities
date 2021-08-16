# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 08:42:44 2021

@author: ctorti
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