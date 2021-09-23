# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:32:37 2021

@author: ctorti
"""


import os
from pathlib import Path
from xnat_tools.sessions import create_session
from xnat_tools.format_pathsDict import create_pathsDict_for_subj_asr
    
    
def download_subj_asrs(
        url, proj_id, subj_label, session=None, 
        export_root_dir='', paths_dict=None, 
        username=None, password=None
        ):
    """
    Download all assessors for a subject of interest from XNAT.
    
    Parameters
    ----------
    url : str
        Address of XNAT (e.g. 'http://10.1.1.20').
    proj_id : str
        The project ID of interest.
    subj_label : str
        The subject label of interest. 
    session : requests session, optional
        If provided a new session request will be avoided.
    export_root_dir : str, optional
        If provided the data will be downloaded to the chosen directory, 
        organised by sub-directories for ProjectLabel, subj_label, Scans or
        Assessors, etc. If not provided, the data will be downloaded to the
        sub-directory "XnatDownloads" within the default downloads directory.
    paths_dict : dict, optional
        Dictionary containing paths of data downloaded.
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = None)
        the user will be prompted to enter a user name.
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = None)
        the user will be prompted to enter a password.
    
    Returns
    -------
    paths_dict : dict
        Dictionary containing paths of data downloaded.
    session : requests session
    """
    
    #import os
    #from pathlib import Path
    #from xnat_requests import create_session
    
    if session == None:
        session = create_session(url, username, password)
    
    if export_root_dir == '':
        export_root_dir = os.path.join(
            Path.home(), "Downloads", "xnat_downloads"
            )
    
    # Get a listing of all resource files for the subject:
    uri = f'{url}/data/projects/{proj_id}/subjects/{subj_label}/files'
                              
    request = session.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    for result in request.json()['ResultSet']['Result']:
        fulluri = f"{url}{result['URI']}"
        
        res_id = result['cat_ID']
        
        name = result['Name']
        
        export_dir = os.path.join(
            export_root_dir, 'projects', proj_id, 'subjects', subj_label,
            'resources', res_id, 'files', name
            )
        
        subj_asr_fpath = os.path.join(export_dir, name)
        
        request = session.get(fulluri)
        
        with open(subj_asr_fpath, 'wb') as file:
            file.write(request.content)
        
        """ Store info in a dictionary: """
        paths_dict, keys = create_pathsDict_for_subj_asr(
            proj_id, subj_label, res_id, name, paths_dict
            )
        
        paths_dict['projects'][proj_id]['subjects'][subj_label]['resources']\
        [res_id]['files'][name].update(
            {'name' : name,
             'res_id' : res_id,
             'subj_asr_dir' : export_dir,
             'subj_asr_fpath' : subj_asr_fpath
             }
            )
    
    return paths_dict, session

def upload_subj_asr(
        subj_asr_fpath, url, proj_id, subj_label, 
        content_label='DRO', session=None, 
        username=None, password=None):
    """
    Upload a subject assessor to a particular subject in XNAT.
    
    Parameters
    ----------
    subj_asr_fpath : str
        Full filepath of assessor to be uploaded to XNAT.
    url : str
        Address of XNAT (e.g. 'http://10.1.1.20').
    proj_id : str
        The project ID of interest.
    subj_label : str
        The subject label of interest. 
    content_label : str, optional
        String to assign to content label. Suggested values are:
            - 'SRO_DRO' for Spatial Registration Object type of DICOM 
            Registration Object
            - 'DSRO_DRO' for Deformable Spatial Registration Object type of 
            DICOM Registration Object
            - 'DRO' generic
        The default value is 'DRO'.
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
    
    #import os
    #from pathlib import Path
    
    if session == None:
        session = create_session(url, username, password)
    
    fname_with_ext = os.path.split(subj_asr_fpath)[1]
    
    if content_label:
        fname, ext = os.path.splitext(
            os.path.split(fname_with_ext)[1]
            )
        
        fname += f'_{content_label}'
        fname_with_ext = fname + ext
    
    # Upload the resource:
    with open(subj_asr_fpath, 'rb') as file:
        uri = f'{url}/data/projects/{proj_id}/subjects/{subj_label}/' +\
            f'resources/DICOM/files/{fname_with_ext}'
        
        request = session.put(
            f'{uri}?inbody=true&file_format=DICOM&overwrite=true',
            data=file
            )
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
        
    print(f'The DRO was uploaded as a subject assessor to: \n{uri}\n')
    
    return session