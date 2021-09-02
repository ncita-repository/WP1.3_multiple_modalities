# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:27:07 2021

@author: ctorti
"""



def get_FOR_uid(url, proj_id, subj_label, exp_label, scan_id,
                session=None, username=None, password=None):
    """
    Get the FrameOfReferenceUID of a single DICOM within a specified scan.
    
    url : str
        Address of XNAT (e.g. 'http://10.1.1.20').
    proj_id : str
        The project ID of interest.
    subj_label : str
        The subject label of interest. 
    exp_label : str
        The DICOM study / XNAT experiment label of interest. 
    scan_id : str
        The DICOM series label / XNAT scan ID of interest.
    session : requests session, optional
        If provided a new session request will be avoided.
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = None)
        the user will be prompted to enter a user name.
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = None)
        the user will be prompted to enter a password.
    
    
    Outputs:
    *******
    
    FORuid : str
        The FrameOfReferenceUID of a DICOM within the DICOM series / XNAT scan
        of interest.
    """
    
    from pydicom.filebase import DicomBytesIO
    from pydicom import dcmread
    
    # Get the Source scan files:
    uri = f'{url}/data/projects/{proj_id}/subjects/{subj_label}/'\
          + f'experiments/{exp_label}/scans/{scan_id}/files'
    
    request = session.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    files = request.json()
    
    # Get the URI of the first DICOM file:
    for file in files['ResultSet']['Result']:
        if file['collection'] == 'DICOM':
            uri = file['URI']
            break
    
    fulluri = f"{url}{uri}"
    
    # Get the file:
    request = session.get(fulluri)
    
    raw = DicomBytesIO(request.content)
            
    dcm = dcmread(raw)
    
    FORuid = dcm.FrameOfReferenceUID
    
    return FORuid


