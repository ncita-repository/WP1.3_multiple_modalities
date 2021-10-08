# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:27:07 2021

@author: ctorti
"""



def get_dicom_metadata(
        url, proj_id, subj_label, exp_label, scan_id,
        session=None, username=None, password=None
        ):
    """
    Get the StudyUID, SeriesUID, FrameOfReferenceUID, ImagePositionPatients and
    ImageOrientationPatient of a DICOM scan.
    
    Parameters
    ----------
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
        If provided a new session request will be avoided. The default value is 
        None.
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = None)
        the user will be prompted to enter a user name. The default value is 
        None.
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = None)
        the user will be prompted to enter a password. The default value is 
        None.
    
    Returns
    -------
    studyUID : str
        Study Instance UID of the DICOM series / XNAT scan.
    seriesUID : str
        Series Instance UID of the DICOM series / XNAT scan.
    FORUID : str
        Frame of reference UID of the DICOM series / XNAT scan.
    SOPUIDs : list of strs
        List of the SOP UIDs of the DICOM series / XNAT scan.
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
    
    IPPs = []
    
    files = request.json()
    
    for file in files['ResultSet']['Result']:
        if file['collection'] == 'DICOM':
            uri = file['URI']
            
            fulluri = f"{url}{uri}"
    
            # Get the file:
            request = session.get(fulluri)
            
            raw = DicomBytesIO(request.content)
                    
            dcm = dcmread(raw)
            
            IPPs.append(dcm.ImagePositionPatient)
    
    # Other metadata can be parsed from the last file:
    studyUID = dcm.StudyInstanceUID
    seriesUID = dcm.SeriesInstanceUID
    FORuid = dcm.FrameOfReferenceUID
    IOP = dcm.ImageOrientationPatient
    
    return studyUID, seriesUID, FORuid, IPPs, IOP


