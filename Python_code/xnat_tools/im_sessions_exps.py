# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 16:48:57 2021

@author: ctorti
"""


def get_im_sessions_by_type(url, session):
    """
    Return a dictionary of image sessions organised by type (e.g. MR, PET).
    

    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : Requests Session
        A valid Requests Session for the XNAT of interest.

    Returns
    -------
    ScanByScanType : dict
        A dictionary with scan type as keys, containing data organised by 
        scan type.
    
    Notes
    -----
    The keys within each sub-dictionary will be `xnat_imagescandata_id', e.g.:
    im_session_by_type['CT'] = {'Project' : 'project_name',
                             'ExpType' : 'experiment_type',
                             'ExpLabel' : 'experiment_label',
                             'ExpId' : 'experiment_id',
                             'Date' : 'experiment_date',
                             'ScanType' : scan_type,
                             'ScanId' : im_scan_data_id,
                             'SeriesDesc' : 'series_description'
                             }
    """
    
    request = session.get(f'{url}/data/experiments')
    
    experiments = request.json()
    
    im_session_by_type = {}
    
    for experiment in experiments['ResultSet']['Result']:
        exp_type = experiment['xsiType']
        project = experiment['project']
        exp_label = experiment['label']
        exp_id = experiment['ID']
        date = experiment['date']
        
        if True:
            print(f' exp_type = {exp_type}')
            print(f' project = {project}')
            print(f' exp_label = {exp_label}')
            print(f' exp_id = {exp_id}')
            print(f' date = {date}')
            print('')
        
        if 'SessionData' in exp_type:
            request = session.get(f'{url}/data/experiments/{exp_id}/scans')
    
            print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' exp_label = {exp_label}')
                print(f' exp_id = {exp_id}')
                print(f' date = {date}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                scan_dict = {}
    
                scan_type = scan['xsiType']
                scan_mod = scan_type.split('xnat:')[1].split('ScanData')[0].upper()
                
                print(f'   scan_type = {scan_type}')
                print(f'   scan_mod = {scan_mod}')
    
                im_scan_data_id = scan['xnat_imagescandata_id']
                series_desc = scan['series_description']
    
                scan_dict = {'Project' : project,
                            'ExpType' : exp_type,
                            'ExpLabel' : exp_label,
                            'ExpId' : exp_id,
                            'Date' : date,
                            'ScanType' : scan_type,
                            'ScanId' : im_scan_data_id,
                            'SeriesDesc' : series_desc}
    
                #if scan_type in list(im_session_by_type.keys()):
                if scan_mod in list(im_session_by_type.keys()):
                    #im_session_by_type[scan_type].append(scan_dict)
                    im_session_by_type[scan_mod].append(scan_dict)
                    
                else:
                    #im_session_by_type.update({scan_type : [scan_dict]})
                    #im_session_by_type.update({scan_mod : [scan_dict]}) # equivalent to below
                    im_session_by_type[scan_mod] = [scan_dict]
    
    return im_session_by_type


def get_im_sessions_by_type_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary containing image sessions organised by image type by
    project. 
    
    If passed a dictionary (with projects as keys), the image session details
    will be added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : Requests Session
        A valid Requests Session for the XNAT of interest.
    data_by_proj : dict, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    The keys within each sub-dictionary will be e.g.:
    data_by_proj['ProjectName]['CT'] = {'ExpType' : 'experiment_type',
                                        'ExpLabel' : 'experiment_label',
                                        'ExpId' : 'experiment_id',
                                        'Date' : 'experiment_date',
                                        'ScanType' : scan_type,
                                        'ScanId' : im_scan_data_id,
                                        'SeriesDesc' : 'series_description'
                                        }
    """
    
    request = session.get(f'{url}/data/experiments')
    
    experiments = request.json()['ResultSet']['Result']
    
    if data_by_proj == None:
        data_by_proj = {}
    
    for experiment in experiments:
        exp_type = experiment['xsiType']
        project = experiment['project']
        exp_label = experiment['label']
        exp_id = experiment['ID']
        exp_date = experiment['date']
        insert_date = experiment['insert_date']
        
        if True:
            print(f' exp_type = {exp_type}')
            print(f' project = {project}')
            print(f' exp_label = {exp_label}')
            print(f' exp_id = {exp_id}')
            print(f' exp_date = {exp_date}')
            print(f' insert_date = {insert_date}')
            print('')
        
        if 'SessionData' in exp_type:
            request = session.get(f'{url}/data/experiments/{exp_id}/scans')
    
            print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' exp_label = {exp_label}')
                print(f' exp_id = {exp_id}')
                print(f' exp_date = {exp_date}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                scan_dict = {}
    
                scan_type = scan['xsiType']
                scan_mod = scan_type.split('xnat:')[1].split('ScanData')[0].upper()
                
                print(f'   scan_type = {scan_type}')
                print(f'   scan_mod = {scan_mod}')
    
                im_scan_data_id = scan['xnat_imagescandata_id']
                series_desc = scan['series_description']
    
                scan_dict = {'ExpType' : exp_type,
                            'ExpLabel' : exp_label,
                            'ExpId' : exp_id,
                            'ExpDate' : exp_date,
                            'ScanType' : scan_type,
                            'ScanId' : im_scan_data_id,
                            'UploadDate' : insert_date,
                            'SeriesDesc' : series_desc}
                
                if project in list(data_by_proj.keys()):
                    if scan_mod in list(data_by_proj[project].keys()):
                        data_by_proj[project][scan_mod].append(scan_dict)
                        
                    else:
                        #data_by_proj[project].update({scan_mod : [scan_dict]}) # equivalent to below
                        data_by_proj[project][scan_mod] = [scan_dict]
                else:
                    data_by_proj[project] = {scan_mod : [scan_dict]}
        
    
    return data_by_proj


def get_num_im_sessions_by_type_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary containing the number of image sessions organised by 
    image type by project. 
    
    If passed a dictionary (with projects as keys), the image session details
    will be added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    url : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : Requests Session
        A valid Requests Session for the XNAT of interest.
    data_by_proj : dictionary or None, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    data_by_proj : dictionary
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    It isn't clear what are the image session types available within XNAT. At
    present that info is not available, e.g. see:
    https://wiki.xnat.org/documentation/xnat-administration/managing-data-types-in-xnat/adding-and-managing-image-scan-data-types
    """
    
    #Alldata_by_proj = get_im_session_by_type_by_proj(url, session, data_by_proj)
    
    
    request = session.get(f'{url}/data/experiments')
    
    experiments = request.json()
    
    if data_by_proj == None:
        data_by_proj = {}
    
    for experiment in experiments['ResultSet']['Result']:
        exp_type = experiment['xsiType']
        project = experiment['project']
        exp_label = experiment['label']
        exp_id = experiment['ID']
        exp_date = experiment['date']
        insert_date = experiment['insert_date']
        
        if False:#True:
            print(f' exp_type = {exp_type}')
            print(f' project = {project}')
            print(f' exp_label = {exp_label}')
            print(f' exp_id = {exp_id}')
            print(f' exp_date = {exp_date}')
            print(f' insert_date = {insert_date}')
            print('')
        
        if 'SessionData' in exp_type:
            request = session.get(f'{url}/data/experiments/{exp_id}/scans')
    
            #print(request, '\n')
    
            if False:#'500' in str(request):
                print(f' project = {project}')
                print(f' exp_label = {exp_label}')
                print(f' exp_id = {exp_id}')
                print(f' exp_date = {exp_date}')
    
    
            scans = request.json()
    
            for scan in scans['ResultSet']['Result']:
                scan_type = scan['xsiType']
                scan_mod = scan_type.split('xnat:')[1].split('ScanData')[0].upper()
                
                key = f"No. of {scan_mod} scans"
                
                #print(f'   scan_type = {scan_type}')
                #print(f'   scan_mod = {scan_mod}')
                #print(f'   key = {key}')
                
                if project in list(data_by_proj.keys()):
                    if key in list(data_by_proj[project].keys()):
                        #print(f'     before = {data_by_proj[project][key]}')
                        data_by_proj[project][key] += 1
                        #print(f'     after = {data_by_proj[project][key]}')
                        
                    else:
                        #print('     before = 0')
                        data_by_proj[project][key] = 1
                        #print(f'     after = {data_by_proj[project][key]}')
                else:
                    #print('     before = 0')
                    #data_by_proj[project][key] = 1
                    data_by_proj[project] = {key : 1}
                    #print(f'     after = {data_by_proj[project][key]}')
        
    
    return data_by_proj


def get_num_exps_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary of the number of experiments organised by project.
    
    If passed a dictionary (with projects as keys), the number of experiments
    will be added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    url : string
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : Requests Session
        A valid Requests Session for the XNAT of interest.
    data_by_proj : dict, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.

    Returns
    -------
    exps_by_type : dict
        A dictionary with scan type as keys, containing data organised by 
        experiment type.
    data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    request = session.get(f'{url}/data/experiments')
    
    experiments = request.json()
    
    if data_by_proj == None:
        data_by_proj = {}
    
    for experiment in experiments['ResultSet']['Result']:
        #exp_type = experiment['xsiType']
        project = experiment['project']
        #exp_label = experiment['label']
        #exp_id = experiment['ID']
        #date = experiment['date']
        #uri = experiment['URI']
        
        #if False:#True:
        #    print(f' exp_type = {exp_type}')
        #    print(f' project = {project}')
        #    print(f' exp_label = {exp_label}')
        #    print(f' exp_id = {exp_id}')
        #    print(f' date = {date}')
        #    print('')
        
        """ The following will return detailed info about the experiment that
        might be useful later:
        
        request = session.get(f'{url}{uri}?format=json')
        
        experiment = request.json()

        print(f"keys in experiment = {list(experiment.keys())}\n")
        
        print(f"len(experiment['items'] = {len(experiment['items'])}\n")
        
        print(f"keys in experiment['items'][0] = {list(experiment['items'][0].keys())}\n")
        
        experiment['items'][0]
        """

        if project in list(data_by_proj.keys()):
            if 'No. of Experiments' in list(data_by_proj[project].keys()):
                #print(f'     before = {data_by_proj[project]['No. of Experiments']}')
                data_by_proj[project]['No. of Experiments'] += 1
                #print(f'     after = {data_by_proj[project]['No. of Experiments']}')
                
            else:
                #print('     before = 0')
                data_by_proj[project]['No. of Experiments'] = 1
                #print(f'     after = {data_by_proj[project]['No. of Experiments']}')
        else:
            #print('     before = 0')
            #data_by_proj[project]['No. of Experiments'] = 1
            data_by_proj[project] = {'No. of Experiments' : 1}
            #print(f'     after = {data_by_proj[project]['No. of Experiments']}')
    
    return data_by_proj
