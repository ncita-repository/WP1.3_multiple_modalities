# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 15:42:06 2021

@author: ctorti
"""


def get_insert_dates_by_proj(url, session):
    """
    Return a dictionary containing a list of insert dates for all image session
    uploads organised by project.

    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : Requests Session
        A valid Requests Session for the XNAT of interest.

    Returns
    -------
    DatesByProj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    """
    
    request = session.get(f'{url}/data/experiments')
    
    experiments = request.json()['ResultSet']['Result']
    
    dates_by_proj = {}
    
    for experiment in experiments:
        exp_type = experiment['xsiType']
        project = experiment['project']
        #exp_label = experiment['label']
        exp_id = experiment['ID']
        #exp_date = experiment['date']
        insert_date = experiment['insert_date']
    
        if False:
            print(f' exp_type = {exp_type}')
            print(f' project = {project}')
            #print(f' exp_label = {exp_label}')
            print(f' exp_id = {exp_id}')
            #print(f' exp_date = {exp_date}')
            print(f' insert_date = {insert_date}')
            print('')
    
        if 'SessionData' in exp_type:
            request = session.get(f'{url}/data/experiments/{exp_id}/scans')
    
            #print(request, '\n')
            
            if False:
                if project in list(dates_by_proj.keys()):
                    if 'insert_dates' in list(dates_by_proj[project].keys()):
                        dates_by_proj[project]['insert_dates'].append(insert_date)
    
                    else:
                        #dates_by_proj[project].update({'insert_dates' : [insert_date]}) # equivalent to below
                        dates_by_proj[project]['insert_dates'] = [insert_date]
                else:
                    dates_by_proj[project] = {'insert_dates' : [insert_date]}
            
            if project in list(dates_by_proj.keys()):
                dates_by_proj[project].append(insert_date)
            else:
                dates_by_proj[project] = [insert_date]
    
    return dates_by_proj


def get_first_last_im_session_uploads_by_proj(
        url, session, data_by_proj=None
        ):
    """
    Return a dictionary containing the first and last dates of image session
    uploads organised by project. 
    
    If passed a dictionary (with projects as keys), the upload details will be
    added to existing data (e.g. investigators organised by project).

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
    """
    
    #from datetime import datetime
    
    # Get insert dates by project:
    dates_by_proj = get_insert_dates_by_proj(url, session)
    
    #print(f'dates_by_proj = {dates_by_proj}')
    
    if data_by_proj == None:
        data_by_proj = {}
    
    for project, insert_dates in dates_by_proj.items():
        #print(f'There are {len(insert_dates)} insert dates for project {project}')
        
        """
        if len(insert_dates) == 1:
            first_date = insert_dates[0]
            last_date = insert_dates[0]
        else:
            # Arbitrarily set initial values:
            first_date = insert_dates[0]
            last_date = insert_dates[0]
            
            for date in insert_dates:
                if date < first_date:
                    first_date = date
                if date > last_date:
                    last_date = date
        """
        
        # Arbitrarily set initial values (or only value if list of length 1):
        first_date = insert_dates[0]
        last_date = insert_dates[0]
        
        # Update first_date and last_date:
        for date in insert_dates:
            if date < first_date:
                first_date = date
            if date > last_date:
                last_date = date
            
        if project in list(data_by_proj.keys()):
            if 'First upload' in list(data_by_proj[project].keys()):
                data_by_proj[project]['First upload'] = first_date
            else:
                #data_by_proj[project].update({'First upload' : first_date}) # equivalent to below
                data_by_proj[project]['First upload'] = first_date
            
            if 'Last upload' in list(data_by_proj[project].keys()):
                data_by_proj[project]['Last upload'] = last_date
            else:
                #data_by_proj[project].update({'Last upload' : last_date}) # equivalent to below
                data_by_proj[project]['Last upload'] = last_date
        else:
            data_by_proj[project] = {'First upload' : first_date}
            #data_by_proj.update({project : {'First upload' : first_date}})
            #data_by_proj[project] = {'Last upload' : last_date}
            data_by_proj[project].update({'Last upload' : last_date})
        
                    
    return data_by_proj


def get_start_date_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary containing the date when project metadata was added to
    XNAT, organised by project. 
    
    If passed a dictionary (with projects as keys), the date will be added to 
    existing data (e.g. investigators organised by project).

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
    It seems that start date within children is when the investigator details
    were addded to XNAT, and the start date within meta is when project-related
    meta data was added (e.g. project description).  So neither dates have
    anything to do with first upload.
    """
    from datetime import datetime
    from xnat_tools.projects import get_project_ids
    
    if data_by_proj == None:
        data_by_proj = {}
    
    # Get a list of project IDs:
    proj_ids = get_project_ids(url, session)
    
    for proj_id in proj_ids:
        project = session.get(f'{url}/data/projects/{proj_id}?format=json').json()
        
        date = project['items'][0]['meta']['start_date']
    
        # Convert from ctime (e.g. Fri May 21 22:14:17 UTC 2021) to datetime 
        # (e.g. 2021-05-21 22:14:17):
        # Note: Both below give same results:
        #date = = datetime.strftime(date, '%Y-%m-%d %H:%M:%S')
        date = str(datetime.strptime(date, '%a %b %d %H:%M:%S UTC %Y'))
        
        data_by_proj[proj_id]['Start date'] = date
    
    return data_by_proj
