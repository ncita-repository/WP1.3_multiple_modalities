# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 16:34:29 2021

@author: ctorti
"""


def get_invs_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary containing the PIs and co-investigators organised by
    project.  
    
    If passed a dictionary (with projects as keys), the investigator details
    will be added to existing data (e.g. subjects organised by project).

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
    Can get list of investigators from:
    
    request = session.get(f'{XnatAddress}/data/investigators')
    investigators = request.json()['ResultSet']['Result']
      
    but resulting list of dictionaries don't include projects.
     
    Instead:
    
    request = session.get(f'{XnatAddress}/xapi/investigators')
    investigators = request.json()
    
    returns list of dictionaries containing details of project, and whether the
    investigator is a PI or co-investigator for the project.
    
    A flat list (not organised by project) of investigators (PIs and CIs):
    
    investigators = [f"{item['firstname']} {item['lastname']}" for item in investigators]
    
    Similarly a flat list of PIs:
    
    request = session.get(f'{url}/data/projects')
    projects = request.json()
    PIs = [f"{item['pi_firstname']} {item['pi_lastname']}" for item in projects['ResultSet']['Result']]
    """
    
    request = session.get(f'{url}/xapi/investigators')
    
    investigators = request.json()
    
    if data_by_proj == None:
        data_by_proj = {}
    
    for investigator in investigators:
        try:
            firstname = investigator['firstname']
        except KeyError:
            firstname = 'N/A'
        try:
            lastname = investigator['lastname']
        except KeyError:
            lastname = 'N/A'
        pri_projects = investigator['primaryProjects']
        co_projects = investigator['investigatorProjects']
        
        if pri_projects:
            for project in pri_projects:
                if project in list(data_by_proj.keys()):
                    data_by_proj[project].update({'PI' : f'{firstname} {lastname}'})
                else:
                    data_by_proj[project] = {'PI' : f'{firstname} {lastname}'}
        
        if co_projects:
            for project in co_projects:
                if project in list(data_by_proj.keys()):
                    if 'Co-investigators' in list(data_by_proj[project].keys()):
                        data_by_proj[project]['Co-investigators'].append(f'{firstname} {lastname}')
                    else:
                        data_by_proj[project].update({'Co-investigators' : [f'{firstname} {lastname}']})
                else:
                    data_by_proj[project] = {'Co-investigators' : [f'{firstname} {lastname}']}
    
    # Add the keys 'PI' and 'CIs' to any project that did not have a PI or CI:
    for project in list(data_by_proj.keys()):
        if not 'PI' in list(data_by_proj[project].keys()):
            data_by_proj[project].update({'PI' : ''})
        
        if not 'Co-investigators' in list(data_by_proj[project].keys()):
            data_by_proj[project].update({'Co-investigators' : []})
    
    return data_by_proj


def get_subjs_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary containing the number of subjects organised by project.
    
    If passed a dictionary (with projects as keys), the subject details will be
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
    
    Notes
    -----
    A flat list (not organised by project) of subjects:
    
    subjects = [item['label'] for item in subjects['ResultSet']['Result']]
    
    The list of subjects will be deleted once the number of subjects has been
    computed.
    """
    
    request = session.get(f'{url}/data/subjects')
    
    subjects = request.json()['ResultSet']['Result']
    
    if data_by_proj == None:
        data_by_proj = {}
        
    for subject in subjects:
        project = subject['project']
        label = subject['label']
        
        if project in list(data_by_proj.keys()):
            if 'Subjects' in list(data_by_proj[project].keys()):
                data_by_proj[project]['Subjects'].append(label)
            else:
                data_by_proj[project].update({'Subjects' : [label]})
        else:
            data_by_proj[project] = {'Subjects' : [label]}
    
    # Add the key 'No. of subjects' to all projects:
    for project in list(data_by_proj.keys()):
        try:
            subjects = data_by_proj[project]['Subjects']
            
            # Delete as no longer required:
            del data_by_proj[project]['Subjects']
        except KeyError:
            subjects = []
        
        data_by_proj[project].update({'No. of subjects' : len(subjects)})
    
    return data_by_proj


def get_users_by_project(url, session, data_by_proj=None):
    """
    Return a dictionary containing the project users organised by project. 
    
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
    """
    
    from xnat_tools.projects import get_project_ids
    
    if data_by_proj == None:
        data_by_proj = {}
    
    # Get a list of project IDs:
    proj_ids = get_project_ids(url, session)
    
    for proj_id in proj_ids:
        request = session.get(f'{url}/data/projects/{proj_id}/users')
        
        project_users = request.json()['ResultSet']['Result']
        
        num_users = len(project_users)
        
        users = []
        
        for i in range(num_users):
            fname = project_users[i]['firstname']
            lname = project_users[i]['lastname']
            
            users.append(f"{fname} {lname}")
        
        if proj_id in list(data_by_proj.keys()):
            data_by_proj[proj_id]['Users'] = users
            data_by_proj[proj_id]['No. of users'] = num_users
        else:
            data_by_proj[proj_id] = {'Users' : users}
            data_by_proj[proj_id] = {'No. of users' : num_users}
    
    return data_by_proj