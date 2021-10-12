# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 16:56:00 2021

@author: ctorti
"""


def get_project_ids(url, session):
    """
    Get a list of project IDs for all projects in an XNAT.
    
    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : requests session
        A valid requests session.
    
    Returns
    -------
    project_ids : list of str
        A list (for each project) of the project name.
    """
    
    request = session.get(f'{url}/data/projects')
    
    projects = request.json()['ResultSet']['Result']
    
    project_ids = [project['ID'] for project in projects]

    return project_ids

def get_project_names(url, session):
    """
    Get a list of project names for all projects in an XNAT.
    
    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : requests session
        A valid requests session.
    
    Returns
    -------
    project_names : list of str
        A list (for each project) of the project name.
    """
    
    request = session.get(f'{url}/data/projects')
    
    projects = request.json()['ResultSet']['Result']
    
    project_names = [project['name'] for project in projects]

    return project_names

def get_proj_desc_by_proj(url, session, data_by_proj=None):
    """
    Return a dictionary containing the project description organised by 
    project. 
    
    If passed a dictionary (with projects as keys), the description will be 
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
    
    if data_by_proj == None:
        data_by_proj = {}
    
    # Get a list of project IDs:
    proj_ids = get_project_ids(url, session)
    
    for proj_id in proj_ids:
        project = session.get(f'{url}/data/projects/{proj_id}?format=json').json()
        
        try:
            # Try to get the description (if it exists):
            desc = project['items'][0]['data_fields']['description']
        except KeyError:
            # There is no description for this project:
            desc = 'None'
        
        data_by_proj[proj_id]['Project description'] = desc
    
    return data_by_proj