# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 16:52:35 2021

@author: ctorti
"""


""" 
Functions that create a "snapshot" of an XNAT and will export it to an
XLSX file.

Example usage as a module (run in a command shell)
--------------------------------------------------

python snapshots.py http://10.1.1.20


Example usage as a function (run in a Python command shell)
-----------------------------------------------------------

from xnat_tools.snapshots import get_xnat_snapshot

get_xnat_snapshot('http://10.1.1.20')
"""

import os
import sys

#code_root = r'C:\Code\WP1.3_multiple_modalities\src'
code_root = os.getcwd()

# Add code_root to the system path so packages can be imported from it:
sys.path.append(code_root)

import time
#from getpass import getpass
#import requests
#from pathlib import Path
import datetime
import argparse

import xnat_tools.sessions
import xnat_tools.invs_subjs_users
import xnat_tools.file_count_size
import xnat_tools.projects
import xnat_tools.dates_times
import xnat_tools.format_pathsDict
import xnat_tools.alias_tokens
import io_tools.exports

from importlib import reload
reload(xnat_tools.sessions)
reload(xnat_tools.invs_subjs_users)
reload(xnat_tools.file_count_size)
reload(xnat_tools.projects)
reload(xnat_tools.dates_times)
reload(xnat_tools.format_pathsDict)
reload(xnat_tools.alias_tokens)
reload(io_tools.exports)

from xnat_tools.sessions import create_session
from xnat_tools.invs_subjs_users import get_invs_by_proj, get_subjs_by_proj
from xnat_tools.invs_subjs_users import get_users_by_project
#from xnat_tools.experiments import get_num_exps_by_proj
#from xnat_tools.im_sessions_exps import get_im_sessions_by_type_by_proj
from xnat_tools.file_count_size import get_file_count_size_by_type_by_proj
from xnat_tools.projects import get_proj_desc_by_proj
#from xnat_tools.dates_times import get_start_date_by_proj
from xnat_tools.dates_times import get_first_last_im_session_uploads_by_proj
from xnat_tools.format_pathsDict import reorder_keys_and_fill_zeros
from xnat_tools.alias_tokens import (
    import_alias_token, is_alias_token_valid, generate_alias_token,
    export_alias_token
    )
from io_tools.exports import (
    export_dict_to_xlsx, export_dict_of_dicts_to_xlsx
    )



def get_xnat_snapshot(
        url, username="", password="", session=None,
        export_xlsx=True, log_to_console=False
        ):
    """
    Get a "snapshot" of info from an XNAT broken down by projects, and 
    XNAT-wide.
    
    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    username : str, optional
        The username for XNAT log-in.  If not provided (i.e. username = "")
        the user will be prompted to enter a user name. The default value is
        "".
    password : str, optional
        The password for XNAT log-in.  If not provided (i.e. password = "")
        the user will be prompted to enter a password. The default value is
        "".
    session : requests session, optional
        If provided a multiple XNAT session requests will be avoided. The 
        default value is None.
    export_xlsx : bool, optional
        If True the "snapshot" will be exported to an XLSX file. The default
        value is True.
    log_to_console : bool, optional
        If True some results will be printed to the console. The default value
        is False.
    
    Returns
    -------
    data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    data_xnat_wide : dict
        A dictionary containing XNAT-wide analytics.
    session : requests session
    
    Notes
    -----
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
    
    If export_xlsx = True, the data dictionary will be exported to the
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
    
    if session == None:
        # Get current working directory:
        cwd = os.getcwd()
        
        # Assumed directory that may contain an XNAT alias token:
        tokenDir = os.path.join(cwd, 'tokens')
        
        print(f'Searching for an XNAT alias token in {tokenDir}')
        # Import an XNAT alias token:
        aliasToken = import_alias_token(tokenDir)
        
        # Is the alias token valid?
        if is_alias_token_valid(aliasToken, url):
            print('Using XNAT alias token to establish XNAT connection.')
            session = create_session(url=url, aliasToken=aliasToken)
        else:
            print('Establishing XNAT connection without XNAT alias token.')
            session = create_session(url=url, username=username)
    
    # Generate an XNAT alias token to avoid making a new authenticated
    # user session:
    aliasToken = generate_alias_token(session)
    
    export_alias_token(aliasToken, tokenDir)
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    #""" Get a list of projects: """
    #request = session.get(f'{url}/data/projects')
    #projects = request.json()
    #Projects = [item['name'] for item in projects['ResultSet']['Result']]
    
    
    #""" Get a list of PIs: """
    #PIs = [f"{item['pi_firstname']} {item['pi_lastname']}" for item in projects['ResultSet']['Result']]
    
    
    """ Get a list of investigators grouped by project: """
    data_by_proj = get_invs_by_proj(url, session)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to fetch investigators.\n')
    
    """ Get a list of subjects grouped by project: """
    data_by_proj = get_subjs_by_proj(url, session, data_by_proj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to fetch subjects.\n')
    
    """ Get the number of experiments by project: """
    #data_by_proj = get_num_exps_by_proj(url, session, data_by_proj)
    
    
    """ Get the number of scans by image session type by project: """
    #data_by_proj = get_im_sessions_by_type_by_proj(url, session, 
    #                                              data_by_proj)
    
    """ Get the number of experiments, image sessions by modality, image scans 
    by modality, and average number of scans per session, all organised by 
    project: """
    #data_by_proj = GetExpsSessionsScansByProject(url, session, data_by_proj)
    
    """ Get the no. of experiments, no. of image sessions by modality, no. of
    image scans by modality, ave. scans/session, ave. image scans/session,
    no. of DICOM/RTSTRUCT/AIM/SEG files, total size of DICOM/RTSTRUCT/SEG files,
    etc, all organised by project: """
    data_by_proj = get_file_count_size_by_type_by_proj(
        url, session, data_by_proj
        )
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to fetch file counts and sizes.\n')
    
    """ Get the project description by project: """
    data_by_proj = get_proj_desc_by_proj(url, session, data_by_proj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to fetch project descriptions.\n')
    
    """ Start date might not be the same as first upload.. """
    #""" Get the start date by project: """
    #data_by_proj = get_start_date_by_proj(url, session, data_by_proj)
    
    
    """ Get the project users by project: """
    data_by_proj = get_users_by_project(url, session, data_by_proj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to fetch users.\n')
    
    """ Get the first and last upload dates by project: """
    data_by_proj = get_first_last_im_session_uploads_by_proj(
        url, session, data_by_proj
        )
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to fetch image session uploads.\n')
    
    """ Re-order the keys at the second tier: """
    data_by_proj = reorder_keys_and_fill_zeros(data_by_proj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to reorder keys.\n')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 2)
    print(f'*Took {Dtime} s to generate data-by-project snapshot.\n')
    
    """ Get XNAT-wide analytics (not broken down by project, i.e.:
        - no. of unique users
        - average no. of users per project
    """
    data_xnat_wide = get_xnat_wide_snapshot(data_by_proj)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 2)
    if log_to_console:
        print(f'*Took {Dtime} s to generate XNAT-wide snapshot.\n')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 2)
    print(f'*Took {Dtime} s to get all snapshots.\n')
    
    
    """ Export the dictionary: """
    if export_xlsx:
        # Directory where snapshots will be exported to:
        exportDir = os.path.join(cwd, 'xnat_snapshots')
        
        time_now = datetime.datetime.now()
        date_time = time_now.strftime('%Y%m%d_%H%M%S')
        
        file_name = date_time + '_XNAT_by_project_snapshot.xlsx'
        
        # Export the by-project snapshot:
        export_dict_of_dicts_to_xlsx(
            data_by_proj, file_name, exportDir=exportDir, 
            nameOfFirstCol='Project'
            )
        
        file_path = os.path.join(exportDir, file_name)
        print(f'XNAT-by-project snapshot exported to {file_path}\n')
        
        file_name = date_time + '_XNAT_wide_snapshot.xlsx'
        
        # Export the XNAT-wide snapshot:
        export_dict_to_xlsx(data_xnat_wide, file_name, exportDir=exportDir)
        
        file_path = os.path.join(exportDir, file_name)
        print(f'XNAT-wide snapshot exported to {file_path}\n')
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'*Took {Dtime} s to export the snapshots to XLSX files.\n')
    
    return data_by_proj, data_xnat_wide, session

def get_xnat_wide_snapshot(data_by_proj):
    """ 
    Get XNAT-wide analytics (not broken down by project).
    
    Parameters
    ----------
    data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Returns
    -------
    data_xnat_wide : dict
        A dictionary containing XNAT-wide analytics.
    """
    
    # Number of decimals for rounding:
    #Decimals = 2
    
    num_projects = len(list(data_by_proj.keys()))
    
    # Initialise:
    unique_users = []
    num_users = 0
    #num_users_by_proj = []
    #num_subjs_by_proj = []
    num_subjs = 0
    num_exps = 0
    num_MR_sessions = 0
    num_CT_sessions = 0
    num_PET_sessions = 0
    num_MR_scans = 0
    num_CT_scans = 0
    num_PET_scans = 0
    num_OT_scans = 0
    num_dcm_im_files = 0
    size_dcm_im_files = 0
    
    
    for proj_id in data_by_proj.keys():
        users = data_by_proj[proj_id]['Users']
        
        num_users += len(users)
        
        #num_users_by_proj.append(len(users))
        
        for user in users:
            if not user in unique_users:
                unique_users.append(user)
        
        num_subjs += data_by_proj[proj_id]['No. of subjects']
        num_exps += data_by_proj[proj_id]['No. of experiments']
        num_MR_sessions += data_by_proj[proj_id]['No. of MR sessions']
        num_CT_sessions += data_by_proj[proj_id]['No. of CT sessions']
        num_PET_sessions += data_by_proj[proj_id]['No. of PET sessions']
        num_MR_scans += data_by_proj[proj_id]['No. of MR scans']
        num_CT_scans += data_by_proj[proj_id]['No. of CT scans']
        num_PET_scans += data_by_proj[proj_id]['No. of PET scans']
        num_OT_scans += data_by_proj[proj_id]['No. of OT scans']
        num_dcm_im_files += data_by_proj[proj_id]['No. of DICOM image files [k]']
        size_dcm_im_files += data_by_proj[proj_id]['Size of DICOM image files [MB]']
        
        #num_subjs_by_proj.append(data_by_proj[proj_id]['No. of subjects'])
    
    # Averages per project:
    ave_num_users_pp = round(num_users/num_projects, 1)
    #ave_num_users_pp = round(num_users/num_projects, 1)
    ave_num_subjs_pp = round(num_subjs/num_projects, 1)
    ave_num_exps_pp = round(num_exps/num_projects, 1)
    ave_num_MR_sessions_pp = round(num_MR_sessions/num_projects, 1)
    ave_num_CT_sessions_pp = round(num_CT_sessions/num_projects, 1)
    ave_num_PET_sessions_pp = round(num_PET_sessions/num_projects, 1)
    ave_num_MR_scans_pp = round(num_MR_scans/num_projects, 1)
    ave_num_CT_scans_pp = round(num_CT_scans/num_projects, 1)
    ave_num_PET_scans_pp = round(num_PET_scans/num_projects, 1)
    ave_num_OT_scans_pp = round(num_OT_scans/num_projects, 1)
    ave_num_dcm_im_files_pp = round(num_dcm_im_files/num_projects, 2)
    ave_size_dcm_im_files_pp = round(size_dcm_im_files/num_projects, 2)
    
    
    data_xnat_wide = {
        'No. of unique users' : len(unique_users),
        'Ave. users/project' : ave_num_users_pp,
        'No. of projects' : num_projects,
        'No. of subjects' : num_subjs,
        'Ave. subjects/project' : ave_num_subjs_pp,
        'No. of experiments' : num_exps,
        'Ave. experiments/project' : ave_num_exps_pp,
        'No. of MR sessions' : num_MR_sessions,
        'Ave. MR sessions/project' : ave_num_MR_sessions_pp,
        'No. of CT sessions' : num_CT_sessions,
        'Ave. CT sessions/project' : ave_num_CT_sessions_pp,
        'No. of PET sessions' : num_PET_sessions,
        'Ave. PET sessions/project' : ave_num_PET_sessions_pp,
        'No. of MR scans' : num_MR_scans,
        'Ave. MR scans/project' : ave_num_MR_scans_pp,
        'No. of CT scans' : num_CT_scans,
        'Ave. CT scans/project' : ave_num_CT_scans_pp,
        'No. of PET scans' : num_PET_scans,
        'Ave. PET scans/project' : ave_num_PET_scans_pp,
        'No. of OT scans' : num_OT_scans,
        'Ave. OT scans/project' : ave_num_OT_scans_pp,
        'No. of DICOM image files [M]' : round(num_dcm_im_files/1000, 2),
        'Ave. DICOM image files/project [k]' : ave_num_dcm_im_files_pp,
        'Size of DICOM image files [GB]' : round(size_dcm_im_files/1000, 2),
        'Ave. size of DICOM image files/project [MB]' : ave_size_dcm_im_files_pp
        }
    
    return data_xnat_wide

if __name__ == '__main__':
    """
    Run snapshots.py as a script.
    
    Example usage in a console:
    
    python snapshots.py http://10.1.1.20
    """
    
    parser = argparse.ArgumentParser(
        description='Arguments for get_xnat_snapshot()'
        )
    
    parser.add_argument(
        'url', 
        help='The XNAT url'
        )
    
    parser.add_argument(
        '--username', '-u',
        nargs='?',
        default='',
        const='',
        help='The XNAT user name'
        )
    
    parser.add_argument(
        '--password', '-p',
        nargs='?',
        default='',
        const='',
        help='The XNAT password'
        )
    
    parser.add_argument(
        '--session', '-s',
        nargs='?',
        default=None,
        const=None,
        help='An existing XNAT Requests session'
        )
    
    parser.add_argument(
        '--export_xlsx', '-e',
        nargs='?',
        default=True,
        const=True,
        help='Export to XLSX files?'
        )
    
    parser.add_argument(
        '--log_to_console', '-l',
        nargs='?',
        default=False,
        const=False,
        help='Log output to the console?'
        )
    
    args = parser.parse_args()
    
    # Run get_xnat_snapshot():
    get_xnat_snapshot(
        args.url, args.username, args.password, args.session, args.export_xlsx,
        args.log_to_console
        )