# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 11:58:55 2021

@author: ctorti
"""


def get_file_count_size_by_type_by_proj(url, session, data_by_proj=None,
                                        log_to_console=False):
    """
    06/07/21 NOTE: 
        This works fine for my XNAT but fails to run on XNAT Central due to 
        missing keys. Use of exception handling isn't sufficient. More robust 
        approach needed to deal with missing data.
        
    
    Return a dictionary containing the number of files and the total size of
    files by modality, all organised by project.
    
    If passed a dictionary (with projects as keys), the above details will be
    added to existing data (e.g. investigators organised by project).

    Parameters
    ----------
    url : str
        URL of XNAT (e.g. 'http://10.1.1.20').
    session : requests session
        A valid Requests Session for the XNAT of interest.
    data_by_proj : dict, optional
        A dictionary with projects as keys, containing data organised by 
        projects. The default is None.
    log_to_console : bool, optional
        If True some info will be printed to the console.

    Returns
    -------
    data_by_proj : dict
        A dictionary with projects as keys, containing data organised by 
        projects.
    
    Notes
    -----
    18/06/21: Expansion of combination of get_num_exps_by_project and 
    get_num_im_session_by_type_by_project, and replaces the latter.
    
    request = session.get(f'{url}{uri}?format=json')
    
    or
    
    request = session.get(f'{url}/data/experiments/{expId}?format=json')
    
    exp = request.json()

    
    exp is a dictionary with one key: 'items'.
        
    exp['items'] is a list of length 1.
    
    exp['items'][0] is a dictionary with keys: 'children', 'meta'
    and 'data_fields'.
    
    exp['items'][0]['meta'] and Exp['items'][0]['data_fields']
    are simple dictionaries.
    
    exp['items'][0]['children'] is a list of variable length. It
    seems that if there are image assessors, 
    exp['items'][0]['children'][0] relate to assessors, and
    exp['items'][0]['children'][1] relate to scans.
    
    exp['items'][0]['children'][i] is a dictionary with keys:
    'field' and 'items'.
    
    exp['items'][0]['children'][i]['field'] is either 
    'assessors/assessor' or 'scans/scan'.
    
    exp['items'][0]['children'][i]['items'] is a list for each
    assessor/scan.
    
    Each item in Exp['items'][0]['children'][i]['items'] is a 
    dictionary with keys: 'children', 'meta' and 'data_fields'.
    
    exp['items'][0]['children'][i]['items'][j]['meta'] and 
    exp['items'][0]['children'][i]['items'][j]['data_fields'] are
    simple dictionaries.
    
    exp['items'][0]['children'][i]['items'][j]['children'] is a
    list of length 2.
    
    exp['items'][0]['children'][i]['items'][j]['children'][0] is a
    dictionary with keys 'field' and 'items'.
    
    exp['items'][0]['children'][i]['items'][j]['children'][0]['field']
    seems to always be 'out/file' (for assessors) or 'file' (for scans).
    
    If the item is an assessor:
        If the assessor is RTSTRUCT:
            exp['items'][0]['children'][i]['items'][j]['children'][0]['items']
            is a list of length 2.
             
            exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][0]
            relates to the AIM, and 
            
            exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][1]
            relates to the RTSTRUCT.
    
        If the assessor is SEG:
            exp['items'][0]['children'][i]['items'][j]['children'][0]['items']
            is a list of length 1.
             
            exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][0]
            relates to the SEG.
        
    
    exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][k] 
    is a dictionary with keys 'children' (= []), 'meta' and 'data_fields'.
    
    The file_size is located in
    exp['items'][0]['children'][i]['items'][j]['children'][0]['items'][k]['data_fields'].
    
    exp['items'][0]['children'][i]['items'][j]['children'][1] is a
    simple dictionary that references the SeriesUID, and has keys 'field' 
    and 'items'.  
    
    exp['items'][0]['children'][i]['items'][j]['children'][1]['field']
    seems to always be 'references/seriesUID'.
    
    exp['items'][0]['children'][i]['items'][j]['children'][1]['items']
    is a list of length 1.
    
    exp['items'][0]['children'][i]['items'][j]['children'][1]['items'][0]
    is a dictionary with keys 'children' (= []), 'meta' and 'data_fields'.
    """
    
    import sys
    from xnat_tools.projects import get_project_ids
    
    if data_by_proj == None:
        data_by_proj = {}
    
    empty_dict = {'No. of experiments' : 0,
                  'No. of MR sessions' : 0,
                  'No. of CT sessions' : 0,
                  'No. of PET sessions' : 0,
                  'Ave. image scans/session' : 0,
                  #'Ave. scans/session' : 0,
                  'No. of MR scans' : 0,
                  'No. of CT scans' : 0,
                  'No. of PET scans' : 0,
                  'No. of OT scans' : 0,
                  'No. of DICOM image files [k]' : 0,
                  'Size of DICOM image files [MB]' : 0,
                  'No. of RTSTRUCT scans' : 0,
                  #'Size of RTSTRUCT scans [MB]' : 0,
                  #'No. of AIM ROI files' : 0,
                  #'Size of AIM ROI files [MB]' : 0,
                  'No. of RTSTRUCT ROI files' : 0,
                  #'Size of RTSTRUCT ROI files [MB]' : 0,
                  'No. of SEG ROI files' : 0,
                  #'Size of SEG ROI files [MB]' : 0,
                  'No. of DICOM files [k]' : 0,
                  'Size of DICOM files [MB]' : 0
                  }
    
    proj_ids = get_project_ids(url, session)
    
    for proj_id in proj_ids:
        if proj_id in data_by_proj.keys():
            data_by_proj[proj_id].update(empty_dict)
        else:
            #data_by_proj[proj_id] = {}
            #data_by_proj[proj_id] = empty_dict
            #data_by_proj[proj_id] = deepcopy(empty_dict)
            data_by_proj[proj_id] = dict(empty_dict) # or empty_dict.copy()
            
            #print('empty_dict = ', empty_dict, '\n')
            
        
        # Get all experiments for this proj_id:
        request = session.get(f'{url}/data/projects/{proj_id}/experiments')
        
        """ 06/07/21: When running snapshot method on XNAT Central I get
        
        HTTPError 403 Client Error: Forbidden for url: 
        https://central.xnat.org/data/projects/{proj_id}/experiments
        
        so try and except for HTTPError.
        """
        
        try:
            exps = request.json()['ResultSet']['Result']
            
            if log_to_console:
                print('****************************************************')
                print(f"Project {proj_id} has {len(exps)} experiments.\n")
            
            key = 'No. of experiments'
            data_by_proj[proj_id][key] = len(exps)
                
            URIs = [exp['URI'] for exp in exps]
            
            num_sessions_this_proj = 0 # initial value of no. of image sessions for this project
            num_scans_this_proj = 0 # initial value of the no. of scans for this project
            num_im_scans_this_proj = 0 # initial value of the no. of image scans for this project
            
            # Loop through each experiment:
            for exp_no in range(len(exps)):
                exp_type = exps[exp_no]['xsiType']
                exp_mod = exp_type.split('xnat:')[1].split('SessionData')[0].upper()
                #project = exps[exp_no]['project']
                exp_label = exps[exp_no]['label']
                exp_id = exps[exp_no]['ID']
                exp_date = exps[exp_no]['date']
                insert_date = exps[exp_no]['insert_date']
            
                if log_to_console:
                    print('----------------------------------------------------')
                    print(f"Experiment {exp_no}:\n")
                    print(f' exp_type = {exp_type}')
                    print(f' exp_mod = {exp_mod}')
                    print(f' Project = {proj_id}')
                    print(f' exp_label = {exp_label}')
                    print(f' exp_id = {exp_id}')
                    print(f' exp_date = {exp_date}')
                    print(f' insert_date = {insert_date}\n')
                
                key = f'No. of {exp_mod} sessions'
                if key in data_by_proj[proj_id].keys():
                    data_by_proj[proj_id][key] += 1
                else:
                    data_by_proj[proj_id][key] = 1
                    
                if 'Session' in exp_type:
                    num_sessions_this_proj += 1
                
                #URI = URIs[exp_no]
                
                # Get more details for this experiment: 
                request = session.get(f'{url}{URIs[exp_no]}?format=json')
                
                exp = request.json()['items'][0]['children']
                
                if log_to_console:
                    print(f" Experiment {exp_no} has {len(exp)} items.\n")
                
                for exp_item_no in range(len(exp)):
                    #exp_item_type = exp[exp_item_no]['field'] # assessors/assessor or scans/scan
                    
                    if 'assessor' in exp[exp_item_no]['field']:
                        exp_item_type = 'assessor'
                    elif 'scan' in exp[exp_item_no]['field']:
                        exp_item_type = 'scan'
                    else:
                        msg = f"Unexpected field value {exp[exp_item_no]['field']}."\
                              + "Was expecting 'assessors/assessors' or 'scans/scan'."
                        raise Exception(msg)
                    
                    #print(f"    Experiment item {exp_item_no} is a {exp_item_type}.\n")
                    
                    if log_to_console:
                        print(f"    Item {exp_item_no} is a {exp_item_type} and has",
                              f"{len(exp[exp_item_no]['items'])} items.\n")
                    
                    # Loop through each assessor/scan sub-item:
                    for sub_item_no in range(len(exp[exp_item_no]['items'])):
                        # Proceed only for items with field 'out/file' (assessors)
                        # or 'file' (scans) (there are also items with field 
                        # 'references/seriesUID':
                        field = exp[exp_item_no]['items'][sub_item_no]['children'][0]['field']
                        
                        if log_to_console:
                            print(f"      Sub-item {sub_item_no} has field {field}.")
                        
                        if 'file' in field:
                            
                            #label = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][0]['data_fields']['label']
                            #item_format = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][0]['data_fields']['format']
                            #
                            #print(f"      Sub-item {sub_item_no} has label {label}.")
                            #print(f"      Sub-item {sub_item_no} has format {item_format}.")
                            
                            # modality doesn't exist for AIM:
                            #modality = exp[exp_item_no]['items'][sub_item_no]['data_fields']['modality']
                            #try:
                            #    modality = exp[exp_item_no]['items'][sub_item_no]['data_fields']['modality']
                            #except KeyError:
                            #    #return exp[exp_item_no]['items'][sub_item_no]
                            #    modality = 'N/A'
                            
                            """ Seems that modality is not accessible from exp[exp_item_no]['items'][sub_item_no]['children']...
                            as are label, item_format, file_size, etc. (as obtained by looping through each sub-sub-item below). """
                            if exp_item_type == 'assessor':
                                modality = exp[exp_item_no]['items'][sub_item_no]['data_fields']['collectionType']
                            else: # exp_item_type = 'scan'
                                modality = exp[exp_item_no]['items'][sub_item_no]['data_fields']['modality']
                            
                            
                            """ Change 'PT' to 'PET' for consistency: """
                            if modality == 'PT':
                                modality = 'PET'
                            
                            N = len(exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'])
                            
                            if log_to_console:
                                print(f"\n      Sub-item {sub_item_no} has {N} items.")
                            
                            # Loop through each sub-sub-item:
                            for sub_sub_item_no in range(N):
                                #if log_to_console:
                                #    print(exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no], '\n')
                                
                                label = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no]['data_fields']['label']
                                item_format = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no]['data_fields']['format']
                                
                                if log_to_console:
                                    print(f"        Sub-sub-item {sub_sub_item_no} has modality {modality}.")
                                    print(f"        Sub-sub-item {sub_sub_item_no} has label {label}.")
                                    print(f"        Sub-sub-item {sub_sub_item_no} has format {item_format}.")
                                
                                # File size is not listed for SNAPSHOTS so skip them:
                                if label == 'SNAPSHOTS':
                                    #return exp[exp_item_no]['items'][sub_item_no]
                                
                                    if log_to_console:
                                        print(f"      Skipping sub_sub_item_no {sub_sub_item_no} since it is type {label}.\n")
                                
                                else:
                                    num_scans_this_proj += 1
                                    
                                    if label == 'DICOM':
                                        # This is an image scan:
                                        num_im_scans_this_proj += 1
                                    
                                    file_count = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no]['data_fields']['file_count']
                                    file_size = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no]['data_fields']['file_size']
                                    
                                    if log_to_console:
                                        print(f"        Sub-sub-item {sub_sub_item_no} has file_count {file_count}.")
                                        print(f"        Sub-sub-item {sub_sub_item_no} has file_size {file_size}.\n")
                                        
                                        if label == 'secondary':
                                            print(exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no], '\n')
                                    
                                    # Proceed only for sub-sub-items with field 'out/file' (ignore 'references/seriesUID'):
                                    if 'file' in exp[exp_item_no]['items'][sub_item_no]['children'][0]['field']:
                                        file_count = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no]['data_fields']['file_count']
                                        file_size = exp[exp_item_no]['items'][sub_item_no]['children'][0]['items'][sub_sub_item_no]['data_fields']['file_size']
                                        file_size = file_size/1000000 # convert to MB
                                        
                                        #print(f"        Sub-sub-item {sub_sub_item_no} has modality {modality}.")
                                        #print(f"        Sub-sub-item {sub_sub_item_no} has file_count {file_count}.")
                                        #print(f"        Sub-sub-item {sub_sub_item_no} has file_size {file_size}.\n")
                                        
                                        #if label == 'AIM':
                                        #    return exp[exp_item_no]['items']#[sub_item_no]
                                        
                                        #if label == 'SEG':
                                        #    return exp[exp_item_no]['items'][sub_item_no]
                                        #    print('')
                                        
                                        #if label == 'secondary':
                                        #    return exp[exp_item_no]['items'][sub_item_no]#['children'][0]['items'][0]
                                        #    return exp[exp_item_no]['items'][sub_item_no - 1]
                                        #    print('')
                                        
                                        #if modality == 'OT':
                                        #    return exp[exp_item_no]['items'][sub_item_no]
                                        
                                        #if label in ['AIM', 'RTSTRUCT', 'SEG']:
                                        if label in ['RTSTRUCT', 'SEG']:
                                            key = f'No. of {label} ROI files'
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += file_count
                                            else:
                                                data_by_proj[proj_id][key] = file_count
                                            
                                            """ Leave out file size for secondary
                                            ROI Collections for now:
                                            key = f'Size of {label} ROI files [MB]'
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += file_size
                                            else:
                                                data_by_proj[proj_id][key] = file_size
                                            """
                                        elif label == 'DICOM':
                                            # This is a DICOM scan (images)
                                            
                                            key = 'No. of DICOM image files [k]'
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += file_count/1000
                                            else:
                                                data_by_proj[proj_id][key] = file_count/1000
                                            
                                            key = 'Size of DICOM image files [MB]'
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += file_size
                                            else:
                                                data_by_proj[proj_id][key] = file_size
                                            
                                            key = f"No. of {modality} scans"
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += 1
                                            else:
                                                data_by_proj[proj_id][key] = 1        
                                        
                                        if item_format == 'DICOM':
                                            # This is a DICOM scan or RTSTRUCT scan or
                                            # ROI Collection of type AIM, RTSTRUCT or SEG:
                                            
                                            key = 'No. of DICOM files [k]'
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += file_count/1000
                                            else:
                                                data_by_proj[proj_id][key] = file_count/1000
                                            
                                            key = 'Size of DICOM files [MB]'
                                            if key in data_by_proj[proj_id].keys():
                                                data_by_proj[proj_id][key] += file_size
                                            else:
                                                data_by_proj[proj_id][key] = file_size
                                            
                                            if label == 'secondary': 
                                                # Secondary DICOMs, e.g. RTSTRUCT or OT scan
                                                """ Replaced file_count with 1 counting scans
                                                not number of files within scans. """
                                                key = f'No. of {modality} scans'
                                                if key in data_by_proj[proj_id].keys():
                                                    #data_by_proj[proj_id][key] += file_count
                                                    data_by_proj[proj_id][key] += 1
                                                else:
                                                    #data_by_proj[proj_id][key] = file_count
                                                    data_by_proj[proj_id][key] = 1
                                                
                                                """ Leave out file size for secondary
                                                DICOMs for now: 
                                                key = f'Size of {modality} scans [MB]'
                                                if key in data_by_proj[proj_id].keys():
                                                    data_by_proj[proj_id][key] += file_size
                                                else:
                                                    data_by_proj[proj_id][key] = file_size
                                                """
                 
                    #key = 'No. of experiments'
                    ##if key in data_by_proj[proj_id].keys():
                    ##    data_by_proj[proj_id][key] += 1
                    ##else:
                    ##    data_by_proj[proj_id][key] = 1
                    #data_by_proj[proj_id][key] += 1
                
                if log_to_console:
                    print('----------------------------------------------------\n')
                
            """ Get the average number of scans by project: """
            #ave_scans = round(num_scans_this_proj/num_sessions_this_proj, 1)
            ave_im_scans = round(num_im_scans_this_proj/num_sessions_this_proj, 1)
            
            #data_by_proj[proj_id]['Ave. scans/session'] = ave_scans
            data_by_proj[proj_id]['Ave. image scans/session'] = ave_im_scans
            
            if log_to_console:
                print('****************************************************\n')
                    
        except:
            if log_to_console:
                print(f'{sys.exc_info()[0]}: \n    {sys.exc_info()[1]}')
            
            data_by_proj[proj_id]['Ave. image scans/session'] = 0
            
        
        #""" Get the average number of scans by project: """
        ##ave_scans = round(num_scans_this_proj/num_sessions_this_proj, 1)
        #ave_im_scans = round(num_im_scans_this_proj/num_sessions_this_proj, 1)
        #
        ##data_by_proj[proj_id]['Ave. scans/session'] = ave_scans
        #data_by_proj[proj_id]['Ave. image scans/session'] = ave_im_scans
        #
        #if log_to_console:
        #    print('****************************************************\n')
    
    
    """ Round the file sizes and file counts expressed in k to 2 decimals: """
    for project in list(data_by_proj.keys()):
        for key, value in data_by_proj[project].items():
            if '[MB]' in key or '[k]' in key:
                data_by_proj[project][key] = round(value, 2)
        
    return data_by_proj