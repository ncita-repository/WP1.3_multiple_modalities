# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:27:07 2021

@author: ctorti
"""


def create_session(xnatCfg):
    """Create a requests session.
    
    Parameters
    ----------
    xnatCfg : dict
        Dictionary containing the 'url', 'username', and 'password' for an
        XNAT instance. If 'username' or 'password' are empty the user will be
        prompted to enter the credentials in the console.
    
    Returns
    -------
    xnatSession : requests.models.Response
    """
    
    from getpass import getpass
    import requests
    #from pathlib import Path
    
    if xnatCfg['url'] == "":
        url = getpass("Enter url: ")
    else:
        url = xnatCfg['url']
    
    if xnatCfg['username'] == "":
        username = getpass("Enter user name: ")
    else:
        username = xnatCfg['username']
    
    if xnatCfg['password'] == "":
        password = getpass("Enter password:")
    else:
        password = xnatCfg['password']
        
    with requests.Session() as xnatSession:
        xnatSession.url = url
        xnatSession.auth = (f'{username}', f'{password}')
    
    # Test connection:
    request = requests.get(xnatSession.url, auth=xnatSession.auth)
    
    # Raise status error if not None:
    if request.raise_for_status() == None:
        print(f'Connection established to XNAT at {url}\n')
    else:
        print(request.raise_for_status())
    
    return xnatSession