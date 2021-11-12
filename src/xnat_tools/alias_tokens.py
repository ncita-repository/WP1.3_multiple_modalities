# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:25:39 2021

@author: ctorti
"""

import requests
import os
from io_tools.exports import export_dict_to_json
from io_tools.imports import import_dict_from_json

def generate_alias_token(xnatSession):
    """
    Generate an XNAT alias token. 
    
    Parameters
    ----------
    xnatSession : requests.models.Response
        A requests session containing the url (xnatSession.url) and 
        authentication details (xnatSession.auth = ('username', 'password')).
    
    Returns
    -------
    aliasToken : dict
        Dictionary containing an alias token.
    
    Note
    ----
    https://wiki.xnat.org/display/XAPI/User+Alias+Token+API
    """
    
    uri = xnatSession.url + "/data/services/tokens/issue"
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    return request.json()

def export_alias_token(aliasToken, exportDir):
    """
    Export an XNAT alias token to the local drive. 
    
    Parameters
    ----------
    aliasToken : dict
        Dictionary containing an alias token.
    exportDir : str
        Directory to export alias token to.
    
    Returns
    -------
    None.
    """
    
    export_dict_to_json(aliasToken, 'XNAT_alias_token', exportDir)
    
    print(f'XNAT alias token exported to {exportDir}\n')

def import_alias_token(exportDir):
    """
    Import an XNAT alias token from the local drive. 
    
    Parameters
    ----------
    exportDir : str
        Directory containing the previously exported alias token.
    
    Returns
    -------
    aliasToken : dict
        Dictionary containing an alias token. Empty if the file was not
        found.
    """
    
    fpath = os.path.join(exportDir, 'XNAT_alias_token.json')
    
    try:
        aliasToken = import_dict_from_json(filepath=fpath)
        
        print('An XNAT alias token was found.')
    
    except FileNotFoundError:
        aliasToken = {}
        
        print('An XNAT alias token was not found.')
    
    return aliasToken

def validate_alias_token(aliasToken, xnatSession):
    """
    Check if an alias token is valid for an existing XNAT session, and return 
    a dictionary containing the user name assigned to the token (if valid), or
    an empty dictionary (if invalid).
    
    Parameters
    ----------
    aliasToken : dict
        Dictionary containing an XNAT alias token.
    xnatSession : requests.models.Response
        A requests session containing the url (xnatSession.url) and 
        authentication details (xnatSession.auth = ('username', 'password')).
    
    Returns
    -------
    validDict : dict
        If valid the dictionary will contain the key:value pair 
        'valid':username, where username is the user name assigned to the
        token.  If invalid, validDict will be empty.
    """
    
    token = aliasToken['alias']
    secret = aliasToken['secret']
    
    uri = xnatSession.url + f"/data/services/tokens/validate/{token}/{secret}"
    
    request = xnatSession.get(uri)
    
    # Raise status error if not None:
    if request.raise_for_status() != None:
        print(request.raise_for_status())
    
    return request.json()

def is_alias_token_valid(aliasToken, url):
    """
    Check if an alias token is valid and return a True if valid.
    
    Parameters
    ----------
    aliasToken : dict
        Dictionary containing an alias token.
    url : str
        XNAT url.
    
    Returns
    -------
    isValid : bool
        True if aliasToken is valid.
    """
    
    if aliasToken == {}:
        return False
    else:
        token = aliasToken['alias']
        secret = aliasToken['secret']
        
        with requests.Session() as xnatSession:
            xnatSession.url = url
            xnatSession.auth = (token, secret)
        
        # Test connection:
        request = requests.get(xnatSession.url, auth=xnatSession.auth)
        
        try:
            result = request.raise_for_status()
            
            if result == None:
                return True
        except:
            return False