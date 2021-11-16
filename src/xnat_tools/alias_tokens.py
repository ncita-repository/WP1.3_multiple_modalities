# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:25:39 2021

@author: ctorti
"""

import requests
import os
import time
from io_tools.exports import export_dict_to_json
from io_tools.imports import import_dict_from_json
from general_tools.general import get_list_of_filePaths

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
    
    # Current datetime
    cdt = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    fname = f'xnat_alias_token_{cdt}'
    
    export_dict_to_json(aliasToken, fname, exportDir)
    
    print(f'XNAT alias token exported to {os.path.join(exportDir, fname)}\n')

def import_alias_token(tokensDir):
    """
    Import an XNAT alias token from the local drive, returning {} if not
    found.
    
    Parameters
    ----------
    tokensDir : str
        Directory containing XNAT alias tokens.
    
    Returns
    -------
    aliasToken : dict
        Dictionary containing an alias token. Empty if no token was found.
    """
    
    fpath = os.path.join(tokensDir, 'XNAT_alias_token.json')
    
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

def is_dict_an_alias_token(possibleToken):
    """
    Check if a dictionary contains keys expected of an XNAT alias token.
    
    Parameters
    ----------
    possibleToken : dict
        A possible XNAT alias token.
    
    Returns
    -------
    is_alias_token : bool
        True if possibleToken is an XNAT alias token, False otherwise.
    """
    
    is_alias_token = False
    
    # Could add more to the list but two should be sufficient:
    keysToMatch = ['alias', 'secret']
    
    keys = possibleToken.keys()
    
    numOfMatches = 0
    
    for keyToMatch in keysToMatch:
        if keyToMatch in keys:
            numOfMatches += 1
    
    if numOfMatches == len(keysToMatch):
        is_alias_token = True
    
    return is_alias_token
    
def import_valid_alias_token(tokensDir, url):
    """
    Import a valid XNAT alias token from the local drive, returning {} if not
    found.
    
    Parameters
    ----------
    tokensDir : str
        Directory containing XNAT alias tokens.
    url : str
        XNAT url.
    
    Returns
    -------
    aliasToken : dict
        Dictionary containing an alias token. Empty if no valid token was 
        found.
    """
    
    aliasToken = {}
    
    # Get list of JSON files in tokensDir:
    fpaths = get_list_of_filePaths(tokensDir, fileExt="json")
    
    if fpaths:
        for fpath in fpaths:
            possibleToken = import_dict_from_json(filepath=fpath)
            
            if is_dict_an_alias_token(possibleToken):
                if is_alias_token_valid(possibleToken, url):
                    print('A valid XNAT alias token was found.')
                    
                    aliasToken = dict(possibleToken)
                    
                    continue
    
    if not aliasToken:
        print('A valid XNAT alias token was not found.')
    
    return aliasToken