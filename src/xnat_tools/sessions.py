# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:27:07 2021

@author: ctorti
"""


def create_session(url='', username='', aliasToken={}):
    """
    Create a requests session to an XNAT.
    
    Parameters
    ----------
    url : str, optional
        XNAT url.  The default value is ''.
    username : str, optional
        XNAT username.  The default value is ''.
    aliasToken : dict, optional
        Dictionary containing the token and secret of an XNAT alias token. 
        The default value is {}.
    
    Returns
    -------
    xnatSession : requests.models.Response
        A requests session containing the url (xnatSession.url) and 
        authentication details (xnatSession.auth = ('username', 'password')).
    """
    
    from getpass import getpass
    import requests
    
    if url == '':
        url = input('Enter XNAT url: ')
    
    # Use alias token if available:
    if aliasToken:
        username = aliasToken['alias']
        password = aliasToken['secret']
    else:
        if username == '':
            username = input('Enter XNAT user name: ')
            
        password = getpass('Enter XNAT password: ')
    
    with requests.Session() as xnatSession:
        xnatSession.url = url
        xnatSession.auth = (username, password)
    
    # Test connection:
    request = requests.get(xnatSession.url, auth=xnatSession.auth)
    
    # Raise status error if not None:
    if request.raise_for_status() == None:
        print(f"Connection established to XNAT at {url}\n")
    else:
        print(request.raise_for_status())
    
    return xnatSession