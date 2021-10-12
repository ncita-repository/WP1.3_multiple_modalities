# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:16:01 2020

@author: ctorti
"""

def ImportDictionaryFromJson(FileName):
    """
    Import a dictionary from a JSON file in the current working directory.
    
    Input:
        FileName   - The file name of the JSON file (which must be in the 
                     current working directory).  If FileName does not include
                     the .json extension it will be added automatically.
        
    Returns:
        Dictionary - A Python dictionary
    """
    
    import json
    
    if not '.json' in FileName:
        FileName += '.json'
    
    with open(FileName, 'r') as JsonFile:
        Dictionary = json.load(JsonFile)
    
    return Dictionary