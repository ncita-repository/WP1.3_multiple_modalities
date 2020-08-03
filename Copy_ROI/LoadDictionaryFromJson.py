# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:16:01 2020

@author: ctorti
"""

def LoadDictionaryFromJson(FileName):
    """
    Export a dictionary to a JSON file.
    
    Input:
        FileName   - The file name of the .json file in the current working 
                     directory
        
    Returns:
        Dictionary - A Python dictionary
    """
    
    import json
    
    with open(FileName + '.json', 'r') as JsonFile:
        Dictionary = json.load(JsonFile)
    
    return Dictionary