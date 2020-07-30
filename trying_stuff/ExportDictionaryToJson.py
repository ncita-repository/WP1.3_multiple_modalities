# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:14:53 2020

@author: ctorti
"""

def ExportDictionaryToJson(Dictionary, FileName):
    """
    Export a dictionary to a JSON file.
    
    Inputs:
        Dictionary - A Python dictionary
        
        FileName   - The file name to assign to the .json file
        
    Returns:
        Nothing returned; a file in the current working directory
    """
    
    import json
    
    with open(FileName + '.json', 'w') as JsonFile:
        json.dump(Dictionary, JsonFile)
    
    return