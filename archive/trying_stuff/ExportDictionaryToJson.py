# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:14:53 2020

@author: ctorti
"""

def ExportDictionaryToJson(Dictionary, FileName):
    """
    Export a dictionary to a JSON file to the current working directory.
    
    Inputs:
        Dictionary - A Python dictionary
        
        FileName   - The file name to assign to the exported file. If FileName
                     does not include the '.json' extension it will be added
                     automatically.
        
    Returns:
        Nothing returned; 
        The file will be saved in the current working directory
    """
    
    import json
    
    if not '.json' in FileName:
        FileName += '.json'
        
    
    with open(FileName, 'w') as JsonFile:
        json.dump(Dictionary, JsonFile)
    
    return