# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:12:04 2021

@author: ctorti
"""

import json


def import_dict_from_json(filepath):
    """
    Import a dictionary from a JSON file from the current working directory.
    
    Parameters
    ----------
    filepath : str
        The file path of the JSON file (which must be in the current working 
        directory).  filepath does not need to include the .json extension.
        
    Returns
    -------
    out : dict
    """
    
    if not '.json' in filepath:
        filepath += '.json'
    
    with open(filepath, 'r') as file:
        out = json.load(file)
    
    return out

def import_list_from_txt(filepath):
    """
    Import a list from a text file from the current working directory.
    
    Parameters
    ----------   
    filepath : str
        The file path of the text file (which must be in the current working 
        directory).  filepath does not need to include the .txt extension.
    
    Returns
    -------  
    out : list
    """
    
    if not '.txt' in filepath:
        filepath += '.txt'
    
    out = []
    
    with open(filepath, 'r') as file:
        for line in file:
            out.append(line.replace('\n', '').split(' '))
    
    return out
