# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:10:00 2021

@author: ctorti
"""


import os
from pathlib import Path
import pandas as pd
import csv
import json


def export_dict_to_csv(dictionary, filename, exportDir='cwd'):
    """
    NOTE:  Not sure this works.
    
    Export a dictionary to a .csv file.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to be exported.
    filename : str
        The file name to assign to the exported file. If filename doesn't 
        include the '.csv' extension it will be added automatically.
    exportDir : str, optional ('cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    if not '.csv' in filename:
        filename += '.csv'
    
    if exportDir == 'cwd':
        filepath = filename
    else:
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        filepath = os.path.join(exportDir, filename)
    
    with open(filepath, 'wb') as output:
        writer = csv.writer(output)
        
        for key, value in dictionary.items():
            print(f"{key} = {value}")
            writer.writerow([key, value])
    
    return

def dict_to_dataframe(dictionary):
    """
    Convert a dictionary to a Pandas dataframe.  
    
    Dictionary assumed to have a simple {key : value} structure.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to be exported.
    
    Returns
    -------
    df : Pandas Dataframe
    """
    
    df = pd.DataFrame.from_dict(data=dictionary, orient='index', columns=[''])
    
    return df

def dict_of_dicts_to_dataframe(dictionary, nameOfFirstCol=''):
    """
    Convert a dictionary of dictionaries to a Pandas dataframe.  
    
    Dictionary assumed to have a nested {key : {key : value}} structure.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to be exported.
    nameOfFirstCol : str, optional ('' by default)
        The name to assign to the first column of the spreadsheet (i.e. the
        first column of the pandas dataframe)
    
    Returns
    -------
    df : Pandas Dataframe
    """
    
    """
    #df = pd.DataFrame(data=dictionary)
    
    # Transpose to swap columns for rows, and remove apostrophes from strings
    # (e.g. lists of strings in the 'Co-investigators' column): 
    #df = (df.T).replace('\'','', regex=True, inplace=True)
    
    df = df.T # <-- pd.DataFrame.from_dict(data=dictionary, orient='index') achieves this
    
    #df = df.replace('\'','', regex=True, inplace=True) # wrong/errors
    
    #df['Co-investigators'].replace('\'','', regex=True, inplace=True) # doesn't work
    
    #df['Co-investigators'] = df['Co-investigators'].replace('\'','', regex=True, inplace=True) # deletes everything in column!
    
    #df['Co-investigators'] = df['Co-investigators'].str.replace(r"[\"\',]", '')
    
    #for key in df['Co-investigators'].keys():
    #    [item.replace(r"[\"\',]", '') for item in df['Co-investigators'][key]]
    """
    
    df = pd.DataFrame.from_dict(data=dictionary, orient='index')
    
    
    """ Add the name nameOfFirstCol to the first column? """
    if nameOfFirstCol != '':
        df.index.name = nameOfFirstCol
    
    return df

def export_dict_to_xlsx(dictionary, filename, exportDir='cwd'):
    """
    Export a dictionary to a .xlsx file.  
    
    Dictionary assumed to have a simple {key : value} structure.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to be exported.
    filename : str
        The file name to assign to the exported file. If filename doesn't 
        include the '.xlsx' extension it will be added automatically.
    exportDir : str, optional ('cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    if not '.xlsx' in filename:
        filename += '.xlsx'
    
    if exportDir == 'cwd':
        filepath = filename
    else:
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        filepath = os.path.join(exportDir, filename)
    
    df = dict_to_dataframe(dictionary)
    
    df.to_excel(filepath)
    
    return

def export_dict_of_dicts_to_xlsx(dictionary, filename, exportDir='cwd', 
                                 nameOfFirstCol=''):
    """
    Export a dictionary of dictionaries to a .xlsx file.  
    
    Dictionary assumed to have a nested {key : {key : value}} structure.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to be exported.
    filename : str
        The file name to assign to the exported file. If filename doesn't 
        include the '.xlsx' extension it will be added automatically.
    exportDir : str, optional ('cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    nameOfFirstCol : str, optional ('' by default)
        The name to assign to the first column of the spreadsheet (i.e. the
        first column of the pandas dataframe)
    """
    
    if not '.xlsx' in filename:
        filename += '.xlsx'
    
    if exportDir == 'cwd':
        filepath = filename
    else:
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        filepath = os.path.join(exportDir, filename)
    
    df = dict_of_dicts_to_dataframe(dictionary, nameOfFirstCol)
    
    df.to_excel(filepath)
    
    return

def export_dict_to_json(dictionary, filename=None, exportDir='cwd'):
    """ 
    Export a dictionary to a .json file.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to be exported.
    filename : str, optional (None by default)
        The file name to assign to the exported file. If filename doesn't 
        include the '.json' extension it will be added automatically. If a
        file name is not provided, a file name will be generated using the
        current date-time followed by '_Dict2Json'.
    exportDir : str, optional ('cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    if filename == None:
        import time
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            
        filename = currentDateTime + '_Dict2Json.json'
    else:
        if not '.json' in filename:
            filename += '.json'
    
    if exportDir == 'cwd':
        filepath = filename
    else:
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        filepath = os.path.join(exportDir, filename)
    
    with open(filepath, 'w') as file:
        json.dump(dictionary, file)
    
    return

def export_list_to_txt(items, filename, exportDir='cwd'):
    """
    Export a list to a .txt file.
    
    Parameters
    ----------
    items : list of strs
        The list to be exported.
    filename : str
        The file name to assign to the exported file. If filename doesn't 
        include the '.txt' extension it will be added automatically.
    exportDir : str, optional ('cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    if not '.txt' in filename:
        filename += '.txt'
    
    if exportDir == 'cwd':
        filepath = filename
    else:
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        filepath = os.path.join(exportDir, filename)
    
    with open(filepath, 'w') as file:
        for line in items:
            file.write(line)
    
    return