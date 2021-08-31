# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:13:20 2021

@author: ctorti
"""

import os
import shutil

def delete_file(filepath):
    try:
        os.remove(filepath)
        print(f'{filepath} deleted')
    except PermissionError as error:
        print(error)
        print(f'{filepath} cannot be deleted\n')

def delete_dir(dirpath):
    try:
        shutil.rmtree(dirpath)
        print(f'{dirpath} deleted')
    except PermissionError as error:
        print(error)
        print(f'{dirpath} cannot be deleted\n')

def get_user_input_as_int(message="Enter an integer", minVal=None, 
                          maxVal=None):
    """
    Get user to input an integer.
    
    Optional minimum and maximum allowable values. If input is not int or not
    within range the user will continue to be prompted to re-enter.
    
    Parameters
    ----------
    message : str, optional
        Message to present to user. The default is "Enter an integer: ".
    minVal : int, optional
        Minimum allowable value. The default is None.
    maxVal : int, optional
        Maximum allowable value. The default is None.
        
    Returns
    -------
    choice : int
        The user input.
    """
    
    if minVal and maxVal:
        message += f' ({minVal} <= integer <= {maxVal}; q to quit): '
    elif minVal:
        message += f' (>= {minVal}; q to quit): '
    elif maxVal:
        message += f' (<= {maxVal}; q to quit): '
    else:
        message += '(q to quit) : '
        
    helpfulMsg = str(message)
    
    while True:
        userInput = input(helpfulMsg)
        
        if userInput == 'q':
            raise Exception("Quit by the user.")
        
        try:
            choice = int(userInput)
        except ValueError:
            helpfulMsg = "\nYou did not enter an integer.\n" + message
            
            continue
        else:
            # userInput was string of int...
            if minVal <= choice <= maxVal:
                # ...and is within the allowed range so exit the loop
                break
            else:
                # ...but is not within the allowed range, so continue looping
                helpfulMsg = f"The integer must be >= {minVal} and <= {maxVal}.\nTry again.\n" + message
                continue
    
    return choice

def get_fpaths(dpath, ext="", includeSubDirs=False):
    """
    NOTE: This function is less efficient than the new one, which uses 
    SimpleITK's ImageSeriesReader.
    
    Get a sorted list of full filepaths for all DICOM files in a directory
    
    Inputs
    ------
    dpath : str
        Path to a file directory.
    ext : str, optional 
        File extension to filter returned list of file paths. The default value
        is "".
    includeSubDirs : bool, optional
        If True the file searching will include sub-directories. Default value
        is False.
    
    Returns
    -------
    fpaths : list of strs
        List of the full filepaths of all files in dpath that meet the criteria.
    """
    
    fpaths = []
    
    if includeSubDirs:
        # Use os.walk to get list of file paths of files at any depth
        for dirName, subDirList, itemList in os.walk(dpath):
            for item in itemList:
                fpath = os.path.join(dirName, item)
                
                isFile = os.path.isfile(fpath)
                
                # Append fpath to fpaths if the item is a file and ext = ""
                # or ext != "" and it is the correct extension:
                if (isFile and (ext == "" or (ext and ext in item.lower()))):
                    fpaths.append(fpath)
    else:
        itemList = os.listdir(dpath)
        
        for item in itemList:
            fpath = os.path.join(dpath, item)
            
            isFile = os.path.isfile(fpath)
            
            if (isFile and (ext == "" or (ext and ext in item.lower()))):
                fpaths.append(fpath) 
    
    return fpaths