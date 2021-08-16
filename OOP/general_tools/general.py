# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:38:32 2021

@author: ctorti
"""

import os
from copy import deepcopy
import numpy as np

def flatten_list(listOfLists):
    """
    Flatten a nested list of lists.
    
    There may be any number of tiers within the nested list.

    Parameters
    ----------
    listOfLists : nested list of int/float/etc.
        Nested list of items.

    Returns
    -------
    flatList : list of int/floats
        A flattened list of items.
    """
    
    if len(listOfLists) == 0:
        return listOfLists
    if isinstance(listOfLists[0], list):
        return flatten_list(listOfLists[0]) + flatten_list(listOfLists[1:])
    return listOfLists[:1] + flatten_list(listOfLists[1:])

def unpack(items):
    """
    Unpack list of 2D or 3D lists (e.g. [x, y] or [x, y, z] coordinates, or 
    [i, j] or [i, j, k] indices) into separate lists for each dimension.
    
    Parameters
    ----------
    items : list of floats/indices
        List of a list of items, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
     
    Returns
    -------
    X : list of floats/indices
        List of all x coordinates in items, e.g. [x0, x1, x2, ...].
    
    Y : list of floats/indices
        List of all y coordinates in items, e.g. [y0, y1, y2, ...].
    
    Z : list of floats/indices
        List of all z coordinates in items, e.g. [z0, z1, z2, ...].
    """
    
    nda = np.array(items)
    
    if len(nda.shape) < 2:
        msg = "The input argument is a list of length 1. There is nothing to" \
              + " unpack."
        
        raise Exception(msg)
        
    dim = nda.shape[1]
    
    # Initialise unpacked lists:
    X = []
    Y = []
    Z = []
    
    if dim == 2:
        # Unpack tuple and store in X, Y:
        for x, y in items:
            X.append(x)
            Y.append(y)
            
        return X, Y
    
    elif dim == 3:
        # Unpack tuple and store in X, Y, Z:
        for x, y, z in items:
            X.append(x)
            Y.append(y)
            Z.append(z)
            
        return X, Y, Z
    
    else:
        msg = f"The input argument has dim = {dim}. Only dim = 2 or dim = 3 "\
              + "is allowed."
        
        raise Exception(msg)

def get_unique_items(items, ignoreZero=False, maintainOrder=False):
    """
    Get list of unique items in a list or Numpy data array.
    
    Parameters
    ----------
    items : list of int/float or a Numpy data array
    ignoreZero : bool, optional (False by default)
        If ignoreZero=True only unique non-zero items will be returned.
    maintainOrder : bool, optional (False by default)
        If maintainOrder=True the unique items will be returned in their
        original order.
        
    Returns
    -------
    uniqueItems : list of int/float or Numpy data array or None
        List or Numpy array of unique items, or None if items is empty.
    """
    
    if not isinstance(items, list) and not isinstance(items, np.ndarray):
        msg = f'The input "items" is data type {type(items)}.  Acceptable '\
              + 'data types are "list" and "numpy.ndarray".'
        raise Exception(msg)
    
    origItems = deepcopy(items)
    
    if isinstance(items, list):
        #items = np.array(items)
        ##items = list(np.concatenate(items).flat)
        #items = np.concatenate(items)
        
        items = np.array(flatten_list(items))
    
    if items.size == 0:
        return None
    
    if ignoreZero:
        inds = np.nonzero(items)
        
        items = [items[ind] for ind in inds[0]]
    
    if maintainOrder:
        uniqueItems = []
        
        for item in items:
            if not item in uniqueItems:
                uniqueItems.append(item)
    else:
        #uniqueItems = list(set(items.tolist()))
        uniqueItems = np.unique(items)
    
    # If items was a list convert back to a list:
    if isinstance(origItems, list):
        uniqueItems = list(uniqueItems) # 18/12/2020
        #uniqueItems = sorted(list(uniqueItems)) # 18/12/2020
    
    return uniqueItems

def get_items_unique_to_within(items, epsilon=1e-5):
    """
    Return a list of items unique to within epsilon.

    Parameters
    ----------
    items : list of float or int
        List of items.
    epsilon : float, optional
        Maximum allowed variation between items not considered to be unique. 
        The default is 1e-5.

    Returns
    -------
    uniqueItems : list of float or int
        List of items that are unique to within epsilon.
    """
    
    uniqueItems = []
    
    n = len(items)
    
    if n < 2:
        print(f'The list has {n} items (it must have at least 2 items).')
        
        return uniqueItems
    
    else:
        uniqueItems.append(items[0])
        
        for i in range(n - 1):
            if abs(items[i] - items[0]) > epsilon:
                uniqueItems.append(items[i])
                
        return uniqueItems

def are_items_equal_to_within_eps_REDUNDANT(item0, item1, epsilon=1e-06):
    """
    07/07/21:
        This is redundant - use are_lists_equal_to_within_eps(), which was
        latter named are_items_equal_to_within_eps.
        
    Are two items equal to within epsilon?
    
    Parameters
    ----------
    item0 : float or int
    item1 : float or int
    epsilon : float, optional
        Maximum allowed variation between items not considered to be unique. 
        The default is 1e-5.

    Returns
    -------
    is_equal : bool
        True if the items are equal to within epsilon.
    """
    
    abs_diff = abs(item0 - item1)
    
    if abs_diff < epsilon:
        is_equal = True
    else:
        is_equal = False
        
    return is_equal

def are_items_equal_to_within_eps(items0, items1, epsilon=1e-06):
    """
    08/07/21 Note: Formerly called are_lists_equal_to_within_eps.
    
    Are two items equal to within epsilon?
    
    The items may be ints, floats, or lists of ints or floats.
    
    Parameters
    ----------
    items0 : list of floats or ints
    items1 : list of floats or ints
    epsilon : float, optional
        Maximum allowed variation between items not considered to be unique. 
        The default is 1e-5.

    Returns
    -------
    are_equal : bool
        True if the items are equal to within epsilon.
    """
    
    # Get the maximum value of their absolute differences:
    abs_diffs = abs(np.array(items0) - np.array(items1))
    max_abs_diffs = max(abs_diffs)
    
    if max_abs_diffs < epsilon:
        are_equal = True
    else:
        are_equal = False
        
    return are_equal

def get_list_of_dtypes(items, unique=True):
    """
    Return a list of the data type of the elements in the list.  
    
    If the parameter unique is True or omitted, the set of data types will be 
    returned.  Hence a single-itemed list containing a single data type is a 
    permitted output.
    
    Parameters
    ----------
    items : list of items
        A list of items (e.g. integers, floats, strings, etc.).
    unique : bool, optional (True by default)
        If True the set of data types will be returned (i.e. only unique data
        types).
        
    Returns
    -------
    dtypes : list of data types
        A list of the data type of each item in ListOfItems.
    """
    
    dtypes = [type(item) for item in items]
    
    if unique:
        dtypes = list(set(dtypes))
    
    return dtypes

def reduce_list_of_str_floats_to_16(origList):
    """ 
    Reduce a list of string floats to 16 characters.
    
    Decimal strings (VR DS in DICOM standard) must be limited to 16 characters. 
    
    Parameters
    ----------  
    origList : list of strs
        The original list of string floats.
    
    Returns
    ------- 
    newList : list of strs
        The list of strings with characters limited to 16.
    """
    
    newList = [item[:16] for item in origList]
    
    return newList

def check_fileExt(filePath, fileExt):
    """
    Check if a filePath has a specified extension.

    Parameters
    ----------
    filePath : str
        Filepath.
    fileExt : str
        File extension.

    Returns
    -------
    isMatch : bool
        True if the specified extension matches.
    """
    
    isMatch = False
    
    if os.path.exists(filePath):
        # is a file
        
        # Split the extension from the path and normalise it to lowercase.
        ext = os.path.splitext(filePath)[-1].lower()

        if ext == fileExt:
            isMatch = False
            
    return isMatch

def does_instance_variable_exist(instanceObj, varName, p2c=False):
    # TODO update docstrings
    """
    Check if an instance variable exists within an instance.
    
    Parameters
    ----------
    instanceObj : instance object
        Instance object.
    varName : str
        Name of the instance variable to check existance in instanceObj.
    p2c : bool
        If True result will be printed to the console.
        
    Returns
    -------
    doesExist : bool
        True if the instance variable exists, False otherwise.
    """
    
    if varName in dir(instanceObj):
        doesExist = True
        
        if p2c:
            print(f'Instance variable {varName} exists.\n')
    else:
        doesExist = False
        
        if p2c:
            print(f'Instance variable {varName} does not exist.\n')
    
    return doesExist