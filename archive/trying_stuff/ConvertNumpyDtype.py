# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:49:17 2020

@author: ctorti
"""


"""
Function:
    ConvertNumpyDtype()
    
Purpose:
    Convert numpy Array from one data type to another, preserving the
    dynamic range.

Input:
    ArrayIn  - Numpy Array whose data type is to be converted; acceptable data
               types include:
                  -- uint8
                  -- uint16
                  -- int8
                  -- int16
                  -- float
                        
    ConvertTo - Data type to convert ArrayIn to
    
    Debug     - If True, some results will be printed

    
Returns:
    ArrayOut  - Numpy Array with data type determined by ConvertTo   
    

Note:

At the moment this function only handles a limited number of possible cases of
"from-to" data types: uint8, uint16, int8, int16 and float.                 
                  
"""

def ConvertNumpyDtype(ArrayIn, ConvertTo, Debug):
    # Import packages:
    import numpy as np
    
    #""" April 8: Temporarily set Debug = True: """
    #Debug = True
    
    # Get info on input Array, first assuming that it is an int dtype:
    try:
        info = np.iinfo(ArrayIn.dtype)
        
    # If ValueError arises, assume it is a float:
    except ValueError:
        info = np.finfo(ArrayIn.dtype)
    
    if Debug:    
        print('\nInput dtype:', ArrayIn.dtype)
        print(f'\nInput: {info}')
        print('\nConvertTo = ', ConvertTo)
      
    
    # Determine value to scale Array up to:
    if 'uint8' in ConvertTo:
        MaxValue = 255
    
    if 'uint16' in ConvertTo:
        MaxValue = 65535
        
    if 'int8' in ConvertTo:
        MaxValue = 127
    
    if 'int16' in ConvertTo:
        MaxValue = 32767
        
    if 'float' in ConvertTo:
        MaxValue = 3.4028235e+38
        
    if Debug:    
        print('\nMax value possible for ' + ConvertTo + f' = {MaxValue}')
    
    # Get min/max values:
    ArrayMin = np.min(ArrayIn)
    ArrayMax = np.max(ArrayIn)
    
    if Debug:
        print(f'\nMin/Max of input Array = [{ArrayMin}, {ArrayMax}]')
    
    #if ArrayMin < 0:
    
    # Deal with negative values if ConvertTo is an unsigned integer
    # (by subtracting the min value):
    if 'u' in ConvertTo:
        Array = ArrayIn - ArrayMin

        # Get min/max values:
        ArrayMin = np.min(Array)
        ArrayMax = np.max(Array)
        
        if Debug:
            print(f'\nMin/Max of Array after subtracting Min = [{ArrayMin}, {ArrayMax}]')
        
    else:
        Array = ArrayIn - 0
    
    # Normalise Array to [0, 1] or [-1, 1] (depends on whether
    # the input Array is signed or unsigned):
    #Array = ArrayIn.astype(np.float64) / info.max # is conversion to float useful?
    #Array = ArrayIn / info.max 
    Array = Array / MaxValue
    
    # Get min/max values:
    ArrayMin = np.min(Array)
    ArrayMax = np.max(Array)
    
    if Debug:
        print(f'\nMin/Max of Array after normalising = [{ArrayMin}, {ArrayMax}]')
        
    
    # Scale Array by MaxValue:
    Array = MaxValue * Array
    
    # Get min/max values:
    ArrayMin = np.min(Array)
    ArrayMax = np.max(Array)
    
    if Debug:
        print(f'\nMin/Max of Array after scaling by {MaxValue} = [{ArrayMin}, {ArrayMax}]')
    
    # Convert Array to desired dtype:
    if 'uint8' in ConvertTo:
        if Debug:
            print('\nConverting to uint8')
        ArrayOut = Array.astype(np.uint8)
        
    if 'uint16' in ConvertTo:
        if Debug:
            print('\nConverting to uint16')
        ArrayOut = Array.astype(np.uint16)
        
    if 'int8' in ConvertTo:
        if Debug:
            print('\nConverting to int8')
        ArrayOut = Array.astype(np.int8)
        
    if 'int16' in ConvertTo:
        if Debug:
            print('\nConverting to int16')
        ArrayOut = Array.astype(np.int16)
        
    if 'float' in ConvertTo:
        if Debug:
            print('\nConverting to float')
        ArrayOut = Array.astype(np.float)
        
        
    # For some reason the attempts to convert to 
        
    
    # Get min/max values:
    ArrayMin = np.min(ArrayOut)
    ArrayMax = np.max(ArrayOut)
    
    if Debug:
            print('\nMin/Max of Array after converting to', ConvertTo, \
                  f'= [{ArrayMin}, {ArrayMax}]')
    
    # Re-check info:
    try:
        info = np.iinfo(ArrayOut.dtype)
        
    # If ValueError arises, assume it is a float:
    except ValueError:
        info = np.finfo(ArrayOut.dtype)
        
    if Debug:
            print('\nOutput dtype:', ArrayOut.dtype)
    if Debug:
            print(f'\nOutput: {info}')
    
    return ArrayOut
