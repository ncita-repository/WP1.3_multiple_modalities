# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:45:41 2020

@author: ctorti
"""


"""
Function:
    ParseTransformParameters()
    
Purpose:
    Parse a SimpleElastix TransformParameters text file for the value(s) for
    a given parameter name.
        
    The text file contains data of the following format:
        
        (FixedImageDimension 2)     <-- single value
        (Size 256 256)              <-- two values
        
Input:
    TxtFilepath                   - full filepath of text file
    
    ParameterName                 - string varible of parameter name to search
    
Returns:
    CroppedTxt or ParameterValues - single string or array of strings
    
"""



def ParseTransformParameters(TxtFilepath, ParameterName):
    # Determine whether the parameters need to be converted to integers or
    # float.  
    
    # Create a list of parameter names that belong to int and float data types:
    """
    Notes: 
        1. Many of the parameters are to remain as strings.
        2. The lists below only cover the default parameters, so any additional
           parameters will need to be added.
        3. If the parameter is not found to be in ints or floats it will be 
           assumed that the parameter is a string.
    """
    
    ints = ['NumberOfParameters', 'FixedImageDimension', 'MovingImageDimension', \
            'Size', 'Index', 'FinalBSplineInterpolationOrder']
    
    floats = ['TransformParameters', 'Spacing', 'Origin', 'Direction', \
              'CenterOfRotationPoint', 'DefaultPixelValue']
    
    if ParameterName in ints:
        dtype = 'int'
    elif ParameterName in floats:
        dtype = 'float'
    else:
        dtype = 'string'
        
    
    # Read in the file:
    txt = open(TxtFilepath, 'r').read()
    
    #print('\ntxt:\n', txt)
    
    # Get index of string variable parameter:
    a = txt.index(ParameterName)
    
    #print(f'\na = {a}')
    
    """ Note:
    The text file has values encoded in the following format:
    (ParameterName ParameterValue1 ParameterValue2 ...)
    
    so the values begin at a+len(ParameterName)+1 and end at the position before ')'.
    """
    StartInd = a + len(ParameterName) + 1 # to account for the space after ParameterName
    
    #print(f'\nStartInd = {StartInd}')
    
    # Find the position of the first ')' after StartInd:
    b = StartInd + txt[StartInd:].index(')')
    
    #print(f'\nb = {b}')
    
    EndInd = b # - 1 <-- no need to subtract 1 since indexing over a range in 
    # Python is not inclusive of the stop index
    
    #print(f'\nEndInd = {EndInd}')
    
    # Crop the text between StartInd and EndInd:
    CroppedTxt = txt[StartInd:EndInd]
    
    #print(f'\nCroppedTxt = {CroppedTxt}')
    
    """ Note:
    CroppedTxt may be single value or may have multiple values. 
    e.g. (FixedImageDimension 2) <-- single value
         (Size 256 256) <-- two values
         
    If CroppedTxt doesn't contain the space character, simply return CroppedTxt.
    If CroppedTxt has a space character, slice and append the elements that preceed
    it to a new array called ParameterValues, then search for another space, repeat 
    until there are no more spaces.
    """
    # Check if there are any space characters in CroppedTxt:
    AreThereSpaces = ' ' in CroppedTxt
    
    #print('\nAreThereSpaces =', AreThereSpaces)
    
    if AreThereSpaces:
        # Initialise an array that will store all ParameterValues:
        ParameterValues = []
        
        # Continue as long as AreThereSpaces is True:
        while AreThereSpaces:
            # Find the first space character:
            ind = CroppedTxt.index(' ')

            # Slice the characters from 0 to ind (don't need to slice to ind-1 
            # since Python doesn't include the stop index):
            ParameterValues.append(CroppedTxt[0:ind])
            
            # Remove from CroppedTxt the characters that were appended to 
            # ParametersValues but not including the space character:
            CroppedTxt = CroppedTxt[ind+1:]
            
            #print(f'\n   CroppedTxt = {CroppedTxt}')
            
            # Check if there are any space characters in CroppedTxt:
            AreThereSpaces = ' ' in CroppedTxt
            
            #print('\n   AreThereSpaces =', AreThereSpaces)
            
        # There are no remaining spaces so simply append what remains of
        # CroppedTxt to ParameterValues:
        if CroppedTxt:
            ParameterValues.append(CroppedTxt)
            
        #return ParameterValues
        if dtype=='int':
            return [int(item) for item in ParameterValues]
        elif dtype=='float':
            return [float(item) for item in ParameterValues]
        else:
            return ParameterValues
        
    else:
        #return CroppedTxt
        if dtype=='int':
            return [int(item) for item in CroppedTxt]
        elif dtype=='float':
            return [float(item) for item in CroppedTxt]
        else:
            return CroppedTxt
    