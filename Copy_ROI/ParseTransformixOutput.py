# -*- coding: utf-8 -*-
"""
Created on Mon May 18 09:58:00 2020

@author: ctorti
"""



"""
Function:
    ParseTransformixOutput()
    
Purpose:
    Parse the output file "outputpoints.txt" from Transformix.
    
    The text file has the format:
    
    Point	0	; InputIndex = [ 252 359 14 ]	; InputPoint = [ 126.344037 171.830275 -54.501356 ]	; OutputIndexFixed = [ 248 353 26 ]	; OutputPoint = [ 126.638100 184.489069 3.193988 ]	; Deformation = [ 0.294064 12.658794 57.695343 ]	; OutputIndexMoving = [ 220 356 4 ]
    Point	1	; InputIndex = [ 251 360 14 ]	; InputPoint = [ 125.463303 172.711009 -54.683935 ]	; OutputIndexFixed = [ 247 354 26 ]	; OutputPoint = [ 125.763112 185.383230 3.057940 ]	; Deformation = [ 0.299810 12.672221 57.741875 ]	; OutputIndexMoving = [ 219 357 4 ]
    ...

    It will be assumed that the file is in the current working directory, and
    that outputpoints.txt was generated for 3D images.

Input:
    None
    
    
Returns:
    PtNos      - array of "Point" numbers
    
    InInds     - array of "InputIndex" indeces
    
    InPts      - array of "InputPoint" points
    
    FixOutInds - array of "OutputIndexFixed" indeces
    
    OutPts     - array of "OutputPoint" points
    
    Defs       - array of "Deformation" values
    
    MovOutInds - array of "OutputIndexMoving" indeces
    
"""


def ParseTransformixOutput():
    # Initialise arrays:
    #PtNos = []
    #InInds = []
    #InPts = []
    #FixOutInds = []
    #OutPts = []
    #Defs = []
    #MovOutInds = []
    
    # Initialise parsed data:
    Data = []
    
    # Define the strings that will be searched for:
    SearchStrings = ['Point', 'InputIndex', 'InputPoint', 
                     'OutputIndexFixed', 'OutputPoint', 
                     'Deformation', 'OutputIndexMoving']
    
    # Open the text file:
    TextFile = open('outputpoints.txt', 'r')

    # Iterate through each line in TextFile:
    for Line in TextFile:
        #print(line)
        
        # Initialise linedata:
        LineData = []
        
        # Get line into a more convenient format by replacing all 
        # '; ' and ' = ' with '',
        # '[ ' and ' ]' with '\t':
        SplitLine = Line.replace('; ', '').replace(' = ', '')\
                        .replace('[ ', '\t').replace(' ]', '\t').split('\t')
        
        
        """ Result is of the form:
        ['Point',
         '598',
         ' InputIndex = ',
         '246 367 14',
         '',
         ' InputPoint = ',
         '121.059633 178.876147 -56.073524',
         '',
         ' OutputIndexFixed = ',
         '242 362 27',
         '',
         ' OutputPoint = ',
         '122.096938 192.876032 5.728234',
         '',
         ' Deformation = ',
         '1.037305 13.999885 61.801758',
         '',
         ' OutputIndexMoving = ',
         '215 364 5',
         '']
        """
        
        #print('\nSplitLine =', SplitLine)
        
        # Append the second row (= value for 'Point'):
        #LineData.append(int(SplitLine[1]))
        
        # Initialise values for split line:
        SplitLineVals = []
        
        # The values for SearchStrings are in the row that follows the row
        # that contains those search strings. 
        # So iterate through each search string and append to LineData the 
        # values in the row beneath them:
        for SearchStr in SearchStrings:
            
            #print('\nSearchStr =', SearchStr)
            
            # The row containing the value of interest is the row containing
            # searchstr + 1:
            ind = SplitLine.index(SearchStr) + 1
            
            # The values with spaces:
            ValsWithSpaces = SplitLine[ind]
            
            # Initialise the array of values without spaces:
            ValsWithoutSpaces = []
            
            # Get the values in ValsWithSpaces by stripping away whitespaces
            # until no characters remain.  
            while ValsWithSpaces and len(ValsWithoutSpaces) < 3:
                # Strip away any whitespace on the left:
                ValsWithSpaces = ValsWithSpaces.lstrip()
                
                # Find the next whitespace:
                SpaceInd = ValsWithSpaces.find(' ')
                
                # If no whitespaces remain, simply append all characters:
                if SpaceInd == -1:
                    # Append all characters as float or int depending on 
                    # SearchStr:
                    if SearchStr in ['InputPoint', 'OutputPoint', 'Deformation']:
                        ValsWithoutSpaces.append(float(ValsWithSpaces))
                    else:
                        ValsWithoutSpaces.append(int(ValsWithSpaces))
                    
                    #print('Value appended to Vals =', ValsWithSpaces)
                    
                    # Set ValsWithSpaces to empty to return from while loop:
                    ValsWithSpaces = ''
                else:
                    # Append the characters from 0 to SpaceInd as a float or 
                    # int depending on SearchStr:
                    if SearchStr in ['InputPoint', 'OutputPoint', 'Deformation']:
                        ValsWithoutSpaces.append(float(ValsWithSpaces[0:SpaceInd]))
                    else:
                        ValsWithoutSpaces.append(int(ValsWithSpaces[0:SpaceInd]))
                        
                    #print('Value appended to Vals =', ValsWithSpaces[0:SpaceInd])
        
                    # Effectively pop ValsWithSpaces from 0 to SpaceInd+1:
                    ValsWithSpaces = ValsWithSpaces[SpaceInd+1::]
                    
                #print('ValsWithoutSpaces =', ValsWithoutSpaces)
                
            #print('ValsWithoutSpaces =', ValsWithoutSpaces)
                    
            # Append ValsWithoutSpaces to SplitLineVals:
            if len(ValsWithoutSpaces) == 1:
                SplitLineVals.append(ValsWithoutSpaces[0])
            else:
                SplitLineVals.append(ValsWithoutSpaces)
            
            #print('\nSplitLineVals =', SplitLineVals)
            
        # Append SplitLineVals to LineData:
        #LineData.append(SplitLineVals)
        LineData.extend(SplitLineVals)
        
        #print('\nLineData =', LineData)
        
        # Append LineData to Data:
        Data.append(LineData)
        
        #print('\nData =', Data)
        
        
    # Close the file:
    TextFile.close()
    
    #return Data 
    #return PtNos, InInds, InPts, FixOutInds, OutPts, Defs, MovOutInds 

        
    # Return the transpose of the array of arrays in Data:
    return [list(item) for item in list(zip(*Data))]  