# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:10:00 2021

@author: ctorti
"""


import os
from pathlib import Path
import time
#from datetime import datetime
import pandas as pd
import csv
import json
import SimpleITK as sitk


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
    items : list of strs or ints or floats
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
            if not isinstance(line, str):
                line = f'{line},\n'
            file.write(line)
    
    return

def export_list_to_csv_NOT_WORKING(items, filename, exportDir='cwd'):
    """
    NOTE:  Not sure this works.
    
    Export a list to a .csv file.
    
    Parameters
    ----------
    items : list of strs or ints or floats
        The list to be exported.
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
        
        for item in items:
            writer.writerow(item)
    
    return

def export_im(im, filename, fileFormat='HDF5ImageIO', exportDir='cwd'):
    """ 
    Export a SimpleITK Image in a desired format. 
    
    Parameters
    ----------
    im : SimpleITK Image 
        The image to be exported.
    filename : str
        The filename (without extension) to assign to the exported file.
    fileFormat : str, optional
        The desired format can be any format specified here:
            https://simpleitk.readthedocs.io/en/master/IO.html#transformations
        e.g. 'HDF5ImageIO', 'NiftiImageIO', 'NrrdImageIO'. The default value is
        'HDF5ImageIO'.
    exportDir : str, optional
        If provided the directory to which the file will be exported to. If the
        directory doesn't exist it will be created. The default value is 'cwd',
        i.e. the current working directory.
    
    Note
    ----
    The file extension will be defined based on fileFormat. Note that this part
    of this function is incomplete.
    """
    
    #if not '.nii' in filename:
    #    filename += '.nii'
    
    # TODO Complete file extension section below
    # Determine the appropriate file extension (THIS IS INCOMPLETE as it does
    # not cover all possible file formats allowed!):
    if fileFormat == 'HDF5ImageIO':
        fileExt = '.hdf'
    elif fileFormat == 'NiftiImageIO':
        fileExt = '.nii'
    elif fileFormat == 'NrrdImageIO':
        fileExt = '.nrrd'
    
    filename += fileExt
    
    if exportDir == 'cwd':
        filepath = filename
    else:
        if not os.path.isdir(exportDir):
            #os.mkdir(exportDir)
            Path(exportDir).mkdir(parents=True)
        
        filepath = os.path.join(exportDir, filename)
    
    writer = sitk.ImageFileWriter()
    writer.SetImageIO(fileFormat)
    writer.SetFileName(filepath)
    writer.Execute(im)
    
    return

def export_newRoicol(newRoicol, srcRoicolFpath, exportDir, fname=''):
    """
    Export ROI Collection (RTS/SEG) to disk.  
    
    Parameters
    ----------
    newRoicol : Pydicom object
        Target ROI Collection (RTS/SEG) to be exported.
    srcRoicolFpath : str
        Full path of the Source RTS/SEG file (used to generate the filename of
        the new RTS/SEG file).
    exportDir : str
        Directory where the new RTS/SEG is to be exported.
    fname : str
        File name to assign.     
    
    Returns
    -------
    newRoicolFpath : str
        Full path of the exported new Target RTS/SEG file.
    """
    
    # Get the filename of the original RTS/SEG file:
    #SrcRoiFname = os.path.split(SrcRoiFpath)[1]
    
    if not os.path.isdir(exportDir):
        #os.mkdir(ExportDir)
        Path(exportDir).mkdir(parents=True)
    
    currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #FnamePrefix = DateTime + '_' + NamePrefix + '_from_'
    #FnamePrefix = DateTime + '_' + TxtToAddToFname.replace(' ', '_')
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = FnamePrefix + SrcRoiFname
    #TrgRoiFname = FnamePrefix
    if fname == '':
        newRoicolFname = currentDateTime + '.dcm'
    else:
        #TrgRoiFname = Fname.replace(' ', '_') + '.dcm'
        #newRoicolFname = f'{fname}_{currentDateTime}.dcm'
        newRoicolFname = f'{fname}.dcm'
    
    newRoicolFpath = os.path.join(exportDir, newRoicolFname)
    
    newRoicol.save_as(newRoicolFpath)
        
    print(f'New Target {newRoicol.Modality} exported to:\n {newRoicolFpath}\n')
    
    return newRoicolFpath