# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:24:41 2021

@author: ctorti
"""


""" Function to download sample DROs. """


import os
from pathlib import Path
import time
import requests
import zipfile
import tarfile
import shutil
from io_tools.general import delete_file, delete_dir



def download_sample_dros():
    """
    Download sample DICOM Registration Objects from:
    https://www.insight-journal.org/browse/publication/923
    """
    
    cwd = os.getcwd()

    DroDir = os.path.join(cwd, 'sample_DRO')
    
    # Create the directory if it doesn't exist
    if not os.path.isdir(DroDir):
        #os.mkdir(ExportDir)
        Path(DroDir).mkdir(parents=True)
    
    # Start timing:
    times = []
    times.append(time.time())
    
    """ Download the zip file:"""
    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    filerequest = requests.get(url)
    
    ZipFname = r'Insight Journal_923_1.zip'
    ZipFpath = os.path.join(DroDir, ZipFname)

    with open(ZipFpath, 'wb') as file:
        file.write(filerequest.content)
        print(f"Downloaded '{ZipFname}' from {url}\n")

    if not file.closed:
        file.close()
        
    """ Extract the zip file. """
    ZipFile = zipfile.ZipFile(ZipFpath, mode='r')
    ZipFile.extractall(DroDir)
    print(f"Extracted '{ZipFname}' to {DroDir}\n")
    
    """ Uncompress the tar file: """
    TarFname = 'dicom-sro.tar.gz'
    TarFpath = os.path.join(DroDir, TarFname)

    with tarfile.open(TarFpath) as file:
        file.extractall(DroDir)
        print(f"Extracted '{TarFname}' to {DroDir}\n")

    if not file.closed:
        file.close()
    
    # The Spatial Registration Storage DRO:
    SpaDroFname = '2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm'
    SpaDroFpath = os.path.join(DroDir, 'dicom-sro', 'rect-offset-sro-rigid', 
                               SpaDroFname)

    # The Deformable Spatial Registration Storage DRO:
    DefDroFname = '2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm'
    DefDroFpath = os.path.join(DroDir, 'dicom-sro', 
                               'sphere-centered-sro-deformable', DefDroFname)

    NewSpaDroFname = 'sample_spatial_dro.dcm'
    NewSpaDroFpath = os.path.join(DroDir, NewSpaDroFname)
    
    NewDefDroFname = 'sample_deformable_dro.dcm'
    NewDefDroFpath = os.path.join(DroDir, NewDefDroFname)
    
    """ Copy the sample DROs to DroDir: """
    shutil.copy(SpaDroFpath, NewSpaDroFpath)
    shutil.copy(DefDroFpath, NewDefDroFpath)

    print(f"Copied sample DROs to {DroDir}\n")
    
    """ Delete unnecessary files and directories: """
    print('Tidying up..\n')
    
    fpath = os.path.join(DroDir, 'dicom-sro.tar.gz')
    delete_file(fpath)
    
    fpath = os.path.join(DroDir, 'Insight Journal_923_1.zip')
    delete_file(fpath)
    
    fpath = os.path.join(DroDir, 'Git_8118756.zip')
    delete_file(fpath)
    
    fpath = os.path.join(DroDir, 'ITKIOTransformDCMTK.pdf')
    delete_file(fpath)
    
    dirpath = os.path.join(DroDir, 'dicom-sro')       
    delete_dir(dirpath)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    print(f'\n*Took {Dtime} s to download the sample DROs.')