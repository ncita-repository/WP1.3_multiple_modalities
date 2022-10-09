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

    droDir = os.path.join(cwd, 'sample_DRO')
    
    # Create the directory if it doesn't exist
    if not os.path.isdir(droDir):
        #os.mkdir(ExportDir)
        Path(droDir).mkdir(parents=True)
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Download the zip file:
    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    filerequest = requests.get(url)
    
    zipFname = r'Insight Journal_923_1.zip'
    zipFpath = os.path.join(droDir, zipFname)

    with open(zipFpath, 'wb') as file:
        file.write(filerequest.content)
        print(f"Downloaded '{zipFname}' from {url}\n")

    if not file.closed:
        file.close()
        
    # Extract the zip file.
    zipFile = zipfile.zipFile(zipFpath, mode='r')
    zipFile.extractall(droDir)
    print(f"Extracted '{zipFname}' to {droDir}\n")
    
    # Uncompress the tar file:
    TarFname = 'dicom-sro.tar.gz'
    TarFpath = os.path.join(droDir, TarFname)

    with tarfile.open(TarFpath) as file:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(file, droDir)
        print(f"Extracted '{TarFname}' to {droDir}\n")

    if not file.closed:
        file.close()
    
    # The Spatial Registration Storage DRO:
    spaDroFname = '2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm'
    spaDroFpath = os.path.join(
        droDir, 'dicom-sro', 'rect-offset-sro-rigid', spaDroFname
        )

    # The Deformable Spatial Registration Storage DRO:
    defDroFname = '2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm'
    defDroFpath = os.path.join(
        droDir, 'dicom-sro', 'sphere-centered-sro-deformable', defDroFname
        )

    newSpaDroFname = 'sample_spatial_dro.dcm'
    newSpaDroFpath = os.path.join(droDir, newSpaDroFname)
    
    newDefDroFname = 'sample_deformable_dro.dcm'
    newDefDroFpath = os.path.join(droDir, newDefDroFname)
    
    # Copy the sample DROs to droDir:
    shutil.copy(spaDroFpath, newSpaDroFpath)
    shutil.copy(defDroFpath, newDefDroFpath)

    print(f"Copied sample DROs to {droDir}\n")
    
    # Delete unnecessary files and directories:
    print('Tidying up..\n')
    
    fpath = os.path.join(droDir, 'dicom-sro.tar.gz')
    delete_file(fpath)
    
    fpath = os.path.join(droDir, 'Insight Journal_923_1.zip')
    delete_file(fpath)
    
    fpath = os.path.join(droDir, 'Git_8118756.zip')
    delete_file(fpath)
    
    fpath = os.path.join(droDir, 'ITKIOTransformDCMTK.pdf')
    delete_file(fpath)
    
    dirpath = os.path.join(droDir, 'dicom-sro')       
    delete_dir(dirpath)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    print(f'\n*Took {Dtime} s to download the sample DROs.')