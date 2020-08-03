# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:01:13 2020

@author: ctorti
"""



def CreateDir(DirPath, Debug):
    # Import packages:
    import os
    
    if Debug:
        print('\nChecking if directory', DirPath, 'exists')
    try:
        os.mkdir(DirPath)
        if Debug:
            print('\nDirectory', DirPath, 'created')
    except FileExistsError:
        if Debug:
            print('\nDirectory', DirPath, 'already exists')