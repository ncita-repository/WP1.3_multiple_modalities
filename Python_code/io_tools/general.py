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
        print(f'{Dirpath} deleted')
    except PermissionError as error:
        print(error)
        print(f'{dirpath} cannot be deleted\n')