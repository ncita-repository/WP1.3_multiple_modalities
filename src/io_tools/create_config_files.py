# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 14:37:02 2021

@author: ctorti
"""


""" Generate new config files """


from importlib import reload
import general_tools.config
reload(general_tools.config)

import sys
import argparse
from general_tools.config import create_config_files


code_root = r'C:\Code\WP1.3_multiple_modalities\OOP\io_tools'

# Add code_root to the system path so packages can be run from it:
sys.path.append(code_root)


def main(cfgDir, password=''):
    """
    Generate new configuration files.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    password : str, optional
        If provided the password will be used when making REST API calls to
        XNAT. The default is "".
    
    Returns
    -------
    None.
    """
    
    create_config_files(cfgDir, password)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments for main()')
    parser.add_argument(
        "--cfgDir", help="The directory containing config files"
        )
    parser.add_argument(
        "--password", action="store_str", str="", 
        help="Password for XNAT connection"
        )
    args = parser.parse_args()
    
    main(args.name)