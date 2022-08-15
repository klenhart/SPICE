#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 15:01:57 2022

@author: chrisbl
"""

import argparse

from ensembl_access import ensembl_access

def parser_setup():
    """
    

    Returns
    -------
    Get species, output direcotry and installation settings for run.

    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="""Specify parent directory of the library. FAS_library folder can already exist in this folder.""")

    parser.add_argument("-s", "--species", type=str, default=None,
                        help="Specify the species.")
    
    parser.add_argument("-l", "--local", action="store_true",
                        help="""Generate the directory and file structure to create a library.
                        Also creates a config tsv file that will be used by the run to quickly access
                        information about the library without having to ping ensembl.
                        This will also download a local version of ensembl.""")

    parser.add_argument("-c", "--config", type=str, default=None,
                        help="Path to a config file of a library. Is required for library creation.")

    args = parser.parse_args()
    config_path = args.config
    output = args.output
    species = args.species
    flag_install_local = args.local

    return output, species, flag_install_local, config_path

def main():
    output_dir, species, flag_install_local, config_path = parser_setup()
    ensembl_access(output_dir, species, flag_install_local, config_path)
    
if __name__ == "__main__":
    main()