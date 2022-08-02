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
    
    parser.add_argument("-o", "--output", type=str,
                        help="""Specify location of the library. FAS_library folder can already exist in this folder.
                        If creating the library, this folder should contain the local ensembl assembly. If it does not,
                        you can download the local ensembl assembly using the -l (--local) argument. It will then be automatically
                        download into the FAS_library folder within the folder given in this argument.""")

    parser.add_argument("-s", "--species", type=str, default=None,
                        help="Specify the species.")
    
    parser.add_argument("-l", "--local", action="store_true",
                        help="Download and unpack a local ensembl assembly into the folder defined in --output.")

    args = parser.parse_args()

    output = args.output
    species = args.species
    flag_install_local = args.local

    return output, species, flag_install_local

def main():
    OUTPUT_DIR, species, flag_install_local = parser_setup()
    ensembl_access(OUTPUT_DIR, species, flag_install_local)
    
if __name__ == "__main__":
    main()