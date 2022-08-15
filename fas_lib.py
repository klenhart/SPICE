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
    
    parser.add_argument("-m", "--movement", action="store_true",
                        help="""Calculate the movement over all transcripts in an expression GTF. 
                        If this argument is used, the name_path and expression_path arguments are required.""")
                        
    parser.add_argument("-n", "--namepath", type=str, default=None,
                        help="""Path to a textfile that contains the names of the expression filepath given in the argument --expressionpath.
                        One name per row.""")

    parser.add_argument("-e", "--expressionpath", type=str, default=None,
                        help="""Path to a textfile that contains the path to the expression gtf files for which the movement shall be calculated. 
                        One sample per row, though several paths can be written into one row seperated by a semicolon. Paths in the same row will be
                        concatenated and considered as several replicates of one sample.""")
                        
    args = parser.parse_args()
    config_path = args.config
    output = args.output
    species = args.species
    name_path = args.namepath
    expression_path = args.expressionpath
    flag_install_local = args.local
    flag_movement = args.movement

    return output, species, flag_install_local, config_path, flag_movement, name_path, expression_path

def main():
    output_dir, species, flag_install_local, config_path,  flag_movement, name_path, expression_path = parser_setup()
    ensembl_access(output_dir, species, flag_install_local, config_path, flag_movement, name_path, expression_path)
    
if __name__ == "__main__":
    main()