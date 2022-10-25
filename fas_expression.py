#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#######################################################################
# Copyright (C) 2022 Christian, Blümel, Julian Dosch
#
# This file is part of grand-trumpet.
#
#  grand-trumpet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  grand-trumpet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with grand-trumpet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
"""
Created on Tue Oct 25 11:00:02 2022

@author: chris
"""

import argparse

import valves.library_class as library_class
import valves.expression_extraction as ee


def parser_setup():
    """
    

    Returns
    -------
    Get paths and flags required to run.

    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config", type=str,
                        help="Path to a config file of a library. Is required for library creation.")
                        
    parser.add_argument("-n", "--namepath", type=str,
                        help="""Path to a textfile that contains the names of the expression filepath given in the argument --expressionpath.
                        One name per row.""")

    parser.add_argument("-e", "--expressionpath", type=str,
                        help="""Path to a textfile that contains the path to the expression gtf files for which the movement shall be calculated. 
                        One sample per row, though several paths can be written into one row seperated by a semicolon. Paths in the same row will be
                        concatenated and considered as several replicates of one sample.""")

                        
    args = parser.parse_args()
    config_path = args.config
    name_path = args.namepath
    expression_path = args.expressionpath

    return config_path, name_path, expression_path

def main():
    config_path, expression_path, name_path = parser_setup()
    
    fas_lib = library_class.Library(config_path, False)
    print("Expression extraction commencing...")
    ee.generate_expression_file(fas_lib, expression_path, name_path)
    
if __name__ == "__main__":
    main()