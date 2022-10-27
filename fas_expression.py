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
import os
import json

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
                        help="Path to a config file of a FAS library.")
                        
    parser.add_argument("-n", "--name", type=str,
                        help="""The name of the condition for which the expression shall be extracted.
                        Output file will be named after this.""")

    parser.add_argument("-e", "--expressionpath", nargs="+", action="append",
                        help="""Path to the expression gtf files for which the expression shall be extracted into the results directory. 
                        If given several paths, they will be considered several replicates of one coondition.""")
                        
    parser.add_argument("-r", "--resultsDir", type=str,
                        help="""Parent directory of the results directory.""")

                        
    args = parser.parse_args()
    config_path = args.config
    name = args.name
    result_path = args.resultsDir + "/result/"
    expression_paths = args.expressionpath[0]

    return config_path, name, expression_paths, result_path

def main():
    config_path, name, expression_paths, result_path = parser_setup()
    
    fas_lib = library_class.Library(config_path, False)
    
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    movement_path = result_path + "movement/"
    if not os.path.exists(movement_path):
        os.makedirs(movement_path)
        
    expression_path = result_path + "expression/"
    if not os.path.exists(expression_path):
        os.makedirs(expression_path)

    main_comparison_path = result_path + "main_comparison/"
    if not os.path.exists(main_comparison_path):
        os.makedirs(main_comparison_path)

    result_config_path = result_path + "result_config.json"
    if not os.path.isfile(result_config_path):
        result_config_dict = dict()
        result_config_dict["movement_dir"] = movement_path
        result_config_dict["expression_dir"] = expression_path
        result_config_dict["main_comparison_dir"] = main_comparison_path
        result_config_dict["conditions"] = dict()
        with open(result_config_path, 'w') as f:
            json.dump(dict(), f,  indent=4)

    print("Expression extraction commencing...")
    ee.generate_expression_file(fas_lib, 
                                result_config_path, 
                                expression_paths, 
                                expression_path, 
                                name)
    
if __name__ == "__main__":
    main()