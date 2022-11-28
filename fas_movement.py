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
Created on Thu Aug 18 11:30:47 2022

@author: chrisbl
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
                        help="Path to a config file of a FAS library.")

    parser.add_argument("-r", "--resultsDir", type=str,
                        help="""Parent directory of the results directory.""")
                        
    parser.add_argument("-n", "--condition", type=str,
                        help="""Name of a condition already present in the results/expression directory
                        of the results. Movement file will be calculated for this condition.""")
    
    parser.add_argument("-t", "--tmhmm", action="store_true",
                        help="tmhmm FAS mode shall be used for the movement calculation.")

    parser.add_argument("-l", "--lcr", action="store_true",
                        help="lcr FAS mode shall be used for the movement calculation.")

    parser.add_argument("-a", "--all", action="store_true",
                        help="all domain type FAS mode shall be used for the movement calculation.")

                        
    args = parser.parse_args()
    config_path = args.config
    condition = args.condition
    result_config_path = args.resultsDir + "/result/result_config.json"
    
    flag_lcr = args.lcr
    flag_tmhmm = args.tmhmm
    flag_all = args.all

    return config_path, result_config_path, condition, flag_lcr, flag_tmhmm, flag_all

def main():
    """
    Returns
    -------
    """
    config_path, result_config_path, condition, flag_lcr, flag_tmhmm, flag_all = parser_setup()

    fas_lib = library_class.Library(config_path, False)
    print("Movement calculation commencing...")
    ee.generate_movement_file(fas_lib, result_config_path, condition, flag_lcr, flag_tmhmm, flag_all)
    
if __name__ == "__main__":
    main()
   