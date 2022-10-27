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
Created on Thu Aug 18 11:31:47 2022

@author: chrisbl
"""

import argparse
import os
import json
import sys

import valves.fas_polygon as poly
import valves.library_class as library_class
import valves.fas_utility as fas_utility

"""
Move the generation of a sample comparison file part of the
"--visualize"
option from fas_handler.py to this
"""

def extract_all_graph(fas_lib, ):
    fas_graphs_dict = dict()
    gene_ids = fas_utility.load_gene_ids_txt(fas_lib.get_config("gene_ids_path"))
    gene_ids = [ gene_id for gene_id in gene_ids if gene_id not in exempt]
    for gene_id in gene_ids:
        fas_graphs_dict[gene_id] = poly.extract_graph(path, gene_id)
    return fas_graphs_dict


def generate_comparison(polygon_dict_list, path_list, fas_lib, name_list):
    title = "x".join(name_list) + "/"
    pictures_path = fas_lib.get_config("root_path") + "pictures/"
    comparison_path = pictures_path + title
    file_path = comparison_path + "polygonFAS_" + "x".join(name_list) + ".tsv"
    
    if not os.path.exists(pictures_path):
        os.makedirs(pictures_path)
    if not os.path.exists(comparison_path):
        os.makedirs(comparison_path)
    
    polygon_dict1, polygon_dict2 = polygon_dict_list
    output = "gene_id\tsample_names\tprot_id\tunscaled_expression\tscaled_expression\tunscaled_rmsd\tscaled_rmsd\tmax_tsl\n"
    for gene_name in polygon_dict1.keys():
        polygon_list = [ polygon_dict1[gene_name], polygon_dict2[gene_name]]
        polygon_dict = poly.prepare_polygon_pair_visual(polygon_list, name_list, gene_name)
        
        max_tsl = fas_utility.find_max_tsl(fas_lib, polygon_dict["categories"])
        
        output_row_list = [polygon_dict["gene_id"],
                           ";".join(polygon_dict["name_list"])]
        unscaled_expr_list = []
        scaled_expr_list = []
        output_row_list.append(";".join(polygon_dict["categories"]))
        unscaled_expr_list = [ ":".join([ str(value) for value in val_list]) for val_list in polygon_dict["unscaled_rel_expr"]]
        scaled_expr_list = [ ":".join([ str(value) for value in val_list]) for val_list in polygon_dict["scaled_rel_expr"]]

        output_row_list.append(";".join(unscaled_expr_list))
        output_row_list.append(";".join(scaled_expr_list))
        output_row_list.append(str(polygon_dict["unscaled_rmsd"]))
        output_row_list.append(str(polygon_dict["scaled_rmsd"]))
        output_row_list.append(max_tsl)

        output_row = "\t".join(output_row_list)
        
        output += output_row + "\n"
    with open(file_path, "w") as f:
        f.write(output)
    return file_path

def sort_by_rmsd(fas_lib, path, flag_more_than_2=True):
    output = "gene_id\tsample_names\tprot_id\tunscaled_expression\tscaled_expression\tunscaled_rmsd\tscaled_rmsd\tmax_tsl\n"
    new_path = path[:-4] + "_sorted.tsv"
    with open(path, "r") as f:
        file = f.read()
        comparisons = file.split("\n")
    comparisons = comparisons[1:]
    comparisons = [ comparison.split("\t") for comparison in comparisons ]
    new_comparisons = []
    for comparison in comparisons:
        if len(comparison) == 8:
            new_comparisons.append(comparison)
    comparisons = new_comparisons
    comparisons = [ comparison[:-3] + [float(comparison[-3]), float(comparison[-2]), int(comparison[-1])] for comparison in comparisons ]
    new_comparisons = []
    for comparison in comparisons:
        if comparison[-3] < 1 and comparison[-3] > 0 and comparison[-3] != comparison[-2]:
            new_comparisons.append(comparison)
    comparisons = new_comparisons
    comparisons = sorted(comparisons, key = lambda x: (x[-1], -x[-3], -x[-2]))
    comparisons = [ comparison[:-3] + [str(comparison[-3]), str(comparison[-2]), str(comparison[-1])] for comparison in comparisons ]
    comparisons = [ "\t".join(comparison) for comparison in comparisons ]
    file = "\n".join(comparisons)
    file = output + file
    with open(new_path, "w") as f:
        f.write(file)


def parser_setup():
    """
    Reads the parser input.

    Returns
    -------
    """  
    
    #Setting up parser:
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--config", type=str,
                        help="Path to a config file of a library.")

    parser.add_argument("-s", "--conditions", nargs="+", action="append",
                        help="""Names of the condition that shall be compared as represented in the result_config.json
                        found in the result directory.""")
                        
    parser.add_argument("-l", "--lcr", type=str,
                        help="lcr FAS mode shall be used for the comparison.")

    parser.add_argument("-t", "--tmhmm", type=str,
                        help="tmhmm FAS mode shall be used for the comparison.")

    parser.add_argument("-a", "--all", type=str,
                        help="all domain type FAS mode shall be used for the comparison.")
    
    args = parser.parse_args()
    
    config_path = args.config
    flag_lcr = args.lcr
    flag_tmhmm = args.tmhmm
    flag_all = args.all
    conditions = args.conditions[0]

    return config_path, conditions, flag_lcr, flag_tmhmm, flag_all
    
def main():
    config_path, conditions, flag_lcr, flag_tmhmm, flag_all = parser_setup()

    conditions = sorted(conditions)

    if flag_lcr:
        fas_mode = "lcr"
    elif flag_tmhmm:
        fas_mode = "tmhmm"
    elif flag_all:
        fas_mode = "all"

    fas_lib = library_class.Library(config_path, False)
    
    result_config_path = fas_lib.get_config("result_config")
    
    with open(result_config_path, "r") as f: 
        result_config_dict = json.load(f)
    available_conditions = list(result_config_dict.keys())
    movement_paths = []
    expression_paths = []
    
    for condition in conditions: 
        if condition not in available_conditions:
            raise Exception("The condition " + condition + " is not available.")
            sys.exit()
        if fas_mode not in result_config_dict[condition]["FAS_modes"]:
            raise Exception("The FAS mode " + fas_mode + " is not available for the condition " + condition + ".")
            sys.exit()
        movement_paths.append(result_config_dict[condition]["movement_path"][fas_mode])
        expression_paths.append(result_config_dict[condition]["expression_path"])
    
    movement_dict_list = []
    expression_dict_list = []
    for path in movement_paths:
        with open(path, "r") as f:
            movement_dict_list.append(json.load(f))
    for path in expression_paths:
        with open(path, "r") as f:
            expression_dict_list.append(json.load(f))
    
    for i, movement_dict in enumerate(movement_dict_list, start=-len(movement_dict_list)):
        if condition[abs(i)] in movement_dict["compared_with"]:
            raise Exception( condition[0] + " and " + condition[1] + " using FAS mode " + fas_mode + " were already compared.")
            sys.exit()
    
    total_expression_dict = 
    
    
            
        
    
        
    
    fas_graphs_dict_list = [ extract_all_graph(fas_lib, path) for path in path_list ]
    name_list = [ prefix + fas_utility.get_name(path) for path in path_list ]
    file_path = generate_comparison(fas_graphs_dict_list, path_list, fas_lib, name_list)
    sort_by_rmsd(fas_lib, file_path, flag_more_than_2=True)


if __name__ == "__main__":
    main()
   