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


def extract_all_graph(fas_lib, movement_dict, exempt=["ENSG00000155657"]):
    # need protein ids and movement scores
    fas_ring_dict = dict()
    gene_ids = fas_utility.load_gene_ids_txt(fas_lib.get_config("gene_ids_path"))
    gene_ids = [ gene_id for gene_id in gene_ids if gene_id not in exempt]
    for gene_id in gene_ids:
        if (gene_id in movement_dict["movement"].keys()):
            fas_ring_dict[gene_id] = [movement_dict["movement"][gene_id]["prot_ids"],
                                      movement_dict["movement"][gene_id]["mean_mov"],
                                      movement_dict["movement"][gene_id]["mean_rel_expr"],
                                      movement_dict["movement"][gene_id]["min_mov"],
                                      movement_dict["movement"][gene_id]["max_mov"],
                                      movement_dict["movement"][gene_id]["plus_std_mov"],
                                      movement_dict["movement"][gene_id]["minus_std_mov"],
                                      movement_dict["movement"][gene_id]["intersample_rmsd_mean"]]
    return fas_ring_dict


def generate_comparison(fas_ring_dict_list, conditions, fas_mode, result_config_dict, fas_lib, movement_dict_list, expression_dict_list, movement_paths, expression_paths):
    title = "@".join(conditions) + "_" + fas_mode + "/"
    comparison_path = result_config_dict["main_comparison_dir"] + title
    file_path = comparison_path + "_".join(["result", "@".join(conditions), fas_mode]) + ".tsv"

    fas_ring_dict1, fas_ring_dict2 = fas_ring_dict_list
    output = "\n".join(["!conditions " + conditions[0] + " " + conditions[1],
                       "!FASmode " + fas_mode,
                       "!ResultOrigin " + result_config_dict["result_dir"],
                       "\t".join(["gene_id",
                                 "prot_id",
                                 "mean_mov",
                                 "min_mov",
                                 "max_mov",
                                 "plus_std_mov",
                                 "minus_std_mov",
                                 "1rmsd_smaller_max",
                                 "1rmsd_smaller_mean_plus_std",
                                 "2rmsd_smaller_max",
                                 "2rmsd_smaller_mean_plus_std",
                                 "rmsd",
                                 "max_tsl"])])
    
    for gene_id in fas_ring_dict1.keys():
        ring_list = [ fas_ring_dict1[gene_id], fas_ring_dict2[gene_id]]
        
        fas_ring_dict = poly.prepare_ring_pair_visual(ring_list, gene_id, conditions)
        max_tsl = fas_utility.find_max_tsl(fas_lib, fas_ring_dict["categories"])
        
        output_row_list = [ fas_ring_dict["gene_id"] ]
        output_row_list.append( ";".join( fas_ring_dict["categories"] ) )
        output_row_list.append( ";".join( [ ":".join(entry) for entry in fas_ring_dict["mean_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(entry) for entry in fas_ring_dict["min_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(entry) for entry in fas_ring_dict["max_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(entry) for entry in fas_ring_dict["plus_std_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(entry) for entry in fas_ring_dict["minus_std_movement"] ] ) )
        output_row_list.append( fas_ring_dict["1rmsd_smaller_max"] ) 
        output_row_list.append( fas_ring_dict["1rmsd_smaller_mean_plus_std"] )
        output_row_list.append( fas_ring_dict["2rmsd_smaller_max"] )
        output_row_list.append( fas_ring_dict["2rmsd_smaller_mean_plus_std"] )
        output_row_list.append( fas_ring_dict["rmsd"] )
        output_row_list.append( max_tsl )

        output_row = "\t".join(output_row_list)
        
        output += "\n" + output_row
    with open(file_path, "w") as f:
        f.write(output)

"""// FOR NOW NOT NECESSARY ANYMORE SINCE JULIAN IS BUILDING THE GUI
def sort_by_rmsd(fas_lib, path, flag_more_than_2=True):
    output = "gene_id\tprot_id\t1_in_2_std\t2_in_1_std\tin_std_sum\trmsd\tmax_tsl\n"
    new_path = path[:-4] + "_sorted.tsv"
    with open(path, "r") as f:
        file = f.read()
        comparisons = file.split("\n")
    comparisons = comparisons[4:]
    output = "\n".join(comparisons[0],
                       comparisons[1],
                       comparisons[2],
                       "\t".join("gene_id",
                                 "prot_id",
                                 "1_in_2_std",
                                 "2_in_1_std",
                                 "in_std_sum",
                                 "rmsd",
                                 "max_tsl"))
    comparisons = [ comparison.split("\t") for comparison in comparisons ]
    new_comparisons = []
    for comparison in comparisons:
        if len(comparison) == 12:
            new_comparisons.append(comparison)
    comparisons = new_comparisons
    comparisons = [ comparison[:2] + comparison[7:9] + [int(comparison[-3]),
                                                        float(comparison[-2]),
                                                        int(comparison[-1])] for comparison in comparisons ]
    new_comparisons = []
    for comparison in comparisons:
        if comparison[-2] < 1:
            new_comparisons.append(comparison)
    comparisons = new_comparisons
    comparisons = sorted(comparisons, key = lambda x: (x[-1], -x[-2], -x[3]))
    comparisons = [ comparison[:2] + [str(comparison[-2]), str(comparison[-1])] for comparison in comparisons ]
    comparisons = [ "\t".join(comparison) for comparison in comparisons ]
    file = "\n".join(comparisons)
    file = output + file
    with open(new_path, "w") as f:
        f.write(file)
"""

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
    
    parser.add_argument("-r", "--resultsDir", type=str,
                        help="""Parent directory of the results directory.""")

    parser.add_argument("-s", "--conditions", nargs="+", action="append",
                        help="""Names of the condition that shall be compared as represented in the result_config.json
                        found in the result directory.""")
                        
    parser.add_argument("-l", "--lcr", action="store_true",
                        help="lcr FAS mode shall be used for the comparison.")

    parser.add_argument("-t", "--tmhmm", action="store_true",
                        help="tmhmm FAS mode shall be used for the comparison.")

    parser.add_argument("-a", "--all", action="store_true",
                        help="all domain type FAS mode shall be used for the comparison.")
    
    args = parser.parse_args()
    
    config_path = args.config
    flag_lcr = args.lcr
    flag_tmhmm = args.tmhmm
    flag_all = args.all
    conditions = args.conditions[0]
    result_config_path = args.resultsDir + "/result/result_config.json"

    return config_path, conditions, result_config_path, flag_lcr, flag_tmhmm, flag_all
    
def main():
    config_path, conditions, result_config_path, flag_lcr, flag_tmhmm, flag_all = parser_setup()

    conditions = sorted(conditions)

    if flag_lcr:
        fas_mode = "lcr"
    elif flag_tmhmm:
        fas_mode = "tmhmm"
    elif flag_all:
        fas_mode = "all"

    fas_lib = library_class.Library(config_path, False)
    
    with open(result_config_path, "r") as f: 
        result_config_dict = json.load(f)
    available_conditions = list(result_config_dict["conditions"].keys())
    movement_paths = []
    expression_paths = []
    
    for condition in conditions: 
        if condition not in available_conditions:
            raise Exception("The condition " + condition + " is not available.")
            sys.exit()
        if fas_mode not in result_config_dict["conditions"][condition]["FAS_modes"]:
            raise Exception("The FAS mode " + fas_mode + " is not available for the condition " + condition + ".")
            sys.exit()
        movement_paths.append(result_config_dict["conditions"][condition]["movement_path"][fas_mode])
        expression_paths.append(result_config_dict["conditions"][condition]["expression_path"])
    
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
            
    fas_ring_dict_list = [ extract_all_graph(fas_lib, movement_dict) for movement_dict in movement_dict_list ]
    
    generate_comparison(fas_ring_dict_list,
                        conditions, fas_mode,
                        result_config_dict,
                        fas_lib,
                        movement_dict_list,
                        expression_dict_list,
                        movement_paths,
                        expression_paths)

if __name__ == "__main__":
    main()
   