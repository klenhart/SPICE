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
import time
import random

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


def prepare_ring_pair_visual(ring_list, gene_id, conditions):
    protein_ids = []
    mean_movs = []
    mean_rel_exprs = []
    min_movs = []
    max_movs = []
    plus_std_movs = []
    minus_std_movs = []
    intersample_rmsd = []

    for prot_ids_list, mean_mov_list, mean_rel_expr_list, min_mov_list, max_mov_list, plus_std_mov_list, minus_std_mov_list, intersample_rmsd_mean in ring_list:
        protein_ids = prot_ids_list
        mean_movs.append(mean_mov_list)
        mean_rel_exprs.append(mean_rel_expr_list)
        min_movs.append(min_mov_list)
        max_movs.append(max_mov_list)
        plus_std_movs.append(plus_std_mov_list)
        minus_std_movs.append(minus_std_mov_list)
        intersample_rmsd.append(intersample_rmsd_mean)

    
    delete_list = []
    for i, mean_rel_expr in enumerate(mean_rel_exprs):
        for j, _ in enumerate(mean_rel_expr):
            if all([ entry[j] == 0 for entry in mean_rel_exprs ]):
                delete_list.append(j)
        break
    
    categories = []
    final_mean_movs = []
    final_min_movs = []
    final_max_movs = []
    final_plus_std_movs = []
    final_minus_std_movs = []
    
    for i in range(2):
        final_mean_movs.append([])
        final_min_movs.append([])
        final_max_movs.append([])
        final_plus_std_movs.append([])
        final_minus_std_movs.append([])
    for i, prot_id in enumerate(protein_ids):
        if i not in delete_list:
            categories.append(prot_id)
            for j in range(2):
                final_mean_movs[j].append(mean_movs[j][i])
                final_min_movs[j].append(min_movs[j][i])
                final_max_movs[j].append(max_movs[j][i])
                final_plus_std_movs[j].append(plus_std_movs[j][i])
                final_minus_std_movs[j].append(minus_std_movs[j][i])

    rmsd = fas_utility.calc_rmsd(final_mean_movs)
    
    rmsd_max = fas_utility.calc_rmsd(final_max_movs)
    rmsd_plus_std = fas_utility.calc_rmsd(final_plus_std_movs)

    rmsd_max_smaller_1 = False if rmsd_max < intersample_rmsd[0] else True
    rmsd_mean_plus_std_smaller_1 = False if rmsd_plus_std < intersample_rmsd[0] else True
    rmsd_max_smaller_2 = False if rmsd_max < intersample_rmsd[1] else True
    rmsd_mean_plus_std_smaller_2 = False if rmsd_plus_std < intersample_rmsd[1] else True
       
    
    output_dict = dict()
    output_dict["gene_id"] = gene_id
    output_dict["categories"] = categories
    output_dict["mean_movement"] = final_mean_movs
    output_dict["min_movement"] = final_min_movs
    output_dict["max_movement"] = final_max_movs
    output_dict["plus_std_movement"] = final_plus_std_movs
    output_dict["minus_std_movement"] = final_minus_std_movs
    
    output_dict["rmsd_max_smaller_1"] = str(int(rmsd_max_smaller_1))
    output_dict["rmsd_mean_plus_std_smaller_1"] = str(int(rmsd_mean_plus_std_smaller_1))
    output_dict["rmsd_max_smaller_2"] = str(int(rmsd_max_smaller_2))
    output_dict["rmsd_mean_plus_std_smaller_2"] = str(int(rmsd_mean_plus_std_smaller_2))
        
    output_dict["rmsd"] = rmsd
    
    return output_dict


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
                                 "rmsd_max_smaller_1",
                                 "rmsd_mean_plus_std_smaller_1",
                                 "rmsd_max_smaller_2",
                                 "rmsd_mean_plus_std_smaller_2",
                                 "rmsd",
                                 "max_tsl"])])
    
    for gene_id in fas_ring_dict1.keys():
        ring_list = [ fas_ring_dict1[gene_id], fas_ring_dict2[gene_id]]
        
        fas_ring_dict = prepare_ring_pair_visual(ring_list, gene_id, conditions)
        max_tsl = fas_utility.find_max_tsl(fas_lib, fas_ring_dict["categories"])
        
        output_row_list = [ fas_ring_dict["gene_id"] ]
        output_row_list.append( ";".join( fas_ring_dict["categories"] ) )
        output_row_list.append( ";".join( [ ":".join(str(entry)) for entry in fas_ring_dict["mean_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(str(entry)) for entry in fas_ring_dict["min_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(str(entry)) for entry in fas_ring_dict["max_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(str(entry)) for entry in fas_ring_dict["plus_std_movement"] ] ) )
        output_row_list.append( ";".join( [ ":".join(str(entry)) for entry in fas_ring_dict["minus_std_movement"] ] ) )
        output_row_list.append( fas_ring_dict["rmsd_max_smaller_1"] ) 
        output_row_list.append( fas_ring_dict["rmsd_mean_plus_std_smaller_1"] )
        output_row_list.append( fas_ring_dict["rmsd_max_smaller_2"] )
        output_row_list.append( fas_ring_dict["rmsd_mean_plus_std_smaller_2"] )
        output_row_list.append( str(fas_ring_dict["rmsd"]) )
        output_row_list.append( str(max_tsl) )

        output_row = "\t".join(output_row_list)
        
        output += "\n" + output_row
    
    if not os.path.exists(comparison_path):
        os.makedirs(comparison_path)    
    with open(file_path, "w") as f:
        f.write(output)
        
    with open(result_config_dict["result_dir"] + "/result.config.json", "r") as f:
        result_dict = json.load(f)
    result_dict["conditions"][conditions[0]]["compared_with"].append(conditions[1])
    result_dict["conditions"][conditions[1]]["compared_with"].append(conditions[0])
    with open(result_config_dict["result_dir"] + "/result.config.json", "w") as f:
        json.dump(result_dict, f, indent=4)

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
    time.sleep(random.randint(10, 50))
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
   