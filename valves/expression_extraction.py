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
Created on Fri Aug  5 12:43:02 2022

@author: chrisbl
"""

import pyranges as pr
import os
import json
import numpy as np

from valves.fas_utility import tsv_to_tuple_list
from valves.fas_utility import longest_common_prefix


def load_expression_gtf(expression_path, flag_filter_unknown=True, flag_remove_transcript_prefix=True):
    gtf = pr.read_gtf(expression_path, as_df=True)
    gtf = [ [gene_id, transcript_id, float(fpkm)] for gene_id,
                                      transcript_id,
                                      fpkm,
                                      feature in zip(gtf['gene_id'],
                                                    gtf['transcript_id'],
                                                    gtf['FPKM'],
                                                    gtf['Feature']) if feature == "transcript"]
    if flag_remove_transcript_prefix:
        gtf = [ [entry[0], entry[1].split(":")[-1], entry[2]] for entry in gtf ]
    if flag_filter_unknown:
        gtf = [ entry for entry in gtf if entry[1].startswith("ENS") ]
    gtf = [ [entry[0].split(".")[0], entry[1].split(".")[0], entry[2]] for entry in gtf ]
    return gtf


def load_protein_coding_ids(path):
    #fas_lib.get_config("protein_coding_ids_path")
    protein_coding_ids = tsv_to_tuple_list(path)
    return protein_coding_ids


def filter_protein_coding(expression_data, protein_coding_ids):
    # FIlter by protein coding
    protein_coding_transcript_ids = [ transcript_id for gene_id, protein_id, transcript_id, tsl, tag in protein_coding_ids ]
    expression_data = [ [gene_id, transcript_id, fpkm] for gene_id, transcript_id, fpkm in expression_data if transcript_id in protein_coding_transcript_ids]
    return expression_data


def fix_expression_ids( expression_data, protein_coding_ids):
    translate_trans_to_gene_dict = dict()
    for gene_id, protein_id, transcript_id, tsl, tag in protein_coding_ids:
        translate_trans_to_gene_dict[transcript_id] = (gene_id, protein_id)
    
    new_expression_data = []
    for gene_id, transcript_id, fpkm in expression_data:
        new_gene_id, new_prot_id = translate_trans_to_gene_dict[transcript_id]
        # If the gene_id is missing.
        if not gene_id.startswith("ENS"):
            new_expression_data.append([new_gene_id, new_prot_id, transcript_id, fpkm])
        # If the gene_id is not missing.
        else:
            # If the gene_id does not belong to the transcript_id
            if new_gene_id != gene_id:
                print(gene_id, "incorrectly associated with", transcript_id)
                print("Dropping from Expression data...")
            # If all is in order
            else:
                new_expression_data.append([gene_id, new_prot_id, transcript_id, fpkm])
    return new_expression_data


def filter_non_expressed(expression_data, fpkm_threshold=0):
    """
    Not in use anymore

    """
    return [ [gene_id, transcript_id, fpkm] for gene_id, transcript_id, fpkm in expression_data if fpkm > fpkm_threshold ]


def make_isoform_protein_id_dict(expression_data):
    isoform_protein_id_dict = dict()
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        isoform_protein_id_dict[gene_id] = []
    for gene_id, protein_id, transcript_id, fpkm in expression_data:
        isoform_protein_id_dict[gene_id].append(protein_id)
    return isoform_protein_id_dict


def isoform_count_distribution_check(isoform_protein_id_dict):
    isoform_count_list= []
    three_or_more_count = 0
    for isoform_list in isoform_protein_id_dict.values():
        isoform_count_list.append(len(isoform_list))
        if len(isoform_list) >= 3:
            three_or_more_count += 1
    isoform_count_list.sort()
    three_or_more_ratio = three_or_more_count / len(isoform_protein_id_dict.values())
    median = isoform_count_list[int(len(isoform_count_list)/2)]
    mean = sum(isoform_count_list) / len(isoform_count_list)
    return median, mean, three_or_more_ratio


def check_for_same_order(m1, m2):
    flag_check = True
    for i, entry in enumerate(m1):
        gene_id_1, transcript_id_1, fpkm_1 = entry
        gene_id_2, transcript_id_2, fpkm_2 = m2[i]
        if gene_id_1 != gene_id_2 or transcript_id_1 != transcript_id_2:
            flag_check = False
            print(i)
            print(gene_id_1, " x ", gene_id_2)
            print(transcript_id_1, " x ", transcript_id_2)
            break
    return flag_check


def filter_exempt_genes(expression_data, exempt_genes):
    if len(exempt_genes) == 0:
        pass
    else:
        expression_data = [ entry for entry in expression_data if entry[0] not in exempt_genes ]
    return expression_data


def join_expression(expression_path_list, protein_coding_path, exempt_genes):
    protein_coding_ids = load_protein_coding_ids(protein_coding_path)
    protein_coding_ids = [ entry for entry in protein_coding_ids if entry[0] not in exempt_genes ]
    expression_dict = dict()
    prot_to_gene_dict = dict()
    for gene_id, prot_id, transcript_id, tsl, tag in protein_coding_ids:
        expression_dict[prot_id] = []
        prot_to_gene_dict[prot_id] = gene_id
    # Extract all the expression from the gtf files and filter out everything unnecessary.
    for expression_path in expression_path_list:
        expression_data = load_expression_gtf(expression_path)
        expression_data = filter_protein_coding(expression_data, protein_coding_ids)
        expression_data = fix_expression_ids(expression_data, protein_coding_ids)
        expression_data = filter_exempt_genes(expression_data, exempt_genes)
        for gene_id, prot_id, transcript_id, fpkm in expression_data:
            expression_dict[prot_id].append(fpkm)
    isoforms_dict = dict()
    for prot_id in expression_dict.keys():
        gene_id = prot_to_gene_dict[prot_id]
        if gene_id in isoforms_dict.keys():
            isoforms_dict[gene_id].append(prot_id)
        else:
            isoforms_dict[gene_id] = [prot_id]
    return expression_dict, isoforms_dict

def load_dist_matrix(fas_lib, flag_lcr, flag_tmhmm, flag_all):
    if flag_lcr:
        distance_master_name = "fas_lcr.json"
        distance_master_path = fas_lib.get_config("fas_lcr_path")
    elif flag_tmhmm:
        distance_master_name = "fas_tmhmm.json"
        distance_master_path = fas_lib.get_config("fas_tmhmm_path")
    elif flag_all:
        distance_master_name = "fas.json"
        distance_master_path = fas_lib.get_config("fas_all_path")
    if distance_master_name in os.listdir(fas_lib.get_config("root_path")):
        with open(fas_lib.get_config("root_path") + distance_master_name, "r") as f: 
            dist_matrix_dict = json.load(f)
    else:
        with open(distance_master_path, "r") as f:
            distance_master = f.read()
        distance_master = distance_master.split("\n")   
        distance_master = distance_master[1:]
        if distance_master[-1] == "":
            distance_master = distance_master[:-1]
        distance_master = [ entry.split("\t") for entry in distance_master ]
        dist_matrix_dict = dict()
        for seed, tax_id, query, fas_1, fas_2 in distance_master:
            #fas_score = (float(fas_1) + float(fas_2)) / 2
            fas_1 = float(fas_1)
            fas_2 = float(fas_2)
            gene_id, prot_id_1, tax_id = seed.split("|")
            gene_id, prot_id_2, tax_id = query.split("|")
            if gene_id not in dist_matrix_dict.keys():
                dist_matrix_dict[gene_id] = {prot_id_1 : { prot_id_2 : fas_2} , prot_id_2 : { prot_id_1 : fas_1} }
            else:
                if prot_id_1 not in dist_matrix_dict[gene_id].keys():
                    dist_matrix_dict[gene_id][prot_id_1] = { prot_id_2 : fas_2}
                else:
                    dist_matrix_dict[gene_id][prot_id_1][prot_id_2] = fas_2
                if prot_id_2 not in dist_matrix_dict[gene_id].keys():
                    dist_matrix_dict[gene_id][prot_id_2] = { prot_id_1 : fas_1}
                else:
                    dist_matrix_dict[gene_id][prot_id_2][prot_id_1] = fas_1
        with open(fas_lib.get_config("root_path") + distance_master_name, 'w') as f:
            json.dump(dist_matrix_dict, f,  indent=4)
    return dist_matrix_dict


def calculate_movement(fas_dist_matrix, expression_vector, gene_id, prot_ids):
    rel_expressions = calculate_relative_expression(expression_vector)
    relative_expression_dict = dict(zip(prot_ids, rel_expressions))
    
    # Extract the distance matrix containing th FAS scores
    dist_matrix = fas_dist_matrix[gene_id]
    movement_dict = dict()
    
    # Calculate the Sigma values.
    for seed_id in prot_ids:
        movement_dict[seed_id] = 0
        for query_id in prot_ids:
            movement_dict[seed_id] += (dist_matrix[seed_id][query_id] * relative_expression_dict[query_id])
    
    movement_list = []
    for prot_id in prot_ids:
        movement_list.append(movement_dict[prot_id])
    return movement_list, rel_expressions


def calculate_relative_expression(expression_vector):
    total = sum(expression_vector)
    if total == 0:
        return expression_vector * 0
    return expression_vector / total
    


def generate_expression_file(fas_lib, result_config_path, expression_paths, expression_path, name):
    release_num = fas_lib.get_config("release_num")
        
    prefix = longest_common_prefix(expression_paths)
    name_cutoff_index = len(prefix)
    expression_names = [ path[name_cutoff_index:] for path in expression_paths ]
        
    expression_dict, isoforms_dict = join_expression(expression_paths,
                                                     fas_lib.get_config("protein_coding_ids_path"),
                                                     ["ENSG00000155657"])
    expression_dict = dict()
    expression_dict["condition"] = name
    expression_dict["replicates"] = expression_names
    expression_dict["normalization"] = "FPKM"
    expression_dict["prefix"] = prefix

    expression_dict["expression"] = dict()

    for gene_id in isoforms_dict.keys():
        expression_dict["expression"][gene_id] = dict()
        for prot_id in isoforms_dict[gene_id]:
            for i, replicate in enumerate(expression_names):
                expression_dict["expression"][gene_id][replicate][prot_id] = expression_dict[prot_id][i]
        
    expression_file_path = expression_path + "_".join(["expression", name, "ENSv" + release_num]) +  ".json"
    with open(expression_file_path, "w") as f:
        json.dump(expression_dict, f,  indent=4)

    with open(result_config_path, "r") as f: 
        result_config_dict = json.load(f)
        
    result_config_dict[name] = dict()
    result_config_dict[name]["prefix"] = prefix
    result_config_dict[name]["replicates"] = expression_names
    result_config_dict[name]["normalization"] = "FPKM"
    result_config_dict[name]["FAS_modes"] = []
    result_config_dict[name]["species"] = fas_lib.get_config("species")
    result_config_dict[name]["release"] = fas_lib.get_config("release")
    result_config_dict[name]["lib_config"] = fas_lib.get_config("self_path")
    result_config_dict[name]["expression_path"] = expression_file_path
    result_config_dict[name]["movement_path"]["all"] = None
    result_config_dict[name]["movement_path"]["lcr"] = None
    result_config_dict[name]["movement_path"]["tmhmm"] = None
    result_config_dict[name]["compared_with"] = []

    with open(result_config_path, 'w') as f:
        json.dump(result_config_dict, f,  indent=4)

def generate_movement_file(fas_lib, name_path, flag_lcr, flag_tmhmm, flag_all):
    result_config_path = fas_lib.get_config("result_config")
    movement_path = fas_lib.get_config("movement")
    fas_dist_matrix = load_dist_matrix(fas_lib, flag_lcr, flag_tmhmm, flag_all)
    
    if flag_lcr:
        fas_mode = "lcr"
    elif flag_tmhmm:
        fas_mode = "tmhmm"
    elif flag_all:
        fas_mode = "all"
     
    with open(name_path, "r") as f:
        name_list = f.read()
    name_list = name_list.split("\n")
    
    with open(result_config_path, "r") as f: 
        result_config_dict = json.load(f)
        
    conditions = list(result_config_dict.keys())
    path_list = []
    for name in name_list:
        if name in conditions:
            if fas_mode in result_config_dict[name]["FAS_modes"]:
                print("""Movement using FAS_""" + fas_mode, """distances has already been calculated.
Check this file:""", result_config_dict[name]["movement_path"][fas_mode])
            else:
                path_list.append(result_config_dict[name]["expression_path"])
        else:
            print("No condition of the name", name, "processed yet. Will skip movement generation for this one.")
    
    for expression_path in path_list:
        with open(expression_path, "r") as f: 
            expression_dict = json.load(f)
        condition = expression_dict["condition"]
        replicates = expression_dict["replicates"]
        
        output_dict = dict()

        output_dict["condition"] = condition
        output_dict["replicates"] = replicates
        output_dict["normalization"] = "FPKM"
        output_dict["movement"] = dict()
        output_dict["compared_with"] = []
        
        for gene_id in list(expression_dict["expression"].keys()):
            expr_matrix = np.array([])
            prot_ids = []
            for replicate in replicates:
                expr_matrix = np.append(expr_matrix, [])
                for prot_id in list(expression_dict["expression"][gene_id][replicate].keys()):
                    prot_ids.append(prot_id)
                    expr_matrix[-1] = np.append(expr_matrix[-1], expression_dict["expression"][gene_id][replicate][prot_id])
            t_expr_matrix = np.transpose(expr_matrix)
            
            min_list = np.array([ min(entry) for entry in t_expr_matrix ])
            min_movement, min_relative_expression = calculate_movement(fas_dist_matrix, min_list, gene_id, prot_ids)
            
            max_list = np.array([ max(entry) for entry in t_expr_matrix ])
            max_movement, max_relative_expression = calculate_movement(fas_dist_matrix, max_list, gene_id, prot_ids)
            
            mean_list = np.array([ np.mean(entry) for entry in t_expr_matrix ])
            mean_movement, mean_relative_expression = calculate_movement(fas_dist_matrix, mean_list, gene_id, prot_ids)
            
            standard_deviation_list = np.array([ np.std(entry) for entry in t_expr_matrix ])
            standard_deviation_movement, standard_deviation_relative_expression = calculate_movement(fas_dist_matrix, standard_deviation_list, gene_id, prot_ids)
            
            plus_std_list = standard_deviation_list + mean_list
            plus_std_movement, plus_std_relative_expression = calculate_movement(fas_dist_matrix, plus_std_list, gene_id, prot_ids)
            
            minus_std_list = standard_deviation_list - mean_list
            minus_std_movement, minus_std_relative_expression = calculate_movement(fas_dist_matrix, minus_std_list, gene_id, prot_ids)
            
            gene_dict = dict()
            gene_dict["prot_ids"] = prot_ids
            gene_dict["min_mov"] = min_movement
            gene_dict["min_rel_expr"] = min_relative_expression
            gene_dict["max_mov"] = max_movement
            gene_dict["max_rel_expr"] = max_relative_expression
            gene_dict["mean_mov"] = mean_movement
            gene_dict["mean_rel_expr"] = mean_relative_expression
            gene_dict["plus_std_mov"] = plus_std_movement
            gene_dict["plus_std_rel_expr"] = plus_std_relative_expression
            gene_dict["minus_std_mov"] = plus_std_movement
            gene_dict["minus_std_rel_expr"] = plus_std_relative_expression
            
            output_dict["movement"][gene_id] = gene_dict

        result_config_dict[condition]["movement_path"][fas_mode] = movement_path + "_".join("movement", condition, fas_mode) + ".json"
        result_config_dict[condition]["FAS_modes"].append(fas_mode)
        
        with open(result_config_dict[condition]["movement_path"][fas_mode], 'w') as f:
            json.dump(output_dict, f,  indent=4)
        
        with open(result_config_path, 'w') as f:
            json.dump(result_config_dict, f,  indent=4)
        

def main():
    pass

if __name__ == "__main__":
    main()
