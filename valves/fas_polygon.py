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
Created on Tue Aug  2 11:13:06 2022
@author: chrisbl
"""

import plotly.graph_objects as go
import math
import os

import valves.fas_utility
    
def calc_rmsd(pair_of_lists):
    if len(pair_of_lists[0]) == 0:
        return 0
    else:
        count = len(pair_of_lists[0])
        difference_list = []
        list_1, list_2 = pair_of_lists
        for i, _ in enumerate(list_1):
            difference_list.append((list_1[i] - list_2[i])**2)
        return math.sqrt(sum(difference_list)/count)

def scale_list(scale_factor, float_list):
    output_list = []
    for value in float_list:
        output_list.append(value * scale_factor)
    return output_list

def extract_graph(path, gene_id):
    with open(path, "r") as f:
        fas_graph = f.read()
    fas_graph = fas_graph.split("\n")
    fas_graph = [ entry for entry in fas_graph if entry.startswith(gene_id)]
    try:
        fas_graph = fas_graph[0].split("\t")
    
        gene_id, isoform_vector, gene_fpkm = fas_graph
    except:
        print(fas_graph)
        print(gene_id)
        print(2.0 * "")
    isoform_vector = isoform_vector.split(";")
    isoform_vector = [ entry.split(":") for entry in isoform_vector ]
    isoform_vector = [ [prot_id, float(sigma_expr), float(rel_expr)] for prot_id, sigma_expr, rel_expr in isoform_vector ]
    prot_ids = [ prot_id for prot_id, sigma_expr, rel_expr in isoform_vector ]
    sigma_exprs = [ sigma_expr for prot_id, sigma_expr, rel_expr in isoform_vector ]
    rel_exprs = [ rel_expr for prot_id, sigma_expr, rel_expr in isoform_vector ]
    return prot_ids, sigma_exprs, rel_exprs, float(gene_fpkm)



def make_graph(fas_lib, gene_id, sample_names, categories, sigma_list, rmsd, filepath):
    fig = go.Figure()
    for i, name in enumerate(sample_names):
        fig.add_trace(go.Scatterpolar(
            r=sigma_list[i] + [sigma_list[i][0]],
            theta=categories + [categories[0]],
            fill="toself",
            name=name
            ))
    fig.add_annotation(x=0, y=-0.2,
                text="RMSD=" + str(rmsd),
                showarrow=False)
    fig.add_annotation(x=0, y=-0.14,
                text=gene_id,
                showarrow=False)
    fig.update_layout(
        autosize=False,
        width=604,
        height=387,
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True,
        paper_bgcolor="LightSteelBlue"
        )
    #fig.show()
    fig.write_image(file=filepath)
    


def prepare_polygon_pair_visual(polygon_list, name_list, gene_id):
    protein_ids = []
    total_fpkms = []
    sigma_exprs = []
    rel_exprs = []
    
    for prot_ids_list, sigma_expr_list, rel_expr_list, total_fpkm in polygon_list:
        protein_ids = prot_ids_list
        total_fpkms.append(total_fpkm)
        sigma_exprs.append(sigma_expr_list)
        rel_exprs.append(rel_expr_list)
    
    delete_list = []
    for i, rel_expr in enumerate(rel_exprs):
        for j, _ in enumerate(rel_expr):
            if all([ entry[j] == 0 for entry in rel_exprs ]):
                delete_list.append(j)
        break
    
    categories = []
    final_sigma_exprs = []
    
    for name in name_list:
        final_sigma_exprs.append([])
    for i, prot_id in enumerate(protein_ids):
        if i not in delete_list:
            categories.append(prot_id)
            for j, name in enumerate(name_list):
                final_sigma_exprs[j].append(sigma_exprs[j][i])
    
    total_fpkm_dict = dict()
    for i, name in enumerate(name_list):
        total_fpkm_dict[name] = total_fpkms[i]
    
    sigma_exprs_dict = dict()
    scaled_sigma_exprs_dict = dict()
    for i, name in enumerate(name_list):
        sigma_exprs_dict[name] = final_sigma_exprs[i]
    
    max_fpkm_name = name_list[0]
    
    for name in total_fpkm_dict.keys():
        if total_fpkm_dict[name] > total_fpkm_dict[max_fpkm_name]:
            max_fpkm_name = name
    
    # Calculate Scales in comparison to largest fpkm.
    scales = dict()
    for name in name_list:
        if total_fpkm_dict[max_fpkm_name] == 0:
            scales[name] = 0
        else:
            scales[name] = total_fpkm_dict[name] / total_fpkm_dict[max_fpkm_name]
    
    for name in name_list:
        scaled_sigma_exprs_dict[name] = scale_list(scales[name], sigma_exprs_dict[name])
        
    scaled_rmsd = calc_rmsd(list(scaled_sigma_exprs_dict.values()))
    unscaled_rmsd = calc_rmsd(list(sigma_exprs_dict.values()))
    
    output_dict = dict()
    output_dict["gene_id"] = gene_id
    output_dict["name_list"] = name_list
    output_dict["categories"] = categories
    output_dict["unscaled_rel_expr"] = list(sigma_exprs_dict.values())
    output_dict["scaled_rel_expr"] =  list(scaled_sigma_exprs_dict.values())
    output_dict["unscaled_rmsd"] = unscaled_rmsd
    output_dict["scaled_rmsd"] = scaled_rmsd
    
    return output_dict

def visualize_fas_polygon(path_list, fas_lib, gene_id, outFormat):    
    path_list = sorted(path_list)
    name_list = [fas_utility.get_name(path_list[0]), fas_utility.get_name(path_list[1])]
    title = "x".join(name_list) + "/"
    pictures_path = fas_lib.get_config("root_path") + "pictures/"
    comparison_path = pictures_path + title
    filepath_scaled = comparison_path + gene_id + "_scaled." + outFormat
    filepath_unscaled = comparison_path + gene_id + "_unscaled." + outFormat
    
    if not os.path.exists(pictures_path):
        os.makedirs(pictures_path)
    if not os.path.exists(comparison_path):
        os.makedirs(comparison_path)
    
    polygon_list = [ extract_graph(path, gene_id) for path in path_list ]   
    polygon_dict = prepare_polygon_pair_visual(polygon_list, name_list, gene_id)
    
    # Unscaled
    make_graph(fas_lib,
               gene_id,
               polygon_dict["name_list"],
               polygon_dict["categories"],
               polygon_dict["unscaled_rel_expr"],
               polygon_dict["unscaled_rmsd"],
               filepath_unscaled)
    
    # Scaled
    make_graph(fas_lib,
               gene_id,
               polygon_dict["name_list"],
               polygon_dict["categories"],
               polygon_dict["scaled_rel_expr"],
               polygon_dict["scaled_rmsd"],
               filepath_scaled)

def main():
    pass


if __name__ == "__main__":
    main()
    