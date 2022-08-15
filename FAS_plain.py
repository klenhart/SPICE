#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:13:06 2022

@author: chrisbl
"""

import plotly.graph_objects as go
import json
import math
import os

from library_class import Library

from ensembl_access import load_gene_ids_txt

########
def test():
    fas_lib = Library("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/config.tsv", False)
    
    with open(fas_lib.get_config("root_path") + "ARNTL_FAS_graph_dict.json", "r") as f: 
        ARNTL_dict = json.load(f)
    
    with open(fas_lib.get_config("root_path") + "DPF2_FAS_graph_dict.json", "r") as f: 
        DPF2_dict = json.load(f)
    
    DPF2_transcript_ids = set()
    
    DPF2_vector_dict = dict()
    DPF2_fpkm_dict = dict()
    
    for sample in DPF2_dict.keys():
        row = DPF2_dict[sample]
        gene_id, transcripts, fpkm = row.split("\t")
    
        DPF2_fpkm_dict[sample] = float(fpkm)
        DPF2_vector_dict[sample] = dict()
        
        protein_id_expression_tuple_list = [entry.split(":") for entry in transcripts.split(";")]
    
        for protein_id, rel_expr in protein_id_expression_tuple_list:
            DPF2_transcript_ids.add(protein_id)
            DPF2_vector_dict[sample][protein_id] = float(rel_expr)
    
    DPF2_transcript_ids = list(DPF2_transcript_ids)
    
    #######
    
    ARNTL_transcript_ids = set()
    
    ARNTL_vector_dict = dict()
    ARNTL_fpkm_dict = dict()
    
    for sample in ARNTL_dict.keys():
        row = ARNTL_dict[sample]
        gene_id, transcripts, fpkm = row.split("\t")
    
        ARNTL_fpkm_dict[sample] = float(fpkm)
        ARNTL_vector_dict[sample] = dict()
        
        protein_id_expression_tuple_list = [entry.split(":") for entry in transcripts.split(";")]
    
        for protein_id, rel_expr in protein_id_expression_tuple_list:
            ARNTL_transcript_ids.add(protein_id)
            ARNTL_vector_dict[sample][protein_id] = float(rel_expr)
    
    ARNTL_transcript_ids = list(ARNTL_transcript_ids)
    
    #########
    
    #ARNTL
    #Example for 
    make_graph_unscaled("ARNTL", ["WT_A2B","R634Q_A2B"], ARNTL_transcript_ids, ARNTL_vector_dict, fas_lib)
    # Example for unscaled while several samples dont work well without scaling
    make_graph_unscaled("ARNTL", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, fas_lib)
    # Same graph but scaled
    make_graph_scaled("ARNTL", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["WTNC_A2", "R634Q_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["WTNC_A2B", "R634Q_A2B"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["WT_A2", "WTNC_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    make_graph_scaled("ARNTL", ["R634Q_A2", "P633L_A2"], ARNTL_transcript_ids, ARNTL_vector_dict, ARNTL_fpkm_dict,fas_lib)
    
    #DPF2
    #Example for 
    make_graph_unscaled("DPF2", ["WT_A2B","R634Q_A2B"], DPF2_transcript_ids, DPF2_vector_dict, fas_lib)
    # Example for unscaled while several samples dont work well without scaling
    make_graph_unscaled("DPF2", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, fas_lib)
    # Same graph but scaled
    make_graph_scaled("DPF2", ["WT_A2", "WT_A2B","R634Q_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)
    make_graph_scaled("DPF2", ["P633L_A2", "P633L_A2B","R634Q_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)
    make_graph_scaled("DPF2", ["WTNC_A2", "R634Q_A2B"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)
    make_graph_scaled("DPF2", ["WTNC_A2B", "R634Q_A2"], DPF2_transcript_ids, DPF2_vector_dict, DPF2_fpkm_dict,fas_lib)

def make_graph_unscaled(gene_name, fas_graph_list, fas_lib, sample_names):
    
    title = gene_name + "_unscaled_" + "x".join(sample_names)
    
    protein_ids = []
    total_fpkms = []
    sigma_exprs = []
    rel_exprs = []
        
    for prot_ids_list, sigma_expr_list, rel_expr_list, total_fpkm in fas_graph_list:
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
    
    for name in sample_names:
        final_sigma_exprs.append([])
    for i, prot_id in enumerate(protein_ids):
        if i not in delete_list:
            categories.append(prot_id)
            for j, name in enumerate(sample_names):
                final_sigma_exprs[j].append(sigma_exprs[j][i])
    
    total_fpkm_dict = dict()
    for i, name in enumerate(sample_names):
        total_fpkm_dict[name] = total_fpkms[i]
    
    sigma_exprs_dict = dict()
    for i, name in enumerate(sample_names):
        sigma_exprs_dict[name] = final_sigma_exprs[i]

    # Calculate Scales in comparison to largest fpkm.
    
    fig = go.Figure()
    
    for name in sample_names:
        fig.add_trace(go.Scatterpolar(
            r=sigma_exprs_dict[name],
            theta=categories,
            fill="toself",
            name=name
            ))
        
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True
        )
    fig.show()
    fig.write_image(file=fas_lib.get_config("root_path") + "pictures/" + title + ".png")
    return calc_rmsd(list(sigma_exprs_dict.values()))
    
#def make_graph_scaled(gene_name, figure_samples, categories, vector_dict,fpkm_dict , fas_lib):
def make_graph_scaled(gene_name, fas_graph_list, fas_lib, sample_names):
    
    title = gene_name + "_scaled_" + "x".join(sample_names)
    
    protein_ids = []
    total_fpkms = []
    sigma_exprs = []
    rel_exprs = []
        
    for prot_ids_list, sigma_expr_list, rel_expr_list, total_fpkm in fas_graph_list:
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
    
    for name in sample_names:
        final_sigma_exprs.append([])
    for i, prot_id in enumerate(protein_ids):
        if i not in delete_list:
            categories.append(prot_id)
            for j, name in enumerate(sample_names):
                final_sigma_exprs[j].append(sigma_exprs[j][i])
    
    total_fpkm_dict = dict()
    for i, name in enumerate(sample_names):
        total_fpkm_dict[name] = total_fpkms[i]
    
    sigma_exprs_dict = dict()
    for i, name in enumerate(sample_names):
        sigma_exprs_dict[name] = final_sigma_exprs[i]
    
    max_fpkm_name = sample_names[0]
    
    for name in total_fpkm_dict.keys():
        if total_fpkm_dict[name] > total_fpkm_dict[max_fpkm_name]:
            max_fpkm_name = name

    # Calculate Scales in comparison to largest fpkm.
    scales = dict()
    for name in sample_names:
            scales[name] = total_fpkm_dict[name] / total_fpkm_dict[max_fpkm_name]
    
    fig = go.Figure()
    
    for name in sample_names:
        #print(sigma_exprs_dict[name])
        sigma_exprs_dict[name] = scale_list(scales[name], sigma_exprs_dict[name])
        fig.add_trace(go.Scatterpolar(
            r=sigma_exprs_dict[name],
            theta=categories,
            fill="toself",
            name=name
            ))
        
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True
        )
    fig.show()
    fig.write_image(file=fas_lib.get_config("root_path") + "pictures/" + title + ".png")

    return calc_rmsd(list(sigma_exprs_dict.values()))
    
def calc_rmsd(pair_of_lists):
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
    fas_graph = fas_graph[0].split("\t")
    gene_id, isoform_vector, gene_fpkm = fas_graph
    isoform_vector = isoform_vector.split(";")
    isoform_vector = [ entry.split(":") for entry in isoform_vector ]
    isoform_vector = [ [prot_id, float(sigma_expr), float(rel_expr)] for prot_id, sigma_expr, rel_expr in isoform_vector ]
    prot_ids = [ prot_id for prot_id, sigma_expr, rel_expr in isoform_vector ]
    sigma_exprs = [ sigma_expr for prot_id, sigma_expr, rel_expr in isoform_vector ]
    rel_exprs = [ rel_expr for prot_id, sigma_expr, rel_expr in isoform_vector ]
    return prot_ids, sigma_exprs, rel_exprs, float(gene_fpkm)



def generate_comparison(fas_graphs_dict_list, fas_lib, sample_names):
    title = "x".join(sample_names)
    fas_graphs_dict1, fas_graphs_dict2 = fas_graphs_dict_list
    output = "gene_id\tsample_names\tprot_id\tunscaled_expression\tscaled_expression\tunscaled_rmsd\tscaled_rmsd\n"
    for gene_name in fas_graphs_dict1.keys():
        fas_graph_list = [fas_graphs_dict1[gene_name], fas_graphs_dict2[gene_name]]
        title = gene_name + "_" + "x".join(sample_names)
        protein_ids = []
        total_fpkms = []
        sigma_exprs = []
        rel_exprs = []
            
        for prot_ids_list, sigma_expr_list, rel_expr_list, total_fpkm in fas_graph_list:
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
        
        for name in sample_names:
            final_sigma_exprs.append([])
        for i, prot_id in enumerate(protein_ids):
            if i not in delete_list:
                categories.append(prot_id)
                for j, name in enumerate(sample_names):
                    final_sigma_exprs[j].append(sigma_exprs[j][i])
        
        total_fpkm_dict = dict()
        for i, name in enumerate(sample_names):
            total_fpkm_dict[name] = total_fpkms[i]
        
        sigma_exprs_dict = dict()
        scaled_sigma_exprs_dict = dict()
        for i, name in enumerate(sample_names):
            sigma_exprs_dict[name] = final_sigma_exprs[i]
        
        max_fpkm_name = sample_names[0]
        
        for name in total_fpkm_dict.keys():
            if total_fpkm_dict[name] > total_fpkm_dict[max_fpkm_name]:
                max_fpkm_name = name

        # Calculate Scales in comparison to largest fpkm.
        scales = dict()
        for name in sample_names:
            if total_fpkm_dict[max_fpkm_name] == 0:
                scales[name] = 0
            else:
                scales[name] = total_fpkm_dict[name] / total_fpkm_dict[max_fpkm_name]
       
        for name in sample_names:
            scaled_sigma_exprs_dict[name] = scale_list(scales[name], sigma_exprs_dict[name])
            
        scaled_rmsd = calc_rmsd(list(scaled_sigma_exprs_dict.values()))
        unscaled_rmsd = calc_rmsd(list(sigma_exprs_dict.values()))
        
        output_row_list = [gene_name,
                           ";".join(sample_names)]
        unscaled_expr_list = []
        scaled_expr_list = []
        output_row_list.append(";".join(categories))
        for name in sample_names:
            unscaled_expr_list.append(":".join([ str(value) for value in sigma_exprs_dict[name]]))
            scaled_expr_list.append(":".join([ str(value) for value in scaled_sigma_exprs_dict[name]]))
        output_row_list.append(";".join(unscaled_expr_list))
        output_row_list.append(";".join(scaled_expr_list))
        output_row_list.append(str(unscaled_rmsd))
        output_row_list.append(str(scaled_rmsd))

        output_row = "\t".join(output_row_list)
        
        output += output_row + "\n"
    with open(fas_lib.get_config("root_path") + "pictures/" + title + ".tsv", "w") as f:
        f.write(output)
        
def extract_all_graph(fas_lib, path, exempt=[]):
    fas_graphs_dict = dict()
    gene_ids = load_gene_ids_txt(fas_lib.get_config("gene_ids_path"))
    gene_ids = [ gene_id for gene_id in gene_ids if gene_id not in exempt]
    for gene_id in gene_ids:
        fas_graphs_dict[gene_id] = extract_graph(path, gene_id)
    return fas_graphs_dict
        
def get_name(path):
    start = path.index("polygonFAS_") + 11
    end = path.index(".tsv")
    return path[start:end]

def make_graph(fas_lib, gene_id, sample_names, categories, sigma_list, polygon_type):
    title = "x".join(sample_names) + "_" + polygon_type + "_" + gene_id
    fig = go.Figure()
    for i, name in enumerate(sample_names):
        fig.add_trace(go.Scatterpolar(
            r=sigma_list[i],
            theta=categories,
            fill="toself",
            name=name
            ))  
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0,1]
                )),
        showlegend=True
        )
    fig.show()
    fig.write_image(file=fas_lib.get_config("root_path") + "pictures/" + title + ".png")

def visualize_fas_polygon(path1 ,path2, fas_lib, gene_id=None, pre_calc_flag=False):
    name1 = get_name(path1)
    name2 = get_name(path2)
    filepath = fas_lib.get_config("root_path") + "pictures/"
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    if pre_calc_flag:
        polygon_dict1 = extract_all_graph(fas_lib, path1, ["ENSG00000155657"])
        polygon_dict2 = extract_all_graph(fas_lib, path2, ["ENSG00000155657"])
        generate_comparison([polygon_dict1, polygon_dict2], fas_lib, [name1, name2])
    else:
        title = "x".join([name1, name2]) + ".tsv"
        if os.path.isfile(filepath + title):
            found = False
            print("Matching precalculated tsv file found. Decoding the polygons found in it.")
            with open(filepath + title, "r") as f:
                all_polygons = f.read().split("\n")
                all_polygons = [ polygon.split("\t") for polygon in all_polygons ]
                all_polygons = [ entry for entry in all_polygons if entry[0] == gene_id]
            if len(all_polygons) > 0:
                found = True
                gene_id, sample_names, categories, unscaled_expression, scaled_expression, unscaled_rmsd, scaled_rmsd = all_polygons[0]
                sample_names = sample_names.split(";")
                categories = categories.split(";")
                unscaled1, unscaled2 = unscaled_expression.split(";")
                scaled1, scaled2 = scaled_expression.split(";")
                unscaled1 = [ int(entry) for entry in unscaled1.split(":") ]
                unscaled2 = [ int(entry) for entry in  unscaled2.split(":") ]
                scaled1 = [ int(entry) for entry in  scaled1.split(":") ]
                scaled2 = [ int(entry) for entry in  scaled2.split(":") ]
                make_graph(fas_lib, gene_id, sample_names, categories, [unscaled1, unscaled2], "unscaled")
                make_graph(fas_lib, gene_id, sample_names, categories, [scaled1, scaled2], "scaled")

        if not found:
            print("No matching precalculated tsv file found. Generating the visualization from the FASpolygons of the individual samples.")
            polygon1 = extract_graph(path1, gene_id)
            polygon2 = extract_graph(path2, gene_id)
            make_graph_scaled(gene_id, [polygon1, polygon2], fas_lib, [name1, name2])
            make_graph_unscaled(gene_id, [polygon1, polygon2], fas_lib, [name1, name2])

def main():
    pass


if __name__ == "__main__":
    main()

