#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:13:06 2022

@author: chrisbl
"""

import plotly.graph_objects as go
import json
import math

from library_class import Library

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

def main():
    fas_lib = Library("/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/config.tsv", False)
    
    r634q = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/FAS_graphs/polygonFAS_R634Q.tsv"
    wtnc = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/FAS_graphs/polygonFAS_WTNC.tsv"
    wt = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/FAS_graphs/polygonFAS_WT.tsv"
    p633l = "/share/project/zarnack/chrisbl/FAS/utility/protein_lib/FAS_library/homo_sapiens/release-107/FAS_graphs/polygonFAS_P633L.tsv"
    
    test_ids_dict = { "RECK" : "ENSG00000122707",
                     "ARNTL" : "ENSG00000133794",
                     "DPF2" : "ENSG00000133884",
                     "PITX2" : "ENSG00000164093",
                     "TPM4" : "ENSG00000167460",
                     "SETD5" : "ENSG00000168137",
                     "UPP1" : "ENSG00000183696",
                     "BCO2" : "ENSG00000197580",     
                     "TPM2" : "ENSG00000198467" }

    r634q_graphs_dict = dict()
    wtnc_graphs_dict = dict()
    
    for key in test_ids_dict.keys():
        r634q_graphs_dict[key] = extract_graph(r634q, test_ids_dict[key])
        wtnc_graphs_dict[key] = extract_graph(wtnc, test_ids_dict[key])

    for gene_name in test_ids_dict.keys():
        rmsd = make_graph_scaled(gene_name,
                                 [r634q_graphs_dict[gene_name],
                                  wtnc_graphs_dict[gene_name]],
                                 fas_lib,
                                 ["R634Q", "WTNC"])
        print(rmsd)
        print(gene_name)
        rmsd = make_graph_unscaled(gene_name,
                                 [r634q_graphs_dict[gene_name],
                                  wtnc_graphs_dict[gene_name]],
                                 fas_lib,
                                 ["R634Q", "WTNC"])
        print(rmsd)
        print(gene_name)
    
    p633l_graphs_dict = dict()
    wt_graphs_dict = dict()
    
    for key in test_ids_dict.keys():
        p633l_graphs_dict[key] = extract_graph(p633l, test_ids_dict[key])
        wt_graphs_dict[key] = extract_graph(wt, test_ids_dict[key])

    for gene_name in test_ids_dict.keys():
        rmsd = make_graph_scaled(gene_name,
                                 [p633l_graphs_dict[gene_name],
                                  wt_graphs_dict[gene_name]],
                                 fas_lib,
                                 ["P633L", "WT"])
        print(rmsd)
        print(gene_name)
        rmsd = make_graph_unscaled(gene_name,
                                 [p633l_graphs_dict[gene_name],
                                  wt_graphs_dict[gene_name]],
                                 fas_lib,
                                 ["P633L", "WT"])
        print(rmsd)
        print(gene_name)
    #test()
    # categories = ["iso1", "iso2", "iso3"]
    
    # fig = go.Figure()
    
    # fig.add_trace(go.Scatterpolar(
    #     r=[0.91, 0.9625, 0.525],
    #     theta=categories,
    #     fill="toself",
    #     name="G_1"))
    
    # fig.update_layout(
    #     polar=dict(
    #         radialaxis=dict(
    #             visible=True,
    #             range=[0,1]
    #             )),
    #     showlegend=True)
    # fig.show()
    # fig.write_image(file="/home/chrisbl/project/FAS_Pipe/diary/exmpl_plot1.png")
    
    # categories = ["iso1", "iso2", "iso3"]
    
    # fig = go.Figure()
    
    # fig.add_trace(go.Scatterpolar(
    #     r=[0.72, 0.825, 0.645],
    #     theta=categories,
    #     fill="toself",
    #     name="F_1"))
    
    # fig.update_layout(
    #     polar=dict(
    #         radialaxis=dict(
    #             visible=True,
    #             range=[0,1]
    #             )),
    #     showlegend=True)
    # fig.show()
    # fig.write_image(file="/home/chrisbl/project/FAS_Pipe/diary/exmpl_plot2.png")
    
    # fig = go.Figure()
    
    # fig.add_trace(go.Scatterpolar(
    #     r=[0.91, 0.9625, 0.525],
    #     theta=categories,
    #     fill="toself",
    #     name="G_1"))
    
    
    # fig.add_trace(go.Scatterpolar(
    #     r=[0.72, 0.825, 0.645],
    #     theta=categories,
    #     fill="toself",
    #     name="F_1"))
    
    # fig.update_layout(
    #     polar=dict(
    #         radialaxis=dict(
    #             visible=True,
    #             range=[0,1]
    #             )),
    #     showlegend=True)
    # fig.show()
    # fig.write_image(file="/home/chrisbl/project/FAS_Pipe/diary/exmpl_plot3.png")

if __name__ == "__main__":
    main()

